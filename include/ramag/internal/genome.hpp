#pragma once
#include "ramag/pairwise_core.hpp"
#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <exception>
#include <map>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <variant>
#ifdef _OPENMP
#include <omp.h>
#endif
namespace ramag::pairwise_detail {
inline constexpr size_t kDpWindow = 5000;
using int_t=int64_t; using uint_t=uint64_t;
using Score_t=int64_t;
using Coord_t=uint64_t; using Length_t=uint64_t; using ChrIndex=uint32_t;
using CigarUnit=uint64_t; using Cigar_t=std::vector<CigarUnit>;
using Strand = ramag::Strand;
constexpr Strand FORWARD = Strand::Forward;
constexpr Strand REVERSE = Strand::Reverse;
using Match = Seed;
extern thread_local PairwiseStatistics* active_statistics;
struct StatisticsScope {
    PairwiseStatistics* previous;
    explicit StatisticsScope(PairwiseStatistics& value) : previous(active_statistics) { active_statistics=&value; }
    ~StatisticsScope() { active_statistics=previous; }
};
struct KswCallTimer {
    bool global;
    std::chrono::steady_clock::time_point begin=std::chrono::steady_clock::now();
    explicit KswCallTimer(bool is_global):global(is_global) {}
    ~KswCallTimer() {
        if(!active_statistics)return;
        const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();
        if(global){++active_statistics->global_ksw_calls;active_statistics->global_ksw_seconds+=seconds;}
        else{++active_statistics->endpoint_ksw_calls;active_statistics->endpoint_ksw_seconds+=seconds;}
    }
};
using MatchVec=std::vector<Match>; using MatchCluster=MatchVec;
using MatchClusterVec=std::vector<MatchCluster>;
using MatchClusterVecPtr=std::shared_ptr<MatchClusterVec>;
inline Coord_t start1(const Match& m){return m.reference_begin;}
inline Coord_t start2(const Match& m){return m.query_begin;}
inline Length_t len1(const Match& m){return m.length;}
inline Length_t len2(const Match& m){return m.length;}
inline int_t diag(const Match& m){return static_cast<int_t>(start2(m))-static_cast<int_t>(start1(m));}
inline int_t diag_reverse(const Match& m){return -static_cast<int_t>(start2(m))-static_cast<int_t>(start1(m));}



inline void releaseCluster(MatchVec& v){MatchVec().swap(v);}
namespace SequenceViews {
struct SequenceManager {
 std::map<ChrIndex,const SequenceRecord*> records;
 explicit SequenceManager(std::span<const SequenceRecord> data){for(const auto& r:data) records.emplace(r.numeric_id,&r);}
 bool tryGetContiguousSubSequence(ChrIndex id,Coord_t begin,Coord_t len,std::string_view& view) const {
  const auto& seq=records.at(id)->bases;
  if(begin>seq.size() || len>seq.size()-begin) throw std::out_of_range("pairwise sequence slice exceeds contig");
  view=std::string_view(seq).substr(begin,len);return true;
 }
 void getSubSequenceInto(ChrIndex id,Coord_t begin,Coord_t len,std::string& out) const {std::string_view v;tryGetContiguousSubSequence(id,begin,len,v);out.assign(v);}
};
using ManagerVariant=std::variant<std::unique_ptr<SequenceManager>>;
}
struct Anchor {
    ChrIndex reference_id;
    Coord_t  reference_begin;
    Length_t ref_len;
    ChrIndex query_id;
    Coord_t  query_begin;
    Length_t qry_len;

    Strand strand;
    uint_t alignment_length{}; // 对齐长度
    uint_t aligned_base{};
    Cigar_t cigar{};                 // 对齐的 CIGAR 字符串

    bool ref_selected = false;
    bool qry_selected = false;
    bool is_linked = false;

    Anchor() = default;

    Anchor(ChrIndex ref_chr, Coord_t reference_begin, Length_t ref_len,
        ChrIndex qry_chr, Coord_t query_begin, Length_t qry_len,
        Strand strand, uint_t align_len, uint_t aligned_base, Cigar_t cigar_str)
        : reference_id(ref_chr),
        reference_begin(reference_begin),
        ref_len(ref_len),
        query_id(qry_chr),
        query_begin(query_begin),
        qry_len(qry_len),
        strand(strand),
        alignment_length(align_len),
        aligned_base(aligned_base),
        cigar(std::move(cigar_str)) {
    }


};

using AnchorVec = std::vector<Anchor>;
class UnionFind
{
public:
    /// 构造并初始化为 n 个独立元素
    explicit UnionFind(std::size_t n = 0);

    /// 清空并重建大小为 n 的并查集
    void reset(std::size_t n);

    /// 查找 x 所在集合的代表（根）
    int_t find(int_t x);

    /// 返回 x 所在集合的元素个数
    int_t set_size(int_t x);

    /// 判断 a 与 b 是否位于同一集合
    bool same(int_t a, int_t b);

    /// 当前独立集合个数
    int_t components() const { return component_cnt_; }

    /**
     * @brief 合并 a、b 所在集合
     * @return true  如果二者原本不连通并成功合并
     *         false 如果二者已在同一集合
     */
    bool unite(int_t a, int_t b);

private:
    std::vector<int_t> parent_;   // 负数：根 & -size；非负：父索引
    int_t component_cnt_{ 0 };      // 当前集合数
};
static const std::array<std::array<int, 5>, 5> HOXD70 = { {
        /*       A      C      G      T     N  */
        /*A*/ {  91,  -114,   -31,  -123, -100 },
        /*C*/ { -114,  100,  -125,   -31, -100 },
        /*G*/ {  -31, -125,   100,  -114, -100 },
        /*T*/ { -123,  -31,  -114,    91, -100 },
        /*N*/ { -100, -100,  -100,  -100, -100 }
    } };

// ------------------------------------------------------------------
// 创建字符到 HOXD70 索引的映射
// 支持大小写 ACGT，其他字符默认映射为 N（索引为 4）
// ------------------------------------------------------------------
static std::array<int8_t, 256> makeCharToIndex() {
    std::array<int8_t, 256> m;
    m.fill(4);  // 默认全部映射为 N
    m['A'] = 0; m['C'] = 1; m['G'] = 2; m['T'] = 3;
    m['a'] = 0; m['c'] = 1; m['g'] = 2; m['t'] = 3;
    return m;
}

// 全局常量，用于快速字符查表打分
static const auto ScoreChar2Idx = makeCharToIndex();

// ------------------------------------------------------------------
// 函数：subsScore
// 功能：返回两个碱基字符 a 和 b 在 HOXD70 中的替换得分
// 注意：大小写不敏感，非 ACGT 一律视为 N（惩罚）
// ------------------------------------------------------------------
inline int subsScore(char a, char b) {
    return HOXD70[ScoreChar2Idx[a]][ScoreChar2Idx[b]];
}
// ------------------------------------------------------------------
// 常量：碱基互补表（大写、小写都支持）
// 使用 C++17 constexpr lambda 构建 256 字节映射表
// ------------------------------------------------------------------
inline constexpr std::array<char, 256> BASE_COMPLEMENT = [] {
    std::array<char, 256> m{};
    for (auto& c : m) c = 'N'; // 默认所有字符都映射为 'N'
    m['A'] = 'T';  m['T'] = 'A';
    m['C'] = 'G';  m['G'] = 'C';
    m['a'] = 't';  m['t'] = 'a';
    m['c'] = 'g';  m['g'] = 'c';
    return m;
    }();

// 互补
inline void baseComplement(std::string& seq) {
    for (char& c : seq) {
        c = BASE_COMPLEMENT[static_cast<unsigned char>(c)];
    }
}

// 反转
inline void reverseSeq(std::string& seq) {
    std::reverse(seq.begin(), seq.end());
}

// 先取反，再反转
inline void reverseComplement(std::string& seq) {
    baseComplement(seq);
    reverseSeq(seq);
}


// ------------------------------------------------------------------
// CIGAR 表示与转换
// ------------------------------------------------------------------

// 单个 CIGAR 操作的压缩编码（uint32_t）
// 高位表示操作长度，低 4 位为操作类型编码
struct CigarSummary {
    uint64_t reference_length = 0;
    uint64_t query_length = 0;
    uint64_t alignment_length = 0;
    uint64_t match_length = 0;
};

struct AlignmentResult {
    Cigar_t cigar;
    CigarSummary summary;
};

struct KswSequenceView {
    std::string_view bases;
    bool reverse_complement = false;
};

CigarSummary summarizeCigar(const Cigar_t& cigar);





struct KSW2AlignConfig {
    const int8_t* mat;                  // 一维替换矩阵 (flattened 5x5)
    int alphabet_size;                 // 通常为 5
    int gap_open;                      // gap open penalty (positive)
    int gap_extend;                    // gap extend penalty (positive)
    int end_bonus;                     // 末端奖励分
    int zdrop = 100;                   // Z-drop 剪枝参数
    int band_width = -1;               // -1 表示全矩阵
    int flag = 0;      // 默认使用全替换矩阵
};


static const std::array<int8_t,25> scoring_matrix=[] {std::array<int8_t,25> m{};for(int i=0;i<5;++i)for(int j=0;j<5;++j){const int x=HOXD70[i][j];m[i*5+j]=static_cast<int8_t>((x+(x>=0?5:-5))/10);}return m;}();
static const int8_t* const dna5_simd_mat=scoring_matrix.data();
inline void init_simd_mat(){}
inline int ksw2GapOpenPenalty(){return 40;}
inline int ksw2GapExtendPenalty(){return 3;}
inline int auto_band(int qlen,int tlen,double rate=0.1,int margin=200){return margin+static_cast<int>(rate*(qlen+tlen/2));}
CigarSummary summarizeCigar(const Cigar_t&);
AlignmentResult globalAlignKSW2Result(KswSequenceView,KswSequenceView);
AlignmentResult extendAlignKSW2Result(KswSequenceView,KswSequenceView,int);

inline void checkKswLengths(KswSequenceView r,KswSequenceView q){if(r.bases.size()>static_cast<size_t>(INT32_MAX/16)||q.bases.size()>static_cast<size_t>(INT32_MAX/16))throw std::overflow_error("KSW2 length/score representation limit");}

uint64_t cigarToInt(char operation,uint64_t length);
void intToCigar(CigarUnit cigar,char& operation,uint64_t& length);
void appendCigar(Cigar_t& destination,const Cigar_t& source);
void prependCigar(Cigar_t& destination,const Cigar_t& source);
void appendCigarOp(Cigar_t& destination,char operation,uint64_t length);
void filterAnchorsByDpTreap(std::span<size_t> indices,AnchorVec& anchors,bool reference_side);
}
namespace ramag::genome_detail {
using KernelClock=std::chrono::steady_clock;
inline double Elapsed(KernelClock::time_point start) {
    return std::chrono::duration<double>(KernelClock::now()-start).count();
}
inline void KernelCheck(const PairwiseCoreOptions& o,std::string_view stage) {
    if(o.interruption_callback) o.interruption_callback(stage);
}
inline void KernelProgress(const PairwiseCoreOptions& o,std::string_view stage,std::uint64_t done,std::uint64_t total) {
    KernelCheck(o,stage);
    if(o.progress_callback) o.progress_callback(stage,done,total);
}
template<class Function> void ParallelKernel(std::size_t count,const PairwiseCoreOptions& o,Function f) {
    std::vector<std::exception_ptr> errors(count);
    if(count>static_cast<std::size_t>(INT64_MAX)) throw AlignmentError("pairwise task count overflow");
#ifdef _OPENMP
    const auto workers=std::max<std::size_t>(1,std::min<std::size_t>(count,o.threads));
#pragma omp parallel for schedule(dynamic) num_threads(workers) if(!omp_in_parallel() && count>1)
#endif
    for(std::int64_t i=0;i<static_cast<std::int64_t>(count);++i) {
        try { f(static_cast<std::size_t>(i)); } catch(...) {errors[static_cast<std::size_t>(i)]=std::current_exception();}
    }
    for(const auto& error:errors) if(error) std::rethrow_exception(error);
}
using Sequences=std::map<SequenceId,const SequenceRecord*>;
inline Sequences CheckSequences(std::span<const SequenceRecord> records) {
    Sequences lookup;
    for(const auto& r:records) {
        // Bounds protect every signed diagonal, linked-gap score and packed
        // CIGAR operation without silently narrowing public coordinates.
        if(r.bases.empty() || r.size()>static_cast<Length>(INT64_MAX/256)) throw AlignmentError("pairwise empty/oversized contig");
        if(!lookup.emplace(r.numeric_id,&r).second) throw AlignmentError("pairwise duplicate numeric sequence ID");
        if(r.bases.find_first_not_of("ACGTN")!=std::string::npos) throw AlignmentError("pairwise requires normalized A/C/G/T/N input");
    }
    return lookup;
}
inline char Complement(char base) {return pairwise_detail::BASE_COMPLEMENT[static_cast<unsigned char>(base)];}
struct NoAnchorVisit { void operator()(unsigned,Length,Length,char,char) const {} };
template<class Visitor=NoAnchorVisit>
AlignmentRecord ConvertAnchor(const pairwise_detail::Anchor& a,const Sequences& refs,const Sequences& queries,bool materialize=true,Visitor visit={}) {
    AlignmentRecord r;
    r.reference_id=a.reference_id; r.query_id=a.query_id;
    r.reference_begin=a.reference_begin; r.reference_end=a.reference_begin+a.ref_len;
    r.query_begin=a.query_begin; r.query_end=a.query_begin+a.qry_len;
    r.strand=a.strand==pairwise_detail::FORWARD?Strand::Forward:Strand::Reverse;
    const auto& ref=refs.at(r.reference_id)->bases; const auto& qry=queries.at(r.query_id)->bases;
    if(r.reference_end>ref.size() || r.query_end>qry.size() || a.ref_len==0 || a.qry_len==0) throw AlignmentError("pairwise alignment span outside sequence");
    Length rp=r.reference_begin, qp=0;
    const auto append=[&](char op,Length n){
        if(n==0)throw AlignmentError("pairwise zero CIGAR operation");
        if(!materialize)return;
        if(!r.cigar.empty()&&r.cigar.back().operation==op)r.cigar.back().length+=n;
        else {

            r.cigar.push_back({op,n});
        }
    };
    for(auto unit:a.cigar) {
        const Length n=unit>>4; const auto op=unit&15;
        if(n==0 || (op!=0&&op!=1&&op!=2&&op!=7&&op!=8)) throw AlignmentError("pairwise invalid packed CIGAR");
        if(op!=1 && n>r.reference_end-rp) throw AlignmentError("pairwise reference CIGAR overrun");
        if(op!=2 && n>a.qry_len-qp) throw AlignmentError("pairwise query CIGAR overrun");
        if(op==1 || op==2) {
            if constexpr(!std::is_same_v<Visitor,NoAnchorVisit>) {
                for(Length j=0;j<n;++j) {
                    const auto rpos=rp+(op==2?j:0),qpos=qp+(op==1?j:0);
                    const char x=op==1?'A':ref[rpos];
                    const char y=op==2?'A':(r.strand==Strand::Forward?qry[r.query_begin+qpos]:Complement(qry[r.query_end-qpos-1]));
                    visit(static_cast<unsigned>(op),rpos,qpos,x,y);
                }
            }
            append(op==1?'I':'D',n); r.edit_distance+=n;
            r.score-=40+3*static_cast<std::int64_t>(n);
            if(op==1)qp+=n;else rp+=n;
        } else {
            for(Length j=0;j<n;++j) {
                const char x=ref[rp++];
                const char y=r.strand==Strand::Forward?qry[r.query_begin+qp]:Complement(qry[r.query_end-qp-1]); ++qp;
                const bool same=x==y&&x!='N';
                if((op==7&&!same)||(op==8&&same))throw AlignmentError("pairwise EQX operation contradicts sequence");
                visit(static_cast<unsigned>(op),rp-1,qp-1,x,y);
                append(same?'=':'X',1); if(!same)++r.edit_distance;
                const auto xi=pairwise_detail::ScoreChar2Idx[static_cast<unsigned char>(x)];
                const auto yi=pairwise_detail::ScoreChar2Idx[static_cast<unsigned char>(y)];
                r.score+=pairwise_detail::dna5_simd_mat[xi*5+yi];
            }
        }
    }
    if(rp!=r.reference_end || qp!=a.qry_len) throw AlignmentError("pairwise incomplete CIGAR consumption");
    return r;
}

}

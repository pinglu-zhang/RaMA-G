// Pairwise alignment kernel. Source attribution and MIT notice are retained
// in third_party/attribution/. Do not substitute unrelated gap backends.
#include "ramag/pairwise_core.hpp"
#include "ksw2.h"
#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <tuple>
#include <unordered_map>
#include <variant>
#ifdef _OPENMP
#include <omp.h>
#endif
namespace ramag::pairwise_detail {
using int_t=int64_t; using uint_t=uint64_t;
using Score_t=int64_t;
using Coord_t=uint64_t; using Length_t=uint64_t; using ChrIndex=uint32_t;
using CigarUnit=uint64_t; using Cigar_t=std::vector<CigarUnit>;
enum Strand { FORWARD, REVERSE };
struct Match {
 ChrIndex ref_chr_index; Coord_t ref_start; ChrIndex qry_chr_index; Coord_t qry_start;
 Length_t length; Strand orientation;
 Length_t match_len() const {return length;}
 Strand strand() const {return orientation;}
 void set_match_len(Length_t n) {length=n;}
};
using MatchVec=std::vector<Match>; using MatchCluster=MatchVec;
using MatchClusterVec=std::vector<MatchCluster>;
using MatchClusterVecPtr=std::shared_ptr<MatchClusterVec>;
inline Coord_t start1(const Match& m){return m.ref_start;}
inline Coord_t start2(const Match& m){return m.qry_start;}
inline Length_t len1(const Match& m){return m.match_len();}
inline Length_t len2(const Match& m){return m.match_len();}
inline int_t diag(const Match& m){return static_cast<int_t>(start2(m))-static_cast<int_t>(start1(m));}
inline int_t diag_reverse(const Match& m){return -static_cast<int_t>(start2(m))-static_cast<int_t>(start1(m));}
inline void releaseCluster(MatchVec& v){MatchVec().swap(v);}
namespace SeqPro {
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
    ChrIndex ref_chr_index;
    Coord_t  ref_start;
    Length_t ref_len;
    ChrIndex qry_chr_index;
    Coord_t  qry_start;
    Length_t qry_len;

    Strand strand;
    uint_t alignment_length{}; // 对齐长度
    uint_t aligned_base{};
    Cigar_t cigar{};                 // 对齐的 CIGAR 字符串

    bool ref_selected = false;
    bool qry_selected = false;
    bool is_linked = false;

    Anchor() = default;

    Anchor(ChrIndex ref_chr, Coord_t ref_start, Length_t ref_len,
        ChrIndex qry_chr, Coord_t qry_start, Length_t qry_len,
        Strand strand, uint_t align_len, uint_t aligned_base, Cigar_t cigar_str)
        : ref_chr_index(ref_chr),
        ref_start(ref_start),
        ref_len(ref_len),
        qry_chr_index(qry_chr),
        qry_start(qry_start),
        qry_len(qry_len),
        strand(strand),
        alignment_length(align_len),
        aligned_base(aligned_base),
        cigar(std::move(cigar_str)) {
    }


};

using AnchorPtr = std::shared_ptr<Anchor>;
using AnchorPtrVec = std::vector<AnchorPtr>;
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
uint64_t cigarToInt(char operation, uint64_t len) {
	uint64_t opCode;

	// 将字符型操作转换为对应的数字编码（4 位）
	switch (operation) {
	case 'M': opCode = 0x0; break; // 对齐匹配（可能包含错配）
	case 'I': opCode = 0x1; break; // 插入
	case 'D': opCode = 0x2; break; // 删除
	case '=': opCode = 0x7; break; // 精确匹配（无错配）
	case 'X': opCode = 0x8; break; // 错配
		// TODO: 可以根据需要添加更多操作类型（如 soft clip, hard clip 等）
	default: opCode = 0xF; break; // 未知操作，使用保留码 0xF
	}

	// 将长度左移 4 位，然后用 bitwise OR 合并 opCode
	return (len << 4) | opCode;
}

void intToCigar(CigarUnit cigar, char& operation, uint64_t& len)
{
	uint64_t opCode = cigar & 0xF;  // 提取低 4 位为操作编码
	len = cigar >> 4;               // 高 28 位表示操作长度

	// 将操作编码转换回对应字符
	switch (opCode) {
	case 0x0: operation = 'M'; break;
	case 0x1: operation = 'I'; break;
	case 0x2: operation = 'D'; break;
	case 0x7: operation = '='; break;
	case 0x8: operation = 'X'; break;
	default:  operation = '?'; break; // 未知操作标记为 '?'
	}
}

void appendCigar(Cigar_t& dst, const Cigar_t& src)
{
	if (src.empty()) return;                  // nothing to do
	size_t idx = 0;

	//// 1) 检查“拼接点”——dst 的最后一个元素 vs. src 的首元素
	if (!dst.empty()) {
		char op_dst; uint64_t len_dst;
		intToCigar(dst.back(), op_dst, len_dst);

		char op_src; uint64_t len_src;
		intToCigar(src.front(), op_src, len_src);

		// 若操作码相同，则二者合并成一个
		if (op_dst == op_src) {
			dst.back() = cigarToInt(op_dst, len_dst + len_src);
			idx = 1;                         // src 第 0 元素已被合并，后续从 1 开始
		}
	}
	// 2) 将 src 剩余元素逐个 push
	for (; idx < src.size(); ++idx)
		dst.push_back(src[idx]);
}

void prependCigar(Cigar_t& dst, const Cigar_t& src)
{
    if (src.empty()) return; // nothing to do
    size_t end = src.size();

    if (!dst.empty()) {
        char op_dst; uint64_t len_dst;
        intToCigar(dst.front(), op_dst, len_dst);

        char op_src; uint64_t len_src;
        intToCigar(src.back(), op_src, len_src);

        if (op_dst == op_src) {
            // 合并 src.back() 和 dst.front()
            dst.front() = cigarToInt(op_dst, len_dst + len_src);
            end = src.size() - 1; // 最后一个元素已被合并，不再插入
        }
    }

    // 把 src[0 ... end-1] 插到 dst 前面
    dst.insert(dst.begin(), src.begin(), src.begin() + end);
}

void appendCigarOp(Cigar_t& dst, char op, uint64_t len)
{
	if (len == 0) return;
	if (!dst.empty()) {
		char op_last; uint64_t len_last;
		intToCigar(dst.back(), op_last, len_last);
		if (op_last == op) {                        // 合并
			dst.back() = cigarToInt(op, len_last + len);
			return;
		}
	}
	dst.push_back(cigarToInt(op, len));
}
UnionFind::UnionFind(std::size_t n) {
    reset(n);
}

// 重置并查集为 n 个独立集合
void UnionFind::reset(std::size_t n) {
    parent_.assign(n, -1);                      // -1 表示单元素集合（根节点，大小=1）
    component_cnt_ = static_cast<int_t>(n);     // 初始有 n 个连通分量
}

// ────────────────────────────────────────────
// Query（查询接口）
// ────────────────────────────────────────────

// 查找：返回元素 x 所在集合的根
int_t UnionFind::find(int_t x) {
    // 路径折半（迭代版）：不断把 x 挂到祖父节点上，加速后续查询
    while (parent_[x] >= 0 && parent_[parent_[x]] >= 0) {
        parent_[x] = parent_[parent_[x]];
        x = parent_[x];
    }
    return (parent_[x] < 0) ? x : parent_[x];
}

// 返回 x 所在集合大小
int_t UnionFind::set_size(int_t x) {
    return -parent_[find(x)];
}

// 判断 a 和 b 是否属于同一集合
bool UnionFind::same(int_t a, int_t b) {
    return find(a) == find(b);
}

// ────────────────────────────────────────────
// Modification（修改接口）
// ────────────────────────────────────────────

// 合并两个集合（按大小合并）
// 返回 true 表示确实发生了合并；false 表示原本就在同一集合
bool UnionFind::unite(int_t a, int_t b) {
    a = find(a);
    b = find(b);
    if (a == b) return false;

    // 按大小合并：parent_ 为负表示大小，数值越小(更负)集合越大
    // 若 a 的集合更小，则交换，让 a 做大根
    if (parent_[a] > parent_[b]) std::swap(a, b);

    parent_[a] += parent_[b];   // 更新根 a 的大小（负值累加）
    parent_[b] = a;             // b 挂到 a
    --component_cnt_;
    return true;
}

// ────────────────────────────────────────────
// filterAndMergeMatches：对 match 列表进行过滤与合并
// 假设：matches 已按 qry_start 排序
// 规则：
// 1) 同一 diagonal 且同向同染色体：合并成更长 match，删除后者
// 2) 同 ref 起点：根据重叠比例过滤较短者；同长时用 tentative 机制处理
// 3) 同 qry 起点：同上
// ────────────────────────────────────────────
void filterAndMergeMatches(MatchVec& matches) {
    if (matches.empty()) return;

    const size_t N = matches.size();
    std::vector<bool> good(N, true);        // good[i]=true 表示保留
    std::vector<bool> tentative(N, false);  // tentative 标记：用于“长度相等且重叠较大”的歧义处理

    for (size_t i = 0; i < N; ++i) {
        if (!good[i]) continue;

        const Match& mi = matches[i];
        int_t  i_diag = mi.qry_start - mi.ref_start;
        Coord_t i_end = mi.qry_start + mi.match_len();  // i 在 query 维度的结束位置（开区间右端）

        // 由于 matches 按 qry_start 排序，只需检查 qry_start <= i_end 的后续元素
        for (size_t j = i + 1; j < N && matches[j].qry_start <= i_end; ++j) {
            if (!good[j]) continue;

            const Match& mj = matches[j];
            int_t j_diag = static_cast<int_t>(mj.qry_start) - static_cast<int_t>(mj.ref_start);

            // --- Case 1: 同一 diagonal（同一条对角线） ---
            // 条件：diag 相等 + strand 相等 + ref/qry 染色体一致
            if (i_diag == j_diag &&
                mi.strand() == mj.strand() &&
                mi.ref_chr_index == mj.ref_chr_index &&
                mi.qry_chr_index == mj.qry_chr_index) {

                // 合并为更长的：计算 mj 相对 mi 的延伸长度
                Coord_t j_extent = mj.match_len() + mj.qry_start - mi.qry_start;
                if (j_extent > matches[i].match_len()) {
                    matches[i].set_match_len(j_extent);
                    i_end = mi.qry_start + j_extent;
                }
                good[j] = false; // 删除 mj
            }

            // --- Case 2: 同一 ref 起点 ---
            else if (mi.ref_start == mj.ref_start &&
                     mi.ref_chr_index == mj.ref_chr_index) {

                int_t overlap = mi.qry_start + mi.match_len() - mj.qry_start;

                if (mi.match_len() < mj.match_len()) {
                    // i 更短：若重叠超过 i 一半，丢弃 i
                    if (overlap >= static_cast<int_t>(mi.match_len() / 2)) {
                        good[i] = false;
                        break;
                    }
                }
                else if (mj.match_len() < mi.match_len()) {
                    // j 更短：若重叠超过 j 一半，丢弃 j
                    if (overlap >= static_cast<int_t>(mj.match_len() / 2))
                        good[j] = false;
                }
                else {
                    // 长度相等：若重叠超过一半，标记 tentative
                    if (overlap >= static_cast<int_t>(mi.match_len() / 2)) {
                        tentative[j] = true;
                        if (tentative[i]) {
                            // 如果 i 也已被 tentative，则丢弃 i
                            good[i] = false;
                            break;
                        }
                    }
                }
            }

            // --- Case 3: 同一 qry 起点 ---
            else if (mi.qry_start == mj.qry_start &&
                     mi.qry_chr_index == mj.qry_chr_index) {

                int64_t overlap = static_cast<int64_t>(mi.ref_start) +
                                  static_cast<int64_t>(mi.match_len()) -
                                  static_cast<int64_t>(mj.ref_start);

                if (mi.match_len() < mj.match_len()) {
                    if (overlap >= static_cast<int64_t>(mi.match_len() / 2)) {
                        good[i] = false;
                        break;
                    }
                }
                else if (mj.match_len() < mi.match_len()) {
                    if (overlap >= static_cast<int64_t>(mj.match_len() / 2))
                        good[j] = false;
                }
                else {
                    if (overlap >= static_cast<int64_t>(mi.match_len() / 2)) {
                        tentative[j] = true;
                        if (tentative[i]) {
                            good[i] = false;
                            break;
                        }
                    }
                }
            }
        }
    }

    // 收集 good 的 matches，生成新数组
    MatchVec filtered;
    filtered.reserve(matches.size());
    for (size_t i = 0; i < N; ++i) {
        if (good[i]) filtered.push_back(matches[i]);
    }
    matches.swap(filtered);
}

// ────────────────────────────────────────────
// buildClusters：将 unique_match 根据 max_gap / diagdiff / diagfactor 聚类
// 逻辑：
// 1) 对 match 排序（按 start2，再按 start1）
// 2) filterAndMergeMatches 做合并压缩
// 3) 使用并查集：若两个 match 在 query 维度距离 sep <= max_gap
//    且 diagonal 差 <= max(diagdiff, diagfactor * sep)，则归为一类
// 4) 最后根据并查集根构建簇列表
// ────────────────────────────────────────────
MatchClusterVec buildClusters(MatchVec& unique_match,
                             int_t  max_gap,
                             int_t  diagdiff,
                             double diagfactor) {
    MatchClusterVec clusters;

    // 特殊情况：0 或 1 个元素，直接返回
    if (unique_match.size() < 2) {
        if (unique_match.size() == 1) {
            clusters.emplace_back();
            clusters.back().push_back(std::move(unique_match[0]));
        }
        return clusters;
    }

    // 判断链方向（假设同一批次一致）
    const bool is_forward = (unique_match.front().strand() == FORWARD);

    // 先排序：按 start2 再按 start1
    std::sort(unique_match.begin(), unique_match.end(),
        [](const Match& a, const Match& b) {
            if (start2(a) < start2(b)) return true;
            if (start2(a) > start2(b)) return false;
            return start1(a) < start1(b);
        });

    // 再合并压缩（可能会删除元素、缩短 vector）
    filterAndMergeMatches(unique_match);

    // 重新获取 N
    const uint_t N = static_cast<uint_t>(unique_match.size());
    if (N < 2) {
        if (N == 1) {
            clusters.emplace_back();
            clusters.back().push_back(std::move(unique_match[0]));
        }
        return clusters;
    }

    UnionFind uf(N);

    // 利用排序后的局部性进行聚类：对每个 i，只检查后续 sep 不超过 max_gap 的 j
    for (uint_t i = 0; i < N; ++i) {
        uint_t i_end = start2(unique_match[i]) + len2(unique_match[i]);

        int_t i_diag = 0;
        if (is_forward) {
            i_diag = diag(unique_match[i]);
        }
        else {
            i_diag = diag_reverse(unique_match[i]);
        }

        for (uint_t j = i + 1; j < N; ++j) {
            int_t sep = static_cast<int_t>(start2(unique_match[j])) - static_cast<int_t>(i_end);

            // 早停：gap 太大则后面的 j 也不可能满足
            if (sep > static_cast<int_t>(max_gap)) break;

            int_t diag_diff = 0;
            if (is_forward) {
                diag_diff = std::abs(diag(unique_match[j]) - i_diag);
            }
            else {
                diag_diff = std::abs(diag_reverse(unique_match[j]) - i_diag);
            }

            // 阈值：max(diagdiff, diagfactor * sep)
            int_t th = std::max(diagdiff, static_cast<int_t>(diagfactor * sep));

            // diagonal 差满足阈值：归并到同一簇
            if (diag_diff <= th) {
                uf.unite(i, j);
            }
        }
    }

    // 根据并查集根构建簇：root -> cluster_id
    std::unordered_map<int_t, int_t> root_to_cluster_id;
    root_to_cluster_id.reserve(N / 4);

    for (uint_t idx = 0; idx < unique_match.size(); idx++) {
        int_t root = uf.find(idx);
        auto it = root_to_cluster_id.find(root);

        int_t cid;
        if (it == root_to_cluster_id.end()) {
            cid = static_cast<int_t>(clusters.size());
            clusters.emplace_back();
            clusters.back().reserve(1);
            root_to_cluster_id[root] = cid;
        }
        else {
            cid = it->second;
        }

        clusters[cid].push_back(std::move(unique_match[idx]));
    }

    return clusters;
}

// ────────────────────────────────────────────
MatchVec bestChainDP(MatchVec& cluster, double diagfactor) {
    if (cluster.empty()) return {};
    if (cluster.size() == 1) return MatchVec{cluster.front()};

    Strand strand = cluster.front().strand();

    std::sort(cluster.begin(), cluster.end(),
        [](const Match& a, const Match& b) { return start2(a) < start2(b); });

    const uint_t N = static_cast<uint_t>(cluster.size());
    std::vector<int_t> score(N), pred(N, -1);
    uint_t best_idx = 0;

    for (uint_t i = 0; i < N; ++i) {
        // 初始分：自身长度
        score[i] = len2(cluster[i]);

        for (uint_t j = 0; j < i; ++j) {
            // query 维度必须不重叠（i 的 start2 需在 j 的末尾之后）
            if (start2(cluster[i]) <= start2(cluster[j]) + len2(cluster[j])) continue;

            int_t d = 0;

            if (strand == FORWARD) {
                // ref 维度也不重叠：i 的 start1 需在 j 的 ref 末尾之后
                int_t prev_endj = start1(cluster[j]) + len1(cluster[j]);
                if (static_cast<int_t>(start1(cluster[i])) <= prev_endj) continue;

                d = std::abs(diag(cluster[i]) - diag(cluster[j]));
            }
            else {
                // 反向链：按原逻辑计算不交叉条件与 sep
                int_t prev_endi = start1(cluster[i]) + len1(cluster[i]);
                if (prev_endi >= static_cast<int_t>(start1(cluster[j]))) continue;

                d = std::abs(diag_reverse(cluster[i]) - diag_reverse(cluster[j]));
            }

            // 候选分数：前链分数 + 当前长度 - diagonal 差惩罚
            int_t cand = score[j] + len2(cluster[i]) - d;
            if (cand > score[i]) {
                score[i] = cand;
                pred[i] = static_cast<int_t>(j);
            }
        }

        if (score[i] > score[best_idx]) best_idx = i;
    }

    // 回溯构建最优链
    MatchVec chain;
    for (int_t k = static_cast<int_t>(best_idx); k != -1; k = pred[k])
        chain.emplace_back(cluster[k]);
    std::reverse(chain.begin(), chain.end());

    return chain;
}
MatchClusterVecPtr clusterChrMatch(MatchVec& unique_match,
                                  uint_t min_cluster_length,
                                  int_t  max_gap,
                                  int_t  diagdiff,
                                  double diagfactor) {

    auto best_chain_clusters = std::make_shared<MatchClusterVec>();

    if (unique_match.size() == 0) {
        return best_chain_clusters;
    }

    // 1) 聚簇
    MatchClusterVec clusters = buildClusters(unique_match, max_gap, diagdiff, diagfactor);

    // 释放 unique_match 的内存（保持原逻辑）
    unique_match.clear();
    unique_match.shrink_to_fit();

    best_chain_clusters->reserve(clusters.size());

    // 2) 每个簇选最佳链
    for (auto& cluster : clusters) {
        if (cluster.empty()) continue;

        MatchVec best_chain = bestChainDP(cluster, diagfactor);

        if (best_chain.empty()) {
            releaseCluster(cluster);
            continue;
        }

        uint_t span = 0;
        // 遍历 best_chain 累加 span（保持原逻辑）
        for (auto& m : best_chain) {
            span += m.match_len();
        }

        // 满足最小簇长度才保留
        if (span >= min_cluster_length) {
            best_chain_clusters->emplace_back(std::move(best_chain));
        }

        // 回收 cluster 剩余元素
        releaseCluster(cluster);
    }

    best_chain_clusters->shrink_to_fit();
    return best_chain_clusters;
}
namespace {
constexpr size_t kMaximumRetainedKswEncoding = 64 * 1024;

struct KswEncodingScratch {
    std::vector<uint8_t> reference;
    std::vector<uint8_t> query;
};

thread_local KswEncodingScratch ksw_encoding_scratch;

char orientedBase(const KswSequenceView sequence, size_t index) {
    if (!sequence.reverse_complement) {
        return sequence.bases[index];
    }
    const size_t reversed_index = sequence.bases.size() - index - 1;
    return BASE_COMPLEMENT[static_cast<uint8_t>(
        sequence.bases[reversed_index])];
}

void encodeSequence(const KswSequenceView sequence,
                    std::vector<uint8_t>& encoded) {
    encoded.resize(sequence.bases.size());
    for (size_t index = 0; index < sequence.bases.size(); ++index) {
        encoded[index] = ScoreChar2Idx[static_cast<uint8_t>(
            orientedBase(sequence, index))];
    }
}

bool isCanonicalExactMatch(const KswSequenceView reference,
                           const KswSequenceView query) {
    if (reference.bases.empty() ||
        reference.bases.size() != query.bases.size()) {
        return false;
    }
    for (size_t index = 0; index < reference.bases.size(); ++index) {
        const char ref_base = orientedBase(reference, index);
        const char query_base = orientedBase(query, index);
        if (ref_base != query_base ||
            (ref_base != 'A' && ref_base != 'C' &&
             ref_base != 'G' && ref_base != 'T')) {
            return false;
        }
    }
    return true;
}

std::vector<uint8_t>& selectEncodingBuffer(
    size_t size, std::vector<uint8_t>& retained,
    std::vector<uint8_t>& temporary) {
    return size <= kMaximumRetainedKswEncoding ? retained : temporary;
}

AlignmentResult copyKswResultAndSummarize(const ksw_extz_t& result) {
    AlignmentResult output;
    output.cigar.reserve(result.n_cigar);
    for (int index = 0; index < result.n_cigar; ++index) {
        const CigarUnit unit = result.cigar[index];
        output.cigar.push_back(unit);
        const uint64_t length = unit >> 4;
        const uint32_t operation = unit & 0xf;
        output.summary.alignment_length += length;
        if (operation != 1) output.summary.reference_length += length;
        if (operation != 2) output.summary.query_length += length;
        if (operation == 0) output.summary.match_length += length;
    }
    return output;
}

struct GlobalKswRun {
    Cigar_t cigar;
    CigarSummary summary;
    int score = std::numeric_limits<int>::min();
};

int maximumSubstitutionScore(const int8_t* matrix, int alphabet_size) {
    int maximum = std::numeric_limits<int>::min();
    for (int index = 0; index < alphabet_size * alphabet_size; ++index) {
        maximum = std::max(maximum, static_cast<int>(matrix[index]));
    }
    return maximum;
}

void appendSimpleGlobalCigar(
    Cigar_t& cigar,
    size_t prefix_length,
    size_t gap_length,
    char gap_operation,
    size_t suffix_length) {
    cigar.clear();
    cigar.reserve(3);
    if (prefix_length > 0) {
        appendCigarOp(
            cigar, 'M', static_cast<uint64_t>(prefix_length));
    }
    if (gap_length > 0) {
        appendCigarOp(
            cigar, gap_operation, static_cast<uint64_t>(gap_length));
    }
    if (suffix_length > 0) {
        appendCigarOp(
            cigar, 'M', static_cast<uint64_t>(suffix_length));
    }
}

// Accept a direct/no-extra-gap path only when it strictly beats an upper
// bound for every path with another gap run. Strict comparison plus a unique
// one-gap placement preserves KSW2's tie-breaking contract.
bool tryProvablyOptimalSimpleGlobalAlignment(
    KswSequenceView reference,
    KswSequenceView query,
    const int8_t* matrix,
    int maximum_substitution_score,
    int gap_open,
    int gap_extend,
    Cigar_t& cigar) {
    const size_t reference_length = reference.bases.size();
    const size_t query_length = query.bases.size();
    const size_t shorter_length =
        std::min(reference_length, query_length);
    const size_t length_difference =
        reference_length > query_length
            ? reference_length - query_length
            : query_length - reference_length;

    if (shorter_length == 0) {
        appendSimpleGlobalCigar(
            cigar, 0, length_difference,
            reference_length > query_length ? 'D' : 'I', 0);
        return true;
    }
    if (gap_open < 0 || gap_extend < 0 ||
        maximum_substitution_score < 0) {
        return false;
    }

    const auto score_pair = [matrix](char reference_base, char query_base) {
        return static_cast<int64_t>(
            matrix[
                ScoreChar2Idx[static_cast<uint8_t>(reference_base)] * 5 +
                ScoreChar2Idx[static_cast<uint8_t>(query_base)]]);
    };

    size_t best_gap_position = 0;
    int64_t best_score = 0;
    bool best_gap_position_is_unique = true;
    if (reference_length == query_length) {
        for (size_t index = 0; index < shorter_length; ++index) {
            best_score += score_pair(
                orientedBase(reference, index),
                orientedBase(query, index));
        }
    } else if (reference_length > query_length) {
        int64_t current_score = 0;
        for (size_t index = 0; index < shorter_length; ++index) {
            current_score += score_pair(
                orientedBase(reference, index + length_difference),
                orientedBase(query, index));
        }
        best_score = current_score;
        for (size_t gap_position = 0;
             gap_position < shorter_length;
             ++gap_position) {
            current_score -= score_pair(
                orientedBase(reference, gap_position + length_difference),
                orientedBase(query, gap_position));
            current_score +=
                score_pair(orientedBase(reference, gap_position),
                           orientedBase(query, gap_position));
            if (current_score > best_score) {
                best_score = current_score;
                best_gap_position = gap_position + 1;
                best_gap_position_is_unique = true;
            } else if (current_score == best_score) {
                best_gap_position_is_unique = false;
            }
        }
    } else {
        int64_t current_score = 0;
        for (size_t index = 0; index < shorter_length; ++index) {
            current_score += score_pair(
                orientedBase(reference, index),
                orientedBase(query, index + length_difference));
        }
        best_score = current_score;
        for (size_t gap_position = 0;
             gap_position < shorter_length;
             ++gap_position) {
            current_score -= score_pair(
                orientedBase(reference, gap_position),
                orientedBase(query, gap_position + length_difference));
            current_score +=
                score_pair(orientedBase(reference, gap_position),
                           orientedBase(query, gap_position));
            if (current_score > best_score) {
                best_score = current_score;
                best_gap_position = gap_position + 1;
                best_gap_position_is_unique = true;
            } else if (current_score == best_score) {
                best_gap_position_is_unique = false;
            }
        }
    }

    if (length_difference > 0) {
        best_score -=
            static_cast<int64_t>(gap_open) +
            static_cast<int64_t>(gap_extend) * length_difference;
    }

    if (length_difference > 0 && !best_gap_position_is_unique) {
        return false;
    }

    int64_t excluded_path_upper_bound;
    if (length_difference == 0) {
        excluded_path_upper_bound =
            static_cast<int64_t>(maximum_substitution_score) *
                (shorter_length - 1) -
            2LL * gap_open - 2LL * gap_extend;
    } else {
        excluded_path_upper_bound =
            static_cast<int64_t>(maximum_substitution_score) *
                shorter_length -
            static_cast<int64_t>(gap_extend) * length_difference -
            2LL * gap_open;
    }
    if (best_score <= excluded_path_upper_bound) {
        return false;
    }

    appendSimpleGlobalCigar(
        cigar,
        best_gap_position,
        length_difference,
        reference_length > query_length ? 'D' : 'I',
        shorter_length - best_gap_position);
    return true;
}

GlobalKswRun runGlobalKsw(
    const std::vector<uint8_t>& reference,
    const std::vector<uint8_t>& query,
    const KSW2AlignConfig& config,
    int band_width) {
    ksw_extz_t result{};
    ksw_extz2_sse(
        nullptr,
        static_cast<int>(query.size()), query.data(),
        static_cast<int>(reference.size()), reference.data(),
        config.alphabet_size, config.mat,
        config.gap_open, config.gap_extend,
        band_width, config.zdrop, config.end_bonus,
        config.flag, &result);

    GlobalKswRun run;
    run.score = result.score;
    run.cigar.reserve(result.n_cigar);
    for (int index = 0; index < result.n_cigar; ++index) {
        const CigarUnit unit = result.cigar[index];
        run.cigar.push_back(unit);
        const uint64_t length = unit >> 4;
        const uint32_t operation = unit & 0xf;
        run.summary.alignment_length += length;
        if (operation != 1) run.summary.reference_length += length;
        if (operation != 2) run.summary.query_length += length;
        if (operation == 0) run.summary.match_length += length;
    }
    free(result.cigar);
    return run;
}

// A path that leaves a width-w diagonal band and finishes at length
// difference d <= w needs at least 2(w + 1) - d gap bases and two gap opens.
// Give every remaining aligned pair the best matrix score; a strict win over
// that upper bound proves every global optimum stays inside the band.
bool bandExcludesEveryGlobalOptimum(
    size_t reference_length,
    size_t query_length,
    int band_width,
    int score,
    int maximum_substitution_score,
    int gap_open,
    int gap_extend) {
    if (band_width < 0) {
        return true;
    }
    const int64_t length_difference = std::abs(
        static_cast<int64_t>(query_length) -
        static_cast<int64_t>(reference_length));
    if (band_width < length_difference) {
        return false;
    }

    const int64_t minimum_gap_bases =
        2LL * (static_cast<int64_t>(band_width) + 1) -
        length_difference;
    const int64_t total_length =
        static_cast<int64_t>(reference_length) +
        static_cast<int64_t>(query_length);
    if (minimum_gap_bases > total_length) {
        return true;
    }

    const int64_t outside_score_upper_bound_times_two =
        static_cast<int64_t>(maximum_substitution_score) *
            (total_length - minimum_gap_bases) -
        2LL * gap_extend * minimum_gap_bases -
        4LL * gap_open;
    return 2LL * score > outside_score_upper_bound_times_two;
}

int bandRequiredToCertifyScore(
    size_t reference_length,
    size_t query_length,
    int score,
    int maximum_substitution_score,
    int gap_open,
    int gap_extend) {
    const int64_t total_length =
        static_cast<int64_t>(reference_length) +
        static_cast<int64_t>(query_length);
    const int64_t length_difference = std::abs(
        static_cast<int64_t>(query_length) -
        static_cast<int64_t>(reference_length));
    const int64_t coefficient =
        static_cast<int64_t>(maximum_substitution_score) +
        2LL * gap_extend;
    const int64_t numerator =
        static_cast<int64_t>(maximum_substitution_score) * total_length -
        4LL * gap_open - 2LL * score;
    const int64_t required_gap_bases =
        numerator < 0 ? 0 : numerator / coefficient + 1;
    const int64_t required_band =
        (required_gap_bases + length_difference + 1) / 2 - 1;
    return static_cast<int>(
        std::max(length_difference, required_band));
}

}  // namespace

Cigar_t globalAlignKSW2(
    const std::string& reference,
    const std::string& query) {
    return globalAlignKSW2Result(
        {reference, false}, {query, false}).cigar;
}

AlignmentResult globalAlignKSW2Result(
    KswSequenceView reference,
    KswSequenceView query) {
    checkKswLengths(reference,query);
    init_simd_mat();

    KSW2AlignConfig config;
    config.mat = dna5_simd_mat;
    config.alphabet_size = 5;
    config.gap_open = ksw2GapOpenPenalty();
    config.gap_extend = ksw2GapExtendPenalty();
    config.end_bonus = 0;
    config.zdrop = -1;
    config.band_width = -1;
    config.flag = KSW_EZ_RIGHT | KSW_EZ_GENERIC_SC;

    Cigar_t simple_cigar;
    const int maximum_substitution_score =
        maximumSubstitutionScore(
            config.mat, config.alphabet_size);
    if (tryProvablyOptimalSimpleGlobalAlignment(
            reference, query, config.mat,
            maximum_substitution_score,
            config.gap_open, config.gap_extend,
            simple_cigar)) {
        AlignmentResult result;
        result.cigar = std::move(simple_cigar);
        result.summary = summarizeCigar(result.cigar);
        return result;
    }

    std::vector<uint8_t> temporary_reference;
    std::vector<uint8_t> temporary_query;
    auto& encoded_reference = selectEncodingBuffer(
        reference.bases.size(), ksw_encoding_scratch.reference,
        temporary_reference);
    auto& encoded_query = selectEncodingBuffer(
        query.bases.size(), ksw_encoding_scratch.query,
        temporary_query);
    encodeSequence(reference, encoded_reference);
    encodeSequence(query, encoded_query);

    const size_t maximum_length =
        std::max(reference.bases.size(), query.bases.size());
    const int64_t length_difference = std::abs(
        static_cast<int64_t>(query.bases.size()) -
        static_cast<int64_t>(reference.bases.size()));
    const bool can_certify_band =
        maximum_substitution_score >= 0 &&
        config.gap_open >= 0 &&
        config.gap_extend >= 0 &&
        maximum_substitution_score + 2 * config.gap_extend > 0;
    if (!can_certify_band || maximum_length <= 64) {
        GlobalKswRun run = runGlobalKsw(
            encoded_reference, encoded_query, config, -1);
        return {std::move(run.cigar), run.summary};
    }

    int band_width = static_cast<int>(
        std::max<int64_t>(64, length_difference));
    if (static_cast<size_t>(band_width) >= maximum_length) {
        GlobalKswRun run = runGlobalKsw(
            encoded_reference, encoded_query, config, -1);
        return {std::move(run.cigar), run.summary};
    }

    GlobalKswRun banded = runGlobalKsw(
        encoded_reference, encoded_query, config, band_width);
    const auto consumes_complete_input = [&reference, &query](
                                             const GlobalKswRun& run) {
        return run.summary.reference_length == reference.bases.size() &&
               run.summary.query_length == query.bases.size();
    };
    if (consumes_complete_input(banded) &&
        bandExcludesEveryGlobalOptimum(
            reference.bases.size(), query.bases.size(), band_width,
            banded.score, maximum_substitution_score,
            config.gap_open, config.gap_extend)) {
        return {std::move(banded.cigar), banded.summary};
    }

    const int required_band = bandRequiredToCertifyScore(
        reference.bases.size(), query.bases.size(), banded.score,
        maximum_substitution_score,
        config.gap_open, config.gap_extend);
    if (required_band > band_width &&
        static_cast<size_t>(required_band) < maximum_length) {
        banded = runGlobalKsw(
            encoded_reference, encoded_query, config, required_band);
        if (consumes_complete_input(banded) &&
            bandExcludesEveryGlobalOptimum(
                reference.bases.size(), query.bases.size(), required_band,
                banded.score, maximum_substitution_score,
                config.gap_open, config.gap_extend)) {
            return {std::move(banded.cigar), banded.summary};
        }
    }

    GlobalKswRun run = runGlobalKsw(
        encoded_reference, encoded_query, config, -1);
    return {std::move(run.cigar), run.summary};
}

Cigar_t extendAlignKSW2(const std::string& ref,
    const std::string& query,
    int zdrop)
{
    return extendAlignKSW2Result(
        {ref, false}, {query, false}, zdrop).cigar;
}

CigarSummary summarizeCigar(const Cigar_t& cigar) {
    CigarSummary summary;
    for (const CigarUnit unit : cigar) {
        const uint64_t length = unit >> 4;
        const uint32_t operation = unit & 0xf;
        summary.alignment_length += length;
        if (operation != 1) summary.reference_length += length;
        if (operation != 2) summary.query_length += length;
        if (operation == 0) summary.match_length += length;
    }
    return summary;
}

static AlignmentResult extendAlignKSW2Impl(
    KswSequenceView ref, KswSequenceView query, int zdrop,
    bool allow_shortcuts)
{
    checkKswLengths(ref,query);
    if (allow_shortcuts && ref.bases.empty() && query.bases.empty()) {
        return {};
    }
    if (allow_shortcuts && isCanonicalExactMatch(ref, query)) {
        AlignmentResult result;
        const auto length = static_cast<uint64_t>(ref.bases.size());
        result.cigar.push_back(cigarToInt('M', length));
        result.summary = {length, length, length, length};
        return result;
    }

    /* ---------- 1. 序列编码 ---------- */
    std::vector<uint8_t> temporary_reference;
    std::vector<uint8_t> temporary_query;
    auto& ref_enc = selectEncodingBuffer(
        ref.bases.size(), ksw_encoding_scratch.reference,
        temporary_reference);
    auto& qry_enc = selectEncodingBuffer(
        query.bases.size(), ksw_encoding_scratch.query,
        temporary_query);
    encodeSequence(ref, ref_enc);
    encodeSequence(query, qry_enc);

    ///* ---------- 2. 配置 ---------- */
    //KSW2AlignConfig cfg = makeTurboKSW2Config(query.size(), ref.size());
    ////KSW2AlignConfig cfg;
    //cfg.zdrop = zdrop;       // 用于提前终止
    //cfg.flag = KSW_EZ_EXTZ_ONLY     // ends-free extension
    //    | KSW_EZ_APPROX_MAX    // 跟踪 ez.max_q/max_t
    //    | KSW_EZ_APPROX_DROP   // 在 approximate 模式下触发 z-drop 就中断
    //    | KSW_EZ_RIGHT;        // （可选）gap 右对齐     // **关键**：启用 extension/ends-free
    //// 若需要右对齐 gaps 建议保留 KSW_EZ_RIGHT
    //cfg.end_bonus = 100;
    //cfg.band_width = -1;
    init_simd_mat();
    KSW2AlignConfig cfg;
	cfg.mat = dna5_simd_mat;
    cfg.zdrop = zdrop;
    cfg.flag = KSW_EZ_EXTZ_ONLY | KSW_EZ_RIGHT |
        KSW_EZ_APPROX_DROP | KSW_EZ_GENERIC_SC;
    cfg.end_bonus = 50;
    cfg.alphabet_size = 5;
    cfg.gap_open = ksw2GapOpenPenalty();
    cfg.gap_extend = ksw2GapExtendPenalty();
    cfg.band_width = auto_band(ref.bases.size(), query.bases.size());


    /* ---------- 3. 调用 KSW2 ---------- */
    ksw_extz_t ez{};
    ksw_extz2_sse(nullptr,
        static_cast<int>(qry_enc.size()), qry_enc.data(),
        static_cast<int>(ref_enc.size()), ref_enc.data(),
        cfg.alphabet_size, cfg.mat,
        cfg.gap_open, cfg.gap_extend,
        cfg.band_width, cfg.zdrop, cfg.end_bonus,
        cfg.flag, &ez);

    // 赋值bool& if_zdrop,int& ref_end,int& qry_end
    /* ---------- 4. 拷贝 & 释放 ---------- */
    AlignmentResult result = copyKswResultAndSummarize(ez);
    free(ez.cigar);                    // ksw2 使用 malloc
    return result;
}

AlignmentResult extendAlignKSW2Result(
    KswSequenceView ref, KswSequenceView query, int zdrop) {
    return extendAlignKSW2Impl(ref, query, zdrop, true);
}

namespace {
constexpr Coord_t kMaximumLinkedGapAlignmentLength = 10000;

constexpr bool linkedGapCanBeAligned(Coord_t ref_gap, Coord_t query_gap) {
    return ref_gap <= kMaximumLinkedGapAlignmentLength &&
           query_gap <= kMaximumLinkedGapAlignmentLength;
}

static_assert(linkedGapCanBeAligned(9999, 9999));
static_assert(linkedGapCanBeAligned(10000, 10000));
static_assert(!linkedGapCanBeAligned(10001, 10000));
static_assert(!linkedGapCanBeAligned(10000, 10001));

struct ManagedSequenceSlice {
    std::string storage;
    KswSequenceView view;
};

const SeqPro::SequenceManager& originalSequenceManager(
    const SeqPro::ManagerVariant& manager) {
    return std::visit([](const auto& pointer) -> const SeqPro::SequenceManager& {
        using Pointer = std::decay_t<decltype(pointer)>;
        if constexpr (std::is_same_v<
                          Pointer,
                          std::unique_ptr<SeqPro::SequenceManager>>) {
            return *pointer;
        } else {
            return pointer->getOriginalManager();
        }
    }, manager);
}

void loadSequenceSlice(const SeqPro::ManagerVariant& manager,
                       ChrIndex chromosome, Coord_t start, Coord_t length,
                       bool reverse_complement,
                       ManagedSequenceSlice& output) {
    output.storage.clear();
    if (length == 0) {
        output.view = {{}, reverse_complement};
        return;
    }
    const auto& original = originalSequenceManager(manager);
    std::string_view contiguous;
    if (original.tryGetContiguousSubSequence(
            chromosome, start, length, contiguous)) {
        output.view = {contiguous, reverse_complement};
        return;
    }
    original.getSubSequenceInto(
        chromosome, start, length, output.storage);
    output.view = {output.storage, reverse_complement};
}

void appendCigarMove(Cigar_t& destination, Cigar_t& source) {
    if (source.empty()) return;
    if (destination.empty()) {
        destination = std::move(source);
        Cigar_t().swap(source);
        return;
    }
    destination.reserve(destination.size() + source.size());
    appendCigar(destination, source);
    Cigar_t().swap(source);
}
}  // namespace

Anchor extendClusterToAnchor(MatchCluster& cluster,
    const SeqPro::ManagerVariant& ref_mgr,
    const SeqPro::ManagerVariant& query_mgr) {
    if (cluster.empty()) return Anchor();
    Anchor anchor;
    const Match& first = cluster.front();

    Strand strand = first.strand();
    bool   fwd = (strand == FORWARD);

    if (!fwd) {
        std::sort(cluster.begin(), cluster.end(),
            [](auto& a, auto& b) { return a.ref_start < b.ref_start; });
    }

    ChrIndex ref_chr = first.ref_chr_index;
    ChrIndex qry_chr = first.qry_chr_index;

    Cigar_t cigar; cigar.reserve(cluster.size() * 2);  // 预估
    Coord_t aln_len = 0;
    Coord_t match_len = 0;
    ManagedSequenceSlice reference_slice;
    ManagedSequenceSlice query_slice;
    for (size_t i = 0;i < cluster.size();++i) {
        const Match& m = cluster[i];

		uint_t len = len1(m);
        uint_t ref_start = m.ref_start + len;
        uint_t qry_start = m.qry_start + len;

        appendCigarOp(cigar, 'M', len);
        aln_len += len;
        match_len += len;

        if (i + 1 == cluster.size()) break;

        const Match& nxt = cluster[i + 1];
        uint_t len2 = len1(nxt);
		uint_t ref_end = nxt.ref_start;
		uint_t qry_end = nxt.qry_start;
        Coord_t query_gap_begin = 0;
        Coord_t query_gap_length = 0;
        if (fwd) {
            query_gap_begin = qry_start;
            query_gap_length = qry_end - qry_start;
        }
        else {
            query_gap_begin = qry_end + len2;
            query_gap_length = m.qry_start - len2 - qry_end;
        }

        loadSequenceSlice(ref_mgr, ref_chr, ref_start,
            ref_end - ref_start, false, reference_slice);
        loadSequenceSlice(query_mgr, qry_chr, query_gap_begin,
            query_gap_length, !fwd, query_slice);
        AlignmentResult gap = globalAlignKSW2Result(
            reference_slice.view, query_slice.view);
        match_len += gap.summary.match_length;
        aln_len += gap.summary.alignment_length;
        appendCigarMove(cigar, gap.cigar);

    }
	const Match& last = cluster.back();
    if (fwd) {
        anchor = Anchor(ref_chr, first.ref_start, last.ref_start + last.match_len() - first.ref_start, qry_chr, first.qry_start, last.qry_start + last.match_len() - first.qry_start, strand, aln_len, match_len, std::move(cigar));
    }
    else {
        anchor = Anchor(ref_chr, first.ref_start, last.ref_start + last.match_len() - first.ref_start, qry_chr, last.qry_start, first.qry_start + first.match_len() - last.qry_start, strand, aln_len, match_len, std::move(cigar));
    }

    return anchor;
}
inline bool intervalOverlap(Coord_t a_lo, Coord_t a_hi, Coord_t b_lo, Coord_t b_hi) {
    // 闭开区间 [lo,hi)，不重叠当且仅当 a_hi <= b_lo 或 b_hi <= a_lo
    return !(a_hi <= b_lo || b_hi <= a_lo);
}
namespace AnchorLinkDetail {

struct Statistics {
    uint64_t candidate_checks{0};
    uint64_t sequence_extractions{0};
    uint64_t direct_ksw_calls{0};
    uint64_t fallback_ksw_calls{0};
    uint64_t long_gap_rejections{0};
    uint64_t maximum_seen_gap{0};
    uint64_t estimated_ksw_cells{0};
};

struct ComponentRange {
    size_t begin{0};
    size_t end{0};
    uint64_t estimated_cost{0};
};

AnchorVec materializeClusterAnchors(
    MatchClusterVec& clusters,
    const SeqPro::ManagerVariant& ref_mgr,
    const SeqPro::ManagerVariant& qry_mgr);

std::vector<ComponentRange> splitAnchorComponents(
    const AnchorVec& anchors,
    Statistics* statistics = nullptr);

AnchorPtrVec linkAnchorRange(
    AnchorVec& anchors,
    size_t begin,
    size_t end,
    const SeqPro::ManagerVariant& ref_mgr,
    const SeqPro::ManagerVariant& qry_mgr,
    Statistics* statistics = nullptr);

}  // namespace AnchorLinkDetail

namespace AnchorLinkDetail {

namespace {

uint64_t saturatingAdd(uint64_t left, uint64_t right) {
    return right > std::numeric_limits<uint64_t>::max() - left
        ? std::numeric_limits<uint64_t>::max()
        : left + right;
}

uint64_t saturatingMultiply(uint64_t left, uint64_t right) {
    if (left != 0 &&
        right > std::numeric_limits<uint64_t>::max() / left) {
        return std::numeric_limits<uint64_t>::max();
    }
    return left * right;
}

uint64_t estimatedGapCells(Coord_t ref_gap, Coord_t query_gap) {
    const uint64_t maximum = std::max<uint64_t>(ref_gap, query_gap);
    const uint64_t band = static_cast<uint64_t>(auto_band(
        static_cast<int>(ref_gap), static_cast<int>(query_gap)));
    return saturatingMultiply(maximum, saturatingAdd(band, band) + 1);
}

void observeGap(Statistics* statistics, int_t ref_gap, int_t query_gap,
                bool rejected) {
    if (!statistics || ref_gap < 0 || query_gap < 0) return;
    statistics->maximum_seen_gap = std::max<uint64_t>(
        statistics->maximum_seen_gap,
        static_cast<uint64_t>(std::max(ref_gap, query_gap)));
    if (rejected) ++statistics->long_gap_rejections;
}

}  // namespace

AnchorVec materializeClusterAnchors(
    MatchClusterVec& clusters,
    const SeqPro::ManagerVariant& ref_mgr,
    const SeqPro::ManagerVariant& qry_mgr) {
    std::sort(clusters.begin(), clusters.end(),
        [](const MatchCluster& left, const MatchCluster& right) {
            if (left.empty() || right.empty()) {
                return left.size() < right.size();
            }
            return left.front().ref_start < right.front().ref_start;
        });

    MatchClusterVec cleaned;
    cleaned.reserve(clusters.size());
    bool have_previous = false;
    ChrIndex previous_ref_chromosome = 0;
    ChrIndex previous_query_chromosome = 0;
    Strand previous_strand = FORWARD;
    Coord_t previous_ref_end = 0;
    Coord_t previous_query_low = 0;
    Coord_t previous_query_high = 0;

    const auto query_bounds = [](const MatchCluster& cluster) {
        Coord_t low = std::numeric_limits<Coord_t>::max();
        Coord_t high = 0;
        for (const auto& match : cluster) {
            const Coord_t first = match.qry_start;
            const Coord_t second = match.qry_start + len2(match);
            low = std::min(low, std::min(first, second));
            high = std::max(high, std::max(first, second));
        }
        if (low == std::numeric_limits<Coord_t>::max()) low = 0;
        return std::pair<Coord_t, Coord_t>{low, high};
    };

    for (auto& cluster : clusters) {
        if (cluster.empty()) continue;
        if (!have_previous ||
            cluster.front().ref_chr_index != previous_ref_chromosome ||
            cluster.front().qry_chr_index != previous_query_chromosome ||
            cluster.front().strand() != previous_strand) {
            MatchCluster kept = std::move(cluster);
            previous_ref_chromosome = kept.front().ref_chr_index;
            previous_query_chromosome = kept.front().qry_chr_index;
            previous_strand = kept.front().strand();
            previous_ref_end = kept.back().ref_start + len1(kept.back());
            const auto [low, high] = query_bounds(kept);
            previous_query_low = low;
            previous_query_high = high;
            cleaned.push_back(std::move(kept));
            have_previous = true;
            continue;
        }

        MatchCluster pruned;
        pruned.reserve(cluster.size());
        for (auto& match : cluster) {
            const bool reference_ok = match.ref_start >= previous_ref_end;
            const Coord_t first = match.qry_start;
            const Coord_t second = match.qry_start + len2(match);
            const Coord_t low = std::min(first, second);
            const Coord_t high = std::max(first, second);
            const bool query_ok = !intervalOverlap(
                low, high, previous_query_low, previous_query_high);
            if (reference_ok && query_ok) {
                pruned.push_back(std::move(match));
            }
        }
        if (pruned.empty()) continue;
        previous_ref_end = pruned.back().ref_start + len1(pruned.back());
        const auto [low, high] = query_bounds(pruned);
        previous_query_low = low;
        previous_query_high = high;
        cleaned.push_back(std::move(pruned));
    }
    clusters.swap(cleaned);
    MatchClusterVec().swap(cleaned);

    AnchorVec anchors;
    anchors.reserve(clusters.size());
    for (auto& cluster : clusters) {
        if (cluster.empty()) continue;
        anchors.push_back(extendClusterToAnchor(cluster, ref_mgr, qry_mgr));
        anchors.back().is_linked = false;
        MatchCluster().swap(cluster);
    }
    MatchClusterVec().swap(clusters);
    return anchors;
}

std::vector<ComponentRange> splitAnchorComponents(
    const AnchorVec& anchors,
    Statistics* statistics) {
    std::vector<ComponentRange> components;
    if (anchors.empty()) return components;
    size_t begin = 0;
    uint64_t prefix_maximum_end =
        static_cast<uint64_t>(anchors.front().ref_start) +
        anchors.front().ref_len;
    const auto append_component = [&](size_t first, size_t last) {
        const uint64_t length = last - first;
        uint64_t cost = saturatingAdd(
            length,
            saturatingMultiply(length, std::min<uint64_t>(length, 2000)));
        for (size_t index = first + 1; index < last; ++index) {
            const Anchor& previous = anchors[index - 1];
            const Anchor& current = anchors[index];
            if (previous.ref_chr_index != current.ref_chr_index ||
                previous.qry_chr_index != current.qry_chr_index ||
                previous.strand != current.strand) {
                continue;
            }
            const int64_t ref_gap = static_cast<int64_t>(current.ref_start) -
                static_cast<int64_t>(previous.ref_start + previous.ref_len);
            const int64_t query_gap = current.strand == FORWARD
                ? static_cast<int64_t>(current.qry_start) -
                    static_cast<int64_t>(
                        previous.qry_start + previous.qry_len)
                : static_cast<int64_t>(previous.qry_start) -
                    static_cast<int64_t>(current.qry_start + current.qry_len);
            if (ref_gap < 0 || query_gap < 0 ||
                !linkedGapCanBeAligned(
                    static_cast<Coord_t>(ref_gap),
                    static_cast<Coord_t>(query_gap))) {
                continue;
            }
            cost = saturatingAdd(cost, estimatedGapCells(
                static_cast<Coord_t>(ref_gap),
                static_cast<Coord_t>(query_gap)));
        }
        components.push_back({first, last, cost});
    };

    for (size_t index = 1; index < anchors.size(); ++index) {
        const Anchor& previous = anchors[index - 1];
        const Anchor& current = anchors[index];
        const bool key_changed =
            previous.ref_chr_index != current.ref_chr_index ||
            previous.qry_chr_index != current.qry_chr_index ||
            previous.strand != current.strand;
        const uint64_t current_start = current.ref_start;
        const bool separated = current_start > prefix_maximum_end &&
            current_start - prefix_maximum_end >
                kMaximumLinkedGapAlignmentLength;
        if (key_changed || separated) {
            if (!key_changed && separated && statistics) {
                ++statistics->long_gap_rejections;
                statistics->maximum_seen_gap = std::max<uint64_t>(
                    statistics->maximum_seen_gap,
                    current_start - prefix_maximum_end);
            }
            append_component(begin, index);
            begin = index;
            prefix_maximum_end =
                static_cast<uint64_t>(current.ref_start) + current.ref_len;
        } else {
            prefix_maximum_end = std::max<uint64_t>(
                prefix_maximum_end,
                static_cast<uint64_t>(current.ref_start) + current.ref_len);
        }
    }
    append_component(begin, anchors.size());
    return components;
}

AnchorPtrVec linkAnchorRange(
    AnchorVec& anchors,
    size_t begin,
    size_t end,
    const SeqPro::ManagerVariant& ref_mgr,
    const SeqPro::ManagerVariant& qry_mgr,
    Statistics* statistics) {
    if (begin > end || end > anchors.size()) {
        throw std::out_of_range("Invalid Anchor linking component range");
    }
    AnchorPtrVec output;
    if (begin == end) return output;
    std::vector<size_t> linked;
    linked.reserve(end - begin);
    constexpr size_t kNoIndex = std::numeric_limits<size_t>::max();
    constexpr int kCandidateLimit = 2000;
    constexpr int kLookBack = 2000;
    constexpr int_t kBreakLength = 200;
    ManagedSequenceSlice reference_slice;
    ManagedSequenceSlice query_slice;

    size_t current_index = begin;
    while (current_index < end) {
        Anchor& current = anchors[current_index];
        if (current.is_linked) {
            ++current_index;
            continue;
        }

        size_t best_index = kNoIndex;
        int_t best_score = std::numeric_limits<int_t>::max();
        int looked = 0;
        for (size_t index = current_index + 1;
             index < end && looked < kCandidateLimit; ++index) {
            Anchor& candidate = anchors[index];
            if (candidate.is_linked) continue;
            ++looked;
            if (statistics) ++statistics->candidate_checks;
            const int_t ref_gap = static_cast<int_t>(candidate.ref_start) -
                static_cast<int_t>(current.ref_start + current.ref_len);
            const int_t query_gap = current.strand == FORWARD
                ? static_cast<int_t>(candidate.qry_start) -
                    static_cast<int_t>(current.qry_start + current.qry_len)
                : static_cast<int_t>(current.qry_start) -
                    static_cast<int_t>(candidate.qry_start + candidate.qry_len);
            if (ref_gap < 0 || query_gap < 0) continue;
            const long greater = std::max(ref_gap, query_gap);
            const long lesser = std::min(ref_gap, query_gap);
            if (greater < kBreakLength ||
                greater - lesser <= kBreakLength) {
                best_index = index;
                break;
            }
            const int_t score = (greater << 1) - lesser;
            if (best_score > score) {
                best_score = score;
                best_index = index;
            }
        }

        bool reached = false;
        if (best_index != kNoIndex) {
            Anchor& best = anchors[best_index];
            const Coord_t ref_gap_begin = current.ref_start + current.ref_len;
            const Coord_t ref_gap_length = best.ref_start - ref_gap_begin;
            Coord_t query_gap_begin = 0;
            Coord_t query_gap_length = 0;
            if (current.strand == FORWARD) {
                query_gap_begin = current.qry_start + current.qry_len;
                query_gap_length = best.qry_start - query_gap_begin;
            } else {
                query_gap_begin = best.qry_start + best.qry_len;
                query_gap_length = current.qry_start - query_gap_begin;
            }

            const bool alignable = linkedGapCanBeAligned(
                ref_gap_length, query_gap_length);
            observeGap(statistics, ref_gap_length, query_gap_length,
                       !alignable);
            if (alignable) {
                if (statistics) {
                    statistics->sequence_extractions += 2;
                    ++statistics->direct_ksw_calls;
                    statistics->estimated_ksw_cells = saturatingAdd(
                        statistics->estimated_ksw_cells,
                        estimatedGapCells(ref_gap_length, query_gap_length));
                }
                loadSequenceSlice(ref_mgr, current.ref_chr_index,
                    ref_gap_begin, ref_gap_length, false, reference_slice);
                loadSequenceSlice(qry_mgr, current.qry_chr_index,
                    query_gap_begin, query_gap_length,
                    current.strand == REVERSE, query_slice);
                AlignmentResult gap = extendAlignKSW2Result(
                    reference_slice.view, query_slice.view,
                    2 * kBreakLength);
                if (gap.summary.reference_length == ref_gap_length &&
                    gap.summary.query_length == query_gap_length) {
                    reached = true;
                    current.ref_len = best.ref_start + best.ref_len -
                        current.ref_start;
                    if (current.strand == FORWARD) {
                        current.qry_len = best.qry_start + best.qry_len -
                            current.qry_start;
                    } else {
                        current.qry_len = current.qry_start + current.qry_len -
                            best.qry_start;
                        current.qry_start = best.qry_start;
                    }
                    current.cigar.reserve(current.cigar.size() +
                        gap.cigar.size() + best.cigar.size());
                    appendCigarMove(current.cigar, gap.cigar);
                    appendCigarMove(current.cigar, best.cigar);
                    current.alignment_length += best.alignment_length +
                        gap.summary.alignment_length;
                    current.aligned_base += gap.summary.match_length +
                        best.aligned_base;
                    best.is_linked = true;
                }
            }
        }

        if (!reached) {
            if (!linked.empty()) {
                size_t best_linked_position = kNoIndex;
                int_t linked_best_score = std::numeric_limits<int_t>::max();
                int checked = 0;
                for (size_t position = linked.size();
                     position > 0 && checked < kLookBack;
                     --position, ++checked) {
                    const Anchor& previous = anchors[linked[position - 1]];
                    if (statistics) ++statistics->candidate_checks;
                    if (previous.strand != current.strand ||
                        previous.ref_chr_index != current.ref_chr_index ||
                        previous.qry_chr_index != current.qry_chr_index) {
                        continue;
                    }
                    const int_t ref_gap = static_cast<int_t>(current.ref_start) -
                        static_cast<int_t>(
                            previous.ref_start + previous.ref_len);
                    const int_t query_gap = current.strand == FORWARD
                        ? static_cast<int_t>(current.qry_start) -
                            static_cast<int_t>(
                                previous.qry_start + previous.qry_len)
                        : static_cast<int_t>(previous.qry_start) -
                            static_cast<int_t>(
                                current.qry_start + current.qry_len);
                    if (ref_gap < 0 || query_gap < 0) continue;
                    const bool alignable = linkedGapCanBeAligned(
                        static_cast<Coord_t>(ref_gap),
                        static_cast<Coord_t>(query_gap));
                    observeGap(statistics, ref_gap, query_gap, !alignable);
                    if (!alignable) continue;
                    const long greater = std::max(ref_gap, query_gap);
                    const long lesser = std::min(ref_gap, query_gap);
                    const int_t score = (greater << 1) - lesser;
                    if (score < linked_best_score) {
                        linked_best_score = score;
                        best_linked_position = position - 1;
                    }
                }

                if (best_linked_position != kNoIndex) {
                    const Anchor& previous =
                        anchors[linked[best_linked_position]];
                    const Coord_t ref_gap_begin =
                        previous.ref_start + previous.ref_len;
                    const Coord_t ref_gap_length =
                        current.ref_start - ref_gap_begin;
                    Coord_t query_gap_begin = 0;
                    Coord_t query_gap_length = 0;
                    if (current.strand == FORWARD) {
                        query_gap_begin = previous.qry_start + previous.qry_len;
                        query_gap_length = current.qry_start - query_gap_begin;
                    } else {
                        query_gap_begin = current.qry_start + current.qry_len;
                        query_gap_length = previous.qry_start - query_gap_begin;
                    }
                    if (statistics) {
                        statistics->sequence_extractions += 2;
                        ++statistics->fallback_ksw_calls;
                        statistics->estimated_ksw_cells = saturatingAdd(
                            statistics->estimated_ksw_cells,
                            estimatedGapCells(
                                ref_gap_length, query_gap_length));
                    }
                    loadSequenceSlice(ref_mgr, current.ref_chr_index,
                        ref_gap_begin, ref_gap_length, false, reference_slice);
                    loadSequenceSlice(qry_mgr, current.qry_chr_index,
                        query_gap_begin, query_gap_length,
                        current.strand == REVERSE, query_slice);
                    AlignmentResult gap = extendAlignKSW2Result(
                        reference_slice.view, query_slice.view,
                        2 * kBreakLength);
                    if (gap.summary.reference_length == ref_gap_length &&
                        gap.summary.query_length == query_gap_length) {
                        current.ref_len += gap.summary.reference_length;
                        current.qry_len += gap.summary.query_length;
                        current.alignment_length +=
                            gap.summary.alignment_length;
                        current.aligned_base += gap.summary.match_length;
                        current.ref_start -= gap.summary.reference_length;
                        if (current.strand == FORWARD) {
                            current.qry_start -= gap.summary.query_length;
                        }
                        prependCigar(current.cigar, gap.cigar);
                    }
                }
            }

            linked.push_back(current_index);
            const Coord_t current_ref_end = current.ref_start + current.ref_len;
            const Coord_t current_query_end = current.strand == FORWARD
                ? current.qry_start + current.qry_len
                : current.qry_start;
            const size_t stop = best_index == kNoIndex ? end : best_index;
            for (size_t index = current_index + 1; index < stop; ++index) {
                Anchor& candidate = anchors[index];
                if (candidate.is_linked) continue;
                const Coord_t candidate_ref_end =
                    candidate.ref_start + candidate.ref_len;
                const Coord_t candidate_query_end = current.strand == FORWARD
                    ? candidate.qry_start + candidate.qry_len
                    : candidate.qry_start;
                if (candidate_ref_end <= current_ref_end &&
                    ((current.strand == FORWARD &&
                      candidate_query_end <= current_query_end) ||
                     (current.strand == REVERSE &&
                      candidate_query_end >= current_query_end))) {
                    candidate.is_linked = true;
                    Cigar_t().swap(candidate.cigar);
                }
            }
            while (current_index < end) {
                const Anchor& candidate = anchors[current_index];
                if (!candidate.is_linked &&
                    candidate.ref_start >= current_ref_end) {
                    break;
                }
                ++current_index;
            }
        }
    }

    output.reserve(linked.size());
    for (const size_t index : linked) {
        anchors[index].is_linked = false;
        output.push_back(std::make_shared<Anchor>(std::move(anchors[index])));
    }
    return output;
}

}  // namespace AnchorLinkDetail

namespace {

constexpr size_t kDpWindow = 5000;
constexpr size_t kMaximumRetainedDpAnchors = 131072;
std::atomic<uint64_t> dp_treap_fallback_count{0};

struct DpTreapBest {
    double value = -std::numeric_limits<double>::infinity();
    size_t index = std::numeric_limits<size_t>::max();

    bool valid() const {
        return index != std::numeric_limits<size_t>::max();
    }
};

DpTreapBest betterDpBest(DpTreapBest left, DpTreapBest right) {
    if (!left.valid()) return right;
    if (!right.valid()) return left;
    if (right.value > left.value ||
        (right.value == left.value && right.index < left.index)) {
        return right;
    }
    return left;
}

uint64_t splitmix64(uint64_t value) {
    value += 0x9e3779b97f4a7c15ULL;
    value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31);
}

struct DpTreapNode {
    long long interval_end = 0;
    size_t original_index = 0;
    uint64_t priority = 0;
    double dp = 0;
    int left = -1;
    int right = -1;
    DpTreapBest subtree_best;
};

class FixedWindowDpTreap {
public:
    explicit FixedWindowDpTreap(std::vector<DpTreapNode>& storage)
        : nodes_(storage) {}

    void reset(size_t node_count) {
        nodes_.resize(node_count);
        root_ = -1;
    }

    void insert(size_t slot, long long interval_end,
                size_t original_index, double dp) {
        auto& node = nodes_[slot];
        node.interval_end = interval_end;
        node.original_index = original_index;
        node.priority = splitmix64(original_index);
        node.dp = dp;
        node.left = -1;
        node.right = -1;
        node.subtree_best = std::isnan(dp)
            ? DpTreapBest{}
            : DpTreapBest{dp, original_index};
        root_ = insertNode(root_, static_cast<int>(slot));
    }

    void erase(long long interval_end, size_t original_index) {
        root_ = eraseNode(root_, interval_end, original_index);
    }

    DpTreapBest bestAll() const {
        return bestForNode(root_);
    }

    DpTreapBest bestEndingAtOrBefore(long long coordinate) const {
        DpTreapBest best;
        int current = root_;
        while (current >= 0) {
            const auto& node = nodes_[current];
            if (node.interval_end <= coordinate) {
                best = betterDpBest(best, bestForNode(node.left));
                if (!std::isnan(node.dp)) {
                    best = betterDpBest(
                        best, {node.dp, node.original_index});
                }
                current = node.right;
            } else {
                current = node.left;
            }
        }
        return best;
    }

private:
    bool keyLess(int left, int right) const {
        const auto& lhs = nodes_[left];
        const auto& rhs = nodes_[right];
        return lhs.interval_end < rhs.interval_end ||
            (lhs.interval_end == rhs.interval_end &&
             lhs.original_index < rhs.original_index);
    }

    DpTreapBest bestForNode(int node) const {
        return node < 0 ? DpTreapBest{} : nodes_[node].subtree_best;
    }

    void update(int node) {
        if (node < 0) return;
        auto best = std::isnan(nodes_[node].dp)
            ? DpTreapBest{}
            : DpTreapBest{nodes_[node].dp,
                          nodes_[node].original_index};
        best = betterDpBest(best, bestForNode(nodes_[node].left));
        best = betterDpBest(best, bestForNode(nodes_[node].right));
        nodes_[node].subtree_best = best;
    }

    void split(int node, int key, int& left, int& right) {
        if (node < 0) {
            left = right = -1;
            return;
        }
        if (keyLess(node, key)) {
            left = node;
            split(nodes_[node].right, key,
                  nodes_[node].right, right);
            update(left);
        } else {
            right = node;
            split(nodes_[node].left, key,
                  left, nodes_[node].left);
            update(right);
        }
    }

    int merge(int left, int right) {
        if (left < 0) return right;
        if (right < 0) return left;
        if (nodes_[left].priority > nodes_[right].priority) {
            nodes_[left].right = merge(nodes_[left].right, right);
            update(left);
            return left;
        }
        nodes_[right].left = merge(left, nodes_[right].left);
        update(right);
        return right;
    }

    int insertNode(int root, int node) {
        if (root < 0) return node;
        if (nodes_[node].priority > nodes_[root].priority) {
            split(root, node, nodes_[node].left, nodes_[node].right);
            update(node);
            return node;
        }
        if (keyLess(node, root)) {
            nodes_[root].left = insertNode(nodes_[root].left, node);
        } else {
            nodes_[root].right = insertNode(nodes_[root].right, node);
        }
        update(root);
        return root;
    }

    int eraseNode(int root, long long interval_end,
                  size_t original_index) {
        if (root < 0) return -1;
        const auto& node = nodes_[root];
        if (node.interval_end == interval_end &&
            node.original_index == original_index) {
            return merge(node.left, node.right);
        }
        if (interval_end < node.interval_end ||
            (interval_end == node.interval_end &&
             original_index < node.original_index)) {
            nodes_[root].left = eraseNode(
                nodes_[root].left, interval_end, original_index);
        } else {
            nodes_[root].right = eraseNode(
                nodes_[root].right, interval_end, original_index);
        }
        update(root);
        return root;
    }

    std::vector<DpTreapNode>& nodes_;
    int root_ = -1;
};

struct DpWorkspace {
    std::vector<double> dp;
    std::vector<int_t> pre;
    std::vector<DpTreapNode> treap_nodes;
};

thread_local DpWorkspace retained_dp_workspace;

void filterAnchorsByDpTreap(AnchorPtrVec result, bool filter_ref) {
    if (result.empty()) return;
    std::sort(result.begin(), result.end(),
        [filter_ref](const AnchorPtr& left, const AnchorPtr& right) {
            return filter_ref ? left->ref_start < right->ref_start
                              : left->qry_start < right->qry_start;
        });

    DpWorkspace temporary_workspace;
    DpWorkspace& workspace = result.size() <= kMaximumRetainedDpAnchors
        ? retained_dp_workspace : temporary_workspace;
    workspace.dp.assign(result.size(), 0);
    workspace.pre.assign(result.size(), -1);
    FixedWindowDpTreap treap(workspace.treap_nodes);
    treap.reset(std::min(result.size(), kDpWindow + 1));

    const auto interval = [&](size_t index) {
        if (filter_ref) {
            return std::pair<long long, long long>{
                static_cast<long long>(result[index]->ref_start),
                static_cast<long long>(result[index]->ref_len)};
        }
        return std::pair<long long, long long>{
            static_cast<long long>(result[index]->qry_start),
            static_cast<long long>(result[index]->qry_len)};
    };
    const auto legacyPredecessor = [&](size_t current, double score) {
        const size_t begin = current > kDpWindow
            ? current - kDpWindow : 0;
        const auto [current_start, current_length] = interval(current);
        const long long current_end = current_start + current_length;
        for (size_t previous = begin; previous < current; ++previous) {
            const auto [previous_start, previous_length] = interval(previous);
            const long long previous_end = previous_start + previous_length;
            const long long overlap = std::max(
                0LL, std::min(previous_end, current_end) -
                         std::max(previous_start, current_start));
            const long long shorter =
                std::min(previous_length, current_length);
            const double overlap_ratio = shorter > 0
                ? static_cast<double>(overlap) /
                    static_cast<double>(shorter)
                : 0.0;
            if (overlap_ratio <= 0.0) {
                const double candidate = workspace.dp[previous] +
                    score - static_cast<double>(overlap);
                if (candidate > workspace.dp[current]) {
                    workspace.dp[current] = candidate;
                    workspace.pre[current] = static_cast<int_t>(previous);
                }
            }
        }
    };

    uint64_t local_fallbacks = 0;
    for (size_t index = 0; index < result.size(); ++index) {
        const double identity = static_cast<float>(
            result[index]->aligned_base) /
            result[index]->alignment_length;
        const double score = result[index]->alignment_length *
            pow(identity, 2);
        workspace.dp[index] = score;

        if (index > kDpWindow) {
            const size_t expired = index - kDpWindow - 1;
            const auto [start, length] = interval(expired);
            treap.erase(start + length, expired);
        }

        const auto [current_start, current_length] = interval(index);
        const DpTreapBest best = current_length == 0
            ? treap.bestAll()
            : treap.bestEndingAtOrBefore(current_start);
        bool fallback = result[index]->alignment_length == 0 ||
            !std::isfinite(score) ||
            (best.valid() && !std::isfinite(best.value));
        if (best.valid() && !fallback) {
            const double candidate = best.value + score;
            const double previous_representable = std::nextafter(
                best.value, -std::numeric_limits<double>::infinity());
            fallback = previous_representable + score == candidate;
            if (!fallback && candidate > workspace.dp[index]) {
                workspace.dp[index] = candidate;
                workspace.pre[index] = static_cast<int_t>(best.index);
            }
        }
        if (fallback) {
            workspace.dp[index] = score;
            workspace.pre[index] = -1;
            legacyPredecessor(index, score);
            ++local_fallbacks;
        }

        treap.insert(index % (kDpWindow + 1),
            current_start + current_length, index,
            workspace.dp[index]);
    }
    dp_treap_fallback_count.fetch_add(
        local_fallbacks, std::memory_order_relaxed);

    double best = 0.0;
    size_t best_index = 0;
    for (size_t index = 0; index < result.size(); ++index) {
        if (workspace.dp[index] > best) {
            best = workspace.dp[index];
            best_index = index;
        }
    }
    for (int index = static_cast<int>(best_index); index >= 0;
         index = workspace.pre[index]) {
        if (filter_ref) result[index]->ref_selected = true;
        else result[index]->qry_selected = true;
        if (workspace.pre[index] == -1) break;
    }
}

void filterAnchorsByDpLegacy(AnchorPtrVec result, bool filter_ref) {
    if (result.empty()) return;
    std::sort(result.begin(), result.end(),
        [filter_ref](const AnchorPtr& left, const AnchorPtr& right) {
            return filter_ref ? left->ref_start < right->ref_start
                              : left->qry_start < right->qry_start;
        });
    std::vector<double> dp(result.size(), 0);
    std::vector<int_t> pre(result.size(), -1);
    const auto interval = [&](size_t index) {
        if (filter_ref) {
            return std::pair<long long, long long>{
                static_cast<long long>(result[index]->ref_start),
                static_cast<long long>(result[index]->ref_len)};
        }
        return std::pair<long long, long long>{
            static_cast<long long>(result[index]->qry_start),
            static_cast<long long>(result[index]->qry_len)};
    };
    for (size_t index = 0; index < result.size(); ++index) {
        const double identity = static_cast<float>(
            result[index]->aligned_base) /
            result[index]->alignment_length;
        const double score = result[index]->alignment_length *
            pow(identity, 2);
        dp[index] = score;
        const size_t begin = index > kDpWindow
            ? index - kDpWindow : 0;
        for (size_t previous = begin; previous < index; ++previous) {
            const auto [previous_start, previous_length] = interval(previous);
            const auto [current_start, current_length] = interval(index);
            const long long previous_end =
                previous_start + previous_length;
            const long long current_end = current_start + current_length;
            const long long overlap = std::max(
                0LL, std::min(previous_end, current_end) -
                         std::max(previous_start, current_start));
            const long long shorter =
                std::min(previous_length, current_length);
            const double overlap_ratio = shorter > 0
                ? static_cast<double>(overlap) /
                    static_cast<double>(shorter)
                : 0.0;
            if (overlap_ratio <= 0.0) {
                const double candidate = dp[previous] + score -
                    static_cast<double>(overlap);
                if (candidate > dp[index]) {
                    dp[index] = candidate;
                    pre[index] = static_cast<int_t>(previous);
                }
            }
        }
    }
    double best = 0.0;
    size_t best_index = 0;
    for (size_t index = 0; index < result.size(); ++index) {
        if (dp[index] > best) {
            best = dp[index];
            best_index = index;
        }
    }
    for (int index = static_cast<int>(best_index); index >= 0;
         index = pre[index]) {
        if (filter_ref) result[index]->ref_selected = true;
        else result[index]->qry_selected = true;
        if (pre[index] == -1) break;
    }
}

}  // namespace


} // namespace ramag::pairwise_detail
// RaMA-G ownership, validation and serialization bridge. Algorithm bodies above
// retain pairwise control flow; this layer neither reconnects nor refilters anchors.
#include <chrono>
namespace ramag {
namespace {
using KernelClock=std::chrono::steady_clock;
double Elapsed(KernelClock::time_point start) {
    return std::chrono::duration<double>(KernelClock::now()-start).count();
}
void KernelCheck(const PairwiseCoreOptions& o,std::string_view stage) {
    if(o.interruption_callback) o.interruption_callback(stage);
}
void KernelProgress(const PairwiseCoreOptions& o,std::string_view stage,std::uint64_t done,std::uint64_t total) {
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
Sequences CheckSequences(std::span<const SequenceRecord> records) {
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
char Complement(char base) {return pairwise_detail::BASE_COMPLEMENT[static_cast<unsigned char>(base)];}
AlignmentRecord ConvertAnchor(const pairwise_detail::Anchor& a,const Sequences& refs,const Sequences& queries) {
    AlignmentRecord r;
    r.reference_id=a.ref_chr_index; r.query_id=a.qry_chr_index;
    r.reference_begin=a.ref_start; r.reference_end=a.ref_start+a.ref_len;
    r.query_begin=a.qry_start; r.query_end=a.qry_start+a.qry_len;
    r.strand=a.strand==pairwise_detail::FORWARD?Strand::Forward:Strand::Reverse;
    const auto& ref=refs.at(r.reference_id)->bases; const auto& qry=queries.at(r.query_id)->bases;
    if(r.reference_end>ref.size() || r.query_end>qry.size() || a.ref_len==0 || a.qry_len==0) throw AlignmentError("pairwise alignment span outside sequence");
    Length rp=r.reference_begin, qp=0;
    const auto append=[&](char op,Length n){if(n==0) throw AlignmentError("pairwise zero CIGAR operation");if(!r.cigar.empty()&&r.cigar.back().operation==op)r.cigar.back().length+=n;else r.cigar.push_back({op,n});};
    for(auto unit:a.cigar) {
        const Length n=unit>>4; const auto op=unit&15;
        if(n==0 || (op!=0&&op!=1&&op!=2&&op!=7&&op!=8)) throw AlignmentError("pairwise invalid packed CIGAR");
        if(op!=1 && n>r.reference_end-rp) throw AlignmentError("pairwise reference CIGAR overrun");
        if(op!=2 && n>a.qry_len-qp) throw AlignmentError("pairwise query CIGAR overrun");
        if(op==1 || op==2) {
            append(op==1?'I':'D',n); r.edit_distance+=n;
            r.score-=40+3*static_cast<std::int64_t>(n);
            if(op==1)qp+=n;else rp+=n;
        } else {
            for(Length j=0;j<n;++j) {
                const char x=ref[rp++];
                const char y=r.strand==Strand::Forward?qry[r.query_begin+qp]:Complement(qry[r.query_end-qp-1]); ++qp;
                const bool same=x==y&&x!='N';
                if((op==7&&!same)||(op==8&&same))throw AlignmentError("pairwise EQX operation contradicts sequence");
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
} // namespace

PairwiseCoreResult AlignPairwiseCore(std::span<const SequenceRecord> references,std::span<const SequenceRecord> queries,std::span<const Seed> seeds,const PairwiseCoreOptions& o) {
    if(o.threads==0 || o.threads>static_cast<unsigned>(INT32_MAX) || o.min_cluster==0 || !std::isfinite(o.diag_factor) || o.diag_factor<0 || o.max_gap>static_cast<Length>(INT64_MAX/256) || o.diag_diff>static_cast<Length>(INT64_MAX/256)) throw AlignmentError("invalid pairwise core configuration");
    auto refs=CheckSequences(references), qrys=CheckSequences(queries);
    for(const auto& q:queries)if(static_cast<long double>(o.diag_factor)*q.size()>static_cast<long double>(INT64_MAX/2))throw AlignmentError("pairwise diagonal threshold cannot be represented");
    using Key=std::tuple<SequenceId,SequenceId,pairwise_detail::Strand>;
    std::map<Key,pairwise_detail::MatchVec> buckets;
    KernelProgress(o,"seed-merge",0,seeds.size());
    for(const auto& s:seeds) {
        KernelCheck(o,"seed-merge");
        if(!refs.contains(s.reference_id)||!qrys.contains(s.query_id))throw AlignmentError("pairwise seed has unknown sequence ID");
        const auto& r=refs.at(s.reference_id)->bases;const auto& q=qrys.at(s.query_id)->bases;
        if(s.length==0||s.reference_begin>r.size()||s.length>r.size()-s.reference_begin||s.query_begin>q.size()||s.length>q.size()-s.query_begin)throw AlignmentError("pairwise seed outside sequence");
        if(s.strand!=Strand::Forward&&s.strand!=Strand::Reverse)throw AlignmentError("pairwise invalid seed strand");
        for(Length i=0;i<s.length;++i) {
            const auto x=r[s.reference_begin+i];const auto y=s.strand==Strand::Forward?q[s.query_begin+i]:Complement(q[s.query_begin+s.length-i-1]);
            if(x=='N'||x!=y)throw AlignmentError("pairwise seed is not an exact canonical match");
        }
        const auto strand=s.strand==Strand::Forward?pairwise_detail::FORWARD:pairwise_detail::REVERSE;
        buckets[{s.reference_id,s.query_id,strand}].push_back({s.reference_id,s.reference_begin,s.query_id,s.query_begin,s.length,strand});
    }
    std::vector<pairwise_detail::MatchVec> groups;
    for(auto& [key,values]:buckets){(void)key;std::sort(values.begin(),values.end(),[](const auto& x,const auto& y){return std::tie(x.qry_start,x.ref_start,x.length)<std::tie(y.qry_start,y.ref_start,y.length);});groups.push_back(std::move(values));}
    PairwiseCoreResult result;
    result.workers=std::min(o.threads,static_cast<std::uint32_t>(std::min<std::size_t>(groups.size(),UINT32_MAX)));
    KernelProgress(o,"chaining",0,groups.size());
    auto start=KernelClock::now();
    std::vector<pairwise_detail::MatchClusterVecPtr> clusters(groups.size());
    ParallelKernel(groups.size(),o,[&](std::size_t i){KernelCheck(o,"chaining");clusters[i]=pairwise_detail::clusterChrMatch(groups[i],o.min_cluster,static_cast<std::int64_t>(o.max_gap),static_cast<std::int64_t>(o.diag_diff),o.diag_factor);});
    result.clustering_seconds=Elapsed(start);
    for(const auto& v:clusters)result.clusters+=v->size();
    KernelProgress(o,"extension",0,groups.size());
    start=KernelClock::now();
    pairwise_detail::SeqPro::ManagerVariant rm{std::make_unique<pairwise_detail::SeqPro::SequenceManager>(references)};
    pairwise_detail::SeqPro::ManagerVariant qm{std::make_unique<pairwise_detail::SeqPro::SequenceManager>(queries)};
    std::vector<pairwise_detail::AnchorPtrVec> extended(groups.size());
    ParallelKernel(groups.size(),o,[&](std::size_t i){
        KernelCheck(o,"extension");
        auto anchors=pairwise_detail::AnchorLinkDetail::materializeClusterAnchors(*clusters[i],rm,qm);
        const auto components=pairwise_detail::AnchorLinkDetail::splitAnchorComponents(anchors);
        for(const auto& component:components){KernelCheck(o,"extension");auto output=pairwise_detail::AnchorLinkDetail::linkAnchorRange(anchors,component.begin,component.end,rm,qm);extended[i].insert(extended[i].end(),std::make_move_iterator(output.begin()),std::make_move_iterator(output.end()));}
    });
    result.extension_seconds=Elapsed(start);
    KernelProgress(o,"conflict-resolution",0,extended.size());
    start=KernelClock::now();
    for(bool reference_side:{true,false}) {
        std::map<SequenceId,pairwise_detail::AnchorPtrVec> dimension;
        for(const auto& group:extended)for(const auto& anchor:group)dimension[reference_side?anchor->ref_chr_index:anchor->qry_chr_index].push_back(anchor);
        std::vector<pairwise_detail::AnchorPtrVec> tasks;
        for(auto& [id,list]:dimension){(void)id;tasks.push_back(std::move(list));}
        ParallelKernel(tasks.size(),o,[&](std::size_t i){KernelCheck(o,"conflict-resolution");pairwise_detail::filterAnchorsByDpTreap(tasks[i],reference_side);});
    }
    result.selection_seconds=Elapsed(start);
    for(const auto& group:extended)for(const auto& a:group){KernelCheck(o,"conflict-resolution");result.records.push_back({ConvertAnchor(*a,refs,qrys),a->ref_selected,a->qry_selected,a->aligned_base,a->alignment_length});}
    std::stable_sort(result.records.begin(),result.records.end(),[](const auto& a,const auto& b){const auto& x=a.record;const auto& y=b.record;return std::tie(x.query_id,x.query_begin,x.reference_id,x.reference_begin,x.strand,x.query_end,x.reference_end)<std::tie(y.query_id,y.query_begin,y.reference_id,y.reference_begin,y.strand,y.query_end,y.reference_end);});
    KernelProgress(o,"conflict-resolution",result.records.size(),result.records.size());
    return result;
}

AlignmentResult AlignPairwiseFromSeeds(const std::vector<SequenceRecord>& refs,const std::vector<SequenceRecord>& queries,const AlignmentOptions& options,std::vector<Seed> seeds,RunStatistics stats) {
    if(options.break_length!=200 || options.max_dp_cells!=4000000 || options.match_score!=2 || options.mismatch_penalty!=4 || options.gap_open_penalty!=4 || options.gap_extend_penalty!=2)throw AlignmentError("legacy extension/scoring overrides are not applicable to the pairwise core");
    PairwiseCoreOptions o{options.max_gap,options.diag_diff,options.diag_factor,options.min_cluster,options.worker_threads,options.interruption_callback,options.progress_callback};
    auto result=AlignPairwiseCore(refs,queries,seeds,o);
    AlignmentResult output;
    output.statistics=std::move(stats);auto& s=output.statistics;
    s.reference_contigs=refs.size();s.query_contigs=queries.size();
    s.reference_bases=0;s.query_bases=0;s.reference_ambiguous_bases=0;s.query_ambiguous_bases=0;
    for(const auto& r:refs){if(r.size()>UINT64_MAX-s.reference_bases)throw AlignmentError("reference statistics overflow");s.reference_bases+=r.size();s.reference_ambiguous_bases+=static_cast<Length>(std::count(r.bases.begin(),r.bases.end(),'N'));}
    for(const auto& q:queries){if(q.size()>UINT64_MAX-s.query_bases)throw AlignmentError("query statistics overflow");s.query_bases+=q.size();s.query_ambiguous_bases+=static_cast<Length>(std::count(q.bases.begin(),q.bases.end(),'N'));}
    s.selected_seed_count=seeds.size();s.chain_count=result.clusters;
    s.candidate_alignment_count=result.records.size();
    s.chaining_route="pairwise-cluster-best-chain+component-link+dual-dp-v1";
    s.chaining_requested_threads=options.worker_threads;s.chaining_worker_threads=result.workers;s.extension_worker_threads=result.workers;
    s.chain_and_extension_seconds=result.clustering_seconds+result.extension_seconds;
    s.conflict_resolution_seconds=result.selection_seconds;
    s.total_seconds=s.seed_seconds+s.chain_and_extension_seconds+s.conflict_resolution_seconds;
    std::map<SequenceId,bool> primary;
    for(auto& item:result.records){if(options.selection==AlignmentSelection::ReciprocalOneToOne&&!(item.reference_selected&&item.query_selected))continue;item.record.primary=!primary[item.record.query_id];primary[item.record.query_id]=true;output.alignments.push_back(std::move(item.record));}
    s.alignment_count=output.alignments.size();s.conflict_rejected_alignment_count=s.candidate_alignment_count-s.alignment_count;
    return output;
}
} // namespace ramag

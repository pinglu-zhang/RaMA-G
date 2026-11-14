#pragma once
#include "ramag/model.hpp"
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <string_view>
#include <vector>

namespace ramag::detail {
struct GapAmbiguity {
    bool ambiguous{};
    Length begin_column{},end_column{};
    std::uint64_t movable_gaps{};
};
// Only normalized run descriptors are allocated. Sequence and CIGAR storage is
// borrowed. The caller reserves this bound in its existing auxiliary ledger.
inline std::uint64_t GapAmbiguityWorkspaceBytes(std::size_t operations) {
    constexpr auto maximum=std::numeric_limits<std::uint64_t>::max();
    if(operations>(maximum-4096)/64)
        throw std::overflow_error("gap ambiguity workspace overflow");
    return static_cast<std::uint64_t>(operations)*64+4096;
}
namespace gap_ambiguity_detail {
struct Block {char operation{};Length length{},reference{},query{},column{};};
static_assert(sizeof(Block)<=64);
inline std::int64_t Add(std::int64_t a,std::int64_t b) {
    constexpr auto hi=std::numeric_limits<std::int64_t>::max();
    constexpr auto lo=std::numeric_limits<std::int64_t>::min();
    if((b>0&&a>hi-b)||(b<0&&a<lo-b))
        throw std::overflow_error("gap ambiguity score overflow");
    return a+b;
}
inline std::int64_t Difference(std::int64_t a,std::int64_t b) {
    constexpr auto hi=std::numeric_limits<std::int64_t>::max();
    constexpr auto lo=std::numeric_limits<std::int64_t>::min();
    if((b>0&&a<lo+b)||(b<0&&a>hi+b))
        throw std::overflow_error("gap ambiguity score difference overflow");
    return a-b;
}
} // namespace gap_ambiguity_detail

// Find equal-score translations of one affine gap within either adjacent paired
// stretch. Score(char,char) supplies the frozen substitution matrix; gap opening
// is 40 and extension is 3. Moving a gap preserves total extension cost. A move
// that reaches a same-kind neighboring gap refunds one opening (40).
//
// Bounds enclose changed PAIRED columns in the ORIGINAL CIGAR, half-open, with
// I/D included in the column coordinate system. If ties occur on both sides or
// around several gaps, their smallest enclosing interval is returned.
// check(columns) is called BEFORE each <=4096 scored-column chunk (two columns
// per displacement comparison), and check(0) checkpoints descriptor processing.
// No alternative alignment CIGAR or per-base state is allocated.
template<class Score,class Check>
GapAmbiguity FindGapTranslationAmbiguity(std::string_view r,std::string_view q,
        std::span<const CigarOp> cigar,Score score,Check check) {
    using namespace gap_ambiguity_detail;
    (void)GapAmbiguityWorkspaceBytes(cigar.size());
    check(std::uint64_t{0});
    std::vector<Block> blocks;blocks.reserve(cigar.size());
    Length rp=0,qp=0,column=0;
    constexpr auto maximum=std::numeric_limits<Length>::max();
    for(std::size_t i=0;i<cigar.size();++i) {
        if((i&4095)==0)check(std::uint64_t{0});
        const auto& op=cigar[i];const auto n=op.length;
        if(!n||(op.operation!='='&&op.operation!='X'&&op.operation!='M'&&
                op.operation!='I'&&op.operation!='D'))
            throw std::invalid_argument("gap ambiguity invalid CIGAR operation");
        if(n>maximum-column||
           (op.operation!='I'&&(rp>r.size()||n>r.size()-rp))||
           (op.operation!='D'&&(qp>q.size()||n>q.size()-qp)))
            throw std::invalid_argument("gap ambiguity CIGAR extent overflow");
        const char kind=(op.operation=='I'||op.operation=='D')?op.operation:'P';
        if(!blocks.empty()&&blocks.back().operation==kind) {
            if(n>maximum-blocks.back().length)
                throw std::overflow_error("gap ambiguity run overflow");
            blocks.back().length+=n;
        } else blocks.push_back({kind,n,rp,qp,column});
        if(op.operation!='I')rp+=n;
        if(op.operation!='D')qp+=n;
        column+=n;
    }
    if(rp!=r.size()||qp!=q.size())
        throw std::invalid_argument("gap ambiguity incomplete CIGAR consumption");
    GapAmbiguity result;
    const auto mark=[&](Length begin,Length end) {
        if(begin>=end)return;
        if(!result.ambiguous){result.begin_column=begin;result.end_column=end;}
        else {result.begin_column=std::min(result.begin_column,begin);
              result.end_column=std::max(result.end_column,end);}
        result.ambiguous=true;
    };
    for(std::size_t bi=0;bi<blocks.size();++bi) {
        if((bi&4095)==0)check(std::uint64_t{0});
        const auto& gap=blocks[bi];
        if(gap.operation=='P')continue;
        bool movable=false;
        for(const bool left:{true,false}) {
            if((left&&bi==0)||(!left&&bi+1==blocks.size()))continue;
            const auto adjacent=left?bi-1:bi+1;
            if(blocks[adjacent].operation!='P')continue;
            const auto limit=blocks[adjacent].length;
            const bool merge=left
                ? (bi>=2&&blocks[bi-2].operation==gap.operation)
                : (bi+2<blocks.size()&&blocks[bi+2].operation==gap.operation);
            std::int64_t delta=0;Length processed=0;
            while(processed<limit) {
                const auto chunk=std::min<Length>(2048,limit-processed);
                check(static_cast<std::uint64_t>(2*chunk));
                for(Length j=0;j<chunk;++j) {
                    const auto k=processed+j+1;
                    Length old_r{},old_q{},new_r{},new_q{};
                    if(left) {
                        old_r=gap.reference-k;old_q=gap.query-k;
                        new_r=gap.operation=='D'?gap.reference+gap.length-k:old_r;
                        new_q=gap.operation=='I'?gap.query+gap.length-k:old_q;
                    } else {
                        old_r=gap.reference+(gap.operation=='D'?gap.length:0)+k-1;
                        old_q=gap.query+(gap.operation=='I'?gap.length:0)+k-1;
                        new_r=gap.reference+k-1;new_q=gap.query+k-1;
                    }
                    const auto before=static_cast<std::int64_t>(score(r[old_r],q[old_q]));
                    const auto after=static_cast<std::int64_t>(score(r[new_r],q[new_q]));
                    delta=Add(delta,Difference(after,before));
                    const auto total=(merge&&k==limit)?Add(delta,40):delta;
                    if(total==0) {
                        movable=true;
                        if(left)mark(gap.column-k,gap.column);
                        else mark(gap.column+gap.length,gap.column+gap.length+k);
                    }
                }
                processed+=chunk;
            }
        }
        if(movable)++result.movable_gaps;
    }
    return result;
}
} // namespace ramag::detail

#pragma once
#include "ramag/centromere_blocks.hpp"
#include <algorithm>
#include <span>

namespace ramag::detail {
struct OrientationEndpoint {
  Position begin{},end{};
  Strand strand{Strand::Forward};
  std::uint64_t block_id{};
};
struct OrientationCoverage {Position begin{},end{};unsigned mask{};};

// Endpoints are frozen direct-support intervals and original known blocks, not
// inferred gaps. Consecutive endpoints partition the remaining unknown space;
// thus a same-direction old block is a pivot, never overwritten or extrapolated.
template<class Emit,class Check>
void VisitConnectedOrientationGaps(std::span<const OrientationEndpoint> endpoints,
    std::span<const OrientationCoverage> labels,std::string_view query,Length limit,
    Emit emit,Check check,bool bridge_unknown=false,
    std::span<const CentromereQueryInterval> forward={},
    std::span<const CentromereQueryInterval> reverse={}) {
  for(std::size_t i=0;i<endpoints.size();++i) {
    check();const auto& right=endpoints[i];
    if(right.begin>=right.end || right.end>query.size() || (i && endpoints[i-1].end>right.begin))
      throw AlignmentError("invalid frozen orientation endpoints");
    if(!i)continue;
    const auto& left=endpoints[i-1];
    if(left.strand!=right.strand || left.end==right.begin)continue;
    CentromereOrientationGap gap;
    gap.begin=left.end;gap.end=right.begin;gap.strand=right.strand;
    gap.left_block_id=left.block_id;gap.right_block_id=right.block_id;
    gap.reason="gap_too_long";
    if(gap.end-gap.begin<=limit) {
      gap.reason="contains_N";
      if(query.substr(gap.begin,gap.end-gap.begin).find('N')==std::string_view::npos) {
        gap.reason="conflicting_or_uncovered_window";
        const unsigned bit=gap.strand==Strand::Forward?1U:2U;
        auto label=std::lower_bound(labels.begin(),labels.end(),gap.begin,
          [](const auto& x,auto p){return x.end<=p;});
        Position covered=gap.begin;
        for(;label!=labels.end() && covered<gap.end;++label) {
          // Require every position, including holes between windows, to have
          // exclusively same-direction evidence. Competition still vetoes.
          if(label->begin>covered || label->mask!=bit)break;
          covered=std::min(gap.end,label->end);
        }
        if(covered==gap.end){gap.accepted=true;gap.reason="continuous_window_gap";}
        else if(bridge_unknown) {
          // Keep the old accepted set. The fallback may cross missing coverage,
          // but a competing/opposite window anywhere in the whole gap vetoes it.
          auto candidate=std::lower_bound(labels.begin(),labels.end(),gap.begin,
            [](const auto& x,auto p){return x.end<=p;});
          bool veto=false;
          for(;candidate!=labels.end() && candidate->begin<gap.end;++candidate)
            if(candidate->mask & ~bit){veto=true;break;}
          if(!veto) {
            // Canonical query unions are sorted and disjoint; alternative copies
            // never add votes. Even one opposing matched base blocks this fallback.
            const auto opposite=gap.strand==Strand::Forward?reverse:forward;
            const auto match=std::lower_bound(opposite.begin(),opposite.end(),gap.begin,
              [](const auto& x,auto p){return x.end<=p;});
            gap.reason="opposite_exact_evidence";
            if(match==opposite.end() || match->begin>=gap.end) {
              gap.accepted=true;gap.reason="same_direction_unknown_gap";
            }
          }
        }
      }
    }
    emit(gap);
  }
}
} // namespace ramag::detail

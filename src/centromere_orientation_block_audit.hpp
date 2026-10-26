#pragma once
#include "centromere_orientation_connections.hpp"
#include <algorithm>
#include <span>

namespace ramag::detail {
// Old chain envelopes remain provenance, not unconditional orientation truth.
// Audit only their internal, anchor-free long gaps. A pure opposite window cell
// retracts the old direction to unknown; it does not assign a new reference copy
// or turn the window's envelope into a new opposite-direction prediction.
// No paths or sequences are copied. Existing chains are reference-ascending;
// reverse chains can therefore be visited backwards in query order.
template<class EmitAudit,class EmitKnown,class Check>
void VisitKnownOrientationBlockAudit(std::span<const CentromereBlock> blocks,
    std::span<const OrientationCoverage> labels,Length query_length,Length gap_limit,
    EmitAudit audit,EmitKnown known,Check check) {
  if(!gap_limit)throw AlignmentError("known-block audit requires a positive gap limit");
  for(std::size_t i=0;i<labels.size();++i) {
    check();const auto& row=labels[i];
    if(row.begin>=row.end || row.end>query_length || row.mask>7U ||
       (i && labels[i-1].end>row.begin))
      throw AlignmentError("invalid known-block audit coverage cells");
  }
  for(const auto& block:blocks) {
    check();if(!block.direction_known)continue;
    if(block.query_begin>block.query_end || block.query_end>query_length ||
       block.reference_begin>block.reference_end)
      throw AlignmentError("invalid known-block audit bounds");
    const auto& anchors=block.planning_anchors;
    for(std::size_t i=0;i<anchors.size();++i) {
      check();const auto& a=anchors[i];
      // Subtraction bounds the additions below and rejects overflow as well.
      if(!a.length || a.query_begin<block.query_begin || a.query_begin>block.query_end ||
         a.length>block.query_end-a.query_begin || a.reference_begin<block.reference_begin ||
         a.reference_begin>block.reference_end || a.length>block.reference_end-a.reference_begin)
        throw AlignmentError("planning anchor lies outside its audited block");
      if(i) {
        const auto& p=anchors[i-1];
        if(p.reference_begin+p.length>a.reference_begin ||
           (block.strand==Strand::Forward ? p.query_begin+p.length>a.query_begin :
             a.query_begin+a.length>p.query_begin))
          throw AlignmentError("planning anchors are not an ordered nonoverlapping chain");
      }
    }
    if(block.query_begin==block.query_end)continue;
    Position retained_begin=block.query_begin;
    const unsigned opposite=block.strand==Strand::Forward?2U:1U;
    auto anchor=[&](std::size_t i)->const CentromereAnchor& {
      return anchors[block.strand==Strand::Forward?i:anchors.size()-1-i];
    };
    for(std::size_t i=1;i<anchors.size();++i) {
      check();const auto& left=anchor(i-1);const auto& right=anchor(i);
      const Position gap_begin=left.query_begin+left.length,gap_end=right.query_begin;
      if(gap_end-gap_begin<=gap_limit)continue;
      auto label=std::lower_bound(labels.begin(),labels.end(),gap_begin,
        [](const auto& row,auto p){return row.end<=p;});
      Position cut_begin=0,cut_end=0;bool pending=false;
      auto emit=[&] {
        if(!pending)return;
        if(retained_begin<cut_begin)
          known(CentromereOrientationKnownInterval{block.id,retained_begin,cut_begin,block.strand});
        audit(CentromereOrientationBlockAudit{block.id,cut_begin,cut_end,gap_begin,gap_end,block.strand});
        retained_begin=cut_end;pending=false;
      };
      for(;label!=labels.end() && label->begin<gap_end;++label) {
        check();if(label->mask!=opposite)continue;
        const auto begin=std::max(gap_begin,label->begin),end=std::min(gap_end,label->end);
        if(begin>=end)continue;
        if(pending && begin==cut_end)cut_end=end;
        else{emit();cut_begin=begin;cut_end=end;pending=true;}
      }
      emit();
    }
    if(retained_begin<block.query_end)
      known(CentromereOrientationKnownInterval{block.id,retained_begin,block.query_end,block.strand});
  }
}
} // namespace ramag::detail

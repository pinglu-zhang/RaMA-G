#pragma once
#include "centromere_direction_island_audit.hpp"

namespace ramag::detail {
// Exact, gap-filling alternatives over frozen raw evidence. No sequence search
// or inferred boundary is promoted to a new root.
struct ExactDirectionGapReview {
  CentromereAnchor core;
  Strand strand{Strand::Forward};
  std::uint64_t source_id{},checks{};
  bool accepted{};
  const char* reason{"ineligible_geometry"};
};
template<class Step>
ExactDirectionGapReview ReviewExactDirectionGap(
    const DirectionIslandEvidence& left,const DirectionIslandEvidence& right,
    std::span<const DirectionIslandEvidence> roots,Length rlength,Length qlength,
    std::string_view reference,std::string_view query,std::pmr::memory_resource* memory,
    Step step,std::uint64_t limit=100000) {
  using namespace direction_island_detail;
  ExactDirectionGapReview out;
  auto tick=[&] {step();if(out.checks==limit){out.reason="budget_exhausted";return false;}++out.checks;return true;};
  Validate(left,rlength,qlength);Validate(right,rlength,qlength);
  if(!memory || reference.size()!=rlength || query.size()!=qlength)
    throw AlignmentError("invalid exact direction gap input");
  if(left.strand!=right.strand)return out;
  const auto& l=left.anchor;const auto& r=right.anchor;
  const auto rb=End(l.reference_begin,l.length),re=r.reference_begin;
  const auto ub=End(Oriented(left,qlength),l.length),ue=Oriented(right,qlength);
  if(rb>=re || ub>=ue || re-rb!=ue-ub)return out;
  const auto n=re-rb;
  if(!left.source_proven || !right.source_proven || !left.source_id || !right.source_id ||
     left.source_id==right.source_id || !left.occurrence_proven || !right.occurrence_proven ||
     left.reference_count!=1 || left.query_count!=1 || right.reference_count!=1 || right.query_count!=1) {
    out.reason="flank_identity_or_occurrence_unproven";return out;
  }
  if(std::min(l.weight,r.weight)<=static_cast<long double>(n)) {
    out.reason="flanks_not_stronger_than_gap";return out;
  }
  out.strand=left.strand==Strand::Forward?Strand::Reverse:Strand::Forward;
  const DirectionIslandEvidence* chosen=nullptr;
  for(const auto& root:roots) {
    if(!tick())return out;
    Validate(root,rlength,qlength);
    if(root.strand!=out.strand || !root.source_proven || !root.occurrence_proven ||
       root.reference_count!=1 || root.query_count!=1 ||
       root.source_id==left.source_id || root.source_id==right.source_id)continue;
    // Express the opposite mapping in the BACKGROUND query orientation: it
    // has slope -1. Cropping the physical query's low end removes the HIGH
    // reference end. This must not become same-end cropping on reverse input.
    const auto& a=root.anchor;
    const auto u=left.strand==Strand::Forward?a.query_begin:qlength-a.query_begin-a.length;
    if(a.length<n || n<=a.length-n)continue; // Retain more of the raw match than is discarded.
    if(a.reference_begin>rb || End(a.reference_begin,a.length)<re ||
       u>ub || End(u,a.length)<ue)continue;
    const auto front=rb-a.reference_begin,back=End(a.reference_begin,a.length)-re;
    if(End(u,back)!=ub || End(u,a.length)-front!=ue)continue;
    if(chosen && !SameMapping(*chosen,root,qlength)) {
      out.reason="tied_opposite_mapping";return out;
    }
    if(!chosen || root.source_id<chosen->source_id)chosen=&root;
  }
  if(!chosen){out.reason="no_unique_exact_gap_filler";return out;}
  const auto qb=left.strand==Strand::Forward?ub:qlength-ue;
  out.core={rb,qb,n,chosen->anchor.weight*static_cast<double>(n)/static_cast<double>(chosen->anchor.length)};
  out.source_id=chosen->source_id;
  if(out.core.weight<=0){out.reason="nonpositive_local_gain";return out;}
  // Conservative competition: if other raw mappings can explain the complete
  // physical gap on either axis, do not force this correspondence. Raw parents
  // and clipped aliases on the chosen mapping never count as new competitors.
  std::pmr::vector<std::pair<Position,Position>> qcover(memory),rcover(memory);
  for(const auto& root:roots) {
    if(!tick())return out;
    if(SameMapping(*chosen,root,qlength))continue;
    const auto& a=root.anchor;
    const auto qa=std::max(qb,a.query_begin),qz=std::min(End(qb,n),End(a.query_begin,a.length));
    const auto ra=std::max(rb,a.reference_begin),rz=std::min(re,End(a.reference_begin,a.length));
    if(qa<qz)qcover.emplace_back(qa,qz);
    if(ra<rz)rcover.emplace_back(ra,rz);
  }
  auto covered=[&](auto& spans) {
    std::sort(spans.begin(),spans.end());Length total=0;Position end=0;
    for(const auto& [a,b]:spans){if(!tick())return n;if(b>end){total+=b-std::max(a,end);end=b;}}
    return total;
  };
  if(covered(qcover)>=n || covered(rcover)>=n) {
    if(std::string_view(out.reason)!="budget_exhausted")out.reason="competing_gap_coverage";
    return out;
  }
  auto complement=[](char c){switch(c){case 'A':return 'T';case 'C':return 'G';case 'G':return 'C';case 'T':return 'A';default:return 'N';}};
  for(Length i=0;i<n;++i) {
    if(!tick())return out;
    const char x=reference[rb+i],y=out.strand==Strand::Forward?query[qb+i]:complement(query[qb+n-1-i]);
    if(std::string_view("ACGT").find(x)==std::string_view::npos || y=='N' || x!=y){out.reason="not_canonical_exact";return out;}
  }
  out.accepted=true;out.reason="unique_exact_core_between_independent_flanks";return out;
}
} // namespace ramag::detail

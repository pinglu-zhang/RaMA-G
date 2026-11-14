#pragma once
#include "shared.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include <limits>
#include <memory_resource>
#include <span>
#include <tuple>
#include <vector>

// Related direction checks share a private module; their admission rules
// and template callback paths remain independent.

// direction_island_audit
namespace ramag::detail {
// source_id identifies the original exact mapping, never a clipped variant.
// occurrence_proven means counts are known for both complete input sequences.
struct DirectionIslandEvidence {
  CentromereAnchor anchor;
  Strand strand{Strand::Forward};
  std::uint64_t source_id{};
  bool source_proven{},occurrence_proven{};
  std::uint32_t reference_count{},query_count{};
  // Nonzero only for anchors in an actual frozen candidate planning block.
  // Index includes anchors that cannot qualify as flanks, so filtering cannot
  // manufacture adjacency. Raw evidence and restored baseline targets are 0.
  std::uint64_t selected_path_id{},selected_index{};
};
struct DirectionIslandTarget {
  std::uint64_t block_id{};
  Position reference_begin{},reference_end{},query_begin{},query_end{};
  Strand strand{Strand::Forward};
};
struct DirectionIslandAuditOptions {
  bool enabled{false};
  std::uint64_t max_checks{100000};
};
struct DirectionIslandReview {
  std::uint64_t block_id{},left_source_id{},right_source_id{},checks{};
  Position begin{},end{}; // Only this original query island can be retracted.
  Position reference_begin{},reference_end{},query_begin{},query_end{}; // Flank frame.
  Strand strand{Strand::Forward},background{Strand::Forward};
  Length reference_covered{},query_covered{},reference_unsupported{},query_unsupported{};
  double weight{},score{};
  bool eligible{},retracted{};
  const char* reason{"disabled"};
  const char* source_role{"candidate"};
};
namespace direction_island_detail {
inline Position End(Position begin,Length length) {
  if(length>std::numeric_limits<Position>::max()-begin)
    throw AlignmentError("direction island coordinate overflow");
  return begin+length;
}
inline void Validate(const DirectionIslandEvidence& e,Length rlength,Length qlength) {
  const auto& a=e.anchor;
  if(!a.length || a.reference_begin>rlength || a.length>rlength-a.reference_begin ||
     a.query_begin>qlength || a.length>qlength-a.query_begin ||
     !std::isfinite(a.weight) || a.weight<0 ||
     (e.strand!=Strand::Forward && e.strand!=Strand::Reverse))
    throw AlignmentError("invalid frozen direction island evidence");
}
inline Position Oriented(const DirectionIslandEvidence& e,Length qlength) {
  return e.strand==Strand::Forward ? e.anchor.query_begin :
    qlength-e.anchor.query_begin-e.anchor.length;
}
inline bool SameMapping(const DirectionIslandEvidence& a,const DirectionIslandEvidence& b,Length qlength) {
  if(a.strand!=b.strand)return false;
  auto diagonal=[&](const auto& e) {
    const auto r=e.anchor.reference_begin,q=Oriented(e,qlength);
    return std::pair{r>=q,r>=q?r-q:q-r};
  };
  return diagonal(a)==diagonal(b);
}
inline bool Overlap(Position a,Position b,Position c,Position d) {return a<d && c<b;}
inline bool SameProvenSource(const DirectionIslandEvidence& a,const DirectionIslandEvidence& b) {
  return a.source_proven && b.source_proven && a.source_id && a.source_id==b.source_id;
}
enum class SelectedNeighborProof {Absent,Proven,BudgetExhausted};
// A selected clipped predecessor/successor is not an independent competing
// copy merely because its discarded raw tail overlaps this flank. Prove that
// relation from the frozen path and canonical raw provenance; no size cutoff.
// This exempts only this flank comparison. The complete raw match must still
// pass the later interior-competition check.
template<class Step>
inline SelectedNeighborProof CompatibleSelectedNeighbor(const DirectionIslandEvidence& raw,
    const DirectionIslandEvidence& flank,std::span<const DirectionIslandEvidence> selected,
    Length qlength,Step step) {
  if(!raw.source_proven || !raw.source_id || !flank.source_proven || !flank.source_id ||
     raw.source_id==flank.source_id || !flank.selected_path_id || raw.strand!=flank.strand ||
     !raw.occurrence_proven || raw.reference_count!=1 || raw.query_count!=1)
    return SelectedNeighborProof::Absent;
  for(const auto& projection:selected) {
    if(!step())return SelectedNeighborProof::BudgetExhausted;
    if(!SameProvenSource(raw,projection) || projection.strand!=raw.strand ||
       projection.selected_path_id!=flank.selected_path_id ||
       !projection.occurrence_proven || projection.reference_count!=1 || projection.query_count!=1)
      continue;
    const auto i=projection.selected_index,j=flank.selected_index;
    if((i<j?j-i:i-j)!=1 || !SameMapping(raw,projection,qlength))continue;
    const auto& r=raw.anchor;const auto& p=projection.anchor;const auto& f=flank.anchor;
    if(p.length>=r.length || p.reference_begin<r.reference_begin || p.query_begin<r.query_begin ||
       End(p.reference_begin,p.length)>End(r.reference_begin,r.length) ||
       End(p.query_begin,p.length)>End(r.query_begin,r.length))continue;
    const auto pu=Oriented(projection,qlength),fu=Oriented(flank,qlength);
    const bool ordered=(End(p.reference_begin,p.length)<=f.reference_begin && End(pu,p.length)<=fu) ||
      (End(f.reference_begin,f.length)<=p.reference_begin && End(fu,f.length)<=pu);
    // Ordered nonoverlap on both axes proves every raw/flank conflict lies
    // outside the exact selected projection, i.e. in its discarded clipping.
    if(ordered)return SelectedNeighborProof::Proven;
  }
  return SelectedNeighborProof::Absent;
}
struct Frame {
  DirectionIslandEvidence left,right;
  Position rb{},re{},qb{},qe{};
};
inline bool Geometry(const DirectionIslandEvidence& left,const DirectionIslandEvidence& right,
    const DirectionIslandTarget& t,Frame& frame) {
  const auto& l=left.anchor;const auto& r=right.anchor;
  if(left.strand!=right.strand || left.strand==t.strand ||
     End(l.query_begin,l.length)>t.query_begin || r.query_begin<t.query_end)return false;
  frame.left=left;frame.right=right;frame.qb=End(l.query_begin,l.length);frame.qe=r.query_begin;
  if(left.strand==Strand::Forward){frame.rb=End(l.reference_begin,l.length);frame.re=r.reference_begin;}
  else {frame.rb=End(r.reference_begin,r.length);frame.re=l.reference_begin;}
  return frame.rb<frame.re && frame.qb<frame.qe && frame.rb<=t.reference_begin &&
    t.reference_end<=frame.re && frame.qb<=t.query_begin && t.query_end<=frame.qe;
}
inline const char* Qualification(const Frame& f) {
  if(!f.left.source_proven || !f.right.source_proven || !f.left.source_id || !f.right.source_id)
    return "independent_sources_unproven";
  if(f.left.source_id==f.right.source_id)return "shared_flank_source";
  if(!f.left.occurrence_proven || !f.right.occurrence_proven)return "flank_occurrences_unproven";
  if(f.left.reference_count!=1 || f.left.query_count!=1 ||
     f.right.reference_count!=1 || f.right.query_count!=1)return "competing_flank_copy";
  const auto width=std::max(f.re-f.rb,f.qe-f.qb);
  if(static_cast<long double>(f.left.anchor.weight)<=width ||
     static_cast<long double>(f.right.anchor.weight)<=width)return "flanks_not_stronger_than_frame";
  return nullptr;
}
inline auto FrameKey(const Frame& f) {
  return std::tuple{std::max(f.re-f.rb,f.qe-f.qb),f.re-f.rb,f.qe-f.qb,
    f.rb,f.qb,f.left.source_id,f.right.source_id};
}
struct WeightedInterval {Position begin{},end{};long double density{};};
// Independent physical coverage and max-density integration avoid duplicate
// votes from repeated records, overlapping raw matches, and clipped variants.
template<class Step>
inline bool Integrate(std::span<const WeightedInterval> intervals,
    std::pmr::memory_resource* memory,Step step,Length& coverage,long double& weight) {
  std::pmr::vector<Position> points(memory);
  for(const auto& x:intervals){if(!step())return false;points.push_back(x.begin);points.push_back(x.end);}
  std::sort(points.begin(),points.end());points.erase(std::unique(points.begin(),points.end()),points.end());
  coverage=0;weight=0;std::pmr::vector<long double> terms(memory);
  for(std::size_t i=1;i<points.size();++i) {
    long double density=0;bool covered=false;
    for(const auto& x:intervals) {
      if(!step())return false;
      if(x.begin<=points[i-1] && points[i]<=x.end){covered=true;density=std::max(density,x.density);}
    }
    if(covered){const auto n=points[i]-points[i-1];coverage+=n;terms.push_back(n*density);}
  }
  // Stable summation also makes physical-axis exchange and RC reflection exact.
  std::sort(terms.begin(),terms.end());for(const auto term:terms)weight+=term;
  return true;
}
} // namespace direction_island_detail

// All spans are frozen caller-owned inputs. The result owns values, never
// references into growing buffers. No search, sequence edit, or truth access.
template<class Check>
DirectionIslandReview AuditDirectionIsland(const DirectionIslandTarget& target,
    std::span<const DirectionIslandEvidence> selected_support,
    std::span<const DirectionIslandEvidence> frozen_flanks,
    std::span<const DirectionIslandEvidence> frozen_evidence,
    Length reference_length,Length query_length,const DirectionIslandAuditOptions& options,
    std::pmr::memory_resource* memory,Check check) {
  using namespace direction_island_detail;
  DirectionIslandReview out;out.block_id=target.block_id;out.begin=target.query_begin;out.end=target.query_end;
  out.strand=target.strand;out.background=target.strand==Strand::Forward?Strand::Reverse:Strand::Forward;
  if(!options.enabled)return out;
  if(!memory || target.reference_begin>=target.reference_end || target.reference_end>reference_length ||
     target.query_begin>=target.query_end || target.query_end>query_length ||
     (target.strand!=Strand::Forward && target.strand!=Strand::Reverse))
    throw AlignmentError("invalid direction island target or allocator");
  auto step=[&] {check();if(out.checks==options.max_checks){out.reason="budget_exhausted";return false;}++out.checks;return true;};
  for(const auto list:{selected_support,frozen_flanks,frozen_evidence})for(const auto& e:list) {
    if(!step())return out;
    Validate(e,reference_length,query_length);
  }
  Frame best;bool found=false;out.reason="no_ordered_flanks";
  for(const auto& left:frozen_flanks) {
    if(!step())return out;
    if(left.strand==target.strand || End(left.anchor.query_begin,left.anchor.length)>target.query_begin)continue;
    for(const auto& right:frozen_flanks) {
    if(!step())return out;
    if(right.strand==target.strand || right.anchor.query_begin<target.query_end)continue;
    Frame f;if(!Geometry(left,right,target,f))continue;
    if(const auto* reason=Qualification(f)){if(!found)out.reason=reason;continue;}
    if(!found || FrameKey(f)<FrameKey(best)){best=f;found=true;}
    }
  }
  if(!found)return out;
  out.left_source_id=best.left.source_id;out.right_source_id=best.right.source_id;
  out.reference_begin=best.rb;out.reference_end=best.re;out.query_begin=best.qb;out.query_end=best.qe;
  if(selected_support.empty()){out.reason="no_exact_support";return out;}
  for(const auto& e:selected_support) {
    const auto& a=e.anchor;
    if(!e.source_proven || !e.source_id){out.reason="support_source_unproven";return out;}
    if(e.source_id==best.left.source_id || e.source_id==best.right.source_id){out.reason="support_reuses_flank_source";return out;}
    if(e.strand!=target.strand || a.reference_begin<target.reference_begin || End(a.reference_begin,a.length)>target.reference_end ||
       a.query_begin<target.query_begin || End(a.query_begin,a.length)>target.query_end)
      throw AlignmentError("direction island support lies outside its target");
  }
  // A target must be one coherent physical chain. Alternative copy mappings
  // must not be accumulated into a single stronger directional explanation.
  for(std::size_t i=0;i<selected_support.size();++i)for(std::size_t j=i+1;j<selected_support.size();++j) {
    if(!step())return out;
    const auto& a=selected_support[i];const auto& b=selected_support[j];
    if(SameMapping(a,b,query_length))continue;
    if(SameProvenSource(a,b)){out.reason="inconsistent_source_mapping";return out;}
    const auto& x=a.anchor;const auto& y=b.anchor;const auto xu=Oriented(a,query_length),yu=Oriented(b,query_length);
    const bool ordered=(End(x.reference_begin,x.length)<=y.reference_begin && End(xu,x.length)<=yu) ||
      (End(y.reference_begin,y.length)<=x.reference_begin && End(yu,y.length)<=xu);
    if(!ordered){out.reason="competing_support_mapping";return out;}
  }
  // Count evidence only at its existing scope. Missing global counts cannot
  // establish a unique competing copy, and cannot justify a rejection either.
  const double strong=std::min(best.left.anchor.weight,best.right.anchor.weight);
  for(const auto& e:frozen_evidence) {
    if(!step())return out;
    const auto& a=e.anchor;
    for(const auto* flank:{&best.left,&best.right}) {
      if(SameProvenSource(e,*flank) && !SameMapping(e,*flank,query_length)) {
        out.reason="inconsistent_flank_source";return out;
      }
      const auto& f=flank->anchor;
      if(a.weight>=flank->anchor.weight && !SameMapping(e,*flank,query_length) &&
         (Overlap(a.query_begin,End(a.query_begin,a.length),f.query_begin,End(f.query_begin,f.length)) ||
          Overlap(a.reference_begin,End(a.reference_begin,a.length),f.reference_begin,End(f.reference_begin,f.length)))) {
        const auto proof=CompatibleSelectedNeighbor(e,*flank,frozen_flanks,query_length,step);
        if(proof==SelectedNeighborProof::BudgetExhausted)return out;
        if(proof!=SelectedNeighborProof::Proven){out.reason="strong_flank_competition";return out;}
      }
    }
    const bool touches=Overlap(a.query_begin,End(a.query_begin,a.length),best.qb,best.qe) ||
      Overlap(a.reference_begin,End(a.reference_begin,a.length),best.rb,best.re);
    if(!touches || a.weight<strong)continue;
    bool own=false;
    for(const auto& support:selected_support) {
      if(!step())return out;
      if(SameMapping(e,support,query_length)){own=true;break;}
    }
    // Even the raw parent of a clipped flank may continue into this frame;
    // that is strong competing evidence, not a second vote for its endpoint.
    if(!own){out.reason="strong_interior_competition";return out;}
  }
  std::pmr::vector<WeightedInterval> reference(memory),query(memory);
  for(const auto& e:selected_support) {
    if(!step())return out;
    const auto& a=e.anchor;const long double density=static_cast<long double>(a.weight)/a.length;
    reference.push_back({a.reference_begin,End(a.reference_begin,a.length),density});
    query.push_back({a.query_begin,End(a.query_begin,a.length),density});
  }
  long double rw=0,qw=0;
  if(!Integrate(reference,memory,step,out.reference_covered,rw) || !Integrate(query,memory,step,out.query_covered,qw))return out;
  // Valid disjoint match chains cover equally many positions on both axes.
  if(out.reference_covered!=out.query_covered || rw!=qw)
    throw AlignmentError("direction island union differs between input axes");
  const auto rlength=best.re-best.rb,qlength=best.qe-best.qb;
  out.reference_unsupported=rlength-out.reference_covered;out.query_unsupported=qlength-out.query_covered;
  out.weight=static_cast<double>(rw);
  if(!std::isfinite(out.weight))throw AlignmentError("direction island weight overflow");
  const auto imbalance=rlength>=qlength?rlength-qlength:qlength-rlength;
  if(imbalance>std::min(out.reference_covered,out.query_covered)) {
    out.reason="length_imbalance_not_explained";return out;
  }
  const long double score=rw-std::max(out.reference_unsupported,out.query_unsupported);
  out.score=static_cast<double>(score);out.eligible=true;out.retracted=score<=0;
  out.reason=out.retracted?"insufficient_absolute_support":"positive_absolute_support";
  return out;
}

// Context is a finite set of actual, nonoverlapping exact pieces. It explains
// only its own physical positions; unexplained frame space still penalizes the
// island. Neither a bounding box nor a duplicate raw record is extra support.
struct DirectionIslandContextDecision {
  CentromereAnchor core;
  Length explained{},reference_unsupported{},query_unsupported{};
  std::uint64_t checks{},components{};
  double score{};
  bool accepted{};
  const char* reason{"no_proven_context"};
};
namespace direction_island_context_detail {
using Range=std::pair<Position,Position>;
inline bool Touch(Position a,Length n,Position b,Length m) {return a<b+m && b<a+n;}
template<class Tick>
bool LargestPiece(const DirectionIslandEvidence& root,const DirectionIslandReview& frame,
    std::span<const Range> rmask,std::span<const Range> qmask,std::pmr::memory_resource* memory,
    Tick tick,CentromereAnchor& result) {
  const auto& a=root.anchor;const auto n=a.length;
  if(a.reference_begin>=frame.reference_end || a.reference_begin+n<=frame.reference_begin ||
     a.query_begin>=frame.query_end || a.query_begin+n<=frame.query_begin)return false;
  Length lo=std::max(a.reference_begin,frame.reference_begin)-a.reference_begin;
  Length hi=std::min(a.reference_begin+n,frame.reference_end)-a.reference_begin;
  const auto qb=std::max(a.query_begin,frame.query_begin),qe=std::min(a.query_begin+n,frame.query_end);
  const auto qlo=root.strand==Strand::Forward?qb-a.query_begin:a.query_begin+n-qe;
  const auto qhi=root.strand==Strand::Forward?qe-a.query_begin:a.query_begin+n-qb;
  lo=std::max(lo,qlo);hi=std::min(hi,qhi);if(lo>=hi)return false;
  std::pmr::vector<Range> blocked(memory);
  auto add=[&](Length b,Length e){b=std::max(b,lo);e=std::min(e,hi);if(b<e)blocked.emplace_back(b,e);};
  for(const auto& [b,e]:rmask) {
    if(!tick())return false;
    const auto x=std::max(b,a.reference_begin),y=std::min(e,a.reference_begin+n);
    if(x<y)add(x-a.reference_begin,y-a.reference_begin);
  }
  for(const auto& [b,e]:qmask) {
    if(!tick())return false;
    const auto x=std::max(b,a.query_begin),y=std::min(e,a.query_begin+n);
    if(x<y) {
      if(root.strand==Strand::Forward)add(x-a.query_begin,y-a.query_begin);
      else add(a.query_begin+n-y,a.query_begin+n-x);
    }
  }
  std::sort(blocked.begin(),blocked.end());Length cursor=lo,best=0,begin=0;
  auto consider=[&](Length b,Length e){if(e>b && e-b>best){begin=b;best=e-b;}};
  for(const auto& [b,e]:blocked){consider(cursor,b);cursor=std::max(cursor,e);}
  consider(cursor,hi);
  // This conservative end repair never turns a tiny surviving substring into
  // an independent unique root. More than half the raw exact mapping survives.
  if(!best || best<=n-best)return false;
  result={a.reference_begin+begin,root.strand==Strand::Forward?a.query_begin+begin:a.query_begin+n-begin-best,
          best,a.weight*static_cast<double>(best)/static_cast<double>(n)};
  return true;
}
}
template<class Step>
DirectionIslandContextDecision ReviewDirectionIslandContext(const DirectionIslandTarget& target,
    std::span<const DirectionIslandEvidence> support,std::span<const DirectionIslandEvidence> roots,
    const DirectionIslandReview& frame,std::string_view reference,std::string_view query,
    std::pmr::memory_resource* memory,Step step,std::uint64_t limit=100000) {
  using namespace direction_island_detail;using namespace direction_island_context_detail;
  DirectionIslandContextDecision out;
  auto tick=[&]{step();if(out.checks==limit){out.reason="budget_exhausted";return false;}++out.checks;return true;};
  if(!frame.eligible || !frame.retracted || support.size()!=1 || !memory)return out;
  const auto& own=support.front();
  if(!own.source_proven || !own.occurrence_proven || own.reference_count!=1 || own.query_count!=1) {
    out.reason="support_occurrences_unproven";return out;
  }
  if(reference.empty() || query.empty())throw AlignmentError("island context requires actual input bases");
  auto exact=[&](const CentromereAnchor& a,Strand sign) {
    if(a.reference_begin>reference.size() || a.length>reference.size()-a.reference_begin ||
       a.query_begin>query.size() || a.length>query.size()-a.query_begin)
      throw AlignmentError("context exact piece outside FASTA");
    auto comp=[](char c){return c=='A'?'T':c=='C'?'G':c=='G'?'C':c=='T'?'A':'N';};
    for(Length i=0;i<a.length;++i) {
      if(!tick())return false;
      const auto x=reference[a.reference_begin+i],y=sign==Strand::Forward?query[a.query_begin+i]:comp(query[a.query_begin+a.length-1-i]);
      if(std::string_view("ACGT").find(x)==std::string_view::npos || y=='N' || x!=y){out.reason="not_canonical_exact";return false;}
    }
    return true;
  };
  std::pmr::vector<Range> rmask(memory),qmask(memory);
  std::pmr::vector<std::size_t> candidates(memory);
  auto usable=[&](const auto& e) {
    return e.source_proven && e.occurrence_proven && e.reference_count==1 && e.query_count==1 &&
      e.source_id!=frame.left_source_id && e.source_id!=frame.right_source_id && !SameMapping(e,own,query.size());
  };
  // First mark shared boundaries against stronger independent mappings. Do
  // not retain the disputed bases merely because the selected island was exact.
  for(std::size_t i=0;i<roots.size();++i) {
    if(!tick())return out;
    const auto& e=roots[i];if(!usable(e))continue;
    Validate(e,reference.size(),query.size());
    CentromereAnchor a;
    const bool located=LargestPiece(e,frame,{}, {},memory,tick,a);
    if(std::string_view(out.reason)=="budget_exhausted")return out;
    if(e.anchor.weight>=own.anchor.weight) {
      const bool touches=Touch(e.anchor.reference_begin,e.anchor.length,own.anchor.reference_begin,own.anchor.length) ||
        Touch(e.anchor.query_begin,e.anchor.length,own.anchor.query_begin,own.anchor.length);
      if(touches && !located){out.reason="strong_context_outside_frame";return out;}
      if(touches) {
        if(!exact(a,e.strand))return out;
        rmask.emplace_back(e.anchor.reference_begin,e.anchor.reference_begin+e.anchor.length);
        qmask.emplace_back(e.anchor.query_begin,e.anchor.query_begin+e.anchor.length);
      }
    }
    if(located)candidates.push_back(i);
  }
  if(!LargestPiece(own,frame,rmask,qmask,memory,tick,out.core)) {
    if(std::string_view(out.reason)!="budget_exhausted")out.reason="no_unambiguous_island_core";
    return out;
  }
  if(!exact(out.core,own.strand))return out;
  rmask.clear();qmask.clear();
  rmask.emplace_back(out.core.reference_begin,out.core.reference_begin+out.core.length);
  qmask.emplace_back(out.core.query_begin,out.core.query_begin+out.core.length);
  std::sort(candidates.begin(),candidates.end(),[&](auto i,auto j) {
    const auto& a=roots[i].anchor;const auto& b=roots[j].anchor;
    return std::tuple{-a.weight,a.reference_begin,a.query_begin,i}<std::tuple{-b.weight,b.reference_begin,b.query_begin,j};
  });
  for(const auto i:candidates) {
    if(!tick())return out;
    const auto& e=roots[i];
    // Equal-weight conflicting explanations cannot be resolved by record order.
    bool tied=false;
    for(const auto j:candidates) {
      if(!tick())return out;
      const auto& other=roots[j];
      if(i==j || e.anchor.weight!=other.anchor.weight || SameMapping(e,other,query.size()))continue;
      if(Touch(e.anchor.reference_begin,e.anchor.length,other.anchor.reference_begin,other.anchor.length) ||
         Touch(e.anchor.query_begin,e.anchor.length,other.anchor.query_begin,other.anchor.length)){tied=true;break;}
    }
    if(tied)continue;
    CentromereAnchor a;if(!LargestPiece(e,frame,rmask,qmask,memory,tick,a)) {
      if(std::string_view(out.reason)=="budget_exhausted")return out;
      continue;
    }
    if(!exact(a,e.strand))return out;
    rmask.emplace_back(a.reference_begin,a.reference_begin+a.length);qmask.emplace_back(a.query_begin,a.query_begin+a.length);
    out.explained+=a.length;++out.components;
  }
  const auto total=out.core.length+out.explained;
  const auto rn=frame.reference_end-frame.reference_begin,qn=frame.query_end-frame.query_begin;
  if(total>rn || total>qn)throw AlignmentError("context pieces overlap a physical frame");
  out.reference_unsupported=rn-total;out.query_unsupported=qn-total;
  out.score=out.core.weight-static_cast<double>(std::max(out.reference_unsupported,out.query_unsupported));
  out.accepted=out.components && out.score>0 && out.score>frame.score;
  out.reason=out.accepted?"positive_owned_support_with_disjoint_context":"unexplained_frame_still_dominates";
  (void)target;return out;
}

namespace direction_island_detail {
class CappedResource final : public std::pmr::memory_resource {
 public:
  CappedResource(std::uint64_t limit,std::uint64_t resident):limit_(limit),live_(resident),peak_(resident) {
    if(resident>limit)throw AlignmentError("direction island resident capacity exceeds planner budget");
  }
  void Account(std::uint64_t bytes) {
    if(bytes>limit_-live_)throw AlignmentError("direction island capacity budget exceeded; final view unchanged");
    live_+=bytes;peak_=std::max(peak_,live_);
  }
  std::uint64_t peak() const{return peak_;}
 private:
  std::uint64_t limit_,live_,peak_;
  void* do_allocate(std::size_t bytes,std::size_t alignment) override {
    Account(bytes);
    try{return std::pmr::new_delete_resource()->allocate(bytes,alignment);}
    catch(...){live_-=bytes;throw;}
  }
  void do_deallocate(void* pointer,std::size_t bytes,std::size_t alignment) override {
    std::pmr::new_delete_resource()->deallocate(pointer,bytes,alignment);live_-=bytes;
  }
  bool do_is_equal(const std::pmr::memory_resource& other) const noexcept override{return this==&other;}
};
struct BatchChecksExceeded {};
inline auto MappingKey(const DirectionIslandEvidence& e,Length qlength) {
  const auto r=e.anchor.reference_begin,q=Oriented(e,qlength);
  return std::tuple{e.strand==Strand::Forward?0U:1U,r>=q,r>=q?r-q:q-r};
}
} // namespace direction_island_detail

// This is called only AFTER RepairInheritedDirections has frozen the actual
// runner-consumed view. Decisions are collected before any interval is removed.
inline void ApplyDirectionIslandAudit(Length reference_length,Length query_length,
    CentromereBlockResult& out,const CentromereBlockOptions& options,
    std::string_view reference={},std::string_view query={}) {
  using namespace direction_island_detail;
  if(!options.direction_island_audit)return;
  if(options.planning_strategy!=CentromereBlockStrategy::OrientationEvidenceV9 ||
     !options.repair_inherited_directions || !out.repair_inherited_directions)
    throw CentromereInputError("direction island audit requires completed inherited-direction repair with orientation-evidence-v9");
  if(out.direction_island_audit_applied)throw AlignmentError("direction island audit is single-pass");
  if(!reference_length || !query_length)throw AlignmentError("empty direction island input");
  const auto started=std::chrono::steady_clock::now();
  CappedResource memory(options.planner_memory_limit_bytes,out.planner_capacity_bytes);
  std::uint64_t checks=0;
  auto check=[&] {
    if(options.interruption_callback)options.interruption_callback("centromere-direction-island-audit");
    if(checks==10000000)throw BatchChecksExceeded{};
    ++checks;
  };
  std::pmr::vector<DirectionIslandReview> reviews(&memory);
  std::pmr::vector<CentromereDirectionIslandContextReview> context_reviews(&memory);
  struct Cut {Position begin,end;Strand strand;};
  std::pmr::vector<Cut> cuts(&memory);
  const auto& frozen_view=out.local_direction_repair_intervals;
  auto fully_known=[&](Position begin,Position end,Strand strand) {
    const auto it=std::lower_bound(frozen_view.begin(),frozen_view.end(),begin,
      [](const auto& row,Position p){return row.end<=p;});
    if(!out.direction_local_protection_applied)
      return it!=frozen_view.end() && it->begin<=begin && end<=it->end && it->strand==strand;
    // Protection preserves old split boundaries. They must not make a fully
    // known equal-strand span disappear from this later absolute-support audit.
    auto cursor=begin;
    for(auto row=it;row!=frozen_view.end() && cursor<end;++row) {
      if(row->begin>cursor || row->strand!=strand)return false;
      cursor=std::max(cursor,row->end);
    }
    return cursor>=end;
  };
  try {
    for(std::size_t i=0;i<frozen_view.size();++i) {
      check();const auto& row=frozen_view[i];
      if(row.begin>=row.end || row.end>query_length ||
         (row.strand!=Strand::Forward && row.strand!=Strand::Reverse) ||
         (i && frozen_view[i-1].end>row.begin))throw AlignmentError("invalid authoritative direction view");
    }
    std::pmr::vector<DirectionIslandEvidence> roots(&memory),flanks(&memory),planning(&memory);
    std::pmr::vector<std::size_t> order(&memory),offsets(&memory);
    // In v9 this store is from InitialEvidence on the two complete inputs;
    // locally searched evidence has separate stores and is never imported here.
    for(std::size_t i=0;i<out.root_candidates.size();++i) {
      check();const auto& r=out.root_candidates[i];
      DirectionIslandEvidence e{r.anchor,r.strand,i+1,true,r.reference_count>0 && r.query_count>0,r.reference_count,r.query_count};
      Validate(e,reference_length,query_length);roots.push_back(e);order.push_back(i);
    }
    std::sort(order.begin(),order.end(),[&](auto a,auto b) {
      return std::tuple{MappingKey(roots[a],query_length),roots[a].anchor.reference_begin,a}<
             std::tuple{MappingKey(roots[b],query_length),roots[b].anchor.reference_begin,b};
    });
    auto origin=[&](const CentromereAnchor& a,Strand strand) {
      DirectionIslandEvidence result{a,strand,0,false,false,0,0};Validate(result,reference_length,query_length);
      const auto key=MappingKey(result,query_length);
      auto first=std::lower_bound(order.begin(),order.end(),key,[&](auto i,const auto& k){return MappingKey(roots[i],query_length)<k;});
      auto last=std::upper_bound(first,order.end(),key,[&](const auto& k,auto i){return k<MappingKey(roots[i],query_length);});
      const DirectionIslandEvidence* match=nullptr;
      for(auto it=first;it!=last;++it) {
        check();const auto& r=roots[*it];const auto& b=r.anchor;
        if(b.reference_begin>a.reference_begin || End(a.reference_begin,a.length)>End(b.reference_begin,b.length))continue;
        if(match) {
          const auto& m=match->anchor;
          // Identical input records share one canonical source. Distinct raw
          // containers cannot be arbitrarily selected as proof of independence.
          if(m.reference_begin!=b.reference_begin || m.query_begin!=b.query_begin || m.length!=b.length ||
             match->reference_count!=r.reference_count || match->query_count!=r.query_count)return result;
          if(r.source_id<match->source_id)match=&r;
        } else match=&r;
      }
      if(match){result.source_id=match->source_id;result.source_proven=true;result.occurrence_proven=match->occurrence_proven;
        result.reference_count=match->reference_count;result.query_count=match->query_count;}
      return result;
    };
    offsets.push_back(0);
    std::uint64_t selected_path_id=0;
    for(const auto& block:out.blocks) {
      check();++selected_path_id;std::uint64_t selected_index=0;
      for(const auto& a:block.planning_anchors) {
        auto e=origin(a,block.strand);e.selected_path_id=selected_path_id;e.selected_index=selected_index++;
        planning.push_back(e);
        if(block.direction_known && fully_known(a.query_begin,End(a.query_begin,a.length),block.strand))flanks.push_back(e);
      }
      offsets.push_back(planning.size());
    }
    // Restored compact-baseline blocks are additional audit targets only.
    // They never supply flank votes or additional global candidate evidence.
    for(const auto& block:out.direction_local_protection_baseline_blocks) {
      check();
      for(const auto& a:block.planning_anchors){check();planning.push_back(origin(a,block.strand));}
      offsets.push_back(planning.size());
    }
    // The original audit has priority over every new context review. Freeze
    // ALL old verdicts before spending any of the remaining batch work budget.
    for(std::size_t i=0;i<out.blocks.size()+out.direction_local_protection_baseline_blocks.size();++i) {
      check();const bool baseline_target=i>=out.blocks.size();
      const auto& b=baseline_target?out.direction_local_protection_baseline_blocks[i-out.blocks.size()]:out.blocks[i];
      if(!b.direction_known || b.query_begin>=b.query_end)continue;
      DirectionIslandReview review;review.block_id=b.id;review.begin=b.query_begin;review.end=b.query_end;review.strand=b.strand;
      review.background=b.strand==Strand::Forward?Strand::Reverse:Strand::Forward;
      if(!fully_known(b.query_begin,b.query_end,b.strand))review.reason="final_view_not_whole_block";
      else {
        const auto support=std::span<const DirectionIslandEvidence>(planning).subspan(offsets[i],offsets[i+1]-offsets[i]);
        review=AuditDirectionIsland({b.id,b.reference_begin,b.reference_end,b.query_begin,b.query_end,b.strand},
          support,flanks,roots,reference_length,query_length,{true,100000},&memory,check);
      }
      review.source_role=baseline_target?"baseline_protection":"candidate";
      reviews.push_back(review);
      if(review.eligible && review.retracted)cuts.push_back({review.begin,review.end,review.strand});
    }
    if(options.direction_island_context_review) {
      struct ContextWorkExhausted {};
      std::uint64_t context_checks=0;
      std::pmr::vector<std::size_t> parents(&memory);
      std::pmr::vector<double> core_weights(&memory);
      auto context_check=[&] {
        if(context_checks==options.direction_island_context_check_limit || checks==10000000)
          throw ContextWorkExhausted{};
        ++context_checks;check();
      };
      bool complete=true;
      try {
        std::size_t review_index=0;
        for(std::size_t i=0;i<out.blocks.size()+out.direction_local_protection_baseline_blocks.size();++i) {
          const bool baseline_target=i>=out.blocks.size();
          const auto& b=baseline_target?out.direction_local_protection_baseline_blocks[i-out.blocks.size()]:out.blocks[i];
          if(!b.direction_known || b.query_begin>=b.query_end)continue;
          const auto current=review_index++;const auto& review=reviews[current];
          if(!review.eligible || !review.retracted)continue;
          const auto support=std::span<const DirectionIslandEvidence>(planning).subspan(offsets[i],offsets[i+1]-offsets[i]);
          const auto decision=ReviewDirectionIslandContext(
            {b.id,b.reference_begin,b.reference_end,b.query_begin,b.query_end,b.strand},
            support,roots,review,reference,query,&memory,context_check);
          parents.push_back(current);core_weights.push_back(decision.core.weight);
          context_reviews.push_back({context_reviews.size()+1,b.id,decision.checks,decision.components,
            b.query_begin,b.query_end,decision.core.reference_begin,decision.core.reference_begin+decision.core.length,
            decision.core.query_begin,decision.core.query_begin+decision.core.length,
            decision.explained,decision.explained,decision.reference_unsupported,decision.query_unsupported,
            review.score,decision.score,decision.accepted,decision.reason});
        }
      } catch(const ContextWorkExhausted&) {
        complete=false;out.direction_island_context_budget_exhausted=true;
        // Context proposals are atomic. The original cuts and verdicts above
        // remain intact, including old rejection of random short inversions.
        for(auto& row:context_reviews){row.accepted=false;row.reason="remaining_batch_checks_exhausted";}
      }
      if(complete) {
        cuts.clear();
        for(std::size_t i=0;i<context_reviews.size();++i) {
          const auto& decision=context_reviews[i];if(!decision.accepted)continue;
          auto& review=reviews[parents[i]];
          const auto begin=decision.core_query_begin,end=decision.core_query_end;
          if(begin>review.begin)cuts.push_back({review.begin,begin,review.strand});
          if(end<review.end)cuts.push_back({end,review.end,review.strand});
          review.begin=begin;review.end=end;
          review.reference_covered=review.query_covered=decision.core_reference_end-decision.core_reference_begin;
          review.reference_unsupported=decision.reference_unsupported;review.query_unsupported=decision.query_unsupported;
          review.weight=core_weights[i];review.score=decision.score;review.retracted=false;
          review.reason="positive_component_owned_support";
        }
        for(const auto& review:reviews)
          if(review.eligible && review.retracted)cuts.push_back({review.begin,review.end,review.strand});
      }
    }
  } catch(const BatchChecksExceeded&) {
    // Only an ORIGINAL audit can reach this handler. New context exhaustion
    // cannot undo previously frozen baseline verdicts.
    reviews.clear();cuts.clear();context_reviews.clear();DirectionIslandReview skip;skip.reason="batch_budget_exhausted";skip.checks=checks;reviews.push_back(skip);
  }

  std::sort(cuts.begin(),cuts.end(),[](const auto& a,const auto& b){return std::tuple{a.strand,a.begin,a.end}<std::tuple{b.strand,b.begin,b.end};});
  std::size_t count=0;
  for(const auto& cut:cuts) {
    if(count && cuts[count-1].strand==cut.strand && cut.begin<=cuts[count-1].end)cuts[count-1].end=std::max(cuts[count-1].end,cut.end);
    else cuts[count++]=cut;
  }
  cuts.resize(count);
  std::pmr::vector<CentromereOrientationKnownInterval> final(&memory);
  for(const auto& row:frozen_view) {
    Position begin=row.begin;
    for(const auto& cut:cuts) {
      if(cut.strand!=row.strand || cut.end<=begin || cut.begin>=row.end)continue;
      if(begin<cut.begin)final.push_back({row.block_id,begin,std::min(row.end,cut.begin),row.strand});
      begin=std::max(begin,cut.end);if(begin>=row.end)break;
    }
    if(begin<row.end)final.push_back({row.block_id,begin,row.end,row.strand});
  }
  // Prepare all public allocations before publishing anything. Capped pmr
  // covers scratch; these retained standard containers are explicitly charged.
  std::vector<CentromereDirectionIslandAudit> public_reviews;
  memory.Account(reviews.size()*sizeof(CentromereDirectionIslandAudit));public_reviews.reserve(reviews.size());
  for(const auto& r:reviews) {
    const auto string_bytes=std::strlen(r.reason)+std::strlen(r.source_role)+2;
    memory.Account(string_bytes);
    public_reviews.push_back({public_reviews.size()+1,r.block_id,r.left_source_id,r.right_source_id,r.checks,
      r.begin,r.end,r.reference_begin,r.reference_end,r.query_begin,r.query_end,r.strand,r.background,
      r.reference_covered,r.query_covered,r.reference_unsupported,r.query_unsupported,r.weight,r.score,r.eligible,r.retracted,r.reason,r.source_role});
    const auto& retained=public_reviews.back();
    const auto actual_string_bytes=retained.reason.capacity()+retained.source_role.capacity()+2;
    if(actual_string_bytes>string_bytes)memory.Account(actual_string_bytes-string_bytes);
  }
  std::vector<CentromereOrientationKnownInterval> public_final;
  std::vector<CentromereQueryInterval> public_unknown;
  memory.Account(final.size()*sizeof(CentromereOrientationKnownInterval));public_final.reserve(final.size());
  memory.Account((final.size()+1)*sizeof(CentromereQueryInterval));public_unknown.reserve(final.size()+1);
  Position cursor=0;
  for(const auto& row:final){public_final.push_back(row);if(cursor<row.begin)public_unknown.push_back({cursor,row.begin});cursor=row.end;}
  if(cursor<query_length)public_unknown.push_back({cursor,query_length});
  const auto old_bytes=out.local_direction_repair_intervals.capacity()*sizeof(CentromereOrientationKnownInterval)+
    out.orientation_unknown.capacity()*sizeof(CentromereQueryInterval);
  std::uint64_t retained=public_reviews.capacity()*sizeof(CentromereDirectionIslandAudit)+
    public_final.capacity()*sizeof(CentromereOrientationKnownInterval)+public_unknown.capacity()*sizeof(CentromereQueryInterval);
  for(const auto& r:public_reviews)retained+=r.reason.capacity()+r.source_role.capacity()+2;
  std::vector<CentromereDirectionIslandContextReview> public_context;
  memory.Account(context_reviews.size()*sizeof(CentromereDirectionIslandContextReview));
  public_context.reserve(context_reviews.size());
  for(auto& row:context_reviews){memory.Account(row.reason.capacity()+1);public_context.push_back(std::move(row));}
  retained+=public_context.capacity()*sizeof(CentromereDirectionIslandContextReview);
  for(const auto& row:public_context)retained+=row.reason.capacity()+1;
  out.direction_island_context_reviews=std::move(public_context);
  out.direction_island_audits=std::move(public_reviews);
  if(!cuts.empty()) {
    out.local_direction_repair_intervals=std::move(public_final);out.orientation_unknown=std::move(public_unknown);
    out.planner_capacity_bytes=out.planner_capacity_bytes>=old_bytes?out.planner_capacity_bytes-old_bytes:0;
    out.planner_capacity_bytes+=retained;
  } else {
    out.planner_capacity_bytes+=out.direction_island_audits.capacity()*sizeof(CentromereDirectionIslandAudit);
    for(const auto& r:out.direction_island_audits)out.planner_capacity_bytes+=r.reason.capacity()+r.source_role.capacity()+2;
    out.planner_capacity_bytes+=out.direction_island_context_reviews.capacity()*sizeof(CentromereDirectionIslandContextReview);
    for(const auto& r:out.direction_island_context_reviews)out.planner_capacity_bytes+=r.reason.capacity()+1;
  }
  out.planner_reserved_bytes=std::max(out.planner_reserved_bytes,memory.peak());
  out.direction_island_checks=checks;out.direction_island_audit_applied=true;
  out.direction_island_audit_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
}
} // namespace ramag::detail

// direction_plan_review
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


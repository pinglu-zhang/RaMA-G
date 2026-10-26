#include "centromere_detail.hpp"
#include <algorithm>
#include <map>
#include <set>
#include <tuple>
#if RAMAG_HAVE_CENTROMERE
#include "anchor.h"
#endif

namespace ramag::detail {
namespace {
using Box=std::tuple<Position,Position,Position,Position>;
using Diagonal=std::tuple<Strand,bool,Length>;
[[maybe_unused]] Box Key(const CentromereRefinementRegion& r){return {r.reference_begin,r.reference_end,r.query_begin,r.query_end};}
// Mapping novelty is an interval union on the SAME oriented diagonal, not
// a new (start,length) tuple. Rediscovering a submatch adds no paired column.
struct Coverage {
  std::map<Diagonal,std::map<Position,Position>> intervals;
  std::size_t segments{};
  bool Add(const CentromerePlanningAnchor& c,Length qlength) {
    const auto& a=c.anchor;const auto u=c.strand==Strand::Forward ? a.query_begin : qlength-a.query_begin-a.length;
    auto& line=intervals[{c.strand,a.reference_begin>=u,a.reference_begin>=u ? a.reference_begin-u : u-a.reference_begin}];
    auto it=line.upper_bound(a.reference_begin);
    if(it!=line.begin() && std::prev(it)->second>=a.reference_begin+a.length)return false;
    auto begin=a.reference_begin,end=a.reference_begin+a.length;
    if(it!=line.begin() && std::prev(it)->second>=begin)it=std::prev(it);
    while(it!=line.end() && it->first<=end) {
      begin=std::min(begin,it->first);end=std::max(end,it->second);it=line.erase(it);--segments;
    }
    line.emplace(begin,end);++segments;return true;
  }
  std::size_t Bytes() const {return intervals.size()*(sizeof(Diagonal)+sizeof(std::map<Position,Position>)+64)+segments*(sizeof(Position)*2+64);}
};
[[maybe_unused]] void Check(const CentromereBlockOptions& o) {
  if(o.interruption_callback)o.interruption_callback("centromere-region-refinement");
}
}
CentromereBlockResult RefineDirectionRegions(const SequenceRecord& reference,const SequenceRecord& query,
    CentromereBlockResult plan,const CentromereBlockOptions& options,Strand preferred) {
#if RAMAG_HAVE_CENTROMERE
  const auto root_count=plan.root_candidates.size();
  std::set<Box> visited;
  Coverage mappings;
  for(const auto& a:plan.root_candidates)(void)mappings.Add(a,query.size());
  std::vector<CentromereRefinementRegion> history;
  std::vector<CentromereRefinementRegion> pending;
  std::vector<CentromereLocalEvidence> provenance;
  std::uint64_t round=0,searches=0,new_mappings=0,peak=plan.planner_reserved_bytes;
  auto capacity=[&] {
    // Conservatively include node allocation overhead of the provenance maps;
    // suffix-index memory is independently guarded by the process RSS limit.
    const auto bytes=plan.root_candidates.capacity()*sizeof(CentromerePlanningAnchor)+
      history.capacity()*sizeof(CentromereRefinementRegion)+provenance.capacity()*sizeof(CentromereLocalEvidence)+
      mappings.Bytes()+visited.size()*(sizeof(Box)+64);
    peak=std::max(peak,static_cast<std::uint64_t>(bytes));
    if(bytes>options.planner_memory_limit_bytes)throw CentromereBudgetExceeded("refinement_memory_budget");
  };
  std::string partial_reason;
  try {
  for(;;) {
    Check(options);
    pending.clear();
    for(auto task:plan.refinement_regions) {
      if(task.reference_begin>=task.reference_end || task.query_begin>=task.query_end)continue;
      if(visited.insert(Key(task)).second) {
        task.id=history.size()+pending.size()+1;task.round=round+1;pending.push_back(task);
      }
    }
    if(pending.empty())break;
    ++round;capacity();
    const auto previous_count=plan.root_candidates.size();
    // One full-input index at a time. Region counts remain local; the separate
    // full-input count query never mistakes local uniqueness for global rarity.
    for(const auto sign:{Strand::Forward,Strand::Reverse}) {
      std::string reverse;
      if(sign==Strand::Reverse) {
        reverse.resize(query.size());
        for(std::size_t i=0;i<query.size();++i) {
          const char c=query.bases[i];
          reverse[query.size()-i-1]=c=='A' ? 'T' : c=='C' ? 'G' : c=='G' ? 'C' : c=='T' ? 'A' : 'N';
        }
      }
      std::vector<rama_source::SequenceInfo> sequences{
        {reference.bases,reference.name,static_cast<uint_t>(reference.size())},
        {sign==Strand::Forward ? std::string_view(query.bases) : std::string_view(reverse),query.name,static_cast<uint_t>(query.size())}};
      rama_source::AnchorFinder finder(sequences,"",options.worker_threads==1 ? 0 : options.worker_threads,
        false,false,100,[&]{Check(options);});
      // Views share the immutable reference; equal patterns at distinct repeat
      // copies reuse one exact full-input count without storing sequence copies.
      std::map<std::string_view,std::pair<std::uint32_t,std::uint32_t>> full_counts;
      for(auto& task:pending) {
        Check(options);++searches;
        if(options.progress_callback)options.progress_callback("centromere-region-refinement",task.id,pending.back().id);
        const auto u=sign==Strand::Forward ? task.query_begin : query.size()-task.query_end;
        finder.visitRegionCandidates({static_cast<uint_t>(task.reference_begin),static_cast<uint_t>(task.reference_end-task.reference_begin),
          static_cast<uint_t>(u),static_cast<uint_t>(task.query_end-task.query_begin)},false,
          [&](const rama_source::RareCandidate& c) {
            const auto& a=c.pair;
            if(!c.first_hit_level || a.second_pos<reference.size()+1)throw AlignmentError("invalid v4 first-hit source");
            const auto qu=a.second_pos-reference.size()-1;
            const auto qb=sign==Strand::Forward ? qu : query.size()-qu-a.match_length;
            if(a.first_pos<task.reference_begin || a.first_pos>task.reference_end ||
              a.match_length>task.reference_end-a.first_pos || qb<task.query_begin || qb>task.query_end ||
              a.match_length>task.query_end-qb)throw AlignmentError("v4 local anchor outside its actual search rectangle");
            CentromerePlanningAnchor candidate{{a.first_pos,qb,a.match_length,a.weight},sign,
              static_cast<std::uint32_t>(c.reference_count),static_cast<std::uint32_t>(c.query_count),static_cast<std::uint32_t>(c.level)};
            // Dedup by coordinate mapping. Original candidates stay immutable;
            // discovery of a local higher weight never rewrites their rarity.
            if(!mappings.Add(candidate,query.size()))return;
            const auto pattern=std::string_view(reference.bases).substr(a.first_pos,a.match_length);
            auto count=full_counts.find(pattern);
            if(count==full_counts.end()) {
              const auto n=finder.countFullOccurrences(a.first_pos,a.match_length);
              count=full_counts.emplace(pattern,std::pair{static_cast<std::uint32_t>(n.first),static_cast<std::uint32_t>(n.second)}).first;
            }
            const auto counts=count->second;
            const auto id=plan.root_candidates.size();
            plan.root_candidates.push_back(candidate);
            provenance.push_back({id,task.id,static_cast<std::uint32_t>(counts.first),static_cast<std::uint32_t>(counts.second)});
            ++task.new_mappings;++new_mappings;capacity();
          });
      }
    }
    for(auto& task:pending){task.status="searched";history.push_back(std::move(task));}
    pending.clear();
    capacity();
    if(plan.root_candidates.size()==previous_count)break;
    // Release the old plan and its path/index storage before rebuilding. The
    // single compact candidate store is moved, not copied, into the new sweep.
    auto candidates=std::move(plan.root_candidates);
    peak=std::max(peak,plan.planner_reserved_bytes);
    auto local_options=options;
    const auto auxiliary=history.capacity()*sizeof(CentromereRefinementRegion)+
      provenance.capacity()*sizeof(CentromereLocalEvidence)+mappings.Bytes()+visited.size()*(sizeof(Box)+64);
    std::uint64_t retained_plan=plan.blocks.capacity()*sizeof(CentromereBlock)+plan.regions.capacity()*sizeof(CentromereDirectionRegion);
    for(const auto& b:plan.blocks)retained_plan+=b.planning_anchors.capacity()*sizeof(CentromereAnchor);
    try {
      if(auxiliary+retained_plan>=options.planner_memory_limit_bytes)throw CentromereBudgetExceeded("refinement_replan_memory_budget");
      local_options.planner_memory_limit_bytes-=auxiliary+retained_plan;
      auto next=ReplanDirectionRegions(reference.size(),query.size(),candidates,local_options,preferred);
      plan=std::move(next);
    } catch(...) {plan.root_candidates=std::move(candidates);throw;}
    peak=std::max(peak,plan.planner_reserved_bytes+auxiliary+retained_plan);
  }
  } catch(const CentromereBudgetExceeded& e) {partial_reason=e.what();}
    catch(const std::bad_alloc&) {partial_reason="refinement_allocation_failed";}
  if(!partial_reason.empty()) {
    for(auto& task:pending){task.status="partial_search";history.push_back(std::move(task));}
    plan.refinement_partial_reason=std::move(partial_reason);
  }
  plan.root_candidate_count=root_count;plan.local_candidates=std::move(provenance);
  plan.refinement_regions=std::move(history);plan.refinement_rounds=round;plan.refinement_searches=searches;
  plan.refinement_new_mappings=new_mappings;plan.planner_reserved_bytes=std::max(peak,plan.planner_reserved_bytes);
  return plan;
#else
  (void)reference;(void)query;(void)plan;(void)options;(void)preferred;
  throw DependencyError("direction refinement requires source RaMA backend");
#endif
}
} // namespace ramag::detail

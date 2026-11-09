#pragma once
#include "centromere_direction_local_protection.hpp"
#include "centromere_direction_island_audit.hpp"
#include <chrono>
#include <cstring>
#include <limits>
#include <tuple>

namespace ramag::detail {
struct FrozenLocalDirectionProtection {
  std::vector<CentromereBlock> blocks;
  std::vector<CentromereOrientationKnownInterval> known;
  std::uint64_t capacity_bytes{};
};
inline std::uint64_t LocalProtectionBlockCapacity(std::span<const CentromereBlock> blocks,std::size_t capacity) {
  std::uint64_t bytes=capacity*sizeof(CentromereBlock);
  for(const auto& b:blocks)bytes+=b.planning_anchors.capacity()*sizeof(CentromereAnchor)+b.status.capacity()+b.reason.capacity()+2;
  return bytes;
}
// Only selected block values and audited known intervals survive; the caller
// moves raw root candidates back before destroying the rest of the baseline.
inline FrozenLocalDirectionProtection FreezeLocalDirectionProtectionBaseline(CentromereBlockResult& baseline) {
  if(baseline.strategy!=CentromereBlockStrategy::OrientationEvidenceV9 || !baseline.orientation_audit_known_blocks)
    throw AlignmentError("local direction protection requires an audited compact v9 baseline");
  FrozenLocalDirectionProtection frozen;
  frozen.blocks=std::move(baseline.blocks);frozen.known=std::move(baseline.orientation_known_intervals);
  frozen.blocks.erase(std::remove_if(frozen.blocks.begin(),frozen.blocks.end(),[](const auto& b){return !b.direction_known;}),frozen.blocks.end());
  std::sort(frozen.blocks.begin(),frozen.blocks.end(),[](const auto& a,const auto& b){return a.id<b.id;});
  std::sort(frozen.known.begin(),frozen.known.end(),[](const auto& a,const auto& b){return std::pair{a.begin,a.end}<std::pair{b.begin,b.end};});
  frozen.capacity_bytes=LocalProtectionBlockCapacity(frozen.blocks,frozen.blocks.capacity())+
    frozen.known.capacity()*sizeof(CentromereOrientationKnownInterval);
  return frozen;
}

inline void ApplyLocalDirectionProtection(Length query_length,FrozenLocalDirectionProtection& baseline,
    CentromereBlockResult& out,const CentromereBlockOptions& options) {
  if(!options.direction_local_protection)return;
  if(!options.direction_island_audit || options.planning_strategy!=CentromereBlockStrategy::OrientationEvidenceV9 || !options.repair_inherited_directions ||
     !out.repair_inherited_directions || out.direction_island_audit_applied)
    throw CentromereInputError("local direction protection requires enabled island audit after inherited-direction repair");
  if(out.direction_local_protection_applied)throw AlignmentError("local direction protection is single-pass");
  const auto started=std::chrono::steady_clock::now();
  if(baseline.capacity_bytes>options.planner_memory_limit_bytes ||
     out.planner_capacity_bytes>options.planner_memory_limit_bytes-baseline.capacity_bytes)
    throw AlignmentError("local direction protection resident capacity exceeds planner budget");
  direction_island_detail::CappedResource memory(options.planner_memory_limit_bytes,out.planner_capacity_bytes+baseline.capacity_bytes);
  std::uint64_t checks=0;
  auto check=[&]{if(options.interruption_callback)options.interruption_callback("centromere-local-direction-protection");++checks;};
  using Geometry=std::tuple<Strand,Position,Position,Position,Position>;
  auto geometry=[](const CentromereBlock& b){return Geometry{b.strand,b.reference_begin,b.reference_end,b.query_begin,b.query_end};};
  std::pmr::vector<Geometry> geometries(&memory);
  std::pmr::vector<LocalDirectionProtectionEvidence> selected(&memory);
  for(const auto& b:out.blocks) {
    check();if(!b.direction_known)continue;geometries.push_back(geometry(b));
    for(const auto& a:b.planning_anchors) {
      check();
      if(!a.length || a.query_begin>query_length || a.length>query_length-a.query_begin)
        throw AlignmentError("invalid selected local direction protection anchor");
      selected.push_back({a.query_begin,a.query_begin+a.length,b.strand});
    }
  }
  std::sort(geometries.begin(),geometries.end());
  auto block_at=[&](std::uint64_t id)->const CentromereBlock& {
    const auto it=std::lower_bound(baseline.blocks.begin(),baseline.blocks.end(),id,[](const auto& b,std::uint64_t value){return b.id<value;});
    if(it==baseline.blocks.end() || it->id!=id)throw AlignmentError("missing compact baseline block for direction protection");
    return *it;
  };
  std::pmr::vector<LocalDirectionProtectionTarget> targets(&memory);
  for(const auto& row:baseline.known) {
    check();const auto& b=block_at(row.block_id);
    const bool changed=!std::binary_search(geometries.begin(),geometries.end(),geometry(b));
    targets.push_back(FreezeLocalDirectionProtectionTarget(row,b,query_length,changed,&memory,check));
  }
  const auto result=ProtectLocalDirectionLabels(targets,selected,out.local_direction_repair_intervals,query_length,true,&memory,check);
  // Allocate all retained public records before changing the authoritative view.
  std::vector<CentromereDirectionProtectionReview> reviews;
  memory.Account(result.reviews.size()*sizeof(CentromereDirectionProtectionReview));reviews.reserve(result.reviews.size());
  std::pmr::vector<std::uint64_t> changed_blocks(&memory);
  for(const auto& r:result.reviews) {
    check();memory.Account(std::strlen(r.reason)+1);
    reviews.push_back({reviews.size()+1,r.block_id,r.begin,r.end,r.strand,r.old_matched,r.counter_matched,r.candidate_same,
      r.candidate_opposite,r.candidate_unknown,r.restored,r.protected_interval,r.region_changed,r.reason});
    if(r.region_changed)changed_blocks.push_back(r.block_id);
  }
  std::sort(changed_blocks.begin(),changed_blocks.end());
  changed_blocks.erase(std::unique(changed_blocks.begin(),changed_blocks.end()),changed_blocks.end());
  const auto regions_changed=static_cast<std::uint64_t>(changed_blocks.size());
  std::vector<CentromereOrientationKnownInterval> final,protected_intervals;
  memory.Account(result.intervals.size()*sizeof(CentromereOrientationKnownInterval));final.reserve(result.intervals.size());
  for(const auto& r:result.intervals){check();final.push_back(r);}
  memory.Account(result.protected_intervals.size()*sizeof(CentromereOrientationKnownInterval));protected_intervals.reserve(result.protected_intervals.size());
  for(const auto& r:result.protected_intervals){check();protected_intervals.push_back(r);}
  std::vector<CentromereQueryInterval> unknown;
  memory.Account((final.size()+1)*sizeof(CentromereQueryInterval));unknown.reserve(final.size()+1);
  Position cursor=0;
  for(const auto& r:final){check();if(cursor<r.begin)unknown.push_back({cursor,r.begin});cursor=r.end;}
  if(cursor<query_length)unknown.push_back({cursor,query_length});
  std::vector<CentromereBlock> audit_blocks;
  memory.Account(protected_intervals.size()*sizeof(CentromereBlock));audit_blocks.reserve(protected_intervals.size());
  for(const auto& row:protected_intervals) {
    check();const auto& old=block_at(row.block_id);CentromereBlock b;
    b.id=old.id;b.query_begin=row.begin;b.query_end=row.end;b.strand=row.strand;b.direction_known=true;
    b.status="baseline_protection";b.reason="frozen_compact_baseline_target";
    memory.Account(b.status.capacity()+b.reason.capacity()+2);
    memory.Account(old.planning_anchors.size()*sizeof(CentromereAnchor));b.planning_anchors.reserve(old.planning_anchors.size());
    b.reference_begin=std::numeric_limits<Position>::max();b.reference_end=0;
    for(const auto& a:old.planning_anchors) {
      check();const auto end=a.query_begin+a.length;
      if(end<=row.begin || a.query_begin>=row.end)continue;
      // The baseline audit only cuts anchor-free gaps. A partial anchor would
      // need density-aware weight integration, so never invent a scaled weight.
      if(a.query_begin<row.begin || end>row.end)throw AlignmentError("compact protection boundary cuts an actual baseline anchor");
      b.planning_anchors.push_back(a);b.reference_begin=std::min(b.reference_begin,a.reference_begin);
      b.reference_end=std::max(b.reference_end,a.reference_begin+a.length);
      (b.strand==Strand::Forward?b.forward_support:b.reverse_support)+=a.length;
    }
    if(b.planning_anchors.empty()){b.reference_begin=old.reference_begin;b.reference_end=old.reference_end;}
    audit_blocks.push_back(std::move(b));
  }
  const auto old_bytes=out.local_direction_repair_intervals.capacity()*sizeof(CentromereOrientationKnownInterval)+
    out.orientation_unknown.capacity()*sizeof(CentromereQueryInterval);
  std::uint64_t retained=reviews.capacity()*sizeof(CentromereDirectionProtectionReview)+
    final.capacity()*sizeof(CentromereOrientationKnownInterval)+protected_intervals.capacity()*sizeof(CentromereOrientationKnownInterval)+
    unknown.capacity()*sizeof(CentromereQueryInterval)+LocalProtectionBlockCapacity(audit_blocks,audit_blocks.capacity());
  for(const auto& r:reviews)retained+=r.reason.capacity()+1;
  out.direction_local_protection_reviews=std::move(reviews);
  out.direction_local_protection_intervals=std::move(protected_intervals);
  out.direction_local_protection_baseline_blocks=std::move(audit_blocks);
  out.local_direction_repair_intervals=std::move(final);out.orientation_unknown=std::move(unknown);
  out.planner_capacity_bytes=out.planner_capacity_bytes>=old_bytes?out.planner_capacity_bytes-old_bytes:0;
  out.planner_capacity_bytes+=retained;out.planner_reserved_bytes=std::max(out.planner_reserved_bytes,memory.peak());
  out.direction_local_protection_regions_changed=regions_changed;out.direction_local_protection_checks=checks;
  out.direction_local_protection_applied=true;
  out.direction_local_protection_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
  // Release the old cache before island auditing; only the necessary old audit
  // targets survive as retained public records already charged above.
  baseline=FrozenLocalDirectionProtection{};
}
} // namespace ramag::detail

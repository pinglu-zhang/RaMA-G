#pragma once
#include "island.hpp"
#include <algorithm>
#include <array>
#include <chrono>
#include <cstring>
#include <limits>
#include <memory_resource>
#include <span>
#include <tuple>
#include <vector>

// Related direction checks share a private module; their admission rules
// and template callback paths remain independent.

// direction_local_protection
namespace ramag::detail {
// This view is the overlap-OFF Replan result on the same frozen root candidates,
// not the result of repeating the complete discovery/local-search pipeline.
struct LocalDirectionProtectionTarget {
  std::uint64_t block_id{};
  Position begin{},end{};
  Strand strand{Strand::Forward};
  Length old_matched{};
  bool region_changed{};
};
struct LocalDirectionProtectionEvidence {
  Position begin{},end{}; // Original physical query coordinates, for either sign.
  Strand strand{Strand::Forward};
};
struct LocalDirectionProtectionReview {
  std::uint64_t block_id{};
  Position begin{},end{};
  Strand strand{Strand::Forward};
  Length old_matched{},counter_matched{},candidate_same{},candidate_opposite{},candidate_unknown{},restored{};
  bool protected_interval{},region_changed{};
  const char* reason{"unchanged_label"};
};
struct LocalDirectionProtectionResult {
  std::pmr::vector<LocalDirectionProtectionReview> reviews;
  std::pmr::vector<CentromereOrientationKnownInterval> intervals,protected_intervals;
  explicit LocalDirectionProtectionResult(std::pmr::memory_resource* memory):
    reviews(memory),intervals(memory),protected_intervals(memory){}
};
namespace local_direction_protection_detail {
inline bool ValidStrand(Strand s){return s==Strand::Forward || s==Strand::Reverse;}
inline unsigned Sign(Strand s){return s==Strand::Forward?0U:1U;}
struct Range {Position begin{},end{};};
template<class Check>
inline void Merge(std::pmr::vector<Range>& ranges,Check check) {
  std::sort(ranges.begin(),ranges.end(),[](const auto& a,const auto& b){return std::pair{a.begin,a.end}<std::pair{b.begin,b.end};});
  std::size_t n=0;
  for(const auto x:ranges) {
    check();
    if(n && x.begin<=ranges[n-1].end)ranges[n-1].end=std::max(ranges[n-1].end,x.end);
    else ranges[n++]=x;
  }
  ranges.resize(n);
}
template<class Check>
inline Length Covered(std::span<const Range> ranges,Position begin,Position end,Check check) {
  Length total=0;
  auto it=std::lower_bound(ranges.begin(),ranges.end(),begin,[](const auto& x,Position p){return x.end<=p;});
  for(;it!=ranges.end() && it->begin<end;++it) {
    check();total+=std::min(it->end,end)-std::max(it->begin,begin);
  }
  return total;
}
} // namespace local_direction_protection_detail

// Compute the old support from actual selected query matches, with no use of
// match weights or envelope spans. Duplicate/overlapping copies count once.
template<class Check>
LocalDirectionProtectionTarget FreezeLocalDirectionProtectionTarget(
    const CentromereOrientationKnownInterval& known,const CentromereBlock& block,
    Length query_length,bool region_changed,std::pmr::memory_resource* memory,Check check) {
  using namespace local_direction_protection_detail;
  if(!memory || !block.direction_known || known.block_id!=block.id || known.strand!=block.strand ||
     !ValidStrand(known.strand) || known.begin>=known.end || known.end>query_length ||
     known.begin<block.query_begin || known.end>block.query_end)
    throw AlignmentError("invalid local direction protection baseline target");
  std::pmr::vector<Range> actual(memory);
  for(const auto& a:block.planning_anchors) {
    check();
    if(!a.length || a.query_begin>query_length || a.length>query_length-a.query_begin ||
       a.query_begin<block.query_begin || a.query_begin+a.length>block.query_end)
      throw AlignmentError("invalid local direction protection baseline anchor");
    const auto begin=std::max(known.begin,a.query_begin),end=std::min(known.end,a.query_begin+a.length);
    if(begin<end)actual.push_back({begin,end});
  }
  Merge(actual,check);
  return {block.id,known.begin,known.end,known.strand,Covered(actual,known.begin,known.end,check),region_changed};
}

// One frozen comparison; no reference positions, sequences, search, truth or
// generated overlay is used as new evidence. The caller applies island audit
// afterwards, and must never invoke this protection again after that audit.
template<class Check>
LocalDirectionProtectionResult ProtectLocalDirectionLabels(
    std::span<const LocalDirectionProtectionTarget> targets,
    std::span<const LocalDirectionProtectionEvidence> selected,
    std::span<const CentromereOrientationKnownInterval> candidate,
    Length query_length,bool enabled,std::pmr::memory_resource* memory,Check check) {
  using namespace local_direction_protection_detail;
  if(!memory)throw AlignmentError("missing local direction protection allocator");
  LocalDirectionProtectionResult out(memory);
  if(!enabled){out.intervals.assign(candidate.begin(),candidate.end());return out;}
  if(!query_length)throw AlignmentError("empty local direction protection input");
  for(std::size_t i=0;i<targets.size();++i) {
    check();const auto& t=targets[i];
    if(!t.block_id || t.begin>=t.end || t.end>query_length || !ValidStrand(t.strand) ||
       t.old_matched>t.end-t.begin || (i && targets[i-1].end>t.begin))
      throw AlignmentError("invalid local direction protection frozen targets");
  }
  for(std::size_t i=0;i<candidate.size();++i) {
    check();const auto& row=candidate[i];
    if(row.begin>=row.end || row.end>query_length || !ValidStrand(row.strand) ||
       (i && candidate[i-1].end>row.begin))
      throw AlignmentError("invalid local direction protection candidate view");
  }
  std::array<std::pmr::vector<Range>,2> actual{std::pmr::vector<Range>(memory),std::pmr::vector<Range>(memory)};
  for(const auto& a:selected) {
    check();
    if(a.begin>=a.end || a.end>query_length || !ValidStrand(a.strand))
      throw AlignmentError("invalid local direction protection selected evidence");
    actual[Sign(a.strand)].push_back({a.begin,a.end});
  }
  for(auto& intervals:actual)Merge(intervals,check);
  for(const auto& t:targets) {
    check();
    LocalDirectionProtectionReview review;review.block_id=t.block_id;review.begin=t.begin;review.end=t.end;
    review.strand=t.strand;review.old_matched=t.old_matched;review.region_changed=t.region_changed;
    // The complete original W is evaluated once, including all its anchors.
    // We never re-evaluate a small unsupported boundary fragment in isolation.
    review.counter_matched=Covered(actual[1U-Sign(t.strand)],t.begin,t.end,check);
    auto it=std::lower_bound(candidate.begin(),candidate.end(),t.begin,[](const auto& row,Position p){return row.end<=p;});
    for(;it!=candidate.end() && it->begin<t.end;++it) {
      check();const auto n=std::min(it->end,t.end)-std::max(it->begin,t.begin);
      (it->strand==t.strand?review.candidate_same:review.candidate_opposite)+=n;
    }
    review.candidate_unknown=t.end-t.begin-review.candidate_same-review.candidate_opposite;
    const auto changed=t.end-t.begin-review.candidate_same;
    if(changed) {
      if(review.counter_matched>review.old_matched)review.reason="strictly_stronger_local_counterevidence";
      else {
        review.protected_interval=true;review.restored=changed;review.reason="counterevidence_not_strictly_stronger";
        out.protected_intervals.push_back({t.block_id,t.begin,t.end,t.strand});
      }
    }
    out.reviews.push_back(review);
  }
  // Subtract the complete protected W from candidate labels, then insert that
  // W exactly once. Adjacent equal-strand W remain distinct original records.
  for(const auto& row:candidate) {
    check();auto begin=row.begin;
    auto it=std::lower_bound(out.protected_intervals.begin(),out.protected_intervals.end(),begin,
      [](const auto& cut,Position p){return cut.end<=p;});
    for(;it!=out.protected_intervals.end() && it->begin<row.end;++it) {
      check();
      if(begin<it->begin)out.intervals.push_back({row.block_id,begin,std::min(row.end,it->begin),row.strand});
      begin=std::max(begin,it->end);if(begin>=row.end)break;
    }
    if(begin<row.end)out.intervals.push_back({row.block_id,begin,row.end,row.strand});
  }
  for(const auto& row:out.protected_intervals) {
    check();
    // Zero means an overlay, not a reference to a candidate macro-block ID.
    // Its actual baseline block provenance remains in the review/protection log.
    out.intervals.push_back({0,row.begin,row.end,row.strand});
  }
  std::sort(out.intervals.begin(),out.intervals.end(),[](const auto& a,const auto& b){return std::pair{a.begin,a.end}<std::pair{b.begin,b.end};});
  return out;
}
} // namespace ramag::detail

// direction_local_protection_integration
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

// orientation_sparse_pivot
namespace ramag::detail {
// One frozen graph per original long gap. Its vertices are actual exact-query
// unions and the two surviving old endpoints, never inferred interval edges.
// This is direction evidence, not a chain of compatible reference-copy mappings.
template<class EmitPiece,class EmitGap,class Check>
void VisitSparseOrientationPivots(std::span<const OrientationEndpoint> endpoints,
    std::span<const OrientationCoverage> labels,std::string_view query,Length limit,
    std::span<const CentromereQueryInterval> known_forward,
    std::span<const CentromereQueryInterval> known_reverse,
    std::span<const CentromereQueryInterval> forward,
    std::span<const CentromereQueryInterval> reverse,
    std::span<const CentromereQueryInterval> occupied,std::pmr::memory_resource* memory,
    EmitPiece emit_piece,EmitGap emit_gap,Check check) {
  if(!limit || !memory)throw AlignmentError("sparse-pivot recovery requires a positive limit and workspace");
  auto validate=[&](auto rows) {
    for(std::size_t i=0;i<rows.size();++i) {
      check();const auto& row=rows[i];
      if(row.begin>=row.end || row.end>query.size() || (i && rows[i-1].end>row.begin))
        throw AlignmentError("invalid sparse-pivot interval index");
    }
  };
  validate(endpoints);validate(labels);validate(known_forward);validate(known_reverse);
  validate(forward);validate(reverse);validate(occupied);
  for(const auto& row:labels)if(row.mask>7U)throw AlignmentError("invalid sparse-pivot window mask");
  auto contains=[](auto rows,Position p) {
    const auto it=std::lower_bound(rows.begin(),rows.end(),p,[](const auto& row,auto x){return row.end<=x;});
    return it!=rows.end() && it->begin<=p;
  };
  using Interval=CentromereQueryInterval;
  std::pmr::vector<Interval> vertices(memory);
  struct Edge {bool pass{};std::string_view reason;};
  std::pmr::vector<Edge> edges(memory);
  for(std::size_t i=1;i<endpoints.size();++i) {
    check();const auto& left=endpoints[i-1];const auto& right=endpoints[i];
    if(left.strand!=right.strand || right.begin-left.end<=limit)continue;
    if(left.strand!=Strand::Forward && left.strand!=Strand::Reverse)
      throw AlignmentError("invalid sparse-pivot endpoint strand");
    const Position begin=left.end,end=right.begin;
    const auto own=left.strand==Strand::Forward?forward:reverse;
    const auto opposite=left.strand==Strand::Forward?reverse:forward;
    const auto known=left.strand==Strand::Forward?known_forward:known_reverse;
    const unsigned bit=left.strand==Strand::Forward?1U:2U;
    auto decision=[&](Position b,Position e,bool accepted,std::string_view reason) {
      emit_gap(CentromereSparseOrientationGap{b,e,begin,end,left.strand,accepted,std::string(reason)});
    };
    if(!contains(known,begin-1) || !contains(own,begin-1) ||
       !contains(known,end) || !contains(own,end)) {
      decision(begin,end,false,"original_pivot_not_supported");continue;
    }
    auto barrier=[&](Position b,Position e)->std::string_view {
      if(b==e)return {};
      auto label=std::lower_bound(labels.begin(),labels.end(),b,[](const auto& row,auto p){return row.end<=p;});
      for(;label!=labels.end() && label->begin<e;++label) {
        check();if(label->mask & ~bit)return "opposing_or_competitive_window";
      }
      const auto other=std::lower_bound(opposite.begin(),opposite.end(),b,[](const auto& row,auto p){return row.end<=p;});
      if(other!=opposite.end() && other->begin<e)return "opposite_exact_evidence";
      if(query.substr(b,e-b).find('N')!=std::string_view::npos)return "contains_N";
      return {};
    };
    vertices.clear();edges.clear();vertices.push_back({begin,begin});
    auto root=std::lower_bound(own.begin(),own.end(),begin,[](const auto& row,auto p){return row.end<=p;});
    for(;root!=own.end() && root->begin<end;++root) {
      check();const auto b=std::max(begin,root->begin),e=std::min(end,root->end);
      const auto reason=barrier(b,e);
      // Never manufacture shorter, newly eligible pivots by removing a vetoed
      // base from a root segment. The normalized query segment is all-or-nothing.
      if(!reason.empty()){decision(b,e,false,reason);continue;}
      vertices.push_back({b,e});
    }
    if(vertices.size()==1){decision(begin,end,false,"no_valid_internal_root");continue;}
    vertices.push_back({end,end});
    for(std::size_t j=1;j<vertices.size();++j) {
      check();const auto b=vertices[j-1].end,e=vertices[j].begin;
      const auto reason=e-b>limit?std::string_view{"gap_too_long"}:barrier(b,e);
      edges.push_back({reason.empty(),reason});
    }
    // A linear graph needs only two reachability frontiers, not path copies or
    // recursive rescans. An isolated middle component cannot claim direction.
    std::size_t l=0,r=vertices.size()-1;
    while(l<edges.size() && edges[l].pass){check();++l;}
    while(r && edges[r-1].pass){check();--r;}
    auto emit_unoccupied=[&](Position b,Position e,bool inferred) {
      auto used=std::lower_bound(occupied.begin(),occupied.end(),b,[](const auto& row,auto p){return row.end<=p;});
      Position cursor=b;
      for(;used!=occupied.end() && used->begin<e;++used) {
        check();
        if(cursor<used->begin)emit_piece(CentromereSparseOrientationInterval{
          0,cursor,std::min(e,used->begin),begin,end,left.strand,inferred});
        cursor=std::max(cursor,used->end);
        if(cursor>=e)return;
      }
      if(cursor<e)emit_piece(CentromereSparseOrientationInterval{0,cursor,e,begin,end,left.strand,inferred});
    };
    for(std::size_t j=1;j<vertices.size();++j) {
      check();const auto b=vertices[j-1].end,e=vertices[j].begin;
      const bool connected=edges[j-1].pass && (j<=l || j-1>=r);
      if(b<e)decision(b,e,connected,connected?"sparse_root_bridge":
        edges[j-1].pass?"isolated_internal_component":edges[j-1].reason);
      // Test the original full edge before subtracting old predictions. Occupied
      // inferred ends cannot turn an over-limit edge into an admissible one.
      if(connected)emit_unoccupied(b,e,true);
      if(j+1<vertices.size() && (j<=l || j>=r))emit_unoccupied(vertices[j].begin,vertices[j].end,false);
    }
  }
}
} // namespace ramag::detail


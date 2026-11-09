#pragma once
#include "ramag/centromere_blocks.hpp"
#include <algorithm>
#include <array>
#include <memory_resource>
#include <span>
#include <vector>

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

#pragma once

#include "ramag/centromere_blocks.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <memory_resource>
#include <numeric>
#include <set>
#include <tuple>

namespace ramag::detail {
// Stable source identities survive node sorting. A cropped end is admissible
// only beside its actual geometric witness; a witness may be cropped at the
// far end, but not at the shared connection. This adds no path-state dimension.
struct AnchorClipWitnesses {
  std::size_t source{std::numeric_limits<std::size_t>::max()};
  std::size_t front{std::numeric_limits<std::size_t>::max()},back{std::numeric_limits<std::size_t>::max()};
  Length trim_front{},trim_back{};
};
inline bool AnchorClipConnectionAllowed(const AnchorClipWitnesses& left,const AnchorClipWitnesses& right) {
  constexpr auto absent=std::numeric_limits<std::size_t>::max();
  return (left.back==absent || (left.back==right.source && !right.trim_front)) &&
    (right.front==absent || (right.front==left.source && !left.trim_back));
}
// End clipping is a conservative local repair. A valid geometric witness is
// insufficient when the remaining unsupported gap is at least as large as
// either side's rarity-weighted matching support. Exact adjacency remains
// legal even for zero weight; original-to-original edges do not use this gate.
inline bool AnchorClipGapHasLocalSupport(Length reference_gap,Length query_gap,
    double left_weight,double right_weight) {
  return (!reference_gap && !query_gap) ||
    static_cast<long double>(std::max(reference_gap,query_gap))<
      static_cast<long double>(std::min(left_weight,right_weight));
}
struct AnchorEndClip {
  Length front{},back{};
  std::size_t front_witness{std::numeric_limits<std::size_t>::max()},back_witness{std::numeric_limits<std::size_t>::max()};
  double front_cost{},back_cost{};
};
struct AnchorEndClipSearch {
  std::pmr::vector<AnchorEndClip> clips;
  std::uint64_t pairs_examined{},local_support_rejections{};
  bool budget_exhausted{};
  explicit AnchorEndClipSearch(std::pmr::memory_resource* memory):clips(memory){}
};

// Generate a linear-sized set of end-clipping alternatives.  Each end retains
// the least clipping required by an ordered, partially overlapping neighbour;
// their combination also permits an anchor between two local indels. This is
// a bounded candidate generator, not an exhaustive optimization over cuts.
// Original anchors remain available to the planner. Clipped alternatives must
// never be counted as additional independent orientation witnesses.
template<class GetAnchor,class Check>
AnchorEndClipSearch FindAnchorEndClips(std::size_t count,GetAnchor get,
    Strand strand,Length query_length,std::pmr::memory_resource* memory,Check check) {
  struct Interval {Position r,u;Length n;double weight;};
  AnchorEndClipSearch result(memory);
  result.clips.resize(count);
  std::pmr::vector<Interval> intervals(memory);
  intervals.reserve(count);
  for(std::size_t i=0;i<count;++i) {
    const auto& a=get(i);
    intervals.push_back({a.reference_begin,
      strand==Strand::Forward ? a.query_begin : query_length-a.query_begin-a.length,a.length,a.weight});
  }
  std::pmr::vector<std::size_t> order(count,0,memory);
  std::iota(order.begin(),order.end(),0);
  // Repetitive input can have quadratically many overlapping rectangles. Stop
  // after a linear amount of geometric work, and reject this entire proposal
  // set on exhaustion rather than silently accepting a prefix of alternatives.
  constexpr std::uint64_t comparisons_per_anchor=64;
  const auto limit=count>std::numeric_limits<std::uint64_t>::max()/comparisons_per_anchor ?
    std::numeric_limits<std::uint64_t>::max() : count*comparisons_per_anchor;
  auto consider=[&](std::size_t ai,std::size_t bi) {
    if(result.pairs_examined==limit){result.budget_exhausted=true;return;}
    ++result.pairs_examined;
    if((result.pairs_examined&1023U)==0)check();
    if(intervals[bi].r<intervals[ai].r)std::swap(ai,bi);
    const auto& a=intervals[ai];const auto& b=intervals[bi];
    if(a.r>=b.r || a.u>=b.u || a.r+a.n>=b.r+b.n || a.u+a.n>=b.u+b.n)return;
    const Length ro=a.r+a.n>b.r ? a.r+a.n-b.r : 0;
    const Length qo=a.u+a.n>b.u ? a.u+a.n-b.u : 0;
    const auto trim=std::max(ro,qo),shorter=std::min(a.n,b.n);
    // Keep local end conflicts distinct from competing repeat-copy mappings:
    // each anchor must retain more sequence than the overlap being removed.
    if(!trim || trim>=shorter || trim>=shorter-trim)return;
    const auto rg=b.r-a.r,qg=b.u-a.u;
    const auto delta=rg>=qg?rg-qg:qg-rg;
    const auto cost=2.0*std::log2(static_cast<double>(delta)+1.0);
    // Geometric compatibility alone is insufficient. The cropped piece must
    // have an actual coherent connection under the planner's existing score;
    // otherwise an unsupported short fragment could enter the macro plan.
    auto keep=[&](Length& old,std::size_t& witness,double& old_cost,
        std::size_t partner,double own_weight,double partner_weight) {
      // Compute the actual unsupported interval after clipping either end.
      // Both alternatives have these gap lengths, but different retained
      // weights, so front and back acceptance remain separate decisions.
      const auto reference_gap=b.r-(a.r+a.n-trim);
      const auto query_gap=b.u-(a.u+a.n-trim);
      if(!AnchorClipGapHasLocalSupport(reference_gap,query_gap,own_weight,partner_weight)) {
        ++result.local_support_rejections;return;
      }
      if(delta && cost>=std::min(own_weight,partner_weight))return;
      if(!old || trim<old || (trim==old && (cost<old_cost ||
         (cost==old_cost && partner<witness)))) {
        old=trim;witness=partner;old_cost=cost;
      }
    };
    const auto aw=a.weight*static_cast<double>(a.n-trim)/static_cast<double>(a.n);
    const auto bw=b.weight*static_cast<double>(b.n-trim)/static_cast<double>(b.n);
    keep(result.clips[ai].back,result.clips[ai].back_witness,result.clips[ai].back_cost,bi,aw,b.weight);
    keep(result.clips[bi].front,result.clips[bi].front_witness,result.clips[bi].front_cost,ai,bw,a.weight);
  };
  for(const bool query_axis:{false,true}) {
    auto begin=[&](auto i){return query_axis ? intervals[i].u : intervals[i].r;};
    std::sort(order.begin(),order.end(),[&](auto a,auto b){
      return std::tuple{begin(a),intervals[a].n,a}<std::tuple{begin(b),intervals[b].n,b};});
    std::pmr::set<std::pair<Position,std::size_t>> active(memory);
    for(const auto i:order) {
      check();const auto start=begin(i);
      while(!active.empty() && active.begin()->first<=start)active.erase(active.begin());
      for(const auto& [end,j]:active) {
        (void)end;
        consider(j,i);
        if(result.budget_exhausted){result.clips.clear();return result;}
      }
      active.emplace(start+intervals[i].n,i);
    }
  }
  return result;
}
} // namespace ramag::detail

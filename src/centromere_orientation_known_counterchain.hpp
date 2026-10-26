#pragma once
#include "ramag/centromere_blocks.hpp"
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <limits>
#include <tuple>

namespace ramag::detail {
// A final, query-only subtraction. No search callback, inferred pivot or
// re-planning path is exposed: every decision reads the frozen upstream result.
inline void ReviewKnownBlockCounterchains(std::string_view query,
    CentromereBlockResult& out,const CentromereBlockOptions& options) {
  if(!options.review_known_block_counterchains)return;
  if(!options.review_gap_interior_opposition || !options.review_unified_gaps ||
     !options.orientation_audit_known_blocks || !options.count_cluster_gap ||
     !options.count_inferred_gap || options.count_inferred_gap>std::numeric_limits<Length>::max()/2)
    throw AlignmentError("known counterchain review requires the frozen gap-interior stack");
  if(out.review_known_block_counterchains)throw AlignmentError("known counterchain review is single-pass");
  const auto started=std::chrono::steady_clock::now();
  out.review_known_block_counterchains=true;
  using Interval=CentromereQueryInterval;
  auto check=[&]{if(options.interruption_callback)options.interruption_callback("centromere-known-counterchain-review");};
  auto merge=[](auto& v) {
    std::sort(v.begin(),v.end(),[](const auto& a,const auto& b){return std::pair{a.begin,a.end}<std::pair{b.begin,b.end};});
    std::size_t n=0;for(const auto& x:v)if(x.begin<x.end) {
      if(n && x.begin<=v[n-1].end)v[n-1].end=std::max(v[n-1].end,x.end);else v[n++]=x;
    }v.resize(n);
  };
  auto subtract=[](Position begin,Position end,const auto& masks,auto emit) {
    auto it=std::lower_bound(masks.begin(),masks.end(),begin,[](const auto& a,auto b){return a.end<=b;});
    for(;it!=masks.end() && it->begin<end;++it) {
      if(begin<it->begin)emit(Interval{begin,std::min(end,it->begin)});
      begin=std::max(begin,it->end);if(begin>=end)return;
    }
    if(begin<end)emit(Interval{begin,end});
  };
  std::array<std::vector<Interval>,2> protected_spans;
  std::vector<Interval> labeled,withdrawn;
  std::array<std::vector<std::size_t>,2> roots;
  std::array<std::vector<Position>,2> prefix;
  struct Gap {std::uint64_t block;Position begin,end;Strand strand;};
  std::vector<Gap> gaps;
  std::vector<std::size_t> opposite,own;
  auto memory_check=[&] {
    std::uint64_t bytes=gaps.capacity()*sizeof(Gap)+(labeled.capacity()+withdrawn.capacity())*sizeof(Interval)+
      (opposite.capacity()+own.capacity())*sizeof(std::size_t)+
      out.known_counterchain_reviews.capacity()*(sizeof(CentromereKnownCounterchainReview)+32)+
      out.known_counterchain_witnesses.capacity()*sizeof(CentromereCollinearWitness)+
      out.known_counterchain_intervals.capacity()*sizeof(CentromereOrientationKnownInterval);
    for(unsigned s=0;s<2;++s)bytes+=protected_spans[s].capacity()*sizeof(Interval)+
      roots[s].capacity()*sizeof(std::size_t)+prefix[s].capacity()*sizeof(Position);
    if(bytes>options.planner_memory_limit_bytes)throw AlignmentError("known counterchain capacity budget exceeded; no truncation");
    out.known_counterchain_capacity_bytes=std::max(out.known_counterchain_capacity_bytes,bytes);
  };
  auto protect=[&](Strand strand,Position begin,Position end) {
    if(begin>end || end>query.size())throw AlignmentError("invalid known counterchain protected anchor");
    const auto limit=options.count_inferred_gap;
    protected_spans[strand==Strand::Forward?0U:1U].push_back({begin-std::min(begin,limit),
      end+std::min<Length>(query.size()-end,limit)});memory_check();
  };
  for(const auto& block:out.blocks)if(block.direction_known)
    for(const auto& a:block.planning_anchors){
      check();if(a.query_begin>query.size() || a.length>query.size()-a.query_begin)
        throw AlignmentError("invalid known counterchain planning anchor");
      protect(block.strand,a.query_begin,a.query_begin+a.length);
    }
  // These are frozen real anchors, including accepted sparse roots. The old
  // 20 kb reach is measured from the original match, never an inferred edge.
  for(const auto& p:out.orientation_unified_pivots){check();protect(p.strand,p.query_begin,p.query_end);}
  for(auto& v:protected_spans)merge(v);
  for(const auto& w:out.orientation_windows)if(w.label!='.')labeled.push_back({w.begin,w.end});
  merge(labeled);
  for(std::size_t i=0;i<out.unified_nodes.size();++i) {
    const auto& n=out.unified_nodes[i];const auto& a=n.source.anchor;
    if(n.search_id || n.neutral || !n.root_id)continue;
    if(!a.length || a.query_begin>query.size() || a.length>query.size()-a.query_begin ||
       a.reference_begin>std::numeric_limits<Position>::max()-a.length ||
       (n.source.strand!=Strand::Forward && n.source.strand!=Strand::Reverse) ||
       !std::isfinite(a.weight) || a.weight<0)throw AlignmentError("invalid frozen counterchain root");
    roots[n.source.strand==Strand::Forward?0U:1U].push_back(i);
  }
  for(unsigned s=0;s<2;++s) {
    auto& ids=roots[s];std::sort(ids.begin(),ids.end(),[&](auto a,auto b){
      const auto& x=out.unified_nodes[a].source.anchor;const auto& y=out.unified_nodes[b].source.anchor;
      return std::tuple{x.query_begin,x.reference_begin,x.length,a}<std::tuple{y.query_begin,y.reference_begin,y.length,b};});
    Position end=0;for(auto i:ids){const auto& a=out.unified_nodes[i].source.anchor;end=std::max(end,a.query_begin+a.length);prefix[s].push_back(end);}
  }
  // Only long internal gaps already challenged by a pure opposite window enter
  // this extra audit. New withdrawals cannot create another eligible gap.
  for(const auto& a:out.orientation_block_audits)
    if(a.gap_begin<a.gap_end && a.gap_end-a.gap_begin>options.count_inferred_gap)
      gaps.push_back({a.block_id,a.gap_begin,a.gap_end,a.strand});
  auto key=[](const auto& g){return std::tuple{g.block,g.begin,g.end,g.strand};};
  std::sort(gaps.begin(),gaps.end(),[&](const auto& a,const auto& b){return key(a)<key(b);});
  gaps.erase(std::unique(gaps.begin(),gaps.end(),[&](const auto& a,const auto& b){return key(a)==key(b);}),gaps.end());
  memory_check();
  auto visit=[&](unsigned sign,Position begin,Position end,auto emit) {
    auto k=static_cast<std::size_t>(std::upper_bound(prefix[sign].begin(),prefix[sign].end(),begin)-prefix[sign].begin());
    for(;k<roots[sign].size();++k) {
      const auto id=roots[sign][k];const auto& a=out.unified_nodes[id].source.anchor;
      if(a.query_begin>=end)break;
      if(a.query_begin+a.length>begin)emit(id);
    }
  };
  // The transition is the existing two-axis collinearity/weight condition.
  // Opposite proposals have the stricter old local-neighbor bound on both axes;
  // retained own chains keep the existing original-anchor 2*20 kb query reach.
  auto edge=[&](const auto& a,const auto& b,unsigned sign,Length bound,bool bound_reference,double& penalty) {
    if(a.query_begin+a.length>b.query_begin)return false;
    const auto qgap=b.query_begin-a.query_begin-a.length;if(qgap>bound)return false;
    Length rgap=0;
    if(!sign){if(a.reference_begin+a.length>b.reference_begin)return false;rgap=b.reference_begin-a.reference_begin-a.length;}
    else {if(b.reference_begin+b.length>a.reference_begin)return false;rgap=a.reference_begin-b.reference_begin-b.length;}
    if(bound_reference && rgap>bound)return false;
    const auto delta=qgap>rgap?qgap-rgap:rgap-qgap;
    penalty=2*std::log2(static_cast<double>(delta)+1);
    return !delta || penalty<std::min(a.weight,b.weight);
  };
  for(const auto& gap:gaps)for(const auto& retained:out.orientation_known_intervals) {
    check();if(retained.block_id!=gap.block || retained.strand!=gap.strand)continue;
    const auto begin=std::max(gap.begin,retained.begin),end=std::min(gap.end,retained.end);
    if(begin>=end)continue;
    const unsigned sign=gap.strand==Strand::Forward?0U:1U;
    subtract(begin,end,labeled,[&](Interval uncovered) {
      subtract(uncovered.begin,uncovered.end,protected_spans[sign],[&](Interval cell) {
        opposite.clear();visit(1-sign,cell.begin,cell.end,[&](auto i) {
          const auto& a=out.unified_nodes[i].source.anchor;
          if(a.query_begin>=cell.begin && a.query_begin+a.length<=cell.end)opposite.push_back(i);
        });memory_check();
        for(std::size_t right=1;right<opposite.size();++right)for(std::size_t left=0;left<right;++left) {
          check();const auto li=opposite[left],ri=opposite[right];
          const auto& a=out.unified_nodes[li].source.anchor;const auto& b=out.unified_nodes[ri].source.anchor;
          ++out.known_counterchain_connection_checks;double penalty=0;
          if(!edge(a,b,1-sign,options.count_cluster_gap,true,penalty))continue;
          CentromereKnownCounterchainReview row;row.id=out.known_counterchain_reviews.size()+1;
          row.block_id=gap.block;row.begin=a.query_begin;row.end=b.query_begin+b.length;
          row.gap_begin=gap.begin;row.gap_end=gap.end;row.strand=gap.strand;
          row.left_node=out.unified_nodes[li].id;row.right_node=out.unified_nodes[ri].id;
          row.opposite_support=a.length+b.length;row.opposite_score=a.weight+b.weight-penalty;row.opposite_nodes=2;
          row.reason="contains_N";
          if(query.substr(row.begin,row.end-row.begin).find_first_not_of("ACGT")==std::string_view::npos) {
            // A valid old-direction edge may span the whole proposal while
            // neither endpoint intersects it. Include its bounded original
            // flanks; independent neighboring chains remain outside the audit.
            const auto reach=2*options.count_inferred_gap;
            auto low=row.begin-std::min(row.begin,reach);
            auto high=row.end+std::min<Length>(query.size()-row.end,reach);
            if(low)--low;
            if(high<query.size())++high;
            own.clear();visit(sign,low,high,[&](auto i){own.push_back(i);});memory_check();
            bool chain=false;
            // Length and weight need not be monotone. The counterchain must
            // strictly beat every intact local root on both axes, not merely
            // the lexicographically selected root. These are separate maxima.
            for(auto i:own){const auto& x=out.unified_nodes[i].source.anchor;
              if(x.query_begin>=row.end || x.query_begin+x.length<=row.begin)continue;
              row.own_support=std::max(row.own_support,x.length);
              row.own_score=std::max(row.own_score,x.weight);row.own_nodes=1;
            }
            Length pair_support=0;double pair_score=0;
            for(std::size_t j=1;j<own.size();++j)for(std::size_t i=0;i<j;++i) {
              check();const auto& x=out.unified_nodes[own[i]].source.anchor;const auto& y=out.unified_nodes[own[j]].source.anchor;
              if(x.query_begin>=row.end || y.query_begin+y.length<=row.begin)continue;
              ++row.connection_checks;++out.known_counterchain_connection_checks;double cost=0;
              if(!edge(x,y,sign,reach,false,cost))continue;
              chain=true;const auto support=x.length+y.length;const auto score=x.weight+y.weight-cost;
              if(support>pair_support || (support==pair_support && score>pair_score)){pair_support=support;pair_score=score;}
            }
            if(chain){row.own_support=pair_support;row.own_score=pair_score;row.own_nodes=2;row.reason="own_bounded_chain";}
            else {
              row.reason="own_support_not_weaker";
              if(row.opposite_support>row.own_support && row.opposite_score>row.own_score) {
                row.retracted=true;row.reason="counterchain_retraction";withdrawn.push_back({row.begin,row.end});
              }
            }
          }
          out.known_counterchain_reviews.push_back(row);
          for(auto [i,ordinal]:{std::pair{li,1U},std::pair{ri,2U}}) {
            const auto& n=out.unified_nodes[i];const auto& w=n.source.anchor;
            out.known_counterchain_witnesses.push_back({row.id,ordinal,n.id,w.reference_begin,w.query_begin,w.length,w.weight});
          }
          memory_check();
        }
      });
    });
  }
  merge(withdrawn);
  for(const auto& old:out.orientation_known_intervals)
    subtract(old.begin,old.end,withdrawn,[&](Interval x){out.known_counterchain_intervals.push_back({old.block_id,x.begin,x.end,old.strand});});
  for(const auto& x:withdrawn)out.orientation_unknown.push_back(x);
  merge(out.orientation_unknown);memory_check();
  out.known_counterchain_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
}
} // namespace ramag::detail

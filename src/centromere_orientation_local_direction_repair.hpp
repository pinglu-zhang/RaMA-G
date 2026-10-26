#pragma once
#include "ramag/centromere_blocks.hpp"
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <limits>
#include <tuple>

namespace ramag::detail {
// One final query-direction view. Candidate roots, challenged gaps, protections
// and unknown space are frozen before any decision; repairs never feed back.
inline void RepairInheritedDirections(std::string_view query,
    CentromereBlockResult& out,const CentromereBlockOptions& options) {
  if(options.repair_drift_consistent_gaps && !options.repair_inherited_directions)
    throw AlignmentError("drift-consistent repair requires inherited direction repair");
  if(!options.repair_inherited_directions)return;
  if(!options.review_known_block_counterchains || !out.review_known_block_counterchains ||
     !options.review_unified_gaps || !options.review_gap_interior_opposition ||
     !options.orientation_audit_known_blocks || !options.count_inferred_gap ||
     options.count_inferred_gap>std::numeric_limits<Length>::max()/2)
    throw AlignmentError("local direction repair requires the frozen known-counterchain stack");
  if(out.repair_inherited_directions)throw AlignmentError("local direction repair is single-pass");
  const auto started=std::chrono::steady_clock::now();
  const auto reach=options.count_inferred_gap,bound=2*reach;
  const bool drift=options.repair_drift_consistent_gaps;
  const Length query_length=query.size();
  auto check=[&]{if(options.interruption_callback)options.interruption_callback("centromere-local-direction-repair");};
  check();out.repair_inherited_directions=true;out.repair_drift_consistent_gaps=drift;
  using Interval=CentromereQueryInterval;
  auto merge=[](auto& v) {
    std::sort(v.begin(),v.end(),[](const auto& a,const auto& b){return std::pair{a.begin,a.end}<std::pair{b.begin,b.end};});
    std::size_t n=0;for(const auto& x:v)if(x.begin<x.end) {
      if(n && x.begin<=v[n-1].end)v[n-1].end=std::max(v[n-1].end,x.end);else v[n++]=x;
    }v.resize(n);
  };
  auto overlap=[](const auto& v,Position begin,Position end) {
    const auto it=std::lower_bound(v.begin(),v.end(),begin,[](const auto& x,auto p){return x.end<=p;});
    return it!=v.end() && it->begin<end;
  };
  auto subtract=[](Position begin,Position end,const auto& masks,auto emit) {
    auto it=std::lower_bound(masks.begin(),masks.end(),begin,[](const auto& x,auto p){return x.end<=p;});
    for(;it!=masks.end() && it->begin<end;++it) {
      if(begin<it->begin)emit(Interval{begin,std::min(end,it->begin)});
      begin=std::max(begin,it->end);if(begin>=end)return;
    }
    if(begin<end)emit(Interval{begin,end});
  };
  auto intersect=[](Position begin,Position end,const auto& masks,auto emit) {
    auto it=std::lower_bound(masks.begin(),masks.end(),begin,[](const auto& x,auto p){return x.end<=p;});
    for(;it!=masks.end() && it->begin<end;++it)emit(Interval{std::max(begin,it->begin),std::min(end,it->end)});
  };
  auto sign_index=[](Strand strand)->unsigned {
    if(strand!=Strand::Forward && strand!=Strand::Reverse)throw AlignmentError("invalid local repair strand");
    return strand==Strand::Forward?0U:1U;
  };
  auto strand_of=[](unsigned sign){return sign?Strand::Reverse:Strand::Forward;};
  auto sign_char=[](unsigned sign){return sign?'-':'+';};
  struct Gap {std::uint64_t block;Position begin,end;Strand strand;};
  std::vector<Gap> gaps;
  std::array<std::vector<Interval>,2> known,protected_spans,assigned,claimed,final;
  std::array<std::vector<std::size_t>,2> roots;
  std::array<std::vector<Position>,2> prefix;
  std::vector<Interval> frozen_unknown,allowed,residual,all_final,read_allowed,ambiguous;
  struct WritePiece {std::uint64_t review;Position begin,end;};
  std::vector<WritePiece> written;
  std::vector<std::size_t> opposite,own;
  auto memory_check=[&] {
    std::uint64_t bytes=gaps.capacity()*sizeof(Gap)+(opposite.capacity()+own.capacity())*sizeof(std::size_t)+
      (frozen_unknown.capacity()+allowed.capacity()+residual.capacity()+all_final.capacity()+
        read_allowed.capacity()+ambiguous.capacity())*sizeof(Interval)+written.capacity()*sizeof(WritePiece)+
      out.local_direction_repair_reviews.capacity()*(sizeof(CentromereLocalDirectionRepairReview)+64)+
      out.local_direction_repair_witnesses.capacity()*sizeof(CentromereCollinearWitness)+
      out.local_direction_repair_changes.capacity()*(sizeof(CentromereLocalDirectionRepairChange)+64)+
      out.local_direction_repair_intervals.capacity()*sizeof(CentromereOrientationKnownInterval);
    for(unsigned s=0;s<2;++s)bytes+=(known[s].capacity()+protected_spans[s].capacity()+assigned[s].capacity()+
      claimed[s].capacity()+final[s].capacity())*sizeof(Interval)+
      roots[s].capacity()*sizeof(std::size_t)+prefix[s].capacity()*sizeof(Position);
    if(bytes>options.planner_memory_limit_bytes)throw AlignmentError("local direction repair capacity budget exceeded; no truncation");
    out.local_direction_repair_capacity_bytes=std::max(out.local_direction_repair_capacity_bytes,bytes);
  };
  // Long candidates reuse one immutable ambiguity index, rather than rescan
  // their possibly overlapping spans. No sequence search or new root is made.
  if(drift)for(Position begin=0;begin<query_length;) {
    if(begin%8192==0)check();
    const auto c=query[begin];
    if(c=='A' || c=='C' || c=='G' || c=='T'){++begin;continue;}
    const auto first=begin++;
    while(begin<query_length) {
      if(begin%8192==0)check();
      const auto d=query[begin];if(d=='A' || d=='C' || d=='G' || d=='T')break;
      ++begin;
    }
    ambiguous.push_back({first,begin});memory_check();
  }
  auto valid_interval=[&](Position begin,Position end) {
    if(begin>end || end>query_length)throw AlignmentError("invalid local direction repair interval");
  };
  auto add_known=[&](const auto& row) {
    valid_interval(row.begin,row.end);known[sign_index(row.strand)].push_back({row.begin,row.end});memory_check();
  };
  auto add=[&](const auto& rows){for(const auto& x:rows){check();add_known(x);}};
  add(out.known_counterchain_intervals);add(out.orientation_intervals);add(out.orientation_long_gap_intervals);
  add(out.orientation_sparse_intervals);add(out.orientation_short_gap_intervals);add(out.orientation_mixed_gap_intervals);
  add(out.orientation_sparse_gap_intervals);add(out.orientation_unified_intervals);add(out.region_intervals);
  add(out.unified_intervals);
  auto add_reviews=[&](const auto& rows){for(const auto& x:rows)if(x.accepted){check();add_known(x);}};
  add_reviews(out.orientation_local_reviews);add_reviews(out.orientation_isolated_reviews);
  add_reviews(out.orientation_collinear_reviews);add_reviews(out.orientation_alternate_reviews);
  for(auto& v:known)merge(v);
  for(const auto& x:known[0])if(overlap(known[1],x.begin,x.end))
    throw AlignmentError("local direction repair input directions overlap");
  frozen_unknown=out.orientation_unknown;memory_check();
  for(const auto& x:frozen_unknown){check();valid_interval(x.begin,x.end);}
  merge(frozen_unknown);
  for(const auto& x:frozen_unknown)if(overlap(known[0],x.begin,x.end) || overlap(known[1],x.begin,x.end))
    throw AlignmentError("local direction repair input unknown overlaps known direction");
  // Validate the frozen ledger, rather than silently filling an upstream hole.
  all_final=frozen_unknown;memory_check();
  for(const auto& v:known)for(const auto& x:v){all_final.push_back(x);memory_check();}
  merge(all_final);
  if(query_length && (all_final.size()!=1 || all_final[0].begin || all_final[0].end!=query_length))
    throw AlignmentError("local direction repair input ledger is incomplete");
  all_final.clear();
  auto protect=[&](Strand strand,Position begin,Position end) {
    valid_interval(begin,end);
    protected_spans[sign_index(strand)].push_back({begin-std::min(begin,reach),end+std::min(query_length-end,reach)});
    memory_check();
  };
  for(const auto& block:out.blocks)if(block.direction_known)for(const auto& a:block.planning_anchors) {
    check();if(a.query_begin>query_length || a.length>query_length-a.query_begin)
      throw AlignmentError("invalid local repair planning anchor");
    protect(block.strand,a.query_begin,a.query_begin+a.length);
  }
  for(const auto& pivot:out.orientation_unified_pivots){check();protect(pivot.strand,pivot.query_begin,pivot.query_end);}
  for(auto& v:protected_spans)merge(v);
  for(std::size_t i=0;i<out.unified_nodes.size();++i) {
    check();const auto& n=out.unified_nodes[i];const auto& a=n.source.anchor;
    if(n.search_id || n.neutral || !n.root_id)continue;
    if(!n.id || !a.length || a.query_begin>query_length || a.length>query_length-a.query_begin ||
       a.reference_begin>std::numeric_limits<Position>::max()-a.length ||
       !std::isfinite(a.weight) || a.weight<0)throw AlignmentError("invalid local repair original root");
    roots[sign_index(n.source.strand)].push_back(i);memory_check();
  }
  for(unsigned s=0;s<2;++s) {
    auto& ids=roots[s];std::sort(ids.begin(),ids.end(),[&](auto a,auto b){
      const auto& x=out.unified_nodes[a];const auto& y=out.unified_nodes[b];
      return std::tuple{x.source.anchor.query_begin,x.source.anchor.reference_begin,x.source.anchor.length,x.id}<
        std::tuple{y.source.anchor.query_begin,y.source.anchor.reference_begin,y.source.anchor.length,y.id};});
    Position end=0;for(auto i:ids){check();const auto& a=out.unified_nodes[i].source.anchor;
      end=std::max(end,a.query_begin+a.length);prefix[s].push_back(end);memory_check();}
  }
  auto visit=[&](unsigned sign,Position begin,Position end,auto emit) {
    auto at=static_cast<std::size_t>(std::upper_bound(prefix[sign].begin(),prefix[sign].end(),begin)-prefix[sign].begin());
    for(;at<roots[sign].size();++at) {
      check();const auto id=roots[sign][at];const auto& a=out.unified_nodes[id].source.anchor;
      if(a.query_begin>=end)break;
      if(a.query_begin+a.length>begin)emit(id);
    }
  };
  // V1 short edges and own-chain vetoes retain their original bounds. V2
  // additionally admits source-LCP 1/1 edges in the same 2*reach diagonal band:
  // long nearly equal gaps may agree, while inconsistent long jumps do not.
  auto edge=[&](const auto& a,const auto& b,unsigned sign,bool bound_reference,double& cost,bool allow_drift=false) {
    if(a.query_begin+a.length>b.query_begin)return false;
    const auto qgap=b.query_begin-a.query_begin-a.length;if(qgap>bound && !allow_drift)return false;
    Length rgap=0;
    if(!sign){if(a.reference_begin+a.length>b.reference_begin)return false;rgap=b.reference_begin-a.reference_begin-a.length;}
    else {if(b.reference_begin+b.length>a.reference_begin)return false;rgap=a.reference_begin-b.reference_begin-b.length;}
    if(bound_reference && rgap>bound && !allow_drift)return false;
    const auto delta=qgap>rgap?qgap-rgap:rgap-qgap;
    if(allow_drift && (qgap>bound || rgap>bound) && delta>bound)return false;
    cost=2*std::log2(static_cast<double>(delta)+1);
    return !delta || cost<std::min(a.weight,b.weight);
  };
  for(const auto& a:out.orientation_block_audits) {
    check();valid_interval(a.gap_begin,a.gap_end);sign_index(a.strand);
    if(a.gap_end-a.gap_begin>reach){gaps.push_back({a.block_id,a.gap_begin,a.gap_end,a.strand});memory_check();}
  }
  auto gap_key=[](const auto& g){return std::tuple{g.block,g.begin,g.end,g.strand};};
  std::sort(gaps.begin(),gaps.end(),[&](const auto& a,const auto& b){return gap_key(a)<gap_key(b);});
  gaps.erase(std::unique(gaps.begin(),gaps.end(),[&](const auto& a,const auto& b){return gap_key(a)==gap_key(b);}),gaps.end());
  for(const auto& gap:gaps) {
    check();const auto oldsign=sign_index(gap.strand),sign=1-oldsign;
    // Independent overlays are deliberately absent. Only this old block's
    // retained pieces and the frozen unknown ledger may change direction.
    allowed=frozen_unknown;memory_check();
    for(const auto& x:out.known_counterchain_intervals)if(x.block_id==gap.block && x.strand==gap.strand) {
      allowed.push_back({x.begin,x.end});memory_check();
    }
    merge(allowed);
    // Already frozen target-direction overlays may supply original evidence,
    // but are never writable or promoted to roots. Old-direction overlays
    // remain absent from the read corridor, so a candidate cannot cross them.
    if(drift) {
      read_allowed=allowed;memory_check();
      for(const auto& x:known[sign]){read_allowed.push_back(x);memory_check();}
      merge(read_allowed);
    }
    intersect(gap.begin,gap.end,drift?read_allowed:allowed,[&](Interval permitted) {
      subtract(permitted.begin,permitted.end,protected_spans[oldsign],[&](Interval cell) {
        opposite.clear();visit(sign,cell.begin,cell.end,[&](auto i) {
          const auto& a=out.unified_nodes[i].source.anchor;
          if(a.query_begin>=cell.begin && a.query_begin+a.length<=cell.end){opposite.push_back(i);memory_check();}
        });
        for(std::size_t left=0;left<opposite.size();++left)for(std::size_t right=left+1;right<opposite.size();++right) {
          check();const auto li=opposite[left],ri=opposite[right];
          const auto& a=out.unified_nodes[li].source.anchor;const auto& b=out.unified_nodes[ri].source.anchor;
          if(!drift && b.query_begin>a.query_begin+a.length && b.query_begin-a.query_begin-a.length>bound)break;
          const auto& left_source=out.unified_nodes[li].source;const auto& right_source=out.unified_nodes[ri].source;
          const bool source_singleton=left_source.reference_count==1 && left_source.query_count==1 &&
            right_source.reference_count==1 && right_source.query_count==1;
          ++out.local_direction_repair_connection_checks;double cost=0;
          if(!edge(a,b,sign,true,cost,drift && source_singleton))continue;
          const auto qgap=b.query_begin-a.query_begin-a.length;
          const auto rgap=!sign?b.reference_begin-a.reference_begin-a.length:a.reference_begin-b.reference_begin-b.length;
          const bool long_edge=qgap>bound || rgap>bound;
          CentromereLocalDirectionRepairReview row;row.id=out.local_direction_repair_reviews.size()+1;
          row.block_id=gap.block;row.begin=a.query_begin;row.end=b.query_begin+b.length;
          row.strand=strand_of(sign);row.inherited_strand=gap.strand;
          row.support=a.length+b.length;row.score=a.weight+b.weight-cost;row.nodes=2;row.reason="contains_N";
          if(!std::isfinite(row.score))throw AlignmentError("nonfinite local direction repair candidate score");
          if(drift?!overlap(ambiguous,row.begin,row.end):
              query.substr(row.begin,row.end-row.begin).find_first_not_of("ACGT")==std::string_view::npos) {
            auto low=row.begin-std::min(row.begin,bound),high=row.end+std::min(query_length-row.end,bound);
            if(low)--low;
            if(high<query_length)++high;
            own.clear();visit(oldsign,low,high,[&](auto i){own.push_back(i);memory_check();});
            for(auto i:own){const auto& x=out.unified_nodes[i].source.anchor;
              if(x.query_begin>=row.end || x.query_begin+x.length<=row.begin)continue;
              row.opposite_support=std::max(row.opposite_support,x.length);
              row.opposite_score=std::max(row.opposite_score,x.weight);row.opposite_nodes=1;
            }
            bool chain=false;
            for(std::size_t i=0;i<own.size() && !chain;++i) {
              const auto& x=out.unified_nodes[own[i]].source.anchor;if(x.query_begin>=row.end)break;
              for(std::size_t j=i+1;j<own.size();++j) {
                check();const auto& y=out.unified_nodes[own[j]].source.anchor;
                if(y.query_begin>x.query_begin+x.length && y.query_begin-x.query_begin-x.length>bound)break;
                if(y.query_begin+y.length<=row.begin)continue;
                ++out.local_direction_repair_connection_checks;double penalty=0;
                if(!edge(x,y,oldsign,false,penalty))continue;
                chain=true;row.opposite_support=x.length+y.length;
                row.opposite_score=x.weight+y.weight-penalty;
                if(!std::isfinite(row.opposite_score))throw AlignmentError("nonfinite local direction repair own-chain score");
                row.opposite_nodes=2;break;
              }
            }
            if(chain)row.reason="own_bounded_chain";
            else if(row.support<=row.opposite_support || row.score<=row.opposite_score)row.reason="own_support_not_weaker";
            else {
              row.accepted=true;row.reason=long_edge?"drift_consistent_root_direction_repair":"bounded_root_direction_repair";
              if(drift)intersect(row.begin,row.end,allowed,[&](Interval piece) {
                assigned[sign].push_back(piece);written.push_back({row.id,piece.begin,piece.end});memory_check();
              });
              else {assigned[sign].push_back({row.begin,row.end});memory_check();}
            }
          }
          out.local_direction_repair_reviews.push_back(row);memory_check();
          for(auto [i,ordinal]:{std::pair{li,1U},std::pair{ri,2U}}) {
            const auto& n=out.unified_nodes[i];const auto& x=n.source.anchor;
            out.local_direction_repair_witnesses.push_back({row.id,ordinal,n.id,x.reference_begin,x.query_begin,x.length,x.weight});
            memory_check();
          }
        }
      });
    });
  }
  for(auto& v:assigned)merge(v);
  for(const auto& x:assigned[0])if(overlap(assigned[1],x.begin,x.end))
    throw AlignmentError("local direction repair proposals conflict");
  // A change is attributed to the first accepted review that covers it, so
  // overlapping witnesses cannot duplicate bases in the public change ledger.
  std::size_t written_index=0;
  for(const auto& row:out.local_direction_repair_reviews)if(row.accepted) {
    check();const auto sign=sign_index(row.strand);
    auto claim=[&](Position begin,Position end) {
      residual.clear();
      subtract(begin,end,claimed[sign],[&](Interval x){residual.push_back(x);memory_check();});
      for(const auto& x:residual)subtract(x.begin,x.end,known[sign],[&](Interval changed) {
        auto emit=[&](Interval piece,char before) {
          out.local_direction_repair_changes.push_back({row.id,piece.begin,piece.end,before,sign_char(sign),row.reason});
          memory_check();
        };
        intersect(changed.begin,changed.end,known[1-sign],[&](Interval piece){emit(piece,sign_char(1-sign));});
        intersect(changed.begin,changed.end,frozen_unknown,[&](Interval piece){emit(piece,'.');});
      });
      claimed[sign].push_back({begin,end});memory_check();merge(claimed[sign]);
    };
    if(drift)while(written_index<written.size() && written[written_index].review==row.id) {
      const auto& x=written[written_index++];claim(x.begin,x.end);
    }
    else claim(row.begin,row.end);
  }
  std::sort(out.local_direction_repair_changes.begin(),out.local_direction_repair_changes.end(),[](const auto& a,const auto& b){
    return std::tuple{a.begin,a.end,a.review_id}<std::tuple{b.begin,b.end,b.review_id};});
  for(unsigned sign=0;sign<2;++sign) {
    for(const auto& x:known[sign]){check();subtract(x.begin,x.end,assigned[1-sign],[&](Interval piece){final[sign].push_back(piece);memory_check();});}
    for(const auto& x:assigned[sign]){final[sign].push_back(x);memory_check();}
    merge(final[sign]);
    for(const auto& x:final[sign]) {
      out.local_direction_repair_intervals.push_back({0,x.begin,x.end,strand_of(sign)});
      all_final.push_back(x);memory_check();
    }
  }
  for(const auto& x:final[0])if(overlap(final[1],x.begin,x.end))throw AlignmentError("local direction repair final directions overlap");
  std::sort(out.local_direction_repair_intervals.begin(),out.local_direction_repair_intervals.end(),[](const auto& a,const auto& b){
    return std::pair{a.begin,a.end}<std::pair{b.begin,b.end};});
  merge(all_final);residual.clear();
  subtract(0,query_length,all_final,[&](Interval x){residual.push_back(x);memory_check();});
  out.orientation_unknown=std::move(residual);
  out.local_direction_repair_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
}
} // namespace ramag::detail

#include "ramag/internal/centromere/gap_review.hpp"
#include "ramag/internal/centromere/evidence_review.hpp"
#include "ramag/internal/centromere/placement_review.hpp"
namespace ramag::detail {
void ReviewKnownBlockCounterchains(std::string_view query,
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

void RepairInheritedDirections(std::string_view query,
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

void ApplyDirectionIslandAudit(Length reference_length,Length query_length,
    CentromereBlockResult& out,const CentromereBlockOptions& options,
    std::string_view reference,std::string_view query) {
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

std::uint64_t LocalProtectionBlockCapacity(std::span<const CentromereBlock> blocks,std::size_t capacity) {
  std::uint64_t bytes=capacity*sizeof(CentromereBlock);
  for(const auto& b:blocks)bytes+=b.planning_anchors.capacity()*sizeof(CentromereAnchor)+b.status.capacity()+b.reason.capacity()+2;
  return bytes;
}

FrozenLocalDirectionProtection FreezeLocalDirectionProtectionBaseline(CentromereBlockResult& baseline) {
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

void ApplyLocalDirectionProtection(Length query_length,FrozenLocalDirectionProtection& baseline,
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
}

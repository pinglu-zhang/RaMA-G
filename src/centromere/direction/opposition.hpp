#pragma once
#include "gap.hpp"
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <limits>
#include <numeric>
#include <span>
#include <tuple>

// Related direction checks share a private module; their admission rules
// and template callback paths remain independent.

// orientation_block_audit
namespace ramag::detail {
// Old chain envelopes remain provenance, not unconditional orientation truth.
// Audit only their internal, anchor-free long gaps. A pure opposite window cell
// retracts the old direction to unknown; it does not assign a new reference copy
// or turn the window's envelope into a new opposite-direction prediction.
// No paths or sequences are copied. Existing chains are reference-ascending;
// reverse chains can therefore be visited backwards in query order.
template<class EmitAudit,class EmitKnown,class Check>
void VisitKnownOrientationBlockAudit(std::span<const CentromereBlock> blocks,
    std::span<const OrientationCoverage> labels,Length query_length,Length gap_limit,
    EmitAudit audit,EmitKnown known,Check check,bool audit_all_supported_conflicts=false) {
  if(!gap_limit)throw AlignmentError("known-block audit requires a positive gap limit");
  for(std::size_t i=0;i<labels.size();++i) {
    check();const auto& row=labels[i];
    if(row.begin>=row.end || row.end>query_length || row.mask>7U ||
       (i && labels[i-1].end>row.begin))
      throw AlignmentError("invalid known-block audit coverage cells");
  }
  for(const auto& block:blocks) {
    check();if(!block.direction_known)continue;
    if(block.query_begin>block.query_end || block.query_end>query_length ||
       block.reference_begin>block.reference_end)
      throw AlignmentError("invalid known-block audit bounds");
    const auto& anchors=block.planning_anchors;
    for(std::size_t i=0;i<anchors.size();++i) {
      check();const auto& a=anchors[i];
      // Subtraction bounds the additions below and rejects overflow as well.
      if(!a.length || a.query_begin<block.query_begin || a.query_begin>block.query_end ||
         a.length>block.query_end-a.query_begin || a.reference_begin<block.reference_begin ||
         a.reference_begin>block.reference_end || a.length>block.reference_end-a.reference_begin)
        throw AlignmentError("planning anchor lies outside its audited block");
      if(i) {
        const auto& p=anchors[i-1];
        if(p.reference_begin+p.length>a.reference_begin ||
           (block.strand==Strand::Forward ? p.query_begin+p.length>a.query_begin :
             a.query_begin+a.length>p.query_begin))
          throw AlignmentError("planning anchors are not an ordered nonoverlapping chain");
      }
    }
    if(block.query_begin==block.query_end)continue;
    Position retained_begin=block.query_begin;
    const unsigned opposite=block.strand==Strand::Forward?2U:1U;
    auto anchor=[&](std::size_t i)->const CentromereAnchor& {
      return anchors[block.strand==Strand::Forward?i:anchors.size()-1-i];
    };
    for(std::size_t i=1;i<anchors.size();++i) {
      check();const auto& left=anchor(i-1);const auto& right=anchor(i);
      const Position gap_begin=left.query_begin+left.length,gap_end=right.query_begin;
      // The opt-in repair admits direct opposing-window evidence at any
      // positive gap size. Planning anchors themselves remain unchanged.
      if(gap_begin==gap_end || (!audit_all_supported_conflicts && gap_end-gap_begin<=gap_limit))continue;
      auto label=std::lower_bound(labels.begin(),labels.end(),gap_begin,
        [](const auto& row,auto p){return row.end<=p;});
      Position cut_begin=0,cut_end=0;bool pending=false;
      auto emit=[&] {
        if(!pending)return;
        if(retained_begin<cut_begin)
          known(CentromereOrientationKnownInterval{block.id,retained_begin,cut_begin,block.strand});
        audit(CentromereOrientationBlockAudit{block.id,cut_begin,cut_end,gap_begin,gap_end,block.strand});
        retained_begin=cut_end;pending=false;
      };
      for(;label!=labels.end() && label->begin<gap_end;++label) {
        check();if(label->mask!=opposite)continue;
        const auto begin=std::max(gap_begin,label->begin),end=std::min(gap_end,label->end);
        if(begin>=end)continue;
        if(pending && begin==cut_end)cut_end=end;
        else{emit();cut_begin=begin;cut_end=end;pending=true;}
      }
      emit();
    }
    if(retained_begin<block.query_end)
      known(CentromereOrientationKnownInterval{block.id,retained_begin,block.query_end,block.strand});
  }
}
} // namespace ramag::detail

// orientation_isolated_review
namespace ramag::detail {
// This does not remove a root or alter its rarity, chain, window or mapping.
// Only a bounded query-only gap can be emitted, from a frozen input snapshot.
template<class Nodes,class Neutral,class Emit,class EmitCounter,class Check>
void VisitIsolatedOppositeGaps(const Nodes& nodes,const Neutral& neutral,
    std::span<const CentromereUnifiedPivot> pivots,std::span<const OrientationCoverage> labels,
    std::string_view query,Length limit,Length neighborhood,
    std::span<const CentromereQueryInterval> known_forward,
    std::span<const CentromereQueryInterval> known_reverse,
    std::pmr::memory_resource* memory,Emit emit,EmitCounter emit_counter,Check check) {
  if(!limit || !neighborhood || !memory || neutral.size()!=nodes.size())
    throw AlignmentError("invalid isolated-opposite review configuration");
  auto intersects=[](auto rows,Position b,Position e) {
    const auto it=std::lower_bound(rows.begin(),rows.end(),b,[](const auto& x,auto p){return x.end<=p;});
    return it!=rows.end() && it->begin<e;
  };
  for(auto rows:{known_forward,known_reverse})for(std::size_t i=0;i<rows.size();++i) {
    check();if(rows[i].begin>=rows[i].end || rows[i].end>query.size() || (i && rows[i-1].end>=rows[i].begin))
      throw AlignmentError("isolated review requires merged known intervals");
  }
  for(const auto& x:known_forward)if(intersects(known_reverse,x.begin,x.end))
    throw AlignmentError("overlapping known directions");
  for(std::size_t i=0;i<labels.size();++i)if(labels[i].begin>=labels[i].end || labels[i].end>query.size() ||
      labels[i].mask>7U || (i && labels[i-1].end>labels[i].begin))throw AlignmentError("invalid isolated review labels");
  using Ids=std::pmr::vector<std::size_t>;
  std::array<Ids,2> roots{Ids(memory),Ids(memory)};
  std::array<std::pmr::vector<Position>,2> prefix{std::pmr::vector<Position>(memory),std::pmr::vector<Position>(memory)};
  for(std::size_t i=0;i<nodes.size();++i) {
    check();const auto& n=nodes[i];const auto& a=n.anchor;
    if(!a.length || a.query_begin>query.size() || a.length>query.size()-a.query_begin ||
       a.reference_begin>std::numeric_limits<Position>::max()-a.length || !std::isfinite(a.weight) || a.weight<0 ||
       (n.strand!=Strand::Forward && n.strand!=Strand::Reverse))throw AlignmentError("invalid normalized review root");
    if(!neutral[i])roots[n.strand==Strand::Forward?0U:1U].push_back(i);
  }
  for(unsigned s=0;s<2;++s) {
    auto& ids=roots[s];std::sort(ids.begin(),ids.end(),[&](auto a,auto b) {
      return std::pair{nodes[a].anchor.query_begin,a}<std::pair{nodes[b].anchor.query_begin,b};});
    Position e=0;for(auto i:ids){e=std::max(e,nodes[i].anchor.query_begin+nodes[i].anchor.length);prefix[s].push_back(e);}
  }
  auto visit=[&](unsigned sign,Position b,Position e,auto fn) {
    const auto& ids=roots[sign];const auto& ends=prefix[sign];
    auto k=static_cast<std::size_t>(std::upper_bound(ends.begin(),ends.end(),b)-ends.begin());
    for(;k<ids.size() && nodes[ids[k]].anchor.query_begin<e;++k) {
      check();const auto i=ids[k];const auto& a=nodes[i].anchor;
      if(a.query_begin+a.length>b && !fn(i))break;
    }
  };
  // Independent, two-sided collinearity, with the existing RaMA diagonal cost.
  auto coherent=[&](std::size_t i,std::size_t j) {
    if(i==j)return false;
    const auto* a=&nodes[i].anchor;const auto* b=&nodes[j].anchor;
    if(a->reference_begin>b->reference_begin)std::swap(a,b);
    const auto sign=nodes[i].strand;
    const Position au=sign==Strand::Forward?a->query_begin:query.size()-a->query_begin-a->length;
    const Position bu=sign==Strand::Forward?b->query_begin:query.size()-b->query_begin-b->length;
    if(a->reference_begin+a->length>b->reference_begin || au+a->length>bu)return false;
    if(b->reference_begin-a->reference_begin-a->length>neighborhood || bu-au-a->length>neighborhood)return false;
    const auto dr=b->reference_begin-a->reference_begin,dq=bu-au,delta=dr>dq?dr-dq:dq-dr;
    return !delta || 2.0*std::log2(static_cast<double>(delta)+1.0)<std::min(a->weight,b->weight);
  };
  using Pointers=std::pmr::vector<const CentromereUnifiedPivot*>;
  std::array<Pointers,4> ends{Pointers(memory),Pointers(memory),Pointers(memory),Pointers(memory)};
  auto position=[](const auto* p){return p->from_left?p->query_end:p->query_begin;};
  for(const auto& p:pivots) {
    check();if(!p.id || !p.source_id || p.source_kind>2 || p.query_begin>=p.query_end || p.query_end>query.size() ||
       p.reference_begin>=p.reference_end || p.reference_end-p.reference_begin!=p.query_end-p.query_begin ||
       (p.strand!=Strand::Forward && p.strand!=Strand::Reverse))throw AlignmentError("invalid isolated review pivot");
    ends[(p.strand==Strand::Forward?0U:2U)+(p.from_left?0U:1U)].push_back(&p);
  }
  for(auto& rows:ends)std::sort(rows.begin(),rows.end(),[&](auto a,auto b){return std::pair{position(a),a->id}<std::pair{position(b),b->id};});
  std::uint64_t next=0;
  for(unsigned sign=0;sign<2;++sign) {
    const auto known=sign==0?known_forward:known_reverse,opposite=sign==0?known_reverse:known_forward;
    const auto& lefts=ends[2*sign];const auto& rights=ends[2*sign+1];const unsigned forbidden=sign==0?6U:5U;
    for(std::size_t k=1;k<known.size();++k) {
      check();const auto& a=known[k-1];const auto& z=known[k];const auto b=a.end,e=z.begin;
      if(intersects(opposite,b,e))continue;
      CentromereLocalGapReview row;row.id=++next;row.begin=b;row.end=e;row.strand=sign==0?Strand::Forward:Strand::Reverse;
      auto decide=[&]() -> const char* {
        auto li=std::upper_bound(lefts.begin(),lefts.end(),b,[&](auto p,auto x){return p<position(x);});
        const auto ri=std::lower_bound(rights.begin(),rights.end(),e,[&](auto x,auto p){return position(x)<p;});
        if(li==lefts.begin() || ri==rights.end())return "missing_original_endpoint";
        --li;const auto l=position(*li),r=position(*ri);
        if(l<=a.begin || r>=z.end)return "missing_original_endpoint";
        li=std::lower_bound(lefts.begin(),lefts.end(),l,[&](auto x,auto p){return position(x)<p;});
        row.left_pivot=(*li)->id;row.right_pivot=(*ri)->id;row.corridor_begin=l;row.corridor_end=r;
        if(r-l>limit && r-l-limit>limit)return "original_endpoint_limit";
        if(intersects(opposite,l,r))return "opposite_known_direction";
        if(query.substr(l,r-l).find_first_not_of("ACGT")!=std::string_view::npos)return "contains_N";
        auto w=std::lower_bound(labels.begin(),labels.end(),l,[](const auto& x,auto p){return x.end<=p;});
        for(;w!=labels.end() && w->begin<r;++w)if(w->mask & forbidden)return "opposing_or_competitive_window";
        bool in_gap=false,clear=true;
        visit(1U-sign,l,r,[&](auto i) {
          const auto& root=nodes[i].anchor;
          in_gap=in_gap || (root.query_begin<e && root.query_begin+root.length>b);
          CentromereOppositeReview counter;counter.review_id=row.id;counter.normalized_id=i+1;
          counter.reference_begin=root.reference_begin;counter.query_begin=root.query_begin;
          counter.length=root.length;counter.weight=root.weight;
          if(root.query_begin<l || root.query_begin+root.length>r)counter.verdict=1;
          else if(intersects(opposite,root.query_begin,root.query_begin+root.length))counter.verdict=2;
          else if(root.length>=(*li)->query_end-(*li)->query_begin || root.length>=(*ri)->query_end-(*ri)->query_begin)counter.verdict=3;
          else {
            auto low=root.query_begin-std::min(neighborhood,root.query_begin);
            const auto end=root.query_begin+root.length;
            auto high=end+std::min(neighborhood,static_cast<Length>(query.size())-end);
            if(low)--low;
            if(high<query.size())++high; // Include exact neighborhood-boundary partners.
            visit(1U-sign,low,high,[&](auto j) {
              if(coherent(i,j)){counter.neighbor_id=j+1;counter.verdict=4;return false;}return true;
            });
          }
          emit_counter(counter);clear=clear && !counter.verdict;return true;
        });
        if(!in_gap)return "no_opposite_root_in_gap";
        if(!clear)return "opposite_evidence_retains_veto";
        row.accepted=true;return "isolated_opposite_locally_outsupported";
      };
      row.reason=decide();emit(row);
    }
  }
}
} // namespace ramag::detail

// orientation_known_counterchain
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

// orientation_local_direction_repair
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


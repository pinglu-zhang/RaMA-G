#pragma once
#include "ramag/centromere_blocks.hpp"
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <tuple>

namespace ramag::detail {
// Frozen unknown gaps are reviewed independently. Accepted inferred intervals
// are never fed back as endpoints, root candidates or search-task generators.
template<class Search>
void ReviewOrientationRegions(std::string_view reference,std::string_view query,
    CentromereBlockResult& out,const CentromereBlockOptions& options,Search search) {
  if(!options.review_region_chains)return;
  if(!options.count_inferred_gap || !options.region_search_axis_limit ||
      options.count_inferred_gap>std::numeric_limits<Length>::max()/2)
    throw AlignmentError("invalid region review bounds");
  const auto started=std::chrono::steady_clock::now();
  out.review_region_chains=true;out.supplement_region_anchors=options.supplement_region_anchors;
  const auto limit=options.count_inferred_gap;
  auto check=[&]{if(options.interruption_callback)options.interruption_callback("centromere-region-chain-review");};
  using Interval=CentromereQueryInterval;
  auto merge=[](auto& v) {
    std::sort(v.begin(),v.end(),[](const auto& a,const auto& b){return std::pair{a.begin,a.end}<std::pair{b.begin,b.end};});
    std::size_t n=0;for(const auto& x:v)if(x.begin<x.end) {
      if(n && x.begin<=v[n-1].end)v[n-1].end=std::max(v[n-1].end,x.end);else v[n++]=x;
    }v.resize(n);
  };
  auto overlap=[](const auto& v,Position b,Position e) {
    const auto p=std::lower_bound(v.begin(),v.end(),b,[](const auto& x,auto y){return x.end<=y;});
    return p!=v.end() && p->begin<e;
  };
  std::array<std::vector<Interval>,2> known;
  auto add=[&](const auto& rows){for(const auto& x:rows)known[x.strand==Strand::Forward?0U:1U].push_back({x.begin,x.end});};
  add(out.orientation_known_intervals);add(out.orientation_intervals);add(out.orientation_long_gap_intervals);
  add(out.orientation_sparse_intervals);add(out.orientation_short_gap_intervals);add(out.orientation_mixed_gap_intervals);
  add(out.orientation_sparse_gap_intervals);add(out.orientation_unified_intervals);
  auto add_reviews=[&](const auto& rows){for(const auto& x:rows)if(x.accepted)known[x.strand==Strand::Forward?0U:1U].push_back({x.begin,x.end});};
  add_reviews(out.orientation_local_reviews);add_reviews(out.orientation_isolated_reviews);
  add_reviews(out.orientation_collinear_reviews);add_reviews(out.orientation_alternate_reviews);
  for(auto& v:known)merge(v);
  for(const auto& x:known[0])if(overlap(known[1],x.begin,x.end))throw AlignmentError("region input known directions overlap");
  const auto frozen_unknown=out.orientation_unknown;
  auto memory_check=[&] {
    const std::uint64_t bytes=out.root_candidates.capacity()*sizeof(CentromerePlanningAnchor)+
      out.region_nodes.capacity()*sizeof(CentromereRegionNode)*6+
      out.region_reviews.capacity()*sizeof(CentromereRegionReview)*3+
      out.region_witnesses.capacity()*sizeof(CentromereCollinearWitness)*2+
      out.region_local_candidates.capacity()*sizeof(CentromereRegionLocalCandidate)*3;
    if(bytes>options.planner_memory_limit_bytes)throw AlignmentError("region review capacity budget exceeded; no truncation");
    out.region_capacity_bytes=std::max(out.region_capacity_bytes,bytes);
  };
  auto normalize=[&](const std::vector<CentromerePlanningAnchor>& input,std::uint64_t task) {
    struct Piece {Position r,u;Length n;double density;std::size_t id;Strand strand;};
    if(input.size()>options.planner_memory_limit_bytes/(sizeof(Piece)*8))throw AlignmentError("region normalization budget exceeded");
    std::vector<Piece> pieces;pieces.reserve(input.size());
    for(std::size_t i=0;i<input.size();++i) {
      check();const auto& x=input[i];const auto& a=x.anchor;
      if(!a.length || a.reference_begin>reference.size() || a.length>reference.size()-a.reference_begin ||
         a.query_begin>query.size() || a.length>query.size()-a.query_begin || !std::isfinite(a.weight) || a.weight<0)
        throw AlignmentError("invalid region source anchor");
      pieces.push_back({a.reference_begin,x.strand==Strand::Forward?a.query_begin:query.size()-a.query_begin-a.length,
        a.length,a.weight/static_cast<double>(a.length),i,x.strand});
    }
    auto diagonal=[](const auto& a){return std::pair{a.r>=a.u,a.r>=a.u?a.r-a.u:a.u-a.r};};
    std::sort(pieces.begin(),pieces.end(),[&](const auto& a,const auto& b){return
      std::tuple{a.strand,diagonal(a),a.r,a.u,a.n,a.id}<std::tuple{b.strand,diagonal(b),b.r,b.u,b.n,b.id};});
    struct Event {Position p;double density;bool add;};
    const auto first=out.region_nodes.size();
    for(std::size_t i=0;i<pieces.size();) {
      check();std::size_t e=i+1;auto finish=pieces[i].r+pieces[i].n;
      while(e<pieces.size() && pieces[e].strand==pieces[i].strand && diagonal(pieces[e])==diagonal(pieces[i]) && pieces[e].r<=finish) {
        finish=std::max(finish,pieces[e].r+pieces[e].n);++e;
      }
      std::vector<Event> events;events.reserve(2*(e-i));
      for(auto k=i;k<e;++k){events.push_back({pieces[k].r,pieces[k].density,true});events.push_back({pieces[k].r+pieces[k].n,pieces[k].density,false});}
      std::sort(events.begin(),events.end(),[](const auto& a,const auto& b){return a.p<b.p;});
      std::map<double,std::size_t> active;auto previous=pieces[i].r;double weight=0;
      for(const auto& ev:events) {
        if(!active.empty())weight+=static_cast<double>(ev.p-previous)*active.rbegin()->first;
        if(ev.add)++active[ev.density];else {auto it=active.find(ev.density);if(it==active.end())throw AlignmentError("region normalization imbalance");if(!--it->second)active.erase(it);}
        previous=ev.p;
      }
      CentromereRegionNode node;node.id=out.region_nodes.size()+1;node.root_id=task?0:pieces[i].id+1;node.search_id=task;
      node.source=input[pieces[i].id];node.source.anchor={pieces[i].r,pieces[i].strand==Strand::Forward?pieces[i].u:query.size()-pieces[i].u-(finish-pieces[i].r),finish-pieces[i].r,weight};
      if(e-i>1){node.source.reference_count=0;node.source.query_count=0;node.source.rare_level=0;}
      out.region_nodes.push_back(node);memory_check();i=e;
    }
    std::vector<std::size_t> order;order.reserve(out.region_nodes.size()-first);
    for(auto i=first;i<out.region_nodes.size();++i)order.push_back(i);
    auto box=[&](auto i){const auto& a=out.region_nodes[i].source.anchor;return std::tuple{a.reference_begin,a.query_begin,a.length};};
    std::sort(order.begin(),order.end(),[&](auto a,auto b){return std::pair{box(a),a}<std::pair{box(b),b};});
    for(std::size_t i=0;i<order.size();) {
      auto e=i;unsigned mask=0;while(e<order.size() && box(order[e])==box(order[i])){mask|=out.region_nodes[order[e]].source.strand==Strand::Forward?1U:2U;++e;}
      if(mask==3)for(auto k=i;k<e;++k)out.region_nodes[order[k]].neutral=true;
      i=e;
    }
  };
  normalize(out.root_candidates,0);
  const auto root_count=out.region_nodes.size();
  // Only exact mapping descriptions are retained. Per-gap vectors contain IDs.
  std::array<std::vector<std::size_t>,2> roots;
  for(std::size_t i=0;i<root_count;++i)if(!out.region_nodes[i].neutral)
    roots[out.region_nodes[i].source.strand==Strand::Forward?0U:1U].push_back(i);
  for(auto& ids:roots)std::sort(ids.begin(),ids.end(),[&](auto a,auto b){const auto& x=out.region_nodes[a].source.anchor;const auto& y=out.region_nodes[b].source.anchor;return std::tuple{x.query_begin,x.reference_begin,a}<std::tuple{y.query_begin,y.reference_begin,b};});
  std::array<std::vector<Position>,2> prefix;
  for(unsigned s=0;s<2;++s){Position end=0;for(auto i:roots[s]){const auto& a=out.region_nodes[i].source.anchor;end=std::max(end,a.query_begin+a.length);prefix[s].push_back(end);}}
  auto visit=[&](unsigned s,Position b,Position e,auto fn) {
    auto k=static_cast<std::size_t>(std::upper_bound(prefix[s].begin(),prefix[s].end(),b)-prefix[s].begin());
    for(;k<roots[s].size();++k){const auto i=roots[s][k];const auto& a=out.region_nodes[i].source.anchor;if(a.query_begin>=e)break;if(a.query_begin+a.length>b && !fn(i))break;}
  };
  struct Gap {Position b,e,lb,le,rb,re;unsigned sign;};
  std::vector<Gap> gaps;
  for(unsigned s=0;s<2;++s)for(std::size_t k=1;k<known[s].size();++k) {
    const auto& l=known[s][k-1];const auto& r=known[s][k];
    if(l.end<r.begin && !overlap(known[1-s],l.end,r.begin))gaps.push_back({l.end,r.begin,l.begin,l.end,r.begin,r.end,s});
  }
  std::sort(gaps.begin(),gaps.end(),[](const auto& a,const auto& b){return std::tuple{a.b,a.e,a.sign}<std::tuple{b.b,b.e,b.sign};});
  auto transition=[&](const auto& a,const auto& b,unsigned s,double& cost) {
    if(a.query_begin+a.length>b.query_begin)return false;
    Length rg=0;
    if(!s){if(a.reference_begin+a.length>b.reference_begin)return false;rg=b.reference_begin-a.reference_begin-a.length;}
    else {if(b.reference_begin+b.length>a.reference_begin)return false;rg=a.reference_begin-b.reference_begin-b.length;}
    const auto qg=b.query_begin-a.query_begin-a.length,d=rg>qg?rg-qg:qg-rg;
    cost=2.0*std::log2(static_cast<double>(d)+1.0);return !d || cost<std::min(a.weight,b.weight);
  };
  auto window_veto=[&](Position b,Position e,unsigned s) -> const char* {
    if(query.substr(b,e-b).find_first_not_of("ACGT")!=std::string_view::npos)return "contains_N";
    for(const auto& w:out.orientation_windows)if(w.begin<e && w.end>b) {
      if(w.label=='?')return "competitive_window";
      if(w.label==(s?'+':'-'))return "opposite_advantage_window";
    }return nullptr;
  };
  auto candidates=[&](const Gap& g,std::uint64_t task) {
    std::vector<std::size_t> ids;
    // The outer accepted anchors may be farther than 20 kb from an inferred
    // edge. Only the eventual real-anchor-to-real-anchor reaches are bounded.
    visit(g.sign,std::max(g.lb,g.b-std::min<Length>(g.b,2*limit)),std::min(g.re,g.e+std::min<Length>(query.size()-g.e,2*limit)),[&](auto i){const auto& a=out.region_nodes[i].source.anchor;
      if(a.query_begin>=g.lb && a.query_begin+a.length<=g.re)ids.push_back(i);
      return true;});
    if(task)for(auto i=root_count;i<out.region_nodes.size();++i) {
      const auto& x=out.region_nodes[i];if(x.search_id==task && !x.neutral && x.source.strand==(g.sign?Strand::Reverse:Strand::Forward))ids.push_back(i);
    }
    std::sort(ids.begin(),ids.end(),[&](auto a,auto b){const auto& x=out.region_nodes[a].source.anchor;const auto& y=out.region_nodes[b].source.anchor;return std::tuple{x.query_begin,x.reference_begin,x.length,a}<std::tuple{y.query_begin,y.reference_begin,y.length,b};});
    return ids;
  };
  auto audit=[&](const Gap& g,std::uint64_t task,std::uint64_t forced_left,std::uint64_t forced_right) {
    CentromereRegionReview row;row.id=out.region_reviews.size()+1;row.begin=g.b;row.end=g.e;row.strand=g.sign?Strand::Reverse:Strand::Forward;row.search_id=task;
    const auto ids=candidates(g,task);row.candidate_nodes=ids.size();
    constexpr auto none=std::numeric_limits<std::size_t>::max();
    struct State {Length support{};double score{};std::size_t count{},previous{none};};
    std::vector<State> states(ids.size());
    auto better=[](const auto& a,const auto& b){if(a.support!=b.support)return a.support>b.support;if(a.score!=b.score)return a.score>b.score;if(a.count!=b.count)return a.count<b.count;return a.previous<b.previous;};
    std::size_t best=none,completed_at=0;
    std::vector<std::size_t> completed(ids.size());std::iota(completed.begin(),completed.end(),0);
    std::sort(completed.begin(),completed.end(),[&](auto a,auto b){const auto& x=out.region_nodes[ids[a]].source.anchor;const auto& y=out.region_nodes[ids[b]].source.anchor;return std::pair{x.query_begin+x.length,a}<std::pair{y.query_begin+y.length,b};});
    std::multimap<Position,std::size_t> active;

    for(std::size_t i=0;i<ids.size();++i) {
      check();const auto& x=out.region_nodes[ids[i]];const auto& a=x.source.anchor;
      const bool start=a.query_begin>=g.lb && a.query_begin+a.length<=g.le && (!forced_left || x.id==forced_left);
      const bool finish=a.query_begin>=g.rb && a.query_begin+a.length<=g.re && (!forced_right || x.id==forced_right);
      if(start && !x.search_id)states[i]={a.length,a.weight,1,none};
      while(completed_at<completed.size()) {
        const auto j=completed[completed_at];const auto& p=out.region_nodes[ids[j]].source.anchor;
        if(p.query_begin+p.length>a.query_begin)break;
        if(states[j].support)active.emplace(p.query_begin+p.length,j);
        ++completed_at;
      }
      const auto oldest=a.query_begin-std::min<Length>(a.query_begin,2*limit);
      active.erase(active.begin(),active.lower_bound(oldest));
      for(const auto& entry:active) {
        const auto j=entry.second;

        const auto& p=out.region_nodes[ids[j]].source.anchor;
        if(!states[j].support || p.query_begin+p.length>a.query_begin)continue;
        const auto gap=a.query_begin-p.query_begin-p.length;
        if(gap>2*limit)continue;
        ++row.connection_checks;double penalty=0;if(!transition(p,a,g.sign,penalty))continue;
        const State next{states[j].support+a.length,states[j].score+a.weight-penalty,states[j].count+1,j};
        if(better(next,states[i]))states[i]=next;
      }
      if(finish && !x.search_id && states[i].support && (best==none || better(states[i],states[best])))best=i;
    }
    out.region_connection_checks+=row.connection_checks;
    std::vector<std::size_t> chain;
    if(best!=none)for(auto i=best;i!=none;i=states[i].previous)chain.push_back(ids[i]);
    std::reverse(chain.begin(),chain.end());
    auto decide=[&]() -> const char* {
      if(chain.empty())return "no_complete_bounded_chain";
      const auto& first=out.region_nodes[chain.front()];const auto& last=out.region_nodes[chain.back()];
      const auto& a=first.source.anchor;const auto& z=last.source.anchor;
      row.left_pivot=first.id;row.right_pivot=last.id;row.corridor_begin=a.query_begin+a.length;row.corridor_end=z.query_begin;
      row.support=states[best].support;row.chain_score=states[best].score;row.chain_nodes=chain.size();
      if(overlap(known[1-g.sign],row.corridor_begin,row.corridor_end))return "opposite_known_direction";
      if(const auto reason=window_veto(row.corridor_begin,row.corridor_end,g.sign))return reason;
      bool clear=true;
      auto counter=[&](std::size_t i) {
        const auto& n=out.region_nodes[i];const auto& r=n.source.anchor;
        if(n.neutral || n.source.strand==row.strand || r.query_begin>=row.corridor_end || r.query_begin+r.length<=row.corridor_begin)return;
        if(r.query_begin<row.corridor_begin || r.query_begin+r.length>row.corridor_end ||
           overlap(known[1-g.sign],r.query_begin,r.query_begin+r.length) || r.length>=a.length || r.length>=z.length){clear=false;return;}
        auto low=r.query_begin-std::min(r.query_begin,options.count_cluster_gap);
        auto high=std::min<Position>(query.size(),r.query_begin+r.length+std::min<Length>(options.count_cluster_gap,query.size()-r.query_begin-r.length));
        if(low)--low;
        if(high<query.size())++high;
        auto partner=[&](std::size_t j) {
          if(i==j)return true;
          const auto& m=out.region_nodes[j];if(m.neutral || m.source.strand!=n.source.strand)return true;
          const auto* l=&r;const auto* rr=&m.source.anchor;if(l->query_begin>rr->query_begin)std::swap(l,rr);
          if(l->query_begin+l->length>rr->query_begin || rr->query_begin-l->query_begin-l->length>options.count_cluster_gap)return true;
          double penalty=0;if(!transition(*l,*rr,1-g.sign,penalty))return true;
          const auto rg=g.sign?rr->reference_begin-l->reference_begin-l->length:l->reference_begin-rr->reference_begin-rr->length;
          if(rg<=options.count_cluster_gap)clear=false;
          return clear;
        };
        visit(1-g.sign,low,high,partner);
        if(task)for(auto j=root_count;j<out.region_nodes.size() && clear;++j)if(out.region_nodes[j].search_id==task)partner(j);
      };
      visit(1-g.sign,row.corridor_begin,row.corridor_end,[&](auto i){counter(i);return clear;});
      if(task)for(auto i=root_count;i<out.region_nodes.size() && clear;++i)if(out.region_nodes[i].search_id==task)counter(i);
      if(!clear)return "opposite_evidence_retains_veto";
      row.accepted=true;return task?"complete_chain_after_one_local_search":"complete_root_chain";
    };
    row.reason=decide();
    if(!chain.empty()) {std::uint64_t ordinal=0;for(auto i:chain){const auto& x=out.region_nodes[i];const auto& a=x.source.anchor;out.region_witnesses.push_back({row.id,++ordinal,x.id,a.reference_begin,a.query_begin,a.length,a.weight});}}
    if(row.accepted)out.region_intervals.push_back({out.region_intervals.size()+1,g.b,g.e,row.strand,row.id});
    out.region_reviews.push_back(std::move(row));memory_check();return out.region_reviews.back().accepted;
  };
  // First pass uses only the original roots. All local tasks are frozen below
  // before any search is executed; new anchors never create another task.
  for(const auto& g:gaps)audit(g,0,0,0);
  if(options.supplement_region_anchors) {
    const auto root_review_count=out.region_reviews.size();
    for(std::size_t k=0;k<gaps.size();++k)if(!out.region_reviews[k].accepted) {
      check();const auto& g=gaps[k];CentromereRegionSearch task;task.id=out.region_searches.size()+1;task.review_id=k+1;task.strand=g.sign?Strand::Reverse:Strand::Forward;
      task.query_begin=g.b;task.query_end=g.e;
      const auto ids=candidates(g,0);std::set<std::tuple<Position,Position,Position,Position>> rectangles;
      for(auto li:ids) {
        const auto& l=out.region_nodes[li].source.anchor;
        if(l.query_begin<g.lb || l.query_begin+l.length>g.le)continue;
        for(auto ri:ids) {
          const auto& r=out.region_nodes[ri].source.anchor;
          if(r.query_begin<g.rb || r.query_begin+r.length>g.re)continue;
          const auto qb=l.query_begin+l.length,qe=r.query_begin;
          if(qe<=qb || qe-qb>options.region_search_axis_limit)continue;
          double penalty=0;if(!transition(l,r,g.sign,penalty))continue;
          const auto rb=g.sign?r.reference_begin+r.length:l.reference_begin+l.length;
          const auto re=g.sign?l.reference_begin:r.reference_begin;
          if(re<=rb || re-rb>options.region_search_axis_limit)continue;
          const auto inserted=rectangles.emplace(rb,re,qb,qe).second;
          if(inserted && rectangles.size()==1){task.reference_begin=rb;task.reference_end=re;task.query_begin=qb;task.query_end=qe;task.left_node=out.region_nodes[li].id;task.right_node=out.region_nodes[ri].id;}
          if(rectangles.size()>1)break;
        }
        if(rectangles.size()>1)break;
      }
      if(rectangles.empty())task.reason="no_bounded_two_sided_rectangle";
      else if(rectangles.size()>1)task.reason="ambiguous_endpoint_correspondence";
      else if(const auto reason=window_veto(task.query_begin,task.query_end,g.sign))task.reason=reason;
      else if(reference.substr(task.reference_begin,task.reference_end-task.reference_begin).find_first_not_of("ACGT")!=std::string_view::npos)task.reason="reference_contains_N";
      else {task.reason="pending";task.input_bases=2*((task.reference_end-task.reference_begin)+(task.query_end-task.query_begin));}
      out.region_searches.push_back(std::move(task));
    }
    std::vector<std::size_t> order(out.region_searches.size());std::iota(order.begin(),order.end(),0);
    std::sort(order.begin(),order.end(),[&](auto a,auto b){const auto& x=out.region_searches[a];const auto& y=out.region_searches[b];return std::tuple{x.input_bases,x.query_begin,x.reference_begin,x.id}<std::tuple{y.input_bases,y.query_begin,y.reference_begin,y.id};});
    for(auto i:order) {
      auto& task=out.region_searches[i];if(task.reason!="pending")continue;
      if(task.input_bases>options.region_search_input_budget-out.region_search_input_bases){task.reason="frozen_input_budget";continue;}
      out.region_search_input_bases+=task.input_bases;
      const auto search_start=std::chrono::steady_clock::now();
      std::vector<CentromerePlanningAnchor> local;
      for(auto strand:{Strand::Forward,Strand::Reverse}) {
        check();auto found=search(task,strand);++out.region_search_calls;
        for(auto& a:found) {
          if(a.strand!=strand || a.anchor.reference_begin<task.reference_begin || a.anchor.reference_begin+a.anchor.length>task.reference_end ||
             a.anchor.query_begin<task.query_begin || a.anchor.query_begin+a.anchor.length>task.query_end)throw AlignmentError("local region evidence outside frozen rectangle");
          out.region_local_candidates.push_back({out.region_local_candidates.size()+1,task.id,a});
          local.push_back(std::move(a));
        }
      }
      task.candidate_count=local.size();normalize(local,task.id);task.searched=true;task.reason="searched_once";
      out.region_search_seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-search_start).count();
    }
    // Exactly one final review per searched task, retaining its original ends.
    for(const auto& task:out.region_searches)if(task.searched) {
      if(task.review_id>root_review_count)throw AlignmentError("region task is not rooted in frozen first pass");
      audit(gaps[task.review_id-1],task.id,task.left_node,task.right_node);
    }
  }
  std::vector<Interval> accepted;for(const auto& x:out.region_intervals)accepted.push_back({x.begin,x.end});merge(accepted);
  out.orientation_unknown.clear();
  for(const auto& u:frozen_unknown) {
    Position at=u.begin;for(const auto& a:accepted){if(a.end<=at)continue;if(a.begin>=u.end)break;if(at<a.begin)out.orientation_unknown.push_back({at,std::min(a.begin,u.end)});at=std::max(at,std::min(a.end,u.end));}
    if(at<u.end)out.orientation_unknown.push_back({at,u.end});
  }
  for(const auto& a:accepted) {
    Length inside=0;for(const auto& u:frozen_unknown)inside+=std::min(a.end,u.end)>std::max(a.begin,u.begin)?std::min(a.end,u.end)-std::max(a.begin,u.begin):0;
    if(inside!=a.end-a.begin)throw AlignmentError("region review changed an existing direction");
  }
  out.region_chain_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
}
} // namespace ramag::detail

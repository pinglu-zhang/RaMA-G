#pragma once
#include "centromere_orientation_isolated_review.hpp"
#include <map>

namespace ramag::detail {
// One frozen query-only overlay. A witness is an alternative exact mapping,
// not a selected reference copy or an alignment of the intervening sequence.
template<class Nodes,class Neutral,class Emit,class EmitCounter,class EmitWitness,class Check,
    class EmitAlternate,class EmitAlternateCounter,class EmitAlternateWitness>
void VisitCollinearSameDirectionGapsImpl(const Nodes& nodes,const Neutral& neutral,
    std::span<const CentromereUnifiedPivot> pivots,std::span<const OrientationCoverage> labels,
    std::string_view query,Length limit,Length neighborhood,
    std::span<const CentromereQueryInterval> known_forward,
    std::span<const CentromereQueryInterval> known_reverse,
    std::pmr::memory_resource* memory,Emit emit,EmitCounter emit_counter,
    EmitWitness emit_witness,Check check,bool review_window_gaps,bool review_alternates,bool review_neutral_windows,
    EmitAlternate emit_alternate,EmitAlternateCounter emit_alternate_counter,
    EmitAlternateWitness emit_alternate_witness) {
  if(!limit || !neighborhood || !memory || neutral.size()!=nodes.size())
    throw AlignmentError("invalid collinear review configuration");
  auto intersects=[](auto rows,Position b,Position e) {
    const auto it=std::lower_bound(rows.begin(),rows.end(),b,[](const auto& x,auto p){return x.end<=p;});
    return it!=rows.end() && it->begin<e;
  };
  for(auto rows:{known_forward,known_reverse})for(std::size_t i=0;i<rows.size();++i) {
    check();if(rows[i].begin>=rows[i].end || rows[i].end>query.size() || (i && rows[i-1].end>=rows[i].begin))
      throw AlignmentError("collinear review requires merged known intervals");
  }
  for(const auto& x:known_forward)if(intersects(known_reverse,x.begin,x.end))
    throw AlignmentError("overlapping known directions");
  for(std::size_t i=0;i<labels.size();++i)if(labels[i].begin>=labels[i].end || labels[i].end>query.size() ||
      labels[i].mask>7U || (i && labels[i-1].end>labels[i].begin))throw AlignmentError("invalid collinear review labels");
  using Ids=std::pmr::vector<std::size_t>;
  std::array<Ids,2> roots{Ids(memory),Ids(memory)};
  std::array<std::pmr::vector<Position>,2> prefix{std::pmr::vector<Position>(memory),std::pmr::vector<Position>(memory)};
  for(std::size_t i=0;i<nodes.size();++i) {
    check();const auto& n=nodes[i];const auto& a=n.anchor;
    if(!a.length || a.query_begin>query.size() || a.length>query.size()-a.query_begin ||
       a.reference_begin>std::numeric_limits<Position>::max()-a.length || !std::isfinite(a.weight) || a.weight<0 ||
       (n.strand!=Strand::Forward && n.strand!=Strand::Reverse))throw AlignmentError("invalid collinear review root");
    if(!neutral[i])roots[n.strand==Strand::Forward?0U:1U].push_back(i);
  }
  for(unsigned s=0;s<2;++s) {
    auto& ids=roots[s];std::sort(ids.begin(),ids.end(),[&](auto a,auto b) {
      return std::pair{nodes[a].anchor.query_begin,a}<std::pair{nodes[b].anchor.query_begin,b};});
    Position end=0;for(auto i:ids){end=std::max(end,nodes[i].anchor.query_begin+nodes[i].anchor.length);prefix[s].push_back(end);}
  }
  auto visit=[&](unsigned sign,Position b,Position e,auto fn) {
    const auto& ids=roots[sign];const auto& ends=prefix[sign];
    auto k=static_cast<std::size_t>(std::upper_bound(ends.begin(),ends.end(),b)-ends.begin());
    for(;k<ids.size() && nodes[ids[k]].anchor.query_begin<e;++k) {
      check();const auto i=ids[k];const auto& a=nodes[i].anchor;
      if(a.query_begin+a.length>b && !fn(i))break;
    }
  };
  auto opposing_partner=[&](std::size_t i,std::size_t j) {
    if(i==j)return false;
    const auto* a=&nodes[i].anchor;const auto* b=&nodes[j].anchor;
    if(a->reference_begin>b->reference_begin)std::swap(a,b);
    const auto sign=nodes[i].strand;
    const Position au=sign==Strand::Forward?a->query_begin:query.size()-a->query_begin-a->length;
    const Position bu=sign==Strand::Forward?b->query_begin:query.size()-b->query_begin-b->length;
    if(a->reference_begin+a->length>b->reference_begin || au+a->length>bu)return false;
    const auto rg=b->reference_begin-a->reference_begin-a->length,qg=bu-au-a->length;
    if(rg>neighborhood || qg>neighborhood)return false;
    const auto delta=rg>qg?rg-qg:qg-rg;
    return !delta || 2.0*std::log2(static_cast<double>(delta)+1.0)<std::min(a->weight,b->weight);
  };
  using Pointers=std::pmr::vector<const CentromereUnifiedPivot*>;
  std::array<Pointers,4> ends{Pointers(memory),Pointers(memory),Pointers(memory),Pointers(memory)};
  auto position=[](const auto* p){return p->from_left?p->query_end:p->query_begin;};
  for(const auto& p:pivots) {
    check();if(!p.id || !p.source_id || p.source_kind>2 || p.query_begin>=p.query_end || p.query_end>query.size() ||
       p.reference_begin>=p.reference_end || p.reference_end-p.reference_begin!=p.query_end-p.query_begin ||
       (p.strand!=Strand::Forward && p.strand!=Strand::Reverse))throw AlignmentError("invalid collinear review pivot");
    ends[(p.strand==Strand::Forward?0U:2U)+(p.from_left?0U:1U)].push_back(&p);
  }
  for(auto& rows:ends)std::sort(rows.begin(),rows.end(),[&](auto a,auto b){return std::pair{position(a),a->id}<std::pair{position(b),b->id};});
  constexpr auto none=std::numeric_limits<std::size_t>::max();
  struct Piece {std::size_t id{};Position rb{},qb{};Length length{};double weight{};bool start{},finish{};};
  struct State {Length support{};double score{};std::size_t count{},previous{none};};
  auto audit=[&](CentromereCollinearReview& row,const CentromereUnifiedPivot& lp,
      const CentromereUnifiedPivot& rp,unsigned sign,auto audit_counter,auto audit_witness) -> const char* {
        const auto b=row.begin,e=row.end,l=lp.query_end,r=rp.query_begin;
        const auto opposite=sign==0?known_reverse:known_forward;
        const unsigned own=sign==0?1U:2U,opposing=sign==0?2U:1U;
        row.left_pivot=lp.id;row.right_pivot=rp.id;row.corridor_begin=l;row.corridor_end=r;
        if(r-l>limit && r-l-limit>limit)return "original_endpoint_limit";
        if(intersects(opposite,l,r))return "opposite_known_direction";
        if(query.substr(l,r-l).find_first_not_of("ACGT")!=std::string_view::npos)return "contains_N";
        Position covered=b;bool competition=false,own_evidence=false;
        auto w=std::lower_bound(labels.begin(),labels.end(),l,[](const auto& x,auto p){return x.end<=p;});
        for(;w!=labels.end() && w->begin<r;++w) {
          if(w->mask & opposing)return "opposite_advantage_window";
          competition=competition || (w->mask & 4U);
          own_evidence=own_evidence || ((w->mask & own) && w->begin<e && w->end>b);
          if((w->mask & own) && w->begin<=covered && w->end>covered)covered=std::min(e,w->end);
        }
        // An uncontested absence of window votes is not a direction claim.
        // The opt-in path only permits the SAME counter-evidence and complete
        // two-flank mapping audit below; inferred edges never become pivots.
        const bool window_gap=covered<e;
        const bool neutral_window_gap=window_gap && !own_evidence;
        if(window_gap) {
          if(!review_window_gaps || (neutral_window_gap && !review_neutral_windows))
            return "noncontinuous_own_window_support";
          if(neutral_window_gap && competition)return "neutral_window_competition";
        }
        row.competitive=competition;
        bool clear=true;
        visit(1U-sign,l,r,[&](auto i) {
          const auto& root=nodes[i].anchor;
          CentromereOppositeReview counter;counter.review_id=row.id;counter.normalized_id=i+1;
          counter.reference_begin=root.reference_begin;counter.query_begin=root.query_begin;counter.length=root.length;counter.weight=root.weight;
          if(root.query_begin<l || root.query_begin+root.length>r)counter.verdict=1;
          else if(intersects(opposite,root.query_begin,root.query_begin+root.length))counter.verdict=2;
          else if(root.length>=lp.query_end-lp.query_begin || root.length>=rp.query_end-rp.query_begin)counter.verdict=3;
          else {
            auto low=root.query_begin-std::min(neighborhood,root.query_begin);
            const auto end=root.query_begin+root.length;
            auto high=end+std::min(neighborhood,static_cast<Length>(query.size())-end);
            if(low)--low;
            if(high<query.size())++high;
            visit(1U-sign,low,high,[&](auto j) {
              if(opposing_partner(i,j)){counter.neighbor_id=j+1;counter.verdict=4;return false;}return true;
            });
          }
          audit_counter(counter);clear=clear && !counter.verdict;return true;
        });
        if(!clear)return "opposite_evidence_retains_veto";
        // Alternative copies are allowed, but both complete original query
        // flank matches must be carried by ONE two-axis nonoverlapping chain.
        std::pmr::vector<Piece> pieces(memory);
        visit(sign,lp.query_begin,rp.query_end,[&](auto i) {
          const auto& root=nodes[i].anchor;const auto end=root.query_begin+root.length;
          const auto qb=std::max(root.query_begin,lp.query_begin),qe=std::min(end,rp.query_end);
          const auto rb=root.reference_begin+(sign==0?qb-root.query_begin:end-qe);
          pieces.push_back({i,rb,qb,qe-qb,root.weight*(static_cast<double>(qe-qb)/static_cast<double>(root.length)),
            root.query_begin<=lp.query_begin && end>=lp.query_end,root.query_begin<=rp.query_begin && end>=rp.query_end});
          return true;
        });
        std::sort(pieces.begin(),pieces.end(),[](const auto& x,const auto& y){return std::tuple{x.qb,x.rb,x.id}<std::tuple{y.qb,y.rb,y.id};});
        row.candidate_nodes=pieces.size();
        std::pmr::vector<State> states(pieces.size(),State{},memory);
        Ids completed(memory);for(std::size_t i=0;i<pieces.size();++i)completed.push_back(i);
        std::sort(completed.begin(),completed.end(),[&](auto x,auto y){return std::pair{pieces[x].qb+pieces[x].length,x}<std::pair{pieces[y].qb+pieces[y].length,y};});
        // Only query-completed states enter the reference range index. No
        // cross-product table, path suffix expansion or candidate truncation.
        std::pmr::multimap<Position,std::size_t> active(memory);
        std::size_t at=0,best=none;
        auto better=[](const State& x,const State& y) {
          if(x.support!=y.support)return x.support>y.support;
          if(x.score!=y.score)return x.score>y.score;
          if(x.count!=y.count)return x.count<y.count;
          return x.previous<y.previous;
        };
        for(std::size_t i=0;i<pieces.size();++i) {
          check();const auto& p=pieces[i];
          while(at<completed.size() && pieces[completed[at]].qb+pieces[completed[at]].length<=p.qb) {
            const auto j=completed[at++];if(states[j].support)active.emplace(sign==0?pieces[j].rb+pieces[j].length:pieces[j].rb,j);
          }
          if(p.start)states[i]={p.length,p.weight,1,none};
          const auto begin=sign==0?active.begin():active.lower_bound(p.rb+p.length);
          const auto finish=sign==0?active.upper_bound(p.rb):active.end();
          for(auto it=begin;it!=finish;++it) {
            check();++row.connection_checks;const auto j=it->second;const auto& prior=pieces[j];
            const auto rg=sign==0?p.rb-prior.rb-prior.length:prior.rb-p.rb-p.length;
            const auto qg=p.qb-prior.qb-prior.length,delta=rg>qg?rg-qg:qg-rg;
            const double cost=2.0*std::log2(static_cast<double>(delta)+1.0);
            if(delta && cost>=std::min(p.weight,prior.weight))continue;
            const State candidate{states[j].support+p.length,states[j].score+p.weight-cost,states[j].count+1,j};
            if(better(candidate,states[i]))states[i]=candidate;
          }
          if(p.finish && states[i].support && (best==none || better(states[i],states[best])))best=i;
        }
        if(best==none)return "no_two_flank_collinear_witness";
        row.support=states[best].support;row.chain_score=states[best].score;row.chain_nodes=states[best].count;
        Ids chain(memory);for(auto i=best;i!=none;i=states[i].previous)chain.push_back(i);
        std::reverse(chain.begin(),chain.end());
        std::uint64_t ordinal=0;
        for(auto i:chain) {
          const auto& p=pieces[i];
          audit_witness(CentromereCollinearWitness{row.id,++ordinal,p.id+1,p.rb,p.qb,p.length,p.weight});
        }
        row.accepted=true;
        if(neutral_window_gap)return "two_flank_collinear_neutral_window_witness";
        return window_gap?"two_flank_collinear_window_gap_witness":"two_flank_collinear_witness";
  };
  struct Pending {
    Position begin{},end{};Strand strand{Strand::Forward};
    std::uint64_t left_pivot{},right_pivot{};std::size_t known_index{};
  };
  std::pmr::vector<Pending> pending(memory);
  std::uint64_t next=0;
  for(unsigned sign=0;sign<2;++sign) {
    const auto known=sign==0?known_forward:known_reverse,opposite=sign==0?known_reverse:known_forward;
    const auto& lefts=ends[2*sign];const auto& rights=ends[2*sign+1];
    for(std::size_t k=1;k<known.size();++k) {
      check();const auto& a=known[k-1];const auto& z=known[k];const auto b=a.end,e=z.begin;
      if(intersects(opposite,b,e))continue;
      CentromereCollinearReview row;row.id=++next;row.begin=b;row.end=e;row.strand=sign==0?Strand::Forward:Strand::Reverse;
      auto decide=[&]() -> const char* {
        auto li=std::upper_bound(lefts.begin(),lefts.end(),b,[&](auto p,auto x){return p<position(x);});
        const auto ri=std::lower_bound(rights.begin(),rights.end(),e,[&](auto x,auto p){return position(x)<p;});
        if(li==lefts.begin() || ri==rights.end())return "missing_original_endpoint";
        --li;const auto l=position(*li),r=position(*ri);
        if(l<=a.begin || r>=z.end)return "missing_original_endpoint";
        li=std::lower_bound(lefts.begin(),lefts.end(),l,[&](auto x,auto p){return position(x)<p;});
        return audit(row,**li,**ri,sign,emit_counter,emit_witness);
      };
      row.reason=decide();emit(row);
      if(review_alternates && !row.accepted)
        pending.push_back({row.begin,row.end,row.strand,row.left_pivot,row.right_pivot,k});
    }
  }
  // The old pass is complete before any new evidence is emitted. Its rows,
  // witnesses, IDs and order are unaffected. Fallbacks see only the frozen
  // known blocks and real pivots, never another fallback's inferred edges.
  std::uint64_t alternate_id=0;
  for(const auto& gap:pending) {
    check();const auto& old=gap;const auto sign=old.strand==Strand::Forward?0U:1U;
    const auto known=sign==0?known_forward:known_reverse;
    const auto& a=known[gap.known_index-1];const auto& z=known[gap.known_index];
    const auto b=old.begin,e=old.end;
    const auto& lefts=ends[2*sign];const auto& rights=ends[2*sign+1];
    // Saturating additions/subtractions avoid overflowing 2*limit. Endpoints
    // outside these ranges cannot possibly cover this gap within two budgets.
    auto maximum_right=b+std::min(limit,static_cast<Length>(query.size())-b);
    maximum_right+=std::min(limit,static_cast<Length>(query.size())-maximum_right);
    auto minimum_left=e-std::min(limit,e);minimum_left-=std::min(limit,minimum_left);
    Pointers eligible_rights(memory);
    auto ri=std::lower_bound(rights.begin(),rights.end(),e,[&](auto p,auto q){return position(p)<q;});
    for(;ri!=rights.end() && position(*ri)<z.end && position(*ri)<=maximum_right;++ri) {
      check();const auto* p=*ri;
      if(p->query_begin>=z.begin && p->query_end<=z.end)eligible_rights.push_back(p);
    }
    if(eligible_rights.empty())continue;
    auto within_limit=[&](Position l,Position r){return r-l<=limit || r-l-limit<=limit;};
    struct Pair {const CentromereUnifiedPivot* left{};std::size_t right{};};
    std::pmr::vector<Pair> heap(memory);
    auto key=[&](const Pair& p) {
      const auto* r=eligible_rights[p.right];const auto l=p.left->query_end;
      return std::tuple{r->query_begin-l,b-l,r->query_begin-e,p.left->id,r->id};
    };
    auto later=[&](const Pair& x,const Pair& y){return key(x)>key(y);};
    auto li=std::upper_bound(lefts.begin(),lefts.end(),a.begin,[&](auto q,auto p){return q<position(p);});
    li=std::lower_bound(li,lefts.end(),minimum_left,[&](auto p,auto q){return position(p)<q;});
    for(;li!=lefts.end() && position(*li)<=b;++li) {
      check();const auto* p=*li;
      if(p->query_begin<a.begin || p->query_end>a.end || !within_limit(p->query_end,eligible_rights.front()->query_begin))continue;
      heap.push_back({p,0});
    }
    std::make_heap(heap.begin(),heap.end(),later);
    while(!heap.empty()) {
      check();std::pop_heap(heap.begin(),heap.end(),later);const auto pair=heap.back();heap.pop_back();
      const auto& lp=*pair.left;const auto& rp=*eligible_rights[pair.right];
      if(lp.id!=old.left_pivot || rp.id!=old.right_pivot) {
        CentromereCollinearReview row;row.id=++alternate_id;row.begin=b;row.end=e;row.strand=old.strand;
        row.reason=audit(row,lp,rp,sign,emit_alternate_counter,emit_alternate_witness);
        emit_alternate(row);if(row.accepted)break;
      }
      const auto following=pair.right+1;
      if(following<eligible_rights.size() && within_limit(lp.query_end,eligible_rights[following]->query_begin)) {
        heap.push_back({pair.left,following});std::push_heap(heap.begin(),heap.end(),later);
      }
    }
  }
}

template<class Nodes,class Neutral,class Emit,class EmitCounter,class EmitWitness,class Check>
void VisitCollinearSameDirectionGaps(const Nodes& nodes,const Neutral& neutral,
    std::span<const CentromereUnifiedPivot> pivots,std::span<const OrientationCoverage> labels,
    std::string_view query,Length limit,Length neighborhood,
    std::span<const CentromereQueryInterval> known_forward,
    std::span<const CentromereQueryInterval> known_reverse,
    std::pmr::memory_resource* memory,Emit emit,EmitCounter emit_counter,
    EmitWitness emit_witness,Check check,bool review_window_gaps=false,bool review_neutral_windows=false) {
  VisitCollinearSameDirectionGapsImpl(nodes,neutral,pivots,labels,query,limit,neighborhood,
    known_forward,known_reverse,memory,emit,emit_counter,emit_witness,check,review_window_gaps,false,review_neutral_windows,
    [](const auto&){},[](const auto&){},[](const auto&){});
}

template<class Nodes,class Neutral,class Emit,class EmitCounter,class EmitWitness,class Check,
    class EmitAlternate,class EmitAlternateCounter,class EmitAlternateWitness>
void VisitCollinearSameDirectionGapsWithAlternateEndpoints(const Nodes& nodes,const Neutral& neutral,
    std::span<const CentromereUnifiedPivot> pivots,std::span<const OrientationCoverage> labels,
    std::string_view query,Length limit,Length neighborhood,
    std::span<const CentromereQueryInterval> known_forward,
    std::span<const CentromereQueryInterval> known_reverse,
    std::pmr::memory_resource* memory,Emit emit,EmitCounter emit_counter,
    EmitWitness emit_witness,Check check,bool review_window_gaps,
    EmitAlternate emit_alternate,EmitAlternateCounter emit_alternate_counter,
    EmitAlternateWitness emit_alternate_witness,bool review_neutral_windows=false) {
  VisitCollinearSameDirectionGapsImpl(nodes,neutral,pivots,labels,query,limit,neighborhood,
    known_forward,known_reverse,memory,emit,emit_counter,emit_witness,check,review_window_gaps,true,review_neutral_windows,
    emit_alternate,emit_alternate_counter,emit_alternate_witness);
}
} // namespace ramag::detail

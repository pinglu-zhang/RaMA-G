#pragma once
#include "ramag/internal/centromere/gap_review.hpp"
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <tuple>

// Related direction checks share a private module; their admission rules
// and template callback paths remain independent.

// orientation_collinear_review
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

// orientation_local_review
namespace ramag::detail {
// Read-only snapshot: an accepted review never supplies another endpoint or
// changes a window label. Nearest facing ORIGINAL pivots dominate farther ones
// for both distance and vetoes, so no endpoint-pair enumeration is necessary.
template<class Emit,class Check>
void VisitLocalCompetitiveGaps(std::span<const CentromereUnifiedPivot> pivots,
    std::span<const OrientationCoverage> labels,std::string_view query,Length limit,
    std::span<const CentromereQueryInterval> forward,std::span<const CentromereQueryInterval> reverse,
    std::span<const CentromereQueryInterval> known_forward,std::span<const CentromereQueryInterval> known_reverse,
    std::pmr::memory_resource* memory,Emit emit,Check check) {
  if(!limit || !memory)throw AlignmentError("local gap review requires positive limit and memory resource");
  auto validate=[&](auto rows) {
    for(std::size_t i=0;i<rows.size();++i) {
      check();const auto& row=rows[i];
      if(row.begin>=row.end || row.end>query.size() || (i && rows[i-1].end>=row.begin))
        throw AlignmentError("local review requires merged interval indices");
    }
  };
  validate(forward);validate(reverse);validate(known_forward);validate(known_reverse);
  for(std::size_t i=0;i<labels.size();++i) {
    const auto& row=labels[i];
    if(row.begin>=row.end || row.end>query.size() || row.mask>7U || (i && labels[i-1].end>row.begin))
      throw AlignmentError("invalid local review window coverage");
  }
  auto intersects=[](auto rows,Position begin,Position end) {
    const auto row=std::lower_bound(rows.begin(),rows.end(),begin,[](const auto& x,auto p){return x.end<=p;});
    return row!=rows.end() && row->begin<end;
  };
  for(const auto& row:known_forward)if(intersects(known_reverse,row.begin,row.end))
    throw AlignmentError("opposite known direction overlap");
  using Index=std::pmr::vector<const CentromereUnifiedPivot*>;
  std::array<Index,4> index{Index(memory),Index(memory),Index(memory),Index(memory)};
  for(const auto& p:pivots) {
    check();
    if(!p.id || !p.source_id || p.source_kind>2 || p.query_begin>=p.query_end || p.query_end>query.size() ||
       p.reference_begin>=p.reference_end || p.reference_end-p.reference_begin!=p.query_end-p.query_begin ||
       (p.strand!=Strand::Forward && p.strand!=Strand::Reverse))
      throw AlignmentError("invalid original local-review endpoint");
    index[(p.strand==Strand::Forward?0U:2U)+(p.from_left?0U:1U)].push_back(&p);
  }
  auto position=[](const auto* p){return p->from_left?p->query_end:p->query_begin;};
  for(auto& ids:index)std::sort(ids.begin(),ids.end(),[&](auto a,auto b){
    return std::pair{position(a),a->id}<std::pair{position(b),b->id};});
  std::uint64_t id=0;
  for(unsigned sign=0;sign<2;++sign) {
    const auto known=sign==0?known_forward:known_reverse;
    const auto opposite_known=sign==0?known_reverse:known_forward,opposite_exact=sign==0?reverse:forward;
    const unsigned own_bit=1U<<sign,opposite_bit=sign==0?2U:1U;
    const auto& lefts=index[2*sign];const auto& rights=index[2*sign+1];
    for(std::size_t i=1;i<known.size();++i) {
      check();const auto& left=known[i-1];const auto& right=known[i];
      const Position begin=left.end,end=right.begin;
      // Adjacent known components must be separated only by unknown positions.
      if(intersects(opposite_known,begin,end))continue;
      CentromereLocalGapReview row;row.id=++id;row.begin=begin;row.end=end;
      row.strand=sign==0?Strand::Forward:Strand::Reverse;
      auto decide=[&]() -> const char* {
        auto li=std::upper_bound(lefts.begin(),lefts.end(),begin,[&](auto p,auto a){return p<position(a);});
        const auto ri=std::lower_bound(rights.begin(),rights.end(),end,[&](auto a,auto p){return position(a)<p;});
        if(li==lefts.begin() || ri==rights.end())return "missing_original_endpoint";
        --li;const Position lb=position(*li),rb=position(*ri);
        // End base must belong to the immediately adjacent same-sign component.
        if(lb<=left.begin || rb>=right.end)return "missing_original_endpoint";
        // Deterministic smallest ID when duplicate positions are supplied.
        li=std::lower_bound(lefts.begin(),lefts.end(),lb,[&](auto a,auto p){return position(a)<p;});
        row.left_pivot=(*li)->id;row.right_pivot=(*ri)->id;
        row.corridor_begin=lb;row.corridor_end=rb;
        const Length distance=rb-lb;
        if(distance>limit && distance-limit>limit)return "original_endpoint_limit";
        if(intersects(opposite_exact,lb,rb))return "opposite_exact_inside_corridor";
        if(intersects(opposite_known,lb,rb))return "opposite_known_inside_corridor";
        if(query.substr(lb,rb-lb).find_first_not_of("ACGT")!=std::string_view::npos)return "contains_N";
        bool competitive=false;Position covered=begin;
        auto w=std::lower_bound(labels.begin(),labels.end(),lb,[](const auto& x,auto p){return x.end<=p;});
        for(;w!=labels.end() && w->begin<rb;++w) {
          check();
          if(w->mask & opposite_bit)return "opposite_advantage_window";
          if(w->end<=begin || w->begin>=end)continue;
          competitive=competitive || (w->mask & 4U);
          if(w->mask & own_bit) {
            if(w->begin>covered)return "same_advantage_not_continuous";
            covered=std::max(covered,std::min(end,w->end));
          }
        }
        if(!competitive)return "not_competitive_gap";
        if(covered!=end)return "same_advantage_not_continuous";
        row.accepted=true;return "local_corridor_clear";
      };
      row.reason=decide();emit(row);
    }
  }
}
} // namespace ramag::detail

// orientation_region_review
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

// orientation_unified_review
namespace ramag::detail {
// Internal source-RaMA adapter shared with differential tests. It returns the
// first-hit pre-chain candidates, in original full-input query coordinates.
std::vector<CentromerePlanningAnchor> SearchCentromereRegionCandidates(
    const SequenceRecord&,const SequenceRecord&,const CentromereRegionSearch&,
    Strand,const CentromereBlockOptions&);
// Frozen unknown gaps are reviewed independently. Accepted inferred intervals
// are never fed back as endpoints, root candidates or search-task generators.
template<class Search>
void ReviewUnifiedOrientationGaps(std::string_view reference,std::string_view query,
    CentromereBlockResult& out,const CentromereBlockOptions& options,Search search) {
  if(!options.review_unified_gaps)return;
  if(!options.count_inferred_gap || !options.unified_search_query_limit || !options.unified_search_reference_limit ||
      options.count_inferred_gap>std::numeric_limits<Length>::max()/2)
    throw AlignmentError("invalid region review bounds");
  const auto started=std::chrono::steady_clock::now();
  out.review_unified_gaps=true;
  out.review_gap_interior_opposition=options.review_gap_interior_opposition;
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
  add(out.region_intervals);
  auto add_reviews=[&](const auto& rows){for(const auto& x:rows)if(x.accepted)known[x.strand==Strand::Forward?0U:1U].push_back({x.begin,x.end});};
  add_reviews(out.orientation_local_reviews);add_reviews(out.orientation_isolated_reviews);
  add_reviews(out.orientation_collinear_reviews);add_reviews(out.orientation_alternate_reviews);
  for(auto& v:known)merge(v);
  for(const auto& x:known[0])if(overlap(known[1],x.begin,x.end))throw AlignmentError("region input known directions overlap");
  const auto frozen_unknown=out.orientation_unknown;
  // normalize() appends one task contiguously. Preserve those bounds so an
  // audit never scans another rectangle's nodes merely to discard their IDs.
  std::vector<std::pair<std::size_t,std::size_t>> local_node_ranges;
  std::vector<bool> legacy_root_accepted;
  auto memory_check=[&](std::uint64_t transient=0) {
    const std::uint64_t bytes=legacy_root_accepted.capacity()/8+1+
      out.root_candidates.capacity()*sizeof(CentromerePlanningAnchor)+
      out.unified_nodes.capacity()*sizeof(CentromereRegionNode)*6+
      out.unified_reviews.capacity()*sizeof(CentromereRegionReview)*3+
      out.unified_witnesses.capacity()*sizeof(CentromereCollinearWitness)*2+
      out.unified_local_candidates.capacity()*sizeof(CentromereRegionLocalCandidate)*3+
      out.unified_gap_decisions.capacity()*sizeof(CentromereUnifiedGapDecision)*2+
      out.unified_hypotheses.capacity()*sizeof(CentromereUnifiedHypothesis)*3+
      out.unified_searches.capacity()*sizeof(CentromereRegionSearch)*3+
      local_node_ranges.capacity()*sizeof(local_node_ranges.front());
    if(transient>options.planner_memory_limit_bytes || bytes>options.planner_memory_limit_bytes-transient)
      throw AlignmentError("unified review capacity budget exceeded; no truncation");
    out.unified_capacity_bytes=std::max(out.unified_capacity_bytes,bytes+transient);
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
    const auto first=out.unified_nodes.size();
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
      CentromereRegionNode node;node.id=out.unified_nodes.size()+1;node.root_id=task?0:pieces[i].id+1;node.search_id=task;
      node.source=input[pieces[i].id];node.source.anchor={pieces[i].r,pieces[i].strand==Strand::Forward?pieces[i].u:query.size()-pieces[i].u-(finish-pieces[i].r),finish-pieces[i].r,weight};
      if(e-i>1){node.source.reference_count=0;node.source.query_count=0;node.source.rare_level=0;}
      out.unified_nodes.push_back(node);memory_check();i=e;
    }
    std::vector<std::size_t> order;order.reserve(out.unified_nodes.size()-first);
    for(auto i=first;i<out.unified_nodes.size();++i)order.push_back(i);
    auto box=[&](auto i){const auto& a=out.unified_nodes[i].source.anchor;return std::tuple{a.reference_begin,a.query_begin,a.length};};
    std::sort(order.begin(),order.end(),[&](auto a,auto b){return std::pair{box(a),a}<std::pair{box(b),b};});
    for(std::size_t i=0;i<order.size();) {
      auto e=i;unsigned mask=0;while(e<order.size() && box(order[e])==box(order[i])){mask|=out.unified_nodes[order[e]].source.strand==Strand::Forward?1U:2U;++e;}
      if(mask==3)for(auto k=i;k<e;++k)out.unified_nodes[order[k]].neutral=true;
      i=e;
    }
    if(task)local_node_ranges.at(task-1)={first,out.unified_nodes.size()};
  };
  normalize(out.root_candidates,0);
  const auto root_count=out.unified_nodes.size();
  // Only exact mapping descriptions are retained. Per-gap vectors contain IDs.
  std::array<std::vector<std::size_t>,2> roots;
  for(std::size_t i=0;i<root_count;++i)if(!out.unified_nodes[i].neutral)
    roots[out.unified_nodes[i].source.strand==Strand::Forward?0U:1U].push_back(i);
  for(auto& ids:roots)std::sort(ids.begin(),ids.end(),[&](auto a,auto b){const auto& x=out.unified_nodes[a].source.anchor;const auto& y=out.unified_nodes[b].source.anchor;return std::tuple{x.query_begin,x.reference_begin,a}<std::tuple{y.query_begin,y.reference_begin,b};});
  std::array<std::vector<Position>,2> prefix;
  for(unsigned s=0;s<2;++s){Position end=0;for(auto i:roots[s]){const auto& a=out.unified_nodes[i].source.anchor;end=std::max(end,a.query_begin+a.length);prefix[s].push_back(end);}}
  auto visit=[&](unsigned s,Position b,Position e,auto fn) {
    auto k=static_cast<std::size_t>(std::upper_bound(prefix[s].begin(),prefix[s].end(),b)-prefix[s].begin());
    for(;k<roots[s].size();++k){const auto i=roots[s][k];const auto& a=out.unified_nodes[i].source.anchor;if(a.query_begin>=e)break;if(a.query_begin+a.length>b && !fn(i))break;}
  };
  struct Gap {Position b,e,lb,le,rb,re;unsigned sign;};
  std::vector<Gap> gaps;
  for(unsigned s=0;s<2;++s)for(std::size_t k=1;k<known[s].size();++k) {
    const auto& l=known[s][k-1];const auto& r=known[s][k];
    if(l.end<r.begin && !overlap(known[1-s],l.end,r.begin))gaps.push_back({l.end,r.begin,l.begin,l.end,r.begin,r.end,s});
  }
  std::sort(gaps.begin(),gaps.end(),[](const auto& a,const auto& b){return std::tuple{a.b,a.e,a.sign}<std::tuple{b.b,b.e,b.sign};});
  // Freeze local-search admission with the previous qualification. A newly
  // accepted root chain must not free budget for a previously unsearched gap.
  legacy_root_accepted.assign(gaps.size(),false);
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
      // A failed count/length dominance test is not itself a competing chain.
      // '?' proceeds to the same actual opposite-anchor audit as clear windows.
      if(w.label==(s?'+':'-'))return "opposite_advantage_window";
    }return nullptr;
  };
  auto candidates=[&](const Gap& g,std::uint64_t task) {
    std::vector<std::size_t> ids;
    // The outer accepted anchors may be farther than 20 kb from an inferred
    // edge. Only the eventual real-anchor-to-real-anchor reaches are bounded.
    visit(g.sign,std::max(g.lb,g.b-std::min<Length>(g.b,2*limit)),std::min(g.re,g.e+std::min<Length>(query.size()-g.e,2*limit)),[&](auto i){const auto& a=out.unified_nodes[i].source.anchor;
      if(a.query_begin>=g.lb && a.query_begin+a.length<=g.re)ids.push_back(i);
      return true;});
    if(task) {
      const auto [first,last]=local_node_ranges.at(task-1);
      for(auto i=first;i<last;++i) {
        const auto& x=out.unified_nodes[i];if(x.search_id==task && !x.neutral && x.source.strand==(g.sign?Strand::Reverse:Strand::Forward))ids.push_back(i);
      }
    }
    std::sort(ids.begin(),ids.end(),[&](auto a,auto b){const auto& x=out.unified_nodes[a].source.anchor;const auto& y=out.unified_nodes[b].source.anchor;return std::tuple{x.query_begin,x.reference_begin,x.length,a}<std::tuple{y.query_begin,y.reference_begin,y.length,b};});
    return ids;
  };
  auto audit=[&](const Gap& g,std::uint64_t task,std::uint64_t forced_left,std::uint64_t forced_right,bool* legacy_accepted=nullptr) {
    CentromereRegionReview row;row.id=out.unified_reviews.size()+1;row.begin=g.b;row.end=g.e;row.strand=g.sign?Strand::Reverse:Strand::Forward;row.search_id=task;
    auto ids=candidates(g,task);
    if(task) {
      const auto& rectangle=out.unified_searches.at(task-1);
      ids.erase(std::remove_if(ids.begin(),ids.end(),[&](auto i) {
        const auto& n=out.unified_nodes[i];const auto& a=n.source.anchor;
        return n.id!=forced_left && n.id!=forced_right &&
          (a.reference_begin<rectangle.reference_begin || a.reference_begin+a.length>rectangle.reference_end ||
           a.query_begin<rectangle.query_begin || a.query_begin+a.length>rectangle.query_end);
      }),ids.end());
    }
    row.candidate_nodes=ids.size();
    constexpr auto none=std::numeric_limits<std::size_t>::max();
    struct State {Length support{};double score{};std::size_t count{},previous{none};};
    std::vector<State> states(ids.size());
    auto better=[](const auto& a,const auto& b){if(a.support!=b.support)return a.support>b.support;if(a.score!=b.score)return a.score>b.score;if(a.count!=b.count)return a.count<b.count;return a.previous<b.previous;};
    std::size_t best=none,completed_at=0;
    std::vector<std::size_t> completed(ids.size());std::iota(completed.begin(),completed.end(),0);
    std::sort(completed.begin(),completed.end(),[&](auto a,auto b){const auto& x=out.unified_nodes[ids[a]].source.anchor;const auto& y=out.unified_nodes[ids[b]].source.anchor;return std::pair{x.query_begin+x.length,a}<std::pair{y.query_begin+y.length,b};});
    std::multimap<Position,std::size_t> active;

    for(std::size_t i=0;i<ids.size();++i) {
      check();const auto& x=out.unified_nodes[ids[i]];const auto& a=x.source.anchor;
      const bool start=a.query_begin>=g.lb && a.query_begin+a.length<=g.le && (!forced_left || x.id==forced_left);
      const bool finish=a.query_begin>=g.rb && a.query_begin+a.length<=g.re && (!forced_right || x.id==forced_right);
      if(start && !x.search_id)states[i]={a.length,a.weight,1,none};
      while(completed_at<completed.size()) {
        const auto j=completed[completed_at];const auto& p=out.unified_nodes[ids[j]].source.anchor;
        if(p.query_begin+p.length>a.query_begin)break;
        if(states[j].support)active.emplace(p.query_begin+p.length,j);
        ++completed_at;
      }
      const auto oldest=a.query_begin-std::min<Length>(a.query_begin,2*limit);
      active.erase(active.begin(),active.lower_bound(oldest));
      for(const auto& entry:active) {
        const auto j=entry.second;

        const auto& p=out.unified_nodes[ids[j]].source.anchor;
        if(!states[j].support || p.query_begin+p.length>a.query_begin)continue;
        const auto gap=a.query_begin-p.query_begin-p.length;
        if(gap>2*limit)continue;
        ++row.connection_checks;double penalty=0;if(!transition(p,a,g.sign,penalty))continue;
        const State next{states[j].support+a.length,states[j].score+a.weight-penalty,states[j].count+1,j};
        if(better(next,states[i]))states[i]=next;
      }
      if(finish && !x.search_id && states[i].support && (best==none || better(states[i],states[best])))best=i;
    }
    out.unified_connection_checks+=row.connection_checks;
    std::vector<std::size_t> chain;
    if(best!=none)for(auto i=best;i!=none;i=states[i].previous)chain.push_back(ids[i]);
    std::reverse(chain.begin(),chain.end());
    auto decide=[&](bool gap_interior) -> const char* {
      // Local opposite evidence must be checked even when its own-direction
      // hypothesis fails to produce a complete chain.
      if(chain.empty() && !task)return "no_complete_bounded_chain";
      const auto& first=out.unified_nodes[chain.empty()?forced_left-1:chain.front()];
      const auto& last=out.unified_nodes[chain.empty()?forced_right-1:chain.back()];
      const auto& a=first.source.anchor;const auto& z=last.source.anchor;
      row.left_pivot=first.id;row.right_pivot=last.id;row.corridor_begin=a.query_begin+a.length;row.corridor_end=z.query_begin;
      if(!chain.empty()){row.support=states[best].support;row.chain_score=states[best].score;row.chain_nodes=chain.size();}
      if(overlap(known[1-g.sign],row.corridor_begin,row.corridor_end))return "opposite_known_direction";
      if(const auto reason=window_veto(row.corridor_begin,row.corridor_end,g.sign))return reason;
      bool clear=true;
      auto counter=[&](std::size_t i) {
        const auto& n=out.unified_nodes[i];const auto& r=n.source.anchor;
        if(n.neutral || n.source.strand==row.strand || r.query_begin>=row.corridor_end || r.query_begin+r.length<=row.corridor_begin)return;
        // A corridor includes already accepted flanks. Context outside the
        // actual unknown gap cannot provide a second independent vote against
        // recovering that gap. Keep the original match intact for all scores.
        const auto touches_gap=[&](const auto& anchor) {
          return anchor.query_begin<g.e && anchor.query_begin+anchor.length>g.b;
        };
        if(gap_interior && !touches_gap(r))return;
        if((!gap_interior &&
            (r.query_begin<row.corridor_begin || r.query_begin+r.length>row.corridor_end)) ||
           overlap(known[1-g.sign],r.query_begin,r.query_begin+r.length) || r.length>=a.length || r.length>=z.length){clear=false;return;}
        auto low=r.query_begin-std::min(r.query_begin,options.count_cluster_gap);
        auto high=std::min<Position>(query.size(),r.query_begin+r.length+std::min<Length>(options.count_cluster_gap,query.size()-r.query_begin-r.length));
        if(low)--low;
        if(high<query.size())++high;
        auto partner=[&](std::size_t j) {
          if(i==j)return true;
          const auto& m=out.unified_nodes[j];if(m.neutral || m.source.strand!=n.source.strand)return true;
          if(gap_interior && !touches_gap(m.source.anchor))return true;
          const auto* l=&r;const auto* rr=&m.source.anchor;if(l->query_begin>rr->query_begin)std::swap(l,rr);
          if(l->query_begin+l->length>rr->query_begin || rr->query_begin-l->query_begin-l->length>options.count_cluster_gap)return true;
          double penalty=0;if(!transition(*l,*rr,1-g.sign,penalty))return true;
          const auto rg=g.sign?rr->reference_begin-l->reference_begin-l->length:l->reference_begin-rr->reference_begin-rr->length;
          if(rg<=options.count_cluster_gap)clear=false;
          return clear;
        };
        visit(1-g.sign,low,high,partner);
        if(task) {
          const auto [first_local,last_local]=local_node_ranges.at(task-1);
          for(auto j=first_local;j<last_local && clear;++j)if(out.unified_nodes[j].search_id==task)partner(j);
        }
      };
      visit(1-g.sign,row.corridor_begin,row.corridor_end,[&](auto i){counter(i);return clear;});
      if(task) {
        const auto [first_local,last_local]=local_node_ranges.at(task-1);
        for(auto i=first_local;i<last_local && clear;++i)if(out.unified_nodes[i].search_id==task)counter(i);
      }
      if(!clear)return "opposite_evidence_retains_veto";
      if(chain.empty())return "no_complete_bounded_chain";
      row.accepted=true;return task?"complete_chain_after_one_local_search":"complete_root_chain";
    };
    if(legacy_accepted && options.review_gap_interior_opposition) {
      // The chain DP and search are not repeated: only qualify the same root
      // evidence twice, retaining old admission and new final evidence.
      (void)decide(false);*legacy_accepted=row.accepted;row.accepted=false;
    }
    row.reason=decide(options.review_gap_interior_opposition);
    if(legacy_accepted && !options.review_gap_interior_opposition)*legacy_accepted=row.accepted;
    if(!chain.empty()) {std::uint64_t ordinal=0;for(auto i:chain){const auto& x=out.unified_nodes[i];const auto& a=x.source.anchor;out.unified_witnesses.push_back({row.id,++ordinal,x.id,a.reference_begin,a.query_begin,a.length,a.weight});}}
    out.unified_reviews.push_back(std::move(row));memory_check();return out.unified_reviews.back().accepted;
  };
  // First pass uses only the original roots. All local tasks are frozen below
  // before any search is executed; new anchors never create another task.
  for(const auto& g:gaps) {
    bool legacy_accepted=false;
    const bool accepted=audit(g,0,0,0,&legacy_accepted);
    legacy_root_accepted[out.unified_gap_decisions.size()]=legacy_accepted;
    const auto review=out.unified_reviews.back().id;
    CentromereUnifiedGapDecision d;d.id=out.unified_gap_decisions.size()+1;d.begin=g.b;d.end=g.e;
    d.strand=g.sign?Strand::Reverse:Strand::Forward;d.root_review_id=review;
    d.accepted=accepted;d.chosen_review_id=accepted?review:0;d.reason=out.unified_reviews.back().reason;
    out.unified_gap_decisions.push_back(std::move(d));
    if(legacy_accepted)out.unified_intervals.push_back({out.unified_intervals.size()+1,g.b,g.e,g.sign?Strand::Reverse:Strand::Forward,review});
  }
  // Exact rectangles share one search, but endpoint hypotheses remain separate.
  // No candidate from a different rectangle can enter a hypothesis chain.
  using Rectangle=std::tuple<Position,Position,Position,Position>;
  std::map<Rectangle,std::uint64_t> cached;
  std::vector<std::vector<std::size_t>> gap_hypotheses(gaps.size());
  for(std::size_t k=0;k<gaps.size();++k)if(!legacy_root_accepted[k]) {
    check();const auto& g=gaps[k];const auto ids=candidates(g,0);
    std::vector<std::size_t> right;
    for(auto i:ids){const auto& a=out.unified_nodes[i].source.anchor;
      if(a.query_begin>=g.rb && a.query_begin+a.length<=g.re)right.push_back(i);}
    auto reference_key=[&](auto i){const auto& a=out.unified_nodes[i].source.anchor;return g.sign?a.reference_begin+a.length:a.reference_begin;};
    std::sort(right.begin(),right.end(),[&](auto a,auto b){return std::pair{reference_key(a),a}<std::pair{reference_key(b),b};});
    for(auto li:ids) {
      const auto& l=out.unified_nodes[li].source.anchor;
      if(l.query_begin<g.lb || l.query_begin+l.length>g.le)continue;
      const auto center=g.sign?l.reference_begin:l.reference_begin+l.length;
      const auto lo=g.sign?center-std::min<Length>(center,options.unified_search_reference_limit):center;
      const auto hi=g.sign?center:center+std::min<Length>(reference.size()-center,options.unified_search_reference_limit);
      auto first=std::lower_bound(right.begin(),right.end(),lo,[&](auto i,auto p){return reference_key(i)<p;});
      for(auto it=first;it!=right.end() && reference_key(*it)<=hi;++it) {
        check();const auto ri=*it;const auto& r=out.unified_nodes[ri].source.anchor;
        const auto qb=l.query_begin+l.length,qe=r.query_begin;
        if(qe<=qb || qe-qb>options.unified_search_query_limit)continue;
        double penalty=0;if(!transition(l,r,g.sign,penalty))continue;
        const auto rb=g.sign?r.reference_begin+r.length:l.reference_begin+l.length;
        const auto re=g.sign?l.reference_begin:r.reference_begin;
        if(re<=rb || re-rb>options.unified_search_reference_limit)continue;
        const Rectangle key{rb,re,qb,qe};auto found=cached.find(key);
        if(found==cached.end()) {
          CentromereRegionSearch task;task.id=out.unified_searches.size()+1;task.review_id=k+1;
          task.left_node=out.unified_nodes[li].id;task.right_node=out.unified_nodes[ri].id;
          task.reference_begin=rb;task.reference_end=re;task.query_begin=qb;task.query_end=qe;
          task.strand=g.sign?Strand::Reverse:Strand::Forward;
          if(const auto reason=window_veto(qb,qe,g.sign))task.reason=reason;
          else if(reference.substr(rb,re-rb).find_first_not_of("ACGT")!=std::string_view::npos)task.reason="reference_contains_N";
          else {task.reason="pending";task.input_bases=2*((re-rb)+(qe-qb));}
          out.unified_searches.push_back(std::move(task));found=cached.emplace(key,out.unified_searches.size()).first;
        }
        CentromereUnifiedHypothesis h;h.id=out.unified_hypotheses.size()+1;h.gap_id=k+1;
        h.left_node=out.unified_nodes[li].id;h.right_node=out.unified_nodes[ri].id;h.search_id=found->second;h.reason="frozen";
        gap_hypotheses[k].push_back(out.unified_hypotheses.size());out.unified_hypotheses.push_back(std::move(h));memory_check();
      }
    }
    out.unified_gap_decisions[k].hypothesis_count=gap_hypotheses[k].size();
  }
  // Reserve complete hypothesis groups. A group that cannot be fully checked
  // remains unknown; cheapest individual copies are not treated as exhaustive.
  local_node_ranges.resize(out.unified_searches.size());memory_check();
  std::vector<std::pair<std::uint64_t,std::size_t>> groups;
  for(std::size_t k=0;k<gaps.size();++k)if(!gap_hypotheses[k].empty()) {
    std::set<std::uint64_t> unique;std::uint64_t cost=0;
    for(auto h:gap_hypotheses[k])unique.insert(out.unified_hypotheses[h].search_id);
    for(auto id:unique){const auto& t=out.unified_searches[id-1];
      if(t.reason=="pending"){
        if(t.input_bases>std::numeric_limits<std::uint64_t>::max()-cost)throw AlignmentError("unified group cost overflow");
        cost+=t.input_bases;}}
    groups.emplace_back(cost,k);
  }
  std::sort(groups.begin(),groups.end());std::vector<bool> reserved(out.unified_searches.size(),false);
  for(const auto& [unused_cost,k]:groups) {
    (void)unused_cost;std::set<std::uint64_t> unique;std::uint64_t cost=0;
    for(auto h:gap_hypotheses[k])unique.insert(out.unified_hypotheses[h].search_id);
    for(auto id:unique)if(!reserved[id-1] && out.unified_searches[id-1].reason=="pending") {
      const auto value=out.unified_searches[id-1].input_bases;
      if(value>std::numeric_limits<std::uint64_t>::max()-cost)throw AlignmentError("unified group cost overflow");
      cost+=value;
    }
    if(cost>options.unified_search_input_budget-out.unified_search_input_bases)continue;
    out.unified_search_input_bases+=cost;for(auto id:unique)if(out.unified_searches[id-1].reason=="pending")reserved[id-1]=true;
  }
  for(std::size_t k=0;k<out.unified_searches.size();++k) {
    auto& task=out.unified_searches[k];if(task.reason!="pending")continue;
    if(!reserved[k]){task.reason="frozen_input_budget";continue;}
    const auto search_start=std::chrono::steady_clock::now();std::vector<CentromerePlanningAnchor> local;
    for(auto strand:{Strand::Forward,Strand::Reverse}) {
      check();auto found=search(task,strand);++out.unified_search_calls;
      // Reserve neither output vector until the incoming raw result and both
      // persistent/normalization copies fit. No candidate is silently dropped.
      constexpr std::uint64_t per_anchor=8*sizeof(CentromerePlanningAnchor);
      if(found.size()>options.planner_memory_limit_bytes/per_anchor ||
         local.size()>options.planner_memory_limit_bytes/per_anchor-found.size())
        throw AlignmentError("unified local candidate capacity budget exceeded; no truncation");
      const auto working_bytes=(found.size()+local.size())*per_anchor;
      auto raw_capacity=out.unified_local_candidates.capacity();
      if(found.size()>out.unified_local_candidates.max_size()-out.unified_local_candidates.size())
        throw AlignmentError("unified local candidate size overflow");
      const auto required=out.unified_local_candidates.size()+found.size();
      if(required>raw_capacity) {
        const auto maximum=out.unified_local_candidates.max_size();
        raw_capacity=std::max(required,raw_capacity>maximum/2?maximum:raw_capacity*2);
      }
      // Account for the proposed capacity before reserve(), including the
      // temporary old/new storage during reallocation. Geometric growth avoids
      // moving every earlier task's raw evidence for each small new task.
      constexpr std::uint64_t raw_bytes=3*sizeof(CentromereRegionLocalCandidate);
      const auto extra_capacity=raw_capacity-out.unified_local_candidates.capacity();
      if(extra_capacity>options.planner_memory_limit_bytes/raw_bytes ||
         working_bytes>options.planner_memory_limit_bytes-extra_capacity*raw_bytes)
        throw AlignmentError("unified local candidate capacity budget exceeded; no truncation");
      memory_check(working_bytes+extra_capacity*raw_bytes);
      local.reserve(local.size()+found.size());
      if(raw_capacity>out.unified_local_candidates.capacity())out.unified_local_candidates.reserve(raw_capacity);
      for(auto& a:found) {
        if(a.strand!=strand || !a.anchor.length || a.anchor.reference_begin<task.reference_begin ||
           a.anchor.reference_begin>task.reference_end || a.anchor.length>task.reference_end-a.anchor.reference_begin ||
           a.anchor.query_begin<task.query_begin || a.anchor.query_begin>task.query_end || a.anchor.length>task.query_end-a.anchor.query_begin)
          throw AlignmentError("unified local evidence outside frozen rectangle");
        out.unified_local_candidates.push_back({out.unified_local_candidates.size()+1,task.id,a});local.push_back(std::move(a));
      }
    }
    task.candidate_count=local.size();normalize(local,task.id);task.searched=true;task.reason="searched_once";
    out.unified_search_seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-search_start).count();
  }
  // One final evaluation per original hypothesis. Never feed inferred output
  // or a newly discovered local anchor back into task generation.
  for(std::size_t k=0;k<gaps.size();++k)if(!legacy_root_accepted[k]) {
    auto& decision=out.unified_gap_decisions[k];bool incomplete=false,missing_audit=false,opposition=false;
    std::uint64_t chosen=decision.accepted?decision.root_review_id:0;
    decision.accepted=false;decision.chosen_review_id=0;
    for(auto hi:gap_hypotheses[k]) {
      auto& h=out.unified_hypotheses[hi];const auto& task=out.unified_searches[h.search_id-1];
      if(!task.searched){
        h.reason=task.reason;
        if(task.reason=="frozen_input_budget")incomplete=true;
        else if(task.reason=="opposite_advantage_window")opposition=true;
        else missing_audit=true;
        continue;
      }
      h.accepted=audit(gaps[k],task.id,h.left_node,h.right_node);h.review_id=out.unified_reviews.back().id;h.reason=out.unified_reviews.back().reason;
      if(h.reason=="opposite_evidence_retains_veto" || h.reason=="opposite_known_direction" || h.reason=="opposite_advantage_window")opposition=true;
      if(h.accepted && !chosen)chosen=h.review_id;
    }
    if(gap_hypotheses[k].empty() && !chosen)decision.reason="no_bounded_endpoint_hypothesis";
    else if(incomplete)decision.reason="incomplete_hypothesis_budget";
    else if(missing_audit)decision.reason="incomplete_hypothesis_audit";
    else if(opposition)decision.reason="competing_hypothesis_evidence";
    else if(!chosen)decision.reason="no_complete_hypothesis_chain";
    else {
      decision.accepted=true;decision.chosen_review_id=chosen;
      decision.reason=gap_hypotheses[k].empty()?"complete_root_chain":"complete_audited_hypotheses";
      out.unified_intervals.push_back({out.unified_intervals.size()+1,decision.begin,decision.end,decision.strand,chosen});
    }
  }
  std::vector<Interval> accepted;for(const auto& x:out.unified_intervals)accepted.push_back({x.begin,x.end});merge(accepted);
  out.orientation_unknown.clear();
  for(const auto& u:frozen_unknown) {
    Position at=u.begin;for(const auto& a:accepted){if(a.end<=at)continue;if(a.begin>=u.end)break;if(at<a.begin)out.orientation_unknown.push_back({at,std::min(a.begin,u.end)});at=std::max(at,std::min(a.end,u.end));}
    if(at<u.end)out.orientation_unknown.push_back({at,u.end});
  }
  for(const auto& a:accepted) {
    Length inside=0;for(const auto& u:frozen_unknown)inside+=std::min(a.end,u.end)>std::max(a.begin,u.begin)?std::min(a.end,u.end)-std::max(a.begin,u.begin):0;
    if(inside!=a.end-a.begin)throw AlignmentError("region review changed an existing direction");
  }
  out.unified_review_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
}
} // namespace ramag::detail


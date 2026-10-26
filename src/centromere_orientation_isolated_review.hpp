#pragma once
#include "centromere_orientation_unified_ends.hpp"
#include <array>
#include <cmath>
#include <numeric>

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

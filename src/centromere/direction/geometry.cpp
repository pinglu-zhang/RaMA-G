#include "ramag/centromere_direction.hpp"
#include <algorithm>
#include <cmath>
#include <map>
#include <tuple>

namespace ramag {
namespace {
struct Piece {
  Position r{},q{};
  Length n{};
  double density{};
  std::uint64_t id{};
};
void ValidateBox(const DirectionBox& b, Length r, Length q) {
  if (b.reference_begin>b.reference_end || b.reference_end>r ||
      b.query_begin>b.query_end || b.query_end>q)
    throw CentromereInputError("direction region outside sequence");
}
}

DirectionChain ScoreDirectionChain(std::span<const DirectionAnchor> anchors,
    const DirectionBox& box, Length reference_length, Length query_length,
    Strand direction, const std::function<void()>& interruption) {
  ValidateBox(box,reference_length,query_length);
  if (direction!=Strand::Forward && direction!=Strand::Reverse)
    throw CentromereInputError("invalid direction");
  const auto qb=direction==Strand::Forward ? box.query_begin : query_length-box.query_end;
  const auto qe=direction==Strand::Forward ? box.query_end : query_length-box.query_begin;
  // Signed diagonals without narrowing uint64 coordinates to int64.
  std::map<std::pair<bool,Length>,std::vector<Piece>> groups;
  for (const auto& a:anchors) {
    if (!a.length || a.reference_begin>reference_length ||
        a.length>reference_length-a.reference_begin || a.query_begin>query_length ||
        a.length>query_length-a.query_begin || !std::isfinite(a.source_weight) ||
        a.source_weight<0 || (a.strand!=Strand::Forward && a.strand!=Strand::Reverse))
      throw CentromereInputError("invalid direction anchor");
    if (a.strand!=direction) continue;
    const auto u=direction==Strand::Forward ? a.query_begin : query_length-a.query_begin-a.length;
    const auto rb=a.reference_begin;
    const Length start=std::max(box.reference_begin>rb ? box.reference_begin-rb : 0,
                               qb>u ? qb-u : 0);
    if (rb>=box.reference_end || u>=qe) continue;
    const Length end=std::min({a.length,box.reference_end-rb,qe-u});
    if (end<=start) continue;
    groups[{rb>=u,rb>=u ? rb-u : u-rb}].push_back(
      {rb+start,u+start,end-start,a.source_weight/static_cast<double>(a.length),a.id});
  }
  DirectionChain out;
  for (auto& [diagonal,pieces]:groups) {
    (void)diagonal;
    if (interruption) interruption();
    std::sort(pieces.begin(),pieces.end(),[](const Piece& a,const Piece& b){
      return std::tie(a.r,a.q,a.n,a.id)<std::tie(b.r,b.q,b.n,b.id);
    });
    for (std::size_t i=0;i<pieces.size();) {
      std::size_t end=i+1;
      Position finish=pieces[i].r+pieces[i].n;
      while(end<pieces.size() && pieces[end].r<=finish) {
        finish=std::max(finish,pieces[end].r+pieces[end].n); ++end;
      }
      DirectionChainNode n;
      n.reference_begin=pieces[i].r; n.oriented_query_begin=pieces[i].q;
      n.length=finish-pieces[i].r;
      // Integrate the maximum per-base source weight, rather than summing
      // duplicate matches or changing the primary length objective.
      std::map<Position,std::vector<std::pair<double,bool>>> events;
      for(std::size_t j=i;j<end;++j) {
        const auto& p=pieces[j];
        events[p.r].push_back({p.density,true});
        events[p.r+p.n].push_back({p.density,false});
        n.source_ids.push_back(p.id);
      }
      std::map<double,std::size_t> active;
      Position previous=n.reference_begin;
      for(const auto& [position,changes]:events) {
        if(!active.empty()) n.weight+=static_cast<double>(position-previous)*active.rbegin()->first;
        for(const auto& [density,add]:changes) {
          if(add) ++active[density];
          else {auto it=active.find(density); if(--it->second==0) active.erase(it);}
        }
        previous=position;
      }
      std::sort(n.source_ids.begin(),n.source_ids.end());
      n.source_ids.erase(std::unique(n.source_ids.begin(),n.source_ids.end()),n.source_ids.end());
      out.nodes.push_back(std::move(n)); i=end;
    }
  }
  std::sort(out.nodes.begin(),out.nodes.end(),[](const auto& a,const auto& b){
    return std::tie(a.reference_begin,a.oriented_query_begin,a.length,a.source_ids)<
           std::tie(b.reference_begin,b.oriented_query_begin,b.length,b.source_ids);
  });
  // Sweep completed reference intervals and query a prefix maximum. Only
  // predecessors with the best PRIMARY length can win; enumerate every tied
  // predecessor to preserve the RaMA secondary score and earliest-index tie.
  // This is exact, not a k-predecessor cap or a beam search.
  std::vector<std::size_t> ends(out.nodes.size());
  std::vector<Position> query_ends;
  for(std::size_t i=0;i<ends.size();++i) {
    ends[i]=i;query_ends.push_back(out.nodes[i].oriented_query_begin+out.nodes[i].length);
  }
  std::sort(ends.begin(),ends.end(),[&](auto a,auto b) {
    return std::pair{out.nodes[a].reference_begin+out.nodes[a].length,a}<
           std::pair{out.nodes[b].reference_begin+out.nodes[b].length,b};
  });
  std::sort(query_ends.begin(),query_ends.end());
  query_ends.erase(std::unique(query_ends.begin(),query_ends.end()),query_ends.end());
  struct Frontier {Length support{};std::vector<std::size_t> indices;};
  std::vector<Frontier> tree(query_ends.size()+1);
  std::size_t completed=0;
  for(std::size_t i=0;i<out.nodes.size();++i) {
    if(interruption) interruption();
    auto& current=out.nodes[i];
    current.support=current.length; current.source_score=current.weight;
    while(completed<ends.size() &&
          out.nodes[ends[completed]].reference_begin+out.nodes[ends[completed]].length<=current.reference_begin) {
      const auto j=ends[completed++];const auto& node=out.nodes[j];
      auto p=static_cast<std::size_t>(std::lower_bound(query_ends.begin(),query_ends.end(),
        node.oriented_query_begin+node.length)-query_ends.begin())+1;
      for(;p<tree.size();p+=p&(~p+1)) {
        auto& entry=tree[p];
        if(node.support>entry.support) {entry.support=node.support;entry.indices.clear();}
        if(node.support==entry.support)entry.indices.push_back(j);
      }
    }
    Length best_support=0;std::vector<std::size_t> eligible;
    auto p=static_cast<std::size_t>(std::upper_bound(query_ends.begin(),query_ends.end(),
      current.oriented_query_begin)-query_ends.begin());
    for(;p;p-=p&(~p+1)) {
      const auto& entry=tree[p];
      if(entry.support>best_support) {best_support=entry.support;eligible.clear();}
      if(entry.support==best_support)eligible.insert(eligible.end(),entry.indices.begin(),entry.indices.end());
    }
    std::sort(eligible.begin(),eligible.end());
    for(const auto j:eligible) {
      const auto& prior=out.nodes[j];
      const auto length=prior.support+current.length;
      const auto rg=current.reference_begin-prior.reference_begin;
      const auto qg=current.oriented_query_begin-prior.oriented_query_begin;
      const auto delta=rg>=qg ? rg-qg : qg-rg;
      const double cost=2.0*std::log2(static_cast<double>(delta)+1.0);
      const double score=prior.source_score+(current.weight>cost ? current.weight-cost : 0.1);
      if(length>current.support || (length==current.support && score>current.source_score)) {
        current.support=length; current.source_score=score; current.predecessor=j;
      }
    }
    if(out.best==kNoDirectionNode || current.support>out.support ||
       (current.support==out.support && current.source_score>out.nodes[out.best].source_score)) {
      out.best=i; out.support=current.support;
    }
  }
  return out;
}
} // namespace ramag

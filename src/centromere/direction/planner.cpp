#include "ramag/centromere_blocks.hpp"
#include "ramag/internal/common/centromere.hpp"
#include "ramag/internal/centromere/direction.hpp"
#include "ramag/internal/centromere/evidence_review.hpp"
#include "ramag/internal/centromere/gap_review.hpp"
#include "ramag/internal/centromere/gap_review.hpp"
#include "ramag/internal/centromere/placement_review.hpp"
#include "ramag/internal/centromere/evidence_review.hpp"
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <deque>
#include <limits>
#include <map>
#include <memory_resource>
#include <memory>
#include <numeric>
#include <set>
#include <tuple>

namespace ramag {
namespace centromere_planning_detail {
struct Item {
  Position rb{},re{},qb{},qe{};
  Length support{};
  double weight{};
  std::size_t endpoint{};
  bool atom{};
};
struct Score {
  Length support{};
  double weight{};
  std::size_t blocks{},predecessor{none};
};
struct Solution {
  Score score;
  Vec<std::size_t> path;
  explicit Solution(Budget& budget):path(&budget){}
};
void Check(const CentromereBlockOptions& options) {
  if(options.interruption_callback)options.interruption_callback("centromere-joint-planning");
}
int Compare(const Score& a,const Score& b) {
  if(a.support!=b.support)return a.support>b.support ? 1 : -1;
  if(a.weight!=b.weight)return a.weight>b.weight ? 1 : -1;
  if(a.blocks!=b.blocks)return a.blocks<b.blocks ? 1 : -1;
  return 0;
}
Position Begin(const Item& item,Strand background,Length qlength) {
  return background==Strand::Forward ? item.qb : qlength-item.qe;
}
Position End(const Item& item,Strand background,Length qlength) {
  return background==Strand::Forward ? item.qe : qlength-item.qb;
}
Solution Solve(const Vec<Item>& items,Strand background,Length qlength,
               const CentromereBlockOptions& options,Budget& budget) {
  Vec<std::size_t> starts(items.size(),0,&budget),ends(items.size(),0,&budget);
  std::iota(starts.begin(),starts.end(),0);ends=starts;
  std::sort(starts.begin(),starts.end(),[&](auto a,auto b){
    return std::tuple{items[a].rb,Begin(items[a],background,qlength),a}<
           std::tuple{items[b].rb,Begin(items[b],background,qlength),b};
  });
  std::sort(ends.begin(),ends.end(),[&](auto a,auto b){return std::pair{items[a].re,a}<std::pair{items[b].re,b};});
  Vec<Position> coordinates(&budget);coordinates.reserve(items.size());
  for(const auto& i:items)coordinates.push_back(End(i,background,qlength));
  std::sort(coordinates.begin(),coordinates.end());
  coordinates.erase(std::unique(coordinates.begin(),coordinates.end()),coordinates.end());
  Vec<Score> scores(items.size(),Score{},&budget);
  Vec<std::size_t> tree(coordinates.size()+1,none,&budget);
  auto better=[&](std::size_t a,std::size_t b) {
    if(a==none)return false;
    if(b==none)return true;
    const int cmp=Compare(scores[a],scores[b]);
    return cmp>0 || (cmp==0 && a<b);
  };
  std::size_t completed=0,best=none;
  for(const auto i:starts) {
    Check(options);
    while(completed<ends.size() && items[ends[completed]].re<=items[i].rb) {
      const auto j=ends[completed++];
      if(!scores[j].blocks)throw AlignmentError("joint sweep used an unfinished predecessor");
      auto p=static_cast<std::size_t>(std::lower_bound(coordinates.begin(),coordinates.end(),
        End(items[j],background,qlength))-coordinates.begin())+1;
      for(;p<tree.size();p+=p&(~p+1))if(better(j,tree[p]))tree[p]=j;
    }
    auto p=static_cast<std::size_t>(std::upper_bound(coordinates.begin(),coordinates.end(),
      Begin(items[i],background,qlength))-coordinates.begin());
    std::size_t prior=none;
    for(;p;p-=p&(~p+1))if(better(tree[p],prior))prior=tree[p];
    auto& score=scores[i];
    if(prior!=none)score=scores[prior];
    score.support+=items[i].support;score.weight+=items[i].weight;
    ++score.blocks;score.predecessor=prior;
    if(better(i,best))best=i;
  }
  Solution result(budget);
  if(best!=none) {
    result.score=scores[best];
    for(auto i=best;i!=none;i=scores[i].predecessor)result.path.push_back(i);
    std::reverse(result.path.begin(),result.path.end());
  }
  return result;
}

// The geometry-only public planner has no sequence payload. The production
// route supplies both sequences and checks every derived exact anchor again.
void ValidateOverlapVariants(const CentromereBlockResult& out,std::string_view reference,
    std::string_view query,const CentromereBlockOptions& options) {
  if(reference.empty() && query.empty())return;
  auto dna=[](char c) {
    switch(c){case 'A':case 'a':return 'A';case 'C':case 'c':return 'C';
      case 'G':case 'g':return 'G';case 'T':case 't':return 'T';default:return '\0';}
  };
  auto complement=[](char c){return c=='A'?'T':c=='C'?'G':c=='G'?'C':c=='T'?'A':'\0';};
  for(const auto& v:out.overlap_variants) {
    const auto& a=v.variant;
    for(Length i=0;i<a.length;++i) {
      if((i&65535U)==0)Check(options);
      const auto r=dna(reference[static_cast<std::size_t>(a.reference_begin+i)]);
      const auto qpos=v.strand==Strand::Forward ? a.query_begin+i : a.query_begin+a.length-1-i;
      const auto q=dna(query[static_cast<std::size_t>(qpos)]);
      if(!r || !q || r!=(v.strand==Strand::Forward?q:complement(q)))
        throw AlignmentError("overlap-clipped candidate is not an exact DNA match");
    }
  }
}

// The same length-first DP and RaMA secondary score as ScoreDirectionChain,
// with compact rows and a bounded allocator for every transient container.
// No source-ID/path copies or per-anchor provenance objects are needed here.
void CompactChains(const std::vector<CentromerePlanningAnchor>& input,Strand strand,Length qlength,
                   Vec<Node>& nodes,Budget& budget,const CentromereBlockOptions& options,CentromereBlockResult& out,bool normalize_only=false) {
  struct Piece {Position r,u;Length n;double density;std::size_t id;};
  Vec<Piece> pieces(&budget);pieces.reserve(static_cast<std::size_t>(std::count_if(input.begin(),input.end(),
    [&](const auto& a){return a.strand==strand;})));
  for(std::size_t i=0;i<input.size();++i)if(input[i].strand==strand) {
    const auto& a=input[i].anchor;
    pieces.push_back({a.reference_begin,strand==Strand::Forward ? a.query_begin : qlength-a.query_begin-a.length,
      a.length,a.weight/static_cast<double>(a.length),i});
  }
  auto diagonal=[](const Piece& a){return std::pair{a.r>=a.u,a.r>=a.u ? a.r-a.u : a.u-a.r};};
  std::sort(pieces.begin(),pieces.end(),[&](const auto& a,const auto& b){
    return std::tuple{diagonal(a),a.r,a.u,a.n,a.id}<std::tuple{diagonal(b),b.r,b.u,b.n,b.id};});
  const auto base=nodes.size();
  struct Event {Position position;double density;bool add;};
  for(std::size_t i=0;i<pieces.size();) {
    Check(options);std::size_t end=i+1;Position finish=pieces[i].r+pieces[i].n;
    while(end<pieces.size() && diagonal(pieces[end])==diagonal(pieces[i]) && pieces[end].r<=finish) {
      finish=std::max(finish,pieces[end].r+pieces[end].n);++end;
    }
    Vec<Event> events(&budget);events.reserve(2*(end-i));
    for(auto j=i;j<end;++j){events.push_back({pieces[j].r,pieces[j].density,true});events.push_back({pieces[j].r+pieces[j].n,pieces[j].density,false});}
    std::sort(events.begin(),events.end(),[](const auto& a,const auto& b){return a.position<b.position;});
    std::pmr::map<double,std::size_t> active(&budget);Position prior=pieces[i].r;double weight=0;
    double first_density=0;bool have_density=false,uniform_density=true;
    for(const auto& e:events) {
      if(!active.empty() && e.position>prior) {
        const auto density=active.rbegin()->first;
        weight+=static_cast<double>(e.position-prior)*density;
        if(have_density && density!=first_density)uniform_density=false;
        if(!have_density){first_density=density;have_density=true;}
      }
      if(e.add)++active[e.density];
      else {auto it=active.find(e.density);if(it==active.end())throw AlignmentError("invalid normalization events");if(!--it->second)active.erase(it);}
      prior=e.position;
    }
    Node n;n.anchor={pieces[i].r,strand==Strand::Forward ? pieces[i].u : qlength-pieces[i].u-(finish-pieces[i].r),finish-pieces[i].r,weight};
    n.strand=strand;n.uniform_density=uniform_density;n.clipping.source=nodes.size();
    nodes.push_back(n);i=end;
  }
  if(options.overlap_compatibility) {
    const auto count=nodes.size()-base;
    auto proposed=detail::FindAnchorEndClips(count,[&](auto i)->const CentromereAnchor& {
      return nodes[base+i].anchor;
    },strand,qlength,&budget,[&](){Check(options);});
    out.overlap_pairs_examined+=proposed.pairs_examined;
    out.overlap_local_support_rejections+=proposed.local_support_rejections;
    out.overlap_budget_exhausted=out.overlap_budget_exhausted || proposed.budget_exhausted;
    if(!proposed.budget_exhausted)for(std::size_t i=0;i<count;++i) {
      const auto original=nodes[base+i];const auto clips=proposed.clips[i];
      // Integral normalized weights cannot be prorated when their maximum raw
      // density changes along the anchor. Preserve such original nodes intact.
      if(!original.uniform_density) {
        if(clips.front || clips.back)++out.overlap_nonuniform_nodes_skipped;
        continue;
      }
      auto append=[&](Length front,Length back) {
        if(!front && !back)return;
        if(front>=original.anchor.length || back>=original.anchor.length-front)return;
        Node n=original;n.overlap_variant=true;
        n.clipping.front=front ? nodes[base+clips.front_witness].clipping.source : none;
        n.clipping.back=back ? nodes[base+clips.back_witness].clipping.source : none;
        n.clipping.trim_front=front;n.clipping.trim_back=back;
        n.anchor.reference_begin+=front;
        n.anchor.query_begin+=strand==Strand::Forward ? front : back;
        n.anchor.length-=front+back;
        n.anchor.weight=original.anchor.weight*static_cast<double>(n.anchor.length)/static_cast<double>(original.anchor.length);
        auto supported=[&](std::size_t witness,double cost,bool before) {
          if(witness>=count)return false;
          const auto& a=before ? nodes[base+witness].anchor : n.anchor;
          const auto& b=before ? n.anchor : nodes[base+witness].anchor;
          auto u=[&](const auto& x){return strand==Strand::Forward ? x.query_begin : qlength-x.query_begin-x.length;};
          if(a.reference_begin+a.length>b.reference_begin || u(a)+a.length>u(b))
            throw AlignmentError("clipping witness overlaps after final trimming");
          if(!detail::AnchorClipGapHasLocalSupport(b.reference_begin-a.reference_begin-a.length,
              u(b)-u(a)-a.length,a.weight,b.weight)) {
            ++out.overlap_local_support_rejections;return false;
          }
          return !cost || cost<std::min(a.weight,b.weight);
        };
        // Combining two individually valid end clips can remove enough weight
        // to invalidate either connection; validate the final piece again.
        if((front && !supported(clips.front_witness,clips.front_cost,true)) ||
           (back && !supported(clips.back_witness,clips.back_cost,false)))return;
        nodes.push_back(n);
        if(out.overlap_variants.size()==out.overlap_variants.capacity())
          Reserve(out.overlap_variants,std::max<std::size_t>(1,out.overlap_variants.capacity()*2),budget);
        out.overlap_variants.push_back({original.anchor,n.anchor,strand,front,back});
        ++out.overlap_variants_added;
      };
      append(clips.front,0);append(0,clips.back);
      if(clips.front && clips.back)append(clips.front,clips.back);
    }
  }
  auto oriented=[&](const Node& n){return strand==Strand::Forward ? n.anchor.query_begin : qlength-n.anchor.query_begin-n.anchor.length;};
  std::sort(nodes.begin()+static_cast<std::ptrdiff_t>(base),nodes.end(),[&](const auto& a,const auto& b){
    return std::tuple{a.anchor.reference_begin,oriented(a),a.anchor.length}<std::tuple{b.anchor.reference_begin,oriented(b),b.anchor.length};});
  if(normalize_only)return;
  Vec<std::size_t> ends(&budget);Vec<Position> coordinates(&budget);
  ends.reserve(nodes.size()-base);coordinates.reserve(nodes.size()-base);
  for(auto i=base;i<nodes.size();++i){ends.push_back(i);coordinates.push_back(oriented(nodes[i])+nodes[i].anchor.length);}
  std::sort(ends.begin(),ends.end(),[&](auto a,auto b){return std::pair{nodes[a].anchor.reference_begin+nodes[a].anchor.length,a}<std::pair{nodes[b].anchor.reference_begin+nodes[b].anchor.length,b};});
  std::sort(coordinates.begin(),coordinates.end());coordinates.erase(std::unique(coordinates.begin(),coordinates.end()),coordinates.end());
  struct Frontier {Length support{};Vec<std::size_t> indices;explicit Frontier(Budget& b):indices(&b){}};
  Vec<Frontier> tree(&budget);tree.reserve(coordinates.size()+1);
  for(std::size_t i=0;i<=coordinates.size();++i)tree.emplace_back(budget);
  Vec<std::size_t> eligible(&budget);std::size_t completed=0;
  for(auto i=base;i<nodes.size();++i) {
    Check(options);auto& n=nodes[i];n.dp_support=n.anchor.length;n.dp_score=n.anchor.weight;
    while(completed<ends.size() && nodes[ends[completed]].anchor.reference_begin+nodes[ends[completed]].anchor.length<=n.anchor.reference_begin) {
      const auto j=ends[completed++];const auto& a=nodes[j];
      auto p=static_cast<std::size_t>(std::lower_bound(coordinates.begin(),coordinates.end(),oriented(a)+a.anchor.length)-coordinates.begin())+1;
      for(;p<tree.size();p+=p&(~p+1)) {
        auto& t=tree[p];if(a.dp_support>t.support){t.support=a.dp_support;t.indices.clear();}
        if(a.dp_support==t.support)t.indices.push_back(j);
      }
    }
    Length best=0;eligible.clear();
    auto p=static_cast<std::size_t>(std::upper_bound(coordinates.begin(),coordinates.end(),oriented(n))-coordinates.begin());
    for(;p;p-=p&(~p+1)) {
      const auto& t=tree[p];if(t.support>best){best=t.support;eligible.clear();}
      if(t.support==best)eligible.insert(eligible.end(),t.indices.begin(),t.indices.end());
    }
    std::sort(eligible.begin(),eligible.end());
    for(const auto j:eligible) {
      const auto& a=nodes[j];const auto rg=n.anchor.reference_begin-a.anchor.reference_begin,qg=oriented(n)-oriented(a);
      const auto delta=rg>=qg ? rg-qg : qg-rg;const auto cost=2.0*std::log2(static_cast<double>(delta)+1.0);
      const auto support=a.dp_support+n.anchor.length;
      const auto score=a.dp_score+(n.anchor.weight>cost ? n.anchor.weight-cost : 0.1);
      if(support>n.dp_support || (support==n.dp_support && score>n.dp_score)) {
        n.dp_support=support;n.dp_score=score;n.dp_predecessor=j;
      }
    }
    n.first=i;n.support=n.anchor.length;n.weight=n.anchor.weight;
    if(n.dp_predecessor!=none) {
      const auto& a=nodes[n.dp_predecessor];const auto rg=n.anchor.reference_begin-a.anchor.reference_begin,qg=oriented(n)-oriented(a);
      const auto delta=rg>=qg ? rg-qg : qg-rg;const auto cost=2.0*std::log2(static_cast<double>(delta)+1.0);
      if(!delta || cost<std::min(n.anchor.weight,a.anchor.weight)) {
        n.predecessor=n.dp_predecessor;n.first=a.first;n.support+=a.support;n.weight+=a.weight;
      }
    }
  }
}

// Exact signed-magnitude diagonals avoid narrowing 64-bit coordinates. The
// coarse query radius is conservative; the original coherence test below is
// always checked on both anchors before ANY path competes for a predecessor.
struct Diagonal {
  bool negative{};
  Length magnitude{};
  bool operator<(const Diagonal& b) const {
    if(negative!=b.negative)return negative;
    return negative ? magnitude>b.magnitude : magnitude<b.magnitude;
  }
};
Diagonal Diag(Position r,Position u){return r<u ? Diagonal{true,u-r} : Diagonal{false,r-u};}
Diagonal Shift(Diagonal d,Length n,bool up) {
  if(d.negative==up) {
    if(n>=d.magnitude)return Diagonal{!up,n-d.magnitude};
    return Diagonal{d.negative,d.magnitude-n};
  }
  const auto maximum=std::numeric_limits<Length>::max();
  return Diagonal{d.negative,n>maximum-d.magnitude ? maximum : d.magnitude+n};
}
Length Radius(double weight) {
  if(weight>=128)return std::numeric_limits<Length>::max();
  const auto value=std::ceil(std::exp2(weight/2.0));
  if(value>=static_cast<double>(std::numeric_limits<Length>::max()))return std::numeric_limits<Length>::max();
  return static_cast<Length>(value);
}
struct Range{std::size_t begin{},end{};};
struct LocalGraph {
  Vec<Range> incoming;
  Vec<std::size_t> edges;
  Vec<double> increments;
  Length query_length{};
  bool hierarchical{};
  bool single_predecessor{};
  LocalGraph(std::size_t size,Budget& b):incoming(size,Range{},&b),edges(&b),increments(&b){}
};
void LocalConnections(const Vec<Node>& nodes,std::size_t base,std::size_t finish,Length qlength,
    LocalGraph& graph,Budget& budget,const CentromereBlockOptions& options) {
  if(base==finish)return;
  const auto strand=nodes[base].strand;
  auto u=[&](std::size_t i){const auto& a=nodes[i].anchor;
    return strand==Strand::Forward ? a.query_begin : qlength-a.query_begin-a.length;};
  Vec<std::size_t> ends(&budget);
  for(auto i=base;i<finish;++i)ends.push_back(i);
  std::sort(ends.begin(),ends.end(),[&](auto a,auto b){return std::pair{nodes[a].anchor.reference_begin+nodes[a].anchor.length,a}<
    std::pair{nodes[b].anchor.reference_begin+nodes[b].anchor.length,b};});
  std::pmr::map<Diagonal,Vec<std::size_t>> active(&budget);
  std::size_t completed=0;
  for(auto i=base;i<finish;++i) {
    Check(options);
    if(options.progress_callback && (i-base)%1024==0)
      options.progress_callback("centromere-local-connections",i-base,finish-base);
    const auto& a=nodes[i].anchor;
    while(completed<ends.size() && nodes[ends[completed]].anchor.reference_begin+nodes[ends[completed]].anchor.length<=a.reference_begin) {
      const auto j=ends[completed++];active.try_emplace(Diag(nodes[j].anchor.reference_begin,u(j))).first->second.push_back(j);
    }
    graph.incoming[i].begin=graph.edges.size();
    const auto d=Diag(a.reference_begin,u(i));const auto radius=Radius(a.weight);
    auto low=radius==std::numeric_limits<Length>::max() ? active.begin() : active.lower_bound(Shift(d,radius,false));
    auto high=radius==std::numeric_limits<Length>::max() ? active.end() : active.upper_bound(Shift(d,radius,true));
    for(auto it=low;it!=high;++it)for(const auto j:it->second) {
      const auto& b=nodes[j].anchor;
      if(!detail::AnchorClipConnectionAllowed(nodes[j].clipping,nodes[i].clipping))continue;
      if(u(j)+b.length>u(i))continue;
      if((nodes[j].overlap_variant || nodes[i].overlap_variant) &&
         !detail::AnchorClipGapHasLocalSupport(a.reference_begin-b.reference_begin-b.length,
           u(i)-u(j)-b.length,a.weight,b.weight))continue;
      const auto rg=a.reference_begin-b.reference_begin,qg=u(i)-u(j);
      const auto delta=rg>=qg ? rg-qg : qg-rg;
      const double cost=2.0*std::log2(static_cast<double>(delta)+1.0);
      if(delta && cost>=std::min(a.weight,b.weight))continue;
      graph.edges.push_back(j);
    }
    graph.incoming[i].end=graph.edges.size();
  }
}

double ConnectionCost(const Node& a,const Node& b,Length qlength) {
  const auto au=a.strand==Strand::Forward ? a.anchor.query_begin : qlength-a.anchor.query_begin-a.anchor.length;
  const auto bu=b.strand==Strand::Forward ? b.anchor.query_begin : qlength-b.anchor.query_begin-b.anchor.length;
  const auto r=b.anchor.reference_begin-a.anchor.reference_begin,q=bu-au;
  return 2.0*std::log2(static_cast<double>(r>=q ? r-q : q-r)+1.0);
}
bool NeedsRefinement(const Node& a,const Node& b,Length qlength) {
  const auto cost=ConnectionCost(a,b,qlength);
  return cost>0 && cost>=std::min(a.anchor.weight,b.anchor.weight);
}

struct Tail {Length support{};double weight{};std::size_t next{none},epoch{};};
struct LocalWorkspace {
  Vec<Tail> tails;
  Vec<std::size_t> heap,visited;
  std::size_t epoch{};
  LocalWorkspace(std::size_t n,Budget& b):tails(n,Tail{},&b),heap(&b),visited(&b){}
};
void TailsTo(std::size_t endpoint,const Vec<Node>& nodes,const LocalGraph& graph,
    LocalWorkspace& w,const CentromereBlockOptions& options,CentromereBlockResult& out) {
  ++w.epoch;w.heap.clear();w.visited.clear();
  // A back-clipped node still owes a right-hand witness and cannot end a path.
  if(nodes[endpoint].clipping.back!=none)return;
  w.tails[endpoint]={nodes[endpoint].anchor.length,nodes[endpoint].anchor.weight,none,w.epoch};
  if(graph.single_predecessor) {
    // A predecessor forest has no branching to schedule. Following its links
    // directly gives the same tails and floating-point evaluation order as
    // the generic heap traversal, without repeated heap operations.
    for(auto j=endpoint;;) {
      w.visited.push_back(j);const auto range=graph.incoming[j];
      if(range.begin==range.end)break;
      if(range.end!=range.begin+1)throw AlignmentError("invalid single-predecessor region graph");
      if((++out.local_chain_transitions&65535U)==0)Check(options);
      const auto p=graph.edges[range.begin];const auto& tail=w.tails[j];
      w.tails[p]={tail.support+nodes[p].anchor.length,
        tail.weight+nodes[p].anchor.weight-nodes[j].anchor.weight+graph.increments[range.begin],j,w.epoch};
      j=p;
    }
    out.local_chain_states+=w.visited.size();return;
  }
  w.heap.push_back(endpoint);
  while(!w.heap.empty()) {
    std::pop_heap(w.heap.begin(),w.heap.end());const auto j=w.heap.back();w.heap.pop_back();w.visited.push_back(j);
    const auto range=graph.incoming[j];const auto tail=w.tails[j];
    for(auto k=range.begin;k<range.end;++k) {
      if((++out.local_chain_transitions&65535U)==0)Check(options);
      const auto p=graph.edges[k];auto& t=w.tails[p];
      const auto support=tail.support+nodes[p].anchor.length;
      const auto weight=graph.hierarchical ? tail.weight+nodes[p].anchor.weight-nodes[j].anchor.weight+graph.increments[k] :
        tail.weight+nodes[p].anchor.weight;
      if(t.epoch!=w.epoch) {
        t={support,weight,j,w.epoch};w.heap.push_back(p);std::push_heap(w.heap.begin(),w.heap.end());
      } else if(support>t.support || (support==t.support && (weight>t.weight || (weight==t.weight && j<t.next))))
        t={support,weight,j,w.epoch};
    }
  }
  out.local_chain_states+=w.visited.size();
}
Item LocalItem(std::size_t first,std::size_t endpoint,const Vec<Node>& nodes,Length support,double weight) {
  const auto& a=nodes[first].anchor;const auto& b=nodes[endpoint].anchor;
  return {a.reference_begin,b.reference_begin+b.length,std::min(a.query_begin,b.query_begin),
    std::max(a.query_begin+a.length,b.query_begin+b.length),support,weight,endpoint,first==endpoint};
}
struct SelectedLocal {std::size_t first{},endpoint{};Length support{};double weight{};std::uint32_t family{};};
struct LocalSolution {
  Score score;
  Vec<SelectedLocal> path;
  explicit LocalSolution(Budget& budget):path(&budget){}
};
// Enumerate all start/end envelopes one endpoint at a time. The reverse DAG
// DP gives the best path for EVERY start without retaining an O(n^2) table.
// A persistent prefix index supplies the macro solution as it stood before
// that start. Only globally useful macro improvements are retained. This is
// the same all-path objective, not a top-k or a sampled start frontier.
LocalSolution SolveLocal(const Vec<Node>& nodes,const LocalGraph& graph,Strand background,Length qlength,
    LocalWorkspace& work,Budget& budget,const CentromereBlockOptions& options,CentromereBlockResult& out,
    const LocalGraph* regional=nullptr) {
  using Index=std::uint32_t;constexpr auto empty_index=std::numeric_limits<Index>::max();
  struct Record {Length support{};Index blocks{},predecessor{empty_index},first{},endpoint{},family{};};
  struct Tree {Index left{empty_index},right{empty_index},best{empty_index};};
  const auto index=[&](std::size_t n)->Index {
    if(n==none)return empty_index;
    if(n>=empty_index)throw AlignmentError("local persistent index exceeds its checked representation range");
    return static_cast<Index>(n);
  };
  const auto wide=[&](Index n)->std::size_t{return n==empty_index ? none : n;};
  struct Version {Position end{};std::size_t root{none};};
  // These append-only persistent stores grow in small independent chunks.
  // A vector's old+new reallocation overlap otherwise exhausts the budget
  // even when the live logical index itself still fits.
  std::pmr::deque<Record> records(&budget);std::pmr::deque<Tree> tree(&budget);Vec<Version> versions(&budget);
  // First-hit unique matches commonly have weight == length. In that exact
  // case the cumulative double is redundant; otherwise retain every actual
  // weight in a parallel bounded store, with no numerical quantization.
  const bool length_weight=!graph.hierarchical && qlength<=9007199254740992ULL &&
    std::all_of(nodes.begin(),nodes.end(),[](const auto& n){return n.anchor.weight==static_cast<double>(n.anchor.length);});
  std::pmr::deque<double> weights(&budget);
  Vec<Position> coordinates(&budget),starts(&budget);Vec<std::size_t> ends(nodes.size(),0,&budget);
  for(const auto& n:nodes) {
    // Every possible block end is an anchor q-end (forward background),
    // or q-begin reflected into the reverse background. The other boundary
    // values are query thresholds, never stored coordinates in this index.
    coordinates.push_back(background==Strand::Forward ? n.anchor.query_begin+n.anchor.length : qlength-n.anchor.query_begin);
    starts.push_back(n.anchor.reference_begin);
  }
  std::sort(coordinates.begin(),coordinates.end());coordinates.erase(std::unique(coordinates.begin(),coordinates.end()),coordinates.end());
  std::sort(starts.begin(),starts.end());starts.erase(std::unique(starts.begin(),starts.end()),starts.end());
  std::iota(ends.begin(),ends.end(),0);
  std::sort(ends.begin(),ends.end(),[&](auto a,auto b){return std::pair{nodes[a].anchor.reference_begin+nodes[a].anchor.length,a}<
    std::pair{nodes[b].anchor.reference_begin+nodes[b].anchor.length,b};});
  auto score=[&](std::size_t id){const auto& r=records[id];return Score{r.support,
    length_weight ? static_cast<double>(r.support) : weights[id],r.blocks,wide(r.predecessor)};};
  auto better_candidate=[&](const Score& a,std::size_t first,std::size_t endpoint,std::uint32_t family,std::size_t b){const auto c=Compare(a,score(b));
    return c>0 || (c==0 && std::tuple{first,endpoint,family,a.predecessor}<
      std::tuple{static_cast<std::size_t>(records[b].first),static_cast<std::size_t>(records[b].endpoint),records[b].family,wide(records[b].predecessor)});};
  auto best=[&](std::size_t a,std::size_t b){return a==none ? b : b==none ? a : better_candidate(score(a),records[a].first,records[a].endpoint,records[a].family,b) ? a : b;};
  auto query=[&](auto&& self,std::size_t root,std::size_t lo,std::size_t hi,std::size_t limit)->std::size_t {
    if(root==none || lo>=limit)return none;
    if(hi<=limit)return wide(tree[root].best);
    const auto mid=lo+(hi-lo)/2;
    return best(self(self,wide(tree[root].left),lo,mid,limit),self(self,wide(tree[root].right),mid,hi,limit));
  };
  std::size_t writable=0;
  auto update=[&](auto&& self,std::size_t root,std::size_t lo,std::size_t hi,std::size_t position,std::size_t record)->std::size_t {
    // Only the completed reference-end versions are snapshots. Nodes born
    // during this SAME end event have never been published to a predecessor
    // query and may be updated in place. Cloning every interim leaf update
    // would retain histories that no possible candidate can ever observe.
    auto id=root;
    if(root==none || root<writable) {
      id=tree.size();(void)index(id);tree.push_back(root==none ? Tree{} : tree[root]);
    }
    if(hi-lo==1){tree[id].best=index(best(wide(tree[id].best),record));return id;}
    const auto mid=lo+(hi-lo)/2;
    if(position<mid){const auto child=self(self,wide(tree[id].left),lo,mid,position,record);tree[id].left=index(child);}
    else {const auto child=self(self,wide(tree[id].right),mid,hi,position,record);tree[id].right=index(child);}
    tree[id].best=index(best(tree[id].left==empty_index ? none : wide(tree[tree[id].left].best),
                            tree[id].right==empty_index ? none : wide(tree[tree[id].right].best)));
    return id;
  };
  std::size_t latest=none;
  std::size_t prior_bucket=none,writable_records=0;
  for(std::size_t begin=0;begin<ends.size();) {
    const auto rend=nodes[ends[begin]].anchor.reference_begin+nodes[ends[begin]].anchor.length;
    const auto bucket=static_cast<std::size_t>(std::lower_bound(starts.begin(),starts.end(),rend)-starts.begin());
    // Only anchor START coordinates ever query a historical version. Several
    // end events with no intervening start are observationally identical to
    // every future candidate, so their unobservable snapshots may coalesce.
    if(bucket!=prior_bucket){writable=tree.size();writable_records=records.size();}
    else versions.pop_back();
    auto finish=begin+1;while(finish<ends.size() && nodes[ends[finish]].anchor.reference_begin+nodes[ends[finish]].anchor.length==rend)++finish;
    for(auto e=begin;e<finish;++e) {
      Check(options);const auto endpoint=ends[e];
      for(std::uint32_t family=0;family<(regional ? 2U : 1U);++family) {
      TailsTo(endpoint,nodes,family ? *regional : graph,work,options,out);
      for(const auto first:work.visited) {
        if(nodes[first].clipping.front!=none)continue;
        const auto& tail=work.tails[first];const auto item=LocalItem(first,endpoint,nodes,tail.support,tail.weight);
        const auto ver=std::upper_bound(versions.begin(),versions.end(),item.rb,[](auto p,const auto& v){return p<v.end;});
        const auto root=ver==versions.begin() ? none : std::prev(ver)->root;
        const auto limit=static_cast<std::size_t>(std::upper_bound(coordinates.begin(),coordinates.end(),Begin(item,background,qlength))-coordinates.begin());
        const auto prior=query(query,root,0,coordinates.size(),limit);
        Score candidate;
        if(prior!=none)candidate=score(prior);
        candidate.support+=tail.support;candidate.weight+=tail.weight;++candidate.blocks;candidate.predecessor=prior;
        const auto pos=static_cast<std::size_t>(std::lower_bound(coordinates.begin(),coordinates.end(),End(item,background,qlength))-coordinates.begin());
        const auto existing=query(query,latest,0,coordinates.size(),pos+1);
        if(existing!=none && !better_candidate(candidate,first,endpoint,family,existing)){++out.local_chain_dominated;continue;}
        auto leaf=latest;std::size_t lo=0,hi=coordinates.size();
        while(leaf!=none && hi-lo>1) {
          const auto mid=lo+(hi-lo)/2;
          if(pos<mid){leaf=wide(tree[leaf].left);hi=mid;}else {leaf=wide(tree[leaf].right);lo=mid;}
        }
        auto id=leaf==none ? none : wide(tree[leaf].best);
        const Record record{candidate.support,index(candidate.blocks),index(prior),index(first),index(endpoint),family};
        // A record born at this end event has no published predecessor
        // references yet. Improving the SAME leaf can reuse it, just as its
        // unobserved tree nodes can be updated in place.
        if(id!=none && id>=writable_records) {
          records[id]=record;if(!length_weight)weights[id]=candidate.weight;
        } else {
          id=records.size();(void)index(id);records.push_back(record);if(!length_weight)weights.push_back(candidate.weight);
        }
        latest=update(update,latest,0,coordinates.size(),pos,id);
      }
      }
    }
    versions.push_back({rend,latest});prior_bucket=bucket;begin=finish;
    if(options.progress_callback && begin%1024==0)options.progress_callback("centromere-local-planning",begin,ends.size());
  }
  LocalSolution result(budget);
  if(latest!=none) {
    auto i=wide(tree[latest].best);result.score=score(i);
    for(;i!=none;i=wide(records[i].predecessor)) {
      const auto& r=records[i];const auto current=score(i);const auto prior=wide(r.predecessor);
      result.path.push_back({r.first,r.endpoint,current.support-(prior==none ? 0 : records[prior].support),
        current.weight-(prior==none ? 0 : score(prior).weight),r.family});
    }
    std::reverse(result.path.begin(),result.path.end());
  }
  return result;
}

// v8 adds finite count-first region candidates to the unchanged v6 macro
// objective. It never discovers matches or expands all suffixes of a path.
struct CountBox {
  Position rb{},re{},qb{},qe{};
  auto Key() const {return std::tuple{rb,re,qb,qe};}
  bool operator<(const CountBox& b) const {return Key()<b.Key();}
};
class CountClusters {
  struct State {std::size_t predecessor{none},first{};std::uint64_t count{};Length support{};double weight{};};
  struct Spatial {CountBox box;std::size_t begin{},end{},left{none},right{none};};
  struct Label {Position begin{},end{};unsigned mask{};};
  const Vec<Node>& nodes_;
  Length qlength_;
  const CentromereBlockOptions& options_;
  CentromereBlockResult& out_;
  Budget& budget_;
  std::string_view reference_,query_;
  Vec<State> states_;
  Vec<bool> neutral_;
  Vec<std::size_t> order_,candidate_by_end_;
  Vec<Spatial> spatial_;
  Vec<Label> labels_;
  std::pmr::map<CountBox,std::pair<std::uint64_t,std::uint64_t>> counts_;
  Position U(std::size_t i) const {
    const auto& n=nodes_[i];return n.strand==Strand::Forward ? n.anchor.query_begin :
      qlength_-n.anchor.query_begin-n.anchor.length;
  }
  CountBox Box(std::size_t i) const {
    const auto& a=nodes_[i].anchor;return {a.reference_begin,a.reference_begin+a.length,a.query_begin,a.query_begin+a.length};
  }
  template<class T> void Append(std::vector<T>& target,T value) {
    if(target.size()==target.capacity())Reserve(target,std::max<std::size_t>(1,target.capacity()*2),budget_);
    target.push_back(std::move(value));
  }
  bool Dominates(std::uint64_t own,std::uint64_t other) const {
    return own>=options_.count_min_anchors && other<=own/options_.count_dominance_ratio;
  }
  static unsigned Bit(Strand s){return s==Strand::Forward ? 1U : 2U;}
  bool Blocked(Position begin,Position end,Strand s,bool inference) const {
    if(begin>=end)return false;
    auto it=std::lower_bound(labels_.begin(),labels_.end(),begin,[](const auto& x,auto p){return x.end<=p;});
    for(;it!=labels_.end() && it->begin<end;++it) {
      if(inference) {if((it->mask&(Bit(s)==1 ? 2U : 1U)) || (it->mask&4U))return true;}
      else if(it->mask==(Bit(s)==1 ? 2U : 1U))return true;
    }
    return false;
  }
  std::size_t BuildSpatial(std::size_t begin,std::size_t end,unsigned depth) {
    Spatial row;row.begin=begin;row.end=end;row.box=Box(order_[begin]);
    for(auto k=begin+1;k<end;++k){const auto b=Box(order_[k]);row.box.rb=std::min(row.box.rb,b.rb);
      row.box.re=std::max(row.box.re,b.re);row.box.qb=std::min(row.box.qb,b.qb);row.box.qe=std::max(row.box.qe,b.qe);}
    const auto id=spatial_.size();spatial_.push_back(row);
    if(end-begin>16) {
      const auto middle=begin+(end-begin)/2;
      std::nth_element(order_.begin()+static_cast<std::ptrdiff_t>(begin),order_.begin()+static_cast<std::ptrdiff_t>(middle),
        order_.begin()+static_cast<std::ptrdiff_t>(end),[&](auto a,auto b) {
          const auto& x=nodes_[a].anchor;const auto& y=nodes_[b].anchor;
          return depth%2 ? std::pair{x.query_begin+x.length/2,a}<std::pair{y.query_begin+y.length/2,b} :
            std::pair{x.reference_begin+x.length/2,a}<std::pair{y.reference_begin+y.length/2,b};});
      const auto left=BuildSpatial(begin,middle,depth+1),right=BuildSpatial(middle,end,depth+1);
      spatial_[id].left=left;spatial_[id].right=right;
    }
    return id;
  }
  void Within(std::size_t tree,const CountBox& b,Vec<std::size_t>& result) const {
    if(tree==none)return;
    const auto& n=spatial_[tree];const auto& x=n.box;
    if(x.re<=b.rb || x.rb>=b.re || x.qe<=b.qb || x.qb>=b.qe)return;
    if(n.left!=none){Within(n.left,b,result);Within(n.right,b,result);return;}
    for(auto k=n.begin;k<n.end;++k) {
      const auto id=order_[k];const auto a=Box(id);
      if(!neutral_[id] && a.rb>=b.rb && a.re<=b.re && a.qb>=b.qb && a.qe<=b.qe)result.push_back(id);
    }
  }
  // Exact maximum cardinality nonoverlapping chain. No gap radius is imposed
  // on this competition audit; different copy chains are never added together.
  std::uint64_t MaxCount(const Vec<std::size_t>& ids,Strand sign) {
    Vec<std::size_t> list(&budget_);
    for(auto i:ids)if(nodes_[i].strand==sign)list.push_back(i);
    std::sort(list.begin(),list.end(),[&](auto a,auto b){return std::tuple{nodes_[a].anchor.reference_begin,U(a),a}<
      std::tuple{nodes_[b].anchor.reference_begin,U(b),b};});
    Vec<std::size_t> ends(list.size(),0,&budget_);Vec<Position> coords(&budget_);
    std::iota(ends.begin(),ends.end(),0);
    for(auto i:list)coords.push_back(U(i)+nodes_[i].anchor.length);
    std::sort(coords.begin(),coords.end());coords.erase(std::unique(coords.begin(),coords.end()),coords.end());
    std::sort(ends.begin(),ends.end(),[&](auto a,auto b){return std::pair{Box(list[a]).re,a}<std::pair{Box(list[b]).re,b};});
    Vec<std::uint64_t> score(list.size(),0,&budget_),tree(coords.size()+1,0,&budget_);
    std::size_t completed=0;std::uint64_t best=0;
    for(std::size_t i=0;i<list.size();++i) {
      while(completed<ends.size() && Box(list[ends[completed]]).re<=Box(list[i]).rb) {
        const auto j=ends[completed++];auto p=static_cast<std::size_t>(std::lower_bound(coords.begin(),coords.end(),
          U(list[j])+nodes_[list[j]].anchor.length)-coords.begin())+1;
        if(!score[j])throw AlignmentError("count sweep consumed an unfinished predecessor");
        for(;p<tree.size();p+=p&(~p+1))tree[p]=std::max(tree[p],score[j]);
      }
      auto p=static_cast<std::size_t>(std::upper_bound(coords.begin(),coords.end(),U(list[i]))-coords.begin());
      std::uint64_t previous=0;for(;p;p-=p&(~p+1))previous=std::max(previous,tree[p]);
      score[i]=previous+1;best=std::max(best,score[i]);
    }
    return best;
  }
  std::pair<std::uint64_t,std::uint64_t> Counts(const CountBox& box) {
    const auto found=counts_.find(box);if(found!=counts_.end())return found->second;
    Check(options_);++out_.count_rectangle_queries;
    Vec<std::size_t> ids(&budget_);Within(spatial_.empty() ? none : 0,box,ids);
    out_.count_rectangle_nodes+=ids.size();
    const auto result=std::pair{MaxCount(ids,Strand::Forward),MaxCount(ids,Strand::Reverse)};
    counts_.emplace(box,result);return result;
  }
  void Windows() {
    Vec<std::size_t> midpoints(&budget_);
    for(std::size_t i=0;i<nodes_.size();++i)if(!neutral_[i])midpoints.push_back(i);
    const auto middle=[&](auto i){const auto& a=nodes_[i].anchor;return a.query_begin+a.length/2;};
    std::sort(midpoints.begin(),midpoints.end(),[&](auto a,auto b){return std::pair{middle(a),a}<std::pair{middle(b),b};});
    const auto width=std::min(qlength_,options_.count_window),last=qlength_-width;
    Vec<Position> starts(&budget_);
    for(Position p=0;;) {
      starts.push_back(p);starts.push_back(last-p);
      if(last-p<options_.count_step)break;
      p+=options_.count_step;
    }
    std::sort(starts.begin(),starts.end());starts.erase(std::unique(starts.begin(),starts.end()),starts.end());
    for(auto begin:starts) {
      Check(options_);CentromereCountWindow w;w.id=out_.count_windows.size()+1;w.begin=begin;w.end=begin+width;
      auto lo=std::lower_bound(midpoints.begin(),midpoints.end(),begin,[&](auto i,auto p){return middle(i)<p;});
      auto hi=std::upper_bound(midpoints.begin(),midpoints.end(),w.end,[&](auto p,auto i){return p<middle(i);});
      for(auto sign:{Strand::Forward,Strand::Reverse}) {
        Vec<std::size_t> intervals(&budget_);
        for(auto it=lo;it!=hi;++it) {
          const auto& a=nodes_[*it].anchor;
          // An even-length match has two central bases. Considering either
          // central base makes membership exactly symmetric under whole RC.
          if(nodes_[*it].strand==sign && a.query_begin+(a.length-1)/2<w.end)intervals.push_back(*it);
        }
        std::sort(intervals.begin(),intervals.end(),[&](auto a,auto b){return std::tuple{Box(a).qe,Box(a).qb,a}<std::tuple{Box(b).qe,Box(b).qb,b};});
        Position end=0;std::uint64_t count=0;
        for(auto i:intervals)if(Box(i).qb>=end){++count;end=Box(i).qe;}
        (sign==Strand::Forward ? w.forward : w.reverse)=count;
      }
      w.label=Dominates(w.forward,w.reverse) ? '+' : Dominates(w.reverse,w.forward) ? '-' :
        std::max(w.forward,w.reverse)>=options_.count_min_anchors ? '?' : '.';
      Append(out_.count_windows,w);
    }
    struct Event {Position p;unsigned label;int increment;};
    Vec<Event> events(&budget_);
    for(const auto& w:out_.count_windows)if(w.label!='.') {
      const unsigned label=w.label=='+' ? 0U : w.label=='-' ? 1U : 2U;
      events.push_back({w.begin,label,1});events.push_back({w.end,label,-1});
    }
    std::sort(events.begin(),events.end(),[](const auto& a,const auto& b){return a.p<b.p;});
    std::array<std::int64_t,3> active{};Position previous=0;
    for(std::size_t i=0;i<events.size();) {
      const auto p=events[i].p;
      const unsigned mask=(active[0]>0 ? 1U:0U)|(active[1]>0 ? 2U:0U)|(active[2]>0 ? 4U:0U);
      if(previous<p && mask)labels_.push_back({previous,p,mask});
      while(i<events.size() && events[i].p==p){active[events[i].label]+=events[i].increment;++i;}
      previous=p;
    }
  }
  bool Better(const State& a,const State& b,std::size_t endpoint,std::size_t prior) const {
    if(a.count!=b.count)return a.count>b.count;
    if(a.support!=b.support)return a.support>b.support;
    if(a.weight!=b.weight)return a.weight>b.weight;
    const auto end=Box(endpoint);
    const auto span=[&](const State& s) {
      const auto start=Box(s.first);
      const auto rs=end.re-start.rb,qs=U(endpoint)+nodes_[endpoint].anchor.length-U(s.first);
      // Compare a 65-bit sum without overflow.
      const auto sum=rs+qs;return std::pair{sum<rs,sum};
    };
    if(span(a)!=span(b))return span(a)<span(b);
    return std::pair{a.first,prior}<std::pair{b.first,b.predecessor};
  }
  void Chains(std::size_t base,std::size_t end) {
    Vec<std::size_t> ended(&budget_);for(auto i=base;i<end;++i)if(!neutral_[i])ended.push_back(i);
    std::sort(ended.begin(),ended.end(),[&](auto a,auto b){return std::pair{Box(a).re,a}<std::pair{Box(b).re,b};});
    std::pmr::set<std::pair<Position,std::size_t>> active(&budget_);
    std::size_t inserted=0,expired=0;
    for(auto i=base;i<end;++i) {
      if(neutral_[i])continue;
      Check(options_);const auto box=Box(i);const auto u=U(i);
      if(options_.progress_callback && (i-base)%8192==0)
        options_.progress_callback("centromere-count-clustering",i-base,end-base);
      const auto lower=box.rb>options_.count_cluster_gap ? box.rb-options_.count_cluster_gap : 0;
      while(inserted<ended.size() && Box(ended[inserted]).re<=box.rb) {
        const auto j=ended[inserted++];active.emplace(U(j)+nodes_[j].anchor.length,j);
      }
      while(expired<inserted && Box(ended[expired]).re<lower) {
        const auto j=ended[expired++];active.erase({U(j)+nodes_[j].anchor.length,j});
      }
      auto& best=states_[i];best={none,i,1,nodes_[i].anchor.length,nodes_[i].anchor.weight};
      const auto qlower=u>options_.count_cluster_gap ? u-options_.count_cluster_gap : 0;
      for(auto it=active.lower_bound({qlower,0});it!=active.end() && it->first<=u;++it) {
        if((++out_.count_neighbor_checks&65535U)==0)Check(options_);
        const auto j=it->second;const auto other=Box(j);
        const auto qb=nodes_[i].strand==Strand::Forward ? other.qe : box.qe;
        const auto qe=nodes_[i].strand==Strand::Forward ? box.qb : other.qb;
        if(Blocked(qb,qe,nodes_[i].strand,false))continue;
        const auto& p=states_[j];
        State candidate{j,p.first,p.count+1,p.support+nodes_[i].anchor.length,p.weight+nodes_[i].anchor.weight};
        if(Better(candidate,best,i,j))best=candidate;
      }
      if(best.count<options_.count_min_anchors)continue;
      const auto item=LocalItem(best.first,i,nodes_,best.support,best.weight);
      CentromereCountCandidate c;c.id=out_.count_candidates.size()+1;c.first=best.first+1;c.endpoint=i+1;
      c.count=best.count;
      c.rb=item.rb;c.re=item.re;c.qb=item.qb;c.qe=item.qe;c.strand=nodes_[i].strand;c.support=best.support;c.weight=best.weight;
      Append(out_.count_candidates,c);
    }
  }
  // Rectangles sharing a first anchor have the same lower corner in their
  // own orientation. Answer all their exact cardinality queries in one sweep;
  // this is the same unrestricted count as MaxCount, not the proposed path's
  // count or a lower bound. Opposite-copy competition is still counted exactly.
  void AuditRegions() {
    Vec<std::size_t> groups(out_.count_candidates.size(),0,&budget_);
    std::iota(groups.begin(),groups.end(),0);
    std::sort(groups.begin(),groups.end(),[&](auto a,auto b) {
      const auto& x=out_.count_candidates[a];const auto& y=out_.count_candidates[b];
      return std::tuple{x.first,x.re,x.qb,x.qe,a}<std::tuple{y.first,y.re,y.qb,y.qe,b};
    });
    for(std::size_t begin=0;begin<groups.size();) {
      Check(options_);auto end=begin+1;const auto& first=out_.count_candidates[groups[begin]];
      CountBox box{first.rb,first.re,first.qb,first.qe};const auto sign=first.strand;
      while(end<groups.size() && out_.count_candidates[groups[end]].first==first.first) {
        const auto& c=out_.count_candidates[groups[end++]];
        box.re=std::max(box.re,c.re);box.qb=std::min(box.qb,c.qb);box.qe=std::max(box.qe,c.qe);
      }
      Vec<std::size_t> all(&budget_),list(&budget_);Within(spatial_.empty() ? none : 0,box,all);
      for(auto i:all)if(nodes_[i].strand==sign)list.push_back(i);
      out_.count_rectangle_nodes+=list.size();
      std::sort(list.begin(),list.end(),[&](auto a,auto b){return std::tuple{Box(a).rb,U(a),a}<std::tuple{Box(b).rb,U(b),b};});
      Vec<std::size_t> ends(list.size(),0,&budget_);std::iota(ends.begin(),ends.end(),0);
      std::sort(ends.begin(),ends.end(),[&](auto a,auto b){return std::pair{Box(list[a]).re,a}<std::pair{Box(list[b]).re,b};});
      Vec<Position> coords(&budget_);for(auto i:list)coords.push_back(U(i)+nodes_[i].anchor.length);
      std::sort(coords.begin(),coords.end());coords.erase(std::unique(coords.begin(),coords.end()),coords.end());
      Vec<std::uint64_t> score(list.size(),0,&budget_),tree(coords.size()+1,0,&budget_);
      std::size_t completed=0,started=0;
      auto query=[&](Position u) {
        auto p=static_cast<std::size_t>(std::upper_bound(coords.begin(),coords.end(),u)-coords.begin());
        std::uint64_t best=0;for(;p;p-=p&(~p+1))best=std::max(best,tree[p]);return best;
      };
      auto promote=[&](Position r) {
        while(completed<ends.size() && Box(list[ends[completed]]).re<=r) {
          const auto j=ends[completed++];
          if(!score[j])throw AlignmentError("batched count consumed an unfinished predecessor");
          auto p=static_cast<std::size_t>(std::lower_bound(coords.begin(),coords.end(),
            U(list[j])+nodes_[list[j]].anchor.length)-coords.begin())+1;
          for(;p<tree.size();p+=p&(~p+1))tree[p]=std::max(tree[p],score[j]);
        }
      };
      for(auto k=begin;k<end;++k) {
        auto& c=out_.count_candidates[groups[k]];
        while(started<list.size() && Box(list[started]).rb<c.re) {
          promote(Box(list[started]).rb);score[started]=query(U(list[started]))+1;++started;
        }
        promote(c.re);const auto count=query(sign==Strand::Forward ? c.qe : qlength_-c.qb);
        (sign==Strand::Forward ? c.forward : c.reverse)=count;
      }
      begin=end;
    }
    for(std::size_t i=0;i<out_.count_candidates.size();++i) {
      Check(options_);auto& c=out_.count_candidates[i];const CountBox box{c.rb,c.re,c.qb,c.qe};
      const auto found=counts_.find(box);
      if(found!=counts_.end()){c.forward=found->second.first;c.reverse=found->second.second;}
      else {
        ++out_.count_rectangle_queries;Vec<std::size_t> ids(&budget_);Within(spatial_.empty() ? none : 0,box,ids);
        out_.count_rectangle_nodes+=ids.size();
        if(c.strand==Strand::Forward)c.reverse=MaxCount(ids,Strand::Reverse);
        else c.forward=MaxCount(ids,Strand::Forward);
        counts_.emplace(box,std::pair{c.forward,c.reverse});
      }
      c.eligible=c.strand==Strand::Forward ? Dominates(c.forward,c.reverse) : Dominates(c.reverse,c.forward);
      if(c.eligible)candidate_by_end_[static_cast<std::size_t>(c.endpoint-1)]=i;
      if(options_.progress_callback && i%8192==0)
        options_.progress_callback("centromere-count-audit",i,out_.count_candidates.size());
    }
  }
  static bool CanonicalGap(std::string_view sequence,Position b,Position e) {
    return e<=sequence.size() && b<=e && sequence.substr(static_cast<std::size_t>(b),static_cast<std::size_t>(e-b)).find_first_not_of("ACGT")==std::string_view::npos;
  }
  bool Infer(std::size_t left,std::size_t right,const Item& region,std::uint64_t region_id) {
    const auto a=Box(left),b=Box(right);const auto sign=nodes_[left].strand;
    CentromereCountGap gap;gap.region_id=region_id;gap.left_node=left+1;gap.right_node=right+1;gap.strand=sign;
    gap.rb=a.re;gap.re=b.rb;gap.qb=sign==Strand::Forward ? a.qe : b.qe;gap.qe=sign==Strand::Forward ? b.qb : a.qb;
    auto finish=[&](const char* reason,bool accepted=false){gap.reason=reason;gap.accepted=accepted;Append(out_.count_gaps,gap);return accepted;};
    if(gap.rb==gap.re || gap.qb==gap.qe)return finish("one_axis_empty");
    if(gap.re-gap.rb>options_.count_inferred_gap || gap.qe-gap.qb>options_.count_inferred_gap)return finish("gap_limit");
    if(reference_.empty() || query_.empty())return finish("sequence_unavailable");
    if(!CanonicalGap(reference_,gap.rb,gap.re) || !CanonicalGap(query_,gap.qb,gap.qe))return finish("noncanonical_gap");
    if(Blocked(gap.qb,gap.qe,sign,true))return finish("window_competition");
    const char label=sign==Strand::Forward ? '+' : '-';
    bool covering=false,qualified=false;
    for(const auto& w:out_.count_windows) {
      if(w.begin>gap.qb)break;
      if(w.label!=label || w.end<gap.qe)continue;
      covering=true;
      // Query context is clipped to this observation window; reference context
      // remains restricted to the selected copy, not all reference repeats.
      const auto votes=Counts({region.rb,region.re,std::max(region.qb,w.begin),std::min(region.qe,w.end)});
      if(!(sign==Strand::Forward ? Dominates(votes.first,votes.second) : Dominates(votes.second,votes.first)))continue;
      gap.window_id=w.id;gap.forward=votes.first;gap.reverse=votes.second;
      if((sign==Strand::Forward ? votes.second : votes.first)>=options_.count_min_anchors)continue;
      qualified=true;break;
    }
    if(!covering)return finish("no_covering_dominant_window");
    if(!qualified)return finish("local_count_or_opposite_chain");
    return finish("count_inferred_gap",true);
  }
public:
  CountClusters(const Vec<Node>& nodes,std::size_t reverse_base,Length qlength,const CentromereBlockOptions& options,
      CentromereBlockResult& out,Budget& budget,std::string_view reference,std::string_view query):
      nodes_(nodes),qlength_(qlength),options_(options),out_(out),budget_(budget),reference_(reference),query_(query),
      states_(nodes.size(),State{},&budget),neutral_(nodes.size(),false,&budget),order_(&budget),
      candidate_by_end_(nodes.size(),none,&budget),spatial_(&budget),labels_(&budget),counts_(&budget) {
    if(!options.count_window || !options.count_step || options.count_step>options.count_window ||
       options.count_min_anchors<2 || options.count_dominance_ratio<2 || !options.count_cluster_gap ||
       !options.count_inferred_gap || options.count_inferred_gap>options.count_cluster_gap)
      throw CentromereInputError("invalid count-cluster parameters");
    if((!reference.empty() || !query.empty()) && (query.size()!=qlength || reference.empty()))
      throw CentromereInputError("count-cluster sequence context is incomplete");
    const auto begin=std::chrono::steady_clock::now();
    Vec<std::size_t> boxes(nodes.size(),0,&budget);std::iota(boxes.begin(),boxes.end(),0);
    std::sort(boxes.begin(),boxes.end(),[&](auto a,auto b){return std::pair{Box(a).Key(),a}<std::pair{Box(b).Key(),b};});
    for(std::size_t i=0;i<boxes.size();) {
      auto j=i;unsigned signs=0;
      while(j<boxes.size() && Box(boxes[i]).Key()==Box(boxes[j]).Key()){signs|=Bit(nodes_[boxes[j]].strand);++j;}
      if(signs==3)for(auto k=i;k<j;++k)neutral_[boxes[k]]=true;
      i=j;
    }
    Windows();order_.resize(nodes.size());std::iota(order_.begin(),order_.end(),0);
    if(!nodes.empty())BuildSpatial(0,nodes.size(),0);
    Chains(0,reverse_base);Chains(reverse_base,nodes.size());AuditRegions();
    Reserve(out_.count_nodes,nodes.size(),budget_);
    for(std::size_t i=0;i<nodes.size();++i) {
      const auto& s=states_[i];out_.count_nodes.push_back({nodes[i].anchor,nodes[i].strand,
        s.predecessor==none ? 0 : s.predecessor+1,s.count ? s.first+1 : 0,s.count,s.support,s.weight,neutral_[i]});
    }
    out_.count_prepare_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();
  }
  const CentromereCountCandidate* Candidate(std::size_t endpoint) const {
    const auto id=candidate_by_end_[endpoint];return id==none ? nullptr : &out_.count_candidates[id];
  }
  void Emit(std::size_t endpoint,bool inversion) {
    const auto candidate=candidate_by_end_[endpoint];
    if(candidate==none)throw AlignmentError("selected count region is not eligible");
    out_.count_candidates[candidate].selected=true;
    const auto& state=states_[endpoint];
    const auto item=LocalItem(state.first,endpoint,nodes_,state.support,state.weight);
    Vec<std::size_t> path(&budget_);
    for(auto j=endpoint;j!=none;j=states_[j].predecessor)path.push_back(j);
    std::reverse(path.begin(),path.end());
    const auto sign=nodes_[endpoint].strand;const auto region_id=out_.regions.size()+1;
    auto supported=[&](std::size_t first,std::size_t last) {
      const auto a=Box(path[first]),z=Box(path[last]);
      Append(out_.supported_intervals,CentromereSupportInterval{region_id,a.rb,z.re,std::min(a.qb,z.qb),std::max(a.qe,z.qe),sign});
    };
    auto known=[&](std::size_t first,std::size_t last,bool inferred) {
      const auto a=Box(path[first]),z=Box(path[last]);CentromereBlock b;
      b.id=out_.blocks.size()+1;b.region_id=region_id;b.reference_begin=a.rb;b.reference_end=z.re;
      b.query_begin=std::min(a.qb,z.qb);b.query_end=std::max(a.qe,z.qe);b.strand=sign;
      b.skeleton=first==last;b.tentative=inferred;b.inversion_candidate=inversion;
      b.reason=inferred ? "count_inferred_gap" : "count_supported_component";
      Reserve(b.planning_anchors,last-first+1,budget_);
      for(auto k=first;k<=last;++k){b.planning_anchors.push_back(nodes_[path[k]].anchor);
        (sign==Strand::Forward ? b.forward_support : b.reverse_support)+=nodes_[path[k]].anchor.length;}
      Append(out_.blocks,std::move(b));
    };
    std::size_t start=0,support_start=0;bool has_inferred=false;const auto first_block=out_.blocks.size();
    for(std::size_t k=0;k+1<path.size();++k)if(NeedsRefinement(nodes_[path[k]],nodes_[path[k+1]],qlength_)) {
      supported(support_start,k);support_start=k+1;
      if(Infer(path[k],path[k+1],item,region_id)){has_inferred=true;continue;}
      known(start,k,has_inferred);start=k+1;has_inferred=false;
      const auto a=Box(path[k]),b=Box(path[k+1]);CentromereBlock gap;
      gap.id=out_.blocks.size()+1;gap.region_id=region_id;gap.reference_begin=a.re;gap.reference_end=b.rb;
      gap.query_begin=sign==Strand::Forward ? a.qe : b.qe;gap.query_end=sign==Strand::Forward ? b.qb : a.qb;
      gap.strand=sign;gap.direction_known=false;gap.tentative=true;gap.status="unresolved";gap.reason="count_gap_unresolved";
      if(gap.reference_begin!=gap.reference_end || gap.query_begin!=gap.query_end)Append(out_.blocks,std::move(gap));
    }
    supported(support_start,path.size()-1);known(start,path.size()-1,has_inferred);
    CentromereDirectionRegion region;region.id=region_id;region.reference_begin=item.rb;region.reference_end=item.re;
    region.query_begin=item.qb;region.query_end=item.qe;region.strand=sign;region.support=item.support;region.secondary_score=item.weight;
    region.block_count=out_.blocks.size()-first_block;region.left_boundary_width=nodes_[path.front()].anchor.length;
    region.right_boundary_width=nodes_[path.back()].anchor.length;
    region.competing_mappings=sign==Strand::Forward ? out_.count_candidates[candidate].reverse : out_.count_candidates[candidate].forward;
    Append(out_.regions,region);
  }
  void ExistingSupport(const CentromereBlock& b) {
    if(b.direction_known)Append(out_.supported_intervals,CentromereSupportInterval{b.region_id,
      b.reference_begin,b.reference_end,b.query_begin,b.query_end,b.strand});
  }
};
// Experimental region-only planner. The two background orders have separate
// persistent macro states, but share each same-direction path traversal.
// No start/end path is dropped and no quadratic table of paths is retained.
class RegionMacroSweep {
  using Index=std::uint32_t;
  static constexpr auto empty_index=std::numeric_limits<Index>::max();
  struct Record {Length support{};Index blocks{},predecessor{empty_index},first{},endpoint{},family{};double weight{};};
  struct Tree {Index left{empty_index},right{empty_index},best{empty_index};};
  struct Version {Position end{};std::size_t root{none};};
  Strand background_;
  Length qlength_;
  std::pmr::deque<Record> records_;
  std::pmr::deque<Tree> tree_;
  Vec<Version> versions_;
  Vec<Position> coordinates_;
  std::size_t latest_{none},writable_{},writable_records_{},prior_bucket_{none};
  static Index Narrow(std::size_t value) {
    if(value==none)return empty_index;
    if(value>=empty_index)throw AlignmentError("region macro index exceeds its checked representation range");
    return static_cast<Index>(value);
  }
  static std::size_t Wide(Index value){return value==empty_index ? none : value;}
  Score Value(std::size_t id) const {
    const auto& r=records_[id];return {r.support,r.weight,r.blocks,Wide(r.predecessor)};
  }
  bool Better(const Score& a,std::size_t first,std::size_t endpoint,std::uint32_t family,std::size_t b) const {
    const auto c=Compare(a,Value(b));const auto& r=records_[b];
    return c>0 || (c==0 && std::tuple{first,endpoint,family,a.predecessor}<
      std::tuple{static_cast<std::size_t>(r.first),static_cast<std::size_t>(r.endpoint),r.family,Wide(r.predecessor)});
  }
  std::size_t Best(std::size_t a,std::size_t b) const {
    return a==none ? b : b==none ? a : Better(Value(a),records_[a].first,records_[a].endpoint,records_[a].family,b) ? a : b;
  }
  std::size_t Query(std::size_t root,std::size_t lo,std::size_t hi,std::size_t limit) const {
    if(root==none || lo>=limit)return none;
    if(hi<=limit)return Wide(tree_[root].best);
    const auto mid=lo+(hi-lo)/2;
    return Best(Query(Wide(tree_[root].left),lo,mid,limit),Query(Wide(tree_[root].right),mid,hi,limit));
  }
  std::size_t Update(std::size_t root,std::size_t lo,std::size_t hi,std::size_t position,std::size_t record) {
    auto id=root;
    if(root==none || root<writable_){id=tree_.size();(void)Narrow(id);tree_.push_back(root==none ? Tree{} : tree_[root]);}
    if(hi-lo==1){tree_[id].best=Narrow(Best(Wide(tree_[id].best),record));return id;}
    const auto mid=lo+(hi-lo)/2;
    if(position<mid){const auto child=Update(Wide(tree_[id].left),lo,mid,position,record);tree_[id].left=Narrow(child);}
    else {const auto child=Update(Wide(tree_[id].right),mid,hi,position,record);tree_[id].right=Narrow(child);}
    tree_[id].best=Narrow(Best(tree_[id].left==empty_index ? none : Wide(tree_[Wide(tree_[id].left)].best),
                              tree_[id].right==empty_index ? none : Wide(tree_[Wide(tree_[id].right)].best)));
    return id;
  }
public:
  RegionMacroSweep(const Vec<Node>& nodes,Strand background,Length qlength,Budget& budget):
      background_(background),qlength_(qlength),records_(&budget),tree_(&budget),versions_(&budget),coordinates_(&budget) {
    for(const auto& n:nodes)coordinates_.push_back(background==Strand::Forward ?
      n.anchor.query_begin+n.anchor.length : qlength-n.anchor.query_begin);
    std::sort(coordinates_.begin(),coordinates_.end());
    coordinates_.erase(std::unique(coordinates_.begin(),coordinates_.end()),coordinates_.end());
  }
  void BeginEvent(std::size_t bucket) {
    if(bucket!=prior_bucket_){writable_=tree_.size();writable_records_=records_.size();}
    else versions_.pop_back();
    prior_bucket_=bucket;
  }
  void Add(const Item& item,std::size_t first,std::size_t endpoint,std::uint32_t family,CentromereBlockResult& out) {
    const auto version=std::upper_bound(versions_.begin(),versions_.end(),item.rb,[](auto p,const auto& v){return p<v.end;});
    const auto root=version==versions_.begin() ? none : std::prev(version)->root;
    const auto limit=static_cast<std::size_t>(std::upper_bound(coordinates_.begin(),coordinates_.end(),Begin(item,background_,qlength_))-coordinates_.begin());
    const auto prior=Query(root,0,coordinates_.size(),limit);
    Score candidate;if(prior!=none)candidate=Value(prior);
    candidate.support+=item.support;candidate.weight+=item.weight;++candidate.blocks;candidate.predecessor=prior;
    const auto position=static_cast<std::size_t>(std::lower_bound(coordinates_.begin(),coordinates_.end(),End(item,background_,qlength_))-coordinates_.begin());
    const auto existing=Query(latest_,0,coordinates_.size(),position+1);
    if(existing!=none && !Better(candidate,first,endpoint,family,existing)){++out.local_chain_dominated;return;}
    auto leaf=latest_;std::size_t lo=0,hi=coordinates_.size();
    while(leaf!=none && hi-lo>1){const auto mid=lo+(hi-lo)/2;
      if(position<mid){leaf=Wide(tree_[leaf].left);hi=mid;}else{leaf=Wide(tree_[leaf].right);lo=mid;}}
    auto id=leaf==none ? none : Wide(tree_[leaf].best);
    const Record record{candidate.support,Narrow(candidate.blocks),Narrow(prior),Narrow(first),Narrow(endpoint),family,candidate.weight};
    if(id!=none && id>=writable_records_)records_[id]=record;
    else{id=records_.size();(void)Narrow(id);records_.push_back(record);}
    latest_=Update(latest_,0,coordinates_.size(),position,id);
  }
  void EndEvent(Position end){versions_.push_back({end,latest_});}
  LocalSolution Result(Budget& budget) const {
    LocalSolution result(budget);
    if(latest_!=none){auto id=Wide(tree_[latest_].best);result.score=Value(id);
      for(;id!=none;id=Wide(records_[id].predecessor)){
        const auto& r=records_[id];const auto prior=Wide(r.predecessor);
        result.path.push_back({r.first,r.endpoint,r.support-(prior==none ? 0 : records_[prior].support),
          r.weight-(prior==none ? 0 : records_[prior].weight),r.family});}
      std::reverse(result.path.begin(),result.path.end());}
    return result;
  }
};
std::pair<LocalSolution,LocalSolution> SolveLocalBoth(const Vec<Node>& nodes,const LocalGraph& graph,const LocalGraph* regional,
    Length qlength,LocalWorkspace& work,Budget& budget,const CentromereBlockOptions& options,CentromereBlockResult& out,
    const CountClusters* clusters=nullptr) {
  RegionMacroSweep forward(nodes,Strand::Forward,qlength,budget),reverse(nodes,Strand::Reverse,qlength,budget);
  Vec<std::size_t> ends(nodes.size(),0,&budget);Vec<Position> starts(&budget);
  std::iota(ends.begin(),ends.end(),0);
  std::sort(ends.begin(),ends.end(),[&](auto a,auto b){return std::pair{nodes[a].anchor.reference_begin+nodes[a].anchor.length,a}<
    std::pair{nodes[b].anchor.reference_begin+nodes[b].anchor.length,b};});
  for(const auto& n:nodes)starts.push_back(n.anchor.reference_begin);
  std::sort(starts.begin(),starts.end());starts.erase(std::unique(starts.begin(),starts.end()),starts.end());
  for(std::size_t begin=0;begin<ends.size();) {
    const auto end=nodes[ends[begin]].anchor.reference_begin+nodes[ends[begin]].anchor.length;
    const auto bucket=static_cast<std::size_t>(std::lower_bound(starts.begin(),starts.end(),end)-starts.begin());
    forward.BeginEvent(bucket);reverse.BeginEvent(bucket);
    auto finish=begin+1;while(finish<ends.size() && nodes[ends[finish]].anchor.reference_begin+nodes[ends[finish]].anchor.length==end)++finish;
    for(auto e=begin;e<finish;++e) {
      Check(options);const auto endpoint=ends[e];
      for(std::uint32_t family=0;family<(regional ? 2U : 1U);++family) {
        TailsTo(endpoint,nodes,family ? *regional : graph,work,options,out);
        for(const auto first:work.visited) {
          // A front-clipped node still owes its preceding witness. Back-end
          // obligations and all actual witness connections are checked in the
          // reverse DAG and local graph, so unrelated nodes cannot bypass it.
          if(nodes[first].clipping.front!=none)continue;
          const auto& tail=work.tails[first];const auto item=LocalItem(first,endpoint,nodes,tail.support,tail.weight);
          forward.Add(item,first,endpoint,family,out);reverse.Add(item,first,endpoint,family,out);
        }
      }
      if(clusters)if(const auto* c=clusters->Candidate(endpoint)) {
        const auto first=static_cast<std::size_t>(c->first-1);
        const auto item=LocalItem(first,endpoint,nodes,c->support,c->weight);
        forward.Add(item,first,endpoint,2,out);reverse.Add(item,first,endpoint,2,out);
      }
    }
    forward.EndEvent(end);reverse.EndEvent(end);begin=finish;
    if(options.progress_callback)options.progress_callback("centromere-local-planning",begin,ends.size());
  }
  return {forward.Result(budget),reverse.Result(budget)};
}

void LocalPlan(CentromereBlockResult& out,const Vec<Node>& nodes,std::size_t reverse_base,Length rlength,Length qlength,
    Strand preferred_background,Budget& budget,const CentromereBlockOptions& options,
    std::string_view reference_bases,std::string_view query_bases) {
  LocalGraph graph(nodes.size(),budget);
  graph.query_length=qlength;
  const bool regions_only=options.planning_strategy==CentromereBlockStrategy::DirectionRegionsOnlyV5;
  const bool fast_direction=options.planning_strategy==CentromereBlockStrategy::FastLocalDirectionV6 ||
    options.planning_strategy==CentromereBlockStrategy::OrientationEvidenceV9;
  const bool count_direction=options.planning_strategy==CentromereBlockStrategy::CountClusterDirectionsV8;
  graph.hierarchical=options.planning_strategy==CentromereBlockStrategy::HierarchicalDirectionRegionsV4 || regions_only;
  LocalConnections(nodes,0,reverse_base,qlength,graph,budget,options);
  LocalConnections(nodes,reverse_base,nodes.size(),qlength,graph,budget,options);
  LocalGraph regional(graph.hierarchical ? nodes.size() : 0,budget);
  regional.query_length=qlength;regional.hierarchical=graph.hierarchical;
  regional.single_predecessor=regions_only;
  if(graph.hierarchical)for(std::size_t i=0;i<nodes.size();++i) {
    regional.incoming[i].begin=regional.edges.size();
    if(nodes[i].dp_predecessor!=none)regional.edges.push_back(nodes[i].dp_predecessor);
    regional.incoming[i].end=regional.edges.size();
  }
  if(graph.hierarchical) {
    for(auto* g:{&graph,&regional}) {
    g->increments.resize(g->edges.size());
    for(std::size_t i=0;i<nodes.size();++i) {
      Check(options);const auto range=g->incoming[i];
      for(auto e=range.begin;e<range.end;++e) {
        const auto cost=ConnectionCost(nodes[g->edges[e]],nodes[i],qlength);
        g->increments[e]=nodes[i].anchor.weight>cost ? nodes[i].anchor.weight-cost : 0.1;
      }
    }
    }
  }
  LocalWorkspace work(nodes.size(),budget);
  std::unique_ptr<CountClusters> clusters;
  if(count_direction)clusters=std::make_unique<CountClusters>(nodes,reverse_base,qlength,options,out,budget,reference_bases,query_bases);
  const auto selection_start=std::chrono::steady_clock::now();
  const auto* regional_graph=graph.hierarchical ? &regional : nullptr;
  auto solutions=(regions_only || fast_direction || count_direction) ? SolveLocalBoth(nodes,graph,regional_graph,qlength,work,budget,options,out,clusters.get()) :
    std::pair{SolveLocal(nodes,graph,Strand::Forward,qlength,work,budget,options,out,regional_graph),
              SolveLocal(nodes,graph,Strand::Reverse,qlength,work,budget,options,out,regional_graph)};
  const auto& forward=solutions.first;const auto& reverse=solutions.second;
  const auto comparison=Compare(forward.score,reverse.score);out.background_tied=comparison==0;
  out.background=comparison>0 ? Strand::Forward : comparison<0 ? Strand::Reverse : preferred_background;
  const auto& selected=out.background==Strand::Forward ? forward.path : reverse.path;
  Vec<std::size_t> by_end(nodes.size(),0,&budget);std::iota(by_end.begin(),by_end.end(),0);
  std::sort(by_end.begin(),by_end.end(),[&](auto a,auto b){return std::pair{nodes[a].anchor.reference_begin+nodes[a].anchor.length,a}<
    std::pair{nodes[b].anchor.reference_begin+nodes[b].anchor.length,b};});
  Reserve(out.blocks,2*selected.size()+1,budget);Position rb=0,ub=0;
  auto unknown=[&](Position re,Position ue){
    if(rb==re && ub==ue)return;
    CentromereBlock b;b.id=out.blocks.size()+1;b.reference_begin=rb;b.reference_end=re;
    b.query_begin=out.background==Strand::Forward ? ub : qlength-ue;b.query_end=out.background==Strand::Forward ? ue : qlength-ub;
    b.strand=out.background;b.direction_known=false;b.tentative=true;b.status="unresolved";b.reason="direction_unresolved";
    if(graph.hierarchical && !regions_only && rb<re && ub<ue) {
      if(out.refinement_regions.size()==out.refinement_regions.capacity())
        Reserve(out.refinement_regions,std::max<std::size_t>(1,out.refinement_regions.capacity()*2),budget);
      out.refinement_regions.push_back({rb,re,b.query_begin,b.query_end});
    }
    if((graph.hierarchical || count_direction) && out.blocks.size()==out.blocks.capacity())
      Reserve(out.blocks,std::max<std::size_t>(1,out.blocks.capacity()*2),budget);
    out.blocks.push_back(std::move(b));
  };
  for(std::size_t i=0;i<selected.size();++i) {
    const auto& s=selected[i];
    if(s.family==2) {
      const auto* c=clusters->Candidate(s.endpoint);
      if(!c)throw AlignmentError("missing selected count candidate");
      const auto item=LocalItem(s.first,s.endpoint,nodes,c->support,c->weight);
      const auto u0=Begin(item,out.background,qlength),u1=End(item,out.background,qlength);
      if(item.rb<rb || u0<ub)throw AlignmentError("count region overlaps on an input axis");
      unknown(item.rb,u0);
      const bool inversion=i>0 && i+1<selected.size() && nodes[s.endpoint].strand!=out.background &&
        nodes[selected[i-1].endpoint].strand==out.background && nodes[selected[i+1].endpoint].strand==out.background;
      clusters->Emit(s.endpoint,inversion);rb=item.re;ub=u1;continue;
    }
    const auto& selected_graph=s.family ? regional : graph;
    TailsTo(s.endpoint,nodes,selected_graph,work,options,out);
    const auto& tail=work.tails[s.first];const auto item=LocalItem(s.first,s.endpoint,nodes,tail.support,tail.weight);
    Vec<std::size_t> saved_path(&budget);
    if(regions_only)for(auto j=s.first;j!=none;j=work.tails[j].next)saved_path.push_back(j);
    const auto u0=Begin(item,out.background,qlength),u1=End(item,out.background,qlength);
    if(item.rb<rb || u0<ub)throw AlignmentError("local frontier plan overlaps on an input axis");
    unknown(item.rb,u0);
    CentromereBlock b;b.id=out.blocks.size()+1;b.reference_begin=item.rb;b.reference_end=item.re;
    b.query_begin=item.qb;b.query_end=item.qe;b.strand=nodes[s.endpoint].strand;b.skeleton=item.atom;
    b.forward_support=b.strand==Strand::Forward ? item.support : 0;b.reverse_support=b.strand==Strand::Reverse ? item.support : 0;
    b.inversion_candidate=i>0 && i+1<selected.size() && b.strand!=out.background &&
      nodes[selected[i-1].endpoint].strand==out.background && nodes[selected[i+1].endpoint].strand==out.background;
    b.reason=item.atom ? "local_exact_anchor" : "local_coherent_frontier_chain";
    std::size_t count=0;for(auto j=s.first;j!=none;j=work.tails[j].next)++count;
    if(!graph.hierarchical) {
      Reserve(b.planning_anchors,count,budget);
      for(auto j=s.first;j!=none;j=work.tails[j].next)b.planning_anchors.push_back(nodes[j].anchor);
    }
    bool ambiguous=false;
    auto alternative=std::lower_bound(by_end.begin(),by_end.end(),item.re,[&](auto j,auto end){return nodes[j].anchor.reference_begin+nodes[j].anchor.length<end;});
    for(;alternative!=by_end.end() && nodes[*alternative].anchor.reference_begin+nodes[*alternative].anchor.length==item.re;++alternative) {
      const auto endpoint=*alternative;const auto& other=nodes[endpoint];
      if(other.strand==b.strand || (other.strand==Strand::Forward ? other.anchor.query_begin+other.anchor.length!=item.qe : other.anchor.query_begin!=item.qb))continue;
      for(std::uint32_t family=0;family<(graph.hierarchical ? 2U : 1U);++family) {
      TailsTo(endpoint,nodes,family ? regional : graph,work,options,out);
      for(const auto first:work.visited)if(nodes[first].anchor.reference_begin==item.rb) {
        if(nodes[first].clipping.front!=none)continue;
        const auto& a=work.tails[first];const auto box=LocalItem(first,endpoint,nodes,a.support,a.weight);
        if(box.qb==item.qb && box.qe==item.qe && a.support==item.support && a.weight==item.weight){ambiguous=true;break;}
      }
      }
      if(ambiguous)break;
    }
    if(ambiguous) {
      b.direction_known=false;b.skeleton=false;b.tentative=true;b.inversion_candidate=false;
      b.status="unresolved";b.reason="direction_competition";budget.Release(b.planning_anchors.capacity()*sizeof(CentromereAnchor));
      std::vector<CentromereAnchor>().swap(b.planning_anchors);
    }
    if(graph.hierarchical && !ambiguous) {
      const auto region_id=out.regions.size()+1;
      CentromereDirectionRegion region;region.id=region_id;region.reference_begin=item.rb;region.reference_end=item.re;
      region.query_begin=item.qb;region.query_end=item.qe;region.strand=b.strand;region.support=item.support;
      region.secondary_score=item.weight;
      region.left_boundary_width=nodes[s.first].anchor.length;region.right_boundary_width=nodes[s.endpoint].anchor.length;
      for(const auto& n:nodes)if(n.strand!=b.strand && n.anchor.reference_begin<item.re &&
        n.anchor.reference_begin+n.anchor.length>item.rb && n.anchor.query_begin<item.qe && n.anchor.query_begin+n.anchor.length>item.qb)
        ++region.competing_mappings;
      // Retain the complete proposed region, while only coherent components
      // count as known-direction blocks. Unsupported internal gaps stay unknown.
      auto emit=[&](std::size_t first,std::size_t last,bool known) {
        const auto& a=nodes[first].anchor;const auto& z=nodes[last].anchor;
        CentromereBlock c;c.id=out.blocks.size()+1;c.region_id=region_id;c.strand=b.strand;
        c.reference_begin=a.reference_begin;c.reference_end=z.reference_begin+z.length;
        c.query_begin=std::min(a.query_begin,z.query_begin);c.query_end=std::max(a.query_begin+a.length,z.query_begin+z.length);
        c.direction_known=known;c.tentative=!known;c.skeleton=known && first==last;
        c.reason=known ? "hierarchical_supported_component" : "direction_unresolved";
        for(auto j=first;;j=work.tails[j].next) {
          if(c.planning_anchors.size()==c.planning_anchors.capacity())
            Reserve(c.planning_anchors,std::max<std::size_t>(1,c.planning_anchors.capacity()*2),budget);
          c.planning_anchors.push_back(nodes[j].anchor);
          (b.strand==Strand::Forward ? c.forward_support : c.reverse_support)+=nodes[j].anchor.length;
          if(j==last)break;
        }
        if(out.blocks.size()==out.blocks.capacity())Reserve(out.blocks,std::max<std::size_t>(1,out.blocks.capacity()*2),budget);
        ++region.block_count;out.blocks.push_back(std::move(c));
      };
      // Ambiguity checks above reuse the workspace. Reconstruct this region's
      // path only once, after those checks have finished.
      if(regions_only) {
        // Only this selected, nonoverlapping path is retained. Restore its
        // links after ambiguity checks rather than solving its entire DAG again.
        for(std::size_t k=0;k<saved_path.size();++k)
          work.tails[saved_path[k]].next=k+1<saved_path.size() ? saved_path[k+1] : none;
      } else TailsTo(s.endpoint,nodes,selected_graph,work,options,out);
      auto first=s.first;
      for(auto j=s.first;j!=none;) {
        const auto next=work.tails[j].next;
        if(next==none || NeedsRefinement(nodes[j],nodes[next],qlength)) {
          emit(first,j,true);
          if(next!=none) {
            CentromereBlock gap;gap.id=out.blocks.size()+1;gap.region_id=region_id;gap.strand=b.strand;
            gap.reference_begin=nodes[j].anchor.reference_begin+nodes[j].anchor.length;gap.reference_end=nodes[next].anchor.reference_begin;
            gap.query_begin=b.strand==Strand::Forward ? nodes[j].anchor.query_begin+nodes[j].anchor.length : nodes[next].anchor.query_begin+nodes[next].anchor.length;
            gap.query_end=b.strand==Strand::Forward ? nodes[next].anchor.query_begin : nodes[j].anchor.query_begin;
            gap.direction_known=false;gap.tentative=true;gap.status="unresolved";gap.reason="direction_unresolved";
            if(gap.reference_begin<gap.reference_end || gap.query_begin<gap.query_end) {
              if(out.blocks.size()==out.blocks.capacity())Reserve(out.blocks,std::max<std::size_t>(1,out.blocks.capacity()*2),budget);
              out.blocks.push_back(std::move(gap));
            }
          }
          first=next;
        }
        j=next;
      }
      if(out.regions.size()==out.regions.capacity())Reserve(out.regions,std::max<std::size_t>(1,out.regions.capacity()*2),budget);
      out.regions.push_back(region);
    } else {
      if(clusters)clusters->ExistingSupport(b);
      if(out.blocks.size()==out.blocks.capacity())Reserve(out.blocks,std::max<std::size_t>(1,out.blocks.capacity()*2),budget);
      out.blocks.push_back(std::move(b));
    }
    rb=item.re;ub=u1;
  }
  unknown(rlength,qlength);
  if(clusters)out.count_selection_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-selection_start).count();
  if(graph.hierarchical && !regions_only) {
    // Evidence enrichment is needed where the provisional solution is unknown.
    // The PROPOSAL's actual endpoints supply the new rectangle; the unknown
    // rectangle only identifies the problem, and never fixes its copy mapping.
    // Refining every losing repeated edge even inside a fully exact solution
    // would manufacture all tiny substrings without any unresolved question.
    Vec<std::pair<Position,Position>> rholes(&budget),qholes(&budget);
    for(const auto& b:out.blocks)if(!b.direction_known) {
      if(b.reference_begin<b.reference_end)rholes.emplace_back(b.reference_begin,b.reference_end);
      if(b.query_begin<b.query_end)qholes.emplace_back(b.query_begin,b.query_end);
      if(b.reference_begin<b.reference_end && b.query_begin<b.query_end) {
        if(out.refinement_regions.size()==out.refinement_regions.capacity())
          Reserve(out.refinement_regions,std::max<std::size_t>(1,out.refinement_regions.capacity()*2),budget);
        out.refinement_regions.push_back({b.reference_begin,b.reference_end,b.query_begin,b.query_end});
      }
    }
    std::sort(rholes.begin(),rholes.end());std::sort(qholes.begin(),qholes.end());
    auto touches=[](const auto& holes,Position b,Position e) {
      const auto p=std::lower_bound(holes.begin(),holes.end(),b,[](const auto& h,auto x){return h.second<=x;});
      return p!=holes.end() && p->first<e;
    };
    for(std::size_t i=0;i<nodes.size();++i) {
      const auto& a=nodes[i];const auto range=regional.incoming[i];
      for(auto e=range.begin;e<range.end;++e) {
        const auto& p=nodes[regional.edges[e]];
        if(!NeedsRefinement(p,a,qlength))continue;
        const auto rbegin=p.anchor.reference_begin+p.anchor.length,re=a.anchor.reference_begin;
        const auto qb=a.strand==Strand::Forward ? p.anchor.query_begin+p.anchor.length : a.anchor.query_begin+a.anchor.length;
        const auto qe=a.strand==Strand::Forward ? a.anchor.query_begin : p.anchor.query_begin;
        if(rbegin<re && qb<qe && touches(rholes,rbegin,re) && touches(qholes,qb,qe)) {
          if(out.refinement_regions.size()==out.refinement_regions.capacity())
            Reserve(out.refinement_regions,std::max<std::size_t>(1,out.refinement_regions.capacity()*2),budget);
          out.refinement_regions.push_back({rbegin,re,qb,qe});
        }
      }
    }
  }
  out.planner_reserved_bytes=budget.peak;
  out.planner_capacity_bytes=out.overlap_variants.capacity()*sizeof(CentromereOverlapVariant)+out.root_candidates.capacity()*sizeof(CentromerePlanningAnchor)+nodes.capacity()*sizeof(Node)+
    (graph.incoming.capacity()+regional.incoming.capacity())*sizeof(Range)+(graph.edges.capacity()+regional.edges.capacity())*sizeof(std::size_t)+
    (graph.increments.capacity()+regional.increments.capacity())*sizeof(double)+out.blocks.capacity()*sizeof(CentromereBlock)+
    work.tails.capacity()*sizeof(Tail)+(work.heap.capacity()+work.visited.capacity()+by_end.capacity())*sizeof(std::size_t);
  for(const auto& b:out.blocks)out.planner_capacity_bytes+=b.planning_anchors.capacity()*sizeof(CentromereAnchor);
}
// One bounded pass. The frozen blocks and roots are read until every review
// completes; only then is the replacement published. No accepted core is a root.
void ReviewLocalDirectionPlan(CentromereBlockResult& out,Length rlength,Length qlength,
    std::string_view reference,std::string_view query,Budget& budget,const CentromereBlockOptions& options) {
  if(!options.direction_plan_review)return;
  if(options.planning_strategy!=CentromereBlockStrategy::OrientationEvidenceV9)
    throw CentromereInputError("direction plan review requires orientation-evidence-v9");
  const auto started=std::chrono::steady_clock::now();
  using Evidence=detail::DirectionIslandEvidence;
  Vec<Evidence> roots(&budget);
  for(std::size_t i=0;i<out.root_candidates.size();++i) {
    const auto& r=out.root_candidates[i];
    roots.push_back({r.anchor,r.strand,i+1,true,r.reference_count>0 && r.query_count>0,r.reference_count,r.query_count});
  }
  auto origin=[&](const CentromereAnchor& a,Strand strand) {
    Evidence answer{a,strand,0,false,false,0,0};
    for(const auto& root:roots) {
      if(++out.direction_plan_review_checks>10000000)throw detail::CentromereBudgetExceeded("direction plan review checks exhausted");
      const auto& b=root.anchor;
      if(root.strand!=strand || !detail::direction_island_detail::SameMapping(answer,root,qlength) ||
         b.reference_begin>a.reference_begin || b.reference_begin+b.length<a.reference_begin+a.length ||
         b.query_begin>a.query_begin || b.query_begin+b.length<a.query_begin+a.length)continue;
      if(answer.source_proven && answer.source_id!=root.source_id) {
        const auto& old=roots[answer.source_id-1];
        if(old.anchor.reference_begin!=b.reference_begin || old.anchor.query_begin!=b.query_begin ||
           old.anchor.length!=b.length || old.reference_count!=root.reference_count || old.query_count!=root.query_count) {
          answer.source_proven=false;answer.source_id=0;return answer;
        }
      } else {
        answer.source_id=root.source_id;answer.source_proven=true;answer.occurrence_proven=root.occurrence_proven;
        answer.reference_count=root.reference_count;answer.query_count=root.query_count;
      }
    }
    return answer;
  };
  Vec<CentromereDirectionPlanReview> reviews(&budget);
  bool exhausted=false;
  try {
    for(const auto& b:out.blocks)if(b.direction_known)for(std::size_t i=1;i<b.planning_anchors.size();++i) {
      Check(options);
      const auto& l=b.planning_anchors[i-1];const auto& r=b.planning_anchors[i];
      const auto u=[&](const auto& a){return b.strand==Strand::Forward?a.query_begin:qlength-a.query_begin-a.length;};
      if(l.reference_begin+l.length>=r.reference_begin || u(l)+l.length>=u(r) ||
         r.reference_begin-l.reference_begin-l.length!=u(r)-u(l)-l.length)continue;
      const auto left=origin(l,b.strand),right=origin(r,b.strand);
      auto d=detail::ReviewExactDirectionGap(left,right,roots,rlength,qlength,reference,query,&budget,[&] {
        if(options.interruption_callback)options.interruption_callback("centromere-direction-plan-review");
        if(++out.direction_plan_review_checks>10000000)throw detail::CentromereBudgetExceeded("direction plan review checks exhausted");
      });
      reviews.push_back({reviews.size()+1,b.id,i-1,d.source_id,d.checks,
        d.core.reference_begin,d.core.reference_begin+d.core.length,d.core.query_begin,d.core.query_begin+d.core.length,
        d.strand,d.core.weight,d.accepted,d.reason});
      if(std::string_view(d.reason)=="budget_exhausted"){exhausted=true;break;}
    }
  } catch(const detail::CentromereBudgetExceeded&) {exhausted=true;}
  // A partial pass is diagnostic only. It cannot silently select a prefix.
  if(exhausted)for(auto& x:reviews){x.accepted=false;x.reason="batch_budget_exhausted";}
  const auto accepted=std::count_if(reviews.begin(),reviews.end(),[](const auto& x){return x.accepted;});
  Reserve(out.direction_plan_reviews,reviews.size(),budget);
  for(const auto& x:reviews)budget.Account(x.reason.capacity()+1);
  std::vector<CentromereBlock> replacement;
  if(accepted) {
    Reserve(replacement,out.blocks.size()+2*static_cast<std::size_t>(accepted),budget);
    auto retain=[&](const CentromereBlock& b,std::size_t begin,std::size_t end) {
      CentromereBlock z;
      z.strand=b.strand;z.direction_known=b.direction_known;z.tentative=b.tentative;z.skeleton=end-begin==1;
      z.inversion_candidate=b.inversion_candidate;z.status=b.status;z.reason=b.reason;z.region_id=b.region_id;
      Reserve(z.planning_anchors,end-begin,budget);
      for(auto i=begin;i<end;++i)z.planning_anchors.push_back(b.planning_anchors[i]);
      const auto& a=z.planning_anchors.front();const auto& t=z.planning_anchors.back();
      z.reference_begin=a.reference_begin;z.reference_end=t.reference_begin+t.length;
      z.query_begin=std::min(a.query_begin,t.query_begin);z.query_end=std::max(a.query_begin+a.length,t.query_begin+t.length);
      for(const auto& x:z.planning_anchors)(z.strand==Strand::Forward?z.forward_support:z.reverse_support)+=x.length;
      replacement.push_back(std::move(z));
    };
    for(const auto& b:out.blocks) {
      std::size_t begin=0;bool split=false;
      for(const auto& x:reviews)if(x.accepted && x.block_id==b.id) {
        retain(b,begin,x.left_index+1);begin=x.left_index+1;split=true;
        CentromereBlock core;core.strand=x.strand;core.reference_begin=x.reference_begin;core.reference_end=x.reference_end;
        core.query_begin=x.query_begin;core.query_end=x.query_end;core.inversion_candidate=true;core.skeleton=true;
        core.reason="exact_gap_opposite_core";
        Reserve(core.planning_anchors,1,budget);core.planning_anchors.push_back({x.reference_begin,x.query_begin,x.reference_end-x.reference_begin,x.weight});
        (core.strand==Strand::Forward?core.forward_support:core.reverse_support)=x.reference_end-x.reference_begin;
        replacement.push_back(std::move(core));
      }
      if(split)retain(b,begin,b.planning_anchors.size());
      else {
        // Account copied retained storage while the frozen source still lives.
        budget.Account(b.planning_anchors.size()*sizeof(CentromereAnchor));
        replacement.push_back(b);
      }
    }
    std::uint64_t oldbytes=out.blocks.capacity()*sizeof(CentromereBlock);
    for(const auto& b:out.blocks)oldbytes+=b.planning_anchors.capacity()*sizeof(CentromereAnchor);
    for(std::size_t i=0;i<replacement.size();++i)replacement[i].id=i+1;
    std::uint64_t newbytes=replacement.capacity()*sizeof(CentromereBlock);
    for(const auto& b:replacement)newbytes+=b.planning_anchors.capacity()*sizeof(CentromereAnchor);
    out.blocks.swap(replacement);
    std::vector<CentromereBlock>().swap(replacement);budget.Release(oldbytes);
    out.planner_capacity_bytes=out.planner_capacity_bytes>=oldbytes?out.planner_capacity_bytes-oldbytes+newbytes:newbytes;
  }
  for(auto& x:reviews)out.direction_plan_reviews.push_back(std::move(x));
  out.planner_capacity_bytes+=out.direction_plan_reviews.capacity()*sizeof(CentromereDirectionPlanReview);
  out.planner_reserved_bytes=std::max(out.planner_reserved_bytes,budget.peak);
  out.direction_plan_review_applied=true;
  out.direction_plan_review_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
}

} // namespace centromere_planning_detail
using namespace centromere_planning_detail;

CentromereBlockResult detail::ReplanDirectionRegions(Length rlength,Length qlength,
    std::vector<CentromerePlanningAnchor>& candidates,const CentromereBlockOptions& options,
    Strand preferred_background,std::string_view reference_bases,std::string_view query_bases) {
  if(options.overlap_compatibility &&
     options.planning_strategy!=CentromereBlockStrategy::LocalCoherentFrontierV3 &&
     options.planning_strategy!=CentromereBlockStrategy::FastLocalDirectionV6 &&
     options.planning_strategy!=CentromereBlockStrategy::OrientationEvidenceV9)
    throw CentromereInputError("overlap compatibility requires local-coherent-frontier-v3, fast-local-direction-v6 or orientation-evidence-v9");
  if(options.orientation_audit_known_blocks &&
     options.planning_strategy!=CentromereBlockStrategy::OrientationEvidenceV9)
    throw CentromereInputError("known-block orientation audit requires orientation-evidence-v9");
  if(options.recover_long_gap_ends &&
     (!options.orientation_connected_gaps || options.planning_strategy!=CentromereBlockStrategy::OrientationEvidenceV9))
    throw CentromereInputError("long-gap endpoint recovery requires connected orientation-evidence-v9");
  if(options.recover_short_gap_ends && !options.recover_sparse_pivots)
    throw CentromereInputError("short-gap endpoint recovery requires the sparse-pivot orientation stack");
  if(options.recover_mixed_gap_ends && !options.recover_short_gap_ends)
    throw CentromereInputError("mixed-gap endpoint recovery requires the short-gap orientation stack");
  if(options.recover_sparse_gap_ends && !options.recover_mixed_gap_ends)
    throw CentromereInputError("sparse-gap endpoint recovery requires the mixed-gap orientation stack");
  if(options.prefer_stronger_sparse_gap_end && !options.recover_sparse_gap_ends)
    throw CentromereInputError("stronger sparse-gap endpoint preference requires sparse-gap endpoint recovery");
  if(options.review_isolated_opposite_anchors && !options.review_local_competitive_gaps)
    throw CentromereInputError("isolated opposite review requires frozen local review layers");
  if(options.review_collinear_window_gaps && !options.review_collinear_fragment_connections)
    throw AlignmentError("collinear window-gap review requires collinear fragment review");
  if(options.review_alternate_original_endpoints && !options.review_collinear_window_gaps)
    throw CentromereInputError("alternate endpoint review requires frozen collinear window-gap review");
  if(options.review_neutral_window_gaps && !options.review_alternate_original_endpoints)
    throw CentromereInputError("neutral window review requires bounded alternate original endpoint review");
  if(options.review_collinear_fragment_connections && !options.review_isolated_opposite_anchors)
    throw CentromereInputError("collinear fragment review requires frozen isolated review layers");
  if(options.review_local_competitive_gaps && !options.unify_real_anchor_ends)
    throw CentromereInputError("local competitive review requires frozen unified endpoints");
  if(options.unify_real_anchor_ends && (!options.prefer_stronger_sparse_gap_end ||
     !options.orientation_bridge_unknown_gaps || !options.orientation_audit_known_blocks ||
     options.planning_strategy!=CentromereBlockStrategy::OrientationEvidenceV9))
    throw CentromereInputError("unified real endpoints require the frozen v9q layers");
  if(options.recover_sparse_pivots && !options.recover_long_gap_ends)
    throw CentromereInputError("sparse-pivot recovery requires long-gap endpoint recovery");
  if(!rlength || !qlength || !options.planner_memory_limit_bytes)
    throw CentromereInputError("joint planner requires nonempty inputs and a memory budget");
  if(preferred_background!=Strand::Forward && preferred_background!=Strand::Reverse)
    throw CentromereInputError("invalid preferred background");
  for(const auto& a:candidates) {
    const auto& x=a.anchor;
    if(!x.length || x.reference_begin>rlength || x.length>rlength-x.reference_begin ||
       x.query_begin>qlength || x.length>qlength-x.query_begin || !std::isfinite(x.weight) || x.weight<0 ||
       (a.strand!=Strand::Forward && a.strand!=Strand::Reverse))
      throw CentromereInputError("invalid joint planning anchor");
  }
  const auto raw_bytes=static_cast<std::uint64_t>(candidates.capacity())*sizeof(CentromerePlanningAnchor);
  Budget budget(options.planner_memory_limit_bytes,raw_bytes);
  const bool hierarchical=options.planning_strategy==CentromereBlockStrategy::HierarchicalDirectionRegionsV4 ||
    options.planning_strategy==CentromereBlockStrategy::DirectionRegionsOnlyV5;
  const bool local=options.planning_strategy==CentromereBlockStrategy::LocalCoherentFrontierV3 || hierarchical ||
    options.planning_strategy==CentromereBlockStrategy::FastLocalDirectionV6 ||
    options.planning_strategy==CentromereBlockStrategy::CountClusterDirectionsV8 ||
    options.planning_strategy==CentromereBlockStrategy::OrientationEvidenceV9;
  if((!reference_bases.empty() || !query_bases.empty()) &&
     (reference_bases.size()!=rlength || query_bases.size()!=qlength))
    throw CentromereInputError("planner sequence lengths do not match");
  CentromereBlockResult out;out.strategy=local ? options.planning_strategy : CentromereBlockStrategy::JointLocalBlocksV2;
  out.root_candidates=std::move(candidates);
  out.root_candidate_count=out.root_candidates.size();
  try {
  Vec<Node> nodes(&budget);
  Vec<Item> items(&budget);
  if(local) {
    CompactChains(out.root_candidates,Strand::Forward,qlength,nodes,budget,options,out,!hierarchical);
    const auto reverse_base=nodes.size();CompactChains(out.root_candidates,Strand::Reverse,qlength,nodes,budget,options,out,!hierarchical);
    if(options.overlap_compatibility)ValidateOverlapVariants(out,reference_bases,query_bases,options);
    LocalPlan(out,nodes,reverse_base,rlength,qlength,preferred_background,budget,options,reference_bases,query_bases);
    const auto review_checkpoint=budget.Live();
    try {
      ReviewLocalDirectionPlan(out,rlength,qlength,reference_bases,query_bases,budget,options);
    } catch(const detail::CentromereBudgetExceeded&) {
      // Every allocation is before the atomic block swap. Frozen science is
      // still intact; release diagnostic capacity and restore the live ledger.
      std::vector<CentromereDirectionPlanReview>().swap(out.direction_plan_reviews);
      budget.Restore(review_checkpoint);
      out.direction_plan_review_applied=true;out.direction_plan_review_budget_exhausted=true;
    } catch(const std::bad_alloc&) {
      std::vector<CentromereDirectionPlanReview>().swap(out.direction_plan_reviews);
      budget.Restore(review_checkpoint);
      out.direction_plan_review_applied=true;out.direction_plan_review_budget_exhausted=true;
    }
    if(options.planning_strategy==CentromereBlockStrategy::OrientationEvidenceV9)
      OrientationEvidence(nodes,qlength,query_bases,out,budget,options);
    return out;
  }
  for(const auto strand:{Strand::Forward,Strand::Reverse}) {
    Check(options);
    const auto base=nodes.size();
    CompactChains(out.root_candidates,strand,qlength,nodes,budget,options,out);
    for(auto i=base;i<nodes.size();++i) {
      const auto& n=nodes[i];
      if(!n.overlap_variant)items.push_back({n.anchor.reference_begin,n.anchor.reference_begin+n.anchor.length,
        n.anchor.query_begin,n.anchor.query_begin+n.anchor.length,n.anchor.length,n.anchor.weight,i,true});
      if(n.predecessor!=none) {
        const auto& first=nodes[n.first].anchor;
        items.push_back({first.reference_begin,n.anchor.reference_begin+n.anchor.length,
          std::min(first.query_begin,n.anchor.query_begin),
          std::max(first.query_begin+first.length,n.anchor.query_begin+n.anchor.length),
          n.support,n.weight,i,false});
      }
    }
  }
  if(options.overlap_compatibility)ValidateOverlapVariants(out,reference_bases,query_bases,options);
  auto forward=Solve(items,Strand::Forward,qlength,options,budget);
  auto reverse=Solve(items,Strand::Reverse,qlength,options,budget);
  // A palindrome can give exactly equal, opposite evidence for the same
  // rectangle. Stable ordering must not turn that tie into a known direction.
  Vec<std::size_t> by_box(items.size(),0,&budget);std::iota(by_box.begin(),by_box.end(),0);
  auto key=[&](std::size_t i){const auto& a=items[i];return std::tuple{a.rb,a.re,a.qb,a.qe,a.support,a.weight};};
  std::sort(by_box.begin(),by_box.end(),[&](auto a,auto b){return std::pair{key(a),a}<std::pair{key(b),b};});
  Vec<bool> ambiguous(items.size(),false,&budget);
  for(std::size_t i=0;i<by_box.size();) {
    auto j=i;unsigned directions=0;
    while(j<by_box.size() && key(by_box[j])==key(by_box[i])) {
      directions|=nodes[items[by_box[j]].endpoint].strand==Strand::Forward ? 1U : 2U;++j;
    }
    if(directions==3)for(auto k=i;k<j;++k)ambiguous[by_box[k]]=true;
    i=j;
  }
  const auto comparison=Compare(forward.score,reverse.score);
  out.background_tied=comparison==0;
  out.background=comparison>0 ? Strand::Forward : comparison<0 ? Strand::Reverse : preferred_background;
  const auto& selected=out.background==Strand::Forward ? forward.path : reverse.path;
  Reserve(out.blocks,2*selected.size()+1,budget);
  Position rb=0,ub=0;
  auto unknown=[&](Position re,Position ue) {
    if(rb==re && ub==ue)return;
    CentromereBlock b;b.id=out.blocks.size()+1;b.reference_begin=rb;b.reference_end=re;
    b.query_begin=out.background==Strand::Forward ? ub : qlength-ue;
    b.query_end=out.background==Strand::Forward ? ue : qlength-ub;
    b.strand=out.background;b.direction_known=false;b.tentative=true;
    b.status="unresolved";b.reason="direction_unresolved";out.blocks.push_back(std::move(b));
  };
  for(std::size_t i=0;i<selected.size();++i) {
    const auto& item=items[selected[i]];const auto& last=nodes[item.endpoint];
    const auto u0=Begin(item,out.background,qlength),u1=End(item,out.background,qlength);
    if(item.rb<rb || u0<ub)throw AlignmentError("joint plan overlaps on an input axis");
    unknown(item.rb,u0);
    CentromereBlock b;b.id=out.blocks.size()+1;b.reference_begin=item.rb;b.reference_end=item.re;
    b.query_begin=item.qb;b.query_end=item.qe;b.strand=last.strand;b.skeleton=item.atom;
    b.forward_support=b.strand==Strand::Forward ? item.support : 0;
    b.reverse_support=b.strand==Strand::Reverse ? item.support : 0;
    b.inversion_candidate=i>0 && i+1<selected.size() && b.strand!=out.background &&
      nodes[items[selected[i-1]].endpoint].strand==out.background &&
      nodes[items[selected[i+1]].endpoint].strand==out.background;
    b.reason=item.atom ? "joint_exact_anchor" : "joint_local_collinear_chain";
    if(ambiguous[selected[i]]) {
      b.direction_known=false;b.skeleton=false;b.tentative=true;b.inversion_candidate=false;
      b.status="unresolved";b.reason="direction_competition";
    }
    std::size_t count=0;
    for(auto n=b.direction_known ? item.endpoint : none;n!=none;n=nodes[n].predecessor){++count;if(item.atom)break;}
    if(count)Reserve(b.planning_anchors,count,budget);
    for(auto n=b.direction_known ? item.endpoint : none;n!=none;n=nodes[n].predecessor) {
      b.planning_anchors.push_back(nodes[n].anchor);
      if(item.atom)break;
    }
    std::reverse(b.planning_anchors.begin(),b.planning_anchors.end());
    out.blocks.push_back(std::move(b));rb=item.re;ub=u1;
  }
  unknown(rlength,qlength);
  out.planner_reserved_bytes=budget.peak;
  out.planner_capacity_bytes=out.overlap_variants.capacity()*sizeof(CentromereOverlapVariant)+raw_bytes+nodes.capacity()*sizeof(Node)+items.capacity()*sizeof(Item)+
    out.blocks.capacity()*sizeof(CentromereBlock)+
    (forward.path.capacity()+reverse.path.capacity()+by_box.capacity())*sizeof(std::size_t)+items.size();
  for(const auto& b:out.blocks)out.planner_capacity_bytes+=b.planning_anchors.capacity()*sizeof(CentromereAnchor);
  if(out.planner_capacity_bytes>options.planner_memory_limit_bytes)
    throw AlignmentError("joint retained planning capacity exceeded memory budget");
  return out;
  } catch(...) {candidates=std::move(out.root_candidates);throw;}
}
CentromereBlockResult PlanJointCentromereBlocks(Length rlength,Length qlength,
    std::vector<CentromerePlanningAnchor> candidates,const CentromereBlockOptions& options,Strand preferred_background) {
  return detail::ReplanDirectionRegions(rlength,qlength,candidates,options,preferred_background);
}
} // namespace ramag

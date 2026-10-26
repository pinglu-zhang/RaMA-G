#include "ramag/centromere_blocks.hpp"
#include "centromere_detail.hpp"
#include "centromere_orientation_connections.hpp"
#include "centromere_orientation_block_audit.hpp"
#include "centromere_orientation_long_gap.hpp"
#include "centromere_orientation_mixed_gap.hpp"
#include "centromere_orientation_sparse_gap.hpp"
#include "centromere_orientation_unified_ends.hpp"
#include "centromere_orientation_local_review.hpp"
#include "centromere_orientation_isolated_review.hpp"
#include "centromere_orientation_collinear_review.hpp"
#include "centromere_orientation_sparse_pivot.hpp"
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
namespace {
constexpr auto none=std::numeric_limits<std::size_t>::max();
struct Node {
  CentromereAnchor anchor;
  Strand strand;
  std::size_t predecessor{none},first{};
  Length support{};
  double weight{};
  Length dp_support{};
  double dp_score{};
  std::size_t dp_predecessor{none};
};
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
class Budget final : public std::pmr::memory_resource {
  std::uint64_t limit_,live_;
  void* do_allocate(std::size_t bytes,std::size_t alignment) override {
    Account(bytes);
    try{return std::pmr::new_delete_resource()->allocate(bytes,alignment);}
    catch(...){live_-=bytes;throw;}
  }
  void do_deallocate(void* p,std::size_t bytes,std::size_t alignment) override {
    std::pmr::new_delete_resource()->deallocate(p,bytes,alignment);live_-=bytes;
  }
  bool do_is_equal(const std::pmr::memory_resource& r) const noexcept override {return this==&r;}
 public:
  std::uint64_t peak{};
  Budget(std::uint64_t limit,std::uint64_t retained):limit_(limit),live_(retained),peak(retained) {
    if(retained>limit)throw detail::CentromereBudgetExceeded("joint root evidence exceeds planner memory budget");
  }
  void Account(std::uint64_t bytes) {
    if(bytes>limit_ || live_>limit_-bytes)
      throw detail::CentromereBudgetExceeded("joint planner live capacity budget exceeded: retained="+std::to_string(live_)+
        " requested="+std::to_string(bytes)+" limit="+std::to_string(limit_)+"; no candidates were truncated");
    live_+=bytes;peak=std::max(peak,live_);
  }
  void Release(std::uint64_t bytes){live_-=bytes;}
};
template<class T>using Vec=std::pmr::vector<T>;
template<class T>void Reserve(std::vector<T>& v,std::size_t size,Budget& budget) {
  const auto old=v.capacity()*sizeof(T);budget.Account(size*sizeof(T));
  v.reserve(size);budget.Release(old);
}
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

// The same length-first DP and RaMA secondary score as ScoreDirectionChain,
// with compact rows and a bounded allocator for every transient container.
// No source-ID/path copies or per-anchor provenance objects are needed here.
void CompactChains(const std::vector<CentromerePlanningAnchor>& input,Strand strand,Length qlength,
                   Vec<Node>& nodes,Budget& budget,const CentromereBlockOptions& options,bool normalize_only=false) {
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
    for(const auto& e:events) {
      if(!active.empty())weight+=static_cast<double>(e.position-prior)*active.rbegin()->first;
      if(e.add)++active[e.density];
      else {auto it=active.find(e.density);if(it==active.end())throw AlignmentError("invalid normalization events");if(!--it->second)active.erase(it);}
      prior=e.position;
    }
    Node n;n.anchor={pieces[i].r,strand==Strand::Forward ? pieces[i].u : qlength-pieces[i].u-(finish-pieces[i].r),finish-pieces[i].r,weight};
    n.strand=strand;nodes.push_back(n);i=end;
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
      if(u(j)+b.length>u(i))continue;
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

#include "centromere_count_clusters.inc"
#include "centromere_macro_sweep.inc"

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
  out.planner_capacity_bytes=out.root_candidates.capacity()*sizeof(CentromerePlanningAnchor)+nodes.capacity()*sizeof(Node)+
    (graph.incoming.capacity()+regional.incoming.capacity())*sizeof(Range)+(graph.edges.capacity()+regional.edges.capacity())*sizeof(std::size_t)+
    (graph.increments.capacity()+regional.increments.capacity())*sizeof(double)+out.blocks.capacity()*sizeof(CentromereBlock)+
    work.tails.capacity()*sizeof(Tail)+(work.heap.capacity()+work.visited.capacity()+by_end.capacity())*sizeof(std::size_t);
  for(const auto& b:out.blocks)out.planner_capacity_bytes+=b.planning_anchors.capacity()*sizeof(CentromereAnchor);
}
#include "centromere_orientation_evidence.inc"
} // namespace

CentromereBlockResult detail::ReplanDirectionRegions(Length rlength,Length qlength,
    std::vector<CentromerePlanningAnchor>& candidates,const CentromereBlockOptions& options,
    Strand preferred_background,std::string_view reference_bases,std::string_view query_bases) {
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
    CompactChains(out.root_candidates,Strand::Forward,qlength,nodes,budget,options,!hierarchical);
    const auto reverse_base=nodes.size();CompactChains(out.root_candidates,Strand::Reverse,qlength,nodes,budget,options,!hierarchical);
    LocalPlan(out,nodes,reverse_base,rlength,qlength,preferred_background,budget,options,reference_bases,query_bases);
    if(options.planning_strategy==CentromereBlockStrategy::OrientationEvidenceV9)
      OrientationEvidence(nodes,qlength,query_bases,out,budget,options);
    return out;
  }
  for(const auto strand:{Strand::Forward,Strand::Reverse}) {
    Check(options);
    const auto base=nodes.size();
    CompactChains(out.root_candidates,strand,qlength,nodes,budget,options);
    for(auto i=base;i<nodes.size();++i) {
      const auto& n=nodes[i];
      items.push_back({n.anchor.reference_begin,n.anchor.reference_begin+n.anchor.length,
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
  out.planner_capacity_bytes=raw_bytes+nodes.capacity()*sizeof(Node)+items.capacity()*sizeof(Item)+
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

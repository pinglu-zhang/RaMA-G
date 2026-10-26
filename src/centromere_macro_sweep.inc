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

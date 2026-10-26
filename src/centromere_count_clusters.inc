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

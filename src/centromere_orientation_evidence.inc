// v9 is deliberately downstream of the unmodified v6 plan. It assigns query
// orientation evidence, not a reference copy and not an alignment.
void OrientationEvidence(const Vec<Node>& nodes,Length qlength,std::string_view query,
    CentromereBlockResult& out,Budget& budget,const CentromereBlockOptions& options) {
  const auto start=std::chrono::steady_clock::now();
  out.orientation_connected_gaps=options.orientation_connected_gaps;
  out.orientation_bridge_unknown_gaps=options.orientation_bridge_unknown_gaps;
  out.orientation_audit_known_blocks=options.orientation_audit_known_blocks;
  out.recover_long_gap_ends=options.recover_long_gap_ends;
  out.recover_sparse_pivots=options.recover_sparse_pivots;
  out.recover_short_gap_ends=options.recover_short_gap_ends;
  out.recover_mixed_gap_ends=options.recover_mixed_gap_ends;
  out.recover_sparse_gap_ends=options.recover_sparse_gap_ends;
  out.prefer_stronger_sparse_gap_end=options.prefer_stronger_sparse_gap_end;
  out.unify_real_anchor_ends=options.unify_real_anchor_ends;
  out.review_local_competitive_gaps=options.review_local_competitive_gaps;
  out.review_isolated_opposite_anchors=options.review_isolated_opposite_anchors;
  out.review_collinear_fragment_connections=options.review_collinear_fragment_connections;
  out.review_collinear_window_gaps=options.review_collinear_window_gaps;
  out.review_alternate_original_endpoints=options.review_alternate_original_endpoints;
  out.review_neutral_window_gaps=options.review_neutral_window_gaps;
  if(options.orientation_bridge_unknown_gaps && !options.orientation_connected_gaps)
    throw CentromereInputError("unknown-gap bridging requires connected orientation gaps");
  if(query.size()!=qlength || !options.count_window || !options.count_step ||
     options.count_dominance_ratio<2 || !options.count_min_anchors || !options.count_inferred_gap)
    throw CentromereInputError("v9 requires query sequence and positive window/count/gap settings (ratio >= 2)");
  using Interval=CentromereQueryInterval;
  std::uint64_t sparse_gap_reason_bytes=0,unified_reason_bytes=0,local_review_reason_bytes=0,isolated_reason_bytes=0,collinear_reason_bytes=0;
  std::uint64_t alternate_reason_bytes=0;
  auto append=[&]<class T>(std::vector<T>& dst,T row) {
    if(dst.size()==dst.capacity())Reserve(dst,std::max<std::size_t>(1,dst.capacity()*2),budget);
    dst.push_back(std::move(row));
  };
  auto merge=[](Vec<Interval>& rows) {
    std::sort(rows.begin(),rows.end(),[](const auto& a,const auto& b){return std::pair{a.begin,a.end}<std::pair{b.begin,b.end};});
    std::size_t size=0;
    for(const auto& row:rows) {
      if(row.begin>=row.end)continue;
      if(size && row.begin<=rows[size-1].end)rows[size-1].end=std::max(rows[size-1].end,row.end);
      else rows[size++]=row;
    }
    rows.resize(size);
  };
  // Neutrality is based on the complete, actual two-axis mapping rectangle.
  Vec<std::size_t> order(nodes.size(),0,&budget);std::iota(order.begin(),order.end(),0);
  auto box=[&](auto i){const auto& a=nodes[i].anchor;return std::tuple{a.reference_begin,a.query_begin,a.length};};
  std::sort(order.begin(),order.end(),[&](auto a,auto b){return std::pair{box(a),a}<std::pair{box(b),b};});
  Vec<bool> neutral(nodes.size(),false,&budget);
  for(std::size_t i=0;i<order.size();) {
    auto j=i;unsigned signs=0;
    while(j<order.size() && box(order[j])==box(order[i])){signs|=nodes[order[j]].strand==Strand::Forward?1U:2U;++j;}
    if(signs==3)for(auto k=i;k<j;++k)neutral[order[k]]=true;
    i=j;
  }
  Vec<Interval> forward(&budget),reverse(&budget),known(&budget),unknown(&budget);
  for(std::size_t i=0;i<nodes.size();++i)if(!neutral[i]) {
    const auto& a=nodes[i].anchor;
    (nodes[i].strand==Strand::Forward?forward:reverse).push_back({a.query_begin,a.query_begin+a.length});
  }
  merge(forward);merge(reverse);
  for(const auto& b:out.blocks)if(b.direction_known && b.query_begin<b.query_end)known.push_back({b.query_begin,b.query_end});
  merge(known);Position cursor=0;
  for(const auto& b:known){if(cursor<b.begin)unknown.push_back({cursor,b.begin});cursor=b.end;}
  if(cursor<qlength)unknown.push_back({cursor,qlength});
  Vec<Length> fp(forward.size()+1,0,&budget),rp(reverse.size()+1,0,&budget);
  for(std::size_t i=0;i<forward.size();++i)fp[i+1]=fp[i]+forward[i].end-forward[i].begin;
  for(std::size_t i=0;i<reverse.size();++i)rp[i+1]=rp[i]+reverse[i].end-reverse[i].begin;
  auto cumulative=[](const auto& intervals,const auto& sums,Position p) {
    const auto it=std::lower_bound(intervals.begin(),intervals.end(),p,[](const auto& a,auto x){return a.end<=x;});
    const auto id=static_cast<std::size_t>(it-intervals.begin());
    return sums[id]+(it!=intervals.end() && it->begin<p ? p-it->begin : 0);
  };
  order.clear();
  for(std::size_t i=0;i<nodes.size();++i)if(!neutral[i])order.push_back(i);
  auto midpoint=[&](auto i){const auto& a=nodes[i].anchor;return a.query_begin+a.length/2;};
  std::sort(order.begin(),order.end(),[&](auto a,auto b){return std::pair{midpoint(a),a}<std::pair{midpoint(b),b};});
  Vec<Position> starts(&budget);
  const auto width=std::min(qlength,options.count_window),last=qlength-width;
  for(Position p=0;;) {
    starts.push_back(p);starts.push_back(last-p);
    if(last-p<options.count_step)break;
    p+=options.count_step;
  }
  std::sort(starts.begin(),starts.end());starts.erase(std::unique(starts.begin(),starts.end()),starts.end());
  Vec<std::size_t> window_ids(&budget);
  for(auto begin:starts) {
    Check(options);CentromereOrientationWindow w;w.id=out.orientation_windows.size()+1;w.begin=begin;w.end=begin+width;
    const auto lo=std::lower_bound(order.begin(),order.end(),begin,[&](auto i,auto p){return midpoint(i)<p;});
    const auto hi=std::upper_bound(order.begin(),order.end(),w.end,[&](auto p,auto i){return p<midpoint(i);});
    for(auto sign:{Strand::Forward,Strand::Reverse}) {
      window_ids.clear();
      for(auto it=lo;it!=hi;++it) {
        const auto& a=nodes[*it].anchor;
        // Preserve v8's two-central-base convention for exact whole-query RC symmetry.
        if(nodes[*it].strand==sign && a.query_begin+(a.length-1)/2<w.end)window_ids.push_back(*it);
      }
      std::sort(window_ids.begin(),window_ids.end(),[&](auto a,auto b) {
        const auto& x=nodes[a].anchor;const auto& y=nodes[b].anchor;
        return std::tuple{x.query_begin+x.length,x.query_begin,a}<std::tuple{y.query_begin+y.length,y.query_begin,b};
      });
      Position end=0;std::uint64_t count=0;
      for(auto id:window_ids){const auto& a=nodes[id].anchor;if(a.query_begin>=end){++count;end=a.query_begin+a.length;}}
      (sign==Strand::Forward?w.forward:w.reverse)=count;
    }
    w.forward_length=cumulative(forward,fp,w.end)-cumulative(forward,fp,w.begin);
    w.reverse_length=cumulative(reverse,rp,w.end)-cumulative(reverse,rp,w.begin);
    auto dominates=[&](auto own,auto other,auto own_length,auto other_length) {
      return own>=options.count_min_anchors && other<=own/options.count_dominance_ratio &&
        other_length<=own_length/options.count_dominance_ratio;
    };
    w.label=dominates(w.forward,w.reverse,w.forward_length,w.reverse_length)?'+':
      dominates(w.reverse,w.forward,w.reverse_length,w.forward_length)?'-':
      std::max(w.forward,w.reverse)>=options.count_min_anchors?'?':'.';
    append(out.orientation_windows,w);
  }
  // A sweep labels disjoint positions; opposing/competitive windows are vetoes.
  // Insufficient windows do not themselves claim either direction.
  struct Event {Position p;unsigned sign;int delta;};
  Vec<Event> events(&budget);
  for(const auto& w:out.orientation_windows)if(w.label!='.') {
    const unsigned sign=w.label=='+'?0U:w.label=='-'?1U:2U;
    events.push_back({w.begin,sign,1});events.push_back({w.end,sign,-1});
  }
  std::sort(events.begin(),events.end(),[](const auto& a,const auto& b){return a.p<b.p;});
  using Label=detail::OrientationCoverage;
  Vec<Label> labels(&budget);std::array<std::int64_t,3> active{};Position previous=0;
  for(std::size_t i=0;i<events.size();) {
    const auto p=events[i].p;const unsigned mask=(active[0]>0?1U:0U)|(active[1]>0?2U:0U)|(active[2]>0?4U:0U);
    if(previous<p)labels.push_back({previous,p,mask});
    do{active[events[i].sign]+=events[i].delta;++i;}while(i<events.size() && events[i].p==p);
    previous=p;
  }
  Vec<CentromereOrientationInterval> support(&budget);
  for(auto sign:{Strand::Forward,Strand::Reverse}) {
    const auto& intervals=sign==Strand::Forward?forward:reverse;const unsigned bit=sign==Strand::Forward?1U:2U;
    std::size_t a=0,b=0,c=0;
    while(a<intervals.size() && b<labels.size() && c<unknown.size()) {
      const auto begin=std::max({intervals[a].begin,labels[b].begin,unknown[c].begin});
      const auto end=std::min({intervals[a].end,labels[b].end,unknown[c].end});
      if(begin<end && labels[b].mask==bit) {
        if(!support.empty() && support.back().strand==sign && support.back().end==begin)support.back().end=end;
        else support.push_back({0,begin,end,sign,false,0});
      }
      if(intervals[a].end==end)++a;
      if(labels[b].end==end)++b;
      if(unknown[c].end==end)++c;
    }
  }
  std::sort(support.begin(),support.end(),[](const auto& a,const auto& b){return std::pair{a.begin,a.end}<std::pair{b.begin,b.end};});
  Vec<detail::OrientationEndpoint> frozen_endpoints(&budget);
  if(options.orientation_connected_gaps) {
    Vec<detail::OrientationEndpoint> endpoints(&budget);
    endpoints.reserve(support.size()+out.blocks.size());
    for(const auto& row:support)endpoints.push_back({row.begin,row.end,row.strand,0});
    for(const auto& block:out.blocks)if(block.direction_known && block.query_begin<block.query_end)
      endpoints.push_back({block.query_begin,block.query_end,block.strand,block.id});
    std::sort(endpoints.begin(),endpoints.end(),[](const auto& a,const auto& b){return a.begin<b.begin;});
    detail::VisitConnectedOrientationGaps(endpoints,labels,query,options.count_inferred_gap,
      [&](const auto& gap) {
        append(out.orientation_gaps,gap);
        if(gap.accepted)append(out.orientation_intervals,CentromereOrientationInterval{
          0,gap.begin,gap.end,gap.strand,true,0});
      },[&]{Check(options);},options.orientation_bridge_unknown_gaps,forward,reverse);
    for(const auto& row:support)append(out.orientation_intervals,row);
    std::sort(out.orientation_intervals.begin(),out.orientation_intervals.end(),
      [](const auto& a,const auto& b){return a.begin<b.begin;});
    for(std::size_t i=0;i<out.orientation_intervals.size();++i)out.orientation_intervals[i].id=i+1;
    if(options.recover_long_gap_ends)frozen_endpoints=std::move(endpoints);
  } else {
  for(std::size_t i=0;i<support.size();++i) {
    Check(options);
    if(i && support[i-1].strand==support[i].strand && support[i-1].end<support[i].begin) {
      const auto begin=support[i-1].end,end=support[i].begin;const auto sign=support[i].strand;
      CentromereOrientationGap gap;gap.begin=begin;gap.end=end;gap.strand=sign;gap.reason="gap_too_long";
      if(end-begin<=options.count_inferred_gap) {
        auto hole=std::upper_bound(unknown.begin(),unknown.end(),begin,[](auto p,const auto& x){return p<x.begin;});
        const bool within=hole!=unknown.begin() && (--hole)->end>=end;
        gap.reason="crosses_known_direction";
        if(within) {
          gap.reason="contains_N";
          if(query.substr(static_cast<std::size_t>(begin),static_cast<std::size_t>(end-begin)).find('N')==std::string_view::npos) {
            gap.reason="conflicting_window";
            auto label=std::lower_bound(labels.begin(),labels.end(),begin,[](const auto& x,auto p){return x.end<=p;});
            bool blocked=false;const unsigned bit=sign==Strand::Forward?1U:2U;
            for(;label!=labels.end() && label->begin<end;++label)if(label->mask!=bit){blocked=true;break;}
            if(!blocked) {
              gap.reason="no_covering_window";
              // Fixed window size => end positions have the same sorted order.
              auto w=std::lower_bound(out.orientation_windows.begin(),out.orientation_windows.end(),end,
                [](const auto& x,auto p){return x.end<p;});
              for(;w!=out.orientation_windows.end() && w->begin<=begin;++w)
                if(w->label==(sign==Strand::Forward?'+':'-')){gap.window_id=w->id;gap.accepted=true;gap.reason="count_inferred_gap";break;}
            }
          }
        }
      }
      if(gap.accepted)append(out.orientation_intervals,CentromereOrientationInterval{
        out.orientation_intervals.size()+1,begin,end,sign,true,gap.window_id});
      append(out.orientation_gaps,std::move(gap));
    }
    auto row=support[i];row.id=out.orientation_intervals.size()+1;append(out.orientation_intervals,row);
  }
  }
  // Map each raw candidate back to its normalized diagonal segment, then retain
  // provenance as IDs only. A repeated copy may be evidence, but never extra votes.
  auto diagonal=[&](const auto& a,Strand sign) {
    const auto u=sign==Strand::Forward?a.query_begin:qlength-a.query_begin-a.length;
    return std::tuple{sign==Strand::Forward?0U:1U,a.reference_begin>=u,
      a.reference_begin>=u?a.reference_begin-u:u-a.reference_begin,a.reference_begin};
  };
  order.resize(nodes.size());std::iota(order.begin(),order.end(),0);
  std::sort(order.begin(),order.end(),[&](auto a,auto b){return diagonal(nodes[a].anchor,nodes[a].strand)<diagonal(nodes[b].anchor,nodes[b].strand);});
  for(std::size_t i=0;i<out.root_candidates.size();++i) {
    if(i%8192==0)Check(options);
    const auto& source=out.root_candidates[i];const auto& a=source.anchor;const auto key=diagonal(a,source.strand);
    auto n=std::upper_bound(order.begin(),order.end(),key,[&](const auto& k,auto j){return k<diagonal(nodes[j].anchor,nodes[j].strand);});
    if(n==order.begin())throw AlignmentError("v9 source is missing its normalized match");
    --n;const auto& normalized=nodes[*n];
    const auto nk=diagonal(normalized.anchor,normalized.strand);
    if(std::get<0>(key)!=std::get<0>(nk) || std::get<1>(key)!=std::get<1>(nk) ||
       std::get<2>(key)!=std::get<2>(nk) || a.reference_begin+a.length>normalized.anchor.reference_begin+normalized.anchor.length)
      throw AlignmentError("v9 normalization provenance mismatch");
    if(neutral[*n])continue;
    auto row=std::lower_bound(out.orientation_intervals.begin(),out.orientation_intervals.end(),a.query_begin,
      [](const auto& x,auto p){return x.end<=p;});
    for(;row!=out.orientation_intervals.end() && row->begin<a.query_begin+a.length;++row)
      if(!row->inferred && row->strand==source.strand)
        append(out.orientation_sources,CentromereOrientationSource{row->id,i+1});
  }
  // A final direction-only audit cannot alter discovery, selected correspondence,
  // old overlays or their endpoints. Retracted old spans stay unknown: do not run
  // support enrichment or gap bridging again on this newly exposed space.
  if(options.orientation_audit_known_blocks) {
    const auto audit_start=std::chrono::steady_clock::now();
    known.clear();
    detail::VisitKnownOrientationBlockAudit(out.blocks,labels,qlength,options.count_inferred_gap,
      [&](const auto& row){append(out.orientation_block_audits,row);},
      [&](const auto& row) {
        append(out.orientation_known_intervals,row);
        known.push_back({row.begin,row.end});
      },[&]{Check(options);});
    merge(known);
    out.orientation_audit_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-audit_start).count();
  }
  // Residual query ledger uses audited known pieces when enabled. Raw blocks
  // remain unchanged for correspondence/provenance and are not the final view.
  Vec<Interval> occupied(known,&budget);
  for(const auto& row:out.orientation_intervals)occupied.push_back({row.begin,row.end});
  merge(occupied);
  if(options.recover_long_gap_ends) {
    const auto recovery_start=std::chrono::steady_clock::now();
    Vec<Interval> known_forward(&budget),known_reverse(&budget);
    auto add_known=[&](const auto& row) {
      (row.strand==Strand::Forward?known_forward:known_reverse).push_back({row.begin,row.end});
    };
    if(options.orientation_audit_known_blocks) {
      for(const auto& row:out.orientation_known_intervals)add_known(row);
    } else {
      for(const auto& b:out.blocks)if(b.direction_known && b.query_begin<b.query_end)
        (b.strand==Strand::Forward?known_forward:known_reverse).push_back({b.query_begin,b.query_end});
    }
    for(const auto& row:out.orientation_intervals)add_known(row);
    merge(known_forward);merge(known_reverse);
    // Only the original long gaps are visited. Auditing may remove a pivot, but
    // cannot create one; current overlays and new recoveries are never endpoints.
    detail::VisitLongOrientationGapEnds(frozen_endpoints,labels,query,options.count_inferred_gap,
      known_forward,known_reverse,forward,reverse,occupied,[&](const auto& gap) {
        append(out.orientation_long_gap_gaps,gap);
        if(gap.accepted)append(out.orientation_long_gap_intervals,CentromereOrientationInterval{
          out.orientation_long_gap_intervals.size()+1,gap.begin,gap.end,gap.strand,true,0});
      },[&]{Check(options);});
    for(const auto& row:out.orientation_long_gap_intervals)occupied.push_back({row.begin,row.end});
    merge(occupied);
    out.orientation_long_gap_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-recovery_start).count();
    if(options.recover_sparse_pivots) {
      const auto sparse_start=std::chrono::steady_clock::now();
      detail::VisitSparseOrientationPivots(frozen_endpoints,labels,query,options.count_inferred_gap,
        known_forward,known_reverse,forward,reverse,occupied,&budget,
        [&](auto row) {
          row.id=out.orientation_sparse_intervals.size()+1;
          append(out.orientation_sparse_intervals,row);
        },[&](const auto& row){append(out.orientation_sparse_gaps,row);},[&]{Check(options);});
      // Provenance uses original candidate IDs. It never constrains a query-only
      // direction to one reference copy, nor adds a second candidate store.
      for(std::size_t i=0;i<out.root_candidates.size();++i) {
        if(i%8192==0)Check(options);
        const auto& source=out.root_candidates[i];const auto& a=source.anchor;
        const auto key=diagonal(a,source.strand);
        auto n=std::upper_bound(order.begin(),order.end(),key,
          [&](const auto& k,auto j){return k<diagonal(nodes[j].anchor,nodes[j].strand);});
        if(n==order.begin())throw AlignmentError("sparse source missing normalized match");
        --n;if(neutral[*n])continue;
        auto row=std::lower_bound(out.orientation_sparse_intervals.begin(),out.orientation_sparse_intervals.end(),a.query_begin,
          [](const auto& x,auto p){return x.end<=p;});
        for(;row!=out.orientation_sparse_intervals.end() && row->begin<a.query_begin+a.length;++row)
          if(!row->inferred && row->strand==source.strand)
            append(out.orientation_sparse_sources,CentromereOrientationSource{row->id,i+1});
      }
      for(const auto& row:out.orientation_sparse_intervals)occupied.push_back({row.begin,row.end});
      merge(occupied);
      out.orientation_sparse_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-sparse_start).count();
    }
    if(options.recover_short_gap_ends) {
      const auto short_start=std::chrono::steady_clock::now();
      // Original endpoints and post-audit pivot support are frozen. Occupancy
      // includes every previous overlay, but those overlays cannot become pivots.
      detail::VisitLongOrientationGapEnds<true>(frozen_endpoints,labels,query,options.count_inferred_gap,
        known_forward,known_reverse,forward,reverse,occupied,[&](const auto& gap) {
          append(out.orientation_short_gap_gaps,gap);
          if(gap.accepted)append(out.orientation_short_gap_intervals,CentromereOrientationInterval{
            out.orientation_short_gap_intervals.size()+1,gap.begin,gap.end,gap.strand,true,0});
        },[&]{Check(options);});
      for(const auto& row:out.orientation_short_gap_intervals)occupied.push_back({row.begin,row.end});
      merge(occupied);
      out.orientation_short_gap_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-short_start).count();
    }
    if(options.recover_mixed_gap_ends) {
      const auto mixed_start=std::chrono::steady_clock::now();
      // Original mixed endpoints retain their own qualification and direction.
      // Occupancy is frozen after all old layers; no new tip becomes an endpoint.
      detail::VisitMixedOrientationGapEnds(frozen_endpoints,labels,query,options.count_inferred_gap,
        known_forward,known_reverse,forward,reverse,occupied,[&](const auto& gap) {
          append(out.orientation_mixed_gap_gaps,gap);
          if(gap.accepted)append(out.orientation_mixed_gap_intervals,CentromereOrientationInterval{
            out.orientation_mixed_gap_intervals.size()+1,gap.begin,gap.end,gap.strand,true,0});
        },[&]{Check(options);});
      for(const auto& row:out.orientation_mixed_gap_intervals)occupied.push_back({row.begin,row.end});
      merge(occupied);
      out.orientation_mixed_gap_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-mixed_start).count();
    }
    if(options.recover_sparse_gap_ends) {
      const auto sparse_gap_start=std::chrono::steady_clock::now();
      auto retain=[&](std::vector<CentromereSparseGapEnd>& rows,CentromereSparseGapEnd row) {
        // Account the retained reason buffers as well as the vector slots. This
        // deliberately includes any inline string storage a library may use;
        // capacity is a conservative estimate, not the measured process RSS.
        const auto bytes=static_cast<std::uint64_t>(row.reason.capacity())+1+
          (row.scan_reason.empty()?0:static_cast<std::uint64_t>(row.scan_reason.capacity())+1);
        budget.Account(bytes);sparse_gap_reason_bytes+=bytes;
        append(rows,std::move(row));
      };
      // The old sparse graph and its exact reachable roots are immutable. No
      // inferred tip is inserted back into the graph or used by a second pass.
      detail::VisitSparseOrientationGapEnds(out.orientation_sparse_gaps,out.orientation_sparse_intervals,
        labels,query,options.count_inferred_gap,forward,reverse,occupied,[&](const auto& gap) {
          retain(out.orientation_sparse_gap_gaps,gap);
          if(gap.accepted) {
            auto row=gap;row.id=out.orientation_sparse_gap_intervals.size()+1;
            retain(out.orientation_sparse_gap_intervals,std::move(row));
          }
        },[&]{Check(options);},options.prefer_stronger_sparse_gap_end);
      for(const auto& row:out.orientation_sparse_gap_intervals)occupied.push_back({row.begin,row.end});
      merge(occupied);
      out.orientation_sparse_gap_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-sparse_gap_start).count();
    }
  }
  if(options.unify_real_anchor_ends) {
    const auto unified_start=std::chrono::steady_clock::now();
    // Qualification is frozen before new evidence is emitted. Every source is
    // an ORIGINAL exact mapping; acceptance clipping never supplies a boundary.
    Vec<CentromereOrientationKnownInterval> audited(out.orientation_known_intervals.begin(),
      out.orientation_known_intervals.end(),&budget);
    std::sort(audited.begin(),audited.end(),[](const auto& a,const auto& b){return a.begin<b.begin;});
    auto nonneutral_mapping=[&](const CentromereAnchor& a,Strand sign) {
      if(!a.length || a.query_begin>qlength || a.length>qlength-a.query_begin)return false;
      const auto key=diagonal(a,sign);
      auto n=std::upper_bound(order.begin(),order.end(),key,
        [&](const auto& k,auto j){return k<diagonal(nodes[j].anchor,nodes[j].strand);});
      if(n==order.begin())return false;
      --n;const auto& normalized=nodes[*n];const auto nk=diagonal(normalized.anchor,normalized.strand);
      return !neutral[*n] && std::get<0>(key)==std::get<0>(nk) && std::get<1>(key)==std::get<1>(nk) &&
        std::get<2>(key)==std::get<2>(nk) && a.reference_begin>=normalized.anchor.reference_begin &&
        a.reference_begin-normalized.anchor.reference_begin<=normalized.anchor.length &&
        a.length<=normalized.anchor.length-(a.reference_begin-normalized.anchor.reference_begin);
    };
    auto pivot=[&](const CentromereAnchor& a,Strand sign,bool from_left,
        std::uint8_t kind,std::uint64_t source_id,std::uint64_t block_id) {
      append(out.orientation_unified_pivots,CentromereUnifiedPivot{0,source_id,block_id,
        a.reference_begin,a.reference_begin+a.length,a.query_begin,a.query_begin+a.length,sign,from_left,kind});
    };
    std::uint64_t planning_id=0;
    for(const auto& block:out.blocks)for(const auto& a:block.planning_anchors) {
      ++planning_id;Check(options);
      if(!block.direction_known || !nonneutral_mapping(a,block.strand))continue;
      for(bool from_left:{false,true})if(detail::QualifyUnifiedPlanningEnd(a,block.id,block.strand,from_left,audited))
        pivot(a,block.strand,from_left,0,planning_id,block.id);
    }
    for(std::size_t i=0;i<out.root_candidates.size();++i) {
      if(i%8192==0)Check(options);
      const auto& root=out.root_candidates[i];const auto& a=root.anchor;
      if(!nonneutral_mapping(a,root.strand))continue;
      const std::span<const Interval> own=root.strand==Strand::Forward?forward:reverse;
      for(bool from_left:{false,true}) {
        // Existing provenance lists all original mappings intersecting accepted
        // query-only exact evidence. No unique reference copy is implied here.
        if(detail::QualifyUnifiedRootEnd(a,root.strand,from_left,false,own,
            std::span<const CentromereOrientationInterval>(out.orientation_intervals)))
          pivot(a,root.strand,from_left,1,i+1,0);
        else if(detail::QualifyUnifiedRootEnd(a,root.strand,from_left,false,own,
            std::span<const CentromereSparseOrientationInterval>(out.orientation_sparse_intervals)))
          pivot(a,root.strand,from_left,2,i+1,0);
      }
    }
    auto endpoint_key=[](const auto& p) {
      return std::tuple{p.strand==Strand::Forward?0U:1U,p.from_left,p.from_left?p.query_end:p.query_begin};
    };
    auto& pivots=out.orientation_unified_pivots;
    std::sort(pivots.begin(),pivots.end(),[&](const auto& a,const auto& b) {
      return std::tuple{endpoint_key(a),a.source_kind,a.source_id,a.block_id}<
        std::tuple{endpoint_key(b),b.source_kind,b.source_id,b.block_id};
    });
    pivots.erase(std::unique(pivots.begin(),pivots.end(),[&](const auto& a,const auto& b) {
      return endpoint_key(a)==endpoint_key(b);
    }),pivots.end());
    for(std::size_t i=0;i<pivots.size();++i)pivots[i].id=i+1;
    Vec<Interval> known_forward(&budget),known_reverse(&budget);
    auto add_known=[&](const auto& row) {
      (row.strand==Strand::Forward?known_forward:known_reverse).push_back({row.begin,row.end});
    };
    for(const auto& row:audited)add_known(row);
    for(const auto& row:out.orientation_intervals)add_known(row);
    for(const auto& row:out.orientation_long_gap_intervals)add_known(row);
    for(const auto& row:out.orientation_sparse_intervals)add_known(row);
    for(const auto& row:out.orientation_short_gap_intervals)add_known(row);
    for(const auto& row:out.orientation_mixed_gap_intervals)add_known(row);
    for(const auto& row:out.orientation_sparse_gap_intervals)add_known(row);
    merge(known_forward);merge(known_reverse);
    detail::VisitUnifiedOrientationEnds(pivots,labels,query,options.count_inferred_gap,
      forward,reverse,known_forward,known_reverse,occupied,&budget,[&](const auto& scan) {
        auto row=scan;const auto bytes=static_cast<std::uint64_t>(row.reason.capacity())+1;
        budget.Account(bytes);unified_reason_bytes+=bytes;append(out.orientation_unified_proposals,std::move(row));
      },[&](const auto& row){append(out.orientation_unified_intervals,row);},[&]{Check(options);});
    // Old accepted tables remain frozen; only the final unknown ledger changes.
    for(const auto& row:out.orientation_unified_intervals)occupied.push_back({row.begin,row.end});
    merge(occupied);
    out.orientation_unified_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-unified_start).count();
    if(options.review_local_competitive_gaps) {
      const auto review_start=std::chrono::steady_clock::now();
      for(const auto& row:out.orientation_unified_intervals)add_known(row);
      merge(known_forward);merge(known_reverse);
      detail::VisitLocalCompetitiveGaps(pivots,labels,query,options.count_inferred_gap,
        forward,reverse,known_forward,known_reverse,&budget,[&](const auto& review) {
          auto row=review;const auto bytes=static_cast<std::uint64_t>(row.reason.capacity())+1;
          budget.Account(bytes);local_review_reason_bytes+=bytes;
          append(out.orientation_local_reviews,std::move(row));
        },[&]{Check(options);});
      // Frozen input snapshots above are never updated by accepted reviews.
      for(const auto& row:out.orientation_local_reviews)if(row.accepted)occupied.push_back({row.begin,row.end});
      merge(occupied);
      out.orientation_local_review_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-review_start).count();
    }
    if(options.review_isolated_opposite_anchors) {
      const auto isolated_start=std::chrono::steady_clock::now();
      for(const auto& row:out.orientation_local_reviews)if(row.accepted)add_known(row);
      merge(known_forward);merge(known_reverse);
      detail::VisitIsolatedOppositeGaps(nodes,neutral,pivots,labels,query,options.count_inferred_gap,
        options.count_cluster_gap,known_forward,known_reverse,&budget,[&](const auto& review) {
          auto row=review;const auto bytes=static_cast<std::uint64_t>(row.reason.capacity())+1;
          budget.Account(bytes);isolated_reason_bytes+=bytes;append(out.orientation_isolated_reviews,std::move(row));
        },[&](const auto& row){append(out.orientation_opposite_reviews,row);},[&]{Check(options);});
      for(const auto& row:out.orientation_isolated_reviews)if(row.accepted)occupied.push_back({row.begin,row.end});
      merge(occupied);
      out.orientation_isolated_review_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-isolated_start).count();
    }
    if(options.review_collinear_fragment_connections) {
      const auto begin=std::chrono::steady_clock::now();
      for(const auto& row:out.orientation_isolated_reviews)if(row.accepted)add_known(row);
      merge(known_forward);merge(known_reverse);
      const auto emit_old=[&](const auto& review) {
          auto row=review;const auto bytes=static_cast<std::uint64_t>(row.reason.capacity())+1;
          budget.Account(bytes);collinear_reason_bytes+=bytes;append(out.orientation_collinear_reviews,std::move(row));
        };
      const auto emit_counter=[&](const auto& row){append(out.orientation_collinear_counters,row);};
      const auto emit_witness=[&](const auto& row){append(out.orientation_collinear_witnesses,row);};
      if(options.review_alternate_original_endpoints) {
        detail::VisitCollinearSameDirectionGapsWithAlternateEndpoints(nodes,neutral,pivots,labels,query,
          options.count_inferred_gap,options.count_cluster_gap,known_forward,known_reverse,&budget,
          emit_old,emit_counter,emit_witness,[&]{Check(options);},options.review_collinear_window_gaps,
          [&](const auto& review) {
            auto row=review;const auto bytes=static_cast<std::uint64_t>(row.reason.capacity())+1;
            budget.Account(bytes);alternate_reason_bytes+=bytes;append(out.orientation_alternate_reviews,std::move(row));
          },[&](const auto& row){append(out.orientation_alternate_counters,row);},
          [&](const auto& row){append(out.orientation_alternate_witnesses,row);},options.review_neutral_window_gaps);
      } else {
        detail::VisitCollinearSameDirectionGaps(nodes,neutral,pivots,labels,query,options.count_inferred_gap,
          options.count_cluster_gap,known_forward,known_reverse,&budget,emit_old,emit_counter,emit_witness,
          [&]{Check(options);},options.review_collinear_window_gaps);
      }
      for(const auto& row:out.orientation_collinear_reviews)if(row.accepted)occupied.push_back({row.begin,row.end});
      for(const auto& row:out.orientation_alternate_reviews)if(row.accepted)occupied.push_back({row.begin,row.end});
      merge(occupied);
      out.orientation_collinear_review_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();
    }
  }
  cursor=0;
  for(const auto& row:occupied){if(cursor<row.begin)append(out.orientation_unknown,Interval{cursor,row.begin});cursor=row.end;}
  if(cursor<qlength)append(out.orientation_unknown,Interval{cursor,qlength});
  out.orientation_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
  out.orientation_capacity_bytes=out.orientation_windows.capacity()*sizeof(CentromereOrientationWindow)+
    out.orientation_intervals.capacity()*sizeof(CentromereOrientationInterval)+out.orientation_sources.capacity()*sizeof(CentromereOrientationSource)+
    out.orientation_gaps.capacity()*sizeof(CentromereOrientationGap)+out.orientation_unknown.capacity()*sizeof(Interval)+
    out.orientation_block_audits.capacity()*sizeof(CentromereOrientationBlockAudit)+
    out.orientation_known_intervals.capacity()*sizeof(CentromereOrientationKnownInterval)+
    out.orientation_long_gap_intervals.capacity()*sizeof(CentromereOrientationInterval)+
    out.orientation_long_gap_gaps.capacity()*sizeof(CentromereOrientationGap)+
    out.orientation_sparse_intervals.capacity()*sizeof(CentromereSparseOrientationInterval)+
    out.orientation_sparse_gaps.capacity()*sizeof(CentromereSparseOrientationGap)+
    out.orientation_sparse_sources.capacity()*sizeof(CentromereOrientationSource)+
    out.orientation_short_gap_intervals.capacity()*sizeof(CentromereOrientationInterval)+
    out.orientation_short_gap_gaps.capacity()*sizeof(CentromereOrientationGap)+
    out.orientation_mixed_gap_intervals.capacity()*sizeof(CentromereOrientationInterval)+
    out.orientation_mixed_gap_gaps.capacity()*sizeof(CentromereOrientationGap)+
    out.orientation_sparse_gap_intervals.capacity()*sizeof(CentromereSparseGapEnd)+
    out.orientation_sparse_gap_gaps.capacity()*sizeof(CentromereSparseGapEnd)+sparse_gap_reason_bytes+
    out.orientation_unified_pivots.capacity()*sizeof(CentromereUnifiedPivot)+
    out.orientation_unified_proposals.capacity()*sizeof(CentromereUnifiedScan)+
    out.orientation_unified_intervals.capacity()*sizeof(CentromereUnifiedInterval)+unified_reason_bytes+
    out.orientation_local_reviews.capacity()*sizeof(CentromereLocalGapReview)+local_review_reason_bytes+
    out.orientation_isolated_reviews.capacity()*sizeof(CentromereLocalGapReview)+isolated_reason_bytes+
    out.orientation_collinear_reviews.capacity()*sizeof(CentromereCollinearReview)+collinear_reason_bytes+
    out.orientation_collinear_counters.capacity()*sizeof(CentromereOppositeReview)+
    out.orientation_collinear_witnesses.capacity()*sizeof(CentromereCollinearWitness)+
    out.orientation_alternate_reviews.capacity()*sizeof(CentromereCollinearReview)+alternate_reason_bytes+
    out.orientation_alternate_counters.capacity()*sizeof(CentromereOppositeReview)+
    out.orientation_alternate_witnesses.capacity()*sizeof(CentromereCollinearWitness)+
    out.orientation_opposite_reviews.capacity()*sizeof(CentromereOppositeReview);
  out.planner_reserved_bytes=budget.peak;out.planner_capacity_bytes+=out.orientation_capacity_bytes;
}

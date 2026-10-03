// Bounded deterministic coverage experiments over existing packed CIGARs.
struct FragmentParentIndex {
    struct Entry {Length b,e,maximum;size_t parent;};
    std::map<SequenceId,std::vector<Entry>> rows;
    static Length build(std::vector<Entry>& v,size_t b,size_t e) {
        if(b==e){return 0;}const size_t m=b+(e-b)/2;
        return v[m].maximum=std::max({v[m].e,build(v,b,m),build(v,m+1,e)});
    }
    FragmentParentIndex(const pairwise_detail::AnchorVec& a,bool reference) {
        for(size_t p=0;p<a.size();++p){const auto& x=a[p];const auto b=reference?x.reference_begin:x.query_begin;
            rows[reference?x.reference_id:x.query_id].push_back({b,b+(reference?x.ref_len:x.qry_len),0,p});}
        for(auto& [id,v]:rows){(void)id;std::sort(v.begin(),v.end(),[](const auto& a0,const auto& b0){return std::tie(a0.b,a0.e,a0.parent)<std::tie(b0.b,b0.e,b0.parent);});build(v,0,v.size());}
    }
    static void search(const std::vector<Entry>& v,size_t b,size_t e,Length lo,Length hi,std::set<size_t>& out) {
        if(b==e||out.size()>256){return;}const auto m=b+(e-b)/2;
        if(v[m].maximum<=lo)return;
        search(v,b,m,lo,hi,out);
        if(out.size()>256)return;
        if(v[m].b<hi){if(v[m].e>lo)out.insert(v[m].parent);search(v,m+1,e,lo,hi,out);}
    }
    void query(SequenceId id,Length b,Length e,std::set<size_t>& out) const {
        const auto it=rows.find(id);if(it!=rows.end())search(it->second,0,it->second.size(),b,e,out);
    }
    std::uint64_t bytes() const {std::uint64_t n=0;for(const auto& [id,v]:rows){(void)id;n+=v.capacity()*sizeof(Entry);}return n;}
};
struct FragmentOwnerIndex {
    struct Owner {Length end;size_t id;};
    std::map<SequenceId,std::map<Length,Owner>> rows;
    void add(SequenceId id,Length b,Length e,size_t owner) {
        auto& m=rows[id];auto it=m.lower_bound(b);
        if((it!=m.end()&&it->first<e)||(it!=m.begin()&&std::prev(it)->second.end>b))
            throw AlignmentError("fragment ownership overlap");
        m.emplace_hint(it,b,Owner{e,owner});
    }
    void remove(SequenceId id,Length b,size_t owner) {
        auto& m=rows.at(id);auto it=m.find(b);
        if(it==m.end()||it->second.id!=owner)throw AlignmentError("fragment ownership removal mismatch");
        m.erase(it);
    }
    void query(SequenceId id,Length b,Length e,std::set<size_t>& out) const {
        auto d=rows.find(id);if(d==rows.end())return;const auto& m=d->second;
        auto it=m.lower_bound(b);if(it!=m.begin()){auto p=std::prev(it);if(p->second.end>b)it=p;}
        for(;it!=m.end()&&it->first<e&&out.size()<=4;++it)if(it->second.end>b)out.insert(it->second.id);
    }
};
struct FragmentTotals {
    Length reference{},query{},paired{},matches{};
    void add(const RecoveryFragment& f) {
        const auto checked=[](Length a,Length b){if(b>UINT64_MAX-a)throw AlignmentError("fragment total overflow");return a+b;};
        const auto r=f.re-f.rb,q=f.qe-f.qb,sum=checked(r,q);
        if(sum<f.columns)throw AlignmentError("fragment invalid column total");
        reference=checked(reference,r);query=checked(query,q);paired=checked(paired,sum-f.columns);matches=checked(matches,f.matches);
    }
};
using FragmentRegions=std::map<SequenceId,std::vector<std::pair<Length,Length>>>;
void MergeFragmentRegions(FragmentRegions& regions) {
    for(auto& [id,v]:regions){(void)id;std::sort(v.begin(),v.end());size_t n=0;
        for(auto x:v){if(n&&x.first<=v[n-1].second)v[n-1].second=std::max(v[n-1].second,x.second);else v[n++]=x;}v.resize(n);}
}
template<class Emit>
void ScanShortExact(const pairwise_detail::Anchor& a,size_t parent,const Sequences& refs,const Sequences& qrys,
                    const OccupiedIntervals& ro,const OccupiedIntervals& qo,const PairwiseCoreOptions& o,Emit emit) {
    const auto& r=refs.at(a.reference_id)->bases;const auto& q=qrys.at(a.query_id)->bases;
    Length rp=a.reference_begin,qp=a.strand==Strand::Forward?a.query_begin:a.query_begin+a.qry_len,col=0;
    RecoveryFragment f{};bool open=false;
    const auto finish=[&] {
        if(open&&f.matches>=o.gap_fill_min_match&&f.matches<o.min_cluster) {
            f.parent=parent;f.column_end=col;f.re=rp;
            if(a.strand==Strand::Forward)f.qe=qp;else{f.qe=f.qb;f.qb=qp;}
            f.columns=f.matches;emit(f);
        }open=false;f={};
    };
    const auto advance=[&](unsigned op,Length count){if(op!=1)rp+=count;if(op!=2){if(a.strand==Strand::Forward)qp+=count;else qp-=count;}col+=count;};
    for(auto u:a.cigar){const unsigned op=u&15;Length left=u>>4;
        if(op==1||op==2){finish();advance(op,left);continue;}
        while(left) {
            KernelCheck(o,"fragment-short-exact");
            const auto rr=ro.run(a.reference_id,rp,left,false),qq=qo.run(a.query_id,qp,left,a.strand==Strand::Reverse);
            const auto count=std::min(rr.second,qq.second);if(!count)throw AlignmentError("short exact cursor failed to advance");
            if(rr.first||qq.first){finish();advance(op,count);left-=count;continue;}
            for(Length k=0;k<count;++k){
                if((k&65535)==0)KernelCheck(o,"fragment-short-exact");
                const char x=r.at(rp),y=a.strand==Strand::Forward?q.at(qp):Complement(q.at(qp-1));
                if(x!='N'&&x==y){if(!open){open=true;f.column_begin=col;f.rb=rp;f.qb=qp;}++f.matches;}else finish();
                advance(op,1);
            }left-=count;
        }
    }finish();
}
void RefineFragmentSelection(const pairwise_detail::AnchorVec& anchors,const Sequences& refs,const Sequences& qrys,
        const PairwiseCoreOptions& o,const std::vector<double>& identities,const std::vector<Length>& parent_matches,
        std::vector<RecoveryFragment>& selected,OccupiedIntervals& ro,OccupiedIntervals& qo,PairwiseStatistics& stats) {
    const auto started=KernelClock::now();
    if(o.fragment_local_exchange) {
        FragmentParentIndex ri(anchors,true),qi(anchors,false);FragmentOwnerIndex owners_r,owners_q;
        std::vector<unsigned char> active(selected.size(),1);
        auto add_owner=[&](size_t id){const auto& f=selected[id];const auto& a=anchors[f.parent];owners_r.add(a.reference_id,f.rb,f.re,id);owners_q.add(a.query_id,f.qb,f.qe,id);};
        for(size_t i=0;i<selected.size();++i)add_owner(i);
        if(o.resident_bytes_callback)stats.fragment_peak_rss_bytes=std::max(stats.fragment_peak_rss_bytes,o.resident_bytes_callback());
        auto remove=[&](const RecoveryFragment& f){const auto& a=anchors[f.parent];ro.remove(a.reference_id,f.rb,f.re);qo.remove(a.query_id,f.qb,f.qe);};
        auto occupy=[&](const RecoveryFragment& f){const auto& a=anchors[f.parent];ro.add(a.reference_id,f.rb,f.re);qo.add(a.query_id,f.qb,f.qe);};
        const OccupiedIntervals empty;
        for(unsigned pass=0;pass<2;++pass) {
            const auto before=stats.exchange_accepted;++stats.exchange_passes;
            for(size_t parent=0;parent<anchors.size();++parent) {
                KernelCheck(o,"fragment-local-exchange");
                const auto& a=anchors[parent];
                // Complete selected parents cannot improve their own placement.
                std::set<size_t> current;owners_r.query(a.reference_id,a.reference_begin,a.reference_begin+a.ref_len,current);
                bool whole=false;for(auto id:current){const auto& f=selected[id];whole|=f.parent==parent&&f.column_begin==0&&f.column_end==a.alignment_length;}
                if(whole)continue;
                std::vector<RecoveryFragment> challenges;
                ScanRecoveryFragments(a,parent,0,a.alignment_length,refs,qrys,empty,empty,o.min_cluster,o,[&](auto f){f.score=static_cast<double>(f.matches)*identities[parent];challenges.push_back(f);});
                for(const auto& challenge:challenges) {
                    std::set<size_t> blocked;owners_r.query(a.reference_id,challenge.rb,challenge.re,blocked);owners_q.query(a.query_id,challenge.qb,challenge.qe,blocked);
                    if(blocked.empty())continue;
                    if(blocked.size()>4){++stats.exchange_blocker_limit;if(pairwise_detail::diagnostic_group_stream)*pairwise_detail::diagnostic_group_stream<<"blocker-cap\t"<<parent<<'\t'<<challenge.column_begin<<'\t'<<challenge.column_end<<'\t'<<challenge.rb<<'\t'<<challenge.re<<'\t'<<challenge.qb<<'\t'<<challenge.qe<<'\n';continue;}
                    if(blocked.size()==1){const auto& b=selected[*blocked.begin()];if(b.parent==parent&&b.column_begin<=challenge.column_begin&&b.column_end>=challenge.column_end)continue;}
                    FragmentRegions rr,qr;auto region=[&](const RecoveryFragment& f){const auto& p=anchors[f.parent];rr[p.reference_id].emplace_back(f.rb,f.re);qr[p.query_id].emplace_back(f.qb,f.qe);};
                    region(challenge);for(auto id:blocked)region(selected[id]);MergeFragmentRegions(rr);MergeFragmentRegions(qr);
                    std::set<size_t> pool{parent};for(auto id:blocked)pool.insert(selected[id].parent);
                    for(const auto& [id,v]:rr)for(auto [b,e]:v)ri.query(id,b,e,pool);
                    for(const auto& [id,v]:qr)for(auto [b,e]:v)qi.query(id,b,e,pool);
                    if(pool.size()>256){++stats.exchange_pool_limit;if(pairwise_detail::diagnostic_group_stream)*pairwise_detail::diagnostic_group_stream<<"pool-cap\t"<<parent<<'\t'<<challenge.column_begin<<'\t'<<challenge.column_end<<'\t'<<challenge.rb<<'\t'<<challenge.re<<'\t'<<challenge.qb<<'\t'<<challenge.qe<<'\n';continue;}
                    ++stats.exchange_trials;FragmentTotals old_total,new_total;
                    for(auto id:blocked){old_total.add(selected[id]);remove(selected[id]);}
                    std::vector<RecoveryFragment> proposal{challenge};occupy(challenge);new_total.add(challenge);
                    std::priority_queue<RecoveryFragment,std::vector<RecoveryFragment>,FragmentHigher> heap;
                    for(auto p:pool){const auto& x=anchors[p];
                        ScanRecoveryFragments(x,p,0,x.alignment_length,refs,qrys,ro,qo,o.min_cluster,o,[&](auto f){f.score=static_cast<double>(f.matches)*identities[p];heap.push(f);});
                    }
                    while(!heap.empty()){KernelCheck(o,"fragment-local-exchange");auto f=heap.top();heap.pop();const auto& x=anchors[f.parent];
                        if(!GapIntervalFree(ro,x.reference_id,f.rb,f.re)||!GapIntervalFree(qo,x.query_id,f.qb,f.qe)) {
                            ScanRecoveryFragments(x,f.parent,f.column_begin,f.column_end,refs,qrys,ro,qo,o.min_cluster,o,[&](auto child){child.score=static_cast<double>(child.matches)*identities[child.parent];heap.push(child);});continue;
                        }occupy(f);proposal.push_back(f);new_total.add(f);
                    }
                    const bool accept=new_total.reference>=old_total.reference&&new_total.query>=old_total.query&&new_total.paired>old_total.paired&&new_total.matches>=old_total.matches;
                    if(pairwise_detail::diagnostic_group_stream) {
                        auto& log=*pairwise_detail::diagnostic_group_stream;
                        log<<"exchange\t"<<pass<<'\t'<<parent<<'\t'<<challenge.column_begin<<'\t'<<challenge.column_end<<'\t'<<accept<<'\t'<<old_total.paired<<'\t'<<new_total.paired;
                        for(auto id:blocked)log<<"\tremoved:"<<selected[id].parent<<':'<<selected[id].column_begin<<':'<<selected[id].column_end;
                        for(const auto& f:proposal)log<<"\tproposed:"<<f.parent<<':'<<f.column_begin<<':'<<f.column_end;
                        log<<'\n';
                    }
                    if(accept){
                        for(auto id:blocked){const auto& f=selected[id];const auto& x=anchors[f.parent];owners_r.remove(x.reference_id,f.rb,id);owners_q.remove(x.query_id,f.qb,id);active[id]=0;}
                        for(const auto& f:proposal){selected.push_back(f);active.push_back(1);add_owner(selected.size()-1);}
                        ++stats.exchange_accepted;stats.exchange_paired_gain+=new_total.paired-old_total.paired;
                    }else{for(const auto& f:proposal)remove(f);for(auto id:blocked)occupy(selected[id]);++stats.exchange_objective_rejected;}
                    stats.refinement_auxiliary_peak_bytes=std::max(stats.refinement_auxiliary_peak_bytes,
                        ri.bytes()+qi.bytes()+selected.capacity()*sizeof(RecoveryFragment)+active.capacity()+
                        proposal.capacity()*sizeof(RecoveryFragment));
                }
            }
            if(stats.exchange_accepted==before)break;
        }
        size_t n=0;for(size_t i=0;i<selected.size();++i)if(active[i])selected[n++]=selected[i];selected.resize(n);
    }
    if(o.fragment_short_exact) {
        std::priority_queue<RecoveryFragment,std::vector<RecoveryFragment>,FragmentHigher> heap;
        for(size_t p=0;p<anchors.size();++p){KernelCheck(o,"fragment-short-exact");if(parent_matches[p]<o.min_cluster)continue;
            ScanShortExact(anchors[p],p,refs,qrys,ro,qo,o,[&](auto f){f.score=static_cast<double>(f.matches)*identities[p];heap.push(f);++stats.short_exact_proposed;});}
        while(!heap.empty()){auto f=heap.top();heap.pop();const auto& a=anchors[f.parent];
            if(!GapIntervalFree(ro,a.reference_id,f.rb,f.re)||!GapIntervalFree(qo,a.query_id,f.qb,f.qe)){
                // Clip only this original exact run; never rediscover the whole parent.
                ScanRecoveryFragments(a,f.parent,f.column_begin,f.column_end,refs,qrys,ro,qo,o.gap_fill_min_match,o,[&](auto c){
                    if(c.matches==c.columns&&c.matches<o.min_cluster){c.score=static_cast<double>(c.matches)*identities[c.parent];heap.push(c);}});continue;}
            ro.add(a.reference_id,f.rb,f.re);qo.add(a.query_id,f.qb,f.qe);selected.push_back(f);
            ++stats.short_exact_accepted;stats.short_exact_bases+=f.matches;
            if(pairwise_detail::diagnostic_group_stream)*pairwise_detail::diagnostic_group_stream<<"short-exact\t"<<f.parent<<'\t'<<f.column_begin<<'\t'<<f.column_end<<'\n';
        }
    }
    stats.fragment_refinement_seconds=Elapsed(started);
}

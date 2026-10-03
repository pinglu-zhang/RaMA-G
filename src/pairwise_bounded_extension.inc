// Bounded anchored extension. No new seed search; budgets are reserved before
// dispatch and never refunded based on runtime or completion order.
struct BoundedEnd {
    size_t parent{};
    SequenceId rid{},qid{};
    Length rb{},qb{},rlen{},qlen{}; // rb/qb are forward/oriented boundary cursors
    Strand strand{};
    bool left{};
    double quality{};
    Length support{};
    std::uint64_t work{};
};
std::uint64_t BoundedBudget(const Sequences& refs,const Sequences& qrys,const PairwiseCoreOptions& o) {
    const auto total=[](const Sequences& seqs){
        std::uint64_t n=0;
        for(const auto& [id,s]:seqs){(void)id;if(s->size()>UINT64_MAX-n)throw AlignmentError("bounded input length overflow");n+=s->size();}
        return n;
    };
    const auto n=std::min(total(refs),total(qrys));
    return std::min(o.bounded_work_budget,n>UINT64_MAX/2?UINT64_MAX:2*n);
}
char BoundedQueryBase(const std::string& q,Length p,Strand strand) {
    return strand==Strand::Forward?q.at(p):Complement(q.at(q.size()-1-p));
}
BoundedEnd MakeBoundedEnd(const pairwise_detail::Anchor& a,size_t parent,bool left,
        const Sequences& refs,const Sequences& qrys,const PairwiseCoreOptions& o,
        const OccupiedIntervals* ro=nullptr,const OccupiedIntervals* qo=nullptr) {
    const auto& r=refs.at(a.reference_id)->bases;const auto& q=qrys.at(a.query_id)->bases;
    if(a.reference_begin>r.size()||a.ref_len>r.size()-a.reference_begin||
       a.query_begin>q.size()||a.qry_len>q.size()-a.query_begin)throw AlignmentError("bounded anchor outside sequence");
    const auto oriented=a.strand==Strand::Forward?a.query_begin:q.size()-a.query_begin-a.qry_len;
    BoundedEnd d;d.parent=parent;d.rid=a.reference_id;d.qid=a.query_id;d.strand=a.strand;d.left=left;
    d.rb=a.reference_begin+(left?0:a.ref_len);d.qb=oriented+(left?0:a.qry_len);
    const auto rmax=std::min(o.bounded_extension_length,left?d.rb:r.size()-d.rb);
    const auto qmax=std::min(o.bounded_extension_length,left?d.qb:q.size()-d.qb);
    while(d.rlen<rmax) {
        const auto p=left?d.rb-1-d.rlen:d.rb+d.rlen;
        if(r[p]=='N'||(ro&&!GapIntervalFree(*ro,d.rid,p,p+1)))break;
        ++d.rlen;
    }
    while(d.qlen<qmax) {
        const auto p=left?d.qb-1-d.qlen:d.qb+d.qlen;
        const auto forward=a.strand==Strand::Forward?p:q.size()-1-p;
        if(BoundedQueryBase(q,p,a.strand)=='N'||(qo&&!GapIntervalFree(*qo,d.qid,forward,forward+1)))break;
        ++d.qlen;
    }
    if(d.rlen&&d.qlen)d.work=3*(d.rlen+1)*(d.qlen+1); // lengths checked <=500
    return d;
}
bool BoundedPriority(const BoundedEnd& a,const BoundedEnd& b) {
    if(a.quality!=b.quality)return a.quality>b.quality;
    if(a.support!=b.support)return a.support>b.support;
    if(a.work!=b.work)return a.work<b.work;
    return std::tie(a.rid,a.rb,a.qid,a.qb,a.strand,a.parent,a.left)<
           std::tie(b.rid,b.rb,b.qid,b.qb,b.strand,b.parent,b.left);
}
void ReserveBounded(std::vector<BoundedEnd>& tasks,std::uint64_t budget,BoundedExtensionStatistics& stats) {
    std::stable_sort(tasks.begin(),tasks.end(),BoundedPriority);
    stats.candidates=tasks.size();size_t kept=0;
    for(const auto& t:tasks) {
        if(t.work>budget-stats.reserved_cells){++stats.budget_skipped;continue;}
        stats.reserved_cells+=t.work;tasks[kept++]=t;
    }
    tasks.resize(kept);stats.scheduled=kept;
    stats.auxiliary_capacity_bytes=tasks.capacity()*sizeof(BoundedEnd);
}
// Returns a partial endpoint CIGAR in increasing-reference orientation. Left
// extension reverses both strings (not their complement) and reverses operations.
std::optional<pairwise_detail::Anchor> RunBoundedEnd(const BoundedEnd& d,
        const Sequences& refs,const Sequences& qrys,const PairwiseCoreOptions& o) {
    KernelCheck(o,"bounded-endpoint");
    const auto& r=refs.at(d.rid)->bases;const auto& q=qrys.at(d.qid)->bases;
    std::string rs,qs;rs.reserve(d.rlen);qs.reserve(d.qlen);
    for(Length k=0;k<d.rlen;++k)rs.push_back(r[d.left?d.rb-1-k:d.rb+k]);
    for(Length k=0;k<d.qlen;++k)qs.push_back(BoundedQueryBase(q,d.left?d.qb-1-k:d.qb+k,d.strand));
    // Same endpoint Z-drop as the production connector (2 * break length 200).
    auto aligned=pairwise_detail::extendAlignKSW2Result({rs,false},{qs,false},400);
    const auto summary=pairwise_detail::summarizeCigar(aligned.cigar);
    if(summary.reference_length>d.rlen||summary.query_length>d.qlen)
        throw AlignmentError("bounded endpoint consumption exceeds window");
    if(!summary.reference_length||!summary.query_length)return std::nullopt;
    if(d.left)std::reverse(aligned.cigar.begin(),aligned.cigar.end());
    const auto rb=d.left?d.rb-summary.reference_length:d.rb;
    const auto oq=d.left?d.qb-summary.query_length:d.qb;
    const auto qb=d.strand==Strand::Forward?oq:q.size()-oq-summary.query_length;
    pairwise_detail::Anchor a(d.rid,rb,summary.reference_length,d.qid,qb,summary.query_length,
        d.strand,summary.alignment_length,summary.match_length,std::move(aligned.cigar));
    (void)ConvertAnchor(a,refs,qrys,false);
    return a;
}
void ObserveBounded(BoundedExtensionStatistics& stats,const PairwiseCoreOptions& o) {
    if(o.resident_bytes_callback)stats.peak_rss_bytes=std::max(stats.peak_rss_bytes,o.resident_bytes_callback());
}
template<class Work,class Reset>
void BoundedBatch(size_t count,PairwiseCoreOptions options,Work work,Reset reset) {
    options.threads=std::min({options.threads,options.bounded_workers,4U});
    for(;;) {
        try {ParallelKernel(count,options,work);return;}
        catch(const std::bad_alloc&) {
            // All worker exceptions have joined before ParallelKernel rethrows.
            // Retry the identical frozen tasks with less simultaneous workspace.
            reset();if(options.threads==1)throw;options.threads=std::max(1U,options.threads/2);
        }
    }
}
// S descriptors reference compact per-group storage. No per-chain CIGAR exists
// until its endpoint work has been admitted to the frozen budget.
pairwise_detail::AnchorVec ExtendShortChains(std::vector<pairwise_detail::ShortChainStore>& stores,
        const pairwise_detail::SequenceViews::ManagerVariant& rm,
        const pairwise_detail::SequenceViews::ManagerVariant& qm,
        const Sequences& refs,const Sequences& qrys,const PairwiseCoreOptions& o,
        BoundedExtensionStatistics& stats) {
    const auto started=KernelClock::now();
    struct Chain {size_t group{},slice{};bool left{},right{};};
    std::vector<Chain> chains;std::vector<BoundedEnd> tasks;
    for(size_t g=0;g<stores.size();++g)for(size_t j=0;j<stores[g].slices.size();++j) {
        const auto& s=stores[g].slices[j];const auto& first=stores[g].seeds[s.begin];
        Length rb=UINT64_MAX,re=0,qb=UINT64_MAX,qe=0;
        for(size_t k=s.begin;k<s.end;++k){const auto& seed=stores[g].seeds[k];rb=std::min(rb,seed.reference_begin);re=std::max(re,seed.reference_begin+seed.length);qb=std::min(qb,seed.query_begin);qe=std::max(qe,seed.query_begin+seed.length);}
        pairwise_detail::Anchor a(first.reference_id,rb,re-rb,first.query_id,qb,qe-qb,first.strand,0,0,{});
        const auto id=chains.size();chains.push_back({g,j,false,false});
        for(bool left:{true,false}){auto t=MakeBoundedEnd(a,id,left,refs,qrys,o);t.quality=1;t.support=s.support;if(t.work)tasks.push_back(t);}
    }
    ReserveBounded(tasks,BoundedBudget(refs,qrys,o),stats);
    stats.auxiliary_capacity_bytes+=chains.capacity()*sizeof(Chain);
    for(const auto& s:stores)stats.auxiliary_capacity_bytes+=s.seeds.capacity()*sizeof(Seed)+s.slices.capacity()*sizeof(pairwise_detail::ShortChainStore::Slice);
    for(const auto& t:tasks){if(t.left)chains[t.parent].left=true;else chains[t.parent].right=true;}
    std::vector<BoundedEnd>().swap(tasks);
    pairwise_detail::AnchorVec output;
    auto worker_options=o;worker_options.threads=std::min({o.threads,o.bounded_workers,4U});
    // Group order allows consumed seed storage to be released immediately.
    for(size_t b=0;b<chains.size();b+=64) {
        const auto n=std::min<size_t>(64,chains.size()-b);
        std::vector<std::optional<pairwise_detail::Anchor>> results(n);
        std::vector<std::uint64_t> calls(n);
        BoundedBatch(n,worker_options,[&](size_t i){
            const auto& c=chains[b+i];if(!c.left&&!c.right)return;
            auto& store=stores[c.group];const auto s=store.slices[c.slice];
            pairwise_detail::MatchCluster chain(store.seeds.begin()+s.begin,store.seeds.begin()+s.end);
            auto a=pairwise_detail::extendClusterToAnchor(chain,rm,qm);
            const auto left_request=MakeBoundedEnd(a,b+i,true,refs,qrys,o);
            const auto right_request=MakeBoundedEnd(a,b+i,false,refs,qrys,o);
            std::optional<pairwise_detail::Anchor> left,right;
            if(c.left){left=RunBoundedEnd(left_request,refs,qrys,o);++calls[i];}
            if(c.right){right=RunBoundedEnd(right_request,refs,qrys,o);++calls[i];}
            pairwise_detail::Cigar_t cigar;
            if(left)pairwise_detail::appendCigarMove(cigar,left->cigar);
            pairwise_detail::appendCigarMove(cigar,a.cigar);
            if(right)pairwise_detail::appendCigarMove(cigar,right->cigar);
            if(left){a.reference_begin-=left->ref_len;a.ref_len+=left->ref_len;a.qry_len+=left->qry_len;if(a.strand==Strand::Forward)a.query_begin-=left->qry_len;}
            if(right){a.ref_len+=right->ref_len;a.qry_len+=right->qry_len;if(a.strand==Strand::Reverse)a.query_begin-=right->qry_len;}
            a.cigar=std::move(cigar);const auto sum=pairwise_detail::summarizeCigar(a.cigar);
            a.alignment_length=sum.alignment_length;a.aligned_base=sum.match_length;a.candidate_source=3;
            const auto record=ConvertAnchor(a,refs,qrys,false);
            const auto matches=a.alignment_length-record.edit_distance;
            if(matches>=o.bounded_output_min_support&&static_cast<double>(matches)/a.alignment_length>=.95)results[i]=std::move(a);
        },[&]{for(auto& a:results)a.reset();std::fill(calls.begin(),calls.end(),0);});
        for(size_t i=0;i<n;++i){stats.endpoint_calls+=calls[i];if(results[i]){
            ++stats.accepted;const auto& a=*results[i];
            const auto record=ConvertAnchor(a,refs,qrys,false);stats.matches+=a.alignment_length-record.edit_distance;
            for(auto u:a.cigar)if((u&15)!=1&&(u&15)!=2)stats.paired_columns+=u>>4;
            output.push_back(std::move(*results[i]));
        }else if(chains[b+i].left||chains[b+i].right)++stats.quality_rejected;}
        const auto last=chains[b+n-1].group;
        const auto next=b+n<chains.size()?chains[b+n].group:stores.size();
        for(size_t g=chains[b].group;g<=last;++g)if(g<next)stores[g]=pairwise_detail::ShortChainStore{};
        ObserveBounded(stats,o);
    }
    stores.clear();stores.shrink_to_fit();stats.seconds=Elapsed(started);return output;
}
pairwise_detail::AnchorVec ExtendSelectedEnds(const pairwise_detail::AnchorVec& anchors,
        const Sequences& refs,const Sequences& qrys,const PairwiseCoreOptions& o,
        BoundedExtensionStatistics& stats) {
    const auto started=KernelClock::now();OccupiedIntervals ro,qo;
    auto flank_options=o;flank_options.gap_fill_min_match=o.bounded_endpoint_min_exact;
    for(const auto& a:anchors)if(a.ref_selected&&a.qry_selected){ro.add(a.reference_id,a.reference_begin,a.reference_begin+a.ref_len);qo.add(a.query_id,a.query_begin,a.query_begin+a.qry_len);}
    std::vector<BoundedEnd> tasks;
    for(size_t i=0;i<anchors.size();++i) {
        const auto& a=anchors[i];if(!a.ref_selected||!a.qry_selected)continue;
        const auto info=InspectGapFlank(a,refs,qrys,flank_options);
        if(info.matches<o.bounded_output_min_support||!info.columns)continue;
        const double identity=static_cast<double>(info.matches)/info.columns;if(identity<.98)continue;
        for(bool left:{true,false})if(left?info.begin_exact:info.end_exact) {
            auto t=MakeBoundedEnd(a,i,left,refs,qrys,o,&ro,&qo);t.quality=identity;t.support=info.matches;if(t.work)tasks.push_back(t);
        }
    }
    ReserveBounded(tasks,BoundedBudget(refs,qrys,o),stats);
    auto worker_options=o;worker_options.threads=std::min({o.threads,o.bounded_workers,4U});
    pairwise_detail::AnchorVec output;
    for(size_t b=0;b<tasks.size();b+=64) {
        const auto n=std::min<size_t>(64,tasks.size()-b);
        std::vector<std::optional<pairwise_detail::Anchor>> results(n);
        BoundedBatch(n,worker_options,[&](size_t i){results[i]=RunBoundedEnd(tasks[b+i],refs,qrys,o);},[&]{for(auto& a:results)a.reset();});
        for(size_t i=0;i<n;++i) {
            ++stats.endpoint_calls;if(!results[i]){++stats.quality_rejected;continue;}
            const auto& a=*results[i];
            // Trim only terminal gaps; the unoccupied window is immutable.
            OccupiedIntervals empty;std::optional<RecoveryFragment> fragment;
            ScanRecoveryFragments(a,0,0,UINT64_MAX,refs,qrys,empty,empty,1,o,[&](const auto& f){if(fragment)throw AlignmentError("bounded endpoint unexpectedly contains N");fragment=f;});
            if(!fragment){++stats.quality_rejected;continue;}
            const auto& f=*fragment;const bool exact=f.matches==f.columns;
            if((exact&&f.matches<o.bounded_endpoint_min_exact)||(!exact&&(f.matches<o.bounded_output_min_support||static_cast<double>(f.matches)/f.columns<.95))){++stats.quality_rejected;continue;}
            if(!GapIntervalFree(ro,a.reference_id,f.rb,f.re)||!GapIntervalFree(qo,a.query_id,f.qb,f.qe)){++stats.conflict_rejected;continue;}
            auto part=MaterializeRecovery(a,f);part.candidate_source=4;
            (void)ConvertAnchor(part,refs,qrys,false);
            ro.add(part.reference_id,f.rb,f.re);qo.add(part.query_id,f.qb,f.qe);
            ++stats.accepted;stats.matches+=f.matches;
            for(auto u:part.cigar)if((u&15)!=1&&(u&15)!=2)stats.paired_columns+=u>>4;
            output.push_back(std::move(part));
        }
        ObserveBounded(stats,o);
    }
    stats.seconds=Elapsed(started);return output;
}

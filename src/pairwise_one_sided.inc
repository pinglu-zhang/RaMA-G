// Additive, non-recursive recovery; the original selected set is immutable.
// Only XOR-selected parents may contribute, and every emitted column must lie
// inside one frozen, sequence-verified pair of original flanks.
pairwise_detail::AnchorVec RecoverOneSidedGuarded(const pairwise_detail::AnchorVec& anchors,
        const Sequences& refs,const Sequences& qrys,const PairwiseCoreOptions& o,
        PairwiseStatistics& stats) {
    const auto begin=KernelClock::now();
    stats.post_selection_strategy="one-sided-guarded-v1";
    stats.recovery_enabled=true;
    for(const auto& a:anchors) {
        if(a.ref_selected&&a.qry_selected)++stats.post_selection_original_records;
        else if(a.ref_selected!=a.qry_selected)++stats.post_selection_one_sided;
        else ++stats.post_selection_both_rejected;
    }
    OccupiedIntervals ro,qo;PairwiseStatistics gs;
    const auto gaps=CollectGuardedGaps(anchors,anchors.size(),refs,qrys,o,gs,ro,qo);
    stats.post_selection_windows=gaps.size();
    stats.post_selection_window_geometry_rejected=gs.gap_fill_geometry_rejected;
    stats.post_selection_window_occupied_rejected=gs.gap_fill_occupied_rejected;
    stats.post_selection_window_quality_rejected=gs.gap_fill_flank_rejected;
    stats.post_selection_window_n_rejected=gs.gap_fill_n_rejected;
    pairwise_detail::AnchorVec output;
    if(gaps.empty()) {
        stats.post_selection_no_window=stats.post_selection_one_sided;
        stats.recovery_seconds=Elapsed(begin);
        return output;
    }
    using Key=std::tuple<SequenceId,SequenceId,Strand>;
    std::map<Key,std::vector<size_t>> groups;
    for(size_t i=0;i<gaps.size();++i)groups[{gaps[i].rid,gaps[i].qid,gaps[i].strand}].push_back(i);
    std::uint64_t aux=gaps.capacity()*sizeof(GuardedGap);
    for(auto& [key,ids]:groups) {
        (void)key;
        std::sort(ids.begin(),ids.end(),[&](size_t a,size_t b){return gaps[a].qb<gaps[b].qb;});
        for(size_t i=1;i<ids.size();++i)
            if(gaps[ids[i-1]].qe>gaps[ids[i]].qb)throw AlignmentError("overlapping frozen recovery windows");
        aux+=ids.capacity()*sizeof(size_t);
    }
    struct Entry {RecoveryFragment f;size_t gap;};
    struct Higher {
        bool operator()(const Entry& a,const Entry& b) const {
            const RecoveryHigher cmp;
            if(cmp(a.f,b.f))return true;
            if(cmp(b.f,a.f))return false;
            return a.gap>b.gap;
        }
    };
    struct Queue:std::priority_queue<Entry,std::vector<Entry>,Higher> {
        std::uint64_t bytes() const {return this->c.capacity()*sizeof(Entry);}
    } queue;
    const auto observe=[&](bool sample) {
        stats.post_selection_queue_peak=std::max<std::uint64_t>(stats.post_selection_queue_peak,queue.size());
        stats.post_selection_auxiliary_peak_bytes=std::max(stats.post_selection_auxiliary_peak_bytes,aux+queue.bytes());
        if(sample&&o.resident_bytes_callback)
            stats.post_selection_peak_rss_bytes=std::max(stats.post_selection_peak_rss_bytes,o.resident_bytes_callback());
    };
    const auto acceptable=[&](const RecoveryFragment& f,size_t index) {
        const auto& g=gaps[index];
        if(f.rb<g.rb||f.re>g.re||f.qb<g.qb||f.qe>g.qe||!f.columns)
            throw AlignmentError("post-selection escaped its frozen window");
        if(static_cast<double>(f.matches)/f.columns<g.minimum_identity) {
            ++stats.post_selection_identity_rejected;
            return false;
        }
        return true;
    };
    const auto scan=[&](size_t parent,size_t window,Length b,Length e,auto emit) {
        const auto& g=gaps[window];
        const RecoveryBounds bounds{g.rb,g.re,g.qb,g.qe};
        ScanRecoveryFragments(anchors[parent],parent,b,e,refs,qrys,ro,qo,o.min_cluster,o,emit,&bounds,&stats);
    };
    observe(true);
    for(size_t parent=0;parent<anchors.size();++parent) {
        const auto& a=anchors[parent];
        if(a.ref_selected==a.qry_selected)continue;
        KernelCheck(o,"post-selection");
        ++stats.recovery_candidates_checked;
        bool intersected=false;
        const auto found=groups.find({a.reference_id,a.query_id,a.strand});
        if(found!=groups.end()) {
            const auto& ids=found->second;
            auto it=std::lower_bound(ids.begin(),ids.end(),a.query_begin,
                [&](size_t id,Length p){return gaps[id].qe<=p;});
            for(;it!=ids.end()&&gaps[*it].qb<a.query_begin+a.qry_len;++it) {
                const auto& g=gaps[*it];
                if(g.re<=a.reference_begin||g.rb>=a.reference_begin+a.ref_len)continue;
                intersected=true;
                scan(parent,*it,0,UINT64_MAX,[&](const RecoveryFragment& f) {
                    if(!acceptable(f,*it))return;
                    queue.push({f,*it});++stats.recovery_proposed_fragments;observe(false);
                });
            }
        }
        if(intersected)++stats.post_selection_parents_in_windows;
        else ++stats.post_selection_no_window;
        if((parent&65535)==0)observe(true);
    }
    observe(true);
    std::uint64_t popped=0;
    while(!queue.empty()) {
        KernelCheck(o,"post-selection");
        const auto entry=queue.top();queue.pop();const auto& f=entry.f;
        bool unchanged=false;
        scan(f.parent,entry.gap,f.column_begin,f.column_end,[&](const RecoveryFragment& g) {
            if(g.column_begin==f.column_begin&&g.column_end==f.column_end&&
               g.matches==f.matches&&g.rb==f.rb&&g.re==f.re&&g.qb==f.qb&&g.qe==f.qe)unchanged=true;
            else if(acceptable(g,entry.gap)) {
                queue.push({g,entry.gap});++stats.recovery_requeued_fragments;observe(false);
            }
        });
        if(!unchanged){++stats.post_selection_conflict_reclips;continue;}
        const auto& a=anchors[f.parent];
        if(!GapIntervalFree(ro,a.reference_id,f.rb,f.re)||!GapIntervalFree(qo,a.query_id,f.qb,f.qe))
            throw AlignmentError("post-selection accepted overlapping fragment");
        auto restored=MaterializeRecovery(a,f);
        const auto verified=ConvertAnchor(restored,refs,qrys,false);
        Length paired=0,ins=0,del=0;
        for(const auto unit:restored.cigar) {
            if((unit&15)==1)ins+=unit>>4;
            else if((unit&15)==2)del+=unit>>4;
            else paired+=unit>>4;
        }
        if(paired<f.matches||verified.edit_distance!=paired-f.matches+ins+del)
            throw AlignmentError("post-selection match summary disagrees with CIGAR");
        ro.add(a.reference_id,f.rb,f.re);qo.add(a.query_id,f.qb,f.qe);
        stats.recovery_reference_bases+=f.re-f.rb;stats.recovery_query_bases+=f.qe-f.qb;
        stats.post_selection_paired_columns+=paired;stats.post_selection_matches+=f.matches;
        stats.post_selection_mismatches+=paired-f.matches;
        stats.post_selection_insertions+=ins;stats.post_selection_deletions+=del;
        ++stats.recovery_accepted_fragments;output.push_back(std::move(restored));
        if((++popped&1023)==0)observe(true);
    }
    observe(true);
    if(stats.recovery_candidates_checked!=stats.post_selection_one_sided||
       stats.post_selection_parents_in_windows+stats.post_selection_no_window!=stats.post_selection_one_sided)
        throw AlignmentError("post-selection candidate accounting mismatch");
    stats.recovery_seconds=Elapsed(begin);
    return output;
}

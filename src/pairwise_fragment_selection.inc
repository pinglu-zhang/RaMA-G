// Deterministic greedy two-axis packing of existing alignment columns.
// Parent identity is frozen; clipping cannot increase a descriptor's priority.
struct FragmentHigher {
    bool operator()(const RecoveryFragment& a,const RecoveryFragment& b) const {
        if(a.score!=b.score)return a.score<b.score;
        if(a.matches!=b.matches)return a.matches<b.matches;
        if(a.parent!=b.parent)return a.parent>b.parent;
        if(a.column_begin!=b.column_begin)return a.column_begin>b.column_begin;
        return a.column_end<b.column_end;
    }
};
#include "pairwise_fragment_refinement.inc"
pairwise_detail::AnchorVec ReselectFragments(pairwise_detail::AnchorVec& anchors,
        const Sequences& refs,const Sequences& qrys,const PairwiseCoreOptions& o,
        PairwiseStatistics& stats,std::vector<RecoveryFragment>* trace=nullptr) {
    const auto begin=KernelClock::now();
    if(!o.min_cluster)throw AlignmentError("fragment selection needs positive minimum support");
    stats.post_selection_strategy="fragment-reselect-v1";
    stats.fragment_parents=anchors.size();
    OccupiedIntervals ro,qo;
    struct Queue:std::priority_queue<RecoveryFragment,std::vector<RecoveryFragment>,FragmentHigher> {
        std::uint64_t bytes() const {return this->c.capacity()*sizeof(RecoveryFragment);}
    } queue;
    std::vector<double> identities(anchors.size());
    const bool refine=o.fragment_local_exchange||o.fragment_short_exact;
    std::vector<Length> parent_matches(refine?anchors.size():0);
    std::vector<RecoveryFragment> selected;
    std::vector<RecoveryFragment> pending;
    pairwise_detail::AnchorVec output;
    const auto observe=[&](bool rss) {
        const auto bytes=queue.bytes()+identities.capacity()*sizeof(double)+
            pending.capacity()*sizeof(RecoveryFragment)+output.capacity()*sizeof(pairwise_detail::Anchor);
        stats.fragment_auxiliary_peak_bytes=std::max(stats.fragment_auxiliary_peak_bytes,bytes);
        if(rss&&o.resident_bytes_callback)stats.fragment_peak_rss_bytes=std::max(stats.fragment_peak_rss_bytes,o.resident_bytes_callback());
    };
    for(size_t parent=0;parent<anchors.size();++parent) {
        KernelCheck(o,"fragment-selection-validate");
        auto& a=anchors[parent];a.ref_selected=a.qry_selected=false;
        const auto& r=refs.at(a.reference_id)->bases;const auto& q=qrys.at(a.query_id)->bases;
        if(a.reference_begin>r.size()||a.ref_len>r.size()-a.reference_begin||
           a.query_begin>q.size()||a.qry_len>q.size()-a.query_begin)
            throw AlignmentError("fragment parent span overflow");
        pending.clear();
        RecoveryFragment f{};bool open=false;Length columns=0,matches=0;
        Length last_col=0,last_r=0,last_q=0,last_matches=0;
        const auto finish=[&] {
            if(open&&last_matches>=o.min_cluster) {
                f.parent=parent;f.column_end=last_col;f.re=last_r;
                if(a.strand==Strand::Forward)f.qe=last_q;
                else {f.qe=f.qb;f.qb=last_q;}
                f.matches=last_matches;f.columns=last_col-f.column_begin;
                pending.push_back(f);
            } else if(open)++stats.fragment_support_rejected;
            open=false;f={};last_matches=0;
        };
        const auto visit=[&](unsigned op,Length rp,Length qp,char x,char y) {
            if((columns&65535)==0)KernelCheck(o,"fragment-selection-validate");
            if(columns==UINT64_MAX)throw AlignmentError("fragment column count overflow");
            const bool paired=op!=1&&op!=2;
            const bool same=paired&&x==y&&x!='N';
            matches+=same;
            if(x=='N'||y=='N'){finish();++columns;++stats.fragment_n_breaks;return;}
            const auto cursor=a.strand==Strand::Forward?a.query_begin+qp:a.query_begin+a.qry_len-qp;
            if(paired&&!open){open=true;f.column_begin=columns;f.rb=rp;f.qb=cursor;}
            if(open&&same)++f.matches;
            ++columns;
            if(open&&paired){last_col=columns;last_r=rp+1;last_q=a.strand==Strand::Forward?cursor+1:cursor-1;last_matches=f.matches;}
        };
        const auto validated=ConvertAnchor(a,refs,qrys,false,visit);
        finish();
        if(!columns||columns!=a.alignment_length||validated.edit_distance>columns||matches!=columns-validated.edit_distance)
            throw AlignmentError("fragment parent summary mismatch");
        identities[parent]=static_cast<double>(matches)/static_cast<double>(columns);
        if(refine)parent_matches[parent]=matches;
        for(auto f0:pending) {
            f0.score=static_cast<double>(f0.matches)*identities[parent];
            queue.push(f0);++stats.fragment_proposed;
        }
        observe((parent&65535)==0);
    }
    std::vector<RecoveryFragment>().swap(pending);
    std::uint64_t popped=0;
    while(!queue.empty()) {
        KernelCheck(o,"fragment-selection");
        const auto f=queue.top();queue.pop();
        auto& a=anchors[f.parent];
        if(!GapIntervalFree(ro,a.reference_id,f.rb,f.re)||!GapIntervalFree(qo,a.query_id,f.qb,f.qe)) {
            ++stats.fragment_reclips;
            PairwiseStatistics scan;
            ScanRecoveryFragments(a,f.parent,f.column_begin,f.column_end,refs,qrys,ro,qo,o.min_cluster,o,[&](RecoveryFragment child) {
                child.score=static_cast<double>(child.matches)*identities[child.parent];
                if(FragmentHigher{}(f,child))throw AlignmentError("fragment clipping increased priority");
                queue.push(child);++stats.fragment_proposed;
            },nullptr,&scan);
            stats.fragment_support_rejected+=scan.post_selection_support_rejected;
            observe(false);continue;
        }
        if(trace&&!refine)trace->push_back(f);
        if(refine)selected.push_back(f);
        ro.add(a.reference_id,f.rb,f.re);qo.add(a.query_id,f.qb,f.qe);
        if(!refine&&f.column_begin==0&&f.column_end==a.alignment_length) {
            a.ref_selected=a.qry_selected=true;++stats.fragment_accepted_whole;
        } else if(!refine) {
            auto fragment=MaterializeRecovery(a,f);
            (void)ConvertAnchor(fragment,refs,qrys,false);
            output.push_back(std::move(fragment));++stats.fragment_accepted_clipped;
        }
        stats.fragment_matches+=f.matches;
        // Both-axis paired count follows from columns and consumed spans.
        const auto rspan=f.re-f.rb,qspan=f.qe-f.qb;
        if(rspan>UINT64_MAX-qspan||rspan+qspan<f.columns)throw AlignmentError("fragment paired count overflow");
        stats.fragment_paired_columns+=rspan+qspan-f.columns;
        stats.fragment_reference_bases+=rspan;stats.fragment_query_bases+=qspan;
        if((++popped&4095)==0)observe(true);
    }
    if(refine) {
        RefineFragmentSelection(anchors,refs,qrys,o,identities,parent_matches,selected,ro,qo,stats);
        stats.fragment_accepted_whole=stats.fragment_accepted_clipped=0;
        stats.fragment_matches=stats.fragment_paired_columns=stats.fragment_reference_bases=stats.fragment_query_bases=0;
        for(const auto& f:selected) {
            auto& a=anchors[f.parent];if(trace)trace->push_back(f);
            if(f.column_begin==0&&f.column_end==a.alignment_length){a.ref_selected=a.qry_selected=true;++stats.fragment_accepted_whole;}
            else{auto part=MaterializeRecovery(a,f);(void)ConvertAnchor(part,refs,qrys,false);output.push_back(std::move(part));++stats.fragment_accepted_clipped;}
            FragmentTotals total;total.add(f);stats.fragment_matches+=total.matches;stats.fragment_paired_columns+=total.paired;
            stats.fragment_reference_bases+=total.reference;stats.fragment_query_bases+=total.query;
        }
    }
    observe(true);stats.fragment_seconds=Elapsed(begin);
    return output;
}

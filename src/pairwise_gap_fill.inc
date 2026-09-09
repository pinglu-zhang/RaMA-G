// Guarded, non-recursive post-selection gap filling. Input anchors are immutable.
struct GapFlankInfo { Length matches{}, columns{}; bool begin_exact{}, end_exact{}; };
struct GuardedGap {
    size_t left{}, right{};
    SequenceId rid{}, qid{};
    Length rb{}, re{}, qb{}, qe{};
    Strand strand{};
    double minimum_identity{};
};
bool GapIntervalFree(const OccupiedIntervals& occupied,SequenceId id,Length b,Length e) {
    return b<e && occupied.run(id,b,e-b)==std::pair<bool,Length>{false,e-b};
}
GapFlankInfo InspectGapFlank(const pairwise_detail::Anchor& a,const Sequences& refs,
                           const Sequences& qrys,const PairwiseCoreOptions& o) {
    KernelCheck(o,"gap-fill-flanks");
    const auto record=ConvertAnchor(a,refs,qrys,false);
    GapFlankInfo info;
    for(auto u:a.cigar) {
        if((u>>4)>UINT64_MAX-info.columns)throw AlignmentError("gap-fill column overflow");
        info.columns+=u>>4;
    }
    if(record.edit_distance>info.columns)throw AlignmentError("gap-fill match accounting");
    info.matches=info.columns-record.edit_distance;
    const auto& r=refs.at(a.reference_id)->bases;
    const auto& q=qrys.at(a.query_id)->bases;
    const auto exact=[&](bool end) {
        const auto unit=end?a.cigar.back():a.cigar.front();
        if((unit&15)==1||(unit&15)==2||(unit>>4)<o.gap_fill_min_match)return false;
        for(Length k=0;k<o.gap_fill_min_match;++k) {
            const auto ri=end?a.reference_begin+a.ref_len-1-k:a.reference_begin+k;
            const auto qi=a.strand==Strand::Forward
                ?(end?a.query_begin+a.qry_len-1-k:a.query_begin+k)
                :(end?a.query_begin+k:a.query_begin+a.qry_len-1-k);
            const char x=r.at(ri),y=a.strand==Strand::Forward?q.at(qi):Complement(q.at(qi));
            if(x=='N'||x!=y)return false;
        }
        return true;
    };
    info.begin_exact=exact(false);info.end_exact=exact(true);return info;
}
pairwise_detail::AnchorVec FillGuardedGaps(const pairwise_detail::AnchorVec& anchors,
        size_t original_count,const Sequences& refs,const Sequences& qrys,
        const PairwiseCoreOptions& o,PairwiseStatistics& stats) {
    using namespace pairwise_detail;
    const auto begin=KernelClock::now();
    if(!o.recover_uncovered_fragments||o.gap_fill_min_match==0||o.threads==0||
       (o.gap_fill!=GapFillMode::ExactGap&&o.gap_fill!=GapFillMode::Ksw2Gap))
        throw AlignmentError("gap-fill requires residual recovery, positive support/threads and a valid mode");
    stats.gap_fill_strategy=o.gap_fill==GapFillMode::ExactGap?"exact-gap":"ksw2-gap";
    OccupiedIntervals ro,qo;
    std::vector<size_t> selected;
    for(size_t i=0;i<anchors.size();++i) {
        const auto& a=anchors[i];if(!(a.ref_selected&&a.qry_selected))continue;
        ro.add(a.reference_id,a.reference_begin,a.reference_begin+a.ref_len);
        qo.add(a.query_id,a.query_begin,a.query_begin+a.qry_len);selected.push_back(i);
    }
    std::stable_sort(selected.begin(),selected.end(),[&](size_t i,size_t j) {
        const auto& a=anchors[i];const auto& b=anchors[j];
        return std::tie(a.query_id,a.query_begin,a.qry_len,a.reference_id,a.reference_begin,i)<
               std::tie(b.query_id,b.query_begin,b.qry_len,b.reference_id,b.reference_begin,j);
    });
    std::vector<GapFlankInfo> flanks(selected.size());
    ParallelKernel(selected.size(),o,[&](size_t i) {
        if(selected[i]<original_count)flanks[i]=InspectGapFlank(anchors[selected[i]],refs,qrys,o);
    });
    std::vector<GuardedGap> gaps;
    for(size_t i=1;i<selected.size();++i) {
        KernelCheck(o,"gap-fill-candidates");
        const auto ai=selected[i-1],bi=selected[i];const auto& a=anchors[ai];const auto& b=anchors[bi];
        if(a.query_id!=b.query_id)continue;
        ++stats.gap_fill_adjacent_pairs;
        if(ai>=original_count||bi>=original_count||a.reference_id!=b.reference_id||a.strand!=b.strand) {
            ++stats.gap_fill_geometry_rejected;continue;
        }
        GuardedGap g;g.left=ai;g.right=bi;g.rid=a.reference_id;g.qid=a.query_id;g.strand=a.strand;
        g.qb=a.query_begin+a.qry_len;g.qe=b.query_begin;
        g.rb=a.strand==Strand::Forward?a.reference_begin+a.ref_len:b.reference_begin+b.ref_len;
        g.re=a.strand==Strand::Forward?b.reference_begin:a.reference_begin;
        if(g.rb>=g.re||g.qb>=g.qe||g.re-g.rb>10000||g.qe-g.qb>10000) {
            ++stats.gap_fill_geometry_rejected;continue;
        }
        if(!GapIntervalFree(ro,g.rid,g.rb,g.re)||!GapIntervalFree(qo,g.qid,g.qb,g.qe)) {
            ++stats.gap_fill_occupied_rejected;continue;
        }
        const auto& x=flanks[i-1];const auto& y=flanks[i];
        const double xi=x.columns?static_cast<double>(x.matches)/x.columns:0;
        const double yi=y.columns?static_cast<double>(y.matches)/y.columns:0;
        const bool inward=a.strand==Strand::Forward?x.end_exact&&y.begin_exact:x.begin_exact&&y.end_exact;
        if(x.matches<o.min_cluster||y.matches<o.min_cluster||xi<0.98||yi<0.98||!inward) {
            ++stats.gap_fill_flank_rejected;continue;
        }
        g.minimum_identity=std::max(0.98,std::min(xi,yi));
        const auto r=std::string_view(refs.at(g.rid)->bases).substr(g.rb,g.re-g.rb);
        const auto q=std::string_view(qrys.at(g.qid)->bases).substr(g.qb,g.qe-g.qb);
        if(r.find('N')!=r.npos||q.find('N')!=q.npos){++stats.gap_fill_n_rejected;continue;}
        gaps.push_back(g);
    }
    std::vector<GapFlankInfo>().swap(flanks);std::vector<size_t>().swap(selected);
    stats.gap_fill_candidates=gaps.size();
    std::stable_sort(gaps.begin(),gaps.end(),[](const auto& a,const auto& b) {
        return std::tie(a.rid,a.rb,a.re,a.qid,a.qb,a.qe,a.strand,a.left,a.right)<
               std::tie(b.rid,b.rb,b.re,b.qid,b.qb,b.qe,b.strand,b.left,b.right);
    });
    struct Proposal { Anchor anchor; bool accepted{},exact{},ksw{},quality{},nonexact{}; Length paired{}; double seconds{}; };
    std::vector<Proposal> results(gaps.size());
    ParallelKernel(gaps.size(),o,[&](size_t i) {
        KernelCheck(o,"gap-fill-alignment");const auto& g=gaps[i];auto& p=results[i];
        const auto r=std::string_view(refs.at(g.rid)->bases).substr(g.rb,g.re-g.rb);
        const auto q=std::string_view(qrys.at(g.qid)->bases).substr(g.qb,g.qe-g.qb);
        bool exact=r.size()==q.size();
        for(size_t k=0;exact&&k<r.size();++k)exact=r[k]==(g.strand==Strand::Forward?q[k]:Complement(q[q.size()-1-k]));
        Cigar_t cigar;
        if(exact){p.exact=true;cigar.push_back((static_cast<uint64_t>(r.size())<<4)|7);}
        else {
            if(o.gap_fill==GapFillMode::ExactGap){p.nonexact=true;return;}
            p.ksw=true;const auto start=KernelClock::now();
            auto aligned=globalAlignKSW2Result({r,false},{q,g.strand==Strand::Reverse});
            p.seconds=Elapsed(start);cigar=std::move(aligned.cigar);
        }
        KernelCheck(o,"gap-fill-alignment");
        const auto summary=summarizeCigar(cigar);
        Anchor full(g.rid,g.rb,r.size(),g.qid,g.qb,q.size(),g.strand,summary.alignment_length,0,std::move(cigar));
        (void)ConvertAnchor(full,refs,qrys,false); // Full global endpoint must close, even when rejected later.
        OccupiedIntervals empty_r,empty_q;std::vector<RecoveryFragment> fragments;
        ScanRecoveryFragments(full,0,0,UINT64_MAX,refs,qrys,empty_r,empty_q,1,o,[&](const auto& f){fragments.push_back(f);});
        if(fragments.empty()){p.quality=true;return;}
        if(fragments.size()!=1)throw AlignmentError("gap-fill unexpected internal break");
        const auto& f=fragments.front();
        if(!exact&&(f.matches<o.min_cluster||static_cast<double>(f.matches)/f.columns<g.minimum_identity)) {
            p.quality=true;return;
        }
        p.anchor=MaterializeRecovery(full,f);
        (void)ConvertAnchor(p.anchor,refs,qrys,false);
        for(auto u:p.anchor.cigar)if((u&15)!=1&&(u&15)!=2)p.paired+=u>>4;
        p.accepted=true;
    });
    AnchorVec output;
    for(size_t i=0;i<results.size();++i) {
        KernelCheck(o,"gap-fill-publication");auto& p=results[i];
        stats.gap_fill_exact+=p.exact;stats.gap_fill_ksw_calls+=p.ksw;stats.gap_fill_ksw_seconds+=p.seconds;
        stats.gap_fill_quality_rejected+=p.quality;stats.gap_fill_nonexact_rejected+=p.nonexact;
        if(!p.accepted)continue;
        auto& a=p.anchor;const auto re=a.reference_begin+a.ref_len,qe=a.query_begin+a.qry_len;
        if(!GapIntervalFree(ro,a.reference_id,a.reference_begin,re)||!GapIntervalFree(qo,a.query_id,a.query_begin,qe)) {
            ++stats.gap_fill_conflict_rejected;continue;
        }
        ro.add(a.reference_id,a.reference_begin,re);qo.add(a.query_id,a.query_begin,qe);
        ++stats.gap_fill_accepted;stats.gap_fill_reference_bases+=a.ref_len;stats.gap_fill_query_bases+=a.qry_len;
        stats.gap_fill_paired_columns+=p.paired;output.push_back(std::move(a));
    }
    stats.gap_fill_seconds=Elapsed(begin);return output;
}

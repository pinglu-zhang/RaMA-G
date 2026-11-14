// Two-flank recovery over the frozen linked main path. Rejected packed CIGARs
// stay borrowed until this function returns; only accepted paths own output.
// Every window is an open interval between mutual neighbours on BOTH axes.
struct VacantWindow {
    size_t left{}, right{};
    RecoveryBounds bounds;
    SequenceId rid{}, qid{};
    Strand strand{};
    double identity{};
    std::vector<RecoveryFragment> alternatives;
    bool exhausted{};
};
std::vector<CigarOp> VacantColumnSlice(const std::vector<CigarOp>& ops,Length begin,Length end) {
    std::vector<CigarOp> out;Length col=0;
    for(const auto& op:ops) {
        const auto b=std::max(col,begin),e=std::min(col+op.length,end);
        if(b<e)AppendLinkOp(out,op.operation,e-b);
        col+=op.length;if(col>=end)break;
    }
    if(col<end)throw AlignmentError("vacant slice outside CIGAR");
    return out;
}
VacantChainStatistics RecoverVacantChains(
        std::span<const SequenceRecord> reference_sequences,
        std::span<const SequenceRecord> query_sequences,
        const pairwise_detail::AnchorVec& parents,
        std::vector<AlignmentRecord>& records,LinkGapMode mode,
        const PairwiseCoreOptions& options,std::uint64_t budget) {
    const auto start=KernelClock::now();
    const auto refs=CheckSequences(reference_sequences),qrys=CheckSequences(query_sequences);
    VacantChainStatistics stats;stats.input_records=records.size();stats.work_budget=budget;
    const auto ro=LinkAxisOrder(records,true),qo=LinkAxisOrder(records,false);
    std::vector<size_t> qrank(records.size());for(size_t i=0;i<qo.size();++i)qrank[qo[i]]=i;
    std::vector<MainPathScan> scans;scans.reserve(records.size());
    for(const auto& a:records)scans.push_back(ScanLinkRecord(a,refs,qrys));
    std::vector<VacantWindow> windows;
    std::map<SequenceId,std::vector<size_t>> query_windows;
    for(size_t i=1;i<ro.size();++i) {
        const auto li=ro[i-1],ri=ro[i];const auto& a=records[li];const auto& b=records[ri];
        const bool fw=a.strand==Strand::Forward;
        if(a.reference_id!=b.reference_id||a.query_id!=b.query_id||a.strand!=b.strand||
           (fw?qrank[li]+1!=qrank[ri]:qrank[ri]+1!=qrank[li]))continue;
        VacantWindow w;w.left=li;w.right=ri;w.rid=a.reference_id;w.qid=a.query_id;w.strand=a.strand;
        w.bounds={a.reference_end,b.reference_begin,fw?a.query_end:b.query_end,fw?b.query_begin:a.query_begin};
        if(w.bounds.rb>=w.bounds.re||w.bounds.qb>=w.bounds.qe)continue;
        const auto& x=scans[li];const auto& y=scans[ri];
        const auto xi=static_cast<double>(x.matches)/x.columns,yi=static_cast<double>(y.matches)/y.columns;
        if(x.matches<options.min_cluster||y.matches<options.min_cluster||x.end_exact<options.gap_fill_min_match||
           y.begin_exact<options.gap_fill_min_match||xi<.98||yi<.98){++stats.flank_rejected;continue;}
        w.identity=std::max(.98,std::min(xi,yi));
        query_windows[w.qid].push_back(windows.size());windows.push_back(std::move(w));
    }
    stats.windows=windows.size();
    for(auto& [id,ids]:query_windows) {
        (void)id;std::sort(ids.begin(),ids.end(),[&](size_t a,size_t b){return windows[a].bounds.qb<windows[b].bounds.qb;});
    }
    // Bound metadata independently of the already-owned candidate CIGARs. A
    // truncated window is rejected, never declared unique from a partial list.
    constexpr size_t per_window_cap=4096;
    constexpr std::uint64_t descriptor_cap=512ULL<<20;
    using Key=std::tuple<SequenceId,Strand,Length,Length,Length,Length,pairwise_detail::Cigar_t>;
    std::vector<std::set<Key>> seen(windows.size());
    std::vector<std::uint64_t> bytes(windows.size());std::uint64_t live_bytes=0;
    const OccupiedIntervals empty;
    for(size_t pi=0;pi<parents.size();++pi) {
        if((pi&4095)==0)KernelCheck(options,"vacant-candidates");
        const auto& a=parents[pi];const auto pos=query_windows.find(a.query_id);
        if(pos==query_windows.end())continue;
        const auto& ids=pos->second;
        auto it=std::lower_bound(ids.begin(),ids.end(),a.query_begin,[&](size_t wi,Length p){return windows[wi].bounds.qe<=p;});
        for(;it!=ids.end()&&windows[*it].bounds.qb<a.query_begin+a.qry_len;++it) {
            const auto wi=*it;auto& w=windows[wi];if(w.exhausted)continue;
            ++stats.candidate_checks;
            // Include competing reference loci, strands and occupied partners.
            const RecoveryBounds bounds{0,refs.at(a.reference_id)->size(),w.bounds.qb,w.bounds.qe};
            ScanRecoveryFragments(a,pi,0,UINT64_MAX,refs,qrys,empty,empty,options.min_cluster,options,[&](const auto& f){
                if(w.exhausted)return;
                const auto packed=MaterializeRecovery(a,f).cigar;
                Key key{a.reference_id,a.strand,f.rb,f.re,f.qb,f.qe,packed};
                if(seen[wi].contains(key)){++stats.duplicates;return;}
                const auto need=sizeof(Key)+sizeof(RecoveryFragment)+128+packed.capacity()*sizeof(packed[0]);
                if(w.alternatives.size()>=per_window_cap||need>descriptor_cap-live_bytes) {
                    w.exhausted=true;++stats.budget_windows;live_bytes-=bytes[wi];bytes[wi]=0;
                    std::vector<RecoveryFragment>().swap(w.alternatives);std::set<Key>().swap(seen[wi]);return;
                }
                seen[wi].insert(std::move(key));w.alternatives.push_back(f);
                live_bytes+=need;bytes[wi]+=need;stats.auxiliary_peak_bytes=std::max(stats.auxiliary_peak_bytes,live_bytes);
                ++stats.fragments;
            },&bounds);
        }
    }
    // No newly accepted record becomes a guard or generates another window.
    struct Join { std::vector<CigarOp> cigar;Length paired{};bool accepted{}; };
    std::vector<Join> joins(records.size());
    for(size_t wi=0;wi<windows.size();++wi) {
        KernelCheck(options,"vacant-chain");
        const auto& w=windows[wi];std::string reason="no-candidate",chosen;
        auto trace=[&]{
            std::ostringstream row;
            row<<wi<<'\t'<<w.left<<'\t'<<w.right<<'\t'<<w.rid<<'\t'<<w.bounds.rb<<'\t'<<w.bounds.re
               <<'\t'<<w.qid<<'\t'<<w.bounds.qb<<'\t'<<w.bounds.qe<<'\t'<<(w.strand==Strand::Forward?'+':'-')
               <<'\t'<<w.alternatives.size()<<'\t'<<reason<<'\t'<<chosen;
            stats.decisions.push_back(row.str());
        };
        if(w.exhausted){reason="candidate-budget";trace();continue;}
        if(w.alternatives.size()*w.alternatives.size()>2000000){++stats.budget_windows;reason="competition-budget";trace();continue;}
        struct Node { RecoveryFragment f;AlignmentRecord record;MainPathScan scan;Length qb{},qe{}; };
        std::vector<Node> nodes;
        const auto qlen=qrys.at(w.qid)->size();
        for(const auto& source:w.alternatives) {
            const auto& a=parents[source.parent];
            if(a.reference_id!=w.rid||a.strand!=w.strand)continue;
            ScanRecoveryFragments(a,source.parent,source.column_begin,source.column_end,refs,qrys,empty,empty,
                options.min_cluster,options,[&](const auto& f){
                auto anchor=MaterializeRecovery(a,f);auto r=ConvertAnchor(anchor,refs,qrys);
                const auto x=ScanLinkRecord(r,refs,qrys);
                if(x.begin_exact<options.gap_fill_min_match||x.end_exact<options.gap_fill_min_match||
                   static_cast<double>(x.matches)/x.columns<w.identity||x.score<=0){++stats.quality_rejected;return;}
                // Compare the SAME complete query interval. Identical pairings
                // are duplicates, not rival loci. Ties at a different locus stay ambiguous.
                bool ambiguous=false;
                for(const auto& other:w.alternatives) {
                    if(other.qb>f.qb||other.qe<f.qe||other.parent==f.parent)continue;
                    const auto& p=parents[other.parent];
                    const RecoveryBounds range{0,refs.at(p.reference_id)->size(),f.qb,f.qe};
                    ScanRecoveryFragments(p,other.parent,other.column_begin,other.column_end,refs,qrys,empty,empty,
                        options.min_cluster,options,[&](const auto& alt){
                        if(alt.qb!=f.qb||alt.qe!=f.qe)return;
                        auto pa=MaterializeRecovery(p,alt);const auto ar=ConvertAnchor(pa,refs,qrys);
                        if(ar.reference_id==r.reference_id&&ar.strand==r.strand&&ar.reference_begin==r.reference_begin&&
                           ar.reference_end==r.reference_end&&ar.cigar==r.cigar)return;
                        const auto ax=ScanLinkRecord(ar,refs,qrys);
                        if(ax.score>=x.score)ambiguous=true;
                    },&range);
                    if(ambiguous)break;
                }
                if(ambiguous){++stats.ambiguous_fragments;return;}
                for(const auto& n:nodes)if(n.record.reference_begin==r.reference_begin&&n.record.reference_end==r.reference_end&&
                    n.record.query_begin==r.query_begin&&n.record.query_end==r.query_end&&n.record.cigar==r.cigar)return;
                const auto qb=w.strand==Strand::Forward?r.query_begin:qlen-r.query_end;
                const auto qe=w.strand==Strand::Forward?r.query_end:qlen-r.query_begin;
                nodes.push_back({f,std::move(r),x,qb,qe});
            },&w.bounds);
        }
        stats.eligible_fragments+=nodes.size();
        if(nodes.empty()){trace();continue;}
        std::sort(nodes.begin(),nodes.end(),[](const auto& a,const auto& b){
            return std::tie(a.record.reference_begin,a.qb,a.record.reference_end,a.qe,a.f.parent)<
                   std::tie(b.record.reference_begin,b.qb,b.record.reference_end,b.qe,b.f.parent);
        });
        if(nodes.size()*nodes.size()>2000000){++stats.budget_windows;reason="chain-budget";trace();continue;}
        const auto& left=records[w.left];const auto& right=records[w.right];
        const auto qend=[&](const AlignmentRecord& a){return w.strand==Strand::Forward?a.query_end:qlen-a.query_begin;};
        const auto qbegin=[&](const AlignmentRecord& a){return w.strand==Strand::Forward?a.query_begin:qlen-a.query_end;};
        const auto edge=[&](Length re,Length qe,Length rb,Length qb)->std::optional<std::int64_t>{
            if(re>rb||qe>qb)return {};
            const auto rn=rb-re,qn=qb-qe;
            if(!rn&&!qn)return 0;
            if(mode==LinkGapMode::None)return {};
            if((rn==qn&&rn<=10000)||(rn<=1000&&qn<=1000)) {
                const auto delta=rn>qn?rn-qn:qn-rn;
                // This is a proposal order, not an alignment score. Every
                // selected edge is subsequently aligned and rescored.
                return delta?40+3*static_cast<std::int64_t>(delta):0;
            }
            return {};
        };
        const auto neg=std::numeric_limits<std::int64_t>::min()/4;
        std::vector<std::int64_t> dp(nodes.size(),neg);std::vector<size_t> prev(nodes.size(),SIZE_MAX);
        std::vector<bool> tied(nodes.size());
        for(size_t j=0;j<nodes.size();++j) {
            const auto& b=nodes[j];
            if(auto c=edge(left.reference_end,qend(left),b.record.reference_begin,b.qb))dp[j]=b.scan.score-*c;
            for(size_t i=0;i<j;++i)if(dp[i]!=neg) {
                const auto& a=nodes[i];auto c=edge(a.record.reference_end,a.qe,b.record.reference_begin,b.qb);
                if(!c)continue;
                const auto score=dp[i]+b.scan.score-*c;
                if(score>dp[j]){dp[j]=score;prev[j]=i;tied[j]=tied[i];}
                else if(score==dp[j])tied[j]=true;
            }
        }
        std::int64_t best=neg;size_t last=SIZE_MAX;bool tie=false;
        for(size_t i=0;i<nodes.size();++i)if(dp[i]!=neg) {
            auto c=edge(nodes[i].record.reference_end,nodes[i].qe,right.reference_begin,qbegin(right));
            if(!c)continue;
            const auto score=dp[i]-*c;
            if(score>best){best=score;last=i;tie=tied[i];}else if(score==best)tie=true;
        }
        if(last==SIZE_MAX){++stats.no_chain;reason="no-complete-chain";trace();continue;}
        if(tie){++stats.tied_chains;reason="tied-chain";trace();continue;}
        std::vector<size_t> path;for(size_t i=last;i!=SIZE_MAX;i=prev[i])path.push_back(i);std::reverse(path.begin(),path.end());
        std::vector<AlignmentRecord> trial{left};
        for(const auto i:path){trial.push_back(nodes[i].record);chosen+=(chosen.empty()?"":",")+std::to_string(nodes[i].f.parent);}
        trial.push_back(right);
        auto local=options;local.threads=1;local.gap_fill_work_budget=budget-stats.work_reserved;
        const auto ls=LinkMainPath(reference_sequences,query_sequences,trial,mode,local,&refs,&qrys);
        stats.work_reserved+=ls.work_reserved;stats.ksw_calls+=ls.ksw_calls;
        if(trial.size()!=1){++stats.failed_joins;reason="unresolved-seam";trace();continue;}
        const auto full=ScanLinkRecord(trial[0],refs,qrys);
        if(full.score<=scans[w.left].score+scans[w.right].score){++stats.quality_rejected;reason="no-score-gain";trace();continue;}
        auto& join=joins[w.left];join.cigar=VacantColumnSlice(trial[0].cigar,scans[w.left].columns,full.columns-scans[w.right].columns);
        join.paired=full.paired-scans[w.left].paired-scans[w.right].paired;join.accepted=true;
        if(!join.paired)throw AlignmentError("vacant bridge contains no pairs");
        ++stats.accepted_chains;stats.accepted_fragments+=path.size();stats.added_paired_columns+=join.paired;
        stats.added_reference_bases+=w.bounds.re-w.bounds.rb;stats.added_query_bases+=w.bounds.qe-w.bounds.qb;
        reason="accepted";trace();
    }
    // Adjacent accepted windows can share a protected flank. Concatenate that
    // flank once, without making a recovered chain into a new search root.
    std::vector<AlignmentRecord> output;output.reserve(records.size());
    for(size_t i=0;i<ro.size();) {
        auto merged=records[ro[i]];bool changed=false;
        while(i+1<ro.size()&&joins[ro[i]].accepted) {
            const auto& join=joins[ro[i]];
            for(const auto& op:join.cigar)AppendLinkOp(merged.cigar,op.operation,op.length);
            const auto& next=records[ro[++i]];
            for(const auto& op:next.cigar)AppendLinkOp(merged.cigar,op.operation,op.length);
            merged.reference_end=next.reference_end;merged.query_begin=std::min(merged.query_begin,next.query_begin);
            merged.query_end=std::max(merged.query_end,next.query_end);merged.primary=merged.primary||next.primary;changed=true;
        }
        if(changed){const auto scan=ScanLinkRecord(merged,refs,qrys);merged.score=scan.score;merged.edit_distance=scan.edits;}
        output.push_back(std::move(merged));++i;
    }
    (void)LinkAxisOrder(output,true);(void)LinkAxisOrder(output,false);
    if(output.size()+stats.accepted_chains!=records.size())throw AlignmentError("vacant fragmentation contract");
    stats.output_records=output.size();stats.seconds=Elapsed(start);
    KernelCheck(options,"vacant-publication");
    records.swap(output);return stats;
}

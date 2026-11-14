// Conservative donor-chain recovery over a fixed B0 path. All original
// candidates participate in rivalry checks; no accepted chain seeds another.
// Resource failure aborts this optional transaction and retains B0 unchanged.
// Copy metadata without duplicating the caller-owned protected CIGAR.
AlignmentRecord CRMetadata(const AlignmentRecord& a) {
    AlignmentRecord b;
    b.reference_id=a.reference_id;b.reference_begin=a.reference_begin;b.reference_end=a.reference_end;
    b.query_id=a.query_id;b.query_begin=a.query_begin;b.query_end=a.query_end;b.strand=a.strand;
    b.score=a.score;b.edit_distance=a.edit_distance;b.primary=a.primary;return b;
}
// Count exactly the normalized operations without constructing a CIGAR.
// This also checks merged-run overflow before any ownership publication.
struct CRCigarPlan {
    size_t operations{};char last{};Length run{};
    void Append(const std::vector<CigarOp>& cigar,const PairwiseCoreOptions& options) {
        if(options.repair_budget)options.repair_budget->Consume(RepairWork::CandidateCigar,cigar.size());
        for(size_t i=0;i<cigar.size();++i) {
            if((i&4095)==0)KernelCheck(options,"chain-cigar-plan");
            const auto& op=cigar[i];if(!op.length)continue;
            if(last==op.operation) {
                if(op.length>UINT64_MAX-run)throw AlignmentError("chain planned CIGAR overflow");
                run+=op.length;
            } else {
                if(operations==SIZE_MAX)throw AlignmentError("chain planned CIGAR count overflow");
                ++operations;last=op.operation;run=op.length;
            }
        }
    }
};
VacantChainStatistics CRRecoverGuardedChains(
        std::span<const SequenceRecord> reference_sequences,
        std::span<const SequenceRecord> query_sequences,
        const pairwise_detail::AnchorVec& parents,
        std::vector<AlignmentRecord>& records,LinkGapMode mode,
        const PairwiseCoreOptions& options,std::uint64_t budget,
        RepairMemoryLease* output_memory=nullptr) {
    const auto start=KernelClock::now();
    const auto shared=options.repair_budget;
    auto check=[&](std::string_view phase){KernelCheck(options,phase);if(shared)shared->Check(phase);};
    auto charge=[&](RepairWork work,std::uint64_t count){if(shared)shared->Consume(work,count);};
    auto memory=shared?shared->Memory(0):RepairMemoryLease{};
    const auto base_bytes=records.size()*(sizeof(MainPathScan)+sizeof(VacantWindow)+4*sizeof(size_t)+128);
    // The caller already owns and accounts the full B0 transaction copy.
    // Reserve only this layer's descriptors and one reusable component workspace.
    constexpr std::uint64_t workspace_bytes=64ULL<<20;
    if(shared)memory.Grow(base_bytes+workspace_bytes);
    const auto refs=RepairSequenceLookup(reference_sequences,options),qrys=RepairSequenceLookup(query_sequences,options);
    VacantChainStatistics stats;stats.input_records=records.size();stats.work_budget=budget;
    const auto ro=LinkAxisOrder(records,true),qo=LinkAxisOrder(records,false);
    std::vector<size_t> qrank(records.size());for(size_t i=0;i<qo.size();++i)qrank[qo[i]]=i;
    std::vector<MainPathScan> scans;scans.reserve(records.size());
    for(const auto& a:records) {
        charge(RepairWork::CandidateCigar,a.cigar.size());
        MainPathScan s;s.score=a.score;s.edits=a.edit_distance;bool leading=true;
        for(const auto& op:a.cigar) {
            s.columns+=op.length;
            if(op.operation=='='||op.operation=='X')s.paired+=op.length;
            if(op.operation=='='){s.matches+=op.length;s.end_exact+=op.length;if(leading)s.begin_exact+=op.length;}
            else {s.end_exact=0;leading=false;}
        }
        scans.push_back(s);
    }
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
    if(windows.empty()){stats.output_records=records.size();stats.seconds=Elapsed(start);return stats;}
    // Select a stable bounded set of the largest actual open components BEFORE
    // expanding donor descriptors. All parents remain in each selected component's
    // rivalry universe, including parents at occupied reference loci.
    std::stable_sort(windows.begin(),windows.end(),[](const auto& a,const auto& b){
        const auto ap=std::min(a.bounds.re-a.bounds.rb,a.bounds.qe-a.bounds.qb);
        const auto bp=std::min(b.bounds.re-b.bounds.rb,b.bounds.qe-b.bounds.qb);
        if(ap!=bp)return ap>bp;
        return std::tie(a.qid,a.bounds.qb,a.rid,a.bounds.rb,a.strand)<
               std::tie(b.qid,b.bounds.qb,b.rid,b.bounds.rb,b.strand);
    });
    // Preserve the frozen conservative task allowance.
    const size_t task_cap=2048;
    if(windows.size()>task_cap){stats.budget_windows+=windows.size()-task_cap;windows.resize(task_cap);}
    query_windows.clear();for(size_t i=0;i<windows.size();++i)query_windows[windows[i].qid].push_back(i);
    for(auto& [id,ids]:query_windows) {
        (void)id;std::sort(ids.begin(),ids.end(),[&](size_t a,size_t b){return windows[a].bounds.qb<windows[b].bounds.qb;});
    }
    // Bound metadata independently of the already-owned candidate CIGARs. A
    // truncated window is rejected, never declared unique from a partial list.
    constexpr size_t per_window_cap=4096;
    const std::uint64_t cap=512ULL<<20;
    const std::uint64_t reserved_base=base_bytes+workspace_bytes;
    const auto ledger=shared?shared->Snapshot():RepairBudgetStatistics{};
    const std::uint64_t descriptor_cap=shared
        ? ledger.memory_limit_bytes-std::min(ledger.memory_limit_bytes,ledger.memory_current_bytes)
        : (reserved_base<cap?cap-reserved_base:0);
    auto descriptor_memory=shared?shared->Memory(0):RepairMemoryLease{};
    auto join_memory=shared?shared->Memory(0):RepairMemoryLease{};
    using Key=std::tuple<SequenceId,Strand,Length,Length,Length,Length,pairwise_detail::Cigar_t>;
    std::vector<std::set<Key>> seen(windows.size());
    std::vector<std::uint64_t> bytes(windows.size());std::uint64_t live_bytes=0;
    const OccupiedIntervals empty;
    for(size_t pi=0;pi<parents.size();++pi) {
        if((pi&4095)==0){check("chain-candidates");charge(RepairWork::CandidateCigar,std::min<size_t>(4096,parents.size()-pi));}
        // Charge the complete parent pass before filtering by query window.
        const auto& a=parents[pi];const auto pos=query_windows.find(a.query_id);
        if(pos==query_windows.end())continue;
        const auto& ids=pos->second;
        auto it=std::lower_bound(ids.begin(),ids.end(),a.query_begin,[&](size_t wi,Length p){return windows[wi].bounds.qe<=p;});
        for(;it!=ids.end()&&windows[*it].bounds.qb<a.query_begin+a.qry_len;++it) {
            const auto wi=*it;auto& w=windows[wi];if(w.exhausted)continue;
            ++stats.candidate_checks;
            charge(RepairWork::CandidateCigar,a.cigar.size());
            if(a.cigar.size()>(8ULL<<20)/sizeof(pairwise_detail::CigarUnit)) {
                w.exhausted=true;++stats.budget_windows;continue;
            }
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
                    std::vector<RecoveryFragment>().swap(w.alternatives);std::set<Key>().swap(seen[wi]);
                    if(shared){descriptor_memory.Reset();descriptor_memory=shared->Memory(live_bytes);}return;
                }
                if(shared)descriptor_memory.Grow(live_bytes+need);
                seen[wi].insert(std::move(key));w.alternatives.push_back(f);
                live_bytes+=need;bytes[wi]+=need;
                stats.auxiliary_peak_bytes=std::max(stats.auxiliary_peak_bytes,reserved_base+live_bytes);
                ++stats.fragments;
            },&bounds);
        }
    }
    // No newly accepted record becomes a guard or generates another window.
    struct Extension { std::vector<CigarOp> cigar;Length reference{},query{},paired{}; };
    struct Join { std::vector<CigarOp> cigar;Extension prefix,suffix;Length paired{};bool full{}; };
    // Only accepted edges own CIGARs. Partial edges do not merge their flanks.
    std::map<size_t,Join> joins;std::uint64_t join_bytes=0;
    const auto full_join=[&](size_t id){const auto p=joins.find(id);return p!=joins.end()&&p->second.full;};
    LocalRealignmentStatistics chain_work;chain_work.work_budget=budget;
    for(size_t wi=0;wi<windows.size();++wi) {
        check("chain-component");
        auto task=shared?shared->Task("chain-component"):RepairTaskGuard{};
        auto component_options=options;
        component_options.interruption_callback=[&](std::string_view phase){
            KernelCheck(options,phase);task.Check();
        };
        const auto& w=windows[wi];std::string reason="no-candidate",chosen;
        auto trace=[&]{
            std::ostringstream row;
            row<<wi<<'\t'<<w.left<<'\t'<<w.right<<'\t'<<w.rid<<'\t'<<w.bounds.rb<<'\t'<<w.bounds.re
               <<'\t'<<w.qid<<'\t'<<w.bounds.qb<<'\t'<<w.bounds.qe<<'\t'<<(w.strand==Strand::Forward?'+':'-')
               <<'\t'<<w.alternatives.size()<<'\t'<<reason<<'\t'<<chosen;
            auto text=row.str();if(!shared||shared->KeepDiagnostic(text.capacity()+sizeof(std::string)))stats.decisions.push_back(std::move(text));
        };
        if(w.exhausted){reason="candidate-budget";trace();continue;}
        if(w.alternatives.size()*w.alternatives.size()>2000000){++stats.budget_windows;reason="competition-budget";trace();continue;}
        std::uint64_t scratch_estimate=0;
        for(const auto& f:w.alternatives) {
            const auto extra=6ULL*sizeof(CigarOp)*f.columns+4ULL*parents[f.parent].cigar.size()*sizeof(pairwise_detail::CigarUnit)+512;
            if(extra>(64ULL<<20)-std::min<std::uint64_t>(scratch_estimate,64ULL<<20)){scratch_estimate=(64ULL<<20)+1;break;}
            scratch_estimate+=extra;
        }
        if(scratch_estimate>(64ULL<<20)){++stats.budget_windows;reason="component-memory-budget";trace();continue;}
        // scratch_estimate is already covered by the single reusable workspace.
        struct Node { RecoveryFragment f;AlignmentRecord record;MainPathScan scan;Length qb{},qe{}; };
        std::vector<Node> nodes;
        const auto qlen=qrys.at(w.qid)->size();
        for(const auto& source:w.alternatives) {
            const auto& a=parents[source.parent];
            if(a.reference_id!=w.rid||a.strand!=w.strand)continue;
            ScanRecoveryFragments(a,source.parent,source.column_begin,source.column_end,refs,qrys,empty,empty,
                options.min_cluster,component_options,[&](const auto& f){
                auto anchor=MaterializeRecovery(a,f);auto r=ConvertAnchor(anchor,refs,qrys,true,NoAnchorVisit{},&component_options);
                const auto x=ScanLinkRecord(r,refs,qrys,&component_options);
                if(x.begin_exact<options.gap_fill_min_match||x.end_exact<options.gap_fill_min_match||
                   static_cast<double>(x.matches)/x.columns<.98||x.score<=0){++stats.quality_rejected;return;}
                // Compare the SAME complete query interval. Identical pairings
                // are duplicates, not rival loci. Ties at a different locus stay ambiguous.
                bool ambiguous=false;
                for(const auto& other:w.alternatives) {
                    task.Check();charge(RepairWork::ChainEdges,1);
                    if(other.qe<=f.qb||other.qb>=f.qe||other.parent==f.parent)continue;
                    const auto& p=parents[other.parent];
                    const RecoveryBounds range{0,refs.at(p.reference_id)->size(),std::max(f.qb,other.qb),std::min(f.qe,other.qe)};
                    ScanRecoveryFragments(p,other.parent,other.column_begin,other.column_end,refs,qrys,empty,empty,
                        options.min_cluster,component_options,[&](const auto& alt){
                        // Partial rivals participate too: otherwise two short rival
                        // pieces could falsely make one long donor look unique.
                        const RecoveryBounds target_range{0,refs.at(a.reference_id)->size(),alt.qb,alt.qe};
                        ScanRecoveryFragments(a,f.parent,f.column_begin,f.column_end,refs,qrys,empty,empty,
                            options.min_cluster,component_options,[&](const auto& target){
                            if(target.qb!=alt.qb||target.qe!=alt.qe)return;
                            auto pa=MaterializeRecovery(p,alt),ta=MaterializeRecovery(a,target);
                            const auto ar=ConvertAnchor(pa,refs,qrys,true,NoAnchorVisit{},&component_options);
                            const auto tr=ConvertAnchor(ta,refs,qrys,true,NoAnchorVisit{},&component_options);
                            if(ar.reference_id==tr.reference_id&&ar.strand==tr.strand&&ar.reference_begin==tr.reference_begin&&
                               ar.reference_end==tr.reference_end&&ar.cigar==tr.cigar)return;
                            const auto ax=ScanLinkRecord(ar,refs,qrys,&component_options);
                            const auto tx=ScanLinkRecord(tr,refs,qrys,&component_options);
                            if(ax.score>=tx.score)ambiguous=true;
                        },&target_range);
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
        charge(RepairWork::ChainEdges,nodes.size()*nodes.size());
        for(size_t j=0;j<nodes.size();++j) {
            task.Check();
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
        // Assemble directly from borrowed flanks/nodes. No vector of copied
        // trial records and no copy of the protected left CIGAR is needed.
        CRCigarPlan left_plan,right_plan;
        left_plan.Append(left.cigar,component_options);right_plan.Append(right.cigar,component_options);
        std::uint64_t merged_ops=left_plan.operations+right_plan.operations;
        for(const auto i:path) {
            CRCigarPlan node_plan;node_plan.Append(nodes[i].record.cigar,component_options);
            merged_ops+=node_plan.operations;
        }
        const AlignmentRecord* previous=&left;
        for(size_t ti=0;ti<=path.size();++ti) {
            const auto& next=ti<path.size()?nodes[path[ti]].record:right;
            const auto rn=next.reference_begin-previous->reference_end;
            const auto qn=qbegin(next)-qend(*previous);
            if(rn>std::numeric_limits<std::uint64_t>::max()-merged_ops||
               qn>std::numeric_limits<std::uint64_t>::max()-merged_ops-rn)
                throw AlignmentError("chain CIGAR capacity overflow");
            merged_ops+=rn+qn;previous=&next;
        }
        if(merged_ops>std::numeric_limits<size_t>::max()/sizeof(CigarOp))
            throw AlignmentError("chain CIGAR allocation overflow");
        const auto merged_bytes=merged_ops*sizeof(CigarOp);
        // Node/rival scratch uses the reserved workspace. Only the portion of
        // this assembly not covered by its remaining space is additional.
        const auto scratch_free=workspace_bytes-scratch_estimate;
        auto assembly_memory=shared?shared->Memory(merged_bytes>scratch_free?merged_bytes-scratch_free:0):RepairMemoryLease{};
        AlignmentRecord merged=CRMetadata(left);merged.cigar.reserve(static_cast<size_t>(merged_ops));
        for(const auto& op:left.cigar)AppendLinkOp(merged.cigar,op.operation,op.length);
        bool complete=true;
        auto task_check=[&]{check("chain-gap");task.Check();};
        for(size_t ti=0;ti<=path.size();++ti) {
            const auto& next=ti<path.size()?nodes[path[ti]].record:right;
            if(ti<path.size()&&chosen.size()<480)chosen+=(chosen.empty()?"":",")+std::to_string(nodes[path[ti]].f.parent);
            const bool fw=merged.strand==Strand::Forward;
            const auto rb=merged.reference_end,re=next.reference_begin;
            const auto qb=fw?merged.query_end:next.query_end;
            const auto qe=fw?next.query_begin:merged.query_begin;
            if(rb>re||qb>qe){complete=false;break;}
            const auto rr=std::string_view(refs.at(merged.reference_id)->bases).substr(rb,re-rb);
            auto qq=qrys.at(merged.query_id)->bases.substr(qb,qe-qb);
            if(!fw){std::reverse(qq.begin(),qq.end());for(auto& c:qq)c=Complement(c);}
            auto gap=AlignRepairGap(rr,qq,task_check,chain_work,shared.get());
            if(!gap){complete=false;break;}
            for(const auto& op:*gap)AppendLinkOp(merged.cigar,op.operation,op.length);
            for(const auto& op:next.cigar)AppendLinkOp(merged.cigar,op.operation,op.length);
            merged.reference_end=next.reference_end;
            merged.query_begin=std::min(merged.query_begin,next.query_begin);
            merged.query_end=std::max(merged.query_end,next.query_end);
            merged.primary=merged.primary||next.primary;
        }
        stats.work_reserved=chain_work.work_reserved;stats.ksw_calls=chain_work.ksw_calls;
        if(!complete){++stats.failed_joins;reason="unresolved-seam";trace();continue;}
        const auto full=ScanLinkRecord(merged,refs,qrys,&component_options);
        // Baseline records have already been replayed by B0. Their protected
        // CIGARs are copied byte-for-byte; only the new joined record is replayed.
        if(full.score<=scans[w.left].score+scans[w.right].score){++stats.quality_rejected;reason="no-score-gain";trace();continue;}
        const auto slice_begin=scans[w.left].columns,slice_end=full.columns-scans[w.right].columns;
        size_t slice_ops=0;Length col=0;
        for(const auto& op:merged.cigar) {
            if(std::max(col,slice_begin)<std::min(col+op.length,slice_end))++slice_ops;
            col+=op.length;if(col>=slice_end)break;
        }
        if(col<slice_end)throw AlignmentError("chain slice outside CIGAR");
        // Check the entire new core, including inherited donor CIGARs. A
        // movable affine gap changes physical pairings even when EQX and the
        // complete score are identical; it must not become confirmed homology.
        auto core_memory=shared?shared->Memory(3ULL*slice_ops*sizeof(CigarOp)+sizeof(Join)+128):RepairMemoryLease{};
        std::vector<CigarOp> core;core.reserve(slice_ops);col=0;
        for(const auto& op:merged.cigar) {
            const auto b=std::max(col,slice_begin),e=std::min(col+op.length,slice_end);
            if(b<e)AppendLinkOp(core,op.operation,e-b);
            col+=op.length;if(col>=slice_end)break;
        }
        const auto rcore=std::string_view(refs.at(w.rid)->bases).substr(w.bounds.rb,w.bounds.re-w.bounds.rb);
        auto query_memory=shared?shared->Memory(w.strand==Strand::Reverse?w.bounds.qe-w.bounds.qb:0):RepairMemoryLease{};
        std::string reverse_core;
        std::string_view qcore=std::string_view(qrys.at(w.qid)->bases).substr(w.bounds.qb,w.bounds.qe-w.bounds.qb);
        if(w.strand==Strand::Reverse) {
            task_check();charge(RepairWork::LocalBases,qcore.size());task_check();
            reverse_core.assign(qcore);std::reverse(reverse_core.begin(),reverse_core.end());
            for(size_t k=0;k<reverse_core.size();++k) {
                if((k&65535)==0)task_check();
                reverse_core[k]=Complement(reverse_core[k]);
            }
            qcore=reverse_core;
        }
        const auto substitution=[](char x,char y) {
            const auto xi=pairwise_detail::ScoreChar2Idx[static_cast<unsigned char>(x)];
            const auto yi=pairwise_detail::ScoreChar2Idx[static_cast<unsigned char>(y)];
            return static_cast<std::int64_t>(pairwise_detail::dna5_simd_mat[xi*5+yi]);
        };
        auto ambiguity_memory=shared?shared->Memory(detail::GapAmbiguityWorkspaceBytes(core.size())):RepairMemoryLease{};
        charge(RepairWork::CandidateCigar,core.size());
        const auto ambiguity_check=[&](std::uint64_t inspected) {
            task_check();charge(RepairWork::LocalBases,inspected);task_check();
        };
        const auto ambiguity=detail::FindGapTranslationAmbiguity(rcore,qcore,core,substitution,ambiguity_check);
        Join join;
        if(!ambiguity.ambiguous) {
            join.cigar=std::move(core);join.full=true;
            join.paired=full.paired-scans[w.left].paired-scans[w.right].paired;
            if(!join.paired)throw AlignmentError("vacant bridge contains no pairs");
        } else {
            const auto core_columns=slice_end-slice_begin;
            if(ambiguity.begin_column>=ambiguity.end_column||ambiguity.end_column>core_columns)
                throw AlignmentError("invalid chain ambiguity bounds");
            stats.ambiguous_columns+=ambiguity.end_column-ambiguity.begin_column;
            ++stats.ambiguous_fragments;
            const auto make_extension=[&](Length begin,Length end,bool prefix) {
                Extension ext;ext.cigar.reserve(core.size());Length at=0;
                for(const auto& op:core) {
                    const auto b=std::max(at,begin),e=std::min(at+op.length,end);
                    if(b<e)AppendLinkOp(ext.cigar,op.operation,e-b);
                    at+=op.length;if(at>=end)break;
                }
                const auto gap=[](char op){return op=='I'||op=='D';};
                if(prefix)while(!ext.cigar.empty()&&gap(ext.cigar.back().operation))ext.cigar.pop_back();
                else {
                    size_t skip=0;while(skip<ext.cigar.size()&&gap(ext.cigar[skip].operation))++skip;
                    if(skip)ext.cigar.erase(ext.cigar.begin(),ext.cigar.begin()+skip);
                }
                for(const auto& op:ext.cigar) {
                    if(op.operation!='I')ext.reference+=op.length;
                    if(op.operation!='D')ext.query+=op.length;
                    if(op.operation=='='||op.operation=='X')ext.paired+=op.length;
                }
                if(!ext.paired){ext={};return ext;}
                AlignmentRecord piece=CRMetadata(left);piece.cigar=std::move(ext.cigar);
                piece.reference_begin=prefix?w.bounds.rb:w.bounds.re-ext.reference;
                piece.reference_end=prefix?w.bounds.rb+ext.reference:w.bounds.re;
                if(w.strand==Strand::Forward) {
                    piece.query_begin=prefix?w.bounds.qb:w.bounds.qe-ext.query;
                    piece.query_end=prefix?w.bounds.qb+ext.query:w.bounds.qe;
                } else {
                    piece.query_begin=prefix?w.bounds.qe-ext.query:w.bounds.qb;
                    piece.query_end=prefix?w.bounds.qe:w.bounds.qb+ext.query;
                }
                const auto part=ScanLinkRecord(piece,refs,qrys,&component_options);
                const auto previous=prefix?left.cigar.back().operation:piece.cigar.back().operation;
                const auto following=prefix?piece.cigar.front().operation:right.cigar.front().operation;
                // Separately scored equal-direction seam gaps paid two opens;
                // concatenation refunds precisely one frozen 40-point open.
                const auto seam_refund=previous==following&&gap(previous)?40:0;
                if(part.score+seam_refund<=0){ext={};return ext;}
                ext.cigar=std::move(piece.cigar);return ext;
            };
            join.prefix=make_extension(0,ambiguity.begin_column,true);
            join.suffix=make_extension(ambiguity.end_column,core_columns,false);
            join.paired=join.prefix.paired+join.suffix.paired;
            if(!join.paired){++stats.quality_rejected;reason="ambiguous-no-positive-outer-core";trace();continue;}
            if(join.prefix.reference+join.suffix.reference>w.bounds.re-w.bounds.rb||
               join.prefix.query+join.suffix.query>w.bounds.qe-w.bounds.qb)
                throw AlignmentError("overlapping stable chain extensions");
            ++stats.partial_chains;
            chosen+=";uncertain-columns="+std::to_string(ambiguity.begin_column)+":"+std::to_string(ambiguity.end_column);
        }
        const auto persistent=sizeof(Join)+128+
            (join.cigar.capacity()+join.prefix.cigar.capacity()+join.suffix.cigar.capacity())*sizeof(CigarOp);
        if(shared)join_memory.Grow(join_bytes+persistent);
        join_bytes+=persistent;
        const auto paired=join.paired;
        const auto new_reference=join.full?w.bounds.re-w.bounds.rb:join.prefix.reference+join.suffix.reference;
        const auto new_query=join.full?w.bounds.qe-w.bounds.qb:join.prefix.query+join.suffix.query;
        const bool complete_join=join.full;
        if(!joins.emplace(w.left,std::move(join)).second)throw AlignmentError("duplicate chain edge");
        if(complete_join)++stats.accepted_chains;
        stats.accepted_fragments+=path.size();stats.added_paired_columns+=paired;
        stats.added_reference_bases+=new_reference;stats.added_query_bases+=new_query;
        reason=complete_join?"accepted":"accepted-partial-ambiguous";trace();
    }
    // With no accepted bridge, retain every record and its owned CIGAR without
    // allocating an output CIGAR pool or replaying unchanged alignments.
    if(!stats.accepted_chains&&!stats.partial_chains){stats.output_records=records.size();stats.seconds=Elapsed(start);return stats;}
    // Donor descriptors are no longer needed during publication planning.
    std::vector<VacantWindow>().swap(windows);std::vector<std::set<Key>>().swap(seen);
    if(shared)descriptor_memory.Reset();
    // Precompute disjoint groups and final operation counts. Unchanged output
    // slots initially contain metadata only, so geometry can be checked before
    // any protected CIGAR ownership leaves the input.
    const auto output_count=records.size()-stats.accepted_chains;
    auto publication_memory=shared?shared->Memory(output_count*(sizeof(AlignmentRecord)+sizeof(size_t))):RepairMemoryLease{};
    std::vector<AlignmentRecord> output(output_count);
    std::vector<size_t> source(output_count,SIZE_MAX);
    std::uint64_t changed_bytes=0;
    size_t oi=0;
    for(size_t i=0;i<ro.size();) {
        check("chain-output-plan");
        const auto first=i;size_t last=i;
        while(last+1<ro.size()&&full_join(ro[last]))++last;
        const Extension* incoming=nullptr;const Extension* outgoing=nullptr;
        if(first) {
            const auto previous=joins.find(ro[first-1]);
            if(previous!=joins.end()&&!previous->second.full&&!previous->second.suffix.cigar.empty())
                incoming=&previous->second.suffix;
        }
        const auto next=joins.find(ro[last]);
        if(next!=joins.end()&&!next->second.full&&!next->second.prefix.cigar.empty())
            outgoing=&next->second.prefix;
        auto& merged=output[oi];merged=CRMetadata(records[ro[first]]);
        if(first==last&&!incoming&&!outgoing)source[oi]=ro[first];
        else {
            CRCigarPlan plan;
            if(incoming)plan.Append(incoming->cigar,options);
            plan.Append(records[ro[first]].cigar,options);
            for(size_t k=first;k<last;++k) {
                plan.Append(joins.at(ro[k]).cigar,options);
                plan.Append(records[ro[k+1]].cigar,options);
            }
            if(outgoing)plan.Append(outgoing->cigar,options);
            const auto ops=plan.operations;
            if(ops>std::numeric_limits<size_t>::max()/sizeof(CigarOp))
                throw AlignmentError("chain output CIGAR capacity overflow");
            changed_bytes+=ops*sizeof(CigarOp);
            if(shared)publication_memory.Grow(output_count*(sizeof(AlignmentRecord)+sizeof(size_t))+changed_bytes);
            merged.cigar.reserve(static_cast<size_t>(ops));
            if(incoming) {
                for(const auto& op:incoming->cigar)AppendLinkOp(merged.cigar,op.operation,op.length);
                merged.reference_begin-=incoming->reference;
                if(merged.strand==Strand::Forward)merged.query_begin-=incoming->query;
                else merged.query_end+=incoming->query;
            }
            for(const auto& op:records[ro[first]].cigar)AppendLinkOp(merged.cigar,op.operation,op.length);
            for(size_t k=first;k<last;++k) {
                for(const auto& op:joins.at(ro[k]).cigar)AppendLinkOp(merged.cigar,op.operation,op.length);
                const auto& following=records[ro[k+1]];
                for(const auto& op:following.cigar)AppendLinkOp(merged.cigar,op.operation,op.length);
                merged.reference_end=following.reference_end;merged.query_begin=std::min(merged.query_begin,following.query_begin);
                merged.query_end=std::max(merged.query_end,following.query_end);merged.primary=merged.primary||following.primary;
            }
            if(outgoing) {
                for(const auto& op:outgoing->cigar)AppendLinkOp(merged.cigar,op.operation,op.length);
                merged.reference_end+=outgoing->reference;
                if(merged.strand==Strand::Forward)merged.query_end+=outgoing->query;
                else merged.query_begin-=outgoing->query;
            }
            if(merged.cigar.size()!=ops)throw AlignmentError("chain output CIGAR count");
            const auto scan=ScanLinkRecord(merged,refs,qrys,&options);
            merged.score=scan.score;merged.edit_distance=scan.edits;
        }
        ++oi;i=last+1;
    }
    if(oi!=output_count)throw AlignmentError("chain output group count");
    (void)LinkAxisOrder(output,true);(void)LinkAxisOrder(output,false);
    stats.output_records=output.size();stats.seconds=Elapsed(start);
    check("chain-publication");
    // Every fallible allocation/replay/geometry/cancellation check is complete.
    // The following ownership transfer cannot allocate or throw.
    static_assert(std::is_nothrow_move_assignable_v<AlignmentRecord>);
    static_assert(std::is_nothrow_move_assignable_v<RepairMemoryLease>);
    static_assert(noexcept(records.swap(output)));
    // Keep this output reservation alive through the caller's final
    // ownership publication, without a second allocation.
    if(output_memory)*output_memory=std::move(publication_memory);
    for(size_t i=0;i<output.size();++i)if(source[i]!=SIZE_MAX)output[i]=std::move(records[source[i]]);
    records.swap(output);return stats;
}


LocalRealignmentStatistics ChainRecoverPaths(
        std::span<const SequenceRecord> references,std::span<const SequenceRecord> queries,
        const pairwise_detail::AnchorVec& parents,std::vector<AlignmentRecord>& fixed,
        const PairwiseCoreOptions& options,std::uint64_t budget,
        RepairMemoryLease* output_memory=nullptr) {
    const auto started=KernelClock::now();
    auto local=options;
    if(!local.repair_budget)local.repair_budget=std::make_shared<RepairBudget>(300.0);
    LocalRealignmentStatistics out;out.input_records=fixed.size();out.output_records=fixed.size();out.work_budget=budget;
    try {
        auto r=CRRecoverGuardedChains(references,queries,parents,fixed,LinkGapMode::Bounded,local,budget,output_memory);
        out.output_records=fixed.size();out.contexts=r.windows;out.candidate_checks=r.candidate_checks;
        out.fragments=r.fragments;out.candidate_budget=r.budget_windows;out.proposals=r.eligible_fragments;
        out.accepted=r.accepted_chains+r.partial_chains;out.merged=r.accepted_chains;
        out.partial_chains=r.partial_chains;out.ambiguous_columns=r.ambiguous_columns;
        out.ambiguous=r.ambiguous_fragments+r.tied_chains;
        out.unsupported=r.no_chain+r.quality_rejected;out.unattached=r.failed_joins;
        out.added_paired_columns=r.added_paired_columns;out.work_reserved=r.work_reserved;
        out.ksw_calls=r.ksw_calls;out.auxiliary_peak_bytes=r.auxiliary_peak_bytes;
        out.decisions=std::move(r.decisions);
    } catch(const RepairBudgetExceeded&) {
        // No publication happened: B0 records retain all their owned CIGARs.
        ++out.candidate_budget;
    }
    out.seconds=Elapsed(started);out.repair_budget=local.repair_budget->Snapshot();
    return out;
}

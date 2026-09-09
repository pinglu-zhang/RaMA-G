// Sequence-aware residual recovery. Original selected anchors remain immutable.
// Descriptors reference packed CIGAR column ranges; only accepted fragments own
// CIGAR storage. Ranking is a deterministic greedy policy, not a global optimum.
class OccupiedIntervals {
    std::map<SequenceId,std::map<Length,Length>> data;
public:
    void add(SequenceId id,Length b,Length e) {
        if(b>=e)throw AlignmentError("recovery empty occupied interval");
        auto& m=data[id];auto it=m.lower_bound(b);
        if(it!=m.begin()){auto p=std::prev(it);if(p->second>=b){b=p->first;e=std::max(e,p->second);it=m.erase(p);}}
        while(it!=m.end()&&it->first<=e){e=std::max(e,it->second);it=m.erase(it);}m.emplace_hint(it,b,e);
    }
    // Boundaries are expressed in consumed-base coordinates, including reverse.
    std::pair<bool,Length> run(SequenceId id,Length p,Length limit,bool reverse=false) const {
        const auto d=data.find(id);if(d==data.end())return {false,limit};
        const auto& m=d->second;
        if(!reverse){
            auto it=m.upper_bound(p);
            if(it!=m.begin()){const auto prev=std::prev(it);if(prev->second>p)return {true,std::min(limit,prev->second-p)};}
            return {false,it==m.end()?limit:std::min(limit,it->first-p)};
        }
        if(p==0)throw AlignmentError("recovery reverse coordinate underflow");
        auto it=m.lower_bound(p);
        if(it==m.begin())return {false,std::min(limit,p)};
        --it;
        if(it->second>=p)return {true,std::min(limit,p-it->first)};
        return {false,std::min(limit,p-it->second)};
    }
};
struct RecoveryFragment {
    size_t parent{};
    Length column_begin{},column_end{},rb{},re{},qb{},qe{},matches{},columns{};
    double score{};
};
struct RecoveryHigher {
    bool operator()(const RecoveryFragment& a,const RecoveryFragment& b) const {
        if(a.score!=b.score)return a.score<b.score;
        if(a.matches!=b.matches)return a.matches<b.matches;
        return std::tie(a.parent,a.rb,a.qb,a.re,a.qe,a.column_begin)>
               std::tie(b.parent,b.rb,b.qb,b.re,b.qe,b.column_begin);
    }
};
template<class Emit>
void ScanRecoveryFragments(const pairwise_detail::Anchor& a,size_t parent,
                          Length column_begin,Length column_end,
                          const Sequences& refs,const Sequences& qrys,
                          const OccupiedIntervals& ro,const OccupiedIntervals& qo,
                          Length minimum,const PairwiseCoreOptions& options,Emit emit) {
    const auto& r=refs.at(a.reference_id)->bases;const auto& q=qrys.at(a.query_id)->bases;
    Length rp=a.reference_begin,qp=a.strand==Strand::Forward?a.query_begin:a.query_begin+a.qry_len,col=0;
    RecoveryFragment f{};bool open=false;Length last_col=0,last_r=0,last_q=0,last_matches=0;
    const auto finish=[&]{
        if(open&&last_matches>=minimum){
            f.parent=parent;f.column_end=last_col;f.re=last_r;
            if(a.strand==Strand::Forward)f.qe=last_q;else{f.qe=f.qb;f.qb=last_q;}
            f.matches=last_matches;f.columns=f.column_end-f.column_begin;
            const double identity=static_cast<float>(f.matches)/f.columns;
            f.score=f.columns*std::pow(identity,2);emit(f);
        }
        open=false;f={};last_matches=0;
    };
    const auto advance=[&](unsigned op,Length n){if(op!=1)rp+=n;if(op!=2){if(a.strand==Strand::Forward)qp+=n;else qp-=n;}col+=n;};
    for(const auto unit:a.cigar){
        const unsigned op=unit&15;Length left=unit>>4;
        if(left==0||(op!=0&&op!=1&&op!=2&&op!=7&&op!=8))throw AlignmentError("recovery invalid packed CIGAR");
        while(left){
            if(col>=column_end){finish();return;}
            if(col<column_begin){const auto n=std::min(left,column_begin-col);advance(op,n);left-=n;continue;}
            if(options.interruption_callback)options.interruption_callback("coverage-recovery");
            Length n=std::min(left,column_end-col);
            const auto rr=op==1?std::pair{false,n}:ro.run(a.reference_id,rp,n);
            const auto qr=op==2?std::pair{false,n}:qo.run(a.query_id,qp,n,a.strand==Strand::Reverse);
            n=std::min(rr.second,qr.second);if(n==0)throw AlignmentError("recovery failed to advance");
            if(rr.first||qr.first){finish();advance(op,n);left-=n;continue;}
            // Occupied stretches were skipped arithmetically. Inspect only free
            // bases, retaining exact N breaks and match support for each slice.
            for(Length k=0;k<n;++k){
                if((k&65535)==0&&options.interruption_callback)options.interruption_callback("coverage-recovery");
                const char x=op==1?'A':r.at(rp);
                const char y=op==2?'A':(a.strand==Strand::Forward?q.at(qp):Complement(q.at(qp-1)));
                const bool paired=op!=1&&op!=2;
                if(x=='N'||y=='N'){finish();advance(op,1);continue;}
                if(paired&&!open){open=true;f.column_begin=col;f.rb=rp;f.qb=qp;f.matches=0;}
                if(open&&paired&&x==y)++f.matches;
                advance(op,1);
                if(open&&paired){last_col=col;last_r=rp;last_q=qp;last_matches=f.matches;}
            }
            left-=n;
        }
    }
    finish();
}
pairwise_detail::Anchor MaterializeRecovery(const pairwise_detail::Anchor& a,const RecoveryFragment& f) {
    pairwise_detail::Cigar_t cigar;Length col=0;
    for(auto unit:a.cigar){
        const Length n=unit>>4;const auto b=std::max(col,f.column_begin),e=std::min(col+n,f.column_end);
        if(b<e){const auto op=unit&15;const auto length=e-b;
            if(length>(UINT64_MAX>>4))throw AlignmentError("recovery packed length overflow");
            if(!cigar.empty()&&(cigar.back()&15)==op){const auto total=(cigar.back()>>4)+length;if(total>(UINT64_MAX>>4))throw AlignmentError("recovery packed sum overflow");cigar.back()=(total<<4)|op;}
            else cigar.push_back((length<<4)|op);
        }
        col+=n;if(col>=f.column_end)break;
    }
    pairwise_detail::Anchor result(a.reference_id,f.rb,f.re-f.rb,a.query_id,f.qb,f.qe-f.qb,a.strand,f.columns,f.matches,std::move(cigar));
    result.ref_selected=result.qry_selected=true;return result;
}
pairwise_detail::AnchorVec RecoverResidualAnchors(const pairwise_detail::AnchorVec& anchors,
        const Sequences& refs,const Sequences& qrys,const PairwiseCoreOptions& options,
        PairwiseStatistics& stats) {
    const auto begin=KernelClock::now();OccupiedIntervals ro,qo;
    for(const auto& a:anchors)if(a.ref_selected&&a.qry_selected){ro.add(a.reference_id,a.reference_begin,a.reference_begin+a.ref_len);qo.add(a.query_id,a.query_begin,a.query_begin+a.qry_len);}
    std::priority_queue<RecoveryFragment,std::vector<RecoveryFragment>,RecoveryHigher> queue;
    for(size_t i=0;i<anchors.size();++i){const auto& a=anchors[i];if(a.ref_selected&&a.qry_selected)continue;
        ++stats.recovery_candidates_checked;
        (void)ConvertAnchor(a,refs,qrys,false); // Validate even candidates with no residual.
        ScanRecoveryFragments(a,i,0,UINT64_MAX,refs,qrys,ro,qo,options.min_cluster,options,[&](const auto& f){queue.push(f);++stats.recovery_proposed_fragments;});
    }
    pairwise_detail::AnchorVec output;
    while(!queue.empty()){
        KernelCheck(options,"coverage-recovery");const auto f=queue.top();queue.pop();const auto& a=anchors[f.parent];
        bool unchanged=false;
        ScanRecoveryFragments(a,f.parent,f.column_begin,f.column_end,refs,qrys,ro,qo,options.min_cluster,options,[&](const auto& g){
            if(g.column_begin==f.column_begin&&g.column_end==f.column_end&&g.matches==f.matches)unchanged=true;
            else{queue.push(g);++stats.recovery_requeued_fragments;}
        });
        if(!unchanged)continue;
        auto recovered=MaterializeRecovery(a,f);(void)ConvertAnchor(recovered,refs,qrys,false);
        ro.add(a.reference_id,f.rb,f.re);qo.add(a.query_id,f.qb,f.qe);
        stats.recovery_reference_bases+=f.re-f.rb;stats.recovery_query_bases+=f.qe-f.qb;
        ++stats.recovery_accepted_fragments;output.push_back(std::move(recovered));
    }
    stats.recovery_seconds=Elapsed(begin);return output;
}

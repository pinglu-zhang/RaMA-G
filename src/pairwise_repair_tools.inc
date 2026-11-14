// Shared bounded-repair primitives. Keep scoring, replay and work accounting
// independent from the candidate-search driver. Sequences remain caller-owned.
void ChargeRepairWork(RepairBudget* budget,RepairWork kind,std::uint64_t amount) {
    if(budget)budget->Consume(kind,amount);
}
RepairMemoryLease ReserveRepairMemory(RepairBudget* budget,std::uint64_t bytes) {
    return budget?budget->Memory(bytes):RepairMemoryLease{};
}

template<class Check>
std::optional<std::vector<CigarOp>> AlignRepairGap(std::string_view r,std::string_view q,
        Check check,LocalRealignmentStatistics& stats,RepairBudget* budget=nullptr) {
    check();ChargeRepairWork(budget,RepairWork::LocalBases,r.size()+q.size());if(r.find('N')!=r.npos||q.find('N')!=q.npos)return {};
    if(r==q){std::vector<CigarOp> c;AppendLinkOp(c,'=',r.size());return c;}
    if(r.empty()||q.empty())return std::vector<CigarOp>{{r.empty()?'I':'D',r.size()+q.size()}};
    // Large windows are first split by searched exact anchors. Never synthesize I+D.
    if(r.size()>1000||q.size()>1000)return {};
    const auto work=3*(r.size()+1)*(q.size()+1);
    if(work>stats.work_budget-stats.work_reserved)return {};
    ChargeRepairWork(budget,RepairWork::DpWork,work);
    auto memory=ReserveRepairMemory(budget,2*work+256*(r.size()+q.size()+2));
    stats.work_reserved+=work;++stats.ksw_calls;
    auto a=pairwise_detail::globalAlignKSW2Result({r,false},{q,false});
    const auto sum=pairwise_detail::summarizeCigar(a.cigar);
    // Local validation below replays every operation against the original FASTA.
    pairwise_detail::Anchor packed(0,0,r.size(),0,0,q.size(),Strand::Forward,
                                  sum.alignment_length,0,std::move(a.cigar));
    std::vector<CigarOp> out;Length rp=0,qp=0,matches=0,paired=0;
    for(const auto unit:packed.cigar) {
        const auto op=unit&15;const Length n=unit>>4;
        ChargeRepairWork(budget,RepairWork::ReplayColumns,n);
        if(op==1){AppendLinkOp(out,'I',n);qp+=n;}
        else if(op==2){AppendLinkOp(out,'D',n);rp+=n;}
        else for(Length k=0;k<n;++k) {
            if(rp>=r.size()||qp>=q.size())throw AlignmentError("local gap CIGAR overrun");
            const bool exact=r[rp++]==q[qp++];matches+=exact;++paired;
            AppendLinkOp(out,exact?'=':'X',1);
        }
    }
    if(rp!=r.size()||qp!=q.size())throw AlignmentError("local gap CIGAR consumption");
    if(!paired||!matches||static_cast<double>(matches)/paired<.98)return {};
    check();return out;
}

// Internal repair paths borrow input that was normalized and validated before B0.
// Rebuilding the ID lookup must not rescan every genome character per repair.
Sequences RepairSequenceLookup(std::span<const SequenceRecord> records,const PairwiseCoreOptions& options) {
    if(!options.repair_budget)return CheckSequences(records);
    Sequences out;
    for(const auto& record:records) {
        KernelCheck(options,"repair-sequence-lookup");
        if(record.bases.empty()||record.size()>static_cast<Length>(INT64_MAX/256))
            throw AlignmentError("pairwise empty/oversized contig");
        if(!out.emplace(record.numeric_id,&record).second)
            throw AlignmentError("pairwise duplicate numeric sequence ID");
    }
    return out;
}

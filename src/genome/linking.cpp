// Pairwise alignment kernel. Source attribution and MIT notice are retained
// in LICENSE. Do not substitute unrelated gap backends.
#include "ramag/pairwise_core.hpp"
#include "ksw2.h"
#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <chrono>
#include <cstdlib>
#include <exception>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <optional>
#include <queue>
#include <set>
#include <type_traits>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <stdexcept>
#include <tuple>
#include <unordered_map>
#include <variant>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "ramag/internal/genome.hpp"
namespace ramag {
namespace {
using namespace genome_detail;
struct GuardedGap {
    size_t left{}, right{};
    SequenceId rid{}, qid{};
    Length rb{}, re{}, qb{}, qe{};
    Strand strand{};
    double minimum_identity{};
};

// Main-path linking consumes the fixed two-sided selection, never candidates.
// CIGAR traversal is reference-forward; query positions are canonicalized before
// constructing edges. Mutual rank adjacency proves both open intervals are free.
struct MainPathScan {
    Length matches{}, columns{}, paired{}, begin_exact{}, end_exact{};
    std::int64_t score{};
    std::uint64_t edits{};
};

void AppendLinkOp(std::vector<CigarOp>& cigar, char op, Length n) {
    if (!n) return;
    if (!cigar.empty() && cigar.back().operation == op) {
        if (n > UINT64_MAX - cigar.back().length) throw AlignmentError("link CIGAR overflow");
        cigar.back().length += n;
    } else cigar.push_back({op, n});
}

MainPathScan ScanLinkRecord(const AlignmentRecord& a, const Sequences& refs,
                           const Sequences& queries) {
    if (!refs.contains(a.reference_id) || !queries.contains(a.query_id) ||
        (a.strand != Strand::Forward && a.strand != Strand::Reverse))
        throw AlignmentError("link invalid sequence ID or strand");
    const auto& r = refs.at(a.reference_id)->bases;
    const auto& q = queries.at(a.query_id)->bases;
    if (a.reference_begin > a.reference_end || a.reference_end > r.size() ||
        a.query_begin > a.query_end || a.query_end > q.size() || a.cigar.empty())
        throw AlignmentError("link invalid record extent");
    MainPathScan s;
    Length rp = a.reference_begin, qp = 0;
    const Length qn = a.query_end - a.query_begin;
    char previous = 0;
    bool leading = true;
    for (const auto& op : a.cigar) {

        if (!op.length || (op.operation != '=' && op.operation != 'X' &&
                          op.operation != 'I' && op.operation != 'D') ||
            op.length > UINT64_MAX - s.columns)
            throw AlignmentError("link invalid extended CIGAR");
        if ((op.operation != 'I' && op.length > a.reference_end - rp) ||
            (op.operation != 'D' && op.length > qn - qp))
            throw AlignmentError("link CIGAR overrun");
        s.columns += op.length;
        if (op.operation == 'I' || op.operation == 'D') {
            // Opening is charged once per maximal same-direction gap, including
            // gaps that used to be split across two input records.
            s.score -= (previous == op.operation ? 0 : 40) +
                       3 * static_cast<std::int64_t>(op.length);
            s.edits += op.length;
            if (op.operation == 'I') qp += op.length;
            else rp += op.length;
            s.end_exact = 0;
            leading = false;
        } else {
            for (Length k = 0; k < op.length; ++k) {
                const char x = r[rp++];
                const char y = a.strand == Strand::Forward ? q[a.query_begin + qp] :
                                      Complement(q[a.query_end - qp - 1]);
                ++qp;
                const bool exact = x == y && x != 'N';
                if ((op.operation == '=') != exact)
                    throw AlignmentError("link EQX contradicts input sequence");
                const auto xi = pairwise_detail::ScoreChar2Idx[static_cast<unsigned char>(x)];
                const auto yi = pairwise_detail::ScoreChar2Idx[static_cast<unsigned char>(y)];
                s.score += pairwise_detail::dna5_simd_mat[xi * 5 + yi];
            }
            s.paired += op.length;
            if (op.operation == '=') {
                s.matches += op.length;
                s.end_exact += op.length;
                if (leading) s.begin_exact += op.length;
            } else {
                s.edits += op.length;
                s.end_exact = 0;
                leading = false;
            }
        }
        previous = op.operation;
    }
    if (rp != a.reference_end || qp != qn)
        throw AlignmentError("link incomplete CIGAR consumption");
    return s;
}

// This checks consumed spans, including I/D, independently on each axis.
// Flags and paired-base coverage are not substitutes for the span contract.
std::vector<size_t> LinkAxisOrder(const std::vector<AlignmentRecord>& records, bool ref) {
    std::vector<size_t> order(records.size());
    std::iota(order.begin(), order.end(), 0);
    const auto key = [&](size_t i) {
        const auto& a = records[i];
        return std::tuple{ref ? a.reference_id : a.query_id,
                          ref ? a.reference_begin : a.query_begin,
                          ref ? a.reference_end : a.query_end, i};
    };
    std::sort(order.begin(), order.end(), [&](size_t i, size_t j) { return key(i) < key(j); });
    for (size_t i = 1; i < order.size(); ++i) {
        const auto& a = records[order[i-1]];
        const auto& b = records[order[i]];
        if ((ref ? a.reference_id == b.reference_id : a.query_id == b.query_id) &&
            (ref ? a.reference_end > b.reference_begin : a.query_end > b.query_begin))
            throw AlignmentError(ref ? "link overlapping reference main path" :
                                       "link overlapping query main path");
    }
    return order;
}

MainPathLinkStatistics LinkMainPath(
        std::span<const SequenceRecord> reference_sequences,
        std::span<const SequenceRecord> query_sequences,
        std::vector<AlignmentRecord>& records, LinkGapMode mode, const PairwiseCoreOptions& o) {
    if ((mode != LinkGapMode::None && mode != LinkGapMode::Bounded) ||
        !o.threads || !o.min_cluster || !o.gap_fill_min_match || o.gap_fill_work_budget > 500000000)
        throw AlignmentError("invalid main-path link options");
    const auto start = KernelClock::now();
    KernelProgress(o, "mainpath-link", 0, records.size());
    const auto refs = CheckSequences(reference_sequences);
    const auto qrys = CheckSequences(query_sequences);
    MainPathLinkStatistics stats;
    stats.input_records = records.size();
    const auto capped_bases = [&](const Sequences& seqs) {
        std::uint64_t n = 0;
        for (const auto& [id, seq] : seqs) {
            (void)id;
            n += std::min(seq->size(), o.gap_fill_work_budget - n);
        }
        return n;
    };
    stats.work_budget = std::min(capped_bases(refs), capped_bases(qrys));
    std::vector<MainPathScan> scans(records.size());
    for (size_t i = 0; i < records.size(); ++i) {
        KernelCheck(o, "mainpath-link-validation");
        const auto& a = records[i];
        if (a.reference_begin >= a.reference_end || a.query_begin >= a.query_end)
            throw AlignmentError("link empty selected span");
        scans[i] = ScanLinkRecord(a, refs, qrys);
    }
    const auto ro = LinkAxisOrder(records, true);
    const auto qo = LinkAxisOrder(records, false);
    std::vector<size_t> qrank(records.size());
    for (size_t i = 0; i < qo.size(); ++i) qrank[qo[i]] = i;

    struct Edge {
        GuardedGap gap;
        unsigned kind{}; // 0 rejected, 1 zero, 2 exact, 3 single-sided, 4 DP
        std::uint64_t work{};
        std::vector<CigarOp> cigar;
    };
    std::vector<Edge> edges(ro.empty() ? 0 : ro.size()-1);
    std::vector<size_t> jobs;
    for (size_t i = 0; i < edges.size(); ++i) {
        KernelCheck(o, "mainpath-link-plan");
        const auto li = ro[i], ri = ro[i+1];
        const auto& a = records[li]; const auto& b = records[ri];
        auto& e = edges[i];
        if (a.reference_id != b.reference_id) continue;
        ++stats.adjacent_pairs;
        const bool forward = a.strand == Strand::Forward;
        if (a.query_id != b.query_id || a.strand != b.strand ||
            (forward ? qrank[li]+1 != qrank[ri] : qrank[ri]+1 != qrank[li])) {
            ++stats.geometry_rejected; continue;
        }
        const Length qlen = qrys.at(a.query_id)->size();
        // Canonical query coordinate increases with reference for either strand.
        const Length a_end = forward ? a.query_end : qlen-a.query_begin;
        const Length b_begin = forward ? b.query_begin : qlen-b.query_end;
        if (a.reference_end > b.reference_begin || a_end > b_begin) {
            ++stats.geometry_rejected; continue;
        }
        auto& g = e.gap;
        g.left=li; g.right=ri; g.rid=a.reference_id; g.qid=a.query_id; g.strand=a.strand;
        g.rb=a.reference_end; g.re=b.reference_begin;
        g.qb=forward ? a_end : qlen-b_begin;
        g.qe=forward ? b_begin : qlen-a_end;
        const Length rn=g.re-g.rb, qn=g.qe-g.qb;
        if (!rn && !qn) { e.kind=1; continue; }
        if (mode == LinkGapMode::None) continue;
        if (rn > 10000 || qn > 10000) { ++stats.length_rejected; continue; }
        const auto& x=scans[li]; const auto& y=scans[ri];
        const double xi=x.columns ? static_cast<double>(x.matches)/x.columns : 0;
        const double yi=y.columns ? static_cast<double>(y.matches)/y.columns : 0;
        if (x.matches < o.min_cluster || y.matches < o.min_cluster ||
            x.end_exact < o.gap_fill_min_match || y.begin_exact < o.gap_fill_min_match ||
            xi < .98 || yi < .98) { ++stats.flank_rejected; continue; }
        g.minimum_identity=std::max(.98, std::min(xi,yi));
        const auto r=std::string_view(refs.at(g.rid)->bases).substr(g.rb,rn);
        const auto q=std::string_view(qrys.at(g.qid)->bases).substr(g.qb,qn);
        if (r.find('N')!=r.npos || q.find('N')!=q.npos) { ++stats.n_rejected; continue; }
        bool exact=rn==qn;
        for (size_t k=0; exact && k<rn; ++k)
            exact=r[k]==(forward ? q[k] : Complement(q[qn-k-1]));
        if (exact) { e.kind=2; e.cigar.push_back({'=',rn}); continue; }
        if (rn > 1000 || qn > 1000) { ++stats.length_rejected; continue; }
        if (!rn || !qn) { e.kind=3; e.cigar.push_back({rn ? 'D' : 'I',rn+qn}); continue; }
        e.work=3*(rn+1)*(qn+1);
        jobs.push_back(i);
    }
    // Reserve the same conservative three-grid work measure as guarded gap-fill.
    // Stable ranking is independent of thread scheduling and no work is refunded.
    std::stable_sort(jobs.begin(),jobs.end(),[&](size_t a,size_t b) {
        const auto& x=edges[a]; const auto& y=edges[b];
        if (x.gap.minimum_identity != y.gap.minimum_identity)
            return x.gap.minimum_identity > y.gap.minimum_identity;
        if (x.work != y.work) return x.work < y.work;
        return a<b; // reference order, including stable original IDs
    });
    size_t admitted=0;
    for (const auto i:jobs) {
        if (edges[i].work > stats.work_budget-stats.work_reserved) { ++stats.budget_skipped; continue; }
        stats.work_reserved+=edges[i].work;
        jobs[admitted++]=i;
    }
    jobs.resize(admitted);

    stats.planning_seconds=Elapsed(start);
    struct Proposal {
        std::vector<CigarOp> cigar;
        bool accepted{}, closure{}, memory{};
        double seconds{};
    };
    constexpr size_t batch_size=64;
    for (size_t offset=0; offset<jobs.size(); offset+=batch_size) {
        const size_t count=std::min(batch_size,jobs.size()-offset);
        std::vector<Proposal> proposals(count);
        ParallelKernel(count,o,[&](size_t j) {
            KernelCheck(o,"mainpath-link-gap");
            const auto& g=edges[jobs[offset+j]].gap;
            auto& p=proposals[j];
            const auto r=std::string_view(refs.at(g.rid)->bases).substr(g.rb,g.re-g.rb);
            const auto q=std::string_view(qrys.at(g.qid)->bases).substr(g.qb,g.qe-g.qb);
            const auto begin=KernelClock::now();
            try {
                auto a=pairwise_detail::globalAlignKSW2Result({r,false},{q,g.strand==Strand::Reverse});
                const auto summary=pairwise_detail::summarizeCigar(a.cigar);
                pairwise_detail::Anchor gap(g.rid,g.rb,r.size(),g.qid,g.qb,q.size(),g.strand,
                                          summary.alignment_length,0,std::move(a.cigar));
                auto aligned=ConvertAnchor(gap,refs,qrys);
                const auto s=ScanLinkRecord(aligned,refs,qrys);
                // Never clip a fixed-endpoint result or synthesize unaligned I+D.
                if (s.paired && s.matches >= o.min_cluster && s.columns &&
                    static_cast<double>(s.matches)/s.columns >= g.minimum_identity) {
                    p.cigar=std::move(aligned.cigar); p.accepted=true;
                }
            } catch (const std::bad_alloc&) { p.memory=true; }
              catch (const AlignmentError&) { p.closure=true; }
            p.seconds=Elapsed(begin);
            KernelCheck(o,"mainpath-link-gap");
        });
        for (size_t j=0; j<count; ++j) {
            KernelCheck(o,"mainpath-link-accept");
            auto& p=proposals[j]; auto& e=edges[jobs[offset+j]];
            ++stats.ksw_calls; stats.ksw_seconds+=p.seconds;
            if (p.accepted) { e.kind=4; e.cigar=std::move(p.cigar); }
            else if (p.memory) ++stats.memory_rejected;
            else if (p.closure) ++stats.closure_rejected;
            else ++stats.quality_rejected;
        }
    }

    const auto assembly_start=KernelClock::now();
    std::vector<AlignmentRecord> output;
    output.reserve(records.size());
    for (size_t begin=0; begin<ro.size();) {
        KernelCheck(o,"mainpath-link-assembly");
        size_t end=begin;
        while (end<edges.size() && edges[end].kind) ++end;
        if (begin==end) { output.push_back(records[ro[begin]]); ++begin; continue; }
        auto merged=records[ro[begin]];
        merged.cigar.clear();
        Length expected_pairs=0, added_r=0, added_q=0, added_pairs=0;
        for (size_t k=begin; k<=end; ++k) {
            const auto& a=records[ro[k]];
            if (k>begin) {
                const auto& e=edges[k-1];
                added_r+=e.gap.re-e.gap.rb; added_q+=e.gap.qe-e.gap.qb;
                for (const auto& op:e.cigar) {
                    AppendLinkOp(merged.cigar,op.operation,op.length);
                    if (op.operation=='=' || op.operation=='X') added_pairs+=op.length;
                }
            }
            for (const auto& op:a.cigar) AppendLinkOp(merged.cigar,op.operation,op.length);
            expected_pairs+=scans[ro[k]].paired;
            merged.primary=merged.primary || a.primary;
        }
        const auto& last=records[ro[end]];
        merged.reference_end=last.reference_end;
        merged.query_begin=std::min(merged.query_begin,last.query_begin);
        merged.query_end=std::max(merged.query_end,last.query_end);
        bool closed=false;
        try {
            const auto scan=ScanLinkRecord(merged,refs,qrys);
            if (scan.paired!=expected_pairs+added_pairs) throw AlignmentError("link lost paired columns");
            merged.score=scan.score; merged.edit_distance=scan.edits;
            closed=true;
        } catch (const AlignmentError&) { ++stats.closure_rejected; }
        if (closed) {
            output.push_back(std::move(merged));
            stats.added_reference_bases+=added_r; stats.added_query_bases+=added_q;
            stats.added_paired_columns+=added_pairs;
            for (size_t k=begin;k<end;++k) {
                switch(edges[k].kind) {
                case 1: ++stats.zero_gap_links; break;
                case 2: ++stats.exact_gap_links; break;
                case 3: ++stats.single_gap_links; break;
                case 4: ++stats.dp_gap_links; break;
                default: throw AlignmentError("link invalid accepted edge");
                }
            }
        } else for (size_t k=begin;k<=end;++k) output.push_back(records[ro[k]]);
        begin=end+1;
    }
    (void)LinkAxisOrder(output,true);
    (void)LinkAxisOrder(output,false);
    // Preserve the primary's path lineage, including reverse chains. There is
    // exactly one primary per query if the fixed selection provided one.
    std::stable_sort(output.begin(),output.end(),[](const auto& a,const auto& b) {
        return std::tie(a.query_id,a.query_begin,a.reference_id,a.reference_begin,a.strand,a.query_end,a.reference_end) <
               std::tie(b.query_id,b.query_begin,b.reference_id,b.reference_begin,b.strand,b.query_end,b.reference_end);
    });
    KernelCheck(o,"mainpath-link-publication");
    stats.output_records=output.size();
    stats.assembly_seconds=Elapsed(assembly_start);
    stats.seconds=Elapsed(start);
    KernelProgress(o,"mainpath-link",records.size(),records.size());
    records=std::move(output);
    return stats;
}

}
MainPathLinkStatistics LinkSelectedAlignmentRecords(
        std::span<const SequenceRecord> references, std::span<const SequenceRecord> queries,
        std::vector<AlignmentRecord>& records, LinkGapMode mode, const PairwiseCoreOptions& options) {
    return LinkMainPath(references, queries, records, mode, options);
}
}

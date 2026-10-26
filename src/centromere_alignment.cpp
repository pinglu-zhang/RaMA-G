#include "ramag/centromere.hpp"
#include "centromere_detail.hpp"
#include <algorithm>
#include <chrono>
#include <climits>
#include <exception>
#include <memory>
#include <string>
#include <string_view>
#include <vector>
#if RAMAG_HAVE_CENTROMERE
#include "anchor.h"
#include "rama_wfa_memory.hpp"
extern "C" {
#include "wavefront/wavefront_align.h"
}
#endif

namespace ramag {
bool CentromereBackendAvailable() noexcept { return RAMAG_HAVE_CENTROMERE != 0; }
bool CentromereCliAnchorsOnly() noexcept { return RAMAG_CENTROMERE_ANCHORS_ONLY != 0; }
#if RAMAG_HAVE_CENTROMERE
namespace {
using Clock = std::chrono::steady_clock;
constexpr bool IsCanonicalBase(char c) noexcept {
  return c=='A' || c=='C' || c=='G' || c=='T';
}
double Seconds(Clock::time_point begin) {
  return std::chrono::duration<double>(Clock::now()-begin).count();
}
void Append(std::vector<CigarOp>& cigar, char op, Length len) {
  if (!len) return;
  if (!cigar.empty() && cigar.back().operation==op) cigar.back().length+=len;
  else cigar.push_back({op,len});
}
void Check(const CentromereAlignmentOptions& options, std::string_view stage) {
  if (options.interruption_callback) options.interruption_callback(stage);
}
// Unknowns have different byte values on the two sides. The pinned WFA uses
// literal equality, so N/N receives the same penalty as any other mismatch.
std::string Encode(std::string_view s, char unknown) {
  std::string result(s);
  for (char& c:result) if (!IsCanonicalBase(c)) c=unknown;
  return result;
}
} // namespace
std::vector<CigarOp> detail::AlignRamaGap(std::string_view r, std::string_view q,
                         const CentromereAlignmentOptions& options,std::uint64_t memory_limit_bytes) {
  Check(options,"centromere-gap");
  std::vector<CigarOp> out;
  if (r.empty()) { Append(out,'I',q.size()); return out; }
  if (q.empty()) { Append(out,'D',r.size()); return out; }
  if ((r.size()<=5 && q.size()>100) || (q.size()<=5 && r.size()>100)) {
    const auto n=std::min(r.size(),q.size());
    for (std::size_t i=0;i<n;++i)
      Append(out,r[i]==q[i] && IsCanonicalBase(r[i]) ? '=' : 'X',1);
    Append(out,'I',q.size()-n); Append(out,'D',r.size()-n);
    return out;
  }
  auto er=Encode(r,'\2'), eq=Encode(q,'\3');
  detail::WfaMemoryScope owned_allocations([&options]{Check(options,"centromere-gap");},memory_limit_bytes);
  auto attributes=wavefront_aligner_attr_default;
  attributes.distance_metric=gap_affine_2p;
  attributes.affine2p_penalties={0,3,4,2,12,1};
  attributes.memory_mode=wavefront_memory_med;
  // Explicitly lock the defaults inherited by the pinned original RaMA.
  attributes.alignment_scope=compute_alignment;
  attributes.alignment_form.span=alignment_end2end;
  attributes.heuristic.strategy=wf_heuristic_wfadaptive;
  attributes.heuristic.min_wavefront_length=10;
  attributes.heuristic.max_distance_threshold=50;
  attributes.heuristic.steps_between_cutoffs=1;
  attributes.system.max_num_threads=1;
  std::unique_ptr<wavefront_aligner_t,decltype(&wavefront_aligner_delete)>
      aligner(wavefront_aligner_new(&attributes),wavefront_aligner_delete);
  if (!aligner) throw std::bad_alloc();
  const int status=wavefront_align(aligner.get(),er.data(),static_cast<int>(er.size()),
                                  eq.data(),static_cast<int>(eq.size()));
  if (status!=0) throw AlignmentError("RaMA WFA2 failed, status="+std::to_string(status));
  std::size_t ri=0, qi=0;
  const auto* cg=aligner->cigar;
  if (!cg || cg->begin_offset<0 || cg->end_offset<cg->begin_offset ||
      cg->end_offset>cg->max_operations)
    throw AlignmentError("invalid WFA2 CIGAR buffer");
  for (int i=cg->begin_offset;i<cg->end_offset;++i) {
    const char op=cg->operations[i];
    if (op=='M' || op=='=' || op=='X') {
      if (ri>=r.size() || qi>=q.size()) throw AlignmentError("WFA2 paired column outside interval");
      Append(out,r[ri]==q[qi] && IsCanonicalBase(r[ri]) ? '=' : 'X',1);
      ++ri; ++qi;
    } else if (op=='I') { Append(out,'I',1); ++qi; }
    else if (op=='D') { Append(out,'D',1); ++ri; }
    else throw AlignmentError("unknown WFA2 operation");
  }
  if (ri!=r.size() || qi!=q.size()) throw AlignmentError("WFA2 incomplete interval consumption");
  Check(options,"centromere-gap");
  return out;
}
namespace {
void Validate(AlignmentRecord& a, const SequenceRecord& r,const SequenceRecord& q) {
  std::size_t ri=0, qi=0;
  std::int64_t cost=0;
  for (const auto& op:a.cigar) {
    if (!op.length) throw AlignmentError("zero-length centromere CIGAR");
    if (op.operation=='=' || op.operation=='X') {
      for (Length i=0;i<op.length;++i) {
        if (ri>=r.size() || qi>=q.size()) throw AlignmentError("centromere paired column out of bounds");
        const bool exact=IsCanonicalBase(r.bases[ri]) && r.bases[ri]==q.bases[qi];
        if (exact!=(op.operation=='=')) throw AlignmentError("centromere CIGAR/base mismatch");
        ++ri; ++qi;
      }
      if (op.operation=='X') cost+=3*static_cast<std::int64_t>(op.length);
    } else {
      const auto n=static_cast<std::int64_t>(op.length);
      cost+=std::min(4+2*n,12+n);
      if (op.operation=='I') qi+=op.length;
      else if (op.operation=='D') ri+=op.length;
      else throw AlignmentError("invalid centromere CIGAR operation");
    }
    if (ri>r.size() || qi>q.size()) throw AlignmentError("centromere CIGAR outside sequence");
  }
  if (ri!=r.size() || qi!=q.size()) throw AlignmentError("centromere CIGAR not end-to-end");
  a.score=-cost;
  a.edit_distance=0;
  for (const auto& op:a.cigar) if (op.operation!='=') a.edit_distance+=op.length;
}
}
#else
std::vector<CigarOp> detail::AlignRamaGap(std::string_view, std::string_view,
    const CentromereAlignmentOptions&,std::uint64_t) {
  throw DependencyError("centromere backend not built; enable RAMAG_ENABLE_CENTROMERE on Linux");
}
#endif

static AlignmentResult RunCentromereCore(const SequenceRecord& reference,
    const SequenceRecord& query,const CentromereAlignmentOptions& options,
    CentromereAnchorResult* anchor_output) {
#if !RAMAG_HAVE_CENTROMERE
  (void)reference; (void)query; (void)options; (void)anchor_output;
  throw DependencyError("centromere backend not built; enable RAMAG_ENABLE_CENTROMERE on Linux");
#else
  if (!options.worker_threads) throw CentromereInputError("centromere threads must be positive");
  if (reference.bases.empty() || query.bases.empty())
    throw CentromereInputError("centromere sequences must be nonempty");
  constexpr auto safe_length=static_cast<Length>((INT_MAX-64)/3);
  if (reference.size()>safe_length || query.size()>safe_length-reference.size())
    throw CentromereInputError("sequence pair exceeds pinned RaMA/WFA2 safe 32-bit coordinate/score range");
  for (const auto* s:{&reference,&query})
    for (char c:s->bases)
      if (!IsCanonicalBase(c) && c!='N')
        throw CentromereInputError("AlignCentromeres requires normalized A/C/G/T/N");
  const auto start=Clock::now();
  AlignmentResult result;
  auto& stats=result.statistics;
  stats.reference_contigs=stats.query_contigs=1;
  stats.reference_bases=reference.size(); stats.query_bases=query.size();
  stats.reference_ambiguous_bases=static_cast<Length>(std::count(reference.bases.begin(),reference.bases.end(),'N'));
  stats.query_ambiguous_bases=static_cast<Length>(std::count(query.bases.begin(),query.bases.end(),'N'));
  stats.actual_seed_route="rama-source-rare-recursive-v1";
  stats.chaining_route="rama-source-collinear-v1";
  stats.seed_requested_threads=options.worker_threads;
  stats.seed_scheduled_threads=stats.seed_worker_threads=options.worker_threads;
  stats.seed_parallel_route="rama-per-call-recursion-pool";
  if (options.progress_callback) options.progress_callback("centromere-anchors",0,0);
  Check(options,"centromere-anchors");
  rama_source::RareMatchPairs anchors;
  std::vector<rama_source::SequenceInfo> sequences{
      {reference.bases,reference.name,static_cast<uint_t>(reference.size())},
      {query.bases,query.name,static_cast<uint_t>(query.size())}};
  {
    rama_source::AnchorFinder finder(sequences,"",
        options.worker_threads==1 ? 0 : options.worker_threads,false,false,100,
        [&options]{Check(options,"centromere-anchors");});
    anchors=finder.lanuchAnchorSearching([&](const rama_source::RareMatchPairs& first) {
      if (!anchor_output) return;
      for (const auto& a:first) {
        if (!a.match_length || a.second_pos<reference.size()+1)
          throw AlignmentError("invalid initial RaMA anchor");
        const Length qb=a.second_pos-reference.size()-1;
        if (a.first_pos+a.match_length>reference.size() || qb+a.match_length>query.size())
          throw AlignmentError("initial RaMA anchor outside sequence");
        for (Length i=0;i<a.match_length;++i)
          if (!IsCanonicalBase(reference.bases[a.first_pos+i]) ||
              reference.bases[a.first_pos+i]!=query.bases[qb+i])
            throw AlignmentError("initial RaMA anchor is not canonical exact match");
        anchor_output->initial.push_back({a.first_pos,qb,a.match_length,a.weight});
      }
    });
  } // Release the index before any WFA workspace is allocated.
  stats.seed_seconds=Seconds(start);
  stats.selected_seed_count=stats.merged_seed_count=anchors.size();
  stats.chain_count=anchors.empty()?0:1;
  Length re=0, qe=0;
  for (const auto& a:anchors) {
    if (!a.match_length || a.second_pos<reference.size()+1)
      throw AlignmentError("invalid RaMA anchor");
    const Length qb=a.second_pos-reference.size()-1;
    if (a.first_pos<re || qb<qe || a.first_pos+a.match_length>reference.size() ||
        qb+a.match_length>query.size()) throw AlignmentError("RaMA anchor order/bounds failure");
    for (Length i=0;i<a.match_length;++i)
      if (!IsCanonicalBase(reference.bases[a.first_pos+i]) ||
          reference.bases[a.first_pos+i]!=query.bases[qb+i])
        throw AlignmentError("RaMA anchor is not canonical exact match");
    re=a.first_pos+a.match_length; qe=qb+a.match_length;
  }
  if (anchor_output) {
    for (const auto& a:anchors)
      anchor_output->final.push_back({a.first_pos,a.second_pos-reference.size()-1,a.match_length,a.weight});
    Check(options,"centromere-anchors");
    stats.total_seconds=Seconds(start);
    anchor_output->statistics=stats;
    return result; // Deliberately before gap partitioning or any WFA allocation.
  }
  const auto gaps=rama_source::AnchorFinder::RareMatchPairs2Intervals(anchors,
      {0,static_cast<uint_t>(reference.size()),0,static_cast<uint_t>(query.size())},
      static_cast<uint_t>(reference.size()));
  AlignmentRecord record;
  record.reference_id=reference.numeric_id; record.query_id=query.numeric_id;
  record.reference_end=reference.size(); record.query_end=query.size();
  record.primary=true;
  const auto extension_start=Clock::now();
  const auto workers=std::min<std::size_t>(options.worker_threads,gaps.size());
  stats.extension_worker_threads=static_cast<std::uint32_t>(workers);
  for (std::size_t first=0;first<gaps.size();first+=64) {
    const auto count=std::min<std::size_t>(64,gaps.size()-first);
    std::vector<std::vector<CigarOp>> batch(count);
    std::vector<std::exception_ptr> errors(count);
    rama_source::ThreadPool pool(workers==1 ? 0 : workers);
    for (std::size_t j=0;j<count;++j) {
      auto task=[&,j] {
        try {
          const auto& g=gaps[first+j];
          if (static_cast<Length>(g.pos1)+g.len1>reference.size() ||
              static_cast<Length>(g.pos2)+g.len2>query.size())
            throw AlignmentError("RaMA gap outside sequence");
          batch[j]=detail::AlignRamaGap(std::string_view(reference.bases).substr(g.pos1,g.len1),
                       std::string_view(query.bases).substr(g.pos2,g.len2),options);
        } catch (...) { errors[j]=std::current_exception(); }
      };
      if (workers==1) task(); else pool.enqueue(task);
    }
    pool.waitAllTasksDone();
    for (std::size_t j=0;j<count;++j) {
      if (errors[j]) std::rethrow_exception(errors[j]);
      for (const auto& c:batch[j]) Append(record.cigar,c.operation,c.length);
      if (first+j<anchors.size()) Append(record.cigar,'=',anchors[first+j].match_length);
    }
    if (options.progress_callback) options.progress_callback("centromere-gap",first+count,gaps.size());
  }
  Validate(record,reference,query);
  result.alignments.push_back(std::move(record));
  stats.candidate_alignment_count=stats.alignment_count=1;
  stats.chain_and_extension_seconds=Seconds(extension_start);
  stats.total_seconds=Seconds(start);
  return result;
#endif
}
AlignmentResult AlignCentromeres(const SequenceRecord& reference,
    const SequenceRecord& query,const CentromereAlignmentOptions& options) {
  return RunCentromereCore(reference,query,options,nullptr);
}
CentromereAnchorResult FindCentromereAnchors(const SequenceRecord& reference,
    const SequenceRecord& query,const CentromereAlignmentOptions& options) {
  CentromereAnchorResult result;
  (void)RunCentromereCore(reference,query,options,&result);
  return result;
}
}

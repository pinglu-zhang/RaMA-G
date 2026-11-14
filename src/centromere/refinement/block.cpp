#include "detail.hpp"
#include "centromere/refinement/block_dp.hpp"
#include "centromere/stitch/kernel.hpp"
#include "centromere/stitch/fixed_gaps.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <limits>
#include <map>
#include <new>
#include <set>
#include <sstream>
#include <tuple>
#include <unordered_map>
#if RAMAG_HAVE_CENTROMERE
#include "anchor.h"
#endif

namespace ramag {
using namespace centromere_refinement_detail;
namespace centromere_refinement_detail {
// Exact 2-bit 32-mers use all 64 bits, so equality has no hash collision
// semantics. Counts saturate at two because only physical uniqueness matters.
int ContextBits(char c) {
  switch (c) {
  case 'A':
    return 0;
  case 'C':
    return 1;
  case 'G':
    return 2;
  case 'T':
    return 3;
  default:
    return -1;
  }
}
ContextCounts CountContexts(const SequenceRecord &sequence,
                            const CentromereEndpointOptions &options) {
  ContextCounts counts;
  std::uint64_t key = 0;
  unsigned run = 0;
  for (std::size_t i = 0; i < sequence.bases.size(); ++i) {
    if ((i % 4096) == 0 && options.kernel.interruption_callback)
      options.kernel.interruption_callback("centromere-block-context");
    const int bits = ContextBits(sequence.bases[i]);
    if (bits < 0) {
      key = 0;
      run = 0;
      continue;
    }
    key = (key << 2U) | static_cast<std::uint64_t>(bits);
    run = std::min(32U, run + 1);
    if (run == 32U) {
      auto &count = counts[key];
      if (count < 2)
        ++count;
    }
  }
  return counts;
}
bool DistinguishableContext(const A &anchor, const SequenceRecord &reference,
                            const ContextCounts &rc, const ContextCounts &qc) {
  if (anchor.length < 64)
    return false;
  std::uint64_t key = 0;
  Length first_end = 0;
  bool first = false;
  for (Length i = 0; i < anchor.length; ++i) {
    const int bits = ContextBits(reference.bases[anchor.reference_begin + i]);
    // Input anchors have already passed canonical exact-match validation.
    key = (key << 2U) | static_cast<std::uint64_t>(bits);
    if (i < 31)
      continue;
    const auto r = rc.find(key), q = qc.find(key);
    if (r == rc.end() || q == qc.end() || r->second != 1 || q->second != 1)
      continue;
    const Length begin = i - 31;
    if (first && begin >= first_end)
      return true;
    if (!first) {
      first = true;
      first_end = i + 1;
    }
  }
  return false;
}
bool CigarOffset(Position r, Position q, Position rb, Position qb, bool dr,
                 bool dq, Length n, Length &offset) {
  if (rb < r || qb < q)
    return false;
  offset = dr ? rb - r : qb - q;
  return offset <= n && rb - r == (dr ? offset : 0) &&
         qb - q == (dq ? offset : 0);
}
bool SliceStitchedCigar(const CentromereStitchResult &stitched,
                        std::uint64_t region, Position rb, Position re,
                        Position qb, Position qe, std::vector<CigarOp> &out) {
  for (const auto &segment : stitched.segments) {
    const auto &a = segment.alignment;
    if (segment.path_role != "primary" || segment.region_id != region ||
        a.strand != Strand::Forward || a.reference_begin > rb ||
        a.reference_end < re || a.query_begin > qb || a.query_end < qe)
      continue;
    Position r = a.reference_begin, q = a.query_begin;
    bool active = false;
    std::vector<CigarOp> sliced;
    for (auto op : a.cigar) {
      if (op.operation != '=' && op.operation != 'X' && op.operation != 'I' &&
          op.operation != 'D')
        throw AlignmentError("block review requires extended stitched CIGAR");
      const bool dr = op.operation != 'I', dq = op.operation != 'D';
      Length begin = 0, end = 0;
      if (!active) {
        if (!CigarOffset(r, q, rb, qb, dr, dq, op.length, begin)) {
          r += dr ? op.length : 0;
          q += dq ? op.length : 0;
          continue;
        }
        active = true;
      }
      if (CigarOffset(r, q, re, qe, dr, dq, op.length, end) && end >= begin) {
        Emit(sliced, op.operation, end - begin);
        out = std::move(sliced);
        return true;
      }
      Emit(sliced, op.operation, op.length - begin);
      r += dr ? op.length : 0;
      q += dq ? op.length : 0;
    }
  }
  return false;
}
} // namespace centromere_refinement_detail

CentromereEndpointResult
RefineCentromereAnchorBlocks(const SequenceRecord &r, const SequenceRecord &q,
                             const CentromereDirectedAnchorResult &input,
                             const CentromereCopyAnchorResult *copies,
                             const CentromereStitchResult &stitched,
                             const CentromereEndpointOptions &options) {
  if (!input.complete || (copies && !copies->complete) || stitched.aborted ||
      stitched.reference_length != r.size() ||
      stitched.query_length != q.size())
    throw AlignmentError("block review requires valid completed anchor stages");
  if (!options.max_reviews || !std::isfinite(options.stage_seconds_limit) ||
      options.stage_seconds_limit <= 0 ||
      !std::isfinite(options.window_seconds_limit) ||
      options.window_seconds_limit <= 0 ||
      options.kernel.score_model !=
          CentromereStitchScoreModel::Centrolign3Piece ||
      stitched.score_model != options.kernel.score_model)
    throw AlignmentError("invalid block review scoring/resources");
  const auto started = Clock::now();
  const auto check = [&] {
    if (options.kernel.interruption_callback)
      options.kernel.interruption_callback("centromere-block-review");
  };
  check();
  CentromereEndpointResult out;
  out.anchors = input;
  std::map<std::uint64_t, std::vector<const A *>> paths;
  std::set<std::uint64_t> ids;
  for (const auto &a : input.final) {
    CheckAnchor(a, r, q);
    if (!ids.insert(a.id).second)
      throw AlignmentError("duplicate block review anchor ID");
    paths[a.region_id].push_back(&a);
  }
  for (auto &[id, path] : paths) {
    (void)id;
    std::sort(path.begin(), path.end(),
              [](auto a, auto b) { return a->chain_index < b->chain_index; });
    for (std::size_t i = 0; i < path.size(); ++i) {
      const auto *a = path[i];
      if (a->chain_index != i + 1)
        throw AlignmentError("block review chain index");
      if (i) {
        const auto *b = path[i - 1];
        if (a->strand != b->strand ||
            b->reference_begin + b->length > a->reference_begin ||
            (a->strand == Strand::Forward
                 ? b->query_begin + b->length > a->query_begin
                 : a->query_begin + a->length > b->query_begin))
          throw AlignmentError("block review chain order");
      }
    }
  }
  const auto rc = CountContexts(r, options), qc = CountContexts(q, options);
  std::vector<Replacement> accepted;
  for (const auto &[region, path] : paths) {
    if (path.empty() || path.front()->strand != Strand::Forward)
      continue;
    std::vector<std::size_t> strong;
    for (std::size_t i = 0; i < path.size(); ++i) {
      check();
      if (DistinguishableContext(*path[i], r, rc, qc))
        strong.push_back(i);
    }
    for (std::size_t k = 1; k < strong.size(); ++k) {
      check();
      const std::size_t li = strong[k - 1], zi = strong[k];
      // A single-root window is already handled by the original stage.
      if (zi - li < 3)
        continue;
      const auto &l = *path[li], &z = *path[zi];
      CentromereEndpointReview review;
      review.id = out.reviews.size() + 1;
      review.region_id = region;
      review.left_anchor_id = l.id;
      review.anchor_id = path[li + 1]->id;
      review.right_anchor_id = z.id;
      review.reference_begin = l.reference_begin + l.length;
      review.reference_end = z.reference_begin;
      review.query_begin = l.query_begin + l.length;
      review.query_end = z.query_begin;
      review.left_gap =
          std::max(path[li + 1]->reference_begin - review.reference_begin,
                   path[li + 1]->query_begin - review.query_begin);
      review.right_gap = std::max(
          review.reference_end - path[zi - 1]->reference_begin -
              path[zi - 1]->length,
          review.query_end - path[zi - 1]->query_begin - path[zi - 1]->length);
      const Length rn = review.reference_end - review.reference_begin,
                   qn = review.query_end - review.query_begin;
      if (zi - li - 1 > 256 || rn > 65536 || qn > 65536)
        review.reason = "block_geometry_budget";
      bool has_root = false;
      for (const auto &a : input.initial)
        if (a.region_id == region && a.strand == Strand::Forward &&
            a.reference_begin >= review.reference_begin &&
            a.reference_begin + a.length <= review.reference_end &&
            a.query_begin >= review.query_begin &&
            a.query_begin + a.length <= review.query_end) {
          has_root = true;
          break;
        }
      if (review.reason.empty() && !has_root)
        review.reason = "no_root_in_block";
      if (review.reason.empty() && copies)
        for (const auto &a : copies->anchors)
          if (Touch(a.reference_begin, a.reference_begin + a.length,
                    review.reference_begin, review.reference_end) ||
              Touch(a.query_begin, a.query_begin + a.length, review.query_begin,
                    review.query_end)) {
            review.reason = "accepted_copy_footprint";
            break;
          }
      if (review.reason.empty())
        for (const auto &a : accepted)
          if (Touch(a.rb, a.re, review.reference_begin, review.reference_end) ||
              Touch(a.qb, a.qe, review.query_begin, review.query_end)) {
            review.reason = "overlapping_accepted_window";
            break;
          }
      if (review.reason.empty() &&
          !SliceStitchedCigar(stitched, region, review.reference_begin,
                              review.reference_end, review.query_begin,
                              review.query_end, review.old_cigar))
        review.reason = "unresolved_or_missing_old_path";
      const auto rv = std::string_view(r.bases).substr(review.reference_begin,
                                                       rn),
                 qv = std::string_view(q.bases).substr(review.query_begin, qn);
      if (review.reason.empty()) {
        review.old_score = detail::ScoreCentromereStitchCigar(
            rv, qv, review.old_cigar, options.kernel.score_model);
        for (auto op : review.old_cigar)
          if (op.operation == '=')
            review.old_exact_bases += op.length;
        if (review.old_exact_bases >= std::min(rn, qn))
          review.reason = "no_exact_headroom";
        else if (rv == qv)
          review.reason = "identical_window_ancestry_ambiguous";
      }
      if (review.reason.empty() &&
          (out.kernel_calls >= options.max_reviews ||
           Seconds(started) >= options.stage_seconds_limit))
        review.reason = "review_budget";
      if (!review.reason.empty()) {
        out.reviews.push_back(std::move(review));
        continue;
      }
      auto kernel = options.kernel;
      kernel.progress_callback = {};
      kernel.observation_callback = {};
      kernel.gap_seconds_limit =
          std::min({kernel.gap_seconds_limit, options.window_seconds_limit,
                    options.stage_seconds_limit - Seconds(started)});
      if (kernel.gap_seconds_limit <= 0) {
        review.reason = "review_budget";
        out.reviews.push_back(std::move(review));
        continue;
      }
      const auto begin = Clock::now();
      ++out.kernel_calls;
      try {
        auto proposal = detail::AlignCentromereBlockDp(rv, qv, kernel);
        review.new_score = detail::ScoreCentromereStitchCigar(
            rv, qv, proposal.cigar, kernel.score_model);
        if (review.new_score != proposal.score)
          throw AlignmentError("block kernel score/replay differs");
        review.scored = true;
        review.engine = proposal.engine;
        review.peak_kernel_bytes = proposal.peak_kernel_bytes;
        review.dp_cells = proposal.dp_cells;
        review.wfa_states = proposal.wfa_states;
        review.new_cigar = std::move(proposal.cigar);
        Replacement p{region,
                      review.reference_begin,
                      review.reference_end,
                      review.query_begin,
                      review.query_end,
                      {}};
        Position rb = p.rb, qb = p.qb;
        for (auto op : review.new_cigar) {
          if (op.operation == 'X')
            review.new_mismatch_bases += op.length;
          if (op.operation == '=') {
            p.exact.push_back({0, region, 0, rb, qb, op.length, Strand::Forward,
                               static_cast<double>(op.length)});
            review.new_exact_bases += op.length;
          }
          if (op.operation != 'I')
            rb += op.length;
          if (op.operation != 'D')
            qb += op.length;
        }
        if (review.new_mismatch_bases)
          review.reason = "new_substitutions_preserved";
        else if (review.new_score <= review.old_score)
          review.reason = review.new_score == review.old_score
                              ? "equal_score_preserved"
                              : "lower_score_preserved";
        else if (review.new_exact_bases <= review.old_exact_bases)
          review.reason = "no_exact_coverage_gain";
        else {
          review.accepted = true;
          review.reason = "strict_score_and_exact_gain";
          accepted.push_back(std::move(p));
          ++out.accepted;
        }
      } catch (const detail::StitchGapResourceLimit &error) {
        review.reason = error.reason;
        review.engine = error.engine;
        review.peak_kernel_bytes = error.peak_kernel_bytes;
        review.dp_cells = error.dp_cells;
        review.wfa_states = error.wfa_states;
      }
      review.seconds = Seconds(begin);
      out.reviews.push_back(std::move(review));
    }
  }
  if (!accepted.empty()) {
    Replace(out.anchors.initial, accepted, input.regions);
    Replace(out.anchors.final, accepted, input.regions);
    Refresh(out.anchors);
  }
  check();
  out.seconds = Seconds(started);
  return out;
}
} // namespace ramag

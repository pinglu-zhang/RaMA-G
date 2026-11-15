#include "ramag/internal/centromere/refinement.hpp"
#include "ramag/internal/centromere/gap_alignment.hpp"
#include "ramag/internal/centromere/gap_alignment.hpp"
#include "ramag/internal/centromere/gap_alignment.hpp"
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
// A cross-region baseline is the sum of actual forward CIGAR fragments already
// emitted inside the rectangle. Unknown space is free, not an invented I+D
// bridge. This favors the input over any complete candidate across that space.
Length ShortDisagreementBurden(const std::vector<CigarOp>&cigar){
  Length n=0;
  for(auto op:cigar)if(op.operation=='X')n+=op.length;
  for(std::size_t i=1;i<cigar.size();++i){
    const auto a=cigar[i-1],b=cigar[i];
    if(((a.operation=='I'&&b.operation=='D')||(a.operation=='D'&&b.operation=='I'))&&a.length<64&&b.length<64)n+=std::min(a.length,b.length);
  }
  return n;
}
CrossFragmentSummary CrossFragments(const CentromereStitchResult &stitched,
                                    const SequenceRecord &reference,
                                    const SequenceRecord &query, Position rb,
                                    Position re, Position qb, Position qe) {
  CrossFragmentSummary summary;
  for (const auto &segment : stitched.segments) {
    const auto &a = segment.alignment;
    if (segment.path_role != "primary" || a.strand != Strand::Forward ||
        a.reference_end < rb || a.reference_begin > re ||
        a.query_end < qb || a.query_begin > qe)
      continue;
    Position r = a.reference_begin, q = a.query_begin;
    Position first_r = 0, first_q = 0, last_r = 0, last_q = 0;
    bool active = false;
    std::vector<CigarOp> cigar;
    for (auto op : a.cigar) {
      if (op.operation != '=' && op.operation != 'X' && op.operation != 'I' &&
          op.operation != 'D')
        throw AlignmentError("cross-region review requires extended CIGAR");
      const bool dr = op.operation != 'I', dq = op.operation != 'D';
      Length begin = 0, end = op.length;
      bool possible = true;
      if (dr) {
        if (r < rb)
          begin = std::max(begin, rb - r);
        if (r > re)
          possible = false;
        else
          end = std::min(end, re - r);
      } else if (r < rb || r > re) {
        possible = false;
      }
      if (dq) {
        if (q < qb)
          begin = std::max(begin, qb - q);
        if (q > qe)
          possible = false;
        else
          end = std::min(end, qe - q);
      } else if (q < qb || q > qe) {
        possible = false;
      }
      if (possible && begin < end) {
        if (!active) {
          first_r = r + (dr ? begin : 0);
          first_q = q + (dq ? begin : 0);
          active = true;
        }
        Emit(cigar, op.operation, end - begin);
        last_r = r + (dr ? end : 0);
        last_q = q + (dq ? end : 0);
      } else if (active) {
        break;
      }
      r += dr ? op.length : 0;
      q += dq ? op.length : 0;
    }
    if (cigar.empty())
      continue;
    summary.short_disagreements+=ShortDisagreementBurden(cigar);
    summary.score += detail::ScoreCentromereStitchCigar(
        std::string_view(reference.bases).substr(first_r, last_r - first_r),
        std::string_view(query.bases).substr(first_q, last_q - first_q), cigar,
        CentromereStitchScoreModel::Centrolign3Piece);
    for (auto op : cigar) {
      if (op.operation == '=')
        summary.exact += op.length;
      if (op.operation == 'X')
        summary.mismatches += op.length;
    }
  }
  return summary;
}
// Outer contexts delimit a problem; they do not uniquely place an unmarked
// repeat core. Internal evidence is an actual 32-base coordinate pairing that
// is unique in both complete inputs. Overlapping starts are retained for exact
// preservation checks, never interpreted as independent observations.
CrossContexts CrossPathContexts(const SequenceRecord &reference,
    const ContextCounts &rc, const ContextCounts &qc,
    const std::vector<CigarOp> &path, Position r, Position q,
    const CentromereEndpointReview &window) {
  CrossContexts out;
  for (auto op : path) {
    if (op.operation == '=' && r < window.reference_end && q < window.query_end) {
      const Length begin = std::max(r < window.reference_begin ? window.reference_begin-r : 0,
                                   q < window.query_begin ? window.query_begin-q : 0);
      const Length end = std::min({op.length,window.reference_end-r,window.query_end-q});
      if (begin <= end && end-begin >= 32) {
        std::uint64_t key = 0;
        for (Length i=begin;i<end;++i) {
          const auto bits=ContextBits(reference.bases[r+i]);
          if(bits<0) throw AlignmentError("noncanonical reviewed exact context");
          key=(key<<2U)|static_cast<std::uint64_t>(bits);
          if(i-begin<31) continue;
          const auto a=rc.find(key),b=qc.find(key);
          if(a!=rc.end() && b!=qc.end() && a->second==1 && b->second==1)
            out.emplace(r+i-31,q+i-31);
        }
      }
    }
    if(op.operation!='I')r+=op.length;
    if(op.operation!='D')q+=op.length;
  }
  return out;
}
CrossContexts CrossOldContexts(const SequenceRecord &r,
    const ContextCounts &rc,const ContextCounts &qc,
    const CentromereStitchResult &old,const CentromereEndpointReview &window) {
  CrossContexts out;
  for(const auto &segment:old.segments) {
    const auto &a=segment.alignment;
    if(segment.path_role!="primary" || a.strand!=Strand::Forward)continue;
    auto contexts=CrossPathContexts(r,rc,qc,a.cigar,a.reference_begin,a.query_begin,window);
    out.insert(contexts.begin(),contexts.end());
  }
  return out;
}
// Cut at exact path vertices, preserving the reviewed operation order. A
// rectangle intersection is insufficient at an I/D corner: both physical
// endpoints must actually occur along this candidate's monotone path.
std::vector<CigarOp> CrossGapSlice(const std::vector<CigarOp> &path,
                                 Position r, Position q, Position rb,
                                 Position re, Position qb, Position qe) {
  bool active = false;
  std::vector<CigarOp> result;
  for (auto op : path) {
    const bool dr = op.operation != 'I', dq = op.operation != 'D';
    auto offset = [&](Position tr, Position tq) -> std::optional<Length> {
      if ((dr && (tr < r || tr - r > op.length)) || (!dr && tr != r) ||
          (dq && (tq < q || tq - q > op.length)) || (!dq && tq != q))
        return std::nullopt;
      const auto n = dr ? tr - r : tq - q;
      if (dr && dq && tr - r != tq - q)
        return std::nullopt;
      return n;
    };
    Length begin = 0;
    if (!active) {
      auto start = offset(rb, qb);
      if (start) {
        active = true;
        begin = *start;
      }
    }
    if (active) {
      auto end = offset(re, qe);
      if (end && *end >= begin) {
        Emit(result, op.operation, *end - begin);
        return result;
      }
      Emit(result, op.operation, op.length - begin);
    }
    r += dr ? op.length : 0;
    q += dq ? op.length : 0;
  }
  throw AlignmentError("reviewed gap endpoints are not on the exact candidate path");
}
std::vector<CentromereStitchGap> CrossFixedGaps(
    const SequenceRecord &r, const SequenceRecord &q,
    const CentromereDirectedAnchorResult &candidate,
    const CentromereEndpointReview &review) {
  std::map<std::uint64_t, std::vector<const A *>> paths;
  for (const auto &region : candidate.regions)
    if (region.strand == Strand::Forward &&
        Touch(region.query_begin, region.query_end,
              review.query_begin, review.query_end))
      paths.emplace(region.id, std::vector<const A *>{});
  for (const auto &a : candidate.final) {
    auto path = paths.find(a.region_id);
    if (a.strand == Strand::Forward && path != paths.end())
      path->second.push_back(&a);
  }
  std::vector<CentromereStitchGap> result;
  for (auto &[region, anchors] : paths) {
    std::sort(anchors.begin(), anchors.end(), [](const A *a, const A *b) {
      return a->chain_index < b->chain_index;
    });
    for (std::size_t i = 1; i < anchors.size(); ++i) {
      const auto &left = *anchors[i - 1], &right = *anchors[i];
      const auto rb = left.reference_begin + left.length;
      const auto qb = left.query_begin + left.length;
      const auto re = right.reference_begin, qe = right.query_begin;
      if (rb > re || qb > qe)
        throw AlignmentError("reviewed exact anchors are not ordered");
      if (rb == re && qb == qe)
        continue;
      const bool contained = rb >= review.reference_begin &&
          re <= review.reference_end && qb >= review.query_begin &&
          qe <= review.query_end;
      if (!contained) {
        // No new gap may silently fall back to the unconstrained solver.
        if (Touch(rb, re, review.reference_begin, review.reference_end) ||
            Touch(qb, qe, review.query_begin, review.query_end))
          throw AlignmentError("reviewed gap crosses the candidate window");
        continue;
      }
      CentromereStitchGap gap;
      gap.path_role = "primary";
      gap.region_id = region;
      gap.strand = Strand::Forward;
      gap.reference_begin = rb; gap.reference_end = re;
      gap.query_begin = qb; gap.query_end = qe;
      gap.left_anchor_id = left.id; gap.right_anchor_id = right.id;
      gap.cigar = CrossGapSlice(review.new_cigar, review.reference_begin,
                               review.query_begin, rb, re, qb, qe);
      gap.score = detail::ScoreCentromereStitchCigar(
          std::string_view(r.bases).substr(rb, re - rb),
          std::string_view(q.bases).substr(qb, qe - qb), gap.cigar,
          CentromereStitchScoreModel::Centrolign3Piece);
      gap.status = "complete"; gap.scored = true;
      result.push_back(std::move(gap));
    }
  }
  return result;
}
std::string CrossReverseComplement(const std::string &sequence) {
  std::string result(sequence.rbegin(), sequence.rend());
  for (auto &c : result) {
    switch (c) {
    case 'A': c = 'T'; break;
    case 'C': c = 'G'; break;
    case 'G': c = 'C'; break;
    case 'T': c = 'A'; break;
    case 'N': break;
    default: throw AlignmentError("cross-region review requires uppercase ACGTN");
    }
  }
  return result;
}
} // namespace centromere_refinement_detail

CentromereEndpointResult RefineCentromereCrossRegionBlocks(
    const SequenceRecord &r, const SequenceRecord &q,
    const CentromereDirectedAnchorResult &input,
    const CentromereCopyAnchorResult *copies,
    const CentromereStitchResult &stitched,
    const CentromereEndpointOptions &options) {
  if (!input.complete || (copies && !copies->complete) || stitched.aborted ||
      stitched.reference_length != r.size() || stitched.query_length != q.size())
    throw AlignmentError("cross-region review requires valid completed stages");
  if (!options.max_reviews || !std::isfinite(options.stage_seconds_limit) ||
      options.stage_seconds_limit <= 0 ||
      !std::isfinite(options.window_seconds_limit) ||
      options.window_seconds_limit <= 0 ||
      options.kernel.score_model != CentromereStitchScoreModel::Centrolign3Piece ||
      stitched.score_model != options.kernel.score_model ||
      !options.kernel.gap_memory_limit_bytes ||
      !std::isfinite(options.kernel.gap_seconds_limit) ||
      options.kernel.gap_seconds_limit <= 0)
    throw AlignmentError("invalid cross-region review scoring/resources");
  const auto started = Clock::now();
  const auto check = [&] {
    if (options.kernel.interruption_callback)
      options.kernel.interruption_callback("centromere-cross-region-review");
  };
  CentromereEndpointResult out;
  out.anchors = input;
  check();
  std::map<std::uint64_t, const CentromereDirectedRegion *> regions;
  std::vector<const A *> ordered;
  std::set<std::uint64_t> ids;
  for (const auto &region : input.regions)
    if (!region.id || region.query_begin >= region.query_end ||
        region.query_end > q.size() ||
        !regions.emplace(region.id, &region).second)
      throw AlignmentError("invalid cross-region input region");
  for (const auto &a : input.final) {
    CheckAnchor(a, r, q);
    auto region = regions.find(a.region_id);
    if (region == regions.end() || a.strand != region->second->strand ||
        a.query_begin < region->second->query_begin ||
        a.query_begin + a.length > region->second->query_end ||
        !ids.insert(a.id).second)
      throw AlignmentError("invalid cross-region anchor region/identity");
    if (a.strand == Strand::Forward)
      ordered.push_back(&a);
  }
  for (const auto &a : input.initial)
    CheckAnchor(a, r, q);
  std::sort(ordered.begin(), ordered.end(), [](const A *a, const A *b) {
    return std::tuple{a->query_begin, a->reference_begin, a->id} <
           std::tuple{b->query_begin, b->reference_begin, b->id};
  });
  const auto rc = CountContexts(r, options), qc = CountContexts(q, options);
  SequenceRecord qrc;
  qrc.bases = CrossReverseComplement(q.bases);
  const auto qrc_counts = CountContexts(qrc, options);
  std::vector<const A *> strong;
  for (const auto *a : ordered) {
    check();
    if (DistinguishableContext(*a, r, rc, qc))
      strong.push_back(a);
  }
  std::vector<Range> accepted_r, accepted_q;
  std::vector<CentromereStitchGap> fixed_gaps;
  accepted_r.reserve(strong.size());
  accepted_q.reserve(strong.size());
  out.reviews.reserve(strong.size());
  for (std::size_t k = 1; k < strong.size(); ++k) {
    check();
    const auto &left = *strong[k - 1], &right = *strong[k];
    if (left.region_id == right.region_id ||
        left.query_begin + left.length >= right.query_begin ||
        left.reference_begin + left.length >= right.reference_begin)
      continue;
    CentromereEndpointReview review;
    review.id = out.reviews.size() + 1;
    review.region_id = left.region_id;
    review.left_anchor_id = left.id;
    review.right_anchor_id = right.id;
    review.reference_begin = left.reference_begin + left.length;
    review.reference_end = right.reference_begin;
    review.query_begin = left.query_begin + left.length;
    review.query_end = right.query_begin;
    const auto rb = review.reference_begin, re = review.reference_end;
    const auto qb = review.query_begin, qe = review.query_end;
    const Length rn = re - rb, qn = qe - qb;
    const auto begin = Clock::now();
    Length inner_count = 0;
    bool has_root = false;
    for (const auto &a : input.final) {
      if (a.strand == Strand::Forward &&
          (Touch(a.query_begin, a.query_begin + a.length, qb, qe) ||
           Touch(a.reference_begin, a.reference_begin + a.length, rb, re))) {
        ++inner_count;
        if (a.query_begin < qb || a.query_begin + a.length > qe ||
            a.reference_begin < rb || a.reference_begin + a.length > re)
          review.reason = "exterior_anchor_crosses_window";
      }
    }
    for (const auto &a : input.initial) {
      if (a.strand != Strand::Forward)
        continue;
      if (!Touch(a.query_begin, a.query_begin + a.length, qb, qe) &&
          !Touch(a.reference_begin, a.reference_begin + a.length, rb, re))
        continue;
      if (a.query_begin < qb || a.query_begin + a.length > qe ||
          a.reference_begin < rb || a.reference_begin + a.length > re)
        review.reason = "initial_anchor_crosses_window";
      else {
        has_root = true;
        if (!review.anchor_id)
          review.anchor_id = a.id;
      }
    }
    if (rn > 65536 || qn > 65536 || inner_count > 256)
      review.reason = "block_geometry_budget";
    if (review.reason.empty() && !has_root)
      review.reason = "no_root_in_block";
    if (review.reason.empty() && copies)
      for (const auto &a : copies->anchors)
        if (Touch(a.reference_begin, a.reference_begin + a.length, rb, re) ||
            Touch(a.query_begin, a.query_begin + a.length, qb, qe)) {
          review.reason = "accepted_copy_footprint";
          break;
        }
    if (review.reason.empty())
      for (const auto &segment : stitched.segments)
        if (segment.path_role == "copy" &&
            (Touch(segment.alignment.reference_begin,
                   segment.alignment.reference_end, rb, re) ||
             Touch(segment.alignment.query_begin,
                   segment.alignment.query_end, qb, qe))) {
          review.reason = "accepted_copy_segment_footprint";
          break;
        }
    if (review.reason.empty())
      for (const auto &a : input.final)
        if (a.strand == Strand::Reverse &&
            (Touch(a.query_begin, a.query_begin + a.length, qb, qe) ||
             Touch(a.reference_begin, a.reference_begin + a.length, rb, re)) &&
            DistinguishableContext(a, r, rc, qrc_counts)) {
          review.reason = "protected_reverse_context";
          break;
        }
    if (review.reason.empty())
      for (std::size_t i = 0; i < accepted_q.size(); ++i)
        if (Touch(rb, re, accepted_r[i].first, accepted_r[i].second) ||
            Touch(qb, qe, accepted_q[i].first, accepted_q[i].second)) {
          review.reason = "overlapping_accepted_window";
          break;
        }
    if (review.reason.empty() && !stitched.complete)
      review.reason = "unresolved_old_paths";
    if (review.reason.empty() &&
        (out.kernel_calls >= options.max_reviews ||
         Seconds(started) >= options.stage_seconds_limit))
      review.reason = "review_budget";
    if (!review.reason.empty()) {
      out.reviews.push_back(std::move(review));
      continue;
    }
    const double window_budget =
        std::min(options.window_seconds_limit,
                 options.stage_seconds_limit - Seconds(started));
    auto kernel = options.kernel;
    kernel.progress_callback = {};
    kernel.observation_callback = {};
    kernel.gap_seconds_limit =
        std::min(kernel.gap_seconds_limit, window_budget);
    kernel.stage_seconds_limit = window_budget;
    const auto caller = options.kernel.interruption_callback;
    kernel.interruption_callback = [&](std::string_view phase) {
      if (caller)
        caller(phase);
      if (Seconds(begin) > window_budget) {
        detail::StitchGapResourceLimit error("gap_time_budget");
        error.engine = "cross_region_review";
        throw error;
      }
    };
    try {
      kernel.interruption_callback("centromere-cross-region-mask");
      const auto old = CrossFragments(stitched, r, q, rb, re, qb, qe);
      review.old_score = old.score;
      review.old_exact_bases = old.exact;
      auto rv = r.bases.substr(rb, rn), qv = q.bases.substr(qb, qn);
      if (old.exact >= std::min(rn, qn)) {
        review.reason = "no_exact_headroom";
      } else if (rv == qv) {
        review.reason = "identical_window_ancestry_ambiguous";
      } else {
        std::vector<bool> allowed(static_cast<std::size_t>(qn), false);
        for (const auto &region : input.regions)
          if (region.strand == Strand::Forward) {
            const auto a = std::max(qb, region.query_begin);
            const auto b = std::min(qe, region.query_end);
            if (a < b)
              std::fill(allowed.begin() + static_cast<std::ptrdiff_t>(a - qb),
                        allowed.begin() + static_cast<std::ptrdiff_t>(b - qb),
                        true);
          }
        for (std::size_t i = 0; i < qv.size(); ++i)
          if (!allowed[i])
            qv[i] = 'N';
        // An opposite path keeps every original anchor and every reference
        // footprint, even if its sequence context is too short to veto review.
        for (const auto &a : input.final)
          if (a.strand == Strand::Reverse) {
            const auto x = std::max(rb, a.reference_begin);
            const auto y = std::min(re, a.reference_begin + a.length);
            if (x < y)
              std::fill(rv.begin() + static_cast<std::ptrdiff_t>(x - rb),
                        rv.begin() + static_cast<std::ptrdiff_t>(y - rb), 'N');
          }
        // The emitted negative span includes paired gap bases between its
        // frozen anchors; protecting anchors alone would miss those relations.
        for (const auto &segment : stitched.segments)
          if (segment.path_role == "primary" &&
              segment.alignment.strand == Strand::Reverse) {
            const auto x = std::max(rb, segment.alignment.reference_begin);
            const auto y = std::min(re, segment.alignment.reference_end);
            if (x < y)
              std::fill(rv.begin() + static_cast<std::ptrdiff_t>(x - rb),
                        rv.begin() + static_cast<std::ptrdiff_t>(y - rb), 'N');
          }
        kernel.interruption_callback("centromere-cross-region-dp");
        ++out.kernel_calls;
        auto proposal = detail::AlignCentromereBlockDpExactGaps(rv, qv, kernel);
        const auto raw_r = std::string_view(r.bases).substr(rb, rn);
        const auto raw_q = std::string_view(q.bases).substr(qb, qn);
        review.new_score = detail::ScoreCentromereStitchCigar(
            raw_r, raw_q, proposal.cigar, kernel.score_model);
        if (review.new_score != proposal.score)
          throw AlignmentError("cross-region candidate raw score differs");
        review.scored = true;
        review.engine = proposal.engine;
        review.peak_kernel_bytes = proposal.peak_kernel_bytes;
        review.dp_cells = proposal.dp_cells;
        review.wfa_states = proposal.wfa_states;
        review.new_cigar = std::move(proposal.cigar);
        std::vector<Replacement> replacements;
        for (const auto &region : input.regions)
          if (region.strand == Strand::Forward &&
              Touch(region.query_begin, region.query_end, qb, qe))
            replacements.push_back({region.id, rb, re,
                                    std::max(qb, region.query_begin),
                                    std::min(qe, region.query_end), {}});
        Position rp = rb, qp = qb;
        for (auto op : review.new_cigar) {
          if (op.operation != '=' && op.operation != 'I' && op.operation != 'D')
            throw AlignmentError("cross-region exact candidate contains mismatch");
          if (op.operation == '=') {
            for (auto &replacement : replacements) {
              const auto x = std::max(qp, replacement.qb);
              const auto y = std::min(qp + op.length, replacement.qe);
              if (x < y) {
                A a{0, replacement.region, 0, rp + x - qp, x, y - x,
                    Strand::Forward, static_cast<double>(y - x)};
                CheckAnchor(A{1, a.region_id, 1, a.reference_begin, a.query_begin,
                              a.length, a.strand, a.weight}, r, q);
                replacement.exact.push_back(a);
                review.new_exact_bases += a.length;
              }
            }
          }
          if (op.operation != 'I')
            rp += op.length;
          if (op.operation != 'D')
            qp += op.length;
        }
        if (rp != re || qp != qe)
          throw AlignmentError("cross-region candidate consumption differs");
        const auto old_contexts=CrossOldContexts(r,rc,qc,stitched,review);
        const auto new_contexts=CrossPathContexts(r,rc,qc,review.new_cigar,rb,qb,review);
        review.old_unique_context_starts=old_contexts.size();
        review.new_unique_context_starts=new_contexts.size();
        if (review.new_score <= review.old_score)
          review.reason = review.new_score == review.old_score
                              ? "equal_score_preserved" : "lower_score_preserved";
        else if (review.new_exact_bases <= review.old_exact_bases)
          review.reason = "no_exact_coverage_gain";
        else if(new_contexts.empty())
          review.reason="no_internal_unique_context";
        else if(!std::includes(new_contexts.begin(),new_contexts.end(),
                              old_contexts.begin(),old_contexts.end()))
          review.reason="existing_internal_context_not_preserved";
        else {
          auto candidate = out.anchors;
          Replace(candidate.initial, replacements, input.regions);
          Replace(candidate.final, replacements, input.regions);
          Refresh(candidate);
          auto candidate_fixed = fixed_gaps;
          auto added_fixed = CrossFixedGaps(r, q, candidate, review);
          candidate_fixed.insert(candidate_fixed.end(),
                                 std::make_move_iterator(added_fixed.begin()),
                                 std::make_move_iterator(added_fixed.end()));
          kernel.interruption_callback("centromere-cross-region-restitch");
          auto restitched = detail::AlignCentromereAnchoredPathsWithFixedGaps(
              r, q, candidate, copies, kernel, candidate_fixed);
          kernel.interruption_callback("centromere-cross-region-verify");
          if (!restitched.complete || restitched.aborted) {
            review.reason = "candidate_restitch_unresolved";
          } else {
            const auto actual = CrossFragments(restitched, r, q, rb, re, qb, qe);
            if (actual.mismatches)
              review.reason = "restitch_substitutions_preserved";
            else if (actual.score <= old.score)
              review.reason = "no_projected_score_gain";
            else if (actual.exact <= old.exact)
              review.reason = "no_projected_exact_gain";
            else {
              // Publish only the fully checked anchor copy; directions,
              // negative roots and accepted copy branches remain unchanged.
              review.reason = "strict_sum_existing_forward_fragments_gain";
              kernel.interruption_callback("centromere-cross-region-publish");
              accepted_r.emplace_back(rb, re);
              accepted_q.emplace_back(qb, qe);
              review.accepted = true;
              ++out.accepted;
              out.anchors = std::move(candidate);
              out.stitched = std::move(restitched);
              fixed_gaps = std::move(candidate_fixed);
            }
          }
        }
      }
    } catch (const detail::StitchGapResourceLimit &error) {
      review.reason = error.reason;
      review.engine = error.engine;
      review.peak_kernel_bytes = std::max(review.peak_kernel_bytes, error.peak_kernel_bytes);
      review.dp_cells = std::max(review.dp_cells, error.dp_cells);
      review.wfa_states = std::max(review.wfa_states, error.wfa_states);
    } catch (const std::bad_alloc &) {
      review.reason = "allocation_failure";
    }
    review.seconds = Seconds(begin);
    out.reviews.push_back(std::move(review));
  }
  check();
  out.seconds = Seconds(started);
  return out;
}
} // namespace ramag

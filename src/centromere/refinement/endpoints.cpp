#include "ramag/internal/centromere/refinement.hpp"
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
double Seconds(Clock::time_point s) {
  return std::chrono::duration<double>(Clock::now() - s).count();
}
std::vector<Range> Union(std::vector<Range> x) {
  std::sort(x.begin(), x.end());
  std::vector<Range> y;
  for (auto p : x)
    if (p.first < p.second) {
      if (y.empty() || y.back().second < p.first)
        y.push_back(p);
      else
        y.back().second = std::max(y.back().second, p.second);
    }
  return y;
}
Length Bases(const std::vector<Range> &x) {
  Length n = 0;
  for (auto [a, b] : x)
    n += b - a;
  return n;
}
bool Touch(Position a, Position b, Position c, Position d) {
  return a < d && c < b;
}
void Emit(std::vector<CigarOp> &v, char c, Length n) {
  if (!n)
    return;
  if (!v.empty() && v.back().operation == c)
    v.back().length += n;
  else
    v.push_back({c, n});
}
void Cigar(std::ostream &o, const std::vector<CigarOp> &c) {
  if (c.empty()) {
    o << '.';
    return;
  }
  for (auto x : c)
    o << x.length << x.operation;
}
void CheckAnchor(const A &a, const SequenceRecord &r, const SequenceRecord &q) {
  if ((a.strand != Strand::Forward && a.strand != Strand::Reverse) || !a.id ||
      !a.chain_index || !a.length || a.reference_begin > r.size() ||
      a.length > r.size() - a.reference_begin || a.query_begin > q.size() ||
      a.length > q.size() - a.query_begin)
    throw AlignmentError("endpoint input anchor bounds/identity");
  for (Length i = 0; i < a.length; ++i) {
    char x = r.bases[a.reference_begin + i];
    char y = a.strand == Strand::Forward
                 ? q.bases[a.query_begin + i]
                 : q.bases[a.query_begin + a.length - 1 - i];
    if (a.strand == Strand::Reverse) {
      switch (y) {
      case 'A':
        y = 'T';
        break;
      case 'C':
        y = 'G';
        break;
      case 'G':
        y = 'C';
        break;
      case 'T':
        y = 'A';
        break;
      default:
        y = 'N';
      }
    }
    if (x != y || !(x == 'A' || x == 'C' || x == 'G' || x == 'T'))
      throw AlignmentError("endpoint input anchor not exact");
  }
}
void Refresh(CentromereDirectedAnchorResult &out) {
  out.uncovered.clear();
  out.reference_overlaps.clear();
  std::map<std::uint64_t, std::vector<Range>> refs;
  for (auto &t : out.tasks) {
    std::vector<Range> ir, iq, fr, fq;
    t.initial_anchor_count = 0;
    t.final_anchor_count = 0;
    for (const auto &a : out.initial)
      if (a.region_id == t.region.id) {
        ++t.initial_anchor_count;
        ir.emplace_back(a.reference_begin, a.reference_begin + a.length);
        iq.emplace_back(a.query_begin, a.query_begin + a.length);
      }
    for (const auto &a : out.final)
      if (a.region_id == t.region.id) {
        ++t.final_anchor_count;
        fr.emplace_back(a.reference_begin, a.reference_begin + a.length);
        fq.emplace_back(a.query_begin, a.query_begin + a.length);
      }
    fr = Union(std::move(fr));
    fq = Union(std::move(fq));
    refs[t.region.id] = fr;
    t.initial_reference_bases = Bases(Union(std::move(ir)));
    t.initial_query_bases = Bases(Union(std::move(iq)));
    t.final_reference_bases = Bases(fr);
    t.final_query_bases = Bases(fq);
    if (t.final_query_bases < t.initial_query_bases)
      throw AlignmentError("endpoint reviewed initial coverage exceeds final");
    t.recursive_new_query_bases = t.final_query_bases - t.initial_query_bases;
    t.unanchored_query_bases = 0;
    t.unanchored_interval_count = 0;
    t.longest_unanchored_query_bases = 0;
    Position cursor = t.region.query_begin;
    auto add = [&](Position b, Position e) {
      if (b == e)
        return;
      out.uncovered.push_back({t.region.id, b, e, t.region.strand});
      t.unanchored_query_bases += e - b;
      ++t.unanchored_interval_count;
      t.longest_unanchored_query_bases =
          std::max(t.longest_unanchored_query_bases, e - b);
    };
    for (auto [b, e] : fq) {
      add(cursor, b);
      cursor = e;
    }
    add(cursor, t.region.query_end);
  }
  for (std::size_t i = 0; i < out.regions.size(); ++i)
    for (std::size_t j = i + 1; j < out.regions.size(); ++j) {
      auto a = out.regions[i].id, b = out.regions[j].id;
      std::size_t x = 0, y = 0;
      while (x < refs[a].size() && y < refs[b].size()) {
        auto lo = std::max(refs[a][x].first, refs[b][y].first),
             hi = std::min(refs[a][x].second, refs[b][y].second);
        if (lo < hi)
          out.reference_overlaps.push_back({a, b, lo, hi});
        if (refs[a][x].second < refs[b][y].second)
          ++x;
        else
          ++y;
      }
    }
  std::sort(out.uncovered.begin(), out.uncovered.end(),
            [](auto a, auto b) { return a.query_begin < b.query_begin; });
}
void Replace(std::vector<A> &list, const std::vector<Replacement> &replacements,
             const std::vector<CentromereDirectedRegion> &regions) {
  std::uint64_t next = 0;
  for (const auto &a : list)
    next = std::max(next, a.id);
  std::vector<A> kept;
  kept.reserve(list.size());
  for (const auto &a : list) {
    bool remove = false;
    for (const auto &p : replacements)
      if (a.region_id == p.region &&
          Touch(a.query_begin, a.query_begin + a.length, p.qb, p.qe)) {
        if (a.query_begin < p.qb || a.query_begin + a.length > p.qe ||
            a.reference_begin < p.rb || a.reference_begin + a.length > p.re)
          throw AlignmentError("endpoint replacement cuts a retained anchor");
        remove = true;
        break;
      }
    if (!remove)
      kept.push_back(a);
  }
  for (const auto &p : replacements)
    for (auto a : p.exact) {
      if (next == std::numeric_limits<std::uint64_t>::max())
        throw AlignmentError("endpoint anchor ID overflow");
      a.id = ++next;
      kept.push_back(a);
    }
  std::map<std::uint64_t, std::size_t> rank;
  for (std::size_t i = 0; i < regions.size(); ++i)
    rank[regions[i].id] = i;
  std::sort(kept.begin(), kept.end(), [&](const A &a, const A &b) {
    return std::tuple{rank.at(a.region_id), a.reference_begin} <
           std::tuple{rank.at(b.region_id), b.reference_begin};
  });
  std::map<std::uint64_t, std::uint64_t> index;
  for (auto &a : kept)
    a.chain_index = ++index[a.region_id];
  list = std::move(kept);
}
} // namespace centromere_refinement_detail

CentromereEndpointResult
RefineCentromereRootEndpoints(const SequenceRecord &r, const SequenceRecord &q,
                              const CentromereDirectedAnchorResult &input,
                              const CentromereCopyAnchorResult *copies,
                              const CentromereEndpointOptions &options) {
  if (!input.complete || (copies && !copies->complete))
    throw AlignmentError("endpoint review requires completed anchor stages");
  if (!options.max_reviews || !std::isfinite(options.stage_seconds_limit) ||
      options.stage_seconds_limit <= 0 ||
      !std::isfinite(options.window_seconds_limit) ||
      options.window_seconds_limit <= 0 ||
      !options.kernel.gap_memory_limit_bytes ||
      !std::isfinite(options.kernel.gap_seconds_limit) ||
      options.kernel.gap_seconds_limit <= 0)
    throw AlignmentError("invalid endpoint resource limits");
  if (options.kernel.score_model !=
      CentromereStitchScoreModel::Centrolign3Piece)
    throw AlignmentError("endpoint review requires centrolign-3piece");
  const auto started = Clock::now();
  auto check = [&] {
    if (options.kernel.interruption_callback)
      options.kernel.interruption_callback("centromere-endpoint-review");
  };
  check();
  std::map<std::uint64_t, std::vector<const A *>> paths;
  std::set<std::uint64_t> ids;
  for (const auto &a : input.final) {
    check();
    CheckAnchor(a, r, q);
    if (!ids.insert(a.id).second)
      throw AlignmentError("duplicate endpoint anchor ID");
    paths[a.region_id].push_back(&a);
  }
  for (auto &[id, path] : paths) {
    (void)id;
    std::sort(path.begin(), path.end(),
              [](auto a, auto b) { return a->chain_index < b->chain_index; });
    for (std::size_t i = 0; i < path.size(); ++i) {
      if (path[i]->chain_index != i + 1)
        throw AlignmentError("endpoint chain index");
      if (i) {
        auto a = path[i - 1], b = path[i];
        if (a->strand != b->strand ||
            a->reference_begin + a->length > b->reference_begin ||
            (a->strand == Strand::Forward
                 ? a->query_begin + a->length > b->query_begin
                 : b->query_begin + b->length > a->query_begin))
          throw AlignmentError("endpoint chain order");
      }
    }
  }
  // Geometry and root membership use increasing oriented query coordinates.
  // Reviews, footprints and published anchors retain physical FASTA coordinates.
  const auto oriented_begin = [&q](const A &a) -> Position {
    if (a.query_begin > q.size() || a.length > q.size() - a.query_begin)
      throw AlignmentError("endpoint initial anchor query bounds");
    return a.strand == Strand::Forward
               ? a.query_begin
               : q.size() - a.query_begin - a.length;
  };
  CentromereEndpointResult out;
  std::vector<Replacement> accepted;
  // Preserve the existing forward work under shared resource and footprint
  // limits. Newly enabled reverse work cannot preempt an old forward review.
  // Both passes retain the same clock, call counter and accepted-window set.
  std::vector<std::pair<std::uint64_t, const std::vector<const A *> *>> order;
  for (const auto strand : {Strand::Forward, Strand::Reverse})
    for (const auto &[id, path] : paths)
      if (!path.empty() && path.front()->strand == strand)
        order.emplace_back(id, &path);
  for (const auto &[id, path_pointer] : order)
    for (std::size_t i = 1; i + 1 < path_pointer->size(); ++i) {
      check();
      const auto &path = *path_pointer;
      const auto &l = *path[i - 1];
      const auto &a = *path[i];
      const auto &z = *path[i + 1];
      const bool reverse = a.strand == Strand::Reverse;
      const auto lu = oriented_begin(l), au = oriented_begin(a),
                 zu = oriented_begin(z);
      const Length lr = a.reference_begin - l.reference_begin - l.length,
                   lq = au - lu - l.length,
                   rr = z.reference_begin - a.reference_begin - a.length,
                   rq = zu - au - a.length;
      if (!((lr == 0 && lq > 0 && rr > 0 && rq == 0) ||
            (lq == 0 && lr > 0 && rq > 0 && rr == 0)))
        continue;
      const auto left = std::max(lr, lq), right = std::max(rr, rq);
      if (l.length <= std::max(left, right) ||
          z.length <= std::max(left, right))
        continue;
      bool root = false;
      for (const auto &b : input.initial) {
        if (b.region_id != id || b.strand != a.strand ||
            b.reference_begin < a.reference_begin)
          continue;
        const auto offset = b.reference_begin - a.reference_begin;
        if (offset > a.length || b.length > a.length - offset)
          continue;
        const auto bu = oriented_begin(b);
        if (bu >= au && bu - au == offset) {
          root = true;
          break;
        }
      }
      if (!root)
        continue;
      CentromereEndpointReview review;
      review.id = out.reviews.size() + 1;
      review.region_id = id;
      review.left_anchor_id = l.id;
      review.anchor_id = a.id;
      review.right_anchor_id = z.id;
      review.reference_begin = l.reference_begin + l.length;
      review.reference_end = z.reference_begin;
      review.query_begin = reverse ? z.query_begin + z.length
                                   : l.query_begin + l.length;
      review.query_end = reverse ? l.query_begin : z.query_begin;
      review.left_gap = left;
      review.right_gap = right;
      review.old_exact_bases = a.length;
      for (const auto &p : accepted)
        if (Touch(p.qb, p.qe, review.query_begin, review.query_end) ||
            Touch(p.rb, p.re, review.reference_begin, review.reference_end)) {
          review.reason = "overlapping_accepted_window";
          break;
        }
      if (review.reason.empty() && copies)
        for (const auto &b : copies->anchors)
          if (Touch(b.query_begin, b.query_begin + b.length, review.query_begin,
                    review.query_end) ||
              Touch(b.reference_begin, b.reference_begin + b.length,
                    review.reference_begin, review.reference_end)) {
            review.reason = "accepted_copy_footprint";
            break;
          }
      if (review.reason.empty() &&
          (out.kernel_calls >= options.max_reviews ||
           Seconds(started) >= options.stage_seconds_limit))
        review.reason = "review_budget";
      if (!review.reason.empty()) {
        out.reviews.push_back(std::move(review));
        continue;
      }
      const auto begin = Clock::now();
      const double window_budget =
          std::min({options.kernel.gap_seconds_limit, options.window_seconds_limit,
                    options.stage_seconds_limit - Seconds(started)});
      if (window_budget <= 0) {
        review.reason = "review_budget";
        out.reviews.push_back(std::move(review));
        continue;
      }
      Emit(review.old_cigar, lr ? 'D' : 'I', left);
      Emit(review.old_cigar, '=', a.length);
      Emit(review.old_cigar, rr ? 'D' : 'I', right);
      const auto rv = std::string_view(r.bases).substr(
          review.reference_begin, review.reference_end - review.reference_begin);
      auto qv = std::string_view(q.bases).substr(
          review.query_begin, review.query_end - review.query_begin);
      std::string oriented_query;
      std::uint64_t orientation_bytes = 0;
      if (reverse) {
        // Only the bounded window is materialized; its storage shares the
        // existing working-memory budget with the candidate kernel.
        if (qv.size() >= options.kernel.gap_memory_limit_bytes) {
          review.reason = "gap_memory_budget";
          review.seconds = Seconds(begin);
          out.reviews.push_back(std::move(review));
          continue;
        }
        try {
          oriented_query.resize(qv.size());
        } catch (const std::bad_alloc &) {
          review.reason = "gap_memory_budget";
          review.seconds = Seconds(begin);
          out.reviews.push_back(std::move(review));
          continue;
        }
        orientation_bytes = oriented_query.capacity() + 1;
        if (orientation_bytes >= options.kernel.gap_memory_limit_bytes) {
          review.reason = "gap_memory_budget";
          review.seconds = Seconds(begin);
          out.reviews.push_back(std::move(review));
          continue;
        }
        for (std::size_t j = 0; j < qv.size(); ++j) {
          if ((j % 4096) == 0) {
            check();
            if (Seconds(begin) >= window_budget) {
              review.reason = "gap_time_budget";
              break;
            }
          }
          switch (qv[qv.size() - 1 - j]) {
          case 'A': oriented_query[j] = 'T'; break;
          case 'C': oriented_query[j] = 'G'; break;
          case 'G': oriented_query[j] = 'C'; break;
          case 'T': oriented_query[j] = 'A'; break;
          case 'N': oriented_query[j] = 'N'; break;
          default: throw AlignmentError("endpoint query requires uppercase ACGTN");
          }
        }
        if (!review.reason.empty()) {
          review.seconds = Seconds(begin);
          out.reviews.push_back(std::move(review));
          continue;
        }
        qv = oriented_query;
      }
      review.old_score = detail::ScoreCentromereStitchCigar(
          rv, qv, review.old_cigar, options.kernel.score_model);
      // Equal strings cannot justify a different ancestral copy choice.
      if (rv == qv) {
        review.reason = "identical_window_ancestry_ambiguous";
        review.seconds = Seconds(begin);
        out.reviews.push_back(std::move(review));
        continue;
      }
      auto kernel = options.kernel;
      kernel.progress_callback = {};
      kernel.observation_callback = {};
      kernel.gap_memory_limit_bytes -= orientation_bytes;
      kernel.gap_seconds_limit = window_budget - Seconds(begin);
      if (kernel.gap_seconds_limit <= 0) {
        review.reason = "gap_time_budget";
        review.seconds = Seconds(begin);
        out.reviews.push_back(std::move(review));
        continue;
      }
      ++out.kernel_calls;
      try {
        // Both strands use the same exact, bounded seven-state kernel and
        // traceback policy. No repeat-length or sample-specific routing.
        auto proposal = detail::AlignCentromereBlockDp(rv, qv, kernel);
        review.new_score = detail::ScoreCentromereStitchCigar(
            rv, qv, proposal.cigar, kernel.score_model);
        if (review.new_score != proposal.score)
          throw AlignmentError("endpoint kernel score/replay differs");
        review.scored = true;
        review.engine = proposal.engine;
        review.peak_kernel_bytes = proposal.peak_kernel_bytes + orientation_bytes;
        review.dp_cells = proposal.dp_cells;
        review.wfa_states = proposal.wfa_states;
        review.new_cigar = std::move(proposal.cigar);
        Replacement p{id,
                      review.reference_begin,
                      review.reference_end,
                      review.query_begin,
                      review.query_end,
                      {}};
        Position rb = p.rb, consumed_query = 0;
        for (auto op : review.new_cigar) {
          if (op.operation == 'X')
            review.new_mismatch_bases += op.length;
          if (op.operation == '=') {
            const Position qb = reverse ? p.qe - consumed_query - op.length
                                        : p.qb + consumed_query;
            p.exact.push_back({0, id, 0, rb, qb, op.length, a.strand,
                               static_cast<double>(op.length)});
            review.new_exact_bases += op.length;
          }
          if (op.operation != 'I')
            rb += op.length;
          if (op.operation != 'D')
            consumed_query += op.length;
        }
        if (rb != p.re || consumed_query != p.qe - p.qb)
          throw AlignmentError("endpoint candidate consumption differs");
        // Endpoint-only refinement must not exchange exact/indel structure
        // for a substitution-rich alignment to another repeat.
        if (Seconds(begin) >= window_budget ||
            Seconds(started) >= options.stage_seconds_limit)
          review.reason = "gap_time_budget";
        else if (review.new_mismatch_bases)
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
      } catch (const detail::StitchGapResourceLimit &e) {
        review.reason = e.reason;
        review.engine = e.engine;
        review.peak_kernel_bytes = e.peak_kernel_bytes + orientation_bytes;
        review.dp_cells = e.dp_cells;
        review.wfa_states = e.wfa_states;
      }
      review.seconds = Seconds(begin);
      if (review.accepted && review.seconds > window_budget) {
        accepted.pop_back();
        --out.accepted;
        review.accepted = false;
        review.reason = "gap_time_budget";
      }
      out.reviews.push_back(std::move(review));
    }
  out.anchors = input;
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

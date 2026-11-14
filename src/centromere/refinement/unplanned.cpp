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
namespace {
constexpr Length UnplannedContextLength = 32;

std::vector<A> ReadUnplannedEvidence(const std::filesystem::path &path) {
  std::ifstream in(path);
  std::string line;
  if (!std::getline(in, line) ||
      line != "anchor_id\tblock_id\treference_begin\treference_end\tquery_begin\tquery_end\tstrand\tmatch_length\tweight")
    throw AlignmentError("invalid unplanned initial evidence header");
  std::vector<A> result;
  while (std::getline(in, line)) {
    A a;
    Position re{}, qe{};
    char sign{};
    std::string extra;
    std::istringstream row(line);
    if (!(row >> a.id >> a.region_id >> a.reference_begin >> re >> a.query_begin >> qe >>
          sign >> a.length >> a.weight) || row >> extra || !a.id || !a.length ||
        (sign != '+' && sign != '-') || re < a.reference_begin || qe < a.query_begin ||
        re - a.reference_begin != a.length || qe - a.query_begin != a.length)
      throw AlignmentError("invalid unplanned initial evidence row");
    a.chain_index = 1;
    a.strand = sign == '+' ? Strand::Forward : Strand::Reverse;
    if (result.size() == 100000) throw detail::StitchGapResourceLimit("candidate_budget");
    result.push_back(a);
  }
  if (in.bad()) throw AlignmentError("cannot read unplanned initial evidence");
  return result;
}

Position OrientedBegin(const A &a, Length query_length) {
  return a.strand == Strand::Forward ? a.query_begin : query_length - a.query_begin - a.length;
}

std::optional<A> LocateUnplannedCore(A a, const SequenceRecord &r, const SequenceRecord &q,
    std::string_view reverse_r, std::string_view reverse_q, const std::function<void()> &check) {
  const auto k = UnplannedContextLength;
  if (a.length < 4 * k) return std::nullopt;
  const auto rb = a.reference_begin + k, qb = OrientedBegin(a, q.size()) + k;
  const auto length = a.length - 2 * k;
  const std::string_view ref = r.bases, query = q.bases;
  const auto qview = a.strand == Strand::Forward ? query : reverse_q;
  const auto qother = a.strand == Strand::Forward ? reverse_q : query;
  const auto word = ref.substr(rb, length);
  auto unique = [](std::string_view sequence, std::string_view reverse, Position pos,
                   std::string_view pattern) {
    return sequence.find(pattern) == pos &&
        sequence.find(pattern, pos + 1) == std::string_view::npos &&
        reverse.find(pattern) == std::string_view::npos;
  };
  auto locates = [&](Length offset, Length n) {
    check();
    const auto pattern = word.substr(offset, n);
    return unique(ref, reverse_r, rb + offset, pattern) &&
           unique(qview, qother, qb + offset, pattern);
  };
  if (!locates(0, length)) return std::nullopt;
  auto shortest = [&](bool suffix) {
    Length lo = k, hi = length;
    while (lo < hi) {
      const auto mid = lo + (hi - lo) / 2;
      if (locates(suffix ? length - mid : 0, mid)) hi = mid;
      else lo = mid + 1;
    }
    return lo;
  };
  const auto left = shortest(false), right = shortest(true);
  // These words establish location. The publication ends inside both of
  // them, so a maximal exact match's ambiguous ends never become output.
  if (left + right + k > length) return std::nullopt;
  a.reference_begin = rb + left;
  a.length = length - left - right;
  a.query_begin = a.strand == Strand::Forward ? qb + left : q.size() - qb - left - a.length;
  a.weight = static_cast<double>(a.length);
  return a;
}

bool UnplannedCoresLink(const A &left, const A &right, Length query_length) {
  if (left.strand != right.strand) return false;
  const auto lq = OrientedBegin(left, query_length), rq = OrientedBegin(right, query_length);
  if (lq + left.length > rq || left.reference_begin + left.length > right.reference_begin)
    return false;
  const auto qgap = rq - lq - left.length;
  const auto rgap = right.reference_begin - left.reference_begin - left.length;
  if (std::max(qgap, rgap) > 65536) return false;
  const auto imbalance = std::max(qgap, rgap) - std::min(qgap, rgap);
  const auto cost = imbalance ? std::min({60 + 30 * imbalance, 800 + 5 * imbalance,
                                         2500 + imbalance}) : 0;
  return cost / 20 < std::min(left.length, right.length);
}
bool UnplannedCoresConflict(const A &a, const A &b) {
  return Touch(a.query_begin, a.query_begin + a.length, b.query_begin, b.query_begin + b.length) ||
      Touch(a.reference_begin, a.reference_begin + a.length, b.reference_begin, b.reference_begin + b.length);
}

void AddUnplannedRegions(CentromereDirectedAnchorResult &state, std::vector<A> &cores,
                         Length reference_length, Length query_length) {
  std::uint64_t region = 0, anchor = 0;
  for (const auto &v : state.regions) region = std::max(region, v.id);
  for (const auto *list : {&state.initial, &state.final})
    for (const auto &v : *list) anchor = std::max(anchor, v.id);
  for (auto &a : cores) {
    if (std::max(region, anchor) == std::numeric_limits<std::uint64_t>::max())
      throw AlignmentError("unplanned core identity overflow");
    a.region_id = ++region;
    a.id = ++anchor;
    CentromereDirectedRegion v{region, a.query_begin, a.query_begin + a.length, a.strand, 0};
    state.regions.push_back(v);
    state.initial.push_back(a);
    state.final.push_back(a);
    CentromereDirectedTask t;
    t.region = v;
    t.reference_length = reference_length;
    t.reference_begin = a.reference_begin;
    t.reference_end = a.reference_begin + a.length;
    t.status = "complete";
    state.tasks.push_back(t);
  }
  std::sort(state.regions.begin(), state.regions.end(), [](const auto &a, const auto &b) {
    return a.query_begin < b.query_begin;
  });
  state.unknown.clear();
  Position end = 0;
  for (const auto &v : state.regions) {
    if (end < v.query_begin) state.unknown.push_back({0, end, v.query_begin, Strand::Forward});
    end = v.query_end;
  }
  if (end < query_length) state.unknown.push_back({0, end, query_length, Strand::Forward});
  Refresh(state);
}
} // namespace

CentromereEndpointResult RecoverCentromereUnplannedCores(const SequenceRecord &r,
    const SequenceRecord &q, const CentromereDirectedAnchorResult &input,
    const CentromereCopyAnchorResult *copies, const CentromereStitchResult &stitched,
    const std::filesystem::path &evidence, const CentromereEndpointOptions &options) {
  CentromereEndpointResult out;
  if (!stitched.complete || stitched.aborted || !input.complete) return out;
  const auto started = Clock::now();
  const auto limit = options.kernel.gap_memory_limit_bytes;
  // Bound the working set before copying anchor/stitch state. Allow room
  // for interval indexes, temporary path validation and container growth.
  std::uint64_t records = input.initial.size() + input.final.size() + input.tasks.size() +
      input.trace.size() + stitched.gaps.size() + stitched.segments.size();
  for (const auto &p : stitched.segments)
    records += p.alignment.cigar.size() + p.anchor_ids.size() + p.gap_ids.size();
  for (const auto &g : stitched.gaps) records += g.cigar.size();
  const auto file_bytes = std::filesystem::file_size(evidence);
  auto report = [&](Position begin, Position end, std::string reason, double elapsed) {
    CentromereEndpointReview row;
    row.id = out.reviews.size() + 1;
    row.query_begin = begin;
    row.query_end = end;
    row.reason = std::move(reason);
    row.seconds = elapsed;
    out.reviews.push_back(std::move(row));
  };
  if (r.size() + q.size() > limit / 16 || records > limit / 4096 || file_bytes > limit / 16) {
    report(0, q.size(), "unplanned_memory_budget", Seconds(started));
    return out;
  }
  auto check_stage = [&] {
    if (options.kernel.interruption_callback)
      options.kernel.interruption_callback("centromere-unplanned-cores");
    if (Seconds(started) >= options.stage_seconds_limit)
      throw detail::StitchGapResourceLimit("unplanned_stage_budget");
  };
  try {
    const auto raw = ReadUnplannedEvidence(evidence);
    if (raw.empty() && !options.unplanned_candidates) return out;
    const auto reverse_r = CrossReverseComplement(r.bases), reverse_q = CrossReverseComplement(q.bases);
    std::vector<Range> ref_used, query_used;
    for (const auto &p : stitched.segments) {
      check_stage();
      const auto &a = p.alignment;
      Position rp = a.reference_begin;
      Position qp = a.strand == Strand::Forward ? a.query_begin : a.query_end;
      for (const auto &op : a.cigar) {
        const bool forward = a.strand == Strand::Forward;
        if (op.operation == '=' || op.operation == 'X') {
          ref_used.emplace_back(rp, rp + op.length);
          query_used.emplace_back(forward ? qp : qp - op.length, forward ? qp + op.length : qp);
        }
        if (op.operation != 'I') rp += op.length;
        if (op.operation != 'D') qp = forward ? qp + op.length : qp - op.length;
      }
    }
    ref_used = Union(std::move(ref_used));
    query_used = Union(std::move(query_used));
    auto overlaps = [](const std::vector<Range> &ranges, Position begin, Position end) {
      const auto it = std::upper_bound(ranges.begin(), ranges.end(), begin,
          [](auto pos, auto range) { return pos < range.second; });
      return it != ranges.end() && it->first < end;
    };
    struct Window {
      Position begin{}, end{};
      double seconds{};
      std::vector<A> cores;
      bool local{};
    };
    std::vector<Window> windows;
    std::size_t attempted = 0, proposals = 0;
    for (const auto &u : input.unknown) {
      const auto begin = u.query_begin, end = u.query_end;
      if (end - begin < 5000 || overlaps(query_used, begin, end)) continue;
      if (++attempted > options.max_reviews) {
        report(begin, end, "unplanned_review_budget", 0);
        break;
      }
      const auto window_started = Clock::now();
      auto check = [&] {
        check_stage();
        if (Seconds(window_started) >= options.window_seconds_limit)
          throw detail::StitchGapResourceLimit("unplanned_window_budget");
      };
      try {
        std::vector<A> candidates;
        std::set<std::tuple<Strand, Position, Position, Length>> seen;
        std::size_t visits = 0;
        for (auto a : raw) {
          check();
          const auto lo = std::max(begin, a.query_begin);
          // Validate bounds before adding untrusted unsigned coordinates.
          if (a.query_begin > q.size() || a.length > q.size() - a.query_begin ||
              a.reference_begin > r.size() || a.length > r.size() - a.reference_begin)
            throw AlignmentError("unplanned initial evidence out of bounds");
          const auto hi = std::min(end, a.query_begin + a.length);
          if (lo >= hi || hi - lo < 4 * UnplannedContextLength) continue;
          if (++visits > 3000) throw detail::StitchGapResourceLimit("candidate_budget");
          a.reference_begin += a.strand == Strand::Forward ? lo - a.query_begin : a.query_begin + a.length - hi;
          a.query_begin = lo;
          a.length = hi - lo;
          CheckAnchor(a, r, q);
          const auto core = LocateUnplannedCore(a, r, q, reverse_r, reverse_q, check);
          if (core && !overlaps(ref_used, core->reference_begin, core->reference_begin + core->length) &&
              seen.emplace(core->strand, core->query_begin, core->reference_begin, core->length).second)
            candidates.push_back(*core);
        }
        Window w{begin, end, 0, {}};
        for (const auto &a : candidates) {
          const bool linked = std::any_of(candidates.begin(), candidates.end(), [&](const auto &b) {
            check();
            return UnplannedCoresLink(a, b, q.size()) || UnplannedCoresLink(b, a, q.size());
          });
          if (linked) w.cores.push_back(a);
        }
        if (w.cores.empty() && options.unplanned_candidates) {
          const auto found = options.unplanned_candidates(begin, end, check);
          if (!found) throw detail::StitchGapResourceLimit("unplanned_local_search_budget");
          w.local = true;
          seen.clear();
          for (const auto &hit : *found) {
            check();
            if (hit.query_begin > end - begin || hit.length > end - begin - hit.query_begin)
              throw AlignmentError("unplanned local hit outside its window");
            A a{1, 0, 1, hit.reference_begin, begin + hit.query_begin, hit.length, hit.strand, 0};
            CheckAnchor(a, r, q);
            const auto core = LocateUnplannedCore(a, r, q, reverse_r, reverse_q, check);
            if (core && !overlaps(ref_used, core->reference_begin, core->reference_begin + core->length) &&
                seen.emplace(core->strand, core->query_begin, core->reference_begin, core->length).second)
              w.cores.push_back(*core);
          }
        }
        check();
        w.seconds = Seconds(window_started);
        if (proposals + w.cores.size() > 3000)
          throw detail::StitchGapResourceLimit("candidate_budget");
        proposals += w.cores.size();
        if (w.cores.empty()) report(begin, end, "unplanned_no_independent_chain", w.seconds);
        else windows.push_back(std::move(w));
      } catch (const detail::StitchGapResourceLimit &error) {
        report(begin, end, error.reason, Seconds(window_started));
        if (Seconds(started) >= options.stage_seconds_limit) break;
      }
    }
    // Competing independently localized hypotheses are withheld on either
    // axis. Disjoint forward and reverse subchains can coexist in one window.
    for (const auto &w : windows) {
      const auto window_started = Clock::now();
      auto check = [&] {
        check_stage();
        if (w.seconds + Seconds(window_started) >= options.window_seconds_limit)
          throw detail::StitchGapResourceLimit("unplanned_window_budget");
      };
      try {
        std::vector<A> kept;
        for (const auto &a : w.cores) {
          bool conflict = false;
          for (const auto &other : windows)
            if ((w.local || !other.local) && std::any_of(other.cores.begin(), other.cores.end(), [&](const auto &b) {
                  check();
                  return &a != &b && UnplannedCoresConflict(a, b);
                })) { conflict = true; break; }
          if (!conflict) kept.push_back(a);
        }
        std::vector<A> chosen;
        for (const auto &a : kept)
          // A discovered core already has two disjoint, globally locating
          // words. It need not bridge a second candidate across an unknown gap.
          if (w.local || std::any_of(kept.begin(), kept.end(), [&](const auto &b) {
                check();
                return UnplannedCoresLink(a, b, q.size()) || UnplannedCoresLink(b, a, q.size());
              })) chosen.push_back(a);
        if (chosen.empty()) {
          report(w.begin, w.end, "unplanned_competing_placement", w.seconds + Seconds(window_started));
          continue;
        }
        auto candidate = out.accepted ? out.anchors : input;
        auto alignment = out.stitched ? *out.stitched : stitched;
        const auto evidence_cores = chosen;
        AddUnplannedRegions(candidate, chosen, r.size(), q.size());
        detail::AppendCentromereExactSingletons(r, q, candidate, copies, chosen, alignment);
        check();
        out.anchors = std::move(candidate);
        out.stitched = std::move(alignment);
        out.accepted += chosen.size();
        const auto elapsed = w.seconds + Seconds(window_started);
        for (std::size_t i = 0; i < chosen.size(); ++i) {
          const auto &a = chosen[i];
          report(a.query_begin, a.query_begin + a.length,
                 w.local ? "unplanned_local_core_gain" : "unplanned_exact_core_gain", elapsed);
          auto &row = out.reviews.back();
          row.anchor_id = evidence_cores[i].id;
          row.region_id = a.region_id;
          row.reference_begin = a.reference_begin;
          row.reference_end = a.reference_begin + a.length;
          row.new_exact_bases = a.length;
          row.accepted = true;
          row.engine = w.local ? (a.strand == Strand::Forward ? "local_core_forward" : "local_core_reverse")
                               : (a.strand == Strand::Forward ? "initial_core_forward" : "initial_core_reverse");
          row.new_cigar = {{'=', a.length}};
        }
      } catch (const detail::StitchGapResourceLimit &error) {
        report(w.begin, w.end, error.reason, w.seconds + Seconds(window_started));
        if (Seconds(started) >= options.stage_seconds_limit) break;
      }
    }
  } catch (const detail::StitchGapResourceLimit &error) {
    report(0, q.size(), error.reason, Seconds(started));
  }
  out.seconds = Seconds(started);
  return out;
}
} // namespace ramag

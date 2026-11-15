#include "ramag/internal/centromere/support.hpp"
#include "ramag/centromere/anchors.hpp"
#include <algorithm>
#include <chrono>
#include <climits>
#include <cmath>
#include <iomanip>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <ostream>
#include <set>
#include <tuple>
#if RAMAG_HAVE_CENTROMERE
#include "anchor.h"
#include "ramag/internal/centromere/copy_cache.hpp"
#endif

namespace ramag {
namespace {
char Sign(Strand s) { return s == Strand::Forward ? '+' : '-'; }
std::string_view EvidenceRoute(const CentromereCopyAnchor &a) {
  if (a.context_evidence == "opposite_gap_competitor" ||
      a.context_evidence == "equivalent_supported_chain_competitor")
    return a.context_evidence;
  if (a.primary_anchor_id && a.self_repeat_id)
    return "primary_self_projection";
  return "unclassified";
}
#if RAMAG_HAVE_CENTROMERE
using Clock = std::chrono::steady_clock;
using Range = std::pair<Position, Position>;
using Key = std::tuple<Position, Position, Length, Strand>;
using Directed = CentromereDirectedAnchor;
struct Budget {};
char Complement(char c) {
  switch (c) {
  case 'A':
    return 'T';
  case 'C':
    return 'G';
  case 'G':
    return 'C';
  case 'T':
    return 'A';
  default:
    return 'N';
  }
}
using ::ramag::detail::IsCanonicalUpper;
std::string Reverse(std::string s) {
  std::reverse(s.begin(), s.end());
  for (char &c : s)
    c = Complement(c);
  return s;
}
void Check(const CentromereAlignmentOptions &o) {
  if (o.interruption_callback)
    o.interruption_callback("centromere-copy-anchors");
}
double Seconds(Clock::time_point t) {
  return std::chrono::duration<double>(Clock::now() - t).count();
}
Key KeyOf(const CentromereCopyAnchor &a) {
  return {a.reference_begin, a.query_begin, a.length, a.strand};
}
Position Oriented(const CentromereCopyAnchor &a, Length q) {
  return a.strand == Strand::Forward ? a.query_begin
                                     : q - a.query_begin - a.length;
}
bool Exact(const SequenceRecord &r, const SequenceRecord &q,
           const CentromereCopyAnchor &a) {
  if (!a.length || a.reference_begin > r.size() ||
      a.length > r.size() - a.reference_begin || a.query_begin > q.size() ||
      a.length > q.size() - a.query_begin)
    return false;
  for (Length i = 0; i < a.length; ++i) {
    char x = r.bases[a.reference_begin + i],
         y = a.strand == Strand::Forward
                 ? q.bases[a.query_begin + i]
                 : Complement(q.bases[a.query_begin + a.length - i - 1]);
    if (!IsCanonicalUpper(x) || x != y)
      return false;
  }
  return true;
}
std::vector<Range> Union(std::vector<Range> in) {
  std::sort(in.begin(), in.end());
  std::vector<Range> out;
  for (auto p : in) {
    if (out.empty() || out.back().second < p.first)
      out.push_back(p);
    else
      out.back().second = std::max(out.back().second, p.second);
  }
  return out;
}
Length Covered(const std::vector<CentromereCopyAnchor> &as, bool query) {
  std::vector<Range> rs;
  for (const auto &a : as) {
    auto b = query ? a.query_begin : a.reference_begin;
    rs.emplace_back(b, b + a.length);
  }
  Length n = 0;
  for (auto [b, e] : Union(std::move(rs)))
    n += e - b;
  return n;
}
std::vector<Range> Gaps(const std::vector<Directed> &as, bool query, Length n) {
  std::vector<Range> rs;
  for (const auto &a : as) {
    auto b = query ? a.query_begin : a.reference_begin;
    rs.emplace_back(b, b + a.length);
  }
  std::vector<Range> gaps;
  Position end = 0;
  for (auto [b, e] : Union(std::move(rs))) {
    if (end < b)
      gaps.emplace_back(end, b);
    end = e;
  }
  if (end < n)
    gaps.emplace_back(end, n);
  return gaps;
}
struct Job {
  CentromereCopyTask task;
  bool query{}, scout{};
  const Directed *left{};
  const Directed *right{};
};
std::vector<Job> Jobs(const SequenceRecord &r, const SequenceRecord &q,
                      const std::vector<Directed> &as) {
  std::vector<Job> out;
  for (bool query : {true, false}) {
    const auto n = query ? q.size() : r.size();
    for (auto [b, e] : Gaps(as, query, n)) {
      const Directed *left = nullptr;
      const Directed *right = nullptr;
      for (const auto &a : as) {
        const auto x = query ? a.query_begin : a.reference_begin;
        if (x + a.length <= b &&
            (!left || x + a.length >
                          (query ? left->query_begin : left->reference_begin) +
                              left->length))
          left = &a;
        if (x >= e && (!right || x < (query ? right->query_begin
                                            : right->reference_begin)))
          right = &a;
      }
      const auto span = e - b;
      Position sb = b > span ? b - span : 0, se = e + std::min(span, n - e);
      if (left)
        sb = std::min(sb, query ? left->query_begin : left->reference_begin);
      if (right)
        se =
            std::max(se, (query ? right->query_begin : right->reference_begin) +
                             right->length);
      CentromereCopyTask t;
      t.side = query ? "query" : "reference";
      t.source_begin = sb;
      t.source_end = se;
      t.target_begin = b;
      t.target_end = e;
      t.unresolved_bases = span;
      out.push_back({t, query, false, left, right});
    }
  }
  std::stable_sort(out.begin(), out.end(), [](const Job &a, const Job &b) {
    const auto x = a.task.target_end - a.task.target_begin,
               y = b.task.target_end - b.task.target_begin;
    if (x != y)
      return x > y;
    return std::tuple(a.query, a.task.target_begin) <
           std::tuple(b.query, b.task.target_begin);
  });
  // Covered-position scouts retain alternatives as ambiguity, never as new
  // accepted history. This also diagnoses identical ancestral copies when the
  // primary happened to cover the entire input in one exact anchor.
  for (bool query : {true, false}) {
    CentromereCopyTask t;
    t.side = query ? "query" : "reference";
    t.source_end = t.target_end = query ? q.size() : r.size();
    out.push_back({t, query, true, nullptr, nullptr});
  }
  for (std::size_t i = 0; i < out.size(); ++i)
    out[i].task.id = i + 1;
  return out;
}
std::uint64_t Count(std::string_view s, std::string_view pattern) {
  std::uint64_t n = 0;
  for (std::size_t p = 0; p <= s.size();) {
    p = s.find(pattern, p);
    if (p == std::string_view::npos)
      break;
    ++n;
    ++p;
  }
  return n;
}
struct Path {
  std::vector<std::size_t> ids;
  double score{};
};
// Copy-specific dynamic program: negative edge gain is never raised to 0.1.
Path Chain(const std::vector<CentromereCopyAnchor> &a, Length q,
           bool query_side, const std::set<std::size_t> &banned = {}) {
  std::vector<double> score(a.size(), -std::numeric_limits<double>::infinity());
  std::vector<std::size_t> prev(a.size(), a.size());
  std::size_t best = a.size();
  for (std::size_t i = 0; i < a.size(); ++i) {
    if (banned.contains(i))
      continue;
    score[i] = a[i].weight;
    const auto qi = Oriented(a[i], q);
    for (std::size_t j = 0; j < i; ++j) {
      if (banned.contains(j) || a[i].strand != a[j].strand ||
          a[i].self_source_traversal != a[j].self_source_traversal)
        continue;
      const auto qj = Oriented(a[j], q);
      if (a[j].reference_begin + a[j].length > a[i].reference_begin ||
          qj + a[j].length > qi ||
          a[j].self_source_oriented_begin + a[j].length >
              a[i].self_source_oriented_begin)
        continue;
      const auto dr = a[i].reference_begin - a[j].reference_begin - a[j].length,
                 dq = qi - qj - a[j].length;
      const auto delta = dr > dq ? dr - dq : dq - dr;
      const auto ds = a[i].self_source_oriented_begin -
                      a[j].self_source_oriented_begin - a[j].length,
                 dt = query_side ? dq : dr;
      const auto self_delta = ds > dt ? ds - dt : dt - ds;
      const double value =
          score[j] + a[i].weight -
          2.0 * std::log2(static_cast<double>(delta) + 1.0) -
          2.0 * std::log2(static_cast<double>(self_delta) + 1.0);
      if (value > score[i]) {
        score[i] = value;
        prev[i] = j;
      }
    }
    if (best == a.size() || score[i] > score[best])
      best = i;
  }
  Path out;
  if (best == a.size())
    return out;
  out.score = score[best];
  for (auto i = best; i != a.size(); i = prev[i])
    out.ids.push_back(i);
  std::reverse(out.ids.begin(), out.ids.end());
  return out;
}
CentromereCopyAnchor ClipTarget(const CentromereCopyAnchor &original,
                                Position begin, Position end, bool query_side) {
  auto x = original;
  const auto b = query_side ? original.query_begin : original.reference_begin;
  const auto lo = std::max(begin, b), hi = std::min(end, b + original.length);
  if (lo >= hi) {
    x.length = 0;
    return x;
  }
  x.length = hi - lo;
  const auto offset = lo - b;
  if (query_side) {
    x.query_begin = lo;
    x.reference_begin =
        original.reference_begin + (x.strand == Strand::Forward
                                        ? offset
                                        : original.length - offset - x.length);
  } else {
    x.reference_begin = lo;
    x.query_begin =
        original.query_begin + (x.strand == Strand::Forward
                                    ? offset
                                    : original.length - offset - x.length);
  }
  x.self_source_oriented_begin = original.self_source_oriented_begin +
                                 x.reference_begin - original.reference_begin;
  return x;
}
Position Counterpart(const CentromereCopyAnchor &a, Position target,
                     bool query_side) {
  const auto begin = query_side ? a.query_begin : a.reference_begin;
  const auto other = query_side ? a.reference_begin : a.query_begin;
  return other + (a.strand == Strand::Forward
                      ? target - begin
                      : a.length - 1 - (target - begin));
}
bool SameLine(const CentromereCopyAnchor &a, const CentromereCopyAnchor &b,
              Length q) {
  return a.strand == b.strand &&
         static_cast<long double>(a.reference_begin) - Oriented(a, q) ==
             static_cast<long double>(b.reference_begin) - Oriented(b, q);
}
std::vector<CentromereCopyAnchor> Merge(std::vector<CentromereCopyAnchor> a,
                                        Length q) {
  std::sort(a.begin(), a.end(), [&](const auto &x, const auto &y) {
    return std::tuple(x.strand,
                      static_cast<long double>(x.reference_begin) -
                          Oriented(x, q),
                      x.reference_begin, x.length) <
           std::tuple(y.strand,
                      static_cast<long double>(y.reference_begin) -
                          Oriented(y, q),
                      y.reference_begin, y.length);
  });
  std::vector<CentromereCopyAnchor> out;
  for (auto x : a) {
    if (!out.empty() && SameLine(out.back(), x, q) &&
        out.back().primary_anchor_id == x.primary_anchor_id &&
        x.reference_begin <= out.back().reference_begin + out.back().length) {
      auto &previous = out.back();
      if (x.reference_begin + x.length <=
          previous.reference_begin + previous.length)
        continue;
      if (x.reference_begin == previous.reference_begin) {
        previous = x;
        continue;
      }
    }
    out.push_back(x);
  }
  std::sort(out.begin(), out.end(), [](const auto &a0, const auto &b0) {
    return std::tuple(a0.reference_begin, a0.query_begin, a0.strand,
                      a0.length) <
           std::tuple(b0.reference_begin, b0.query_begin, b0.strand, b0.length);
  });
  return out;
}
// Enumerate the actual bounded rectangle, including mismatches between seeds.
// A caller must supply source proof before these candidates can be accepted.
void VisitCross(
    const SequenceRecord &r, const SequenceRecord &q, Position rb, Position re,
    Position qb, Position qe, Strand strand,
    const CentromereCopyAnchorOptions &co,
    const CentromereAlignmentOptions &options,
    const std::function<void(CentromereCopyAnchor, std::uint64_t)> &sink,
    bool supplemental = true) {
  if (rb >= re || qb >= qe)
    return;
  std::string rs = r.bases.substr(rb, re - rb),
              qs = q.bases.substr(qb, qe - qb);
  if (strand == Strand::Reverse)
    qs = Reverse(std::move(qs));
  std::vector<rama_source::SequenceInfo> data{
      {rs, "copy-window-reference", static_cast<uint_t>(rs.size())},
      {qs, "copy-window-query", static_cast<uint_t>(qs.size())}};
  rama_source::AnchorFinder finder(
      data, "", options.worker_threads, false, false, co.max_rare_level,
      [&] { Check(options); }, true);
  finder.visitRegionCandidates(
      {0, static_cast<uint_t>(rs.size()), 0, static_cast<uint_t>(qs.size())},
      supplemental, [&](const rama_source::RareCandidate &c) {
        const auto local = c.pair.second_pos - rs.size() - 1;
        CentromereCopyAnchor x;
        x.reference_begin = rb + c.pair.first_pos;
        x.length = c.pair.match_length;
        x.query_begin =
            strand == Strand::Forward ? qb + local : qe - local - x.length;
        x.strand = strand;
        if (Exact(r, q, x))
          sink(x, c.level);
      });
}
Range OppositeGap(const Job &job) {
  if (!job.left || !job.right || job.left->strand != job.right->strand)
    return {};
  const auto &l = *job.left;
  const auto &r = *job.right;
  const auto lb = job.query ? l.reference_begin : l.query_begin;
  const auto rb = job.query ? r.reference_begin : r.query_begin;
  return l.strand == Strand::Forward ? Range{lb + l.length, rb}
                                     : Range{rb + r.length, lb};
}
bool Flanks(const Job &job, const SequenceRecord &r, const SequenceRecord &q,
            std::map<Key, bool> &cache) {
  if (!job.left || !job.right || job.left == job.right ||
      job.left->strand != job.right->strand)
    return false;
  const auto &l = *job.left;
  const auto &h = *job.right;
  if (job.query) {
    if (l.strand == Strand::Forward
            ? l.reference_begin + l.length > h.reference_begin
            : h.reference_begin + h.length > l.reference_begin)
      return false;
  } else {
    if (l.strand == Strand::Forward ? l.query_begin + l.length > h.query_begin
                                    : h.query_begin + h.length > l.query_begin)
      return false;
  }
  auto unique = [&](const Directed &a) {
    Key k{a.reference_begin, a.query_begin, a.length, a.strand};
    auto it = cache.find(k);
    if (it != cache.end())
      return it->second;
    std::string pattern = r.bases.substr(a.reference_begin, a.length);
    const auto other = a.strand == Strand::Forward ? pattern : Reverse(pattern);
    bool ok = Count(r.bases, pattern) == 1 && Count(q.bases, other) == 1;
    cache.emplace(k, ok);
    return ok;
  };
  return unique(l) && unique(h);
}
#endif
} // namespace

CentromereCopyAnchorResult
FindCentromereCopyAnchors(const SequenceRecord &reference,
                          const SequenceRecord &query,
                          const CentromereDirectedAnchorResult &primary,
                          const CentromereCopyAnchorOptions &co,
                          const CentromereAlignmentOptions &options) {
#if RAMAG_HAVE_CENTROMERE
  if (!options.worker_threads || !co.max_rare_level ||
      !co.local_candidate_limit || !co.total_candidate_limit)
    throw CentromereInputError(
        "copy anchor threads, rarity and budgets must be positive");
  if (!primary.complete)
    throw CentromereInputError("copy anchors require a complete fixed primary");
  if (reference.bases.empty() || query.bases.empty())
    throw CentromereInputError("copy anchor sequences must be nonempty");
  constexpr auto safe_length = static_cast<Length>((INT_MAX - 64) / 3);
  if (reference.size() > safe_length / 2 || query.size() > safe_length / 2)
    throw CentromereInputError(
        "copy self-search exceeds pinned source coordinate range");
  for (const auto *s : {&reference, &query})
    for (char c : s->bases)
      if (!IsCanonicalUpper(c) && c != 'N')
        throw CentromereInputError("copy anchors require normalized A/C/G/T/N");
  for (const auto &a : primary.final) {
    if (a.strand != Strand::Forward && a.strand != Strand::Reverse)
      throw CentromereInputError("copy primary has invalid strand");
    CentromereCopyAnchor x;
    x.reference_begin = a.reference_begin;
    x.query_begin = a.query_begin;
    x.length = a.length;
    x.strand = a.strand;
    if (!Exact(reference, query, x))
      throw CentromereInputError("copy primary has invalid or inexact anchors");
  }
  Check(options);
  const auto start = Clock::now();
  CentromereCopyAnchorResult result;
  std::vector<CentromereCopyAmbiguity> pending_runner_witnesses;
  std::uint64_t branch_id = 0, group_id = 0;
  std::map<Key, bool> unique_cache;
  std::map<std::tuple<Position, Length, Strand>,
           std::pair<std::uint64_t, std::uint64_t>>
      counts;
  detail::CopyCountCache<> count_cache(reference.bases);
  auto weight = [&](CentromereCopyAnchor &a) {
    const auto key = std::tuple(a.reference_begin, a.length, a.strand);
    auto it = counts.find(key);
    if (it == counts.end()) {
      const auto reference_count =
          count_cache.Get(detail::CopyCountTarget::Reference, a.reference_begin,
                          a.length, false, [&](std::string_view pattern) {
                            return Count(reference.bases, pattern);
                          });
      const auto query_count = count_cache.Get(
          detail::CopyCountTarget::Query, a.reference_begin, a.length,
          a.strand == Strand::Reverse, [&](std::string_view pattern) {
            return Count(query.bases, pattern);
          });
      it = counts.emplace(key, std::pair(reference_count, query_count)).first;
    }
    a.weight = static_cast<double>(a.length) /
               std::sqrt(static_cast<double>(it->second.first) *
                         static_cast<double>(it->second.second));
  };
  std::set<Key> emitted;
  std::set<Key> global_projected_keys;
  std::set<Position> query_starts, query_ends, reference_starts, reference_ends;
  for (const auto &a : primary.final) {
    query_starts.insert(a.query_begin);
    query_ends.insert(a.query_begin + a.length);
    reference_starts.insert(a.reference_begin);
    reference_ends.insert(a.reference_begin + a.length);
  }
  std::set<std::tuple<std::string, Position, Position, Length, Strand>>
      global_self_keys;
  for (auto job : Jobs(reference, query, primary.final)) {
    Check(options);
    auto &t = job.task;
    const auto task_start = Clock::now();
    std::vector<CentromereCopyAnchor> projected;
    std::set<Key> projected_keys;
    std::set<std::tuple<Position, Position, Length, Strand>> self_keys;
    std::vector<CentromereCopySelfRepeat> repeats;
    bool exhausted = false;
    auto add_budget = [&](bool global_new = true) {
      if (t.candidate_count >= co.local_candidate_limit ||
          (global_new && result.total_candidates >= co.total_candidate_limit))
        throw Budget{};
      ++t.candidate_count;
      if (global_new)
        ++result.total_candidates;
    };
    const auto &seq = job.query ? query : reference;
    if (result.total_candidates >= co.total_candidate_limit)
      exhausted = true;
    const auto self_start = Clock::now();
    try {
      if (exhausted)
        throw Budget{};
      for (auto strand : {Strand::Forward, Strand::Reverse}) {
        std::string source = seq.bases.substr(t.source_begin,
                                              t.source_end - t.source_begin),
                    target = seq.bases.substr(t.target_begin,
                                              t.target_end - t.target_begin);
        if (strand == Strand::Reverse)
          target = Reverse(std::move(target));
        std::vector<rama_source::SequenceInfo> data{
            {source, "copy-source", static_cast<uint_t>(source.size())},
            {target, "copy-target", static_cast<uint_t>(target.size())}};
        rama_source::AnchorFinder finder(
            data, "", options.worker_threads, false, false, co.max_rare_level,
            [&] { Check(options); }, true);
        finder.visitRegionCandidates(
            {0, static_cast<uint_t>(source.size()), 0,
             static_cast<uint_t>(target.size())},
            true, [&](const rama_source::RareCandidate &c) {
              const Position s = t.source_begin + c.pair.first_pos;
              const Position local = c.pair.second_pos - source.size() - 1;
              Length len = c.pair.match_length;
              if (local > target.size() || len > target.size() - local)
                return;
              const Position d = strand == Strand::Forward
                                     ? t.target_begin + local
                                     : t.target_end - local - len;
              if (job.scout && strand == Strand::Forward && s != d &&
                  len >= std::max(s, d) - std::min(s, d))
                len = std::max(s, d) - std::min(s, d);
              if (job.scout && s + len != d && d + len != s)
                return;
              if (s == d || (s < d + len && d < s + len))
                return; // physical diagonal and overlapping periodic self
                        // matches
              if (!job.scout && s < t.target_end && t.target_begin < s + len)
                return;
              // A local copy seed must be attached to observed task structure.
              // Interior accidental short matches do not create new tasks
              // solely by equality.
              if (!job.scout) {
                const auto &starts =
                    job.query ? query_starts : reference_starts;
                const auto &ends = job.query ? query_ends : reference_ends;
                const bool boundary =
                    d == t.target_begin || d + len == t.target_end ||
                    starts.contains(s) || ends.contains(s + len);
                if (!boundary)
                  return;
              }

              const auto key =
                  std::tuple(std::min(s, d), std::max(s, d), len, strand);
              if (self_keys.contains(key))
                return;
              const auto global_key = std::tuple(t.side, std::min(s, d),
                                                 std::max(s, d), len, strand);
              add_budget(!global_self_keys.contains(global_key));
              self_keys.insert(key);
              global_self_keys.insert(global_key);
              CentromereCopySelfRepeat repeat;
              repeat.id = result.self_repeats.size() + repeats.size() + 1;
              repeat.task_id = t.id;
              repeat.side = t.side;
              repeat.first_begin = std::min(s, d);
              repeat.second_begin = std::max(s, d);
              repeat.length = len;
              repeat.strand = strand;
              repeat.first_count = s < d ? c.reference_count : c.query_count;
              repeat.second_count = s < d ? c.query_count : c.reference_count;
              repeat.rare_level = c.level;
              if (len != c.pair.match_length)
                repeat.count_scope = "local_lcp_before_period_clip";
              repeats.push_back(repeat);
              for (const auto &a : primary.final) {
                const auto ab = job.query ? a.query_begin : a.reference_begin;
                const auto b = std::max(s, ab),
                           e = std::min(s + len, ab + a.length);
                if (b >= e)
                  continue;
                CentromereCopyAnchor x;
                x.task_id = t.id;
                x.self_repeat_id = repeat.id;
                x.primary_anchor_id = a.id;
                x.length = e - b;
                x.strand =
                    a.strand == strand ? Strand::Forward : Strand::Reverse;
                const auto target_begin =
                    strand == Strand::Forward ? d + b - s : d + len - (e - s);
                const auto other_begin =
                    job.query ? a.reference_begin + (a.strand == Strand::Forward
                                                         ? b - ab
                                                         : ab + a.length - e)
                              : a.query_begin + (a.strand == Strand::Forward
                                                     ? b - ab
                                                     : ab + a.length - e);
                x.reference_begin = job.query ? other_begin : target_begin;
                x.query_begin = job.query ? target_begin : other_begin;
                x.self_source_traversal = job.query ? a.strand : strand;
                x.self_source_oriented_begin =
                    x.self_source_traversal == Strand::Forward ? b
                                                               : seq.size() - e;
                if (Exact(reference, query, x) &&
                    !projected_keys.contains(KeyOf(x))) {
                  add_budget(!global_projected_keys.contains(KeyOf(x)));
                  projected_keys.insert(KeyOf(x));
                  global_projected_keys.insert(KeyOf(x));
                  projected.push_back(x);
                }
              }
            });
      }
      t.self_seconds = Seconds(self_start);
    } catch (const Budget &) {
      exhausted = true;
    }
    t.self_seconds = Seconds(self_start);
    result.self_repeats.insert(result.self_repeats.end(), repeats.begin(),
                               repeats.end());
    t.projected_count = projected.size();
    t.retained_evidence_count = projected.size();
    if (exhausted) {
      t.status = "budget_exhausted";
      if (job.scout)
        t.unresolved_bases = seq.size();
      ++result.budget_exhausted_tasks;
      if (!projected.empty()) {
        ++group_id;
        for (auto a : projected) {
          a.context_evidence = "truncated_candidate_evidence";
          result.ambiguities.push_back(
              {group_id,
               static_cast<std::uint64_t>(result.ambiguities.size() + 1), a, 0,
               0, "budget_exhausted"});
          ++t.ambiguity_count;
        }
      }
      t.total_seconds = Seconds(task_start);
      result.tasks.push_back(t);
      continue;
    }
    if (projected.empty()) {
      t.status = repeats.empty() ? "no_repeat" : "no_projected_repeat";
      t.total_seconds = Seconds(task_start);
      result.tasks.push_back(t);
      continue;
    }
    std::vector<CentromereCopyAnchor> accepted;
    std::vector<CentromereCopyAnchor> task_evidence = projected;
    std::set<Key> evidence_keys;
    for (const auto &x : task_evidence)
      evidence_keys.insert(KeyOf(x));
    t.retained_evidence_count = task_evidence.size();
    auto retain_evidence = [&](const CentromereCopyAnchor &x) {
      if (evidence_keys.insert(KeyOf(x)).second)
        task_evidence.push_back(x);
      t.retained_evidence_count = task_evidence.size();
      if (t.retained_evidence_count > t.candidate_count)
        throw AlignmentError("copy evidence exceeded charged candidate count");
    };
    enum class Phase { Projection, Chain, Recursion, Finished };
    auto phase = Phase::Projection;
    auto phase_start = Clock::now();
    auto finish_phase = [&] {
      const auto elapsed = Seconds(phase_start);
      switch (phase) {
      case Phase::Projection:
        t.projection_seconds = elapsed;
        break;
      case Phase::Chain:
        t.chain_seconds = elapsed;
        break;
      case Phase::Recursion:
        t.recursive_seconds = elapsed;
        break;
      case Phase::Finished:
        break;
      }
    };
    auto next_phase = [&](Phase next) {
      finish_phase();
      phase = next;
      phase_start = Clock::now();
    };

    try {

      projected = Merge(std::move(projected), query.size());
      for (auto &a : projected)
        weight(a);
      next_phase(Phase::Chain);
      const bool context =
          !job.scout && Flanks(job, reference, query, unique_cache);
      std::vector<CentromereCopyAnchor> rivals;
      std::set<Key> rival_keys;
      auto charge_evidence = [&](CentromereCopyAnchor &x) {
        if (!projected_keys.contains(KeyOf(x))) {
          add_budget(!global_projected_keys.contains(KeyOf(x)));
          projected_keys.insert(KeyOf(x));
          global_projected_keys.insert(KeyOf(x));
        }
        x.task_id = t.id;
        x.self_source_oriented_begin = x.reference_begin;
        x.self_source_traversal = Strand::Forward;
        weight(x);
        retain_evidence(x);
      };
      auto retain_rival = [&](CentromereCopyAnchor x) {
        if (rival_keys.contains(KeyOf(x)))
          return;
        charge_evidence(x);
        rival_keys.insert(KeyOf(x));
        rivals.push_back(x);
        t.rival_candidate_count = rivals.size();
      };
      const auto opposite = OppositeGap(job);
      if (!job.scout && opposite.first < opposite.second) {
        const auto rb = job.query ? opposite.first : t.target_begin,
                   re = job.query ? opposite.second : t.target_end;
        const auto qb = job.query ? t.target_begin : opposite.first,
                   qe = job.query ? t.target_end : opposite.second;
        for (auto sign : {Strand::Forward, Strand::Reverse})
          VisitCross(reference, query, rb, re, qb, qe, sign, co, options,
                     [&](CentromereCopyAnchor x, std::uint64_t) {
                       // Rivals need no primary proof to refute a unique claim.
                       // They can never become accepted relations from this
                       // search.
                       // Score all bounded competing candidates, including
                       // interior fragments. Boundary-only evidence can falsely
                       // create a strict winner.
                       x.context_evidence = "opposite_gap_competitor";
                       retain_rival(x);
                     });
      }

      std::map<std::tuple<Position, Position, Length, Strand>, std::uint64_t>
          self_ids;
      for (const auto &x : repeats)
        self_ids[{x.first_begin, x.second_begin, x.length, x.strand}] = x.id;
      auto search_proved = [&](Position rb, Position re, Position qb,
                               Position qe, Strand sign, Position sb,
                               Position se, Strand traversal,
                               bool supplemental) {
        std::vector<CentromereCopyAnchor> found;
        std::set<Key> seen;
        VisitCross(
            reference, query, rb, re, qb, qe, sign, co, options,
            [&](CentromereCopyAnchor raw, std::uint64_t level) {
              for (const auto &root : primary.final) {
                const auto ab =
                    job.query ? root.reference_begin : root.query_begin;
                const auto xb =
                    job.query ? raw.reference_begin : raw.query_begin;
                const auto lo = std::max(ab, xb),
                           hi = std::min(ab + root.length, xb + raw.length);
                if (lo >= hi)
                  continue;
                auto x = raw;
                x.length = hi - lo;
                const auto offset = lo - xb;
                if (job.query) {
                  x.reference_begin = lo;
                  x.query_begin =
                      raw.query_begin + (sign == Strand::Forward
                                             ? offset
                                             : raw.length - offset - x.length);
                } else {
                  x.query_begin = lo;
                  x.reference_begin = raw.reference_begin +
                                      (sign == Strand::Forward
                                           ? offset
                                           : raw.length - offset - x.length);
                }
                const Position physical_source =
                    (job.query ? root.query_begin : root.reference_begin) +
                    (root.strand == Strand::Forward ? lo - ab
                                                    : ab + root.length - hi);
                const auto self_sign =
                    root.strand == sign ? Strand::Forward : Strand::Reverse;
                x.self_source_traversal = job.query ? root.strand : self_sign;
                x.self_source_oriented_begin =
                    x.self_source_traversal == Strand::Forward
                        ? physical_source
                        : seq.size() - physical_source - x.length;
                if (x.self_source_traversal != traversal ||
                    x.self_source_oriented_begin < sb ||
                    x.self_source_oriented_begin + x.length > se ||
                    !Exact(reference, query, x) ||
                    !seen.insert(KeyOf(x)).second)
                  continue;
                const auto target =
                    job.query ? x.query_begin : x.reference_begin;
                if (physical_source < target + x.length &&
                    target < physical_source + x.length)
                  continue;
                for (Length k = 0; k < x.length; ++k) {
                  const auto c = seq.bases[physical_source + k],
                             d = self_sign == Strand::Forward
                                     ? seq.bases[target + k]
                                     : Complement(seq.bases[target + x.length -
                                                            k - 1]);
                  if (!IsCanonicalUpper(c) || c != d)
                    throw AlignmentError(
                        "copy recursive self proof is inconsistent");
                }
                const auto self_key = std::tuple(
                    std::min(physical_source, target),
                    std::max(physical_source, target), x.length, self_sign);
                auto proof = self_ids.find(self_key);
                if (proof == self_ids.end()) {
                  const auto global = std::tuple(
                      t.side, std::min(physical_source, target),
                      std::max(physical_source, target), x.length, self_sign);
                  add_budget(!global_self_keys.contains(global));
                  global_self_keys.insert(global);
                  self_keys.insert(self_key);
                  CentromereCopySelfRepeat ev;
                  ev.id = result.self_repeats.size() + 1;
                  ev.task_id = t.id;
                  ev.side = t.side;
                  ev.first_begin = std::min(physical_source, target);
                  ev.second_begin = std::max(physical_source, target);
                  ev.length = x.length;
                  ev.strand = self_sign;
                  ev.rare_level = level;
                  ev.count_scope = "full_self";
                  const auto pattern =
                      seq.bases.substr(physical_source, x.length);
                  ev.first_count = Count(seq.bases, pattern);
                  ev.second_count =
                      Count(seq.bases, self_sign == Strand::Forward
                                           ? pattern
                                           : Reverse(pattern));
                  if (physical_source > target)
                    std::swap(ev.first_count, ev.second_count);
                  result.self_repeats.push_back(ev);
                  proof = self_ids.emplace(self_key, ev.id).first;
                }
                if (!projected_keys.contains(KeyOf(x))) {
                  add_budget(!global_projected_keys.contains(KeyOf(x)));
                  projected_keys.insert(KeyOf(x));
                  global_projected_keys.insert(KeyOf(x));
                  ++t.projected_count;
                }
                x.task_id = t.id;
                x.self_repeat_id = proof->second;
                x.primary_anchor_id = root.id;
                weight(x);
                found.push_back(x);
                retain_evidence(x);
              }
            },
            supplemental);
        std::sort(found.begin(), found.end(), [](const auto &a, const auto &b) {
          return KeyOf(a) < KeyOf(b);
        });
        return found;
      };
      auto covers = [&](const std::vector<CentromereCopyAnchor> &candidates,
                        const CentromereCopyAnchor &anchor) {
        std::vector<Range> ranges;
        for (const auto &c : candidates)
          if (SameLine(c, anchor, query.size()))
            ranges.emplace_back(c.reference_begin,
                                c.reference_begin + c.length);
        for (auto [b, e] : Union(std::move(ranges)))
          if (b <= anchor.reference_begin &&
              e >= anchor.reference_begin + anchor.length)
            return true;
        return false;
      };
      // A task may contain several target copies. A branch's self-source must
      // remain monotone; reuse of that source starts a separate branch. Allow a
      // single branch to bridge mismatches between exact anchors without
      // filling those gaps in the exported coordinates.
      std::vector<std::vector<CentromereCopyAnchor>> groups;
      groups.push_back(std::move(projected));
      std::vector<Range> supported_copy_spans;
      for (auto group : groups) {

        while (!group.empty()) {
          Check(options);
          ++t.branch_selection_count;
          std::sort(
              group.begin(), group.end(),
              [](const auto &a, const auto &b) { return KeyOf(a) < KeyOf(b); });
          const auto best = Chain(group, query.size(), job.query);
          if (best.ids.empty())
            break;
          Position begin = std::numeric_limits<Position>::max(), end = 0;
          auto target_position = [&](const CentromereCopyAnchor &a) {
            return job.query ? a.query_begin : a.reference_begin;
          };
          for (auto id : best.ids) {
            const auto &a = group[id];
            begin = std::min(begin, target_position(a));
            end = std::max(end, target_position(a) + a.length);
          }
          std::set<std::size_t> outside;
          for (std::size_t i = 0; i < group.size(); ++i) {
            const auto b = target_position(group[i]);
            if (b >= end || b + group[i].length <= begin)
              outside.insert(i);
          }
          // A primary-covered counterpart is still a competing physical copy.
          // Enumerate placements from the longest actual supported anchor, then
          // verify every supported anchor under the same
          // translation/reflection. Inter-anchor gap characters are
          // deliberately not used as evidence. These placements can only veto a
          // unique claim, never seed inference.
          if (!job.scout) {
            const auto &counterpart = job.query ? reference : query;
            Position cb = std::numeric_limits<Position>::max(), ce = 0;
            std::size_t seed_id = best.ids.front();
            for (auto id : best.ids) {
              const auto &a = group[id];
              const auto b = job.query ? a.reference_begin : a.query_begin;
              cb = std::min(cb, b);
              ce = std::max(ce, b + a.length);
              if (a.length > group[seed_id].length ||
                  (a.length == group[seed_id].length &&
                   KeyOf(a) < KeyOf(group[seed_id])))
                seed_id = id;
            }
            const auto &seed = group[seed_id];
            const auto seed_begin =
                job.query ? seed.reference_begin : seed.query_begin;
            const auto span = ce - cb;
            const auto original_pattern =
                counterpart.bases.substr(seed_begin, seed.length);
            for (bool reflected : {false, true}) {
              const auto pattern =
                  reflected ? Reverse(original_pattern) : original_pattern;
              const auto seed_offset =
                  reflected ? span - (seed_begin - cb) - seed.length
                            : seed_begin - cb;
              for (std::size_t at = 0; at <= counterpart.size();) {
                Check(options);
                at = counterpart.bases.find(pattern, at);
                if (at == std::string::npos)
                  break;
                const auto occurrence = static_cast<Position>(at);
                ++at;
                if (occurrence < seed_offset)
                  continue;
                const auto placement = occurrence - seed_offset;
                if (span > counterpart.size() - placement ||
                    (!reflected && placement == cb))
                  continue;
                auto project = [&](const CentromereCopyAnchor &a) {
                  auto x = a;
                  const auto old =
                      job.query ? a.reference_begin : a.query_begin;
                  const auto mapped =
                      placement +
                      (reflected ? span - (old - cb) - a.length : old - cb);
                  if (job.query)
                    x.reference_begin = mapped;
                  else
                    x.query_begin = mapped;
                  if (reflected)
                    x.strand = a.strand == Strand::Forward ? Strand::Reverse
                                                           : Strand::Forward;
                  x.id = x.branch_id = x.chain_index = 0;
                  x.self_repeat_id = x.primary_anchor_id = 0;
                  x.context_evidence = "equivalent_supported_chain_competitor";
                  return x;
                };
                // Even a placement rejected by a later anchor consumes the
                // inspected exact seed relation budget. Never retain before
                // charging, and never allow a partial task to survive
                // exhaustion.
                auto candidate_seed = project(seed);
                if (!Exact(reference, query, candidate_seed))
                  throw AlignmentError(
                      "equivalent copy seed projection failed");
                charge_evidence(candidate_seed);
                bool supported = true;
                for (auto id : best.ids) {
                  auto candidate = project(group[id]);
                  if (!Exact(reference, query, candidate)) {
                    supported = false;
                    break;
                  }
                  // A later failed anchor must not make its already verified
                  // exact prefix relations free. Full-chain acceptance and
                  // unique-candidate accounting are separate decisions.
                  charge_evidence(candidate);
                }
                if (!supported)
                  continue;
                for (auto id : best.ids)
                  retain_rival(project(group[id]));
              }
            }
          }
          double runner = 0;
          std::vector<CentromereCopyAnchor> runner_witness;
          for (auto id : best.ids) {
            auto excluded = outside;
            excluded.insert(id);
            auto alternative = Chain(group, query.size(), job.query, excluded);
            runner = std::max(runner, alternative.score);
          }
          for (auto sign : {Strand::Forward, Strand::Reverse}) {
            std::vector<CentromereCopyAnchor> local_rivals;
            for (const auto &original : rivals) {
              if (original.strand != sign)
                continue;
              const auto b = job.query ? original.query_begin
                                       : original.reference_begin,
                         e = b + original.length;
              const auto lo = std::max(b, begin), hi = std::min(e, end);
              if (lo >= hi)
                continue;
              auto x = original;
              x.length = hi - lo;
              const auto offset = lo - b;
              if (job.query) {
                x.query_begin = lo;
                x.reference_begin = original.reference_begin +
                                    (sign == Strand::Forward
                                         ? offset
                                         : original.length - offset - x.length);
              } else {
                x.reference_begin = lo;
                x.query_begin = original.query_begin +
                                (sign == Strand::Forward
                                     ? offset
                                     : original.length - offset - x.length);
              }
              x.self_source_oriented_begin = x.reference_begin;
              weight(x);
              local_rivals.push_back(x);
            }
            std::sort(local_rivals.begin(), local_rivals.end(),
                      [](const auto &a, const auto &b) {
                        return KeyOf(a) < KeyOf(b);
                      });
            const auto rival_path = Chain(local_rivals, query.size(), false);
            if (!rival_path.ids.empty() && rival_path.score >= runner) {
              runner_witness.clear();
              for (auto id : rival_path.ids)
                runner_witness.push_back(local_rivals[id]);
            }
            runner = std::max(runner, rival_path.score);
          }
          const double margin = best.score - runner;
          bool strict = margin > 1e-12 * std::max(1.0, std::abs(best.score));
          bool local_competition = false;
          Position local_begin = 0, local_end = 0;
          double local_selected_score = 0, local_rival_score = 0;
          std::vector<CentromereCopyAnchor> local_witness;
          // A long uncontested prefix must not buy certainty for a disputed
          // tail. Compare existing mutually exclusive evidence on the same
          // physical target span; all shortened weights are recomputed.
          auto compare_overlap =
              [&](const std::vector<CentromereCopyAnchor> &path,
                  bool projected_proof) {
                if (path.empty() || local_competition)
                  return;
                Position lo = std::numeric_limits<Position>::max(), hi = 0;
                for (const auto &x : path) {
                  const auto b = target_position(x);
                  lo = std::min(lo, b);
                  hi = std::max(hi, b + x.length);
                }
                lo = std::max(lo, begin);
                hi = std::min(hi, end);
                if (lo >= hi)
                  return;
                std::vector<CentromereCopyAnchor> proposed, alternative;
                auto score_clipped = [&](CentromereCopyAnchor &x) {
                  if (!projected_keys.contains(KeyOf(x))) {
                    add_budget(!global_projected_keys.contains(KeyOf(x)));
                    projected_keys.insert(KeyOf(x));
                    global_projected_keys.insert(KeyOf(x));
                  }
                  weight(x);
                  retain_evidence(x);
                };
                for (auto id : best.ids) {
                  auto x = ClipTarget(group[id], lo, hi, job.query);
                  if (x.length) {
                    score_clipped(x);
                    proposed.push_back(x);
                  }
                }
                for (const auto &original : path) {
                  auto x = ClipTarget(original, lo, hi, job.query);
                  if (x.length) {
                    score_clipped(x);
                    if (!projected_proof)
                      x.self_source_oriented_begin = x.reference_begin;
                    alternative.push_back(x);
                  }
                }
                if (proposed.empty() || alternative.empty())
                  return;
                // Use the existing displacement penalty as a contextual
                // evidence requirement. An incidental few-base collision cannot
                // locate a different copy merely because both clipped strings
                // are equal.
                double placement_cost = std::numeric_limits<double>::infinity();
                bool conflict = false;
                for (const auto &a : alternative)
                  for (const auto &b : proposed) {
                    const auto ab = target_position(a), bb = target_position(b);
                    const auto overlap_begin = std::max(ab, bb),
                               finish = std::min(ab + a.length, bb + b.length);
                    if (overlap_begin >= finish || SameLine(a, b, query.size()))
                      continue;
                    conflict = true;
                    const auto x = Counterpart(a, overlap_begin, job.query),
                               y = Counterpart(b, overlap_begin, job.query);
                    const auto shift = x > y ? x - y : y - x;
                    placement_cost = std::min(
                        placement_cost,
                        2.0 * std::log2(static_cast<double>(shift) + 1.0));
                  }
                if (!conflict)
                  return;
                auto order = [](const auto &a, const auto &b) {
                  return KeyOf(a) < KeyOf(b);
                };
                std::sort(proposed.begin(), proposed.end(), order);
                std::sort(alternative.begin(), alternative.end(), order);
                const double ps =
                    Chain(proposed, query.size(), job.query).score;
                const double as = Chain(alternative, query.size(),
                                        projected_proof ? job.query : false)
                                      .score;
                if (as <= placement_cost)
                  return;
                const double eps =
                    1e-12 * std::max({1.0, std::abs(ps), std::abs(as)});
                if (as + eps >= ps) {
                  local_competition = true;
                  local_begin = lo;
                  local_end = hi;
                  local_selected_score = ps;
                  local_rival_score = as;
                  local_witness = std::move(alternative);
                }
              };
          if (strict) {
            std::vector<CentromereCopyAnchor> alternatives;
            for (const auto &candidate : group) {
              // A selected path may change diagonal across a gap. A candidate
              // on the prefix diagonal can still conflict with its tail: only
              // actual target overlap determines whether it is an alternative.
              bool conflicts = false;
              for (auto id : best.ids) {
                const auto &selected = group[id];
                const auto cb = target_position(candidate),
                           sb = target_position(selected);
                if (std::max(cb, sb) <
                        std::min(cb + candidate.length, sb + selected.length) &&
                    !SameLine(candidate, selected, query.size())) {
                  conflicts = true;
                  break;
                }
              }
              if (conflicts)
                alternatives.push_back(candidate);
            }
            // Single supported exact blocks and multi-anchor competing paths
            // are both eligible; disjoint target copies never enter a window.
            for (const auto &candidate : alternatives) {
              compare_overlap({candidate}, true);
              if (local_competition)
                break;
            }
            if (!local_competition) {
              const auto path = Chain(alternatives, query.size(), job.query);
              std::vector<CentromereCopyAnchor> competing;
              for (auto id : path.ids)
                competing.push_back(alternatives[id]);
              compare_overlap(competing, true);
            }
            for (const auto &rival : rivals) {
              compare_overlap({rival}, false);
              if (local_competition)
                break;
            }
            if (!local_competition) {
              // Veto-only opposite-gap evidence has the same local chaining
              // opportunity as projected evidence. Several fragments can
              // jointly locate a competing copy even when no single fragment
              // does.
              for (auto sign : {Strand::Forward, Strand::Reverse}) {
                std::vector<CentromereCopyAnchor> signed_rivals;
                for (const auto &rival : rivals)
                  if (rival.strand == sign) {
                    auto x = rival;
                    x.self_source_oriented_begin = x.reference_begin;
                    signed_rivals.push_back(x);
                  }
                std::sort(signed_rivals.begin(), signed_rivals.end(),
                          [](const auto &a, const auto &b) {
                            return KeyOf(a) < KeyOf(b);
                          });
                const auto chain = Chain(signed_rivals, query.size(), false);
                std::vector<CentromereCopyAnchor> competing;
                for (auto id : chain.ids)
                  competing.push_back(signed_rivals[id]);
                compare_overlap(competing, false);
                if (local_competition)
                  break;
              }
            }
            if (local_competition) {
              strict = false;
              ++group_id;
              std::uint64_t witness_order = 0;
              for (auto x : local_witness) {
                const std::string route(EvidenceRoute(x));
                x.context_evidence = "same_target_overlap_" +
                                     std::to_string(local_begin) + "_" +
                                     std::to_string(local_end);
                result.ambiguities.push_back(
                    {group_id,
                     static_cast<std::uint64_t>(result.ambiguities.size() + 1),
                     x, x.weight, local_selected_score - local_rival_score,
                     "local_overlap_competition", local_selected_score,
                     local_rival_score, true, local_begin, local_end});
                auto &diagnostic = result.ambiguities.back();
                diagnostic.evidence_role = "local_runner_witness";
                diagnostic.evidence_route = route;
                diagnostic.witness_order = ++witness_order;
                diagnostic.witness_size = local_witness.size();
                ++t.ambiguity_count;
              }
            }
          }
          bool local_verified = true;
          Position source_begin = std::numeric_limits<Position>::max(),
                   source_end = 0;
          for (auto id : best.ids) {
            const auto &a = group[id];
            const auto sb =
                a.self_source_traversal == Strand::Forward
                    ? a.self_source_oriented_begin
                    : seq.size() - a.self_source_oriented_begin - a.length;
            source_begin = std::min(source_begin, sb);
            source_end = std::max(source_end, sb + a.length);
          }
          const bool tandem_boundary =
              source_end == t.target_begin || source_begin == t.target_end;
          if (context && strict && tandem_boundary) {
            Position rb = std::numeric_limits<Position>::max(), re = 0, qb = rb,
                     qe = 0, sb = rb, se = 0;
            for (auto id : best.ids) {
              const auto &a = group[id];
              rb = std::min(rb, a.reference_begin);
              re = std::max(re, a.reference_begin + a.length);
              qb = std::min(qb, a.query_begin);
              qe = std::max(qe, a.query_begin + a.length);
              sb = std::min(sb, a.self_source_oriented_begin);
              se = std::max(se, a.self_source_oriented_begin + a.length);
            }
            ++t.recursion_node_count;
            const auto window =
                search_proved(rb, re, qb, qe, group[best.ids[0]].strand, sb, se,
                              group[best.ids[0]].self_source_traversal, false);
            const auto local_path = Chain(window, query.size(), job.query);
            std::vector<CentromereCopyAnchor> local_chain;
            for (auto id : local_path.ids)
              local_chain.push_back(window[id]);
            for (auto id : best.ids)
              local_verified = local_verified && covers(local_chain, group[id]);
          }

          if (context && strict && local_verified && tandem_boundary) {
            supported_copy_spans.emplace_back(begin, end);
            ++branch_id;
            std::uint64_t index = 0;
            for (auto id : best.ids) {
              auto a = group[id];
              a.branch_id = branch_id;
              a.chain_index = ++index;
              a.context_evidence = "two_unique_primary_flanks";
              if (emitted.insert(KeyOf(a)).second)
                accepted.push_back(a);
            }
          } else {
            ++group_id;
            const auto reason =
                job.scout           ? "competing_main_context"
                : !context          ? "insufficient_context"
                : local_competition ? "local_overlap_branch_ambiguous"
                : !strict           ? "equal_support"
                : !local_verified   ? "no_local_support"
                                    : "unresolved_tandem_boundary";
            for (std::size_t i = 0; i < group.size(); ++i) {
              if (outside.contains(i))
                continue;
              auto a = group[i];
              a.context_evidence =
                  context ? "two_unique_primary_flanks" : "unresolved";
              result.ambiguities.push_back(
                  {group_id,
                   static_cast<std::uint64_t>(result.ambiguities.size() + 1), a,
                   a.weight, margin, reason, best.score, runner, true});
              ++t.ambiguity_count;
            }
            // Preserve the already scored rival path as diagnostics only.
            // Append at the end so existing candidate IDs remain unchanged.
            std::uint64_t witness_order = 0;
            for (const auto &a : runner_witness) {
              CentromereCopyAmbiguity diagnostic{
                  group_id,   0,      a,    a.weight, margin, reason,
                  best.score, runner, true, begin,    end};
              diagnostic.evidence_role = "runner_witness";
              diagnostic.evidence_route = EvidenceRoute(a);
              diagnostic.witness_order = ++witness_order;
              diagnostic.witness_size = runner_witness.size();
              pending_runner_witnesses.push_back(std::move(diagnostic));
              ++t.ambiguity_count;
            }
          }
          // Recurse only on not-yet-processed target pieces. No accepted branch
          // is used as a new projection root; every retained piece keeps
          // original proof.
          std::vector<CentromereCopyAnchor> remaining;
          std::set<Key> remaining_keys;
          for (const auto &a : group) {
            const auto b = target_position(a), e = b + a.length;
            for (const auto &part :
                 {Range{b, std::min(e, begin)}, Range{std::max(b, end), e}}) {
              if (part.first >= part.second)
                continue;
              auto x = a;
              const auto offset = part.first - b;
              x.length = part.second - part.first;
              if (job.query) {
                x.query_begin = part.first;
                x.reference_begin =
                    a.reference_begin + (a.strand == Strand::Forward
                                             ? offset
                                             : a.length - offset - x.length);
              } else {
                x.reference_begin = part.first;
                x.query_begin =
                    a.query_begin + (a.strand == Strand::Forward
                                         ? offset
                                         : a.length - offset - x.length);
              }
              x.self_source_oriented_begin =
                  a.self_source_oriented_begin +
                  (x.reference_begin - a.reference_begin);
              if (!projected_keys.contains(KeyOf(x))) {
                add_budget(!global_projected_keys.contains(KeyOf(x)));
                projected_keys.insert(KeyOf(x));
                global_projected_keys.insert(KeyOf(x));
              }
              weight(x);
              if (remaining_keys.insert(KeyOf(x)).second)
                remaining.push_back(x);
              retain_evidence(x);
            }
          }
          group = std::move(remaining);
        }
      }
      // Outer primary flanks locate the entire gap, not arbitrary short matches
      // inside it. Require independently supported copy-chain spans to explain
      // that structural gap. Internal mismatches remain unanchored, not filled.
      const auto spans = Union(supported_copy_spans);
      const bool gap_supported = spans.size() == 1 &&
                                 spans[0].first == t.target_begin &&
                                 spans[0].second == t.target_end;
      if (!gap_supported && !accepted.empty()) {
        ++group_id;
        for (auto a : accepted) {
          emitted.erase(KeyOf(a));
          result.ambiguities.push_back(
              {group_id,
               static_cast<std::uint64_t>(result.ambiguities.size() + 1), a,
               a.weight, 0, "incomplete_tandem_context"});
          ++t.ambiguity_count;
        }
        accepted.clear();
      }
      next_phase(Phase::Recursion);
      struct GapNode {
        Position rb, re, ub, ue, sb, se;
        Strand sign, traversal;
        std::uint64_t branch, depth;
      };
      std::map<std::uint64_t, std::vector<CentromereCopyAnchor>> branches;
      for (const auto &a : accepted)
        branches[a.branch_id].push_back(a);
      std::vector<GapNode> pending;
      for (auto &[id, anchors] : branches) {
        std::sort(anchors.begin(), anchors.end(),
                  [](const auto &a, const auto &b) {
                    return a.reference_begin < b.reference_begin;
                  });
        for (std::size_t i = 1; i < anchors.size(); ++i) {
          const auto &l = anchors[i - 1];
          const auto &r = anchors[i];
          pending.push_back({l.reference_begin + l.length, r.reference_begin,
                             Oriented(l, query.size()) + l.length,
                             Oriented(r, query.size()),
                             l.self_source_oriented_begin + l.length,
                             r.self_source_oriented_begin, l.strand,
                             l.self_source_traversal, id, 1});
        }
      }
      while (!pending.empty()) {
        Check(options);
        const auto node = pending.back();
        pending.pop_back();
        if (node.rb >= node.re || node.ub >= node.ue || node.sb >= node.se)
          continue;
        ++t.recursion_node_count;
        t.max_recursion_depth = std::max(t.max_recursion_depth, node.depth);
        const auto qb = node.sign == Strand::Forward ? node.ub
                                                     : query.size() - node.ue,
                   qe = node.sign == Strand::Forward ? node.ue
                                                     : query.size() - node.ub;
        auto candidates = search_proved(node.rb, node.re, qb, qe, node.sign,
                                        node.sb, node.se, node.traversal, true);
        const auto best = Chain(candidates, query.size(), job.query);
        if (best.ids.empty())
          continue;
        double runner = 0;
        for (auto id : best.ids)
          runner = std::max(
              runner, Chain(candidates, query.size(), job.query, {id}).score);
        if (best.score - runner <=
            1e-12 * std::max(1.0, std::abs(best.score))) {
          ++group_id;
          for (const auto &x : candidates) {
            result.ambiguities.push_back(
                {group_id,
                 static_cast<std::uint64_t>(result.ambiguities.size() + 1), x,
                 x.weight, best.score - runner, "recursive_equal_support",
                 best.score, runner, true});
            ++t.ambiguity_count;
          }
          continue;
        }
        Position rb = node.rb, ub = node.ub, sb = node.sb;
        for (auto id : best.ids) {
          auto x = candidates[id];
          pending.push_back({rb, x.reference_begin, ub,
                             Oriented(x, query.size()), sb,
                             x.self_source_oriented_begin, node.sign,
                             node.traversal, node.branch, node.depth + 1});
          rb = x.reference_begin + x.length;
          ub = Oriented(x, query.size()) + x.length;
          sb = x.self_source_oriented_begin + x.length;
          x.branch_id = node.branch;
          x.context_evidence = "original_primary_self_repeat_bounded_recursion";
          if (emitted.insert(KeyOf(x)).second) {
            accepted.push_back(x);
            ++t.recursive_added_count;
          }
        }
        pending.push_back({rb, node.re, ub, node.ue, sb, node.se, node.sign,
                           node.traversal, node.branch, node.depth + 1});
      }
      next_phase(Phase::Finished);
      std::sort(
          accepted.begin(), accepted.end(), [](const auto &a, const auto &b) {
            return std::tuple(a.branch_id, a.reference_begin, a.query_begin) <
                   std::tuple(b.branch_id, b.reference_begin, b.query_begin);
          });
      std::uint64_t previous_branch = 0, index = 0;
      for (auto &a : accepted) {
        if (a.branch_id != previous_branch) {
          previous_branch = a.branch_id;
          index = 0;
        }
        a.chain_index = ++index;
      }
      for (auto a : accepted) {
        a.id = result.anchors.size() + 1;
        result.anchors.push_back(a);
        ++t.accepted_count;
      }
      if (!job.scout)
        t.unresolved_bases = (t.target_end - t.target_begin) -
                             std::min<Length>(t.target_end - t.target_begin,
                                              Covered(accepted, job.query));

      t.total_seconds = Seconds(task_start);
      t.status = t.accepted_count    ? "accepted"
                 : t.ambiguity_count ? "ambiguous"
                                     : "no_local_support";
      result.tasks.push_back(t);
    } catch (const Budget &) {
      finish_phase();
      for (const auto &a : accepted)
        emitted.erase(KeyOf(a));
      accepted.clear();
      ++result.budget_exhausted_tasks;
      ++group_id;
      for (auto x : task_evidence) {
        x.context_evidence = "truncated_candidate_evidence";
        result.ambiguities.push_back(
            {group_id,
             static_cast<std::uint64_t>(result.ambiguities.size() + 1), x, 0, 0,
             "budget_exhausted"});
        ++t.ambiguity_count;
      }
      t.accepted_count = 0;
      t.recursive_added_count = 0;
      t.unresolved_bases =
          job.scout ? seq.size() : t.target_end - t.target_begin;
      t.status = "budget_exhausted";
      t.total_seconds = Seconds(task_start);
      result.tasks.push_back(t);
    }
  }
  if (result.total_candidates !=
      global_self_keys.size() + global_projected_keys.size())
    throw AlignmentError("copy candidate accounting is not transactional");
  for (auto &diagnostic : pending_runner_witnesses) {
    diagnostic.candidate_id = result.ambiguities.size() + 1;
    result.ambiguities.push_back(std::move(diagnostic));
  }
  result.total_seconds = Seconds(start);
  result.complete = true;
  return result;
#else
  (void)reference;
  (void)query;
  (void)primary;
  (void)co;
  (void)options;
  throw DependencyError("centromere backend not built; copy anchors require "
                        "RAMAG_ENABLE_CENTROMERE");
#endif
}
void WriteCentromereCopyAnchors(std::ostream &o,
                                const CentromereCopyAnchorResult &r) {
  o << "anchor_id\tbranch_id\tchain_index\treference_begin\treference_"
       "end\tquery_begin\tquery_end\tstrand\tlength\tweight\ttask_id\tself_"
       "repeat_id\tprimary_anchor_id\tcontext_evidence\n"
    << std::setprecision(17);
  for (const auto &a : r.anchors)
    o << a.id << '\t' << a.branch_id << '\t' << a.chain_index << '\t'
      << a.reference_begin << '\t' << a.reference_begin + a.length << '\t'
      << a.query_begin << '\t' << a.query_begin + a.length << '\t'
      << Sign(a.strand) << '\t' << a.length << '\t' << a.weight << '\t'
      << a.task_id << '\t' << a.self_repeat_id << '\t' << a.primary_anchor_id
      << '\t' << a.context_evidence << '\n';
}
void WriteCentromereCopyAmbiguities(std::ostream &o,
                                    const CentromereCopyAnchorResult &r) {
  o << "group_id\tcandidate_id\treference_begin\treference_end\tquery_"
       "begin\tquery_end\tstrand\tlength\tweight\ttask_id\tself_repeat_"
       "id\tprimary_anchor_id\tcandidate_anchor_score\tgroup_score_"
       "margin\treason\tgroup_best_chain_"
       "score\trunner_chain_score\tscored\tevaluation_target_begin\tevaluation_"
       "target_end\tevidence_role\tevidence_route\twitness_order\twitness_"
       "size\n"
    << std::setprecision(17);
  for (const auto &x : r.ambiguities) {
    const auto &a = x.anchor;
    o << x.group_id << '\t' << x.candidate_id << '\t' << a.reference_begin
      << '\t' << a.reference_begin + a.length << '\t' << a.query_begin << '\t'
      << a.query_begin + a.length << '\t' << Sign(a.strand) << '\t' << a.length
      << '\t' << a.weight << '\t' << a.task_id << '\t' << a.self_repeat_id
      << '\t' << a.primary_anchor_id << '\t';
    if (x.scored)
      o << x.score << '\t' << x.score_margin;
    else
      o << "NA\tNA";
    o << '\t' << x.reason << '\t';
    if (x.scored)
      o << x.group_best_chain_score << '\t' << x.runner_chain_score;
    else
      o << "NA\tNA";
    o << '\t' << x.scored << '\t';
    if (x.evaluation_target_begin < x.evaluation_target_end)
      o << x.evaluation_target_begin << '\t' << x.evaluation_target_end;
    else
      o << "NA\tNA";
    o << '\t' << x.evidence_role << '\t'
      << (x.evidence_route.empty() ? EvidenceRoute(a)
                                   : std::string_view(x.evidence_route))
      << '\t';
    if (x.witness_order && x.witness_size)
      o << x.witness_order << '\t' << x.witness_size;
    else
      o << "NA\tNA";
    o << '\n';
  }
}
void WriteCentromereCopySelfRepeats(std::ostream &o,
                                    const CentromereCopyAnchorResult &r) {
  o << "self_repeat_id\ttask_id\tside\tfirst_begin\tfirst_end\tsecond_"
       "begin\tsecond_end\tstrand\tlength\tfirst_count\tsecond_count\trare_"
       "level\tcount_scope\n";
  for (const auto &a : r.self_repeats)
    o << a.id << '\t' << a.task_id << '\t' << a.side << '\t' << a.first_begin
      << '\t' << a.first_begin + a.length << '\t' << a.second_begin << '\t'
      << a.second_begin + a.length << '\t' << Sign(a.strand) << '\t' << a.length
      << '\t' << a.first_count << '\t' << a.second_count << '\t' << a.rare_level
      << '\t' << a.count_scope << '\n';
}
void WriteCentromereCopyTasks(std::ostream &o,
                              const CentromereCopyAnchorResult &r) {
  o << "task_id\tside\tsource_begin\tsource_end\ttarget_begin\ttarget_"
       "end\tcandidate_count\tprojected_count\taccepted_count\tambiguity_"
       "count\trecursion_node_count\tmax_recursion_depth\tself_"
       "seconds\tprojection_seconds\tchain_seconds\ttotal_"
       "seconds\tstatus\tunresolved_bases\tbranch_selection_count\trecursive_"
       "added_count\trecursive_seconds\tretained_evidence_count\trival_"
       "candidate_count\n"
    << std::setprecision(17);
  for (const auto &t : r.tasks)
    o << t.id << '\t' << t.side << '\t' << t.source_begin << '\t'
      << t.source_end << '\t' << t.target_begin << '\t' << t.target_end << '\t'
      << t.candidate_count << '\t' << t.projected_count << '\t'
      << t.accepted_count << '\t' << t.ambiguity_count << '\t'
      << t.recursion_node_count << '\t' << t.max_recursion_depth << '\t'
      << t.self_seconds << '\t' << t.projection_seconds << '\t'
      << t.chain_seconds << '\t' << t.total_seconds << '\t' << t.status << '\t'
      << t.unresolved_bases << '\t' << t.branch_selection_count << '\t'
      << t.recursive_added_count << '\t' << t.recursive_seconds << '\t'
      << t.retained_evidence_count << '\t' << t.rival_candidate_count << '\n';
}
void WriteCentromereCopySummary(std::ostream &o,
                                const CentromereCopyAnchorResult &r) {
  o << "key\tvalue\ncomplete\t" << r.complete << "\naccepted_anchor_count\t"
    << r.anchors.size() << "\nambiguity_count\t" << r.ambiguities.size()
    << "\nself_repeat_count\t" << r.self_repeats.size() << "\ntask_count\t"
    << r.tasks.size() << "\ntotal_candidates\t" << r.total_candidates
    << "\nbudget_exhausted_tasks\t" << r.budget_exhausted_tasks
    << "\nwfa_calls\t" << r.wfa_calls << "\ntotal_seconds\t"
    << std::setprecision(17) << r.total_seconds << '\n';
}
} // namespace ramag

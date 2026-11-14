#include "common/centromere.hpp"
#include "ramag/centromere_stitch.hpp"
#include "centromere/stitch/kernel.hpp"
#include "centromere/stitch/fixed_gaps.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <exception>
#include <iomanip>
#include <limits>
#include <map>
#include <ostream>
#include <set>
#include <tuple>
#include <utility>

namespace ramag {
namespace {
using Clock = std::chrono::steady_clock;
using Interval = std::pair<Position, Position>;

struct Anchor {
  std::uint64_t id{}, index{};
  Position rb{}, qb{};
  Length length{};
  Strand strand{Strand::Forward};
};
struct Path {
  std::uint64_t id{}, region_id{}, branch_id{};
  std::string role;
  std::vector<Anchor> anchors;
  std::vector<std::uint64_t> gaps;  // zero means adjacent anchors with no gap
};

double Seconds(Clock::time_point start) {
  return std::chrono::duration<double>(Clock::now() - start).count();
}
using ::ramag::detail::IsCanonicalUpper;
char Complement(char c) {
  switch (c) {
    case 'A': return 'T';
    case 'C': return 'G';
    case 'G': return 'C';
    case 'T': return 'A';
    default: return 'N';
  }
}
char Sign(Strand strand) { return strand == Strand::Forward ? '+' : '-'; }
void Append(std::vector<CigarOp>& cigar, char op, Length length) {
  if (!length) return;
  if (!cigar.empty() && cigar.back().operation == op) {
    if (length > std::numeric_limits<Length>::max() - cigar.back().length)
      throw AlignmentError("stitch CIGAR length overflow");
    cigar.back().length += length;
  } else {
    cigar.push_back({op, length});
  }
}
Length UnionLength(std::vector<Interval> intervals) {
  std::sort(intervals.begin(), intervals.end());
  Length result = 0;
  Position begin = 0, end = 0;
  for (const auto& interval : intervals) {
    if (interval.first >= interval.second) continue;
    if (interval.first > end) {
      result += end - begin;
      begin = interval.first;
      end = interval.second;
    } else {
      end = std::max(end, interval.second);
    }
  }
  return result + end - begin;
}
void ValidateSequence(const SequenceRecord& sequence) {
  for (char c : sequence.bases)
    if (!IsCanonicalUpper(c) && c != 'N')
      throw AlignmentError("stitch input must contain normalized uppercase ACGTN");
}
void ValidateAnchor(const Anchor& anchor, const SequenceRecord& reference,
                    const SequenceRecord& query) {
  if (!anchor.id || !anchor.index || !anchor.length ||
      anchor.rb > reference.size() || anchor.length > reference.size() - anchor.rb ||
      anchor.qb > query.size() || anchor.length > query.size() - anchor.qb ||
      (anchor.strand != Strand::Forward && anchor.strand != Strand::Reverse))
    throw AlignmentError("invalid frozen stitch anchor identity, bounds or strand");
  for (Length i = 0; i < anchor.length; ++i) {
    const char r = reference.bases[anchor.rb + i];
    const char q = anchor.strand == Strand::Forward
        ? query.bases[anchor.qb + i]
        : Complement(query.bases[anchor.qb + anchor.length - i - 1]);
    if (!IsCanonicalUpper(r) || r != q)
      throw AlignmentError("frozen stitch anchor is not a canonical exact match");
  }
}
void ValidatePath(Path& path, const SequenceRecord& reference,
                  const SequenceRecord& query) {
  std::sort(path.anchors.begin(), path.anchors.end(),
            [](const Anchor& a, const Anchor& b) { return a.index < b.index; });
  for (std::size_t i = 0; i < path.anchors.size(); ++i) {
    const auto& anchor = path.anchors[i];
    ValidateAnchor(anchor, reference, query);
    if (anchor.index != i + 1)
      throw AlignmentError("stitch path chain_index must be unique and consecutive from 1");
    if (!i) continue;
    const auto& previous = path.anchors[i - 1];
    const bool query_order = anchor.strand == Strand::Forward
        ? previous.qb + previous.length <= anchor.qb
        : anchor.qb + anchor.length <= previous.qb;
    if (anchor.strand != previous.strand || previous.rb + previous.length > anchor.rb ||
        !query_order)
      throw AlignmentError("frozen stitch path is not strand-consistent and nonoverlapping");
  }
}
std::vector<Path> MakePaths(const SequenceRecord& reference, const SequenceRecord& query,
                           const CentromereDirectedAnchorResult& primary,
                           const CentromereCopyAnchorResult* copies,
                           CentromereStitchResult& out) {
  if (!primary.complete || (copies && !copies->complete))
    throw AlignmentError("stitch requires completed frozen anchor stages");
  std::map<std::uint64_t, CentromereDirectedRegion> regions;
  std::vector<Interval> known;
  for (const auto& region : primary.regions) {
    if (!region.id || region.query_begin >= region.query_end ||
        region.query_end > query.size() ||
        (region.strand != Strand::Forward && region.strand != Strand::Reverse) ||
        !regions.emplace(region.id, region).second)
      throw AlignmentError("invalid frozen stitch direction region");
    known.emplace_back(region.query_begin, region.query_end);
  }
  std::sort(known.begin(), known.end());
  for (std::size_t i = 1; i < known.size(); ++i)
    if (known[i].first < known[i - 1].second)
      throw AlignmentError("frozen stitch query direction regions overlap");
  out.unknown_query_bases = query.size() - UnionLength(known);

  std::map<std::uint64_t, Path> main_paths, copy_paths;
  std::set<std::uint64_t> primary_ids, copy_ids;
  for (const auto& anchor : primary.final) {
    const auto found = regions.find(anchor.region_id);
    if (found == regions.end() || !primary_ids.insert(anchor.id).second)
      throw AlignmentError("stitch primary anchor has unknown region or duplicate ID");
    const auto& region = found->second;
    if (anchor.strand != region.strand || anchor.query_begin < region.query_begin ||
        anchor.query_begin > region.query_end ||
        anchor.length > region.query_end - anchor.query_begin)
      throw AlignmentError("stitch primary anchor falls outside its frozen region");
    auto& path = main_paths[anchor.region_id];
    path.role = "primary";
    path.region_id = anchor.region_id;
    path.anchors.push_back({anchor.id, anchor.chain_index, anchor.reference_begin,
                            anchor.query_begin, anchor.length, anchor.strand});
  }
  if (copies) {
    for (const auto& anchor : copies->anchors) {
      if (!anchor.branch_id || !copy_ids.insert(anchor.id).second)
        throw AlignmentError("stitch copy anchor has invalid branch or duplicate ID");
      auto& path = copy_paths[anchor.branch_id];
      path.role = "copy";
      path.branch_id = anchor.branch_id;
      path.anchors.push_back({anchor.id, anchor.chain_index, anchor.reference_begin,
                              anchor.query_begin, anchor.length, anchor.strand});
    }
  }
  out.primary_anchor_count = primary.final.size();
  out.copy_anchor_count = copies ? copies->anchors.size() : 0;
  out.primary_path_count = main_paths.size();
  out.copy_path_count = copy_paths.size();
  out.no_anchor_region_count = regions.size() - main_paths.size();
  std::vector<Path> paths;
  for (auto& item : main_paths) paths.push_back(std::move(item.second));
  for (auto& item : copy_paths) paths.push_back(std::move(item.second));
  std::vector<Interval> reference_support, query_support, primary_query_support;
  for (auto& path : paths) {
    path.id = &path - paths.data() + 1;
    ValidatePath(path, reference, query);
    const auto& first = path.anchors.front();
    const auto& last = path.anchors.back();
    const Interval qi{std::min(first.qb, last.qb),
                      std::max(first.qb + first.length, last.qb + last.length)};
    reference_support.emplace_back(first.rb, last.rb + last.length);
    query_support.push_back(qi);
    if (path.role == "primary") primary_query_support.push_back(qi);
    if (path.anchors.size() == 1) ++out.singleton_path_count;
    for (std::size_t i = 1; i < path.anchors.size(); ++i) {
      const auto& left = path.anchors[i - 1];
      const auto& right = path.anchors[i];
      const Position rb = left.rb + left.length, re = right.rb;
      const Position qb = left.strand == Strand::Forward ? left.qb + left.length
                                                        : right.qb + right.length;
      const Position qe = left.strand == Strand::Forward ? right.qb : left.qb;
      if (rb == re && qb == qe) {
        path.gaps.push_back(0);
        continue;
      }
      CentromereStitchGap gap;
      gap.gap_id = out.gaps.size() + 1;
      gap.path_id = path.id;
      gap.path_role = path.role;
      gap.region_id = path.region_id;
      gap.branch_id = path.branch_id;
      gap.left_anchor_id = left.id;
      gap.right_anchor_id = right.id;
      gap.reference_begin = rb;
      gap.reference_end = re;
      gap.query_begin = qb;
      gap.query_end = qe;
      gap.strand = left.strand;
      gap.status = "pending";
      path.gaps.push_back(gap.gap_id);
      out.gaps.push_back(std::move(gap));
    }
  }
  out.reference_support_bases = UnionLength(reference_support);
  out.query_support_bases = UnionLength(query_support);
  out.known_query_outside_primary_support_bases =
      UnionLength(known) - UnionLength(primary_query_support);
  return paths;
}
std::string_view ReferenceView(const SequenceRecord& reference, Position begin, Position end) {
  return std::string_view(reference.bases).substr(begin, end - begin);
}
std::string_view QueryView(const SequenceRecord& query, const std::string& reverse_query,
                          Position begin, Position end, Strand strand) {
  if (strand == Strand::Forward)
    return std::string_view(query.bases).substr(begin, end - begin);
  return std::string_view(reverse_query).substr(query.size() - end, end - begin);
}
void Refresh(CentromereStitchResult& out, Clock::time_point started) {
  std::vector<Interval> reference_aligned, query_aligned, reference_paired, query_paired;
  for (const auto& segment : out.segments) {
    const auto& alignment = segment.alignment;
    reference_aligned.emplace_back(alignment.reference_begin, alignment.reference_end);
    query_aligned.emplace_back(alignment.query_begin, alignment.query_end);
    Position r = alignment.reference_begin;
    Position q = alignment.strand == Strand::Forward ? alignment.query_begin : alignment.query_end;
    for (const auto& op : alignment.cigar) {
      if (op.operation == '=' || op.operation == 'X') {
        reference_paired.emplace_back(r, r + op.length);
        query_paired.emplace_back(alignment.strand == Strand::Forward ? q : q - op.length,
                                  alignment.strand == Strand::Forward ? q + op.length : q);
      }
      if (op.operation != 'I') r += op.length;
      if (op.operation != 'D') {
        if (alignment.strand == Strand::Forward) q += op.length;
        else q -= op.length;
      }
    }
  }
  out.reference_aligned_bases = UnionLength(reference_aligned);
  out.query_aligned_bases = UnionLength(query_aligned);
  out.reference_paired_bases = UnionLength(reference_paired);
  out.query_paired_bases = UnionLength(query_paired);
  out.reference_unaligned_bases = out.reference_length - out.reference_aligned_bases;
  out.query_unaligned_bases = out.query_length - out.query_aligned_bases;
  out.total_seconds = Seconds(started);
}
void ValidateSegment(CentromereStitchSegment& segment, const SequenceRecord& reference,
                     const SequenceRecord& query, const std::string& reverse_query,
                     CentromereStitchScoreModel model) {
  auto& alignment = segment.alignment;
  alignment.score = detail::ScoreCentromereStitchCigar(
      ReferenceView(reference, alignment.reference_begin, alignment.reference_end),
      QueryView(query, reverse_query, alignment.query_begin, alignment.query_end, alignment.strand),
      alignment.cigar, model);
  alignment.edit_distance = 0;
  for (const auto& op : alignment.cigar)
    if (op.operation != '=') alignment.edit_distance += op.length;
}
void StartSegment(const Path& path, const Anchor& anchor, const SequenceRecord& reference,
                  const SequenceRecord& query, CentromereStitchResult& out) {
  CentromereStitchSegment segment;
  segment.segment_id = out.segments.size() + 1;
  segment.path_id = path.id;
  segment.path_role = path.role;
  segment.region_id = path.region_id;
  segment.branch_id = path.branch_id;
  auto& alignment = segment.alignment;
  alignment.reference_id = reference.numeric_id;
  alignment.query_id = query.numeric_id;
  alignment.reference_begin = anchor.rb;
  alignment.reference_end = anchor.rb + anchor.length;
  alignment.query_begin = anchor.qb;
  alignment.query_end = anchor.qb + anchor.length;
  alignment.strand = anchor.strand;
  // The path role is authoritative; no SAM/PAF selection is performed here.
  alignment.primary = path.role == "primary";
  alignment.cigar.push_back({'=', anchor.length});
  segment.anchor_ids.push_back(anchor.id);
  out.segments.push_back(std::move(segment));
}
void ExtendSegment(CentromereStitchSegment& segment, const Anchor& anchor,
                   const CentromereStitchGap* gap) {
  if (gap) {
    for (const auto& op : gap->cigar) Append(segment.alignment.cigar, op.operation, op.length);
    segment.gap_ids.push_back(gap->gap_id);
  }
  Append(segment.alignment.cigar, '=', anchor.length);
  segment.anchor_ids.push_back(anchor.id);
  segment.alignment.reference_end = anchor.rb + anchor.length;
  segment.alignment.query_begin = std::min(segment.alignment.query_begin, anchor.qb);
  segment.alignment.query_end = std::max(segment.alignment.query_end, anchor.qb + anchor.length);
}
void CountEngine(CentromereStitchResult& out, std::string_view engine) {
  if (engine == "popoa_dp") ++out.dp_calls;
  if (engine == "ggwfa") ++out.wavefront_calls;
}
// IDs are deliberately excluded: replacing reviewed roots refreshes their IDs
// and chain indices, while physical endpoints and path ownership are stable.
using FixedGapKey = std::tuple<std::string, std::uint64_t, std::uint64_t, Strand,
                               Position, Position, Position, Position>;
FixedGapKey FixedKey(const CentromereStitchGap& gap) {
  return {gap.path_role, gap.region_id, gap.branch_id, gap.strand,
          gap.reference_begin, gap.reference_end, gap.query_begin, gap.query_end};
}
std::map<FixedGapKey, const CentromereStitchGap*> FixedGapIndex(
    const std::vector<CentromereStitchGap>& fixed,
    const std::vector<CentromereStitchGap>& actual) {
  if (fixed.empty()) return {};
  std::set<FixedGapKey> eligible;
  for (const auto& gap : actual) {
    if (!eligible.insert(FixedKey(gap)).second)
      throw AlignmentError("duplicate physical gap in frozen stitch paths");
  }
  std::map<FixedGapKey, const CentromereStitchGap*> index;
  for (const auto& gap : fixed) {
    const auto key = FixedKey(gap);
    if (!index.emplace(key, &gap).second)
      throw AlignmentError("duplicate reviewed stitch gap key");
    if (!eligible.count(key))
      throw AlignmentError("reviewed stitch gap does not match a current adjacent-anchor gap");
    if (!gap.scored || gap.status != "complete")
      throw AlignmentError("reviewed stitch gap must contain a completed scored CIGAR");
    for (const auto& op : gap.cigar) {
      if (!op.length || (op.operation != '=' && op.operation != 'I' && op.operation != 'D'))
        throw AlignmentError("reviewed stitch gap permits only positive =/I/D operations");
    }
  }
  return index;
}
void IdList(std::ostream& output, const std::vector<std::uint64_t>& ids) {
  if (ids.empty()) { output << '.'; return; }
  for (std::size_t i = 0; i < ids.size(); ++i) {
    if (i) output << ',';
    output << ids[i];
  }
}
void Cigar(std::ostream& output, const std::vector<CigarOp>& cigar) {
  if (cigar.empty()) { output << '*'; return; }
  for (const auto& op : cigar) output << op.length << op.operation;
}
void Field(std::ostream& output, std::string_view value) {
  // Reasons may originate in a host callback; keep the formal TSV rectangular.
  if (value.empty()) { output << '.'; return; }
  for (char c : value) output << ((c == '\t' || c == '\n' || c == '\r') ? ' ' : c);
}
}  // namespace

std::string_view CentromereStitchScoreModelName(CentromereStitchScoreModel model) {
  switch (model) {
    case CentromereStitchScoreModel::Centrolign3Piece: return "centrolign-3piece";
    case CentromereStitchScoreModel::Rama2PieceControl: return "rama-2piece-control";
  }
  throw AlignmentError("unknown centromere stitch score model");
}

CentromereStitchResult detail::AlignCentromereAnchoredPathsWithFixedGaps(
    const SequenceRecord& reference, const SequenceRecord& query,
    const CentromereDirectedAnchorResult& primary,
    const CentromereCopyAnchorResult* accepted_copies,
    const CentromereStitchOptions& options,
    const std::vector<CentromereStitchGap>& fixed_gaps) {
  if (!CentromereBackendAvailable())
    throw DependencyError("centromere anchored alignment requires the centromere backend");
  (void)CentromereStitchScoreModelName(options.score_model);
  if (!options.max_dp_cells || !options.gap_memory_limit_bytes ||
      !std::isfinite(options.gap_seconds_limit) || options.gap_seconds_limit <= 0 ||
      !std::isfinite(options.stage_seconds_limit) || options.stage_seconds_limit <= 0)
    throw AlignmentError("stitch resource limits must be positive and finite");
  const auto started = Clock::now();
  CentromereStitchResult out;
  out.score_model = options.score_model;
  out.reference_length = reference.size();
  out.query_length = query.size();
  bool in_observation = false;
  auto observe = [&] {
    Refresh(out, started);
    if (options.observation_callback) {
      in_observation = true;
      options.observation_callback(out);
      in_observation = false;
    }
  };
  try {
    ValidateSequence(reference);
    ValidateSequence(query);
    auto paths = MakePaths(reference, query, primary, accepted_copies, out);
    const auto fixed_index = FixedGapIndex(fixed_gaps, out.gaps);
    std::set<FixedGapKey> consumed_fixed;
    std::string reverse_query;
    if (std::any_of(paths.begin(), paths.end(), [](const Path& path) {
          return path.anchors.front().strand == Strand::Reverse;
        })) {
      reverse_query.reserve(query.bases.size());
      for (auto i = query.bases.rbegin(); i != query.bases.rend(); ++i)
        reverse_query.push_back(Complement(*i));
    }
    bool stage_exhausted = false;
    auto check = [&](std::string_view phase) {
      if (options.interruption_callback) options.interruption_callback(phase);
      if (stage_exhausted || Seconds(started) >= options.stage_seconds_limit)
        throw detail::StitchGapResourceLimit("stage_time_budget");
    };
    for (const auto& path : paths) {
      if (options.interruption_callback) options.interruption_callback("stitch-path");
      StartSegment(path, path.anchors.front(), reference, query, out);
      ValidateSegment(out.segments.back(), reference, query, reverse_query, options.score_model);
      for (std::size_t i = 1; i < path.anchors.size(); ++i) {
        const auto gap_id = path.gaps[i - 1];
        if (!gap_id) {
          ExtendSegment(out.segments.back(), path.anchors[i], nullptr);
          ValidateSegment(out.segments.back(), reference, query, reverse_query, options.score_model);
          continue;
        }
        auto& gap = out.gaps[gap_id - 1];
        const auto gap_started = Clock::now();
        try {
          check("stitch-gap");
          const auto r = ReferenceView(reference, gap.reference_begin, gap.reference_end);
          const auto q = QueryView(query, reverse_query, gap.query_begin, gap.query_end, gap.strand);
          const auto fixed = fixed_index.empty() ? fixed_index.end() : fixed_index.find(FixedKey(gap));
          if (fixed != fixed_index.end()) {
            const auto& reviewed = *fixed->second;
            const auto score = detail::ScoreCentromereStitchCigar(
                r, q, reviewed.cigar, options.score_model);
            if (score != reviewed.score)
              throw AlignmentError("reviewed stitch gap score disagrees with original FASTA replay");
            check("stitch-reviewed-gap");
            if (Seconds(gap_started) >= options.gap_seconds_limit)
              throw detail::StitchGapResourceLimit("gap_time_budget");
            gap.cigar = reviewed.cigar;
            gap.score = score;
            gap.engine = "reviewed_exact_gap_replay";
            consumed_fixed.insert(fixed->first);
          } else if (r.empty() || q.empty()) {
            gap.engine = "geometry";
            if (r.empty()) Append(gap.cigar, 'I', q.size());
            else Append(gap.cigar, 'D', r.size());
            gap.score = detail::ScoreCentromereStitchCigar(r, q, gap.cigar, options.score_model);
          } else {
            auto local = options;
            local.interruption_callback = check;
            local.progress_callback = {};
            local.observation_callback = {};
            auto aligned = detail::AlignCentromereStitchGap(r, q, local);
            gap.engine = std::move(aligned.engine);
            gap.score = aligned.score;
            gap.cigar = std::move(aligned.cigar);
            gap.dp_cells = aligned.dp_cells;
            gap.wfa_states = aligned.wfa_states;
            gap.peak_kernel_bytes = aligned.peak_kernel_bytes;
            CountEngine(out, gap.engine);
            if (detail::ScoreCentromereStitchCigar(r, q, gap.cigar, options.score_model) != gap.score)
              throw AlignmentError("stitch kernel score disagrees with independent CIGAR replay");
          }
          gap.scored = true;
          gap.status = "complete";
          ++out.completed_gap_count;
        } catch (const detail::StitchGapResourceLimit& error) {
          gap.status = "unresolved";
          gap.reason = error.reason;
          gap.engine = error.engine;
          gap.dp_cells = error.dp_cells;
          gap.wfa_states = error.wfa_states;
          gap.peak_kernel_bytes = error.peak_kernel_bytes;
          CountEngine(out, gap.engine);
          gap.cigar.clear();
          gap.scored = false;
          ++out.unresolved_gap_count;
          if (error.reason == "stage_time_budget") stage_exhausted = true;
        } catch (const std::bad_alloc&) {
          gap.status = "unresolved";
          gap.reason = "gap_allocation_failed";
          gap.cigar.clear();
          gap.scored = false;
          ++out.unresolved_gap_count;
        }
        gap.seconds = Seconds(gap_started);
        out.gap_seconds += gap.seconds;
        if (gap.status == "complete") ExtendSegment(out.segments.back(), path.anchors[i], &gap);
        else StartSegment(path, path.anchors[i], reference, query, out);
        ValidateSegment(out.segments.back(), reference, query, reverse_query, options.score_model);
        if (options.progress_callback)
          options.progress_callback("stitch-gaps", out.completed_gap_count + out.unresolved_gap_count, out.gaps.size());
        observe();
      }
      ValidateSegment(out.segments.back(), reference, query, reverse_query, options.score_model);
      observe();
    }
    // All keys were checked against current adjacency before traversal. A
    // resource refusal leaves a valid gap unresolved, rather than converting
    // the documented partial result into a stale-record integrity error.
    if (!out.unresolved_gap_count && consumed_fixed.size() != fixed_index.size())
      throw AlignmentError("unused reviewed stitch gap");
    out.complete = out.unresolved_gap_count == 0;
    observe();
    return out;
  } catch (...) {
    const auto original = std::current_exception();
    out.complete = false;
    out.aborted = true;
    // Failure metadata is best effort: allocation failure here must not mask
    // the caller's original interruption, bad allocation, or integrity error.
    try {
      try { std::rethrow_exception(original); }
      catch (const std::exception& error) { out.abort_reason = error.what(); }
      catch (...) { out.abort_reason = "non-standard exception"; }
    } catch (...) {}
    // A failing host observer must not recursively invoke itself or replace
    // an original interruption/error while exposing the previously kept data.
    if (!in_observation) {
      try { observe(); } catch (...) {}
    }
    std::rethrow_exception(original);
  }
}

CentromereStitchResult AlignCentromereAnchoredPaths(
    const SequenceRecord& reference, const SequenceRecord& query,
    const CentromereDirectedAnchorResult& primary,
    const CentromereCopyAnchorResult* accepted_copies,
    const CentromereStitchOptions& options) {
  return detail::AlignCentromereAnchoredPathsWithFixedGaps(
      reference, query, primary, accepted_copies, options, {});
}

void detail::AppendCentromereExactSingletons(
    const SequenceRecord& reference, const SequenceRecord& query,
    const CentromereDirectedAnchorResult& primary, const CentromereCopyAnchorResult* copies,
    std::span<const CentromereDirectedAnchor> added, CentromereStitchResult& out) {
  const auto started = Clock::now();
  const auto previous_seconds = out.total_seconds;
  CentromereStitchResult summary;
  (void)MakePaths(reference, query, primary, copies, summary);
  std::uint64_t path_id = 0;
  for (const auto& segment : out.segments) path_id = std::max(path_id, segment.path_id);
  std::string reverse_query;
  if (std::any_of(added.begin(), added.end(), [](const auto &a) { return a.strand == Strand::Reverse; })) {
    reverse_query.reserve(query.size());
    for (auto i = query.bases.rbegin(); i != query.bases.rend(); ++i)
      reverse_query.push_back(Complement(*i));
  }
  for (const auto& a : added) {
    if (path_id == std::numeric_limits<std::uint64_t>::max())
      throw AlignmentError("unplanned stitch path ID overflow");
    Path path;
    path.id = ++path_id;
    path.role = "primary";
    path.region_id = a.region_id;
    StartSegment(path, {a.id, a.chain_index, a.reference_begin, a.query_begin, a.length, a.strand},
                 reference, query, out);
    ValidateSegment(out.segments.back(), reference, query, reverse_query, out.score_model);
  }
  out.primary_anchor_count = summary.primary_anchor_count;
  out.primary_path_count = summary.primary_path_count;
  out.no_anchor_region_count = summary.no_anchor_region_count;
  out.singleton_path_count = summary.singleton_path_count;
  out.reference_support_bases = summary.reference_support_bases;
  out.query_support_bases = summary.query_support_bases;
  out.unknown_query_bases = summary.unknown_query_bases;
  out.known_query_outside_primary_support_bases = summary.known_query_outside_primary_support_bases;
  Refresh(out, started);
  out.total_seconds += previous_seconds;
}

void WriteCentromereStitchSegments(std::ostream& output, const CentromereStitchResult& result) {
  output << "segment_id\tpath_id\tpath_role\tregion_id\tbranch_id\treference_begin\treference_end"
            "\tquery_begin\tquery_end\tstrand\tcigar\tscore\tscore_model\tedit_distance\tanchor_ids\tgap_ids\n";
  for (const auto& segment : result.segments) {
    const auto& alignment = segment.alignment;
    output << segment.segment_id << '\t' << segment.path_id << '\t' << segment.path_role << '\t'
           << segment.region_id << '\t' << segment.branch_id << '\t'
           << alignment.reference_begin << '\t' << alignment.reference_end << '\t'
           << alignment.query_begin << '\t' << alignment.query_end << '\t' << Sign(alignment.strand) << '\t';
    Cigar(output, alignment.cigar);
    output << '\t' << alignment.score << '\t' << CentromereStitchScoreModelName(result.score_model)
           << '\t' << alignment.edit_distance << '\t';
    IdList(output, segment.anchor_ids);
    output << '\t';
    IdList(output, segment.gap_ids);
    output << '\n';
  }
  if (!output) throw AlignmentError("failed to write stitch segments");
}
void WriteCentromereStitchGaps(std::ostream& output, const CentromereStitchResult& result) {
  output << std::setprecision(std::numeric_limits<double>::max_digits10)
         << "gap_id\tpath_id\tpath_role\tregion_id\tbranch_id\tleft_anchor_id\tright_anchor_id"
            "\treference_begin\treference_end\tquery_begin\tquery_end\tstrand\tstatus\treason\tengine"
            "\tcigar\tscore\tscore_model\tdp_cells\twfa_states\tpeak_kernel_bytes\tseconds\n";
  for (const auto& gap : result.gaps) {
    output << gap.gap_id << '\t' << gap.path_id << '\t' << gap.path_role << '\t' << gap.region_id
           << '\t' << gap.branch_id << '\t' << gap.left_anchor_id << '\t' << gap.right_anchor_id
           << '\t' << gap.reference_begin << '\t' << gap.reference_end << '\t' << gap.query_begin
           << '\t' << gap.query_end << '\t' << Sign(gap.strand) << '\t' << gap.status << '\t';
    Field(output, gap.reason); output << '\t';
    Field(output, gap.engine); output << '\t';
    Cigar(output, gap.cigar); output << '\t';
    if (gap.scored) output << gap.score; else output << "NA";
    output << '\t' << CentromereStitchScoreModelName(result.score_model) << '\t' << gap.dp_cells
           << '\t' << gap.wfa_states << '\t' << gap.peak_kernel_bytes << '\t' << gap.seconds << '\n';
  }
  if (!output) throw AlignmentError("failed to write stitch gaps");
}
void WriteCentromereStitchSummary(std::ostream& output, const CentromereStitchResult& result) {
  output << std::setprecision(std::numeric_limits<double>::max_digits10)
         << "field\tvalue\nscore_model\t" << CentromereStitchScoreModelName(result.score_model)
         << "\nstatus\t" << (result.aborted ? "aborted" : result.complete ? "complete" : "partial")
         << "\ncomplete_scope\teligible_internal_gaps\n";
#define STITCH_FIELD(name) output << #name << '\t' << result.name << '\n'
  STITCH_FIELD(reference_length); STITCH_FIELD(query_length);
  STITCH_FIELD(primary_path_count); STITCH_FIELD(copy_path_count);
  STITCH_FIELD(primary_anchor_count); STITCH_FIELD(copy_anchor_count);
  STITCH_FIELD(no_anchor_region_count); STITCH_FIELD(singleton_path_count);
  STITCH_FIELD(completed_gap_count); STITCH_FIELD(unresolved_gap_count);
  STITCH_FIELD(dp_calls); STITCH_FIELD(wavefront_calls);
  STITCH_FIELD(reference_support_bases); STITCH_FIELD(query_support_bases);
  STITCH_FIELD(reference_aligned_bases); STITCH_FIELD(query_aligned_bases);
  STITCH_FIELD(reference_paired_bases); STITCH_FIELD(query_paired_bases);
  STITCH_FIELD(unknown_query_bases); STITCH_FIELD(known_query_outside_primary_support_bases);
  STITCH_FIELD(reference_unaligned_bases); STITCH_FIELD(query_unaligned_bases);
  STITCH_FIELD(total_seconds); STITCH_FIELD(gap_seconds);
  STITCH_FIELD(aborted);
  output << "abort_reason\t"; Field(output, result.abort_reason); output << '\n';
#undef STITCH_FIELD
  output << "segment_count\t" << result.segments.size() << "\ngap_count\t" << result.gaps.size() << '\n';
  if (!output) throw AlignmentError("failed to write stitch summary");
}
}  // namespace ramag

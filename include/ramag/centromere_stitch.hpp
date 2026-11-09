#pragma once

#include "ramag/centromere_copy_anchors.hpp"
#include <functional>
#include <iosfwd>
#include <string_view>

namespace ramag {

enum class CentromereStitchScoreModel { Centrolign3Piece, Rama2PieceControl };
struct CentromereStitchResult;

struct CentromereStitchOptions {
  CentromereStitchScoreModel score_model{CentromereStitchScoreModel::Centrolign3Piece};
  std::uint64_t max_dp_cells{4'000'000};
  std::uint64_t gap_memory_limit_bytes{512ULL * 1024ULL * 1024ULL};
  double gap_seconds_limit{120.0}, stage_seconds_limit{1800.0};
  // Gap kernels and the assembler propagate caller exceptions unchanged.
  std::function<void(std::string_view)> interruption_callback;
  std::function<void(std::string_view, std::uint64_t, std::uint64_t)> progress_callback;
  // Snapshot after each completed gap/path and immediately before propagating
  // a caller exception. Borrowed result is valid only during this callback.
  std::function<void(const CentromereStitchResult&)> observation_callback;
};

struct CentromereStitchSegment {
  std::uint64_t segment_id{}, path_id{}, region_id{}, branch_id{};
  std::string path_role;  // primary or copy; independent of SAM/PAF primary flag
  AlignmentRecord alignment;
  std::vector<std::uint64_t> anchor_ids, gap_ids;
};

struct CentromereStitchGap {
  std::uint64_t gap_id{}, path_id{}, region_id{}, branch_id{};
  std::string path_role;
  std::uint64_t left_anchor_id{}, right_anchor_id{};
  Position reference_begin{}, reference_end{}, query_begin{}, query_end{};
  Strand strand{Strand::Forward};
  std::string status, reason, engine;  // complete or unresolved; no heuristic success
  std::vector<CigarOp> cigar;
  std::int64_t score{};
  bool scored{};
  std::uint64_t dp_cells{}, wfa_states{}, peak_kernel_bytes{};
  double seconds{};
};

struct CentromereStitchResult {
  CentromereStitchScoreModel score_model{CentromereStitchScoreModel::Centrolign3Piece};
  std::vector<CentromereStitchSegment> segments;
  std::vector<CentromereStitchGap> gaps;
  Length reference_length{}, query_length{};
  std::uint64_t primary_path_count{}, copy_path_count{};
  std::uint64_t primary_anchor_count{}, copy_anchor_count{};
  std::uint64_t no_anchor_region_count{}, singleton_path_count{};
  std::uint64_t completed_gap_count{}, unresolved_gap_count{};
  std::uint64_t dp_calls{}, wavefront_calls{};
  // All coverage is a physical interval union; copy paths may share either axis.
  Length reference_support_bases{}, query_support_bases{};
  Length reference_aligned_bases{}, query_aligned_bases{};
  Length reference_paired_bases{}, query_paired_bases{};
  Length unknown_query_bases{}, known_query_outside_primary_support_bases{};
  Length reference_unaligned_bases{}, query_unaligned_bases{};
  double total_seconds{}, gap_seconds{};
  // Complete means every eligible internal gap succeeded, not whole-input
  // coverage. Tails, zero-anchor regions and ambiguous copies are not extended.
  bool complete{};
  // Diagnostic snapshots preserve completed segments on cancellation/hard error.
  bool aborted{};
  std::string abort_reason;
};

[[nodiscard]] std::string_view CentromereStitchScoreModelName(
    CentromereStitchScoreModel model);
// Pure consumer of frozen anchors. Paths are independent: reference/query
// overlap between paths is valid. Only consecutive anchors in one path bound
// a gap; no tail extension, new anchor search, or ambiguity resolution occurs.
// Coordinates stay on the original forward inputs. Negative CIGARs traverse
// reference increasingly and query decreasingly; I consumes query, D reference.
[[nodiscard]] CentromereStitchResult AlignCentromereAnchoredPaths(
    const SequenceRecord& reference, const SequenceRecord& query,
    const CentromereDirectedAnchorResult& primary,
    const CentromereCopyAnchorResult* accepted_copies = nullptr,
    const CentromereStitchOptions& options = {});

void WriteCentromereStitchSegments(std::ostream&, const CentromereStitchResult&);
void WriteCentromereStitchGaps(std::ostream&, const CentromereStitchResult&);
void WriteCentromereStitchSummary(std::ostream&, const CentromereStitchResult&);

}  // namespace ramag

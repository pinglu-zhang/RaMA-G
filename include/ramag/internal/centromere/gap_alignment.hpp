#pragma once
#include "ramag/centromere/stitch.hpp"
#include "ramag/alignment.hpp"
#include <span>
#include <string_view>

namespace ramag::detail {
struct StitchGapResult {
  std::vector<CigarOp> cigar;
  std::int64_t score{};
  std::string engine;
  std::uint64_t peak_kernel_bytes{}, dp_cells{}, wfa_states{};
  double elapsed_seconds{};
};
class StitchGapResourceLimit : public AlignmentError {
 public:
  explicit StitchGapResourceLimit(std::string why) : AlignmentError(why), reason(std::move(why)) {}
  std::string reason, engine;
  std::uint64_t peak_kernel_bytes{}, dp_cells{}, wfa_states{};
  double elapsed_seconds{};
};
[[nodiscard]] StitchGapResult AlignCentromereStitchGap(
    std::string_view reference, std::string_view oriented_query,
    const CentromereStitchOptions& options);
[[nodiscard]] std::int64_t ScoreCentromereStitchCigar(
    std::string_view reference, std::string_view oriented_query,
    std::span<const CigarOp> cigar, CentromereStitchScoreModel model);
}  // namespace ramag::detail


namespace ramag::detail {
// Private consumer for a previously reviewed exact gap path. Each record must
// match a current adjacent-anchor gap by path role/region/branch, strand and
// physical endpoints. IDs are intentionally ignored. All records must be
// complete and scored; =/I/D consumption and score are replayed against the
// original FASTA. Invalid, duplicated or stale records are hard errors.
[[nodiscard]] CentromereStitchResult AlignCentromereAnchoredPathsWithFixedGaps(
    const SequenceRecord& reference, const SequenceRecord& query,
    const CentromereDirectedAnchorResult& primary,
    const CentromereCopyAnchorResult* accepted_copies,
    const CentromereStitchOptions& options,
    const std::vector<CentromereStitchGap>& fixed_gaps);
// Add validated exact singleton paths to an existing assembly, preserving
// every existing segment and gap record, including its path identity.
void AppendCentromereExactSingletons(
    const SequenceRecord &, const SequenceRecord &, const CentromereDirectedAnchorResult &,
    const CentromereCopyAnchorResult *, std::span<const CentromereDirectedAnchor>,
    CentromereStitchResult &);

}  // namespace ramag::detail

namespace ramag::detail {
// Exact linear seven-state three-piece DP. Rolling score rows and nine
// traceback bits per cell are budgeted before allocation. This private
// candidate kernel is used by explicit endpoint and multi-anchor block review.
[[nodiscard]] StitchGapResult
AlignCentromereBlockDp(std::string_view reference,
                       std::string_view oriented_query,
                       const CentromereStitchOptions &options);
// Candidate-only constrained optimum under the same three-piece score.
// Only canonical exact diagonals and I/D gaps are permitted; no substitution
// diagonal, including N/N, can become part of its traceback. The unconstrained
// API above and its deterministic tie order remain unchanged.
[[nodiscard]] StitchGapResult
AlignCentromereBlockDpExactGaps(std::string_view reference,
                                std::string_view oriented_query,
                                const CentromereStitchOptions &options);
// Private forced-split entry for exact small-data differential checks. The
// production API selects its split only after the full trace exceeds budget.
// This entry still applies the same fixed 16 MiB safety reserve; no CLI knob.
[[nodiscard]] StitchGapResult AlignCentromereBlockDpTwoStrips(
    std::string_view reference, std::string_view oriented_query,
    const CentromereStitchOptions &options, std::uint64_t retained_suffix_rows);
// Forced storage split for independent exact-gap candidate tests only.
[[nodiscard]] StitchGapResult AlignCentromereBlockDpExactGapsTwoStrips(
    std::string_view reference, std::string_view oriented_query,
    const CentromereStitchOptions &options, std::uint64_t retained_suffix_rows);
} // namespace ramag::detail

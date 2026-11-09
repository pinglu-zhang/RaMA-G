#pragma once
#include "centromere_stitch_kernel.hpp"
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

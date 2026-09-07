#pragma once
#include "ramag/alignment.hpp"
#include <span>

namespace ramag {
// Immutable invocation configuration. Scoring is frozen pairwise HOXD70/10,
// open=40, extend=3; it never reads PAIRWISE_* environment variables.
struct PairwiseCoreOptions {
    Length max_gap{90};
    Length diag_diff{5};
    double diag_factor{0.12};
    Length min_cluster{65};
    std::uint32_t threads{1};
    std::function<void(std::string_view)> interruption_callback;
    std::function<void(std::string_view,std::uint64_t,std::uint64_t)> progress_callback;
};
struct PairwiseAlignment {
    AlignmentRecord record;
    bool reference_selected{};
    bool query_selected{};
    // Frozen pairwise selection accounting, not a calibrated identity estimate.
    Length selection_matching_columns{};
    Length selection_alignment_columns{};
};
struct PairwiseCoreResult {
    std::vector<PairwiseAlignment> records;
    std::uint64_t clusters{};
    std::uint32_t workers{};
    double clustering_seconds{};
    double extension_seconds{};
    double selection_seconds{};
};
// Graph-free, file-free, Sufkit-free entry. Inputs live until this call returns;
// returned records own their memory. Hooks may be called from worker threads.
// Coordinates are original-forward, zero-based, half-open. Seeds must be exact
// canonical A/C/G/T matches. Independent calls may run concurrently.
[[nodiscard]] PairwiseCoreResult AlignPairwiseCore(
    std::span<const SequenceRecord> references,
    std::span<const SequenceRecord> queries,
    std::span<const Seed> seeds,
    const PairwiseCoreOptions& options = {});

// CLI/Sufkit bridge: retains all or the pairwise reference/query DP intersection.
// No legacy reciprocal-interval resolver is applied to these records.
[[nodiscard]] AlignmentResult AlignPairwiseFromSeeds(
    const std::vector<SequenceRecord>& references,
    const std::vector<SequenceRecord>& queries,
    const AlignmentOptions& options,
    std::vector<Seed> seeds,
    RunStatistics seed_statistics = {});
} // namespace ramag

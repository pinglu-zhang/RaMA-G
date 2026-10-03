#pragma once
#include "ramag/alignment.hpp"
#include <span>

namespace ramag {
enum class GapFillMode { Off, ExactGap, Ksw2Gap, BudgetedKsw2GapV1 };
enum class PostSelectionPolicy { Baseline, OneSidedGuardedV1, FragmentReselectV1 };
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
    std::function<std::uint64_t()> resident_bytes_callback;
    // Experimental, sequence-aware post-selection recovery. Default off.
    // Returns original candidates plus recovered fragments carrying both flags.
    bool recover_uncovered_fragments{false};
    GapFillMode gap_fill{GapFillMode::Off};
    // Minimum inward exact support at each original flank (CLI: min_match).
    Length gap_fill_min_match{20};
    // Budgeted mode reserves up to three complete DP grids before dispatch.
    // Hosts may reduce this cap; actual budget is also bounded by input lengths.
    std::uint64_t gap_fill_work_budget{500000000};
    // Independently selectable coverage experiments; ordinary defaults remain off.
    bool preserve_link_candidates{false};
    bool link_precheck{false};
    bool recovery_require_flanks{false};
    // Experimental post-selection recovery or fragment-level re-selection.
    // Mutually exclusive with legacy recovery and gap filling.
    PostSelectionPolicy post_selection{PostSelectionPolicy::Baseline};
    // Independent coverage experiments, valid only with fragment re-selection.
    bool fragment_local_exchange{false};
    bool fragment_short_exact{false};
    bool supplement_pruned_clusters{false};
    bool supplement_alternative_chain{false};
    bool half_open_chaining{false};
    bool strand_aware_diagonal{false};
    bool strand_aware_merge{false};
    // Independent bounded experiments; thresholds do not change main/D chains.
    bool supplement_short_chain{false};
    bool extend_selected_endpoints{false};
    Length short_chain_min_support{50};
    Length bounded_output_min_support{65};
    Length bounded_extension_length{500};
    std::uint64_t bounded_work_budget{5000000000ULL};
    std::uint32_t bounded_workers{4};
    Length bounded_endpoint_min_exact{20};
};
struct PairwiseAlignment {
    AlignmentRecord record;
    bool reference_selected{};
    bool query_selected{};
    // Frozen pairwise selection accounting, not a calibrated identity estimate.
    Length selection_matching_columns{};
    Length selection_alignment_columns{};
    bool recovered_fragment{};
    bool gap_filled_fragment{};
    bool endpoint_extended_fragment{};
    // 0 main, 1 pre-cleanup chain, 2 alternative chain, 3 short chain, 4 endpoint.
    std::uint8_t candidate_source{};
};
struct PairwiseCoreResult {
    std::vector<PairwiseAlignment> records;
    std::uint64_t clusters{};
    std::uint32_t workers{};
    double clustering_seconds{};
    double extension_seconds{};
    double selection_seconds{};
    std::vector<MemoryObservation> memory_observations;
    PairwiseStatistics statistics;
};
// Resolve this binary's experimental build policy and explicit alignment options.
// CLI and library clients can use the same conversion; no process-global state.
[[nodiscard]] PairwiseCoreOptions ConfiguredPairwiseCoreOptions(const AlignmentOptions& options);
[[nodiscard]] std::string PairwiseConfigurationText(const AlignmentOptions& options);
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

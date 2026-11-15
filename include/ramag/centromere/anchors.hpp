#pragma once
#include "ramag/centromere.hpp"
#include <iosfwd>


namespace ramag {

// Frozen query-only direction evidence. IDs are independent of old block IDs.
struct CentromereDirectedRegion {
  std::uint64_t id{};
  Position query_begin{}, query_end{};
  Strand strand{Strand::Forward};
  std::size_t source_line{};
  friend bool operator==(const CentromereDirectedRegion&,
                         const CentromereDirectedRegion&) = default;
};
struct CentromereDirectedAnchor {
  std::uint64_t id{}, region_id{}, chain_index{};
  Position reference_begin{}, query_begin{};
  Length length{};
  Strand strand{Strand::Forward};
  double weight{};
  friend bool operator==(const CentromereDirectedAnchor&,
                         const CentromereDirectedAnchor&) = default;
};
struct CentromereDirectedTrace {
  std::uint64_t region_id{};
  CentromereRecursionNode node;
};
struct CentromereDirectedQueryInterval {
  std::uint64_t region_id{};
  Position query_begin{}, query_end{};
  Strand strand{Strand::Forward};
};
struct CentromereDirectedReferenceOverlap {
  std::uint64_t first_region_id{}, second_region_id{};
  Position reference_begin{}, reference_end{};
};
struct CentromereDirectedTask {
  CentromereDirectedRegion region;
  Length reference_length{};
  Position reference_begin{}, reference_end{};
  std::uint64_t initial_anchor_count{}, final_anchor_count{};
  Length initial_reference_bases{}, initial_query_bases{};
  Length final_reference_bases{}, final_query_bases{}, recursive_new_query_bases{};
  Length unanchored_query_bases{}, longest_unanchored_query_bases{};
  std::uint64_t unanchored_interval_count{}, recursion_node_count{};
  std::uint32_t max_recursion_depth{};
  double total_seconds{};
  std::string status;
  CentromereCoverageGuardStatistics coverage_guard,window_coverage_guard;
};
struct CentromereDirectedAnchorResult {
  std::vector<CentromereDirectedRegion> regions;
  std::vector<CentromereDirectedAnchor> initial, final;
  std::vector<CentromereDirectedAnchor> full_reference_initial;
  std::vector<CentromereDirectedTrace> trace;
  // Known direction without an exact anchor is not orientation-unknown.
  std::vector<CentromereDirectedQueryInterval> uncovered, unknown;
  std::vector<CentromereDirectedReferenceOverlap> reference_overlaps;
  std::vector<CentromereDirectedTask> tasks;
  double total_seconds{};
  std::uint64_t wfa_calls{};
  bool complete{};
};

// Read the final direction TSV (block_id/query_begin/query_end/strand), or
// this module's canonical five-column format. Legacy block IDs are metadata.
[[nodiscard]] std::vector<CentromereDirectedRegion> ReadCentromereDirectedRegions(
    const std::filesystem::path& path);
// Forward segments search the full reference. Reverse segments use an ordered
// reference window large enough for the query span when both adjacent forward
// segments provide ordered flanking anchors and the unconstrained initial chain
// escapes that window;
// otherwise they retain full-reference search. No truth is used to choose windows.
// Uses corrected half-open recursion. Never recomputes directions or calls WFA.
// No global one-to-one reference assignment is implied. Exceptions propagate;
// a partially completed run is never returned as successful.
[[nodiscard]] CentromereDirectedAnchorResult FindCentromereAnchorsInDirections(
    const SequenceRecord& reference, const SequenceRecord& query,
    std::span<const CentromereDirectedRegion> regions,
    const CentromereAlignmentOptions& options = {});

void WriteCentromereDirectedRegions(std::ostream&, const CentromereDirectedAnchorResult&);
void WriteCentromereDirectedInitialAnchors(std::ostream&, const CentromereDirectedAnchorResult&);
void WriteCentromereDirectedFullReferenceInitialAnchors(std::ostream&, const CentromereDirectedAnchorResult&);
void WriteCentromereDirectedFinalAnchors(std::ostream&, const CentromereDirectedAnchorResult&);
void WriteCentromereDirectedTrace(std::ostream&, const CentromereDirectedAnchorResult&);
void WriteCentromereDirectedUncovered(std::ostream&, const CentromereDirectedAnchorResult&);
void WriteCentromereDirectedReferenceOverlaps(std::ostream&, const CentromereDirectedAnchorResult&);
void WriteCentromereDirectedSummary(std::ostream&, const CentromereDirectedAnchorResult&);

} // namespace ramag


namespace ramag {
struct CentromereCopyAnchorOptions {
  std::uint32_t max_rare_level{100};
  std::uint64_t local_candidate_limit{3000}, total_candidate_limit{100000};
};
struct CentromereCopyAnchor {
  std::uint64_t id{}, branch_id{}, chain_index{};
  Position reference_begin{}, query_begin{};
  Length length{};
  Strand strand{Strand::Forward};
  double weight{};
  std::uint64_t task_id{}, self_repeat_id{}, primary_anchor_id{};
  std::string context_evidence;
  // Physical self-source in the direction traversed by increasing reference.
  Position self_source_oriented_begin{};
  Strand self_source_traversal{Strand::Forward};
};
struct CentromereCopyAmbiguity {
  std::uint64_t group_id{}, candidate_id{};
  CentromereCopyAnchor anchor;
  double score{}, score_margin{};
  std::string reason;
  double group_best_chain_score{}, runner_chain_score{};
  bool scored{};
  Position evaluation_target_begin{}, evaluation_target_end{};
  // Diagnostics only: one actual scored runner path, not all possible rivals.
  std::string evidence_role{"candidate"}, evidence_route{};
  std::uint64_t witness_order{}, witness_size{};
};
struct CentromereCopySelfRepeat {
  std::uint64_t id{}, task_id{};
  std::string side;
  Position first_begin{}, second_begin{};
  Length length{};
  Strand strand{Strand::Forward};
  std::uint64_t first_count{}, second_count{}, rare_level{};
  std::string count_scope{"local_lcp"};
};
struct CentromereCopyTask {
  std::uint64_t id{};
  std::string side;
  Position source_begin{}, source_end{}, target_begin{}, target_end{};
  std::uint64_t candidate_count{}, projected_count{}, accepted_count{},
      ambiguity_count{};
  std::uint64_t recursion_node_count{}, max_recursion_depth{};
  std::uint64_t branch_selection_count{}, recursive_added_count{};
  std::uint64_t retained_evidence_count{}, rival_candidate_count{};
  double recursive_seconds{};
  double self_seconds{}, projection_seconds{}, chain_seconds{}, total_seconds{};
  std::string status;
  Length unresolved_bases{};
};
struct CentromereCopyAnchorResult {
  std::vector<CentromereCopyAnchor> anchors;
  std::vector<CentromereCopyAmbiguity> ambiguities;
  std::vector<CentromereCopySelfRepeat> self_repeats;
  std::vector<CentromereCopyTask> tasks;
  std::uint64_t total_candidates{}, budget_exhausted_tasks{}, wfa_calls{};
  double total_seconds{};
  bool complete{};
};
// Sparse exact relations only. The primary is immutable and is the only
// cross-sequence projection root. Ambiguities never enter accepted output.
[[nodiscard]] CentromereCopyAnchorResult
FindCentromereCopyAnchors(const SequenceRecord &reference,
                          const SequenceRecord &query,
                          const CentromereDirectedAnchorResult &primary,
                          const CentromereCopyAnchorOptions &copy_options = {},
                          const CentromereAlignmentOptions &options = {});
void WriteCentromereCopyAnchors(std::ostream &,
                                const CentromereCopyAnchorResult &);
void WriteCentromereCopyAmbiguities(std::ostream &,
                                    const CentromereCopyAnchorResult &);
void WriteCentromereCopySelfRepeats(std::ostream &,
                                    const CentromereCopyAnchorResult &);
void WriteCentromereCopyTasks(std::ostream &,
                              const CentromereCopyAnchorResult &);
void WriteCentromereCopySummary(std::ostream &,
                                const CentromereCopyAnchorResult &);
} // namespace ramag

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

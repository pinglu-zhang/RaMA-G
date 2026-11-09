#pragma once

#include "ramag/model.hpp"
#include "ramag/alignment.hpp"

#include <filesystem>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

namespace ramag {

enum class CentromereMode { Genome, WholeInput, AnnotatedGenome };

// One recursive search node in the coordinates of the supplied, oriented input.
// IDs are deterministic preorder IDs (root 1, parent 0) emitted after workers drain.
struct CentromereRecursionNode {
  std::uint64_t node_id{}, parent_id{};
  std::uint32_t depth{}, rare_level{};
  Position reference_begin{}, reference_end{}, query_begin{}, query_end{};
  std::uint64_t candidate_count{}, selected_anchor_count{};
  std::string stop_reason;
  std::string guard_status{"disabled"};
  std::uint64_t query_added{}, reference_added{}, cache_bytes{};
};

struct CentromereAlignmentOptions {
  std::uint32_t worker_threads{1};
  std::function<void(std::string_view, std::uint64_t, std::uint64_t)> progress_callback;
  // May be called by computation workers; must be thread-safe.
  std::function<void(std::string_view)> interruption_callback;
  // Opt in to mathematically complete half-open recursive gap boundaries.
  // Legacy APIs preserve the pinned source boundary behavior by default.
  bool corrected_recursive_boundaries{false};
  std::function<void(const CentromereRecursionNode&)> recursion_observer;
};

[[nodiscard]] bool CentromereBackendAvailable() noexcept;
[[nodiscard]] bool CentromereCliAnchorsOnly() noexcept;
struct CentromereAnchor {
  Position reference_begin{}, query_begin{};
  Length length{};
  double weight{};
};
// Proposal timing is summed across workers; cache bytes are retained proposal
// capacity for this search, not process RSS or total solver workspace.
struct CentromereCoverageGuardStatistics {
  bool enabled{};
  double proposal_seconds_sum{},selection_seconds{},refill_seconds{};
  std::uint64_t baseline_nodes{},cached_nodes{},cache_exhausted_nodes{},proposal_budget_exhausted_nodes{};
  std::uint64_t accepted_nodes{},reused_subtrees{},refill_searches{},refill_nodes{},refill_nonempty_nodes{},cache_bytes_retained{};
};
struct CentromereAnchorResult {
  std::vector<CentromereAnchor> initial, final;
  RunStatistics statistics;
  CentromereCoverageGuardStatistics coverage_guard;
};
// Source initial/final anchors, not a gapped alignment. No WFA call.
[[nodiscard]] CentromereAnchorResult FindCentromereAnchors(
    const SequenceRecord& reference, const SequenceRecord& query,
    const CentromereAlignmentOptions& options = {});
// Exactly one nonempty sequence per side, supplied orientation. Ambiguous
// bases are N, never exact matches. This is not an inversion caller.
[[nodiscard]] AlignmentResult AlignCentromeres(
    const SequenceRecord& reference, const SequenceRecord& query,
    const CentromereAlignmentOptions& options = {});

struct CentromereSpec {
  bool whole_input{false};
  std::filesystem::path reference_bed;
  std::filesystem::path query_bed;

  [[nodiscard]] CentromereMode Mode() const noexcept;
};

struct CentromereRegion {
  std::string contig;
  Position begin{};
  Position end{};
  std::string name;
  std::size_t source_line{};
};

struct CentromereAnnotations {
  std::vector<CentromereRegion> reference;
  std::vector<CentromereRegion> query;
};

class CentromereInputError : public std::runtime_error {
 public:
  using std::runtime_error::runtime_error;
};

[[nodiscard]] std::string CentromereModeName(CentromereMode mode);
// BED order and optional names are metadata, never a pairing instruction.
[[nodiscard]] std::vector<CentromereRegion> ReadCentromereBed(
    const std::filesystem::path& path);
[[nodiscard]] CentromereAnnotations LoadCentromereAnnotations(
    const CentromereSpec& spec);
// Bounds use the same first-token FASTA identifiers as the genome reader.
void ValidateCentromereBounds(std::span<const CentromereRegion> regions,
                             std::span<const SequenceRecord> sequences);
// Refuse unsupported builds and BED mixed alignment before ordinary alignment.
void RequireCentromereBackend(const CentromereSpec& spec);

}  // namespace ramag

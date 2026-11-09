#pragma once
#include "ramag/centromere_directed_anchors.hpp"
#include <iosfwd>

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

#pragma once
#include "ramag/centromere.hpp"
#include <array>
#include <compare>
#include <iosfwd>
#include <limits>
#include <memory>


namespace ramag {

// Experimental evidence API. All query positions use the original forward
// sequence. No WFA, CIGAR, or confirmed biological event is produced.
struct DirectionBox {
  Position reference_begin{}, reference_end{}, query_begin{}, query_end{};
  friend auto operator<=>(const DirectionBox&, const DirectionBox&) = default;
};
struct DirectionOrigin {
  std::uint64_t task_id{};
  std::uint32_t depth{}, rare_level{}, reference_interval_count{}, query_interval_count{};
  bool first_hit_level{};
};
struct DirectionOriginNode {
  DirectionOrigin value;
  mutable std::shared_ptr<const DirectionOriginNode> previous;
  std::size_t count{};
  ~DirectionOriginNode() {
    // A long unshared history must not recursively destroy one stack frame
    // per origin during cancellation or an evidence-capacity failure.
    while(previous && previous.use_count()==1) {
      auto node=std::move(previous);previous=std::move(node->previous);
    }
  }
};
class DirectionOriginInterner;
// Immutable shared histories. Occurrence expansion often gives millions of
// anchors the same provenance; each anchor retains its complete ordered view.
class DirectionOrigins {
  friend class DirectionOriginInterner;
  std::shared_ptr<const DirectionOriginNode> tail_;
 public:
  class const_iterator {
    std::shared_ptr<std::vector<const DirectionOriginNode*>> path_;
    std::size_t position_{};
   public:
    const_iterator() = default;
    explicit const_iterator(const std::shared_ptr<const DirectionOriginNode>& tail) {
      if(!tail) return;
      path_=std::make_shared<std::vector<const DirectionOriginNode*>>();
      path_->reserve(tail->count);
      for(auto p=tail.get();p;p=p->previous.get()) path_->push_back(p);
    }
    const DirectionOrigin& operator*() const {return (*path_)[path_->size()-1-position_]->value;}
    const_iterator& operator++() {if(++position_==path_->size())path_.reset();return *this;}
    friend bool operator==(const const_iterator& a,const const_iterator& b) {
      if(!a.path_ || !b.path_)return !a.path_ && !b.path_;
      return (*a.path_)[a.path_->size()-1-a.position_]==(*b.path_)[b.path_->size()-1-b.position_];
    }
  };
  bool empty() const noexcept {return !tail_;}
  std::size_t size() const noexcept {return tail_ ? tail_->count : 0;}
  const DirectionOrigin& back() const {return tail_->value;}
  const_iterator begin() const {return const_iterator(tail_);}
  const_iterator end() const {return {};}
};
struct DirectionAnchor {
  std::uint64_t id{};
  Position reference_begin{}, query_begin{};
  Length length{};
  Strand strand{Strand::Forward};
  double source_weight{};
  Length full_reference_count{}, full_query_count{};
  DirectionOrigins origins;
};
inline constexpr std::size_t kNoDirectionNode = std::numeric_limits<std::size_t>::max();
struct DirectionChainNode {
  Position reference_begin{}, oriented_query_begin{};
  Length length{}, support{};
  double weight{}, source_score{};
  std::size_t predecessor{kNoDirectionNode};
  std::vector<std::uint64_t> source_ids;
};
struct DirectionChain {
  // All endpoint optima are retained as predecessor links, not copied paths.
  std::vector<DirectionChainNode> nodes;
  std::size_t best{kNoDirectionNode};
  Length support{};
};
struct DirectionRegion {
  std::uint64_t id{}, task_id{};
  std::uint32_t depth{};
  DirectionBox box;
  Length scale{};
  std::array<DirectionChain,2> chains;
  double advantage{}; // signed descriptive ratio, never a probability
};
struct DirectionTask {
  std::uint64_t id{}, parent_id{};
  std::uint32_t depth{};
  DirectionBox box;
  std::string stop_reason;
};
struct InversionEvidence {
  std::uint64_t id{}, region_id{};
  DirectionBox block;
  Strand strand{Strand::Reverse}, background{Strand::Forward};
  std::string classification, reason;
  Length support{}, opposite_support{}, context_background_support{};
  Length longest_match{};
  std::size_t independent_segments{};
  std::vector<std::uint64_t> anchors;
  std::vector<std::uint64_t> enclosing_opposite_candidates;
  std::uint64_t left_flank{}, right_flank{}; // zero means unavailable
  std::array<Position,2> reference_left{}, reference_right{}, query_left{}, query_right{};
  bool left_bounded{}, right_bounded{}, competing{};
};
enum class DirectionSearchStrategy {
  RamaFirstHitFocused,
  LegacyMultiscale
};
struct CentromereDirectionOptions : CentromereAlignmentOptions {
  DirectionSearchStrategy search_strategy{DirectionSearchStrategy::RamaFirstHitFocused};
  // Only the legacy comparator may supplement beyond the first hit tier.
  bool supplemental_levels{true};
  bool recurse{true};
  // A failure guard, not a pruning budget. Includes retained evidence capacity;
  // is not a claim that allocator/index/process RSS is limited to this value.
  std::uint64_t evidence_memory_limit_bytes{100000000000ULL};
  // Optional low-frequency accounting diagnostics. Does not affect selection.
  std::function<void(std::string_view)> memory_diagnostic_callback;
};
struct CentromereDirectionResult {
  Length reference_length{}, query_length{};
  std::vector<DirectionAnchor> anchors;
  std::vector<DirectionRegion> regions;
  std::vector<DirectionTask> tasks;
  std::vector<InversionEvidence> candidates;
  RunStatistics statistics;
};

// Pure geometry entry point, also used by exhaustive independent test oracles.
// Input coordinates and finite nonnegative weights are validated; the sequence
// entry point additionally validates every exact match against input bases.
[[nodiscard]] DirectionChain ScoreDirectionChain(
    std::span<const DirectionAnchor> anchors, const DirectionBox& region,
    Length reference_length, Length query_length, Strand direction,
    const std::function<void()>& interruption = {});

[[nodiscard]] CentromereDirectionResult AnalyzeCentromereDirections(
    const SequenceRecord& reference, const SequenceRecord& query,
    const CentromereDirectionOptions& options = {});
[[nodiscard]] bool CentromereCliDirections() noexcept;
[[nodiscard]] DirectionSearchStrategy CentromereCliDirectionStrategy() noexcept;
[[nodiscard]] std::string_view DirectionStrategyName(DirectionSearchStrategy) noexcept;

// Writers only serialize evidence; caller owns transactional file publication.
void WriteDirectionAnchors(std::ostream&, const CentromereDirectionResult&);
void WriteDirectionRegions(std::ostream&, const CentromereDirectionResult&);
void WriteDirectionChains(std::ostream&, const CentromereDirectionResult&);
void WriteDirectionTasks(std::ostream&, const CentromereDirectionResult&);
void WriteInversionEvidence(std::ostream&, const CentromereDirectionResult&);
} // namespace ramag


namespace ramag {
enum class CentromereBlockStrategy { InitialDirectionBlocksV1, JointLocalBlocksV2, LocalCoherentFrontierV3,
  HierarchicalDirectionRegionsV4, DirectionRegionsOnlyV5, FastLocalDirectionV6, CountClusterDirectionsV8,
  OrientationEvidenceV9 };
struct CentromereOrientationWindow {
  std::uint64_t id{};
  Position begin{},end{};
  std::uint64_t forward{},reverse{};
  Length forward_length{},reverse_length{};
  char label{'.'};
};
// Query-only evidence. No reference box exists until copy correspondence is resolved.
struct CentromereOrientationInterval {
  std::uint64_t id{};
  Position begin{},end{};
  Strand strand{Strand::Forward};
  bool inferred{};
  std::uint64_t window_id{}; // Single covering window; zero also denotes a continuous-window union.
};
struct CentromereOrientationSource {
  std::uint64_t interval_id{},root_id{}; // One-based references, no second anchor copy.
};
struct CentromereOrientationGap {
  Position begin{},end{};
  Strand strand{Strand::Forward};
  std::uint64_t window_id{};
  std::string reason;
  bool accepted{};
  std::uint64_t left_block_id{},right_block_id{}; // Zero: new exact support; otherwise an unchanged v6 block.
};
struct CentromereQueryInterval {Position begin{},end{};};
// Original accepted exact endpoints, never inferred or clipped boundaries.
struct CentromereUnifiedPivot {
  std::uint64_t id{},source_id{},block_id{};
  Position reference_begin{},reference_end{},query_begin{},query_end{};
  Strand strand{Strand::Forward};
  bool from_left{}; // Grow right from query_end; otherwise left from query_begin.
  std::uint8_t source_kind{}; // 0 planning anchor, 1 root exact, 2 sparse root exact.
};
struct CentromereUnifiedScan {
  std::uint64_t id{},pivot_id{};
  Position begin{},end{};
  Strand strand{Strand::Forward};
  std::string reason;
};
struct CentromereUnifiedInterval {
  std::uint64_t id{};
  Position begin{},end{};
  Strand strand{Strand::Forward};
  std::uint64_t scan_id{};
};
struct CentromereLocalGapReview {
  std::uint64_t id{},left_pivot{},right_pivot{};
  Position begin{},end{},corridor_begin{},corridor_end{};
  Strand strand{Strand::Forward};
  bool accepted{};
  std::string reason;
};
// The root remains in all upstream evidence. Verdict 0 only removes its hard
// veto for this review: 1 boundary-crossing, 2 known direction, 3 strong root,
// 4 locally coherent independent opposite neighbor all retain the veto.
struct CentromereOppositeReview {
  std::uint64_t review_id{},normalized_id{},neighbor_id{};
  Position reference_begin{},query_begin{};
  Length length{};
  double weight{};
  unsigned verdict{};
};
struct CentromereCollinearReview : CentromereLocalGapReview {
  Length support{};
  double chain_score{};
  std::uint64_t candidate_nodes{},connection_checks{},chain_nodes{};
  bool competitive{};
};
// Exact witness coordinates; do not imply a selected unique reference copy.
struct CentromereCollinearWitness {
  std::uint64_t review_id{},ordinal{},normalized_id{};
  Position reference_begin{},query_begin{};
  Length length{};
  double weight{};
};
// Sparse roots can support query direction without selecting a reference copy.
struct CentromereSparseOrientationInterval {
  std::uint64_t id{};
  Position begin{},end{},parent_begin{},parent_end{};
  Strand strand{Strand::Forward};
  bool inferred{};
};
struct CentromereSparseOrientationGap {
  Position begin{},end{},parent_begin{},parent_end{};
  Strand strand{Strand::Forward};
  bool accepted{};
  std::string reason;
};
// One bounded tip of an ORIGINAL rejected sparse edge. The source interval ID
// joins orientation_sparse_sources; root bounds are the complete exact-query
// union, never an occupied/inferred cut edge or a selected reference copy.
struct CentromereSparseGapEnd {
  std::uint64_t id{},source_interval_id{};
  Position begin{},end{},parent_begin{},parent_end{},edge_begin{},edge_end{};
  Position root_begin{},root_end{};
  Strand strand{Strand::Forward};
  bool from_left{},accepted{};
  std::string reason;
  // Enabled only for the independent meeting-tip policy. Support is accepted
  // real query union inside this original root, never inferred gap span.
  bool meeting_compared{};
  Length accepted_root_support{},competing_root_support{};
  Position scan_begin{},scan_end{};
  std::string scan_reason;
};
// Retracted query-only direction, not a new reference correspondence or alignment.
struct CentromereOrientationBlockAudit {
  std::uint64_t block_id{};
  Position begin{},end{},gap_begin{},gap_end{};
  Strand strand{Strand::Forward};
};
// Final one-pass withdrawal evidence; strand is the inherited old direction.
struct CentromereKnownCounterchainReview {
  std::uint64_t id{},block_id{},left_node{},right_node{};
  Position begin{},end{},gap_begin{},gap_end{};
  Strand strand{Strand::Forward};
  Length opposite_support{},own_support{};
  double opposite_score{},own_score{};
  std::uint64_t opposite_nodes{},own_nodes{},connection_checks{};
  bool retracted{};
  std::string reason;
};
// Final direction changes use existing original roots; upstream evidence is immutable.
struct CentromereLocalDirectionRepairReview {
  std::uint64_t id{},block_id{};
  Position begin{},end{};
  Strand strand{Strand::Forward},inherited_strand{Strand::Forward};
  Length support{},opposite_support{};
  double score{},opposite_score{};
  std::uint64_t nodes{},opposite_nodes{};
  bool accepted{};
  std::string reason;
};
struct CentromereLocalDirectionRepairChange {
  std::uint64_t review_id{};
  Position begin{},end{};
  char before_strand{'.'},after_strand{'.'};
  std::string reason;
};
struct CentromereOrientationKnownInterval {
  std::uint64_t block_id{};
  Position begin{},end{};
  Strand strand{Strand::Forward};
};
struct CentromereCountWindow {
  std::uint64_t id{};
  Position begin{},end{};
  std::uint64_t forward{},reverse{};
  char label{'.'}; // +, -, ? (competition), . (insufficient).
};
struct CentromereCountNode {
  CentromereAnchor anchor;
  Strand strand{Strand::Forward};
  std::uint64_t predecessor{},first{},count{};
  Length support{};
  double weight{};
  bool neutral{};
};
struct CentromereCountCandidate {
  std::uint64_t id{},first{},endpoint{},count{},forward{},reverse{};
  Position rb{},re{},qb{},qe{};
  Strand strand{Strand::Forward};
  Length support{};
  double weight{};
  bool eligible{},selected{};
};
struct CentromereCountGap {
  std::uint64_t region_id{},left_node{},right_node{},window_id{},forward{},reverse{};
  Position rb{},re{},qb{},qe{};
  Strand strand{Strand::Forward};
  bool accepted{};
  std::string reason;
};
struct CentromereSupportInterval {
  std::uint64_t region_id{};
  Position rb{},re{},qb{},qe{};
  Strand strand{Strand::Forward};
};
struct CentromerePlanningAnchor {
  CentromereAnchor anchor;
  Strand strand{Strand::Forward};
  std::uint32_t reference_count{},query_count{},rare_level{};
};
struct CentromereBlock {
  std::uint64_t id{};
  Position reference_begin{},reference_end{},query_begin{},query_end{};
  Strand strand{Strand::Forward};
  bool skeleton{},tentative{},inversion_candidate{};
  Length forward_support{},reverse_support{};
  std::string status{"planned"},reason;
  bool direction_known{true}; // Unknown rows own positions but are never predictions.
  std::vector<CentromereAnchor> planning_anchors; // Original-forward query positions.
  std::uint64_t region_id{};
};
struct CentromereDirectionRegion {
  std::uint64_t id{};
  Position reference_begin{},reference_end{},query_begin{},query_end{};
  Strand strand{Strand::Forward};
  Length support{};
  double secondary_score{};
  std::uint64_t competing_mappings{},block_count{};
  Length paired_columns{},insertions{},deletions{};
  // Exact terminal matches bound the region; they are not resolved breakpoints.
  Length left_boundary_width{},right_boundary_width{};
};
struct CentromereRefinementRegion {
  Position reference_begin{},reference_end{},query_begin{},query_end{};
  std::uint64_t id{},round{},new_mappings{};
  std::string status{"pending"};
};
struct CentromereLocalEvidence {
  std::uint64_t candidate_index{}; // Shared candidate store, no second anchor copy.
  std::uint64_t region_id{};
  std::uint32_t full_reference_count{},full_query_count{};
};
struct CentromereBlockAnchor {
  CentromereAnchor anchor;
  Strand strand{Strand::Forward};
  std::uint64_t block_id{};
};
struct CentromereUnresolved {
  CentromereBlock interval;
  bool unpaired{};
};
// Complete-chain query direction review. Local evidence retains local rarity.
struct CentromereRegionNode {
  std::uint64_t id{},root_id{},search_id{};
  CentromerePlanningAnchor source;
  bool neutral{};
};
struct CentromereRegionReview : CentromereCollinearReview {
  std::uint64_t search_id{};
};
struct CentromereRegionSearch {
  std::uint64_t id{},review_id{},left_node{},right_node{};
  Position reference_begin{},reference_end{},query_begin{},query_end{};
  Strand strand{Strand::Forward};
  std::uint64_t input_bases{},candidate_count{};
  bool searched{};
  std::string reason;
};
struct CentromereRegionLocalCandidate {
  std::uint64_t id{},search_id{};
  CentromerePlanningAnchor source;
};
// A frozen unknown query gap and its independently audited endpoint hypotheses.
// A direction verdict does not assert a unique reference-copy correspondence.
struct CentromereUnifiedGapDecision {
  std::uint64_t id{},root_review_id{},chosen_review_id{},hypothesis_count{};
  Position begin{},end{};
  Strand strand{Strand::Forward};
  bool accepted{};
  std::string reason;
};
struct CentromereUnifiedHypothesis {
  std::uint64_t id{},gap_id{},left_node{},right_node{},search_id{},review_id{};
  bool accepted{};
  std::string reason;
};
struct CentromereOverlapVariant {
  CentromereAnchor original,variant;
  Strand strand{Strand::Forward};
  Length trim_front{},trim_back{};
};
struct CentromereDirectionProtectionReview {
  std::uint64_t id{},block_id{}; // block_id belongs to the compact baseline view.
  Position begin{},end{};
  Strand strand{Strand::Forward};
  Length old_matched{},counter_matched{},candidate_same{},candidate_opposite{},candidate_unknown{},restored{};
  bool protected_interval{},region_changed{};
  std::string reason;
};
struct CentromereDirectionIslandAudit {
  std::uint64_t id{},block_id{},left_source_id{},right_source_id{},checks{};
  Position begin{},end{},reference_begin{},reference_end{},query_begin{},query_end{};
  Strand strand{Strand::Forward},background{Strand::Forward};
  Length reference_covered{},query_covered{},reference_unsupported{},query_unsupported{};
  double weight{},score{};
  bool eligible{},retracted{};
  std::string reason;
  std::string source_role{"candidate"}; // Baseline IDs never alias candidate block IDs.
};
struct CentromereDirectionIslandContextReview {
  std::uint64_t id{},block_id{},checks{},components{};
  Position original_begin{},original_end{},core_reference_begin{},core_reference_end{},core_query_begin{},core_query_end{};
  Length reference_explained{},query_explained{},reference_unsupported{},query_unsupported{};
  double old_score{},score{};
  bool accepted{};
  std::string reason;
};
struct CentromereDirectionPlanReview {
  std::uint64_t id{},block_id{},left_index{},source_id{},checks{};
  Position reference_begin{},reference_end{},query_begin{},query_end{};
  Strand strand{Strand::Forward};
  double weight{};
  bool accepted{};
  std::string reason;
};
struct CentromereBlockResult {
  bool direction_plan_review_applied{},direction_plan_review_budget_exhausted{};
  bool direction_island_context_budget_exhausted{};
  std::vector<CentromereDirectionIslandContextReview> direction_island_context_reviews;
  double direction_plan_review_seconds{};
  std::uint64_t direction_plan_review_checks{};
  std::vector<CentromereDirectionPlanReview> direction_plan_reviews;
  bool direction_local_protection_applied{};
  double direction_local_protection_baseline_seconds{},direction_local_protection_seconds{};
  std::uint64_t direction_local_protection_checks{},direction_local_protection_regions_changed{};
  std::vector<CentromereDirectionProtectionReview> direction_local_protection_reviews;
  std::vector<CentromereOrientationKnownInterval> direction_local_protection_intervals;
  // Frozen old targets for the post-protection island audit, not alignment output,
  // candidate roots, macro blocks, window votes or independent flank evidence.
  std::vector<CentromereBlock> direction_local_protection_baseline_blocks;
  bool direction_island_audit_applied{};
  double direction_island_audit_seconds{};
  std::uint64_t direction_island_checks{};
  std::vector<CentromereDirectionIslandAudit> direction_island_audits;
  std::uint64_t overlap_pairs_examined{},overlap_variants_added{};
  std::uint64_t overlap_nonuniform_nodes_skipped{};
  std::uint64_t overlap_local_support_rejections{}; // Proposal checks, not distinct events.
  bool overlap_budget_exhausted{};
  std::vector<CentromereOverlapVariant> overlap_variants;
  AlignmentResult alignment;
  std::vector<CentromereAnchor> initial_forward,initial_reverse; // Original-forward query coordinates.
  std::vector<CentromereBlock> blocks;
  std::vector<CentromereBlockAnchor> refined;
  std::vector<CentromereUnresolved> unresolved;
  Strand background{Strand::Forward};
  bool background_tied{},complete{true};
  std::uint64_t wfa_calls{};
  double initial_seconds{},planning_seconds{},recursion_seconds{},wfa_seconds{};
  CentromereBlockStrategy strategy{CentromereBlockStrategy::InitialDirectionBlocksV1};
  std::vector<CentromerePlanningAnchor> root_candidates;
  std::uint64_t planner_reserved_bytes{},planner_capacity_bytes{};
  std::uint64_t local_chain_states{},local_chain_transitions{},local_chain_dominated{};
  std::vector<CentromereDirectionRegion> regions;
  std::vector<CentromereRefinementRegion> refinement_regions;
  std::vector<CentromereLocalEvidence> local_candidates;
  std::uint64_t refinement_rounds{},refinement_searches{},refinement_new_mappings{};
  std::uint64_t root_candidate_count{};
  std::string refinement_partial_reason;
  std::vector<CentromereCountWindow> count_windows;
  std::vector<CentromereCountNode> count_nodes;
  std::vector<CentromereCountCandidate> count_candidates;
  std::vector<CentromereCountGap> count_gaps;
  std::vector<CentromereSupportInterval> supported_intervals;
  std::uint64_t count_neighbor_checks{},count_rectangle_queries{},count_rectangle_nodes{};
  double count_prepare_seconds{},count_selection_seconds{};
  std::vector<CentromereOrientationWindow> orientation_windows;
  std::vector<CentromereOrientationInterval> orientation_intervals;
  std::vector<CentromereOrientationSource> orientation_sources;
  std::vector<CentromereOrientationGap> orientation_gaps;
  std::vector<CentromereQueryInterval> orientation_unknown;
  double orientation_seconds{};
  std::uint64_t orientation_capacity_bytes{};
  bool orientation_connected_gaps{};
  bool orientation_bridge_unknown_gaps{};
  bool orientation_audit_known_blocks{};
  double orientation_audit_seconds{};
  std::vector<CentromereOrientationBlockAudit> orientation_block_audits;
  // When audit is enabled, replaces raw block envelopes in query direction predictions.
  std::vector<CentromereOrientationKnownInterval> orientation_known_intervals;
  // Independent query-only overlay. Old evidence IDs and audit outputs are unchanged.
  bool recover_long_gap_ends{};
  double orientation_long_gap_seconds{};
  std::vector<CentromereOrientationInterval> orientation_long_gap_intervals;
  std::vector<CentromereOrientationGap> orientation_long_gap_gaps;
  bool recover_sparse_pivots{};
  double orientation_sparse_seconds{};
  std::vector<CentromereSparseOrientationInterval> orientation_sparse_intervals;
  std::vector<CentromereSparseOrientationGap> orientation_sparse_gaps;
  std::vector<CentromereOrientationSource> orientation_sparse_sources;
  bool recover_short_gap_ends{};
  double orientation_short_gap_seconds{};
  std::vector<CentromereOrientationInterval> orientation_short_gap_intervals;
  std::vector<CentromereOrientationGap> orientation_short_gap_gaps;
  bool recover_mixed_gap_ends{};
  double orientation_mixed_gap_seconds{};
  std::vector<CentromereOrientationInterval> orientation_mixed_gap_intervals;
  std::vector<CentromereOrientationGap> orientation_mixed_gap_gaps;
  bool recover_sparse_gap_ends{};
  bool prefer_stronger_sparse_gap_end{};
  double orientation_sparse_gap_seconds{};
  std::vector<CentromereSparseGapEnd> orientation_sparse_gap_intervals;
  std::vector<CentromereSparseGapEnd> orientation_sparse_gap_gaps;
  bool unify_real_anchor_ends{};
  double orientation_unified_seconds{};
  std::vector<CentromereUnifiedPivot> orientation_unified_pivots;
  std::vector<CentromereUnifiedScan> orientation_unified_proposals;
  std::vector<CentromereUnifiedInterval> orientation_unified_intervals;
  bool review_local_competitive_gaps{};
  double orientation_local_review_seconds{};
  std::vector<CentromereLocalGapReview> orientation_local_reviews;
  bool review_isolated_opposite_anchors{};
  double orientation_isolated_review_seconds{};
  std::vector<CentromereLocalGapReview> orientation_isolated_reviews;
  std::vector<CentromereOppositeReview> orientation_opposite_reviews;
  bool review_collinear_fragment_connections{};
  bool review_collinear_window_gaps{};
  double orientation_collinear_review_seconds{};
  std::vector<CentromereCollinearReview> orientation_collinear_reviews;
  std::vector<CentromereOppositeReview> orientation_collinear_counters;
  std::vector<CentromereCollinearWitness> orientation_collinear_witnesses;
  bool review_alternate_original_endpoints{};
  bool review_neutral_window_gaps{};
  bool review_region_chains{},supplement_region_anchors{};
  double region_chain_seconds{},region_search_seconds{};
  std::uint64_t region_search_calls{},region_search_input_bases{},region_connection_checks{},region_capacity_bytes{};
  std::vector<CentromereRegionNode> region_nodes;
  std::vector<CentromereRegionReview> region_reviews;
  std::vector<CentromereCollinearWitness> region_witnesses;
  std::vector<CentromereUnifiedInterval> region_intervals;
  std::vector<CentromereRegionSearch> region_searches;
  std::vector<CentromereRegionLocalCandidate> region_local_candidates;
  bool repair_inherited_directions{},repair_drift_consistent_gaps{};
  double local_direction_repair_seconds{};
  std::uint64_t local_direction_repair_connection_checks{},local_direction_repair_capacity_bytes{};
  std::vector<CentromereLocalDirectionRepairReview> local_direction_repair_reviews;
  std::vector<CentromereCollinearWitness> local_direction_repair_witnesses;
  std::vector<CentromereLocalDirectionRepairChange> local_direction_repair_changes;
  // Complete final query direction view; old overlays remain provenance.
  std::vector<CentromereOrientationKnownInterval> local_direction_repair_intervals;
  bool review_known_block_counterchains{};
  double known_counterchain_seconds{};
  std::uint64_t known_counterchain_connection_checks{},known_counterchain_capacity_bytes{};
  std::vector<CentromereKnownCounterchainReview> known_counterchain_reviews;
  std::vector<CentromereCollinearWitness> known_counterchain_witnesses;
  // Final retained old-block pieces; the original audited pieces stay immutable.
  std::vector<CentromereOrientationKnownInterval> known_counterchain_intervals;
  bool review_unified_gaps{},review_gap_interior_opposition{};
  double unified_review_seconds{},unified_search_seconds{};
  std::uint64_t unified_search_calls{},unified_search_input_bases{},unified_connection_checks{},unified_capacity_bytes{};
  std::vector<CentromereRegionNode> unified_nodes;
  std::vector<CentromereRegionReview> unified_reviews;
  std::vector<CentromereCollinearWitness> unified_witnesses;
  std::vector<CentromereUnifiedInterval> unified_intervals;
  std::vector<CentromereRegionSearch> unified_searches;
  std::vector<CentromereRegionLocalCandidate> unified_local_candidates;
  std::vector<CentromereUnifiedGapDecision> unified_gap_decisions;
  std::vector<CentromereUnifiedHypothesis> unified_hypotheses;
  std::vector<CentromereCollinearReview> orientation_alternate_reviews;
  std::vector<CentromereOppositeReview> orientation_alternate_counters;
  std::vector<CentromereCollinearWitness> orientation_alternate_witnesses;
};
struct CentromereBlockOptions : CentromereAlignmentOptions {
  // Opt-in candidate repairs; disabled preserves the frozen direction route.
  bool direction_plan_review{false}; // Frozen-root local competition; no additional searches.
  bool overlap_compatibility{false};
  bool direction_evidence_repair{false};
  bool direction_island_audit{false};
  bool direction_island_context_review{false}; // Credit only proven disjoint components in the same flank frame.
  std::uint64_t direction_island_context_check_limit{10000000}; // Also bounded by remaining original batch checks.
  bool direction_local_protection{false}; // FindCentromereDirectionRegions only; compact frozen-root baseline.
  CentromereBlockStrategy planning_strategy{CentromereBlockStrategy::InitialDirectionBlocksV1};
  std::uint64_t planner_memory_limit_bytes{10000000000ULL};
  std::uint64_t wfa_memory_limit_bytes{512ULL*1024*1024};
  double wfa_seconds_limit{120},alignment_seconds_limit{1800};
  Length count_window{10000},count_step{5000},count_cluster_gap{10000},count_inferred_gap{20000};
  std::uint32_t count_min_anchors{3},count_dominance_ratio{2};
  // Experimental connection-only policy. Does not change v6 blocks or window votes.
  bool orientation_connected_gaps{false};
  // Only bridges insufficient coverage, never opposing/competitive evidence.
  bool orientation_bridge_unknown_gaps{false};
  // Retract pure opposing-window evidence only inside long, anchor-free old block gaps.
  bool orientation_audit_known_blocks{false};
  // Recover only frozen long-gap ends with contiguous pure direction evidence.
  bool recover_long_gap_ends{false};
  // Existing exact roots only; inferred output never becomes a new pivot.
  bool recover_sparse_pivots{false};
  // Continuous-window terminal recovery in frozen <=limit gaps; never bridges a veto.
  bool recover_short_gap_ends{false};
  // Opposite original endpoints recover independently; the center stays unknown.
  bool recover_mixed_gap_ends{false};
  // Frozen accepted exact sparse roots only; '.' windows may be traversed once.
  bool recover_sparse_gap_ends{false};
  // Meeting tips compete by frozen accepted real root support; ties stay unknown.
  bool prefer_stronger_sparse_gap_end{false};
  // Frozen real endpoints share one bounded scan rule; same-direction tips may meet.
  bool unify_real_anchor_ends{false};
  // Same-direction frozen endpoints only. Competition is reviewed in the actual
  // connection corridor; opposite windows and in-corridor evidence still veto.
  bool review_local_competitive_gaps{false};
  bool review_isolated_opposite_anchors{false};
  bool review_collinear_fragment_connections{false};
  // Permit an incomplete own-window cover to reach the unchanged full audit.
  bool review_collinear_window_gaps{false};
  // Review bounded original endpoint alternatives only after the old pair fails.
  bool review_alternate_original_endpoints{false};
  // Zero own-window votes may reach full audit only without competitive votes.
  bool review_neutral_window_gaps{false};
  // One frozen region review and at most one search wave; never runs WFA.
  bool review_region_chains{false},supplement_region_anchors{false};
  Length region_search_axis_limit{100000};
  std::uint64_t region_search_input_budget{1000000};
  bool review_unified_gaps{false};
  // Only evidence supporting the actual unknown gap may veto its recovery.
  // Preserve full anchor lengths for strong evidence; never crop a real match.
  bool review_gap_interior_opposition{false};
  // Existing root counterchains may retract old inherited directions once, last.
  bool review_known_block_counterchains{false};
  // One final local competition over frozen original roots; no new searches.
  bool repair_inherited_directions{false};
  // Optional unique-root long edges with bounded two-axis drift; v1 is unchanged.
  bool repair_drift_consistent_gaps{false};
  Length unified_search_query_limit{100000},unified_search_reference_limit{512000};
  std::uint64_t unified_search_input_budget{4000000};
  // Synchronous, read-only evidence snapshots; never influences the plan.
  std::function<void(const CentromereBlockResult&)> observation_callback;
};

void WriteCentromereDirectionIslandContextReviews(std::ostream&,const CentromereBlockResult&);
void WriteCentromereDirectionPlanReviews(std::ostream&,const CentromereBlockResult&);
void WriteCentromereDirectionLocalProtectionReviews(std::ostream&,const CentromereBlockResult&);
void WriteCentromereDirectionIslandAudits(std::ostream&,const CentromereBlockResult&);
void WriteCentromereOverlapVariants(std::ostream&,const CentromereBlockResult&);
void WriteCentromereOverlapSummary(std::ostream&,const CentromereBlockResult&);
[[nodiscard]] bool CentromereCliBlocks() noexcept;
[[nodiscard]] bool CentromereCliFullAlignment() noexcept;
[[nodiscard]] CentromereBlockOptions CentromereCliBlockOptions();
[[nodiscard]] std::string_view CentromereBlockStrategyName(CentromereBlockStrategy) noexcept;
[[nodiscard]] std::vector<CentromereAnchor> FindInitialCentromereAnchors(
    const SequenceRecord&,const SequenceRecord&,const CentromereAlignmentOptions& = {});
// Pure planner. Inputs must each be a validated, nonoverlapping source chain.
[[nodiscard]] CentromereBlockResult PlanCentromereBlocks(
    Length reference_length,Length query_length,
    std::vector<CentromereAnchor> forward,std::vector<CentromereAnchor> reverse);
// Sparse block competition. Candidate discovery remains the original first-hit RaMA.
[[nodiscard]] CentromereBlockResult PlanJointCentromereBlocks(
    Length reference_length,Length query_length,
    std::vector<CentromerePlanningAnchor> candidates,
    const CentromereBlockOptions& options = {},Strand preferred_background = Strand::Forward);
[[nodiscard]] CentromereBlockResult AlignCentromeresWithDirections(
    const SequenceRecord&,const SequenceRecord&,const CentromereBlockOptions& = {});
// Direction planning only: no gap recursion, WFA, CIGAR or alignment records.
// Unknown intervals stay explicit. The caller selects the experimental strategy.
[[nodiscard]] CentromereBlockResult FindCentromereDirectionRegions(
    const SequenceRecord&,const SequenceRecord&,const CentromereBlockOptions&);
void WriteCentromereBlocks(std::ostream&,const CentromereBlockResult&);
void WriteCentromereInitialAnchors(std::ostream&,const CentromereBlockResult&);
void WriteCentromereRefinedAnchors(std::ostream&,const CentromereBlockResult&);
void WriteCentromereUnresolved(std::ostream&,const CentromereBlockResult&);
void WriteCentromereRootCandidates(std::ostream&,const CentromereBlockResult&);
void WriteCentromerePlanningAnchors(std::ostream&,const CentromereBlockResult&);
void WriteCentromereDirectionRegions(std::ostream&,const CentromereBlockResult&);
void WriteCentromereRefinementRegions(std::ostream&,const CentromereBlockResult&);
void WriteCentromereLocalCandidates(std::ostream&,const CentromereBlockResult&);
void WriteCentromereCountWindows(std::ostream&,const CentromereBlockResult&);
void WriteCentromereCountNodes(std::ostream&,const CentromereBlockResult&);
void WriteCentromereCountCandidates(std::ostream&,const CentromereBlockResult&);
void WriteCentromereCountGaps(std::ostream&,const CentromereBlockResult&);
void WriteCentromereSupportedIntervals(std::ostream&,const CentromereBlockResult&);
void WriteCentromereOrientationWindows(std::ostream&,const CentromereBlockResult&);
void WriteCentromereOrientationIntervals(std::ostream&,const CentromereBlockResult&);
void WriteCentromereOrientationSources(std::ostream&,const CentromereBlockResult&);
void WriteCentromereOrientationGaps(std::ostream&,const CentromereBlockResult&);
void WriteCentromereOrientationUnknown(std::ostream&,const CentromereBlockResult&);
void WriteCentromereOrientationBlockAudit(std::ostream&,const CentromereBlockResult&);
void WriteCentromereOrientationKnownIntervals(std::ostream&,const CentromereBlockResult&);
void WriteCentromereOrientationLongGapEnds(std::ostream&,const CentromereBlockResult&);
void WriteCentromereOrientationLongGapGaps(std::ostream&,const CentromereBlockResult&);
void WriteCentromereOrientationSparsePivots(std::ostream&,const CentromereBlockResult&);
void WriteCentromereOrientationSparseGaps(std::ostream&,const CentromereBlockResult&);
void WriteCentromereOrientationSparseSources(std::ostream&,const CentromereBlockResult&);
void WriteCentromereOrientationShortGapEnds(std::ostream&,const CentromereBlockResult&);
void WriteCentromereOrientationShortGapGaps(std::ostream&,const CentromereBlockResult&);
void WriteCentromereOrientationMixedGapEnds(std::ostream&,const CentromereBlockResult&);
void WriteCentromereOrientationMixedGapGaps(std::ostream&,const CentromereBlockResult&);
void WriteCentromereOrientationSparseGapEnds(std::ostream&,const CentromereBlockResult&);
void WriteCentromereOrientationSparseGapGaps(std::ostream&,const CentromereBlockResult&);
void WriteCentromereOrientationUnifiedPivots(std::ostream&,const CentromereBlockResult&);
void WriteCentromereOrientationUnifiedScans(std::ostream&,const CentromereBlockResult&);
void WriteCentromereOrientationUnifiedEnds(std::ostream&,const CentromereBlockResult&);
void WriteCentromereOrientationLocalReviews(std::ostream&,const CentromereBlockResult&);
void WriteCentromereOrientationIsolatedReviews(std::ostream&,const CentromereBlockResult&);
void WriteCentromereOrientationOppositeReviews(std::ostream&,const CentromereBlockResult&);
void WriteCentromereOrientationIsolatedGaps(std::ostream&,const CentromereBlockResult&);
void WriteCentromereOrientationCollinearReviews(std::ostream&,const CentromereBlockResult&);
void WriteCentromereOrientationCollinearCounters(std::ostream&,const CentromereBlockResult&);
void WriteCentromereOrientationCollinearWitnesses(std::ostream&,const CentromereBlockResult&);
void WriteCentromereOrientationCollinearGaps(std::ostream&,const CentromereBlockResult&);
void WriteCentromereOrientationAlternateReviews(std::ostream&,const CentromereBlockResult&);
void WriteCentromereOrientationAlternateCounters(std::ostream&,const CentromereBlockResult&);
void WriteCentromereOrientationAlternateWitnesses(std::ostream&,const CentromereBlockResult&);
void WriteCentromereOrientationAlternateGaps(std::ostream&,const CentromereBlockResult&);
void WriteCentromereRegionReviews(std::ostream&,const CentromereBlockResult&);
void WriteCentromereRegionNodes(std::ostream&,const CentromereBlockResult&);
void WriteCentromereRegionWitnesses(std::ostream&,const CentromereBlockResult&);
void WriteCentromereRegionGaps(std::ostream&,const CentromereBlockResult&);
void WriteCentromereRegionSearches(std::ostream&,const CentromereBlockResult&);
void WriteCentromereRegionLocalCandidates(std::ostream&,const CentromereBlockResult&);
void WriteCentromereUnifiedReviewReviews(std::ostream&,const CentromereBlockResult&);
void WriteCentromereUnifiedReviewNodes(std::ostream&,const CentromereBlockResult&);
void WriteCentromereUnifiedReviewWitnesses(std::ostream&,const CentromereBlockResult&);
void WriteCentromereUnifiedReviewIntervals(std::ostream&,const CentromereBlockResult&);
void WriteCentromereUnifiedReviewSearches(std::ostream&,const CentromereBlockResult&);
void WriteCentromereUnifiedReviewLocalCandidates(std::ostream&,const CentromereBlockResult&);
void WriteCentromereUnifiedReviewGaps(std::ostream&,const CentromereBlockResult&);
void WriteCentromereUnifiedReviewHypotheses(std::ostream&,const CentromereBlockResult&);
void WriteCentromereOrientationLocalGaps(std::ostream&,const CentromereBlockResult&);
void WriteCentromereKnownCounterchainReviews(std::ostream&,const CentromereBlockResult&);
void WriteCentromereKnownCounterchainWitnesses(std::ostream&,const CentromereBlockResult&);
void WriteCentromereKnownCounterchainIntervals(std::ostream&,const CentromereBlockResult&);
void WriteCentromereLocalDirectionRepairReviews(std::ostream&,const CentromereBlockResult&);
void WriteCentromereLocalDirectionRepairWitnesses(std::ostream&,const CentromereBlockResult&);
void WriteCentromereLocalDirectionRepairIntervals(std::ostream&,const CentromereBlockResult&);
void WriteCentromereLocalDirectionRepairChanges(std::ostream&,const CentromereBlockResult&);
} // namespace ramag

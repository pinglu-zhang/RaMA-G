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

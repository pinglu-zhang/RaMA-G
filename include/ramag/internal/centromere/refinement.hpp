#pragma once
#include "ramag/centromere/stitch.hpp"
#include <optional>
#include <chrono>
#include <set>
#include <unordered_map>

namespace ramag {
// An explicit post-search, pre-stitch stage. Raw search results remain
// auditable.
struct CentromereEndpointOptions {
  CentromereStitchOptions kernel;
  std::uint64_t max_reviews{64};
  double stage_seconds_limit{60.0}, window_seconds_limit{10.0};
  // Optional reference-unique, bidirectional seed discovery. Returned
  // coordinates are physical, with query coordinates relative to the window.
  std::function<std::optional<std::vector<Seed>>(Position, Position,
      const std::function<void()> &)> unplanned_candidates;
};
struct CentromereEndpointReview {
  std::uint64_t id{}, region_id{}, left_anchor_id{}, anchor_id{},
      right_anchor_id{};
  Position reference_begin{}, reference_end{}, query_begin{}, query_end{};
  Length left_gap{}, right_gap{}, old_exact_bases{}, new_exact_bases{},
      new_mismatch_bases{};
  std::int64_t old_score{}, new_score{};
  Length old_unique_context_starts{}, new_unique_context_starts{};
  bool scored{}, accepted{};
  std::string reason, engine;
  std::vector<CigarOp> old_cigar, new_cigar;
  std::uint64_t peak_kernel_bytes{}, dp_cells{}, wfa_states{};
  double seconds{};
};
struct CentromereEndpointResult {
  CentromereDirectedAnchorResult anchors;
  std::vector<CentromereEndpointReview> reviews;
  // Cross-region reviews return the actual validated assembly. Consumers must
  // use it: a new unconstrained solve could change the reviewed exact edits.
  std::optional<CentromereStitchResult> stitched;
  std::uint64_t accepted{}, kernel_calls{};
  double seconds{};
};
// Reviews primary roots in either strand between compensating one-sided gaps.
// Geometry uses oriented query coordinates; outputs retain physical coordinates.
// Tails and accepted copy footprints stay unchanged. Existing forward work
// retains priority under the shared resource and footprint limits.
// Review CIGARs advance reference and follow the region strand on query;
// review coordinates always describe original FASTA half-open intervals.
// Only strict score AND exact-coverage improvements without new substitutions
// are accepted. Identical windows retain the old copy choice. No truth input.
[[nodiscard]] CentromereEndpointResult
RefineCentromereRootEndpoints(const SequenceRecord &, const SequenceRecord &,
                              const CentromereDirectedAnchorResult &,
                              const CentromereCopyAnchorResult * = nullptr,
                              const CentromereEndpointOptions & = {});
// One non-iterative review of bounded multi-anchor blocks after the original
// stitch. Outer anchors require two nonoverlapping 32-mers unique in both
// input strings. The actual stitched CIGAR is the old-score baseline; missing
// or unresolved paths, ambiguous contexts and equal scores keep the input.
[[nodiscard]] CentromereEndpointResult RefineCentromereAnchorBlocks(
    const SequenceRecord &, const SequenceRecord &,
    const CentromereDirectedAnchorResult &, const CentromereCopyAnchorResult *,
    const CentromereStitchResult &, const CentromereEndpointOptions & = {});
// Conservative forward cross-region review. Existing regions and every reverse
// anchor stay fixed. Unknown/negative query space, reverse anchor footprints
// and emitted negative reference spans cannot receive candidate exact pairs.
// Accepted copy segment spans also remain protected. Reviews compare against
// sum-existing-forward-fragments: old_score sums actual clipped CIGAR fragment
// scores, old_cigar stays empty because no full old cross-region path exists.
// new_score/new_cigar describe the complete exact-gap candidate; its projected
// Internal globally unique 32-base pairings must exist on the candidate path;
// every such old internal coordinate pairing must remain. Purely repetitive
// cores keep the old copy assignment. The projected
// and actually restitched fragments must independently improve score and exact
// coverage with zero substitutions before initial/final are updated together.
[[nodiscard]] CentromereEndpointResult RefineCentromereCrossRegionBlocks(
    const SequenceRecord &, const SequenceRecord &,
    const CentromereDirectedAnchorResult &, const CentromereCopyAnchorResult *,
    const CentromereStitchResult &, const CentromereEndpointOptions & = {});
// Bounded joint competition of local rare chains before effective regions and
// roots are published. Opens unknown/weak-opposition space with two-sided
// sequence context, retains strong reverse/copy evidence and old unique pairs.
[[nodiscard]] CentromereEndpointResult RefineCentromereJointPaths(
    const SequenceRecord &, const SequenceRecord &,
    const CentromereDirectedAnchorResult &, const CentromereCopyAnchorResult *,
    const CentromereStitchResult &, const CentromereEndpointOptions & = {});
// Extend only empty, direction-unknown query space from the direction
// predictor's existing exact-anchor evidence. Never reopen an old path.
[[nodiscard]] CentromereEndpointResult RecoverCentromereUnplannedCores(
    const SequenceRecord &, const SequenceRecord &, const CentromereDirectedAnchorResult &,
    const CentromereCopyAnchorResult *, const CentromereStitchResult &,
    const std::filesystem::path &initial_evidence, const CentromereEndpointOptions & = {});

void WriteCentromereEndpointReviews(std::ostream &,
                                    const CentromereEndpointResult &);
void WriteCentromereCrossRegionReviews(std::ostream &,
                                       const CentromereEndpointResult &);
} // namespace ramag

// Private, non-owning interfaces shared by the five refinement stages.
// All returned containers own their storage; input sequence/anchor references
// remain valid for the duration of the synchronous stage call.

namespace ramag::centromere_refinement_detail {
using Clock = std::chrono::steady_clock;
using Range = std::pair<Position, Position>;
using A = CentromereDirectedAnchor;
struct Replacement {
  std::uint64_t region{};
  Position rb{}, re{}, qb{}, qe{};
  std::vector<A> exact;
};
struct CrossFragmentSummary {
  std::int64_t score{};
  Length exact{}, mismatches{}, short_disagreements{};
};
using ContextCounts = std::unordered_map<std::uint64_t, unsigned char>;
using CrossContexts = std::set<std::pair<Position, Position>>;

double Seconds(Clock::time_point s);
std::vector<Range> Union(std::vector<Range> x);
Length Bases(const std::vector<Range> &x);
bool Touch(Position a, Position b, Position c, Position d);
void Emit(std::vector<CigarOp> &v, char c, Length n);
void CheckAnchor(const A &a, const SequenceRecord &r, const SequenceRecord &q);
void Refresh(CentromereDirectedAnchorResult &out);
void Replace(std::vector<A> &list, const std::vector<Replacement> &replacements,
             const std::vector<CentromereDirectedRegion> &regions);
int ContextBits(char c);
ContextCounts CountContexts(const SequenceRecord &sequence,
                            const CentromereEndpointOptions &options);
bool DistinguishableContext(const A &anchor, const SequenceRecord &reference,
                            const ContextCounts &rc, const ContextCounts &qc);
bool CigarOffset(Position r, Position q, Position rb, Position qb, bool dr,
                 bool dq, Length n, Length &offset);
bool SliceStitchedCigar(const CentromereStitchResult &stitched,
                        std::uint64_t region, Position rb, Position re,
                        Position qb, Position qe, std::vector<CigarOp> &out);
Length ShortDisagreementBurden(const std::vector<CigarOp>&cigar);
CrossFragmentSummary CrossFragments(const CentromereStitchResult &stitched,
                                    const SequenceRecord &reference,
                                    const SequenceRecord &query, Position rb,
                                    Position re, Position qb, Position qe);
CrossContexts CrossPathContexts(const SequenceRecord &reference,
    const ContextCounts &rc, const ContextCounts &qc,
    const std::vector<CigarOp> &path, Position r, Position q,
    const CentromereEndpointReview &window);
CrossContexts CrossOldContexts(const SequenceRecord &r,
    const ContextCounts &rc,const ContextCounts &qc,
    const CentromereStitchResult &old,const CentromereEndpointReview &window);
std::vector<CigarOp> CrossGapSlice(const std::vector<CigarOp> &path,
                                 Position r, Position q, Position rb,
                                 Position re, Position qb, Position qe);
std::vector<CentromereStitchGap> CrossFixedGaps(
    const SequenceRecord &r, const SequenceRecord &q,
    const CentromereDirectedAnchorResult &candidate,
    const CentromereEndpointReview &review);
std::string CrossReverseComplement(const std::string &sequence);
} // namespace ramag::centromere_refinement_detail

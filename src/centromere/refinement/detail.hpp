#pragma once
// Private, non-owning interfaces shared by the five refinement stages.
// All returned containers own their storage; input sequence/anchor references
// remain valid for the duration of the synchronous stage call.
#include "refinement.hpp"
#include <chrono>
#include <set>
#include <unordered_map>

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

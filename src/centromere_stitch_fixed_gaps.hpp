#pragma once

#include "ramag/centromere_stitch.hpp"

namespace ramag::detail {
// Private consumer for a previously reviewed exact gap path. Each record must
// match a current adjacent-anchor gap by path role/region/branch, strand and
// physical endpoints. IDs are intentionally ignored. All records must be
// complete and scored; =/I/D consumption and score are replayed against the
// original FASTA. Invalid, duplicated or stale records are hard errors.
[[nodiscard]] CentromereStitchResult AlignCentromereAnchoredPathsWithFixedGaps(
    const SequenceRecord& reference, const SequenceRecord& query,
    const CentromereDirectedAnchorResult& primary,
    const CentromereCopyAnchorResult* accepted_copies,
    const CentromereStitchOptions& options,
    const std::vector<CentromereStitchGap>& fixed_gaps);
// Add validated exact singleton paths to an existing assembly, preserving
// every existing segment and gap record, including its path identity.
void AppendCentromereExactSingletons(
    const SequenceRecord &, const SequenceRecord &, const CentromereDirectedAnchorResult &,
    const CentromereCopyAnchorResult *, std::span<const CentromereDirectedAnchor>,
    CentromereStitchResult &);

}  // namespace ramag::detail

#pragma once

#include "ramag/centromere_stitch.hpp"
#include <iosfwd>

namespace ramag {

// Export each completed segment independently. Partial results retain their
// completed segments; aborted snapshots and malformed input are rejected before
// writing. No new pairs, tail extension, or cross-path stitching is performed.
// PAF and SAM tags preserve role/path identity; SAM secondary records can be
// accepted additional copy relations and must not be filtered for full scoring.
// All SAM records use hard clipping and contain their oriented query fragment.
// MAF uses ref.<name> and query.<name> to distinguish equal input names.
void WriteCentromereStitchPaf(std::ostream&, const SequenceRecord& reference,
                             const SequenceRecord& query, const CentromereStitchResult&);
void WriteCentromereStitchSam(std::ostream&, const SequenceRecord& reference,
                             const SequenceRecord& query, const CentromereStitchResult&);
void WriteCentromereStitchMaf(std::ostream&, const SequenceRecord& reference,
                             const SequenceRecord& query, const CentromereStitchResult&);

}  // namespace ramag

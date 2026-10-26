#pragma once
#include "ramag/centromere.hpp"
#include "ramag/centromere_blocks.hpp"
namespace ramag::detail {
class CentromereBudgetExceeded : public AlignmentError {
 public: using AlignmentError::AlignmentError;
};
// Private interval kernel, also exercised against an independent DP oracle.
std::vector<CigarOp> AlignRamaGap(std::string_view reference, std::string_view query,
                                const CentromereAlignmentOptions& options,
                                std::uint64_t memory_limit_bytes = 0);
CentromereBlockResult RefineDirectionRegions(const SequenceRecord&,const SequenceRecord&,
    CentromereBlockResult,const CentromereBlockOptions&,Strand);
// Restores the shared evidence store to the caller if planning fails.
CentromereBlockResult ReplanDirectionRegions(Length,Length,
    std::vector<CentromerePlanningAnchor>&,const CentromereBlockOptions&,Strand,
    std::string_view reference_bases={},std::string_view query_bases={});
}

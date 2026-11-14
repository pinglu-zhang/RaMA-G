#pragma once
#include "ramag/centromere_stitch.hpp"
#include "ramag/alignment.hpp"
#include <span>
#include <string_view>
namespace ramag::detail {
struct StitchGapResult {
  std::vector<CigarOp> cigar;
  std::int64_t score{};
  std::string engine;
  std::uint64_t peak_kernel_bytes{}, dp_cells{}, wfa_states{};
  double elapsed_seconds{};
};
class StitchGapResourceLimit : public AlignmentError {
 public:
  explicit StitchGapResourceLimit(std::string why) : AlignmentError(why), reason(std::move(why)) {}
  std::string reason, engine;
  std::uint64_t peak_kernel_bytes{}, dp_cells{}, wfa_states{};
  double elapsed_seconds{};
};
[[nodiscard]] StitchGapResult AlignCentromereStitchGap(
    std::string_view reference, std::string_view oriented_query,
    const CentromereStitchOptions& options);
[[nodiscard]] std::int64_t ScoreCentromereStitchCigar(
    std::string_view reference, std::string_view oriented_query,
    std::span<const CigarOp> cigar, CentromereStitchScoreModel model);
}  // namespace ramag::detail

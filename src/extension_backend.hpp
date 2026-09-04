#pragma once

#include "ramag/alignment.hpp"

#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

namespace ramag::internal {

struct GapAlignment {
    std::vector<CigarOp> cigar;
    std::int64_t score{};
};

struct ExtensionCallMetrics {
    std::uint64_t estimated_cells{};
    std::uint64_t full_matrix_cells{};
    std::uint64_t effective_band_width{};
    std::uint64_t effective_block_size{};
};

[[nodiscard]] GapAlignment AlignKsw2Gap(
    std::string_view reference,
    std::string_view query,
    const AlignmentOptions& options,
    bool automatic_band,
    ExtensionCallMetrics* metrics);

[[nodiscard]] GapAlignment AlignBlockGap(
    std::string_view reference,
    std::string_view query,
    const AlignmentOptions& options,
    bool exact,
    ExtensionCallMetrics* metrics);

[[nodiscard]] Length Ksw2AutomaticBandWidth(
    Length reference_length,
    Length query_length,
    const AlignmentOptions& options);

[[nodiscard]] std::string_view ConfiguredExtensionBackend();
[[nodiscard]] std::string ConfiguredExtensionRoute();

}  // namespace ramag::internal

#include "extension_backend.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace ramag::internal {

std::string_view ConfiguredExtensionBackend() {
    return RAMAG_EXTENSION_BACKEND;
}

std::string ConfiguredExtensionRoute() {
    const std::string_view backend = ConfiguredExtensionBackend();
    if (backend == "scalar") {
        return "exact+ungapped+bounded-scalar-dp";
    }
    if (backend == "ksw2-exact") {
        return "exact+ungapped+ksw2-global-exact";
    }
    if (backend == "ksw2-band-auto") {
        return "exact+ungapped+ksw2-green-global+ksw2-auto-band-v1";
    }
    if (backend == "block-exact") {
        return "exact+ungapped+block-aligner-full-block";
    }
    if (backend == "block-adaptive") {
        return "exact+ungapped+block-aligner-adaptive-32-128";
    }
    throw AlignmentError("unknown compiled extension backend");
}

Length Ksw2AutomaticBandWidth(Length reference_length,
                              Length query_length,
                              const AlignmentOptions& options) {
    const Length separation = std::max(reference_length, query_length);
    if (separation == 0) {
        return 0;
    }
    const Length endpoint_delta = reference_length > query_length
                                      ? reference_length - query_length
                                      : query_length - reference_length;
    const long double slack =
        static_cast<long double>(options.diag_diff) +
        static_cast<long double>(options.diag_factor) *
            static_cast<long double>(separation);
    if (!std::isfinite(slack) || slack < 0.0L ||
        slack > static_cast<long double>(std::numeric_limits<Length>::max())) {
        throw AlignmentError("KSW2 automatic band width is outside the 64-bit range");
    }
    const Length chain_slack = static_cast<Length>(std::ceil(slack));
    return std::min(separation, std::max(endpoint_delta, chain_slack));
}

}  // namespace ramag::internal

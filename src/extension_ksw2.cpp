#include "extension_backend.hpp"

#include <ksw2.h>

#include <algorithm>
#include <array>
#include <cstdlib>
#include <limits>
#include <memory>
#include <string>
#include <vector>

namespace ramag::internal {
namespace {

void Append(std::vector<CigarOp>& cigar, char operation, Length length) {
    if (length == 0) return;
    if (!cigar.empty() && cigar.back().operation == operation) {
        if (length > std::numeric_limits<Length>::max() - cigar.back().length) {
            throw AlignmentError("KSW2 CIGAR length exceeds the 64-bit range");
        }
        cigar.back().length += length;
    } else {
        cigar.push_back(CigarOp{operation, length});
    }
}

std::vector<std::uint8_t> Encode(std::string_view sequence) {
    std::vector<std::uint8_t> encoded;
    encoded.reserve(sequence.size());
    for (const char base : sequence) {
        switch (base) {
            case 'A': encoded.push_back(0); break;
            case 'C': encoded.push_back(1); break;
            case 'G': encoded.push_back(2); break;
            case 'T': encoded.push_back(3); break;
            case 'N': encoded.push_back(4); break;
            default:
                throw AlignmentError("KSW2 received a non-normalized base");
        }
    }
    return encoded;
}

void ValidateScores(const AlignmentOptions& options) {
    const auto fits_positive_i8 = [](std::int32_t value) {
        return value > 0 && value <= std::numeric_limits<std::int8_t>::max();
    };
    if (!fits_positive_i8(options.match_score) ||
        !fits_positive_i8(options.mismatch_penalty) ||
        !fits_positive_i8(options.gap_open_penalty) ||
        !fits_positive_i8(options.gap_extend_penalty)) {
        throw AlignmentError(
            "KSW2 requires positive scoring values representable as int8_t");
    }
    const std::int64_t vector_limit =
        static_cast<std::int64_t>(options.match_score) +
        2LL * (static_cast<std::int64_t>(options.gap_open_penalty) +
               static_cast<std::int64_t>(options.gap_extend_penalty));
    if (vector_limit > std::numeric_limits<std::int8_t>::max()) {
        throw AlignmentError("KSW2 scoring values exceed its SSE byte range");
    }
}

std::uint64_t FullCells(std::size_t reference_length,
                        std::size_t query_length) {
    const std::uint64_t rows = static_cast<std::uint64_t>(reference_length) + 1;
    const std::uint64_t columns = static_cast<std::uint64_t>(query_length) + 1;
    if (rows > std::numeric_limits<std::uint64_t>::max() / columns) {
        throw AlignmentError("KSW2 matrix cell estimate exceeds the 64-bit range");
    }
    return rows * columns;
}

std::uint64_t BandedCells(std::size_t reference_length,
                          std::size_t query_length,
                          std::uint64_t width) {
    std::uint64_t cells = 0;
    for (std::uint64_t row = 0;
         row <= static_cast<std::uint64_t>(reference_length); ++row) {
        const std::uint64_t begin = row > width ? row - width : 0;
        const std::uint64_t end = std::min<std::uint64_t>(
            static_cast<std::uint64_t>(query_length), row + width);
        if (begin <= end) {
            const std::uint64_t count = end - begin + 1;
            if (count > std::numeric_limits<std::uint64_t>::max() - cells) {
                throw AlignmentError("KSW2 band cell estimate exceeds the 64-bit range");
            }
            cells += count;
        }
    }
    return cells;
}

struct FreeCigar {
    void operator()(std::uint32_t* pointer) const noexcept {
        std::free(pointer);
    }
};

}  // namespace

GapAlignment AlignKsw2Gap(std::string_view reference,
                          std::string_view query,
                          const AlignmentOptions& options,
                          bool automatic_band,
                          ExtensionCallMetrics* metrics) {
    if (reference.empty() || query.empty()) {
        throw AlignmentError("KSW2 adapter requires non-empty sequences");
    }
    if (reference.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
        query.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw AlignmentError("KSW2 sequence length exceeds its int interface");
    }
    ValidateScores(options);

    const std::vector<std::uint8_t> encoded_reference = Encode(reference);
    const std::vector<std::uint8_t> encoded_query = Encode(query);
    std::array<std::int8_t, 25> matrix{};
    for (std::size_t row = 0; row < 5; ++row) {
        for (std::size_t column = 0; column < 5; ++column) {
            matrix[row * 5 + column] =
                row < 4 && row == column
                    ? static_cast<std::int8_t>(options.match_score)
                    : static_cast<std::int8_t>(-options.mismatch_penalty);
        }
    }

    int width = -1;
    Length automatic_width = 0;
    if (automatic_band) {
        automatic_width = Ksw2AutomaticBandWidth(
            static_cast<Length>(reference.size()),
            static_cast<Length>(query.size()), options);
        if (automatic_width > static_cast<Length>(std::numeric_limits<int>::max())) {
        throw AlignmentError("KSW2 automatic band width exceeds its int interface");
        }
        width = static_cast<int>(automatic_width);
    }

    int score = KSW_NEG_INF;
    int cigar_count = 0;
    std::uint32_t* raw_cigar = nullptr;
#if defined(RAMAG_EXTENSION_KSW2_EXACT)
    if (automatic_band) {
        throw AlignmentError("KSW2 exact binary requested the banded route");
    }
    ksw_extz_t extension{};
    ksw_extz2_sse(
        nullptr,
        static_cast<int>(encoded_query.size()), encoded_query.data(),
        static_cast<int>(encoded_reference.size()), encoded_reference.data(),
        5, matrix.data(),
        static_cast<std::int8_t>(options.gap_open_penalty),
        static_cast<std::int8_t>(options.gap_extend_penalty),
        width, -1, 0, KSW_EZ_GENERIC_SC, &extension);
    if (extension.zdropped != 0) {
        throw AlignmentError("KSW2 unexpectedly stopped at Z-drop");
    }
    score = extension.score;
    cigar_count = extension.n_cigar;
    raw_cigar = extension.cigar;
#elif defined(RAMAG_EXTENSION_KSW2_BAND_AUTO)
    if (!automatic_band) {
        throw AlignmentError("KSW2 band binary requested the exact route");
    }
    int cigar_capacity = 0;
    score = ksw_gg(
        nullptr,
        static_cast<int>(encoded_query.size()), encoded_query.data(),
        static_cast<int>(encoded_reference.size()), encoded_reference.data(),
        5, matrix.data(),
        static_cast<std::int8_t>(options.gap_open_penalty),
        static_cast<std::int8_t>(options.gap_extend_penalty),
        width, &cigar_capacity, &cigar_count, &raw_cigar);
#else
#error "extension_ksw2.cpp requires a KSW2 backend definition"
#endif
    std::unique_ptr<std::uint32_t, FreeCigar> cigar_owner(raw_cigar);
    if (score == KSW_NEG_INF || cigar_count <= 0 || raw_cigar == nullptr) {
        throw AlignmentError("KSW2 did not produce a complete global alignment");
    }

    GapAlignment converted;
    converted.score = score;
    std::size_t reference_offset = 0;
    std::size_t query_offset = 0;
    for (int index = 0; index < cigar_count; ++index) {
        const std::uint32_t packed = raw_cigar[index];
        const std::uint32_t operation = packed & 0x0fU;
        const Length length = static_cast<Length>(packed >> 4U);
        if (length == 0) {
            throw AlignmentError("KSW2 returned a zero-length CIGAR operation");
        }
        if (operation == KSW_CIGAR_MATCH) {
            if (length > reference.size() - reference_offset ||
                length > query.size() - query_offset) {
                throw AlignmentError("KSW2 match exceeds its input sequences");
            }
            for (Length offset = 0; offset < length; ++offset) {
                const char reference_base =
                    reference[reference_offset + static_cast<std::size_t>(offset)];
                const char query_base =
                    query[query_offset + static_cast<std::size_t>(offset)];
                Append(converted.cigar,
                       reference_base != 'N' && reference_base == query_base ? '=' : 'X',
                       1);
            }
            reference_offset += static_cast<std::size_t>(length);
            query_offset += static_cast<std::size_t>(length);
        } else if (operation == KSW_CIGAR_INS) {
            if (length > query.size() - query_offset) {
                throw AlignmentError("KSW2 insertion exceeds the query");
            }
            Append(converted.cigar, 'I', length);
            query_offset += static_cast<std::size_t>(length);
        } else if (operation == KSW_CIGAR_DEL) {
            if (length > reference.size() - reference_offset) {
                throw AlignmentError("KSW2 deletion exceeds the reference");
            }
            Append(converted.cigar, 'D', length);
            reference_offset += static_cast<std::size_t>(length);
        } else {
            throw AlignmentError("KSW2 returned an unsupported CIGAR operation");
        }
    }
    if (reference_offset != reference.size() || query_offset != query.size()) {
        throw AlignmentError("KSW2 CIGAR does not reach the global endpoint");
    }

    if (metrics != nullptr) {
        metrics->effective_band_width = automatic_width;
        metrics->full_matrix_cells = FullCells(reference.size(), query.size());
        metrics->estimated_cells = automatic_band
                                       ? BandedCells(reference.size(), query.size(),
                                                     automatic_width)
                                       : metrics->full_matrix_cells;
    }
    return converted;
}

}  // namespace ramag::internal

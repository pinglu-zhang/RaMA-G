#include "extension_backend.hpp"

#include <block_aligner.h>

#include <algorithm>
#include <array>
#include <limits>
#include <string_view>

namespace ramag::internal {
namespace {

constexpr std::size_t kMinimumBlockSize = 32;
constexpr std::size_t kAdaptiveMaximumBlockSize = 128;
constexpr std::size_t kExactMaximumBlockSize = 16'384;

void Append(std::vector<CigarOp>& cigar, char operation, Length length) {
    if (length == 0) return;
    if (!cigar.empty() && cigar.back().operation == operation) {
        if (length > std::numeric_limits<Length>::max() - cigar.back().length) {
            throw AlignmentError("Block Aligner CIGAR exceeds the 64-bit range");
        }
        cigar.back().length += length;
    } else {
        cigar.push_back(CigarOp{operation, length});
    }
}

void ValidateSequence(std::string_view sequence) {
    for (const char base : sequence) {
        if (base != 'A' && base != 'C' && base != 'G' &&
            base != 'T' && base != 'N') {
            throw AlignmentError("Block Aligner received a non-normalized base");
        }
    }
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
            "Block Aligner requires positive scoring values representable as int8_t");
    }
    const std::int64_t first_gap_cost =
        static_cast<std::int64_t>(options.gap_open_penalty) +
        static_cast<std::int64_t>(options.gap_extend_penalty);
    if (first_gap_cost > std::numeric_limits<std::int8_t>::max()) {
        throw AlignmentError("Block Aligner gap-open value exceeds int8_t");
    }
}

std::size_t ExactBlockSize(std::size_t reference_length,
                           std::size_t query_length) {
    const std::size_t maximum = std::max(reference_length, query_length);
    if (maximum == std::numeric_limits<std::size_t>::max()) {
        throw AlignmentError("Block Aligner full-block dimension overflows size_t");
    }
    const std::size_t required = std::max(kMinimumBlockSize, maximum + 1);
    std::size_t block = 1;
    while (block < required) {
        if (block > kExactMaximumBlockSize / 2) {
            throw AlignmentError(
                "Block Aligner exact block would exceed the 16384-cell limit");
        }
        block *= 2;
    }
    return block;
}

std::uint64_t FullCells(std::size_t reference_length,
                        std::size_t query_length) {
    const std::uint64_t rows = static_cast<std::uint64_t>(reference_length) + 1;
    const std::uint64_t columns = static_cast<std::uint64_t>(query_length) + 1;
    if (rows > std::numeric_limits<std::uint64_t>::max() / columns) {
        throw AlignmentError("Block Aligner cell estimate exceeds the 64-bit range");
    }
    return rows * columns;
}

class BlockContext {
public:
    ~BlockContext() { ResetStorage(); ResetMatrix(); }

    BlockContext(const BlockContext&) = delete;
    BlockContext& operator=(const BlockContext&) = delete;
    BlockContext() = default;

    void Prepare(std::size_t query_length,
                 std::size_t reference_length,
                 std::size_t block_size,
                 const AlignmentOptions& options) {
        PrepareMatrix(options);
        if (handle_ == nullptr || query_length > query_capacity_ ||
            reference_length > reference_capacity_ ||
            block_size > block_capacity_) {
            ResetStorage();
            query_capacity_ = std::max(query_capacity_, query_length);
            reference_capacity_ = std::max(reference_capacity_, reference_length);
            block_capacity_ = std::max(block_capacity_, block_size);
            query_ = block_new_padded_aa(query_capacity_, block_capacity_);
            reference_ = block_new_padded_aa(reference_capacity_, block_capacity_);
            handle_ = block_new_aa_trace(
                query_capacity_, reference_capacity_, block_capacity_);
            if (query_ == nullptr || reference_ == nullptr || handle_ == nullptr) {
                ResetStorage();
                throw AlignmentError("Block Aligner failed to allocate its worker context");
            }
        }
    }

    void SetSequences(std::string_view query, std::string_view reference) {
        block_set_bytes_padded_aa(
            query_, reinterpret_cast<const std::uint8_t*>(query.data()),
            query.size(), block_capacity_);
        block_set_bytes_padded_aa(
            reference_, reinterpret_cast<const std::uint8_t*>(reference.data()),
            reference.size(), block_capacity_);
    }

    [[nodiscard]] BlockHandle handle() const noexcept { return handle_; }
    [[nodiscard]] const PaddedBytes* query() const noexcept { return query_; }
    [[nodiscard]] const PaddedBytes* reference() const noexcept { return reference_; }
    [[nodiscard]] const AAMatrix* matrix() const noexcept { return matrix_; }

private:
    void PrepareMatrix(const AlignmentOptions& options) {
        if (matrix_ != nullptr && match_ == options.match_score &&
            mismatch_ == options.mismatch_penalty) {
            return;
        }
        ResetMatrix();
        matrix_ = block_new_simple_aamatrix(
            static_cast<std::int8_t>(options.match_score),
            static_cast<std::int8_t>(-options.mismatch_penalty));
        if (matrix_ == nullptr) {
            throw AlignmentError("Block Aligner failed to allocate its scoring matrix");
        }
        constexpr std::array<std::uint8_t, 5> bases{'A', 'C', 'G', 'T', 'N'};
        for (const std::uint8_t left : bases) {
            for (const std::uint8_t right : bases) {
                const std::int8_t score =
                    left != 'N' && left == right
                        ? static_cast<std::int8_t>(options.match_score)
                        : static_cast<std::int8_t>(-options.mismatch_penalty);
                block_set_aamatrix(matrix_, left, right, score);
            }
        }
        match_ = options.match_score;
        mismatch_ = options.mismatch_penalty;
    }

    void ResetStorage() noexcept {
        if (handle_ != nullptr) block_free_aa_trace(handle_);
        if (query_ != nullptr) block_free_padded_aa(query_);
        if (reference_ != nullptr) block_free_padded_aa(reference_);
        handle_ = nullptr;
        query_ = nullptr;
        reference_ = nullptr;
    }

    void ResetMatrix() noexcept {
        if (matrix_ != nullptr) block_free_aamatrix(matrix_);
        matrix_ = nullptr;
    }

    BlockHandle handle_{};
    PaddedBytes* query_{};
    PaddedBytes* reference_{};
    AAMatrix* matrix_{};
    std::size_t query_capacity_{};
    std::size_t reference_capacity_{};
    std::size_t block_capacity_{};
    std::int32_t match_{};
    std::int32_t mismatch_{};
};

thread_local BlockContext context;

}  // namespace

GapAlignment AlignBlockGap(std::string_view reference,
                           std::string_view query,
                           const AlignmentOptions& options,
                           bool exact,
                           ExtensionCallMetrics* metrics) {
    if (reference.empty() || query.empty()) {
        throw AlignmentError("Block Aligner adapter requires non-empty sequences");
    }
    ValidateSequence(reference);
    ValidateSequence(query);
    ValidateScores(options);

    const std::size_t block_size = exact
                                       ? ExactBlockSize(reference.size(), query.size())
                                       : kAdaptiveMaximumBlockSize;
    context.Prepare(query.size(), reference.size(), block_size, options);
    context.SetSequences(query, reference);

    const Gaps gaps{
        static_cast<std::int8_t>(
            -(options.gap_open_penalty + options.gap_extend_penalty)),
        static_cast<std::int8_t>(-options.gap_extend_penalty)};
    const SizeRange sizes{
        exact ? block_size : kMinimumBlockSize,
        block_size};
    block_align_aa_trace(context.handle(), context.query(), context.reference(),
                         context.matrix(), gaps, sizes, 0);
    const AlignResult result = block_res_aa_trace(context.handle());
    if (result.query_idx != query.size() ||
        result.reference_idx != reference.size()) {
        throw AlignmentError("Block Aligner did not reach the global endpoint");
    }

    Cigar* raw_cigar = block_new_cigar(result.query_idx, result.reference_idx);
    if (raw_cigar == nullptr) {
        throw AlignmentError("Block Aligner failed to allocate its CIGAR");
    }
    struct CigarGuard {
        Cigar* value;
        ~CigarGuard() { block_free_cigar(value); }
    } cigar_guard{raw_cigar};
    block_cigar_aa_trace(context.handle(), result.query_idx,
                         result.reference_idx, raw_cigar);

    GapAlignment converted;
    converted.score = result.score;
    std::size_t reference_offset = 0;
    std::size_t query_offset = 0;
    const std::size_t operation_count = block_len_cigar(raw_cigar);
    for (std::size_t index = 0; index < operation_count; ++index) {
        const OpLen operation = block_get_cigar(raw_cigar, index);
        const Length length = static_cast<Length>(operation.len);
        if (length == 0) {
            throw AlignmentError("Block Aligner returned a zero-length CIGAR operation");
        }
        if (operation.op == M || operation.op == Eq || operation.op == X) {
            if (length > reference.size() - reference_offset ||
                length > query.size() - query_offset) {
                throw AlignmentError("Block Aligner match exceeds its input sequences");
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
        } else if (operation.op == I) {
            if (length > query.size() - query_offset) {
                throw AlignmentError("Block Aligner insertion exceeds the query");
            }
            Append(converted.cigar, 'I', length);
            query_offset += static_cast<std::size_t>(length);
        } else if (operation.op == D) {
            if (length > reference.size() - reference_offset) {
                throw AlignmentError("Block Aligner deletion exceeds the reference");
            }
            Append(converted.cigar, 'D', length);
            reference_offset += static_cast<std::size_t>(length);
        } else {
            throw AlignmentError("Block Aligner returned an unsupported CIGAR operation");
        }
    }
    if (reference_offset != reference.size() || query_offset != query.size()) {
        throw AlignmentError("Block Aligner CIGAR does not consume both sequences");
    }

    if (metrics != nullptr) {
        metrics->effective_block_size = static_cast<std::uint64_t>(block_size);
        metrics->full_matrix_cells = FullCells(reference.size(), query.size());
        metrics->estimated_cells = exact
                                       ? metrics->full_matrix_cells
                                       : 0;
    }
    return converted;
}

}  // namespace ramag::internal

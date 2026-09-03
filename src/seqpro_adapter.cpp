#include "ramag/seqpro_adapter.hpp"

#include "ramag/runtime.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cctype>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <system_error>
#include <unordered_set>
#include <utility>

#ifndef RAMAG_HAVE_SEQPRO
#define RAMAG_HAVE_SEQPRO 0
#endif

#if RAMAG_HAVE_SEQPRO
#include <seqpro/error.h>
#include <seqpro/fasta_index.h>
#include <seqpro/indexed_fasta.h>
#endif

namespace ramag {
namespace {

[[maybe_unused, nodiscard]] std::string PathForMessage(
    const std::filesystem::path& path) {
    return path.empty() ? std::string{"<empty path>"} : path.string();
}

struct FileSha256 {
    std::uint64_t size_bytes{0};
    std::string hex_digest;
};

class Sha256State {
  public:
    void Update(const unsigned char* data, std::size_t size_bytes) {
        if (size_bytes >
            static_cast<std::size_t>(std::numeric_limits<std::uint64_t>::max() -
                                     total_size_bytes_)) {
            throw FastaError("file is too large for SHA-256 size accounting");
        }
        total_size_bytes_ += static_cast<std::uint64_t>(size_bytes);
        std::size_t consumed = 0;
        while (consumed < size_bytes) {
            const std::size_t available = block_.size() - buffered_size_bytes_;
            const std::size_t copied = std::min(available, size_bytes - consumed);
            std::copy_n(data + consumed, copied,
                        block_.begin() +
                            static_cast<std::ptrdiff_t>(buffered_size_bytes_));
            buffered_size_bytes_ += copied;
            consumed += copied;
            if (buffered_size_bytes_ == block_.size()) {
                Transform(block_.data());
                buffered_size_bytes_ = 0;
            }
        }
    }

    [[nodiscard]] std::string Finalize() {
        if (total_size_bytes_ >
            std::numeric_limits<std::uint64_t>::max() / 8U) {
            throw FastaError("file is too large for SHA-256 bit-length encoding");
        }
        const std::uint64_t bit_length = total_size_bytes_ * 8U;
        block_[buffered_size_bytes_++] = 0x80U;
        if (buffered_size_bytes_ > 56U) {
            std::fill(block_.begin() +
                          static_cast<std::ptrdiff_t>(buffered_size_bytes_),
                      block_.end(), 0U);
            Transform(block_.data());
            buffered_size_bytes_ = 0;
        }
        std::fill(block_.begin() +
                      static_cast<std::ptrdiff_t>(buffered_size_bytes_),
                  block_.begin() + 56, 0U);
        for (std::size_t byte_index = 0; byte_index < 8U; ++byte_index) {
            const unsigned int shift =
                static_cast<unsigned int>((7U - byte_index) * 8U);
            block_[56U + byte_index] =
                static_cast<unsigned char>((bit_length >> shift) & 0xffU);
        }
        Transform(block_.data());

        std::ostringstream digest;
        digest << std::hex << std::setfill('0');
        for (const std::uint32_t word : state_) {
            digest << std::setw(8) << word;
        }
        return digest.str();
    }

  private:
    static constexpr std::array<std::uint32_t, 64> kRoundConstants{
        0x428a2f98U, 0x71374491U, 0xb5c0fbcfU, 0xe9b5dba5U,
        0x3956c25bU, 0x59f111f1U, 0x923f82a4U, 0xab1c5ed5U,
        0xd807aa98U, 0x12835b01U, 0x243185beU, 0x550c7dc3U,
        0x72be5d74U, 0x80deb1feU, 0x9bdc06a7U, 0xc19bf174U,
        0xe49b69c1U, 0xefbe4786U, 0x0fc19dc6U, 0x240ca1ccU,
        0x2de92c6fU, 0x4a7484aaU, 0x5cb0a9dcU, 0x76f988daU,
        0x983e5152U, 0xa831c66dU, 0xb00327c8U, 0xbf597fc7U,
        0xc6e00bf3U, 0xd5a79147U, 0x06ca6351U, 0x14292967U,
        0x27b70a85U, 0x2e1b2138U, 0x4d2c6dfcU, 0x53380d13U,
        0x650a7354U, 0x766a0abbU, 0x81c2c92eU, 0x92722c85U,
        0xa2bfe8a1U, 0xa81a664bU, 0xc24b8b70U, 0xc76c51a3U,
        0xd192e819U, 0xd6990624U, 0xf40e3585U, 0x106aa070U,
        0x19a4c116U, 0x1e376c08U, 0x2748774cU, 0x34b0bcb5U,
        0x391c0cb3U, 0x4ed8aa4aU, 0x5b9cca4fU, 0x682e6ff3U,
        0x748f82eeU, 0x78a5636fU, 0x84c87814U, 0x8cc70208U,
        0x90befffaU, 0xa4506cebU, 0xbef9a3f7U, 0xc67178f2U};

    static std::uint32_t Choose(std::uint32_t first,
                                std::uint32_t second,
                                std::uint32_t third) noexcept {
        return (first & second) ^ (~first & third);
    }

    static std::uint32_t Majority(std::uint32_t first,
                                  std::uint32_t second,
                                  std::uint32_t third) noexcept {
        return (first & second) ^ (first & third) ^ (second & third);
    }

    static std::uint32_t BigSigmaZero(std::uint32_t value) noexcept {
        return std::rotr(value, 2) ^ std::rotr(value, 13) ^
               std::rotr(value, 22);
    }

    static std::uint32_t BigSigmaOne(std::uint32_t value) noexcept {
        return std::rotr(value, 6) ^ std::rotr(value, 11) ^
               std::rotr(value, 25);
    }

    static std::uint32_t SmallSigmaZero(std::uint32_t value) noexcept {
        return std::rotr(value, 7) ^ std::rotr(value, 18) ^ (value >> 3U);
    }

    static std::uint32_t SmallSigmaOne(std::uint32_t value) noexcept {
        return std::rotr(value, 17) ^ std::rotr(value, 19) ^ (value >> 10U);
    }

    void Transform(const unsigned char* block) {
        std::array<std::uint32_t, 64> schedule{};
        for (std::size_t word_index = 0; word_index < 16U; ++word_index) {
            const std::size_t offset = word_index * 4U;
            schedule[word_index] =
                (static_cast<std::uint32_t>(block[offset]) << 24U) |
                (static_cast<std::uint32_t>(block[offset + 1U]) << 16U) |
                (static_cast<std::uint32_t>(block[offset + 2U]) << 8U) |
                static_cast<std::uint32_t>(block[offset + 3U]);
        }
        for (std::size_t word_index = 16U; word_index < schedule.size();
             ++word_index) {
            schedule[word_index] =
                SmallSigmaOne(schedule[word_index - 2U]) +
                schedule[word_index - 7U] +
                SmallSigmaZero(schedule[word_index - 15U]) +
                schedule[word_index - 16U];
        }

        std::uint32_t first = state_[0];
        std::uint32_t second = state_[1];
        std::uint32_t third = state_[2];
        std::uint32_t fourth = state_[3];
        std::uint32_t fifth = state_[4];
        std::uint32_t sixth = state_[5];
        std::uint32_t seventh = state_[6];
        std::uint32_t eighth = state_[7];
        for (std::size_t round = 0; round < schedule.size(); ++round) {
            const std::uint32_t temporary_one =
                eighth + BigSigmaOne(fifth) + Choose(fifth, sixth, seventh) +
                kRoundConstants[round] + schedule[round];
            const std::uint32_t temporary_two =
                BigSigmaZero(first) + Majority(first, second, third);
            eighth = seventh;
            seventh = sixth;
            sixth = fifth;
            fifth = fourth + temporary_one;
            fourth = third;
            third = second;
            second = first;
            first = temporary_one + temporary_two;
        }
        state_[0] += first;
        state_[1] += second;
        state_[2] += third;
        state_[3] += fourth;
        state_[4] += fifth;
        state_[5] += sixth;
        state_[6] += seventh;
        state_[7] += eighth;
    }

    std::array<std::uint32_t, 8> state_{
        0x6a09e667U, 0xbb67ae85U, 0x3c6ef372U, 0xa54ff53aU,
        0x510e527fU, 0x9b05688cU, 0x1f83d9abU, 0x5be0cd19U};
    std::array<unsigned char, 64> block_{};
    std::size_t buffered_size_bytes_{0};
    std::uint64_t total_size_bytes_{0};
};

[[maybe_unused, nodiscard]] FileSha256 Sha256File(
    const std::filesystem::path& path) {
    std::ifstream input(path, std::ios::binary);
    if (!input) {
        throw FastaError("cannot open FASTA index for SHA-256: '" +
                         PathForMessage(path) + "'");
    }
    Sha256State state;
    std::array<char, 64U * 1024U> buffer{};
    std::uint64_t size_bytes = 0;
    while (input) {
        input.read(buffer.data(),
                   static_cast<std::streamsize>(buffer.size()));
        const std::streamsize available = input.gcount();
        if (available > 0) {
            const auto available_size = static_cast<std::size_t>(available);
            if (available_size >
                static_cast<std::size_t>(
                    std::numeric_limits<std::uint64_t>::max() - size_bytes)) {
                throw FastaError("FASTA index is too large for 64-bit provenance: '" +
                                 PathForMessage(path) + "'");
            }
            state.Update(
                reinterpret_cast<const unsigned char*>(buffer.data()),
                available_size);
            size_bytes += static_cast<std::uint64_t>(available_size);
        }
    }
    if (input.bad()) {
        throw FastaError("cannot read FASTA index for SHA-256: '" +
                         PathForMessage(path) + "'");
    }
    return FileSha256{size_bytes, state.Finalize()};
}

[[maybe_unused, nodiscard]] std::filesystem::path AbsoluteOrOriginal(
    const std::filesystem::path& path) {
    std::error_code error;
    auto absolute = std::filesystem::absolute(path, error);
    return error ? path : absolute.lexically_normal();
}

[[maybe_unused]] void RequireWorkDirectory(const std::filesystem::path& fasta_path,
                                           const std::filesystem::path& work_dir) {
    if (work_dir.empty()) {
        throw FastaError(
            "SeqPro work directory is required; refusing to place an index "
            "beside the FASTA implicitly");
    }

    std::error_code error;
    std::filesystem::create_directories(work_dir, error);
    if (error) {
        throw FastaError("cannot create SeqPro work directory '" +
                         PathForMessage(work_dir) + "': " + error.message());
    }
    if (!std::filesystem::is_directory(work_dir, error) || error) {
        throw FastaError("SeqPro work path is not a directory: '" +
                         PathForMessage(work_dir) + "'");
    }

    // The adapter must never publish <input>.fai or another sidecar directly
    // in the input directory.  Resolve symlinks when possible so an apparently
    // different spelling cannot bypass that contract.
    const auto input_parent =
        AbsoluteOrOriginal(fasta_path).parent_path().lexically_normal();
    const auto absolute_work = AbsoluteOrOriginal(work_dir).lexically_normal();
    std::error_code equivalent_error;
    const bool equivalent =
        std::filesystem::equivalent(input_parent, absolute_work, equivalent_error);
    if ((!equivalent_error && equivalent) || input_parent == absolute_work) {
        throw FastaError(
            "SeqPro work directory must differ from the FASTA input directory; "
            "refusing to write index sidecars beside '" +
            PathForMessage(fasta_path) + "'");
    }
}

[[maybe_unused, nodiscard]] std::filesystem::path ResolveIndexPath(
    const SeqProInputOptions& options) {
    if (options.index_basename.empty()) {
        throw FastaError(
            "SeqPro index basename is required (use distinct reference/query roles)");
    }
    const std::filesystem::path basename{options.index_basename};
    if (basename.has_parent_path() || basename.filename() != basename ||
        basename == "." || basename == "..") {
        throw FastaError("SeqPro index basename must be a file name, not a path: '" +
                         options.index_basename + "'");
    }
    for (const unsigned char byte : options.index_basename) {
        const bool alpha_numeric =
            (byte >= static_cast<unsigned char>('A') &&
             byte <= static_cast<unsigned char>('Z')) ||
            (byte >= static_cast<unsigned char>('a') &&
             byte <= static_cast<unsigned char>('z')) ||
            (byte >= static_cast<unsigned char>('0') &&
             byte <= static_cast<unsigned char>('9'));
        if (!alpha_numeric && byte != static_cast<unsigned char>('.') &&
            byte != static_cast<unsigned char>('_') &&
            byte != static_cast<unsigned char>('-')) {
            throw FastaError(
                "SeqPro index basename may contain only ASCII letters, digits, '.', "
                "'_' and '-': '" +
                options.index_basename + "'");
        }
    }

    std::string file_name = options.index_basename;
    if (file_name.size() < 4 || file_name.substr(file_name.size() - 4) != ".fai") {
        file_name += ".fai";
    }
    return options.work_dir / file_name;
}

struct ExternalFastaIndexCopy {
    std::filesystem::path source_path;
    FileSha256 source_digest;
    bool copied{false};
};

[[maybe_unused, nodiscard]] ExternalFastaIndexCopy
CopyExternalFastaIndexIfPresent(const std::filesystem::path& fasta_path,
                                const std::filesystem::path& destination_path) {
    const std::filesystem::path source_path{fasta_path.string() + ".fai"};
    std::error_code path_error;
    const bool source_exists = std::filesystem::exists(source_path, path_error);
    if (path_error) {
        throw FastaError("cannot inspect input-side FASTA index '" +
                         PathForMessage(source_path) + "': " +
                         path_error.message());
    }
    if (!source_exists) {
        return {};
    }
    if (!std::filesystem::is_regular_file(source_path, path_error) || path_error) {
        throw FastaError("input-side FASTA index is not a regular file: '" +
                         PathForMessage(source_path) +
                         (path_error ? "' (" + path_error.message() + ")" : "'"));
    }

    const FileSha256 digest_before_copy = Sha256File(source_path);
    if (!std::filesystem::copy_file(source_path, destination_path,
                                    std::filesystem::copy_options::none,
                                    path_error)) {
        throw FastaError("cannot copy input-side FASTA index '" +
                         PathForMessage(source_path) + "' to isolated work path '" +
                         PathForMessage(destination_path) + "'" +
                         (path_error ? ": " + path_error.message() : ""));
    }
    const FileSha256 copied_digest = Sha256File(destination_path);
    const FileSha256 source_digest_after_copy = Sha256File(source_path);
    if (copied_digest.size_bytes != digest_before_copy.size_bytes ||
        copied_digest.hex_digest != digest_before_copy.hex_digest ||
        source_digest_after_copy.size_bytes != digest_before_copy.size_bytes ||
        source_digest_after_copy.hex_digest != digest_before_copy.hex_digest) {
        throw FastaError(
            "input-side FASTA index changed while it was copied to the run work "
            "directory: '" +
            PathForMessage(source_path) + "'");
    }
    return ExternalFastaIndexCopy{AbsoluteOrOriginal(source_path),
                                  digest_before_copy, true};
}

[[maybe_unused, nodiscard]] bool IsAsciiAlphabetic(unsigned char byte) noexcept {
    return (byte >= static_cast<unsigned char>('A') &&
            byte <= static_cast<unsigned char>('Z')) ||
           (byte >= static_cast<unsigned char>('a') &&
            byte <= static_cast<unsigned char>('z'));
}

[[maybe_unused, nodiscard]] char NormalizeBase(unsigned char byte) noexcept {
    switch (byte) {
        case 'A':
        case 'a':
            return 'A';
        case 'C':
        case 'c':
            return 'C';
        case 'G':
        case 'g':
            return 'G';
        case 'T':
        case 't':
            return 'T';
        default:
            return 'N';
    }
}

[[maybe_unused, nodiscard]] std::vector<std::string> ReadCompleteHeaders(
    const std::filesystem::path& fasta_path) {
    std::ifstream input(fasta_path, std::ios::binary);
    if (!input) {
        throw FastaError("cannot open FASTA header stream '" +
                         PathForMessage(fasta_path) + "'");
    }
    std::vector<std::string> headers;
    std::array<char, 64U * 1024U> buffer{};
    std::string header;
    bool at_line_start = true;
    bool collecting_header = false;
    const auto finish_header = [&]() {
        std::size_t begin = 0;
        while (begin < header.size() &&
               std::isspace(static_cast<unsigned char>(header[begin])) != 0) {
            ++begin;
        }
        std::size_t end = header.size();
        while (end > begin &&
               std::isspace(static_cast<unsigned char>(header[end - 1])) != 0) {
            --end;
        }
        headers.emplace_back(header.substr(begin, end - begin));
        header.clear();
    };
    while (input) {
        CheckInterruption("input");
        input.read(buffer.data(), static_cast<std::streamsize>(buffer.size()));
        const auto available = input.gcount();
        for (std::streamsize index = 0; index < available; ++index) {
            const char value = buffer[static_cast<std::size_t>(index)];
            if (value == '\n') {
                if (collecting_header) finish_header();
                collecting_header = false;
                at_line_start = true;
            } else if (at_line_start && value == '>') {
                collecting_header = true;
                at_line_start = false;
            } else {
                if (collecting_header) header.push_back(value);
                at_line_start = false;
            }
        }
    }
    if (input.bad()) {
        throw FastaError("I/O failure while reading FASTA headers '" +
                         PathForMessage(fasta_path) + "'");
    }
    if (collecting_header) finish_header();
    return headers;
}

[[maybe_unused]] void CheckedAdd(Length& value,
                                 Length increment,
                                 std::string_view description) {
    if (increment > std::numeric_limits<Length>::max() - value) {
        throw FastaError(std::string{description} +
                         " exceeds the 64-bit RaMA-G coordinate range");
    }
    value += increment;
}

[[maybe_unused, noreturn]] void ThrowUnavailable() {
    throw FastaError(
        "the SeqPro FASTA adapter is unavailable because RaMA-G was built "
        "without the pinned SeqPro dependency");
}

#if RAMAG_HAVE_SEQPRO

[[nodiscard]] std::string ErrorCodeName(seqpro::ErrorCode code) {
    switch (code) {
        case seqpro::ErrorCode::kInvalidArgument:
            return "invalid-argument";
        case seqpro::ErrorCode::kIoError:
            return "io-error";
        case seqpro::ErrorCode::kInvalidFasta:
            return "invalid-fasta";
        case seqpro::ErrorCode::kInvalidFastaIndex:
            return "invalid-fasta-index";
        case seqpro::ErrorCode::kStaleFastaIndex:
            return "stale-fasta-index";
        case seqpro::ErrorCode::kDuplicateSequenceName:
            return "duplicate-sequence-name";
        case seqpro::ErrorCode::kSequenceNotFound:
            return "sequence-not-found";
        case seqpro::ErrorCode::kSequenceRangeOutOfBounds:
            return "sequence-range-out-of-bounds";
        case seqpro::ErrorCode::kIntegerOverflow:
            return "integer-overflow";
        case seqpro::ErrorCode::kUnsupportedFileFormat:
            return "unsupported-file-format";
    }
    return "unknown";
}

[[nodiscard]] std::string BuildActionName(seqpro::FastaIndexBuildAction action) {
    switch (action) {
        case seqpro::FastaIndexBuildAction::kCreated:
            return "created";
        case seqpro::FastaIndexBuildAction::kReused:
            return "reused";
        case seqpro::FastaIndexBuildAction::kAdoptedExternalIndex:
            return "adopted-external-index";
        case seqpro::FastaIndexBuildAction::kRebuilt:
            return "rebuilt";
    }
    return "unknown";
}

[[nodiscard]] std::string IndexOriginName(seqpro::FastaIndexOrigin origin) {
    switch (origin) {
        case seqpro::FastaIndexOrigin::kSeqProVerified:
            return "seqpro-verified";
        case seqpro::FastaIndexOrigin::kExternalStandardFai:
            return "external-standard-fai";
    }
    return "unknown";
}

[[nodiscard]] std::string VerificationStatusName(
    seqpro::IndexVerificationStatus status) {
    switch (status) {
        case seqpro::IndexVerificationStatus::kStructureValidated:
            return "structure-validated";
        case seqpro::IndexVerificationStatus::kMetadataValidated:
            return "metadata-validated";
        case seqpro::IndexVerificationStatus::kFullContentValidated:
            return "full-content-validated";
    }
    return "unknown";
}

[[noreturn]] void RethrowSeqProFailure(std::string_view operation,
                                       const seqpro::SeqProError& error) {
    const std::string message = "SeqPro " + std::string{operation} + " failed [" +
                                ErrorCodeName(error.error_code()) + "]: " +
                                error.what();
    if (error.error_code() == seqpro::ErrorCode::kUnsupportedFileFormat) {
        throw UnsupportedFastaFormat(message);
    }
    throw FastaError(message);
}

#endif

}  // namespace

bool SeqProAdapterAvailable() noexcept {
#if RAMAG_HAVE_SEQPRO
    return true;
#else
    return false;
#endif
}

SeqProInputResult ReadFastaWithSeqPro(const std::filesystem::path& fasta_path,
                                     const SeqProInputOptions& options) {
#if RAMAG_HAVE_SEQPRO
    if (fasta_path.empty()) {
        throw FastaError("FASTA path is empty");
    }
    ValidateUncompressedFasta(fasta_path);
    RequireWorkDirectory(fasta_path, options.work_dir);
    const std::filesystem::path fasta_index_path = ResolveIndexPath(options);

    try {
        const ExternalFastaIndexCopy external_index =
            CopyExternalFastaIndexIfPresent(fasta_path, fasta_index_path);
        seqpro::FastaIndexBuildOptions build_options;
        build_options.fasta_index_path = fasta_index_path;
        build_options.force_rebuild =
            external_index.copied ? false : options.force_rebuild;
        build_options.write_seqpro_metadata = !external_index.copied;
        const seqpro::FastaIndexBuildReport build_report =
            seqpro::BuildFastaIndex(fasta_path, build_options);

        seqpro::IndexedFastaOptions open_options;
        open_options.fasta_index_path = fasta_index_path;
        open_options.file_access_pattern = seqpro::FileAccessPattern::kSequential;
        open_options.index_verification_mode =
            seqpro::IndexVerificationMode::kFast;
        open_options.require_seqpro_metadata = !external_index.copied;
        const seqpro::IndexedFasta indexed =
            seqpro::IndexedFasta::Open(fasta_path, open_options);

        if (external_index.copied) {
            if (build_report.build_action !=
                    seqpro::FastaIndexBuildAction::kReused ||
                !build_report.metadata_path.empty() ||
                indexed.fasta_index_origin() !=
                    seqpro::FastaIndexOrigin::kExternalStandardFai ||
                indexed.index_verification_status() !=
                    seqpro::IndexVerificationStatus::kStructureValidated) {
                throw FastaError(
                    "SeqPro external-FAI route did not preserve the required "
                    "structure-only validation provenance");
            }
        } else if (build_report.metadata_path.empty() ||
                   indexed.fasta_index_origin() !=
                       seqpro::FastaIndexOrigin::kSeqProVerified ||
                   indexed.index_verification_status() !=
                       seqpro::IndexVerificationStatus::kMetadataValidated) {
            throw FastaError(
                "SeqPro full-scan route did not produce validated metadata");
        }

        if (indexed.sequence_count() == 0) {
            throw FastaError("FASTA contains no records: '" +
                             PathForMessage(fasta_path) + "'");
        }
        const std::uint64_t record_count =
            static_cast<std::uint64_t>(indexed.sequence_count());
        const std::uint64_t available_ids =
            static_cast<std::uint64_t>(
                std::numeric_limits<SequenceId>::max() - options.first_id) +
            1U;
        if (record_count > available_ids) {
            throw FastaError(
                "FASTA contains too many records for 32-bit RaMA-G sequence ids");
        }
        const auto complete_headers = ReadCompleteHeaders(fasta_path);
        if (complete_headers.size() != indexed.sequence_count()) {
            throw FastaError(
                "SeqPro record count differs from the FASTA header catalog for '" +
                PathForMessage(fasta_path) + "'");
        }

        SeqProInputResult result;
        result.fasta.source = AbsoluteOrOriginal(fasta_path);
        result.fasta_index_path = build_report.fasta_index_path;
        result.metadata_path = build_report.metadata_path;
        result.source_fasta_index_path = external_index.source_path;
        result.source_fasta_index_size_bytes =
            external_index.source_digest.size_bytes;
        result.source_fasta_index_sha256 =
            external_index.source_digest.hex_digest;
        result.source_fasta_index_copied = external_index.copied;
        result.external_fasta_index_adopted = external_index.copied;
        result.build_action = BuildActionName(build_report.build_action);
        result.index_origin = IndexOriginName(indexed.fasta_index_origin());
        result.verification_status =
            VerificationStatusName(indexed.index_verification_status());
        result.fasta.sequences.reserve(indexed.sequence_count());

        std::unordered_set<std::string> seen_names;
        for (std::size_t ordinal = 0; ordinal < indexed.sequence_count(); ++ordinal) {
            if (ordinal > static_cast<std::size_t>(
                              std::numeric_limits<seqpro::SequenceId>::max())) {
                throw FastaError(
                    "SeqPro record ordinal exceeds its 32-bit sequence-id range");
            }
            const auto view = indexed.SequenceById(
                static_cast<seqpro::SequenceId>(ordinal));
            const std::uint64_t sequence_length = view.sequence_length();
            if (sequence_length == 0) {
                throw FastaError("FASTA record '" +
                                 std::string{view.sequence_name()} + "' in '" +
                                 PathForMessage(fasta_path) +
                                 "' has an empty sequence");
            }
            if (sequence_length >
                static_cast<std::uint64_t>(
                    std::numeric_limits<std::size_t>::max())) {
                throw FastaError("FASTA record '" +
                                 std::string{view.sequence_name()} +
                                 "' cannot be materialized in this process");
            }

            SequenceRecord record;
            record.numeric_id = static_cast<SequenceId>(
                static_cast<std::uint64_t>(options.first_id) + ordinal);
            record.name = std::string{view.sequence_name()};
            record.header = complete_headers[ordinal];
            if (record.name.empty()) {
                throw FastaError("SeqPro returned an empty FASTA record name for '" +
                                 PathForMessage(fasta_path) + "'");
            }
            const auto header_separator = record.header.find_first_of(
                " \t\v\f\r\n");
            const auto header_name = record.header.substr(0, header_separator);
            if (header_name != record.name) {
                throw FastaError(
                    "SeqPro record name differs from FASTA header at ordinal " +
                    std::to_string(ordinal) + " in '" +
                    PathForMessage(fasta_path) + "'");
            }
            if (!seen_names.insert(record.name).second) {
                throw FastaError("duplicate FASTA record identifier '" +
                                 record.name + "' in '" +
                                 PathForMessage(fasta_path) + "'");
            }

            record.bases.resize(static_cast<std::size_t>(sequence_length));
            view.CopySubsequenceTo(0, record.bases.data(), record.bases.size());
            for (std::size_t position = 0; position < record.bases.size();
                 ++position) {
                const unsigned char byte =
                    static_cast<unsigned char>(record.bases[position]);
                if (!IsAsciiAlphabetic(byte)) {
                    std::ostringstream diagnostic;
                    diagnostic << "invalid non-alphabetic FASTA byte 0x" << std::hex
                               << std::uppercase << std::setw(2) << std::setfill('0')
                               << static_cast<unsigned int>(byte) << std::dec
                               << " at one-based position " << (position + 1)
                               << " of record '" << record.name << "' in '"
                               << PathForMessage(fasta_path) << "'";
                    throw FastaError(diagnostic.str());
                }
                const char normalized = NormalizeBase(byte);
                record.bases[position] = normalized;
                if (normalized == 'N') {
                    CheckedAdd(result.fasta.ambiguous_bases, 1,
                               "FASTA ambiguous-base count");
                }
            }
            CheckedAdd(result.fasta.total_bases, sequence_length,
                       "FASTA base count");
            result.fasta.sequences.push_back(std::move(record));
        }

        return result;
    } catch (const UnsupportedFastaFormat&) {
        throw;
    } catch (const FastaError&) {
        throw;
    } catch (const seqpro::SeqProError& error) {
        RethrowSeqProFailure("FASTA index/open/read", error);
    }
#else
    static_cast<void>(fasta_path);
    static_cast<void>(options);
    ThrowUnavailable();
#endif
}

}  // namespace ramag

#pragma once

#include "ramag/model.hpp"

#include <cstdint>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace ramag {

class FastaError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class UnsupportedFastaFormat : public FastaError {
public:
    using FastaError::FastaError;
};

struct FastaData {
    std::filesystem::path source;
    std::vector<SequenceRecord> sequences;
    Length total_bases{};
    Length ambiguous_bases{};
};

enum class FastaCompression : std::uint8_t { Plain, Gzip };

// Detect plain versus gzip/BGZF input from magic bytes. A gzip-looking suffix
// must contain a valid gzip header; other compression formats are rejected.
[[nodiscard]] FastaCompression DetectFastaCompression(
    const std::filesystem::path& path);

// Read a strict plain or gzip FASTA file. IDs are the first header token and
// must be unique.  ASCII letters other than A/C/G/T are normalized to N; other
// sequence characters are rejected.  Records with no bases are rejected.
[[nodiscard]] FastaData ReadFasta(const std::filesystem::path& path,
                                  SequenceId first_id = 0);

[[nodiscard]] bool IsCanonicalBase(char base) noexcept;

// Reverse-complement a normalized sequence.  Lower-case A/C/G/T are accepted;
// every non-canonical character is returned as N.
[[nodiscard]] std::string ReverseComplement(std::string_view bases);

// Reverse-complement an ordered collection with one bounded OpenMP team. The
// output order is identical to the input order, each output byte has exactly
// one writer, and an already-active OpenMP region forces the scalar path so a
// caller cannot accidentally multiply the shared worker budget. Small
// collections stay scalar to avoid making thread startup the dominant cost.
[[nodiscard]] std::vector<std::string> ReverseComplements(
    const std::vector<SequenceRecord>& sequences,
    std::uint32_t worker_threads = 1);

}  // namespace ramag

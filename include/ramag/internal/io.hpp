#pragma once
#include "ramag/model.hpp"
#include <cstdint>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>
#include <iosfwd>
#include <span>



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



namespace ramag {

struct OutputPaths {
  std::filesystem::path sam;
  std::filesystem::path paf;
  std::filesystem::path delta;
  std::filesystem::path maf;
  std::filesystem::path chain;
};

class WriterError : public std::runtime_error {
 public:
  using std::runtime_error::runtime_error;
};

[[nodiscard]] OutputPaths MakeOutputPaths(const std::filesystem::path& prefix);
void WriteSam(std::ostream& output,
              const FastaData& references,
              const FastaData& queries,
              std::span<const AlignmentRecord> alignments,
              std::string_view command_line);
void WritePaf(std::ostream& output,
              const FastaData& references,
              const FastaData& queries,
              std::span<const AlignmentRecord> alignments);
void WriteDelta(std::ostream& output,
                const FastaData& references,
                const FastaData& queries,
                std::span<const AlignmentRecord> alignments);
void WriteMaf(std::ostream& output,
              const FastaData& references,
              const FastaData& queries,
              std::span<const AlignmentRecord> alignments);
void WriteChain(std::ostream& output,
                const FastaData& references,
                const FastaData& queries,
                std::span<const AlignmentRecord> alignments);

void ValidateSamFile(const std::filesystem::path& path);
void ValidatePafFile(const std::filesystem::path& path);
void ValidateDeltaFile(const std::filesystem::path& path);
void ValidateMafFile(const std::filesystem::path& path);
void ValidateChainFile(const std::filesystem::path& path);

[[nodiscard]] std::string SamCigar(const AlignmentRecord& alignment,
                                   Length query_length);
[[nodiscard]] std::string BuildMdTag(const AlignmentRecord& alignment,
                                     const SequenceRecord& reference,
                                     const SequenceRecord& query);
[[nodiscard]] std::vector<std::int64_t> CigarToDelta(
    std::span<const CigarOp> cigar);

}  // namespace ramag

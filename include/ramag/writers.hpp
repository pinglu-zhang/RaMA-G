#pragma once

#include "ramag/cli.hpp"
#include "ramag/fasta.hpp"
#include "ramag/model.hpp"

#include <filesystem>
#include <iosfwd>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

namespace ramag {

struct OutputPaths {
  std::filesystem::path sam;
  std::filesystem::path paf;
  std::filesystem::path delta;
  std::filesystem::path maf;
  std::filesystem::path chain;
  std::filesystem::path manifest;
  std::filesystem::path complete;
};

class WriterError : public std::runtime_error {
 public:
  using std::runtime_error::runtime_error;
};

[[nodiscard]] OutputPaths MakeOutputPaths(const std::filesystem::path& prefix);
[[nodiscard]] OutputPaths MakeOutputPaths(const RunSpec& spec);
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

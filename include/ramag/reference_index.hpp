#pragma once

#include "ramag/cli.hpp"
#include "ramag/fasta.hpp"

#include <filesystem>
#include <string>

namespace ramag {

struct ReferenceIndexPaths {
  std::filesystem::path index;
  std::filesystem::path manifest;
  std::filesystem::path complete;
};

[[nodiscard]] ReferenceIndexPaths MakeReferenceIndexPaths(
    const std::filesystem::path& index);

[[nodiscard]] ReferenceIndexPaths RunReferenceIndexPipeline(
    const IndexSpec& spec,
    std::string invocation,
    const std::filesystem::path& binary_path);

// Require the RaMA-G companion manifest/marker and bind them to the supplied
// normalized reference before Sufkit loads the index payload.
void ValidateReferenceIndexBundle(const std::filesystem::path& index,
                                  const FastaData& reference);

}  // namespace ramag

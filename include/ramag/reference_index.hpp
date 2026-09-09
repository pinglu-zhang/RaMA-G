#pragma once

#include "ramag/cli.hpp"
#include "ramag/fasta.hpp"

#include <filesystem>
#include <string>
#include <map>

namespace ramag {
class SufkitSeedIndex;
class ProgressSession;
class RunLogger;

struct ReferenceIndexPaths {
  std::filesystem::path index;
};

[[nodiscard]] ReferenceIndexPaths MakeReferenceIndexPaths(
    const std::filesystem::path& index);

[[nodiscard]] ReferenceIndexPaths RunReferenceIndexPipeline(
    const IndexSpec& spec,
    std::string invocation,
    const std::filesystem::path& binary_path,
    const CpuAffinityInfo& launch_affinity, RunLogger* logger = nullptr);

// Publish the existing in-memory object; never rebuild or reparse the reference.
// The published index has an independent lifetime from alignment outputs.
[[nodiscard]] ReferenceIndexPaths SaveReferenceIndex(
    const IndexSpec& spec, const FastaData& reference,
    const SufkitSeedIndex& index, std::string invocation,
    const std::filesystem::path& binary_path, ProgressSession& progress,
    double build_seconds, std::map<std::string, double>& timings, RunLogger* logger = nullptr);

}  // namespace ramag

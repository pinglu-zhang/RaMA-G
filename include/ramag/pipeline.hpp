#pragma once

#include "ramag/cli.hpp"
#include "ramag/model.hpp"
#include "ramag/writers.hpp"

#include <filesystem>
#include <string>

namespace ramag {

class RunLogger;

struct RunOutcome {
  OutputPaths paths;
  RunStatistics statistics;
};

struct BatchQueryOutcome {
  std::string name;
  std::string status{"not-run"};
  int exit_code{};
  double elapsed_seconds{};
  std::uint64_t alignment_count{};
  std::filesystem::path log_path;
  std::string message;
};

struct BatchOutcome {
  std::vector<BatchQueryOutcome> queries;
  int exit_code{};
  std::uint64_t index_build_calls{};
  std::uint64_t index_load_calls{};
};

[[nodiscard]] BatchOutcome RunBatchPipeline(const BatchSpec& spec,
    std::string invocation, std::filesystem::path binary_path,
    RunLogger* logger = nullptr);

[[nodiscard]] RunOutcome RunAlignmentPipeline(const RunSpec& spec,
                                              std::string invocation,
                                              std::filesystem::path binary_path, RunLogger* logger = nullptr);

}  // namespace ramag

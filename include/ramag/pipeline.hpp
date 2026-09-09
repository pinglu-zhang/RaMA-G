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

[[nodiscard]] RunOutcome RunAlignmentPipeline(const RunSpec& spec,
                                              std::string invocation,
                                              std::filesystem::path binary_path, RunLogger* logger = nullptr);

}  // namespace ramag

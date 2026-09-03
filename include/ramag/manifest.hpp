#pragma once

#include "ramag/cli.hpp"
#include "ramag/fasta.hpp"
#include "ramag/model.hpp"
#include "ramag/writers.hpp"

#include <cstdint>
#include <filesystem>
#include <map>
#include <span>
#include <string>
#include <vector>

namespace ramag {

struct ArtifactReport {
  std::string format;
  std::filesystem::path path;
  std::uintmax_t bytes{0};
  std::string state{"validated"};
};

struct ManifestData {
  std::string run_id;
  std::string started_utc;
  std::string finished_utc;
  RunSpec run_spec;
  const FastaData* reference{nullptr};
  const FastaData* query{nullptr};
  RunStatistics statistics;
  std::vector<ArtifactReport> artifacts;
  std::map<std::string, double> stage_wall_seconds;
  std::map<std::string, std::string> adapter_provenance;
  std::string invocation;
  std::filesystem::path binary_path;
  std::string input_route;
  std::string index_route;
  std::string seeding_route;
  std::string chaining_route{"sparse-exact-edge-components-v1"};
  std::string extension_route{"exact+ungapped+bounded-scalar-dp"};
  std::string input_parallel_route{"serial-reference-query-v1"};
  std::uint32_t input_requested_workers{1};
  std::uint32_t input_actual_workers{1};
  std::uint32_t actual_threads{1};
  std::string status{"success"};
  int exit_code{0};
};

[[nodiscard]] std::string ManifestJson(const ManifestData& data);

}  // namespace ramag

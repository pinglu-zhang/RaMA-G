#pragma once

#include "ramag/alignment.hpp"

#include <cstdint>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace ramag {

struct FormatSelection {
  bool sam{true};
  bool paf{true};
  bool delta{true};
  bool maf{false};
  bool chain{false};

  [[nodiscard]] bool Any() const noexcept {
    return sam || paf || delta || maf || chain;
  }
  [[nodiscard]] std::string ToString() const;
};

enum class OutputFormat : std::uint8_t { Sam, Paf, Delta, Maf, Chain };

struct OutputRequest {
  OutputFormat format{OutputFormat::Paf};
  std::filesystem::path path;
};

enum class ProgressMode : std::uint8_t { Auto, On, Off };

struct ProgressOptions {
  ProgressMode mode{ProgressMode::Auto};
  std::uint32_t interval_seconds{10};
  bool interval_explicit{false};
};

struct RunSpec {
  std::filesystem::path reference_path;
  std::filesystem::path reference_index_path;
  std::filesystem::path query_path;
  std::filesystem::path output_prefix;
  std::vector<OutputRequest> outputs;
  std::filesystem::path work_dir;
  std::uint32_t threads{1};
  FormatSelection formats{};
  ProgressOptions progress{};
  AlignmentOptions alignment{};
};

struct IndexSpec {
  std::filesystem::path reference_path;
  std::filesystem::path output_path;
  std::filesystem::path work_dir;
  std::uint32_t threads{1};
  ProgressOptions progress{};
};

enum class CommandKind : std::uint8_t { Align, Index };

struct OpenMpRuntimeInfo {
  bool enabled{false};
  std::string runtime{"unavailable"};
  std::int32_t specification_date{0};
  std::int32_t max_threads{1};
  std::int32_t thread_limit{1};
  std::int32_t dynamic{0};
  std::int32_t max_active_levels{1};
  std::uint32_t configured_requested_threads{0};
};

struct CliParseResult {
  CommandKind command{CommandKind::Align};
  RunSpec run_spec;
  IndexSpec index_spec;
  bool show_help{false};
  bool show_version{false};
  bool print_effective_config{false};
};

class CliError : public std::runtime_error {
 public:
  using std::runtime_error::runtime_error;
};

[[nodiscard]] CliParseResult ParseCommandLine(int argc, const char* const* argv);
[[nodiscard]] std::string HelpText();
[[nodiscard]] std::string VersionText();
[[nodiscard]] std::string SeedModeName(SeedMode mode);
[[nodiscard]] std::string AlignmentSelectionName(AlignmentSelection selection);
[[nodiscard]] std::string OutputFormatName(OutputFormat format);
[[nodiscard]] std::string ProgressModeName(ProgressMode mode);
[[nodiscard]] std::string EffectiveConfigText(const RunSpec& spec);
[[nodiscard]] std::string EffectiveIndexConfigText(const IndexSpec& spec);
[[nodiscard]] OpenMpRuntimeInfo CurrentOpenMpRuntimeInfo();
void ConfigureOpenMpRuntime(const RunSpec& spec);
void ConfigureOpenMpRuntime(const IndexSpec& spec);
void ValidateRunSpec(const RunSpec& spec);
void ValidateIndexSpec(const IndexSpec& spec);

}  // namespace ramag

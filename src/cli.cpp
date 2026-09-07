#include "ramag/cli.hpp"

#include <algorithm>
#include <atomic>
#include <cctype>
#include <charconv>
#include <cmath>
#include <cerrno>
#include <cstring>
#include <limits>
#include <sstream>
#include <system_error>

#if defined(__linux__)
#include <sched.h>
#endif

#if RAMAG_OPENMP_ENABLED
#include <omp.h>
#endif

namespace ramag {
namespace {

std::atomic<std::uint32_t> configured_openmp_requested_threads{0};

template <typename Integer>
Integer ParseInteger(std::string_view text, std::string_view option) {
  Integer value{};
  const auto result = std::from_chars(text.data(), text.data() + text.size(), value);
  if (result.ec != std::errc{} || result.ptr != text.data() + text.size()) {
    throw CliError("invalid integer for " + std::string(option) + ": " +
                   std::string(text));
  }
  return value;
}

double ParseDouble(std::string_view text, std::string_view option) {
  std::size_t consumed = 0;
  double value = 0.0;
  try {
    value = std::stod(std::string(text), &consumed);
  } catch (const std::exception&) {
    throw CliError("invalid number for " + std::string(option) + ": " +
                   std::string(text));
  }
  if (consumed != text.size() || !std::isfinite(value)) {
    throw CliError("invalid number for " + std::string(option) + ": " +
                   std::string(text));
  }
  return value;
}

std::string RequireValue(int argc, const char* const* argv, int& index,
                         std::string_view option) {
  if (index + 1 >= argc) throw CliError("missing value for " + std::string(option));
  return argv[++index];
}

SeedMode ParseSeedMode(std::string_view text) {
  if (text == "fast") return SeedMode::FastHierarchical;
  if (text == "mumreference") return SeedMode::MumReference;
  if (text == "mum") return SeedMode::Mum;
  if (text == "smem") return SeedMode::Smem;
  if (text == "maxmatch") return SeedMode::MaxMatch;
  throw CliError("unknown --seed-mode value: " + std::string(text));
}

AlignmentSelection ParseAlignmentSelection(std::string_view text) {
  if (text == "all") return AlignmentSelection::All;
  if (text == "one-to-one") return AlignmentSelection::ReciprocalOneToOne;
  throw CliError("unknown --selection-mode value: " + std::string(text));
}

ProgressMode ParseProgressMode(std::string_view text) {
  if (text == "auto") return ProgressMode::Auto;
  if (text == "on") return ProgressMode::On;
  if (text == "off") return ProgressMode::Off;
  throw CliError("unknown --progress value: " + std::string(text));
}

void SelectFormat(FormatSelection& formats, std::string_view token,
                  std::string_view option) {
  bool* selected = nullptr;
  if (token == "sam") selected = &formats.sam;
  else if (token == "paf") selected = &formats.paf;
  else if (token == "delta") selected = &formats.delta;
  else if (token == "maf") selected = &formats.maf;
  else if (token == "chain") selected = &formats.chain;
  else throw CliError("unknown " + std::string(option) + " entry: " +
                      std::string(token));
  if (*selected) throw CliError("duplicate " + std::string(option) + " entry: " +
                                std::string(token));
  *selected = true;
}

FormatSelection ParseFormats(std::string_view text) {
  FormatSelection formats{false, false, false, false, false};
  std::size_t begin = 0;
  while (begin <= text.size()) {
    const auto comma = text.find(',', begin);
    const auto end = comma == std::string_view::npos ? text.size() : comma;
    SelectFormat(formats, text.substr(begin, end - begin), "--formats");
    if (comma == std::string_view::npos) break;
    begin = comma + 1;
  }
  if (!formats.Any()) throw CliError("--formats must select at least one output format");
  return formats;
}

OutputFormat InferOutputFormat(const std::filesystem::path& path) {
  const auto extension = path.extension().string();
  if (extension == ".sam") return OutputFormat::Sam;
  if (extension == ".paf") return OutputFormat::Paf;
  if (extension == ".delta") return OutputFormat::Delta;
  if (extension == ".maf") return OutputFormat::Maf;
  if (extension == ".chain") return OutputFormat::Chain;
  throw CliError("--output must end in .sam, .paf, .delta, .maf, or .chain: " +
                 path.string());
}

std::filesystem::path OutputBase(const std::filesystem::path& path) {
  auto result = path;
  result.replace_extension();
  return result.lexically_normal();
}

void AddExplicitOutput(RunSpec& spec, const std::filesystem::path& path) {
  if (path.empty()) throw CliError("--output path must not be empty");
  const auto format = InferOutputFormat(path);
  if (std::any_of(spec.outputs.begin(), spec.outputs.end(),
                  [format](const OutputRequest& value) { return value.format == format; })) {
    throw CliError("duplicate --output format: " + OutputFormatName(format));
  }
  const auto normalized = path.lexically_normal();
  if (std::any_of(spec.outputs.begin(), spec.outputs.end(),
                  [&](const OutputRequest& value) {
                    return value.path.lexically_normal() == normalized;
                  })) {
    throw CliError("duplicate --output path: " + path.string());
  }
  const auto base = OutputBase(path);
  if (spec.outputs.empty()) {
    spec.output_prefix = base;
    spec.formats = FormatSelection{false, false, false, false, false};
  } else if (base != spec.output_prefix.lexically_normal()) {
    throw CliError("all --output paths must share one directory and base prefix");
  }
  switch (format) {
    case OutputFormat::Sam: spec.formats.sam = true; break;
    case OutputFormat::Paf: spec.formats.paf = true; break;
    case OutputFormat::Delta: spec.formats.delta = true; break;
    case OutputFormat::Maf: spec.formats.maf = true; break;
    case OutputFormat::Chain: spec.formats.chain = true; break;
  }
  spec.outputs.push_back(OutputRequest{format, path});
}

std::string AbsoluteForDisplay(const std::filesystem::path& path) {
  if (path.empty()) return {};
  std::error_code error;
  const auto absolute = std::filesystem::absolute(path, error);
  return error ? path.string() : absolute.lexically_normal().string();
}

void ValidateInput(const std::filesystem::path& path, std::string_view label) {
  if (path.empty()) throw CliError(std::string(label) + " path is empty");
  std::error_code error;
  if (!std::filesystem::is_regular_file(path, error) || error) {
    throw CliError(std::string(label) + " is not a readable regular file: " + path.string());
  }
}

void ValidateWorkDir(const std::filesystem::path& path) {
  if (path.empty()) throw CliError("--work-dir must not be empty");
  std::error_code error;
  const bool exists = std::filesystem::exists(path, error);
  if (error) throw CliError("cannot inspect --work-dir: " + error.message());
  if (exists && !std::filesystem::is_directory(path, error)) {
    throw CliError("--work-dir exists but is not a directory: " + path.string());
  }
  if (error) throw CliError("cannot inspect --work-dir type: " + error.message());
}

void ValidateOutputParent(const std::filesystem::path& path) {
  if (path.parent_path().empty()) return;
  std::error_code error;
  const bool exists = std::filesystem::exists(path.parent_path(), error);
  if (error) throw CliError("cannot inspect output parent: " + error.message());
  if (exists && !std::filesystem::is_directory(path.parent_path(), error)) {
    throw CliError("output parent exists but is not a directory: " +
                   path.parent_path().string());
  }
}

void ValidateProgress(const ProgressOptions& options) {
  if (options.interval_seconds == 0) {
    throw CliError("--progress-interval must be at least 1");
  }
  if (options.mode == ProgressMode::Off && options.interval_explicit) {
    throw CliError("--progress-interval cannot be used with --progress off");
  }
}

template <typename Spec>
void ConfigureOpenMpRuntimeImpl(const Spec& spec) {
#if RAMAG_OPENMP_ENABLED
  if (spec.threads > static_cast<std::uint32_t>(std::numeric_limits<int>::max())) {
    throw CliError("--threads exceeds the OpenMP runtime integer limit");
  }
  omp_set_dynamic(0);
  omp_set_max_active_levels(1);
  omp_set_num_threads(static_cast<int>(spec.threads));
  configured_openmp_requested_threads.store(spec.threads, std::memory_order_relaxed);
  const auto info = CurrentOpenMpRuntimeInfo();
  if (info.dynamic != 0 || info.max_active_levels != 1 ||
      info.max_threads != static_cast<std::int32_t>(spec.threads) ||
      info.thread_limit < static_cast<std::int32_t>(spec.threads)) {
    throw CliError("OpenMP runtime cannot honor the requested --threads budget");
  }
#else
  static_cast<void>(spec);
  throw CliError("this RaMA-G binary was built without required OpenMP support");
#endif
}

}  // namespace

OpenMpRuntimeInfo CurrentOpenMpRuntimeInfo() {
  OpenMpRuntimeInfo info;
  info.runtime = RAMAG_OPENMP_RUNTIME;
  info.specification_date = RAMAG_OPENMP_SPEC_DATE;
  info.configured_requested_threads =
      configured_openmp_requested_threads.load(std::memory_order_relaxed);
#if RAMAG_OPENMP_ENABLED
  info.enabled = true;
  info.max_threads = omp_get_max_threads();
  info.thread_limit = omp_get_thread_limit();
  info.dynamic = omp_get_dynamic();
  info.max_active_levels = omp_get_max_active_levels();
#endif
  return info;
}

std::string CpuAffinityInfo::CpuList() const {
  if (!supported) return "unavailable";
  if (logical_cpus.empty()) return "empty";
  std::ostringstream output;
  for (std::size_t begin = 0; begin < logical_cpus.size();) {
    std::size_t end = begin;
    while (end + 1 < logical_cpus.size() &&
           logical_cpus[end + 1] == logical_cpus[end] + 1U) {
      ++end;
    }
    if (begin != 0) output << ',';
    output << logical_cpus[begin];
    if (end != begin) output << '-' << logical_cpus[end];
    begin = end + 1;
  }
  return output.str();
}

CpuAffinityInfo CurrentCpuAffinity() {
  CpuAffinityInfo info;
#if defined(__linux__)
  cpu_set_t mask;
  CPU_ZERO(&mask);
  if (sched_getaffinity(0, sizeof(mask), &mask) != 0) {
    throw CliError("cannot inspect CPU affinity: " +
                   std::string(std::strerror(errno)));
  }
  info.supported = true;
  info.source = "sched-affinity";
  for (std::uint32_t cpu = 0; cpu < CPU_SETSIZE; ++cpu) {
    if (CPU_ISSET(static_cast<int>(cpu), &mask)) {
      info.logical_cpus.push_back(cpu);
    }
  }
  if (info.logical_cpus.empty()) {
    throw CliError("CPU affinity contains no available logical CPUs");
  }
#endif
  return info;
}

void ValidateCpuAffinityBudget(const CpuAffinityInfo& affinity,
                               std::uint32_t requested_threads,
                               std::string_view context) {
  if (!affinity.supported) return;
  if (affinity.logical_cpus.size() <
      static_cast<std::size_t>(requested_threads)) {
    throw CliError(std::string(context) + " requested " +
                   std::to_string(requested_threads) +
                   " threads but current CPU affinity permits only " +
                   std::to_string(affinity.logical_cpus.size()) +
                   " logical CPUs (allowed=" + affinity.CpuList() + ")");
  }
}

namespace {

void RestoreCpuAffinity(const CpuAffinityInfo& affinity) {
#if defined(__linux__)
  if (!affinity.supported) return;
  cpu_set_t mask;
  CPU_ZERO(&mask);
  for (const auto cpu : affinity.logical_cpus) {
    if (cpu >= CPU_SETSIZE) {
      throw CliError("launch CPU affinity contains an unsupported CPU id: " +
                     std::to_string(cpu));
    }
    CPU_SET(static_cast<int>(cpu), &mask);
  }
  if (sched_setaffinity(0, sizeof(mask), &mask) != 0) {
    throw CliError("cannot restore launch CPU affinity: " +
                   std::string(std::strerror(errno)));
  }
#else
  static_cast<void>(affinity);
#endif
}

}  // namespace

void ConfigureOpenMpRuntime(const RunSpec& spec) { ConfigureOpenMpRuntimeImpl(spec); }

void ConfigureIndexRuntime(const IndexSpec& spec,
                           const CpuAffinityInfo& launch_affinity) {
  ValidateCpuAffinityBudget(launch_affinity, spec.threads,
                            "reference-index construction");
  ConfigureOpenMpRuntimeImpl(spec);
  // libgomp may bind the calling thread to the first OpenMP place while its
  // runtime is initialized. CaPS uses a separate std::thread/parlay scheduler,
  // whose workers inherit the caller's affinity. Restore exactly the CPU set
  // granted at process launch and never widen beyond taskset/cgroup limits.
  RestoreCpuAffinity(launch_affinity);
  const auto restored = CurrentCpuAffinity();
  if (launch_affinity.supported &&
      restored.logical_cpus != launch_affinity.logical_cpus) {
    throw CliError("restored CPU affinity differs from the process launch mask "
                   "(launch=" + launch_affinity.CpuList() +
                   ", restored=" + restored.CpuList() + ")");
  }
  ValidateCpuAffinityBudget(restored, spec.threads,
                            "reference-index construction");
}

std::string FormatSelection::ToString() const {
  std::string result;
  const auto append = [&](std::string_view value) {
    if (!result.empty()) result.push_back(',');
    result.append(value);
  };
  if (sam) append("sam");
  if (paf) append("paf");
  if (delta) append("delta");
  if (maf) append("maf");
  if (chain) append("chain");
  return result;
}

std::string OutputFormatName(OutputFormat format) {
  switch (format) {
    case OutputFormat::Sam: return "sam";
    case OutputFormat::Paf: return "paf";
    case OutputFormat::Delta: return "delta";
    case OutputFormat::Maf: return "maf";
    case OutputFormat::Chain: return "chain";
  }
  throw CliError("invalid internal output format");
}

std::string ProgressModeName(ProgressMode mode) {
  switch (mode) {
    case ProgressMode::Auto: return "auto";
    case ProgressMode::On: return "on";
    case ProgressMode::Off: return "off";
  }
  throw CliError("invalid internal progress mode");
}

std::string SeedModeName(SeedMode mode) {
  switch (mode) {
    case SeedMode::FastHierarchical: return "fast";
    case SeedMode::MumReference: return "mumreference";
    case SeedMode::Mum: return "mum";
    case SeedMode::Smem: return "smem";
    case SeedMode::MaxMatch: return "maxmatch";
  }
  throw CliError("invalid internal seed mode");
}

std::string AlignmentSelectionName(AlignmentSelection selection) {
  switch (selection) {
    case AlignmentSelection::All: return "all";
    case AlignmentSelection::ReciprocalOneToOne: return "one-to-one";
  }
  throw CliError("invalid internal alignment selection mode");
}

CliParseResult ParseCommandLine(int argc, const char* const* argv) {
  CliParseResult parsed;
  if (argc <= 1) { parsed.show_help = true; return parsed; }
  const std::string_view first = argv[1];
  if (first == "--help" || first == "-h" || first == "help") {
    parsed.show_help = true; return parsed;
  }
  if (first == "--version" || first == "version") {
    parsed.show_version = true; return parsed;
  }

  if (first == "index") {
    parsed.command = CommandKind::Index;
    bool have_reference = false;
    bool have_output = false;
    bool have_work_dir = false;
    for (int index = 2; index < argc; ++index) {
      const std::string_view option = argv[index];
      if (option == "--help" || option == "-h") parsed.show_help = true;
      else if (option == "--version") parsed.show_version = true;
      else if (option == "--print-effective-config") parsed.print_effective_config = true;
      else if (option == "--reference") {
        parsed.index_spec.reference_path = RequireValue(argc, argv, index, option);
        have_reference = true;
      } else if (option == "--output") {
        parsed.index_spec.output_path = RequireValue(argc, argv, index, option);
        have_output = true;
      } else if (option == "--work-dir") {
        parsed.index_spec.work_dir = RequireValue(argc, argv, index, option);
        have_work_dir = true;
      } else if (option == "--threads") {
        parsed.index_spec.threads = ParseInteger<std::uint32_t>(RequireValue(argc, argv, index, option), option);
      } else if (option == "--progress") {
        parsed.index_spec.progress.mode = ParseProgressMode(RequireValue(argc, argv, index, option));
      } else if (option == "--progress-interval") {
        parsed.index_spec.progress.interval_seconds = ParseInteger<std::uint32_t>(RequireValue(argc, argv, index, option), option);
        parsed.index_spec.progress.interval_explicit = true;
      } else throw CliError("unknown index option: " + std::string(option));
    }
    if (parsed.show_help || parsed.show_version) return parsed;
    if (!have_reference) throw CliError("missing required option --reference");
    if (!have_output) throw CliError("missing required option --output");
    if (!have_work_dir) throw CliError("missing required option --work-dir");
    ValidateIndexSpec(parsed.index_spec);
    return parsed;
  }

  if (first != "align") {
    throw CliError("expected subcommand 'align' or 'index'; got: " + std::string(first));
  }
  parsed.command = CommandKind::Align;
  bool have_reference = false;
  bool have_query = false;
  bool have_prefix = false;
  bool have_work_dir = false;
  bool have_formats = false;
  bool have_smem_min_occurrences = false;
  for (int index = 2; index < argc; ++index) {
    const std::string_view option = argv[index];
    if (option == "--help" || option == "-h") parsed.show_help = true;
    else if (option == "--version") parsed.show_version = true;
    else if (option == "--print-effective-config") parsed.print_effective_config = true;
    else if (option == "--reference") {
      parsed.run_spec.reference_path = RequireValue(argc, argv, index, option);
      have_reference = true;
    } else if (option == "--reference-index") {
      parsed.run_spec.reference_index_path = RequireValue(argc, argv, index, option);
    } else if (option == "--query") {
      parsed.run_spec.query_path = RequireValue(argc, argv, index, option);
      have_query = true;
    } else if (option == "--output-prefix") {
      if (!parsed.run_spec.outputs.empty()) throw CliError("--output-prefix cannot be combined with --output");
      parsed.run_spec.output_prefix = RequireValue(argc, argv, index, option);
      have_prefix = true;
    } else if (option == "--output") {
      if (have_prefix || have_formats) throw CliError("--output cannot be combined with --output-prefix or --formats");
      AddExplicitOutput(parsed.run_spec, RequireValue(argc, argv, index, option));
    } else if (option == "--work-dir") {
      parsed.run_spec.work_dir = RequireValue(argc, argv, index, option);
      have_work_dir = true;
    } else if (option == "--threads") {
      parsed.run_spec.threads = ParseInteger<std::uint32_t>(RequireValue(argc, argv, index, option), option);
    } else if (option == "--formats") {
      if (!parsed.run_spec.outputs.empty()) throw CliError("--formats cannot be combined with --output");
      parsed.run_spec.formats = ParseFormats(RequireValue(argc, argv, index, option));
      have_formats = true;
    } else if (option == "--progress") {
      parsed.run_spec.progress.mode = ParseProgressMode(RequireValue(argc, argv, index, option));
    } else if (option == "--progress-interval") {
      parsed.run_spec.progress.interval_seconds = ParseInteger<std::uint32_t>(RequireValue(argc, argv, index, option), option);
      parsed.run_spec.progress.interval_explicit = true;
    } else if (option == "--seed-mode") {
      parsed.run_spec.alignment.seed_mode = ParseSeedMode(RequireValue(argc, argv, index, option));
    } else if (option == "--selection-mode") {
      parsed.run_spec.alignment.selection = ParseAlignmentSelection(RequireValue(argc, argv, index, option));
    } else if (option == "--min-match") {
      parsed.run_spec.alignment.min_match = ParseInteger<Length>(RequireValue(argc, argv, index, option), option);
    } else if (option == "--smem-min-occurrences") {
      parsed.run_spec.alignment.smem_min_occurrences = ParseInteger<std::uint64_t>(RequireValue(argc, argv, index, option), option);
      have_smem_min_occurrences = true;
    } else if (option == "--max-gap") {
      parsed.run_spec.alignment.max_gap = ParseInteger<Length>(RequireValue(argc, argv, index, option), option);
    } else if (option == "--diag-diff") {
      parsed.run_spec.alignment.diag_diff = ParseInteger<Length>(RequireValue(argc, argv, index, option), option);
    } else if (option == "--diag-factor") {
      parsed.run_spec.alignment.diag_factor = ParseDouble(RequireValue(argc, argv, index, option), option);
    } else if (option == "--min-cluster") {
      parsed.run_spec.alignment.min_cluster = ParseInteger<Length>(RequireValue(argc, argv, index, option), option);
    } else if (option == "--break-length") {
      parsed.run_spec.alignment.break_length = ParseInteger<Length>(RequireValue(argc, argv, index, option), option);
    } else if (option == "--max-dp-cells") {
      parsed.run_spec.alignment.max_dp_cells = ParseInteger<std::uint64_t>(RequireValue(argc, argv, index, option), option);
    } else throw CliError("unknown option: " + std::string(option));
  }
  if (parsed.show_help || parsed.show_version) return parsed;
  if (!have_reference) throw CliError("missing required option --reference");
  if (!have_query) throw CliError("missing required option --query");
  if (!have_prefix && parsed.run_spec.outputs.empty()) {
    throw CliError("missing required option --output-prefix or --output");
  }
  if (!have_work_dir) throw CliError("missing required option --work-dir");
  if (have_smem_min_occurrences && parsed.run_spec.alignment.seed_mode != SeedMode::Smem) {
    throw CliError("--smem-min-occurrences requires --seed-mode smem");
  }
  ValidateRunSpec(parsed.run_spec);
  return parsed;
}

void ValidateRunSpec(const RunSpec& spec) {
#if RAMAG_USE_PAIRWISE_CORE
  if (spec.alignment.break_length != 200 || spec.alignment.max_dp_cells != 4000000) {
    throw CliError("pairwise core: legacy --break-length/--max-dp-cells overrides are not applicable");
  }
#endif
  ValidateInput(spec.reference_path, "reference");
  ValidateInput(spec.query_path, "query");
  if (!spec.reference_index_path.empty()) ValidateInput(spec.reference_index_path, "reference index");
  if (spec.output_prefix.empty()) throw CliError("output prefix must not be empty");
  ValidateWorkDir(spec.work_dir);
  ValidateOutputParent(spec.output_prefix);
  for (const auto& output : spec.outputs) ValidateOutputParent(output.path);
  if (spec.threads == 0) throw CliError("--threads must be at least 1");
  if (!spec.formats.Any()) throw CliError("at least one output format is required");
  ValidateProgress(spec.progress);
  if (spec.alignment.min_match == 0) throw CliError("--min-match must be at least 1");
  if (spec.alignment.smem_min_occurrences == 0) throw CliError("--smem-min-occurrences must be at least 1");
  if (spec.alignment.seed_mode != SeedMode::Smem && spec.alignment.smem_min_occurrences != 1) {
    throw CliError("non-default smem_min_occurrences requires seed mode smem");
  }
  if (spec.alignment.min_cluster == 0) throw CliError("--min-cluster must be at least 1");
  if (spec.alignment.break_length == 0) throw CliError("--break-length must be at least 1");
  if (spec.alignment.diag_factor < 0.0) throw CliError("--diag-factor must be non-negative");
  if (spec.alignment.max_dp_cells == 0) throw CliError("--max-dp-cells must be at least 1");
  if (spec.formats.delta) {
    const auto contains_whitespace = [](const std::filesystem::path& path) {
      const auto text = AbsoluteForDisplay(path);
      return std::any_of(text.begin(), text.end(), [](unsigned char value) {
        return std::isspace(value) != 0;
      });
    };
    if (contains_whitespace(spec.reference_path) || contains_whitespace(spec.query_path)) {
      throw CliError("delta format cannot encode whitespace in input paths; choose paths without whitespace or omit delta");
    }
  }
}

void ValidateIndexSpec(const IndexSpec& spec) {
  ValidateInput(spec.reference_path, "reference");
  if (spec.output_path.empty() || spec.output_path.extension() != ".sufidx") {
    throw CliError("index --output must end in .sufidx");
  }
  ValidateOutputParent(spec.output_path);
  ValidateWorkDir(spec.work_dir);
  if (spec.threads == 0) throw CliError("--threads must be at least 1");
  ValidateProgress(spec.progress);
}

std::string HelpText() {
  return R"(RaMA-G pairwise whole-genome aligner

Usage:
  ramag align --reference REF.fa[.gz] --query QUERY.fa[.gz]
              --output-prefix PREFIX --work-dir DIR [options]
  ramag align --reference REF.fa[.gz] --query QUERY.fa[.gz]
              --output RESULT.paf [--output RESULT.maf ...] --work-dir DIR
  ramag index --reference REF.fa[.gz] --output REF.sufidx
              --work-dir DIR [--threads N]

Alignment I/O:
  --reference PATH          Reference FASTA or gzip FASTA
  --reference-index PATH    Reusable index created by 'ramag index'
  --query PATH              Query FASTA or gzip FASTA
  --output PATH             Repeatable; suffix selects sam|paf|delta|maf|chain
  --output-prefix PATH      Legacy output prefix
  --formats LIST            Legacy comma list (default: sam,paf,delta)
  --work-dir PATH           Run-local workspace

Core options:
  --threads N               Worker budget (default: 1)
  --progress auto|on|off    Terminal progress (default: auto)
  --progress-interval N     Progress seconds, N >= 1 (default: 10)
  --seed-mode MODE          fast|mumreference|mum|smem|maxmatch (default: fast)
  --selection-mode MODE     all|one-to-one (default: all)
  --min-match N             Minimum exact seed length (default: 20)
  --smem-min-occurrences N  Minimum SMEM reference occurrences (default: 1)
  --print-effective-config  Validate, print effective configuration, and exit

Auditable baseline tuning:
  --max-gap N               Maximum seed-chain gap (default: 90)
  --diag-diff N             Diagonal tolerance (default: 5)
  --diag-factor X           Relative diagonal tolerance (default: 0.12)
  --min-cluster N           Minimum cluster span (default: 65)
  --break-length N          Chain split threshold (default: 200)
  --max-dp-cells N          Hard bound for scalar DP (default: 4000000)

Signals on POSIX systems:
  SIGUSR1 prints a progress snapshot. SIGINT/SIGTERM request safe interruption;
  interrupted runs publish no final output or completion marker and are not resumable.
)";
}

std::string VersionText() {
  const auto openmp = CurrentOpenMpRuntimeInfo();
  std::ostringstream output;
  output << "ramag " << RAMAG_VERSION << " (commit " << RAMAG_GIT_COMMIT
         << "; build " << RAMAG_BUILD_TYPE << "; compiler " << RAMAG_CXX_COMPILER
         << "; openmp " << (openmp.enabled ? "enabled" : "disabled")
         << "; runtime " << openmp.runtime << "; specification_date "
         << openmp.specification_date << "; max_threads " << openmp.max_threads
         << "; max_active_levels " << openmp.max_active_levels
         << "; configured_requested_threads " << openmp.configured_requested_threads
         << "; sufkit_divsufsort_openmp "
         << (RAMAG_SUFKIT_DIVSUFSORT_OPENMP ? "enabled" : "disabled") << ')';
  return output.str();
}

std::string EffectiveConfigText(const RunSpec& spec) {
  const auto openmp = CurrentOpenMpRuntimeInfo();
  std::ostringstream output;
#if RAMAG_USE_PAIRWISE_CORE
  output << "alignment_core=pairwise\n"
         << "pairwise_source_commit=7d08359e0df7f7e6ffcfe67217c3399761cb2129\n"
         << "extension_scoring=scaled-HOXD70;gap-open=40;gap-extend=3\n"
         << "break_length_applicability=fixed-pairwise-200\n"
         << "max_dp_cells_applicability=legacy-only\n"
         << "selection_contract=pairwise-reference-query-dp-intersection\n";
#endif
  output << "command=align\n"
         << "reference=" << AbsoluteForDisplay(spec.reference_path) << '\n'
         << "reference_index=" << AbsoluteForDisplay(spec.reference_index_path) << '\n'
         << "index_action="
         << (spec.reference_index_path.empty() ? "built" : "loaded") << '\n'
         << "query=" << AbsoluteForDisplay(spec.query_path) << '\n'
         << "output_prefix=" << AbsoluteForDisplay(spec.output_prefix) << '\n';
  for (const auto& request : spec.outputs) output << "output=" << AbsoluteForDisplay(request.path) << '\n';
  output << "work_dir=" << AbsoluteForDisplay(spec.work_dir) << '\n'
         << "threads=" << spec.threads << '\n'
         << "openmp_enabled=" << (openmp.enabled ? "true" : "false") << '\n'
         << "openmp_runtime=" << openmp.runtime << '\n'
         << "openmp_specification_date=" << openmp.specification_date << '\n'
         << "openmp_runtime_max_threads=" << openmp.max_threads << '\n'
         << "openmp_thread_limit=" << openmp.thread_limit << '\n'
         << "openmp_dynamic=" << openmp.dynamic << '\n'
         << "openmp_max_active_levels=" << openmp.max_active_levels << '\n'
         << "openmp_requested_threads=" << spec.threads << '\n'
         << "openmp_configured_requested_threads=" << openmp.configured_requested_threads << '\n'
         << "sufkit_divsufsort_openmp=" << (RAMAG_SUFKIT_DIVSUFSORT_OPENMP ? "true" : "false") << '\n'
         << "formats=" << spec.formats.ToString() << '\n'
         << "progress=" << ProgressModeName(spec.progress.mode) << '\n'
         << "progress_interval=" << spec.progress.interval_seconds << '\n'
         << "seed_mode=" << SeedModeName(spec.alignment.seed_mode) << '\n'
         << "selection_mode=" << AlignmentSelectionName(spec.alignment.selection) << '\n'
         << "min_match=" << spec.alignment.min_match << '\n'
         << "smem_min_occurrences=" << spec.alignment.smem_min_occurrences << '\n'
         << "max_gap=" << spec.alignment.max_gap << '\n'
         << "diag_diff=" << spec.alignment.diag_diff << '\n'
         << "diag_factor=" << spec.alignment.diag_factor << '\n'
         << "min_cluster=" << spec.alignment.min_cluster << '\n'
         << "break_length=" << spec.alignment.break_length << '\n'
         << "max_dp_cells=" << spec.alignment.max_dp_cells << '\n';
  return output.str();
}

std::string EffectiveIndexConfigText(const IndexSpec& spec) {
  std::ostringstream output;
  output << "command=index\n"
         << "reference=" << AbsoluteForDisplay(spec.reference_path) << '\n'
         << "output=" << AbsoluteForDisplay(spec.output_path) << '\n'
         << "work_dir=" << AbsoluteForDisplay(spec.work_dir) << '\n'
         << "threads=" << spec.threads << '\n'
         << "progress=" << ProgressModeName(spec.progress.mode) << '\n'
         << "progress_interval=" << spec.progress.interval_seconds << '\n';
  return output.str();
}

}  // namespace ramag

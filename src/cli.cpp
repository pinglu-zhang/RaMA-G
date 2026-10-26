#include "ramag/cli.hpp"
#include "ramag/centromere_direction.hpp"
#include "ramag/centromere_blocks.hpp"
#include "ramag/pairwise_core.hpp"

#include <algorithm>
#include <atomic>
#include <cctype>
#include <charconv>
#include <cmath>
#include <cerrno>
#include <cstring>
#include <limits>
#include <fstream>
#include <set>
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

  const bool batch = first == "batch";
  if (first != "align" && !batch) {
    throw CliError("expected subcommand 'align', 'index' or 'batch'; got: " + std::string(first));
  }
  parsed.command = batch ? CommandKind::Batch : CommandKind::Align;
  if (batch) parsed.run_spec.formats = {false, true, false, false, false};
  bool have_reference = false;
  bool have_query = false;
  bool have_prefix = false;
  bool have_work_dir = false;
  bool have_formats = false;
  bool have_smem_min_occurrences = false;
  std::string centromere_incompatible;
  for (int index = 2; index < argc; ++index) {
    const std::string_view option = argv[index];
    if (option == "--reference-index" || option == "--save" || option == "--seed-mode" ||
        option == "--selection-mode" || option == "--min-match" || option == "--smem-min-occurrences" ||
        option == "--max-gap" || option == "--diag-diff" || option == "--diag-factor" ||
        option == "--min-cluster" || option == "--break-length" || option == "--max-dp-cells") {
      if (centromere_incompatible.empty()) centromere_incompatible = option;
    }
    if (option == "--help" || option == "-h") parsed.show_help = true;
    else if (option == "--version") parsed.show_version = true;
    else if (option == "--print-effective-config") parsed.print_effective_config = true;
    else if (option == "--reference") {
      parsed.run_spec.reference_path = RequireValue(argc, argv, index, option);
      have_reference = true;
    } else if (option == "--reference-index") {
      parsed.run_spec.reference_index_path = RequireValue(argc, argv, index, option);
    } else if (option == "--save") {
      if (!parsed.run_spec.save_index_path.empty()) throw CliError("--save may appear only once");
      parsed.run_spec.save_index_path = RequireValue(argc, argv, index, option);
      if (parsed.run_spec.save_index_path.empty()) throw CliError("--save requires an index file path");
    } else if (option == "--query") {
      parsed.run_spec.query_path = RequireValue(argc, argv, index, option);
      if (batch) {
        auto stem = parsed.run_spec.query_path.filename();
        if (stem.extension() == ".gz") stem = stem.stem();
        if (stem.extension() == ".fa" || stem.extension() == ".fna" || stem.extension() == ".fasta") stem = stem.stem();
        parsed.batch_spec.queries.push_back({stem.string(), parsed.run_spec.query_path});
      }
      have_query = true;
    } else if (option == "-c" || option == "--centromere") {
      if (batch) throw CliError("centromere options are supported only by align");
      if (parsed.run_spec.centromere.whole_input) throw CliError("--centromere may appear only once");
      parsed.run_spec.centromere.whole_input = true;
    } else if (option == "--ref-cen" || option == "--reference-centromeres" ||
               option == "--query-cen" || option == "--query-centromeres") {
      if (batch) throw CliError("centromere options are supported only by align");
      auto& path = (option == "--ref-cen" || option == "--reference-centromeres")
                       ? parsed.run_spec.centromere.reference_bed
                       : parsed.run_spec.centromere.query_bed;
      if (!path.empty()) throw CliError(std::string(option) + " may appear only once (including its alias)");
      path = RequireValue(argc, argv, index, option);
      if (path.empty()) throw CliError(std::string(option) + " requires a BED path");
    } else if (batch && option == "--seqfile") {
      if (!parsed.batch_spec.seqfile.empty()) throw CliError("--seqfile may appear only once");
      parsed.batch_spec.seqfile = RequireValue(argc, argv, index, option);
    } else if (batch && option == "--output-dir") {
      if (!parsed.batch_spec.output_dir.empty()) throw CliError("--output-dir may appear only once");
      parsed.batch_spec.output_dir = RequireValue(argc, argv, index, option);
    } else if (option == "--output-prefix") {
      if (batch) throw CliError("batch uses --output-dir, not --output-prefix");
      if (!parsed.run_spec.outputs.empty()) throw CliError("--output-prefix cannot be combined with --output");
      parsed.run_spec.output_prefix = RequireValue(argc, argv, index, option);
      have_prefix = true;
    } else if (option == "--output") {
      if (batch) throw CliError("batch uses --output-dir and --formats, not --output");
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
  if (parsed.run_spec.centromere.whole_input && !centromere_incompatible.empty())
    throw CliError("--centromere does not accept genome-only option " + centromere_incompatible);
  if (parsed.run_spec.centromere.whole_input && CentromereCliAnchorsOnly() &&
      !parsed.print_effective_config && (have_formats || !parsed.run_spec.outputs.empty()))
    throw CliError("this build exports centromere anchors only; use --output-prefix, not --output/--formats");
  if (!have_reference) throw CliError("missing required option --reference");
  if (have_smem_min_occurrences && parsed.run_spec.alignment.seed_mode != SeedMode::Smem) {
    throw CliError("--smem-min-occurrences requires --seed-mode smem");
  }
  if (batch) {
    auto& spec = parsed.batch_spec;
    if (!spec.seqfile.empty()) {
      if (have_query) throw CliError("--query and --seqfile are mutually exclusive");
      std::ifstream input(spec.seqfile);
      if (!input) throw CliError("cannot read seqfile: " + spec.seqfile.string());
      std::string line;
      std::size_t line_number = 0;
      while (std::getline(input, line)) {
        ++line_number;
        const auto begin = line.find_first_not_of(" \t\r");
        if (begin == std::string::npos || line[begin] == '#') continue;
        const auto split = line.find_first_of(" \t", begin);
        const auto path_begin = split == std::string::npos ? split : line.find_first_not_of(" \t", split);
        if (path_begin == std::string::npos || path_begin > line.find_last_not_of(" \t\r"))
          throw CliError("seqfile requires name and FASTA path at line " + std::to_string(line_number));
        std::filesystem::path path(line.substr(path_begin, line.find_last_not_of(" \t\r") - path_begin + 1));
        if (path.is_relative()) path = spec.seqfile.parent_path() / path;
        spec.queries.push_back({line.substr(begin, split - begin), path.lexically_normal()});
      }
      if (input.bad()) throw CliError("error reading seqfile");
    }
    if (!have_work_dir) parsed.run_spec.work_dir = spec.output_dir / ".ramag-work";
    parsed.run_spec.query_path.clear();
    spec.common = parsed.run_spec;
    ValidateBatchSpec(spec);
    return parsed;
  }
  if (!have_query) throw CliError("missing required option --query");
  if (!have_prefix && parsed.run_spec.outputs.empty()) {
    throw CliError("missing required option --output-prefix or --output");
  }
  if (!have_work_dir) {
    if (parsed.run_spec.centromere.whole_input)
      parsed.run_spec.work_dir = parsed.run_spec.output_prefix.parent_path() / ".ramag-work";
    else throw CliError("missing required option --work-dir");
  }
  if (have_smem_min_occurrences && parsed.run_spec.alignment.seed_mode != SeedMode::Smem) {
    throw CliError("--smem-min-occurrences requires --seed-mode smem");
  }
  ValidateRunSpec(parsed.run_spec);
  return parsed;
}

void ValidateRunSpec(const RunSpec& spec) {
  if (spec.centromere.whole_input &&
      (!spec.reference_index_path.empty() || !spec.save_index_path.empty()))
    throw CliError("--centromere does not support index loading or saving");
  try {
    (void)LoadCentromereAnnotations(spec.centromere);
  } catch (const CentromereInputError& error) {
    throw CliError(error.what());
  }
  if (!spec.save_index_path.empty()) {
    if (!spec.reference_index_path.empty()) throw CliError("--save and --reference-index are mutually exclusive");
    ValidateOutputParent(spec.save_index_path);
    std::vector<std::filesystem::path> protected_paths{spec.reference_path, spec.query_path};
    if (!spec.centromere.reference_bed.empty()) protected_paths.push_back(spec.centromere.reference_bed);
    if (!spec.centromere.query_bed.empty()) protected_paths.push_back(spec.centromere.query_bed);
    for (const auto* suffix : {".sam", ".paf", ".delta", ".maf", ".chain"}) {
      protected_paths.emplace_back(spec.output_prefix.string() + suffix);
    }
    for (const auto& output : spec.outputs) protected_paths.push_back(output.path);
    for (const auto* suffix : {""}) {
      const std::filesystem::path target(spec.save_index_path.string() + suffix);
      if (std::filesystem::symlink_status(target).type() != std::filesystem::file_type::not_found) {
        throw CliError("refusing to overwrite index artifact: " + target.string());
      }
      const auto normalized = std::filesystem::weakly_canonical(target);
      for (const auto& path : protected_paths) {
        if (normalized == std::filesystem::weakly_canonical(path)) {
          throw CliError("saved index path conflicts with input or alignment output: " + target.string());
        }
      }
    }
  }
  if (spec.alignment.break_length != 200 || spec.alignment.max_dp_cells != 4000000) {
    throw CliError("pairwise core: legacy --break-length/--max-dp-cells overrides are not applicable");
  }
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
  if (spec.formats.delta && !(spec.centromere.whole_input && CentromereCliAnchorsOnly())) {
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

void ValidateIndexSpec(const IndexSpec& spec, bool require_suffix) {
  ValidateInput(spec.reference_path, "reference");
  if (spec.output_path.empty() || (require_suffix && spec.output_path.extension() != ".sufidx")) {
    throw CliError("index --output must end in .sufidx");
  }
  ValidateOutputParent(spec.output_path);
  ValidateWorkDir(spec.work_dir);
  if (spec.threads == 0) throw CliError("--threads must be at least 1");
  ValidateProgress(spec.progress);
}

RunSpec BatchQuerySpec(const BatchSpec& spec, const BatchQuery& query) {
  auto run = spec.common;
  run.query_path = query.path;
  run.save_index_path.clear();
  run.output_prefix = spec.output_dir / query.name / "alignment";
  run.work_dir = spec.common.work_dir / "queries" / query.name;
  run.outputs.clear();
  return run;
}

void ValidateBatchSpec(const BatchSpec& spec) {
  if (spec.common.centromere.Mode() != CentromereMode::Genome) {
    throw CliError("centromere options are supported only by align");
  }
  if (spec.output_dir.empty()) throw CliError("batch requires --output-dir");
  if (spec.queries.empty()) throw CliError("batch requires at least one query");
  if (!spec.common.outputs.empty() || !spec.common.output_prefix.empty()) throw CliError("batch output paths are generated from --output-dir");
  std::set<std::string> names;
  std::vector<std::filesystem::path> targets{spec.output_dir / "batch.tsv"};
  const auto ancestor = [](const std::filesystem::path& a, const std::filesystem::path& b) {
    auto i = a.begin(); auto j = b.begin();
    for (; i != a.end() && j != b.end() && *i == *j; ++i, ++j) {}
    return i == a.end();
  };
  for (const auto& query : spec.queries) {
    const auto& name = query.name;
    if (name.empty() || name == "." || name == ".." || name == ".ramag-work" || name == "batch.tsv" ||
        name.find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789._-") != std::string::npos)
      throw CliError("invalid or reserved batch query name: " + name);
    if (!names.insert(name).second) throw CliError("duplicate batch query name: " + name + "; use --seqfile with unique names");
    if (query.path.empty()) throw CliError("empty batch query path");
    auto run = BatchQuerySpec(spec, query);
    // Validate shared settings without opening an individual query: missing or
    // malformed query input belongs to that item's result, not batch startup.
    run.query_path = spec.common.reference_path;
    run.save_index_path = spec.common.save_index_path;
    ValidateRunSpec(run);
    if (spec.common.formats.delta && query.path.string().find_first_of(" \t\n\r") != std::string::npos)
      throw CliError("delta format cannot encode whitespace in batch query paths");
    const auto prefix = run.output_prefix.string();
    if (run.formats.sam) targets.emplace_back(prefix + ".sam");
    if (run.formats.paf) targets.emplace_back(prefix + ".paf");
    if (run.formats.delta) targets.emplace_back(prefix + ".delta");
    if (run.formats.maf) targets.emplace_back(prefix + ".maf");
    if (run.formats.chain) targets.emplace_back(prefix + ".chain");
  }
  if (!spec.common.save_index_path.empty()) targets.push_back(spec.common.save_index_path);
  std::set<std::filesystem::path> unique;
  std::vector<std::filesystem::path> inputs{spec.common.reference_path};
  if (!spec.common.reference_index_path.empty()) inputs.push_back(spec.common.reference_index_path);
  if (!spec.seqfile.empty()) inputs.push_back(spec.seqfile);
  for (const auto& query : spec.queries) inputs.push_back(query.path);
  for (const auto& target : targets) {
    ValidateOutputParent(target);
    if (std::filesystem::symlink_status(target).type() != std::filesystem::file_type::not_found)
      throw CliError("batch output already exists: " + target.string());
    const auto normalized = std::filesystem::weakly_canonical(std::filesystem::absolute(target));
    for (const auto& previous : unique)
      if (ancestor(previous, normalized) || ancestor(normalized, previous)) throw CliError("batch output paths conflict");
    unique.insert(normalized);
    for (const auto& input : inputs) {
      const auto in = std::filesystem::weakly_canonical(std::filesystem::absolute(input));
      if (ancestor(normalized, in) || ancestor(in, normalized)) throw CliError("batch output conflicts with input: " + target.string());
    }
    const auto work = std::filesystem::weakly_canonical(std::filesystem::absolute(spec.common.work_dir));
    if (ancestor(normalized, work) || ancestor(work, normalized)) throw CliError("batch output conflicts with work directory");
  }
}

std::string EffectiveBatchConfigText(const BatchSpec& spec) {
  auto run = BatchQuerySpec(spec, spec.queries.front());
  run.save_index_path = spec.common.save_index_path;
  std::string common = EffectiveConfigText(run);
  const auto command = common.find("command=align\n");
  if (command != std::string::npos) common.replace(command, 14, "command=batch\n");
  std::ostringstream out;
  out << common << "batch_output_dir=" << AbsoluteForDisplay(spec.output_dir)
      << "\nbatch_query_count=" << spec.queries.size() << "\nbatch_query_execution=serial\nbatch_failure_policy=continue\n";
  for (std::size_t i = 0; i < spec.queries.size(); ++i) {
    const auto& q = spec.queries[i];
    out << "batch_query." << i + 1 << ".name=" << q.name << '\n'
        << "batch_query." << i + 1 << ".path=" << AbsoluteForDisplay(q.path) << '\n'
        << "batch_query." << i + 1 << ".output_prefix=" << AbsoluteForDisplay(BatchQuerySpec(spec, q).output_prefix) << '\n';
  }
  return out.str();
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
  ramag batch --reference REF.fa[.gz] [--reference-index REF.sufidx]
              --query QUERY.fa[.gz] [--query OTHER.fa[.gz] ...]
              --output-dir DIR [--threads N]
  ramag batch --reference REF.fa[.gz] --seqfile queries.txt --output-dir DIR

Batch options:
  --query PATH              Repeatable; mutually exclusive with --seqfile
  --seqfile PATH            Query name and FASTA path per line; relative to seqfile
  --output-dir PATH         Per-query NAME/alignment outputs and batch.tsv
  --work-dir PATH           Optional (default: OUTPUT_DIR/.ramag-work)
  --formats LIST            sam,paf,delta,maf,chain (batch default: paf)
  Shared alignment options below also apply. Queries run serially, reusing one
  resident index. Query failures are recorded and later queries continue.

Alignment I/O:
  --reference PATH          Reference FASTA or gzip FASTA
  --reference-index PATH    Reusable index created by 'ramag index'
  --save PATH               Save the newly built index, then align; no overwrite
  --query PATH              Query FASTA or gzip FASTA
  --output PATH             Repeatable; suffix selects sam|paf|delta|maf|chain
  --output-prefix PATH      Legacy output prefix
  --formats LIST            Legacy comma list (default: sam,paf,delta)
  --work-dir PATH           Run-local workspace

Centromere input (align only; source-derived RaMA + WFA2 on Linux):
  -c, --centromere          Process exactly one nonempty sequence per input
  --ref-cen BED             Centromere intervals in the full reference genome
  --query-cen BED           Centromere intervals in the full query genome
  Long aliases: --reference-centromeres, --query-centromeres
  Both BEDs are required together; cannot combine with -c. BED is zero-based,
  half-open, with >=3 columns and optional names. No automatic region detection.
  --print-effective-config validates this interface without running alignment.
  -c uses supplied orientation, fixed RaMA rules, and deterministic N mismatches.
  Default development build: anchors only; use --output-prefix for two TSV files.
  Check --print-effective-config for centromere_stage; no WFA in anchors-only mode.
  With -c only, omitted --work-dir defaults to .ramag-work beside the output.
  Experimental direction builds export candidate evidence, not confirmed inversions.
  Genome tuning/index options are incompatible with -c.
  BED mixed alignment remains unsupported; it does not fall back to genome mode.

Core options:
  --threads N               Worker budget (default: 1)
  --progress auto|on|off    Terminal progress (default: auto)
  --progress-interval N     Progress seconds, N >= 1 (default: 10)
  --seed-mode MODE          fast|mumreference|mum|smem|maxmatch (default: fast)
  --selection-mode MODE     all|one-to-one (default: all)
  --min-match N             Minimum exact seed length (default: 20)
  --smem-min-occurrences N  Minimum SMEM reference occurrences (default: 1)
  --print-effective-config  Validate, print effective configuration, and exit

Pairwise + KSW2 tuning:
  --max-gap N               Maximum seed-chain gap (default: 90)
  --diag-diff N             Diagonal tolerance (default: 5)
  --diag-factor X           Relative diagonal tolerance (default: 0.12)
  --min-cluster N           Minimum cluster span (default: 65)
  --break-length N          Compatibility only: accepts the default 200
  --max-dp-cells N          Compatibility only: accepts the default 4000000

Signals on POSIX systems:
  SIGUSR1 prints a progress snapshot. SIGINT/SIGTERM request safe interruption;
  interrupted alignments publish no successful alignment marker. Independently
  completed saved indexes are retained. There is no alignment resume interface.
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
  if (spec.centromere.whole_input) {
    const auto count_options=CentromereCliBlockOptions();
    const bool count_mode=count_options.planning_strategy==CentromereBlockStrategy::CountClusterDirectionsV8;
    const bool orientation_mode=count_options.planning_strategy==CentromereBlockStrategy::OrientationEvidenceV9;
    const bool evidence_mode=count_mode || orientation_mode;
    output << "command=align\nalignment_core=rama-centromere\ncentromere_mode=centromere-only\n"
           << "centromere_backend=" << (CentromereBackendAvailable() ? "available" : "not-built") << '\n'
           << "rama_source_commit=7bc89c4902bc9254f18cc4b05c4bf491a10bd7cb\n"
           << "wfa2_commit=99a5879a3a2972eff168cbcc6610b6cf02844ed0\n"
           << "centromere_stage=" << (orientation_mode ? "orientation-evidence-v9" : count_mode ? "count-direction-evidence" : CentromereCliBlocks() ? CentromereBlockStrategyName(CentromereCliBlockOptions().planning_strategy) : CentromereCliDirections() ? "direction-evidence" : CentromereCliAnchorsOnly() ? "anchors-only" : "complete-alignment") << '\n'
           << "wfa_execution=" << (CentromereCliAnchorsOnly() ? "disabled" : "enabled") << '\n'
           << "alignment_scope=" << (evidence_mode ? "direction-evidence-no-alignment" : CentromereCliBlocks() ? "fixed-direction-partitioned-pair" : CentromereCliAnchorsOnly() ? "sparse-exact-anchors" : "single-pair-end-to-end") << '\n'
           << "orientation=" << (orientation_mode ? "v6-plus-query-only-evidence-unresolved-reference-copy" : count_mode ? "joint-local-count-clusters" : CentromereCliBlocks() ?
              (CentromereCliBlockOptions().planning_strategy!=CentromereBlockStrategy::InitialDirectionBlocksV1 ? "joint-local-direction-and-correspondence" : "initial-both-then-fixed-blocks") :
              CentromereCliDirections() ? "independent-forward-and-reverse" : "supplied-forward") << '\n'
           << "centromere_inversion_detection=" << (evidence_mode ? "experimental-direction-evidence" : CentromereCliBlocks() ? "preliminary-fixed-blocks" : CentromereCliDirections() ? "experimental-anchor-evidence" : "disabled") << '\n'
           << "rare_match_max_count=100\n"
           << "extension_scoring=dual-affine;match=0;mismatch=3;open1=4;extend1=2;open2=12;extend2=1\n"
           << "wfa_memory=medium\nwfa_heuristic=wfadaptive:10:50:1\n"
           << "ambiguous_bases=N;hard-break-anchors;paired-N-is-mismatch\n"
           << "reference=" << AbsoluteForDisplay(spec.reference_path) << '\n'
           << "query=" << AbsoluteForDisplay(spec.query_path) << '\n'
           << "output_prefix=" << AbsoluteForDisplay(spec.output_prefix) << '\n'
           << "threads=" << spec.threads << "\nformats="
           << (orientation_mode ? "initial-anchors.tsv,root-candidates.tsv,blocks.tsv,planning-anchors.tsv,orientation-windows.tsv,orientation-evidence.tsv,orientation-sources.tsv,orientation-gaps.tsv,orientation-unknown.tsv" : count_mode ? "initial-anchors.tsv,root-candidates.tsv,blocks.tsv,planning-anchors.tsv,count-windows.tsv,count-nodes.tsv,count-candidates.tsv,count-gaps.tsv,supported-intervals.tsv,direction-regions.tsv" : CentromereCliDirections() ? "direction-anchors.tsv,direction-regions.tsv,direction-chains.tsv,direction-tasks.tsv,inversion-candidates.tsv" : CentromereCliAnchorsOnly() ? "initial-anchors.tsv,final-anchors.tsv" : spec.formats.ToString()) << '\n'
           << "work_dir=" << AbsoluteForDisplay(spec.work_dir) << '\n';
    if(orientation_mode && count_options.orientation_audit_known_blocks)
      output<<"orientation_audit_outputs=orientation-block-audit.tsv,orientation-known-blocks.tsv\n"
        <<"orientation_known_view="<<(count_options.repair_inherited_directions?
          "orientation-local-direction-repair-intervals.tsv;prior-direction-views-and-raw-blocks-are-provenance":count_options.review_known_block_counterchains?
          "orientation-known-counterchain-intervals.tsv;orientation-known-blocks.tsv-and-raw-blocks-are-provenance":
          "orientation-known-blocks.tsv;raw-blocks-are-provenance")<<'\n';
    if(orientation_mode && count_options.repair_inherited_directions)
      output<<"orientation_local_direction_repair_policy="<<(count_options.repair_drift_consistent_gaps?
          "single-pass-drift-consistent-root-competition-v2":"single-pass-frozen-root-competition-v1")<<'\n'
        <<"orientation_local_direction_repair_outputs=orientation-local-direction-repair-reviews.tsv,orientation-local-direction-repair-witnesses.tsv,orientation-local-direction-repair-intervals.tsv,orientation-local-direction-repair-changes.tsv\n";
    if(orientation_mode && count_options.review_known_block_counterchains)
      output<<"orientation_known_counterchain_policy=final-local-opposite-chain-retraction-v1\n"
        <<"orientation_known_counterchain_outputs=orientation-known-counterchain-reviews.tsv,orientation-known-counterchain-witnesses.tsv,orientation-known-counterchain-intervals.tsv\n";
    if(orientation_mode && count_options.recover_long_gap_ends)
      output<<"orientation_long_gap_outputs=orientation-long-gap-ends.tsv,orientation-long-gap-gaps.tsv\n";
    if(orientation_mode && count_options.recover_sparse_pivots)
      output<<"orientation_sparse_outputs=orientation-sparse-pivots.tsv,orientation-sparse-gaps.tsv,orientation-sparse-sources.tsv\n";
    if(orientation_mode && count_options.recover_short_gap_ends)
      output<<"orientation_short_gap_outputs=orientation-short-gap-ends.tsv,orientation-short-gap-gaps.tsv\n";
    if(orientation_mode && count_options.recover_mixed_gap_ends)
      output<<"orientation_mixed_gap_outputs=orientation-mixed-gap-ends.tsv,orientation-mixed-gap-gaps.tsv\n";
    if(orientation_mode && count_options.recover_sparse_gap_ends)
      output<<"orientation_sparse_gap_outputs=orientation-sparse-gap-ends.tsv,orientation-sparse-gap-gaps.tsv\n";
    if(orientation_mode && count_options.unify_real_anchor_ends)
      output<<"orientation_unified_outputs=orientation-unified-pivots.tsv,orientation-unified-proposals.tsv,orientation-unified-ends.tsv\n";
    if(evidence_mode)output<<"block_planning_strategy="<<CentromereBlockStrategyName(count_options.planning_strategy)<<"\ndirection_recursion=disabled\n"
      <<"count_window="<<count_options.count_window<<"\ncount_step="<<count_options.count_step
      <<"\ncount_min_anchors="<<count_options.count_min_anchors<<"\ncount_dominance_ratio="<<count_options.count_dominance_ratio
      <<"\ncount_cluster_gap="<<count_options.count_cluster_gap<<"\ncount_inferred_gap="<<count_options.count_inferred_gap
      <<"\nblock_planning_memory_bytes="<<count_options.planner_memory_limit_bytes<<'\n';
    if(orientation_mode)output<<"orientation_length_ratio="<<count_options.count_dominance_ratio
      <<"\norientation_connection_policy="<<(count_options.orientation_bridge_unknown_gaps?"same-direction-unknown-gap-v1":count_options.orientation_connected_gaps?"continuous-windows-known-pivots-v1":"single-window-v9")
      <<"\norientation_audit_policy="<<(count_options.orientation_audit_known_blocks?"opposite-window-old-gap-v1":"disabled")
      <<"\norientation_long_gap_policy="<<(count_options.recover_long_gap_ends?"bounded-pure-window-endpoints-v1":"disabled")
      <<"\norientation_sparse_policy="<<(count_options.recover_sparse_pivots?"frozen-sparse-root-pivots-v1":"disabled")
      <<"\norientation_short_gap_policy="<<(count_options.recover_short_gap_ends?"bounded-short-pure-window-ends-v1":"disabled")
      <<"\norientation_mixed_gap_policy="<<(count_options.recover_mixed_gap_ends?"bounded-opposite-direction-pure-window-ends-v1":"disabled")
      <<"\norientation_sparse_gap_policy="<<(count_options.recover_sparse_gap_ends?"bounded-frozen-sparse-root-gap-ends-v1":"disabled")
      <<"\norientation_scope="<<(count_options.orientation_audit_known_blocks?"v6-query-evidence-with-old-block-audit":"v6-unknown-query-only")
      <<"\nreference_copy=unresolved\n";
    if(orientation_mode && count_options.prefer_stronger_sparse_gap_end)
      output<<"orientation_sparse_gap_meeting_policy=stronger-frozen-accepted-root-v1\n";
    if(orientation_mode && count_options.unify_real_anchor_ends)
      output<<"orientation_unified_policy=frozen-real-anchor-ends-v1\n";
    if(orientation_mode && count_options.review_local_competitive_gaps)
      output<<"orientation_local_review_policy=same-direction-original-corridor-v1\n"
        <<"orientation_local_review_outputs=orientation-local-reviews.tsv,orientation-local-gaps.tsv\n";
    if(orientation_mode && count_options.review_collinear_fragment_connections)
      output<<"orientation_collinear_review_policy=two-flank-collinear-witness-v1\n"
        <<"orientation_collinear_review_outputs=orientation-collinear-reviews.tsv,orientation-collinear-counters.tsv,orientation-collinear-witnesses.tsv,orientation-collinear-gaps.tsv\n";
    if(orientation_mode && count_options.review_collinear_window_gaps)
      output<<"orientation_collinear_window_policy=positive-window-break-full-audit-v1\n";
    if(orientation_mode && count_options.review_neutral_window_gaps)
      output<<"orientation_neutral_window_policy=bounded-uncontested-window-full-review-v1\n";
    if(orientation_mode && count_options.review_alternate_original_endpoints)
      output<<"orientation_alternate_endpoint_policy=bounded-original-endpoint-full-review-v1\n"
        <<"orientation_alternate_endpoint_outputs=orientation-alternate-reviews.tsv,orientation-alternate-counters.tsv,orientation-alternate-witnesses.tsv,orientation-alternate-gaps.tsv\n";
    if(orientation_mode && count_options.review_isolated_opposite_anchors)
      output<<"orientation_isolated_review_policy=bounded-contextual-opposite-v1\n"
        <<"orientation_isolated_review_outputs=orientation-isolated-reviews.tsv,orientation-opposite-reviews.tsv,orientation-isolated-gaps.tsv\n";
    if (CentromereCliBlocks()) {
      const auto budgets=CentromereCliBlockOptions();
      output << "wfa_gap_memory_bytes=" << budgets.wfa_memory_limit_bytes
             << "\nblock_planning_strategy=" << CentromereBlockStrategyName(budgets.planning_strategy)
             << "\nblock_planning_memory_bytes=" << budgets.planner_memory_limit_bytes
             << "\nunknown_direction=" << (budgets.planning_strategy!=CentromereBlockStrategy::InitialDirectionBlocksV1 ? "unresolved-no-wfa" : "background-fallback")
             << "\nwfa_gap_seconds=" << budgets.wfa_seconds_limit
             << "\nblock_alignment_seconds=" << budgets.alignment_seconds_limit
             << "\nwfa_workers=1\npartial_exit_code=9\npartial_suffix=.partial\n";
    }
    if (CentromereCliDirections())
      output << "direction_objective=effective-collinear-length-v1\n"
             << "direction_search_strategy=" << DirectionStrategyName(CentromereCliDirectionStrategy()) << '\n'
             << "direction_supplemental_rare_levels="
             << (CentromereCliDirectionStrategy()==DirectionSearchStrategy::LegacyMultiscale ? "true" : "false") << '\n'
             << "direction_windows="
             << (CentromereCliDirectionStrategy()==DirectionSearchStrategy::LegacyMultiscale ?
               "80bp-doubling-half-step-reflection-symmetric" : "paired-chain-gaps-and-opposite-support-only") << '\n'
             << "direction_recursion=both-directions\n"
             << "direction_actual_workers=1\n"
             << "direction_evidence_memory_guard_bytes=100000000000\n"
             << "direction_breakpoints=uncertainty-intervals-not-base-resolved\n";
    for (const auto& request:spec.outputs)
      output << "output=" << AbsoluteForDisplay(request.path) << '\n';
    return output.str();
  }
  output << PairwiseConfigurationText(spec.alignment);
  output << "alignment_core="
         << (spec.centromere.Mode() == CentromereMode::Genome ? "pairwise" : "centromere-backend-pending") << '\n'
         << "pairwise_source_commit=7d08359e0df7f7e6ffcfe67217c3399761cb2129\n"
         << "extension_scoring=scaled-HOXD70;gap-open=40;gap-extend=3\n"
         << "break_length_applicability=fixed-pairwise-200\n"
         << "max_dp_cells_applicability=legacy-only\n"
         << "selection_contract=pairwise-reference-query-dp-intersection\n";
  output << "command=align\n"
         << "reference=" << AbsoluteForDisplay(spec.reference_path) << '\n'
         << "reference_index=" << AbsoluteForDisplay(spec.reference_index_path) << '\n'
         << "save=" << AbsoluteForDisplay(spec.save_index_path) << '\n'
         << "index_persistence=" << (spec.save_index_path.empty() ? "not-requested" : "requested") << '\n'
         << "lcp_encoding_policy=" << (spec.reference_index_path.empty() ? RAMAG_INTERNAL_LCP_ENCODING : "from-index-file") << '\n'
         << "index_action="
         << (spec.reference_index_path.empty() ? "built" : "loaded") << '\n'
         << "query=" << AbsoluteForDisplay(spec.query_path) << '\n'
         << "output_prefix=" << AbsoluteForDisplay(spec.output_prefix) << '\n';
  for (const auto& request : spec.outputs) output << "output=" << AbsoluteForDisplay(request.path) << '\n';
  if (spec.centromere.Mode() != CentromereMode::Genome) {
    const auto annotations = LoadCentromereAnnotations(spec.centromere);
    output << "centromere_mode=" << CentromereModeName(spec.centromere.Mode()) << '\n'
           << "alignment_scope=" << (spec.centromere.whole_input ? "all-input-records" : "whole-genome-with-specialised-regions") << '\n'
           << "reference_centromeres=" << AbsoluteForDisplay(spec.centromere.reference_bed) << '\n'
           << "query_centromeres=" << AbsoluteForDisplay(spec.centromere.query_bed) << '\n'
           << "reference_centromere_regions=" << annotations.reference.size() << '\n'
           << "query_centromere_regions=" << annotations.query.size() << '\n'
           << "centromere_pairing=not-inferred-from-BED-order-or-names\n"
           << "centromere_backend=not-integrated\n"
           << "centromere_inversion_detection=not-integrated\n";
  }
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

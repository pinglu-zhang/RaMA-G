#include "ramag/pipeline.hpp"
#include "ramag/pairwise_core.hpp"

#include "ramag/alignment.hpp"
#include "ramag/fasta.hpp"
#include "ramag/reference_index.hpp"
#include "ramag/runtime.hpp"
#include "ramag/logging.hpp"
#include <map>
#include "ramag/sufkit_adapter.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdlib>
#include <ctime>
#include <exception>
#include <fstream>
#include <functional>
#include <iterator>
#include <new>
#include <random>
#include <sstream>
#include <string_view>
#include <system_error>
#include <utility>
#include <vector>

#if RAMAG_OPENMP_ENABLED
#include <omp.h>
#endif

namespace ramag {
namespace {

using Clock = std::chrono::steady_clock;

double SecondsBetween(Clock::time_point begin, Clock::time_point end) {
  return std::chrono::duration<double>(end - begin).count();
}

std::filesystem::path TemporarySibling(const std::filesystem::path& final_path,
                                       std::string_view nonce) {
  return final_path.string() + ".tmp." + std::string(nonce);
}

std::string MakeNonce() {
  const auto ticks = Clock::now().time_since_epoch().count();
  std::random_device random;
  std::ostringstream value;
  value << std::hex << ticks << '-' << random();
  return value.str();
}

std::string AbsolutePathString(const std::filesystem::path& path) {
  if (path.empty()) {
    return {};
  }
  std::error_code error;
  const auto absolute = std::filesystem::absolute(path, error);
  return (error ? path : absolute.lexically_normal()).string();
}

void EnsureDirectory(const std::filesystem::path& path,
                     std::string_view description) {
  if (path.empty()) {
    return;
  }
  std::error_code error;
  std::filesystem::create_directories(path, error);
  if (error || !std::filesystem::is_directory(path)) {
    throw WriterError("cannot create " + std::string(description) + ": " +
                      path.string() + (error ? " (" + error.message() + ")" : ""));
  }
}

void RemoveIfPresent(const std::filesystem::path& path) noexcept {
  std::error_code ignored;
  std::filesystem::remove(path, ignored);
}

void PublishNewFile(const std::filesystem::path& temporary,
                    const std::filesystem::path& final_path) {
  std::error_code error;
  // temporary is a sibling of final_path, so an atomic hard-link publication
  // makes the already-closed inode visible without POSIX rename's overwrite
  // race. The temporary link remains until the transaction commits.
  std::filesystem::create_hard_link(temporary, final_path, error);
  if (error) {
    throw WriterError("cannot atomically publish new artifact " +
                      final_path.string() + " without overwriting: " +
                      error.message());
  }
}

class TemporaryFiles {
 public:
  void Add(std::filesystem::path path) { paths_.push_back(std::move(path)); }
  ~TemporaryFiles() {
    for (const auto& path : paths_) {
      RemoveIfPresent(path);
    }
  }

 private:
  std::vector<std::filesystem::path> paths_;
};

class PublishedFiles {
 public:
  void Add(std::filesystem::path temporary, std::filesystem::path final_path) {
    paths_.emplace_back(std::move(temporary), std::move(final_path));
  }
  void Commit() noexcept { committed_ = true; }
  ~PublishedFiles() {
    if (!committed_) {
      for (const auto& [temporary, final_path] : paths_) {
        std::error_code error;
        // Only remove a final name that still identifies the inode published
        // by this run. A concurrently replaced path is not ours to delete.
        if (std::filesystem::equivalent(temporary, final_path, error) && !error) {
          RemoveIfPresent(final_path);
        }
      }
    }
  }

 private:
  std::vector<std::pair<std::filesystem::path, std::filesystem::path>> paths_;
  bool committed_{false};
};

void WriteFile(const std::filesystem::path& path,
               const std::function<void(std::ostream&)>& writer) {
  std::ofstream output(path, std::ios::binary | std::ios::trunc);
  if (!output) {
    throw WriterError("cannot open temporary output: " + path.string());
  }
  writer(output);
  output.flush();
  if (!output) {
    throw WriterError("failed to flush temporary output: " + path.string());
  }
  output.close();
  if (!output) {
    throw WriterError("failed to close temporary output: " + path.string());
  }
}

void MaybeInjectWriterFailure(std::string_view format) {
  const char* requested = std::getenv("RAMAG_TEST_FAIL_WRITER");
  if (requested != nullptr && std::string_view(requested) == format) {
    throw WriterError("injected " + std::string(format) + " writer failure");
  }
}

void MaybeInjectValidatorFailure(std::string_view format) {
  const char* requested = std::getenv("RAMAG_TEST_FAIL_VALIDATOR");
  if (requested != nullptr && std::string_view(requested) == format) {
    throw WriterError("injected " + std::string(format) +
                      " validator failure");
  }
}

void MaybeInjectPrecomputeFailure() {
  if (std::getenv("RAMAG_TEST_FAIL_DEPENDENCY") != nullptr) {
    throw DependencyError("injected dependency/index failure");
  }
  if (std::getenv("RAMAG_TEST_FAIL_RESOURCE") != nullptr) {
    throw std::bad_alloc{};
  }
}

void MaybeInjectInputFailure(std::string_view role) {
  const char* requested = std::getenv("RAMAG_TEST_FAIL_INPUT");
  if (requested == nullptr) {
    return;
  }
  const std::string_view selection{requested};
  if (selection == role || selection == "both") {
    throw FastaError("injected " + std::string(role) + " input failure");
  }
}

struct ParallelFastaInputs {
  std::array<FastaData, 2> results;
  std::array<std::exception_ptr, 2> failures;
  std::array<int, 2> worker_ids{-1, -1};
  std::uint32_t requested_workers{1};
  std::uint32_t actual_workers{1};
  std::string route{"serial-reference-query-v1"};
};

ParallelFastaInputs ReadFastaInputs(const RunSpec& spec) {
  ParallelFastaInputs inputs;
  inputs.requested_workers = spec.threads > 1 ? 2U : 1U;

  const auto read_role = [&](std::size_t role_index) noexcept {
    try {
#if RAMAG_OPENMP_ENABLED
      inputs.worker_ids[role_index] =
          omp_in_parallel() != 0 ? omp_get_thread_num() : 0;
#else
      inputs.worker_ids[role_index] = 0;
#endif
      const bool is_reference = role_index == 0;
      const std::string_view role = is_reference ? "reference" : "query";
      MaybeInjectInputFailure(role);
      inputs.results[role_index] = ReadFasta(is_reference ? spec.reference_path : spec.query_path);
    } catch (...) {
      // No exception may cross an OpenMP structured block. Each role owns a
      // distinct result and failure slot; failures are rethrown below in the
      // stable reference-before-query order.
      inputs.failures[role_index] = std::current_exception();
    }
  };

#if RAMAG_OPENMP_ENABLED
  if (inputs.requested_workers == 2U && omp_in_parallel() == 0) {
    int observed_workers = 1;
#pragma omp parallel num_threads(2) default(none) \
    shared(inputs, observed_workers, read_role)
    {
#pragma omp single
      { observed_workers = omp_get_num_threads(); }
#pragma omp sections
      {
#pragma omp section
        { read_role(0); }
#pragma omp section
        { read_role(1); }
      }
    }
    const bool two_active_workers =
        observed_workers >= 2 && inputs.worker_ids[0] >= 0 &&
        inputs.worker_ids[1] >= 0 &&
        inputs.worker_ids[0] != inputs.worker_ids[1];
    inputs.actual_workers = two_active_workers ? 2U : 1U;
    inputs.route = two_active_workers
                       ? "openmp-sections-reference-query-v1"
                       : "openmp-sections-single-active-worker-v1";
  } else {
    if (inputs.requested_workers == 2U) {
      inputs.route = "serial-nested-openmp-guard-v1";
    }
    read_role(0);
    read_role(1);
  }
#else
  read_role(0);
  read_role(1);
#endif

  for (const auto& failure : inputs.failures) {
    if (failure != nullptr) {
      std::rethrow_exception(failure);
    }
  }
  return inputs;
}

void MaybeInjectPublicationRace(const std::filesystem::path& final_path) {
  const char* requested = std::getenv("RAMAG_TEST_RACE_TARGET");
  if (requested == nullptr || final_path.extension() !=
                                  ("." + std::string(requested))) {
    return;
  }
  std::error_code exists_error;
  if (std::filesystem::exists(final_path, exists_error) || exists_error) {
    return;
  }
  std::ofstream competitor(final_path, std::ios::binary | std::ios::out);
  if (!competitor) {
    throw WriterError("cannot create injected publication-race target");
  }
  competitor << "race-owner\n";
  competitor.close();
  if (!competitor) {
    throw WriterError("cannot close injected publication-race target");
  }
}

struct ArtifactReport {
  std::string format;
  std::filesystem::path path;
  std::uintmax_t bytes{};
};

ArtifactReport Artifact(std::string format,
                        const std::filesystem::path& temporary,
                        const std::filesystem::path& final_path) {
  std::error_code error;
  const auto bytes = std::filesystem::file_size(temporary, error);
  if (error) {
    throw WriterError("cannot inspect temporary artifact " + temporary.string() +
                      ": " + error.message());
  }
  return {std::move(format), final_path, bytes};
}

std::string SeedRoute(SeedMode mode) {
  switch (mode) {
    case SeedMode::FastHierarchical:
      return "baseline-formal-mam-skeleton+whole-query-mem-filter";
    case SeedMode::MumReference:
      return "baseline-formal-mam";
    case SeedMode::MaxMatch:
      return "baseline-formal-mem";
    case SeedMode::Mum:
      return "baseline-formal-mum";
    case SeedMode::Smem:
      return "baseline-generalized-smem";
  }
  return "invalid";
}

void EnsureTargetsAbsent(const RunSpec& spec, const OutputPaths& paths) {
  std::vector<std::filesystem::path> targets;
  if (spec.formats.sam) {
    targets.push_back(paths.sam);
  }
  if (spec.formats.paf) {
    targets.push_back(paths.paf);
  }
  if (spec.formats.delta) {
    targets.push_back(paths.delta);
  }
  if (spec.formats.maf) {
    targets.push_back(paths.maf);
  }
  if (spec.formats.chain) {
    targets.push_back(paths.chain);
  }
  for (const auto& target : targets) {
    std::error_code error;
    if (std::filesystem::exists(target, error)) {
      throw CliError("output target already exists (no --force in v1): " +
                     target.string());
    }
    if (error) {
      throw CliError("cannot inspect output target " + target.string() + ": " +
                     error.message());
    }
  }
}

void VerifyPublishedArtifact(const ArtifactReport& artifact) {
  std::error_code error;
  const auto bytes = std::filesystem::file_size(artifact.path, error);
  if (error || bytes != artifact.bytes) {
    throw WriterError("published artifact size mismatch: " + artifact.path.string());
  }
}

}  // namespace

RunOutcome RunAlignmentPipeline(const RunSpec& spec,
                                std::string invocation,
                                std::filesystem::path binary_path, RunLogger* logger) {
  ValidateRunSpec(spec);
  if (!spec.save_index_path.empty() && !SufkitAdapterAvailable()) {
    throw DependencyError("--save requires the pinned Sufkit production adapter");
  }
  const auto total_begin = Clock::now();
  const auto paths = MakeOutputPaths(spec);
  const auto run_id = logger ? logger->Id() : MakeNonce();
  const auto run_directory = spec.work_dir / "runs" / run_id;
  std::string stage = "configuration";
  std::filesystem::path retained_index;
  EnsureDirectory(spec.work_dir, "work directory");
  EnsureDirectory(run_directory, "run diagnostic directory");
  ProgressSession progress(spec.progress, run_id, spec.threads, logger);

  try {
    progress.Stage("configuration");
    EnsureDirectory(spec.output_prefix.parent_path(), "output directory");
    EnsureTargetsAbsent(spec, paths);
    std::map<std::string, double> timings;
    stage = "input";
    progress.Stage(stage);
    const auto input_begin = Clock::now();
    FastaData reference;
    FastaData query;
    std::string input_route;
    std::string input_parallel_route = "serial-reference-query-v1";
    std::uint32_t input_requested_workers = 1;
    std::uint32_t input_actual_workers = 1;
    std::map<std::string, std::string> adapter_provenance;
    const auto reference_compression =
        DetectFastaCompression(spec.reference_path);
    const auto query_compression = DetectFastaCompression(spec.query_path);
    adapter_provenance["input.reference.compression"] =
        reference_compression == FastaCompression::Gzip ? "gzip" : "plain";
    adapter_provenance["input.query.compression"] =
        query_compression == FastaCompression::Gzip ? "gzip" : "plain";
    auto inputs = ReadFastaInputs(spec);
    input_requested_workers = inputs.requested_workers;
    input_actual_workers = inputs.actual_workers;
    input_parallel_route = std::move(inputs.route);
    reference = std::move(inputs.results[0]);
    query = std::move(inputs.results[1]);
    input_route = "kseq-zlib-strict-fasta";
    if (logger) logger->Info("input reference_bases=" + std::to_string(reference.total_bases) +
        " query_bases=" + std::to_string(query.total_bases) + " actual_workers=" + std::to_string(input_actual_workers));
    timings["input"] = SecondsBetween(input_begin, Clock::now());
    progress.Update(static_cast<std::uint64_t>(reference.sequences.size() +
                                               query.sequences.size()),
                    static_cast<std::uint64_t>(reference.sequences.size() +
                                               query.sequences.size()),
                    "FASTA records");

    stage = "index-build";
    MaybeInjectPrecomputeFailure();
    const auto alignment_begin = Clock::now();
    AlignmentResult result;
    std::string index_route;
    std::string seeding_route;
    std::uint32_t actual_threads = 1;
    AlignmentOptions alignment_options = spec.alignment;
    alignment_options.worker_threads = spec.threads;
    alignment_options.resident_bytes_callback = CurrentResidentBytes;
    alignment_options.progress_callback =
        [&](std::string_view phase, std::uint64_t completed,
            std::uint64_t total) {
          stage.assign(phase);
          std::string detail;
          if (phase == "seed-merge") {
            detail = "selected_seeds=" + std::to_string(total);
          } else if (phase == "chaining") {
            detail = "merged_seeds=" + std::to_string(total);
          } else if (phase == "extension") {
            detail = "chains=" + std::to_string(total);
          } else if (phase == "conflict-resolution") {
            detail = "candidate_alignments=" + std::to_string(total);
          }
          progress.Stage(stage, completed, total, std::move(detail));
        };
    alignment_options.interruption_callback =
        [](std::string_view phase) { CheckInterruption(phase); };
    if (SufkitAdapterAvailable()) {
      SufkitSeedResult seed_result;
      RunStatistics seed_statistics;
      const auto memory_begin = Clock::now();
      const auto observe_index = [&](std::string phase, std::uint64_t index_bytes) {
        MemoryObservation observation;
        observation.stage = std::move(phase);
        observation.elapsed_seconds = SecondsBetween(memory_begin, Clock::now());
        observation.rss_bytes = CurrentResidentBytes();
        observation.index_estimated_bytes = index_bytes;
        observation.seed_capacity_bytes = seed_result.seeds.capacity() * sizeof(Seed);
        seed_statistics.memory_observations.push_back(std::move(observation));
      };
      // Seeds and their statistics own all data needed by the alignment core.
      // Release the full SA/ISA/LCP before allocating pairwise workspaces.
      {
      const auto index_begin = Clock::now();
      const bool load_persistent_index = !spec.reference_index_path.empty();
      stage = load_persistent_index ? "index-load" : "index-build";
      progress.Stage(stage, 0, 0,
                     load_persistent_index ? spec.reference_index_path.string()
                                           : "ephemeral");
      SufkitIndexOptions core_index_options{spec.threads};
      core_index_options.mam_worker_cap = spec.threads;
      core_index_options.mam_workspace_limit_bytes = UINT64_MAX;
      core_index_options.boundary_mem_occurrence_limit = UINT64_MAX;
      auto index = [&]() {
        if (load_persistent_index) {
          return SufkitSeedIndex::Load(
              spec.reference_index_path, reference.sequences,
              core_index_options);
        }
        const auto build_affinity = CurrentCpuAffinity();
        ValidateCpuAffinityBudget(build_affinity, spec.threads, "Sufkit alignment build");
        adapter_provenance["sufkit.index.pre_build_cpu_list"] = build_affinity.CpuList();
        adapter_provenance["sufkit.index.pre_build_cpu_count"] = std::to_string(build_affinity.logical_cpus.size());
        return SufkitSeedIndex::Build(
            reference.sequences, core_index_options);
      }();
      const double index_seconds = SecondsBetween(index_begin, Clock::now());
      timings[load_persistent_index ? "index_load" : "index_build"] =
          index_seconds;
      adapter_provenance["sufkit.index.save_path"] = spec.save_index_path.string();
      adapter_provenance["sufkit.index.persisted"] = "false";
      if (!spec.save_index_path.empty()) {
        stage = "index-save";
        IndexSpec save_spec{spec.reference_path, spec.save_index_path,
                            spec.work_dir, spec.threads, spec.progress};
        const auto saved = SaveReferenceIndex(save_spec, reference, index,
            invocation, binary_path, progress, index_seconds, timings, logger);
        retained_index = saved.index;
        adapter_provenance["sufkit.index.persisted"] = "true";
      }
      adapter_provenance["sufkit.index.action"] =
          load_persistent_index ? "loaded" : "built";
      adapter_provenance["sufkit.index.loaded_with_commit"] =
          std::string(kRequiredSufkitCommit);
      adapter_provenance["sufkit.index.path"] =
          load_persistent_index
              ? AbsolutePathString(spec.reference_index_path)
              : (spec.save_index_path.empty() ? "ephemeral" : AbsolutePathString(spec.save_index_path));
      adapter_provenance["sufkit.index.reference_validation"] =
          load_persistent_index ? "sufkit-crc+structure+metadata+fingerprint"
                                : "built-from-current-reference";
      stage = "seed-enumeration";
      progress.Stage(stage, 0,
                     static_cast<std::uint64_t>(query.sequences.size()));
      observe_index("index-ready", index.BuildStatistics().resident_core_bytes);
      seed_result = index.Enumerate(reference.sequences, query.sequences,
                                         alignment_options);
      seed_statistics.memory_observations.insert(seed_statistics.memory_observations.end(),
          std::make_move_iterator(seed_result.statistics.memory_observations.begin()),
          std::make_move_iterator(seed_result.statistics.memory_observations.end()));
      observe_index("seeds-enumerated", index.BuildStatistics().resident_core_bytes);
      seed_statistics.mem_seed_count =
          seed_result.statistics.mem_occurrence_count;
      seed_statistics.mam_seed_count =
          seed_result.statistics.mam_occurrence_count;
      seed_statistics.mum_seed_count =
          seed_result.statistics.mum_occurrence_count;
      seed_statistics.smem_interval_count =
          seed_result.statistics.smem_interval_count;
      seed_statistics.smem_coordinate_seed_count =
          seed_result.statistics.smem_occurrence_count;
      seed_statistics.selected_seed_count =
          seed_result.statistics.selected_seed_count;
      seed_statistics.seed_seconds = seed_result.statistics.enumeration_seconds;
      seed_statistics.actual_seed_route = seed_result.statistics.actual_route;
      seed_statistics.seed_requested_threads =
          seed_result.statistics.query_requested_threads;
      seed_statistics.seed_scheduled_threads =
          seed_result.statistics.query_scheduled_threads;
      seed_statistics.seed_worker_threads =
          seed_result.statistics.query_actual_threads;
      seed_statistics.seed_task_count =
          seed_result.statistics.query_strand_task_count;
      seed_statistics.seed_tasks_completed =
          seed_result.statistics.query_strand_tasks_completed;
      seed_statistics.seed_parallel_route =
          seed_result.statistics.query_parallel_route;
      seed_statistics.seed_mam_worker_cap =
          seed_result.statistics.mam_worker_cap;
      seed_statistics.seed_mam_tile_bases =
          seed_result.statistics.mam_tile_bases;
      seed_statistics.seed_oriented_query_count =
          seed_result.statistics.oriented_query_count;
      seed_statistics.seed_mam_tile_task_count =
          seed_result.statistics.mam_tile_task_count;
      seed_statistics.seed_mam_tile_tasks_completed =
          seed_result.statistics.mam_tile_tasks_completed;
      seed_statistics.seed_mam_short_query_task_count =
          seed_result.statistics.mam_short_query_task_count;
      seed_statistics.seed_mam_boundary_task_count =
          seed_result.statistics.mam_boundary_task_count;
      seed_statistics.seed_mam_boundary_tasks_completed =
          seed_result.statistics.mam_boundary_tasks_completed;
      seed_statistics.seed_mam_tile_raw_count =
          seed_result.statistics.mam_tile_raw_count;
      seed_statistics.seed_mam_tile_globally_maximal_count =
          seed_result.statistics.mam_tile_globally_maximal_count;
      seed_statistics.seed_mam_boundary_mem_raw_count =
          seed_result.statistics.mam_boundary_mem_raw_count;
      seed_statistics.seed_mam_boundary_pattern_count =
          seed_result.statistics.mam_boundary_pattern_count;
      seed_statistics.seed_mam_boundary_reference_unique_count =
          seed_result.statistics.mam_boundary_reference_unique_count;
      seed_statistics.seed_mam_boundary_recovered_count =
          seed_result.statistics.mam_boundary_recovered_count;
      seed_statistics.seed_mam_boundary_mem_occurrence_limit =
          seed_result.statistics.mam_boundary_mem_occurrence_limit;
      seed_statistics.seed_mam_workspace_baseline_bytes =
          seed_result.statistics.mam_workspace_baseline_bytes;
      seed_statistics.seed_mam_workspace_peak_bytes =
          seed_result.statistics.mam_workspace_peak_bytes;
      seed_statistics.seed_mam_workspace_limit_bytes =
          seed_result.statistics.mam_workspace_limit_bytes;
      timings["seed_enumeration"] =
          seed_result.statistics.enumeration_seconds;
      const auto& index_stats = seed_result.statistics.index;
      adapter_provenance["sufkit.library_version"] =
          index_stats.library_version;
      adapter_provenance["sufkit.commit"] =
          std::string{kRequiredSufkitCommit};
      adapter_provenance["sufkit.source_state"] =
          "exact-commit-clean-at-configure";
      adapter_provenance["sufkit.format_version"] =
          index_stats.format_version;
      adapter_provenance["sufkit.backend"] = index_stats.backend;
      adapter_provenance["sufkit.backend_signature"] =
          index_stats.backend_signature;
      adapter_provenance["sufkit.resource_profile"] =
          index_stats.resource_profile;
      adapter_provenance["sufkit.sampling_rate"] =
          std::to_string(index_stats.sampling_rate);
      adapter_provenance["sufkit.coordinate_width_bits"] =
          std::to_string(index_stats.coordinate_width);
      adapter_provenance["sufkit.stored_coordinate_width_bits"] =
          std::to_string(index_stats.stored_coordinate_width);
      adapter_provenance["sufkit.lcp_encoding"] =
          index_stats.lcp_encoding;
      adapter_provenance["sufkit.acceleration"] = index_stats.acceleration;
      adapter_provenance["sufkit.lookup_acceleration"] =
          index_stats.lookup_acceleration;
      adapter_provenance["sufkit.requested_threads"] =
          std::to_string(index_stats.requested_threads);
      adapter_provenance["sufkit.build.caps_min_reference_bases"] =
          std::to_string(index_stats.parallel_caps_min_reference_bases);
      adapter_provenance["sufkit.query.requested_threads"] =
          std::to_string(seed_result.statistics.query_requested_threads);
      adapter_provenance["sufkit.query.scheduled_threads"] =
          std::to_string(seed_result.statistics.query_scheduled_threads);
      adapter_provenance["sufkit.query.actual_threads"] =
          std::to_string(seed_result.statistics.query_actual_threads);
      adapter_provenance["sufkit.query.strand_tasks"] =
          std::to_string(seed_result.statistics.query_strand_task_count);
      adapter_provenance["sufkit.query.strand_tasks_completed"] =
          std::to_string(seed_result.statistics.query_strand_tasks_completed);
      adapter_provenance["sufkit.query.parallel_route"] =
          seed_result.statistics.query_parallel_route;
      adapter_provenance["sufkit.query.mam_worker_cap"] =
          std::to_string(seed_result.statistics.mam_worker_cap);
      adapter_provenance["sufkit.query.mam_tile_bases"] =
          std::to_string(seed_result.statistics.mam_tile_bases);
      adapter_provenance["sufkit.query.oriented_queries"] =
          std::to_string(seed_result.statistics.oriented_query_count);
      adapter_provenance["sufkit.query.mam_tile_tasks"] =
          std::to_string(seed_result.statistics.mam_tile_task_count);
      adapter_provenance["sufkit.query.mam_boundary_tasks"] =
          std::to_string(seed_result.statistics.mam_boundary_task_count);
      adapter_provenance["sufkit.query.mam_boundary_mem_raw"] =
          std::to_string(seed_result.statistics.mam_boundary_mem_raw_count);
      adapter_provenance["sufkit.query.mam_boundary_recovered"] =
          std::to_string(seed_result.statistics.mam_boundary_recovered_count);
      adapter_provenance["sufkit.query.mam_workspace_peak_bytes"] =
          std::to_string(seed_result.statistics.mam_workspace_peak_bytes);
      adapter_provenance["sufkit.reference_contigs"] =
          std::to_string(index_stats.reference_contigs);
      adapter_provenance["sufkit.reference_bases"] =
          std::to_string(index_stats.reference_bases);
      adapter_provenance["sufkit.reference_ambiguous_bases"] =
          std::to_string(index_stats.reference_ambiguous_bases);
      adapter_provenance["sufkit.text_symbols"] =
          std::to_string(index_stats.text_symbols);
      adapter_provenance["sufkit.suffix_rows"] =
          std::to_string(index_stats.suffix_count);
      adapter_provenance["sufkit.serialized_bytes"] =
          std::to_string(index_stats.serialized_bytes);
      adapter_provenance["sufkit.text_bytes"] =
          std::to_string(index_stats.text_bytes);
      adapter_provenance["sufkit.suffix_array_bytes"] =
          std::to_string(index_stats.suffix_array_bytes);
      adapter_provenance["sufkit.inverse_suffix_array_bytes"] =
          std::to_string(index_stats.inverse_suffix_array_bytes);
      adapter_provenance["sufkit.lcp_bytes"] =
          std::to_string(index_stats.lcp_bytes);
      adapter_provenance["sufkit.lcp_primary_bytes"] =
          std::to_string(index_stats.lcp_primary_bytes);
      adapter_provenance["sufkit.lcp_overflow_anchors"] =
          std::to_string(index_stats.lcp_overflow_anchors);
      adapter_provenance["sufkit.lcp_overflow_bytes"] =
          std::to_string(index_stats.lcp_overflow_bytes);
      adapter_provenance["sufkit.lcp_guide_bytes"] =
          std::to_string(index_stats.lcp_guide_bytes);
      adapter_provenance["sufkit.child_bytes"] =
          std::to_string(index_stats.child_bytes);
      adapter_provenance["sufkit.auxiliary_bytes"] =
          std::to_string(index_stats.auxiliary_bytes);
      adapter_provenance["sufkit.resident_core_bytes"] =
          std::to_string(index_stats.resident_core_bytes);
      adapter_provenance["sufkit.learned_index_bytes"] =
          std::to_string(index_stats.learned_index_bytes);
      adapter_provenance["sufkit.prefix_directory.enabled"] =
          index_stats.prefix_directory_enabled ? "true" : "false";
      adapter_provenance["sufkit.prefix_directory.bytes"] =
          std::to_string(index_stats.prefix_directory_bytes);
      adapter_provenance["sufkit.build.total_seconds"] =
          std::to_string(index_stats.total_build_seconds);
      adapter_provenance["sufkit.build.suffix_array_seconds"] =
          std::to_string(index_stats.suffix_array_seconds);
      adapter_provenance["sufkit.build.storage_compaction_seconds"] =
          std::to_string(index_stats.storage_compaction_seconds);
      adapter_provenance["sufkit.build.isa_seconds"] =
          std::to_string(index_stats.isa_seconds);
      adapter_provenance["sufkit.build.lcp_seconds"] =
          std::to_string(index_stats.lcp_seconds);
      adapter_provenance["sufkit.build.child_seconds"] =
          std::to_string(index_stats.child_seconds);
      adapter_provenance["sufkit.build.learned_index_seconds"] =
          std::to_string(index_stats.learned_index_seconds);
      index_route = "sufkit-full-sa:" + index_stats.backend +
                    ":sampling=" + std::to_string(index_stats.sampling_rate) +
                    ":acceleration=" + index_stats.acceleration;
      seeding_route = seed_result.statistics.actual_route;
      }
      observe_index("index-released", 0);
      result = AlignPairwiseFromSeeds(
          reference.sequences, query.sequences, alignment_options,
          std::move(seed_result.seeds), std::move(seed_statistics));
    } else {
      stage = "seed-enumeration";
      progress.Stage(stage);
      RunStatistics oracle_statistics;
      auto oracle_seeds = EnumerateSeeds(reference.sequences, query.sequences,
                                        alignment_options, &oracle_statistics);
      result = AlignPairwiseFromSeeds(reference.sequences, query.sequences,
                                   alignment_options, std::move(oracle_seeds),
                                   std::move(oracle_statistics));
      index_route = "bundled-canonical-kmer-oracle-index";
      seeding_route = SeedRoute(spec.alignment.seed_mode);
    }
    actual_threads = std::max(
        {input_actual_workers, result.statistics.seed_worker_threads,
         result.statistics.chaining_worker_threads,
         result.statistics.extension_worker_threads});
    const double seed_chain_extend_seconds =
        SecondsBetween(alignment_begin, Clock::now());
    timings["seed_chain_extend"] = std::max(
        0.0, seed_chain_extend_seconds -
                 result.statistics.conflict_resolution_seconds);

    // Alignment-core conflict resolution is the sole authority for
    // de-duplication, primary assignment, and deterministic record order.
    // Writers consume this immutable result and never select independently.
    timings["conflict_resolution"] =
        result.statistics.conflict_resolution_seconds;

    const auto temporary_sam = TemporarySibling(paths.sam, run_id);
    const auto temporary_paf = TemporarySibling(paths.paf, run_id);
    const auto temporary_delta = TemporarySibling(paths.delta, run_id);
    const auto temporary_maf = TemporarySibling(paths.maf, run_id);
    const auto temporary_chain = TemporarySibling(paths.chain, run_id);
    TemporaryFiles temporary_files;
    PublishedFiles published_files;
    std::vector<std::pair<std::filesystem::path, std::filesystem::path>> publish;
    std::vector<ArtifactReport> artifacts;

    stage = "writers";
    progress.Stage(stage, 0,
                   static_cast<std::uint64_t>(result.alignments.size()),
                   "alignment_records=" +
                       std::to_string(result.alignments.size()));
    const auto writers_begin = Clock::now();
    if (spec.formats.sam) {
      stage = "writer-sam";
      progress.Stage(stage, 0,
                     static_cast<std::uint64_t>(result.alignments.size()),
                     "output_records=" +
                         std::to_string(result.alignments.size()));
      temporary_files.Add(temporary_sam);
      const auto begin = Clock::now();
      WriteFile(temporary_sam, [&](std::ostream& output) {
        WriteSam(output, reference, query, result.alignments, invocation);
      });
      MaybeInjectWriterFailure("sam");
      stage = "validation-sam";
      progress.Stage(stage, 0, 1);
      MaybeInjectValidatorFailure("sam");
      ValidateSamFile(temporary_sam);
      timings["writer_sam"] = SecondsBetween(begin, Clock::now());
      artifacts.push_back(Artifact("sam", temporary_sam, paths.sam));
      publish.emplace_back(temporary_sam, paths.sam);
    }
    if (spec.formats.paf) {
      stage = "writer-paf";
      progress.Stage(stage, 0,
                     static_cast<std::uint64_t>(result.alignments.size()),
                     "output_records=" +
                         std::to_string(result.alignments.size()));
      temporary_files.Add(temporary_paf);
      const auto begin = Clock::now();
      WriteFile(temporary_paf, [&](std::ostream& output) {
        WritePaf(output, reference, query, result.alignments);
      });
      MaybeInjectWriterFailure("paf");
      stage = "validation-paf";
      progress.Stage(stage, 0, 1);
      MaybeInjectValidatorFailure("paf");
      ValidatePafFile(temporary_paf);
      timings["writer_paf"] = SecondsBetween(begin, Clock::now());
      artifacts.push_back(Artifact("paf", temporary_paf, paths.paf));
      publish.emplace_back(temporary_paf, paths.paf);
    }
    if (spec.formats.delta) {
      stage = "writer-delta";
      progress.Stage(stage, 0,
                     static_cast<std::uint64_t>(result.alignments.size()),
                     "output_records=" +
                         std::to_string(result.alignments.size()));
      temporary_files.Add(temporary_delta);
      const auto begin = Clock::now();
      WriteFile(temporary_delta, [&](std::ostream& output) {
        WriteDelta(output, reference, query, result.alignments);
      });
      MaybeInjectWriterFailure("delta");
      stage = "validation-delta";
      progress.Stage(stage, 0, 1);
      MaybeInjectValidatorFailure("delta");
      ValidateDeltaFile(temporary_delta);
      timings["writer_delta"] = SecondsBetween(begin, Clock::now());
      artifacts.push_back(Artifact("delta", temporary_delta, paths.delta));
      publish.emplace_back(temporary_delta, paths.delta);
    }
    if (spec.formats.maf) {
      stage = "writer-maf";
      progress.Stage(stage, 0,
                     static_cast<std::uint64_t>(result.alignments.size()),
                     "output_records=" +
                         std::to_string(result.alignments.size()));
      temporary_files.Add(temporary_maf);
      const auto begin = Clock::now();
      WriteFile(temporary_maf, [&](std::ostream& output) {
        WriteMaf(output, reference, query, result.alignments);
      });
      MaybeInjectWriterFailure("maf");
      stage = "validation-maf";
      progress.Stage(stage, 0, 1);
      MaybeInjectValidatorFailure("maf");
      ValidateMafFile(temporary_maf);
      timings["writer_maf"] = SecondsBetween(begin, Clock::now());
      artifacts.push_back(Artifact("maf", temporary_maf, paths.maf));
      publish.emplace_back(temporary_maf, paths.maf);
    }
    if (spec.formats.chain) {
      stage = "writer-chain";
      progress.Stage(stage, 0,
                     static_cast<std::uint64_t>(result.alignments.size()),
                     "output_records=" +
                         std::to_string(result.alignments.size()));
      temporary_files.Add(temporary_chain);
      const auto begin = Clock::now();
      WriteFile(temporary_chain, [&](std::ostream& output) {
        WriteChain(output, reference, query, result.alignments);
      });
      MaybeInjectWriterFailure("chain");
      stage = "validation-chain";
      progress.Stage(stage, 0, 1);
      MaybeInjectValidatorFailure("chain");
      ValidateChainFile(temporary_chain);
      timings["writer_chain"] = SecondsBetween(begin, Clock::now());
      artifacts.push_back(Artifact("chain", temporary_chain, paths.chain));
      publish.emplace_back(temporary_chain, paths.chain);
    }
    timings["writers_total"] = SecondsBetween(writers_begin, Clock::now());

    stage = "publication";
    progress.Stage(stage, 0, static_cast<std::uint64_t>(publish.size()));
    std::uint64_t published_count = 0;
    for (const auto& [temporary, final_path] : publish) {
      CheckInterruption(stage);
      published_files.Add(temporary, final_path);
      MaybeInjectPublicationRace(final_path);
      PublishNewFile(temporary, final_path);
      progress.Update(++published_count,
                      static_cast<std::uint64_t>(publish.size()));
    }
    for (const auto& artifact : artifacts) {
      VerifyPublishedArtifact(artifact);
      if (artifact.format == "sam") {
        ValidateSamFile(artifact.path);
      } else if (artifact.format == "paf") {
        ValidatePafFile(artifact.path);
      } else if (artifact.format == "delta") {
        ValidateDeltaFile(artifact.path);
      } else if (artifact.format == "maf") {
        ValidateMafFile(artifact.path);
      } else if (artifact.format == "chain") {
        ValidateChainFile(artifact.path);
      }
    }
    progress.Finish(std::to_string(result.statistics.alignment_count) +
                    " alignments");
    if (logger) {
      logger->Info("statistics.reference_contigs=" + std::to_string(result.statistics.reference_contigs), false);
      logger->Info("statistics.query_contigs=" + std::to_string(result.statistics.query_contigs), false);
      logger->Info("statistics.reference_bases=" + std::to_string(result.statistics.reference_bases), false);
      logger->Info("statistics.query_bases=" + std::to_string(result.statistics.query_bases), false);
      logger->Info("statistics.reference_ambiguous_bases=" + std::to_string(result.statistics.reference_ambiguous_bases), false);
      logger->Info("statistics.query_ambiguous_bases=" + std::to_string(result.statistics.query_ambiguous_bases), false);
      logger->Info("statistics.mem_seed_count=" + std::to_string(result.statistics.mem_seed_count), false);
      logger->Info("statistics.mam_seed_count=" + std::to_string(result.statistics.mam_seed_count), false);
      logger->Info("statistics.mum_seed_count=" + std::to_string(result.statistics.mum_seed_count), false);
      logger->Info("statistics.smem_interval_count=" + std::to_string(result.statistics.smem_interval_count), false);
      logger->Info("statistics.smem_coordinate_seed_count=" + std::to_string(result.statistics.smem_coordinate_seed_count), false);
      logger->Info("statistics.selected_seed_count=" + std::to_string(result.statistics.selected_seed_count), false);
      logger->Info("statistics.merged_seed_count=" + std::to_string(result.statistics.merged_seed_count), false);
      logger->Info("statistics.chain_count=" + std::to_string(result.statistics.chain_count), false);
      logger->Info("statistics.candidate_alignment_count=" + std::to_string(result.statistics.candidate_alignment_count), false);
      logger->Info("statistics.conflict_rejected_alignment_count=" + std::to_string(result.statistics.conflict_rejected_alignment_count), false);
      logger->Info("statistics.alignment_count=" + std::to_string(result.statistics.alignment_count), false);
      logger->Info("statistics.exact_gap_count=" + std::to_string(result.statistics.exact_gap_count), false);
      logger->Info("statistics.ungapped_gap_count=" + std::to_string(result.statistics.ungapped_gap_count), false);
      logger->Info("statistics.dp_gap_count=" + std::to_string(result.statistics.dp_gap_count), false);
      logger->Info("statistics.seed_requested_threads=" + std::to_string(result.statistics.seed_requested_threads), false);
      logger->Info("statistics.seed_scheduled_threads=" + std::to_string(result.statistics.seed_scheduled_threads), false);
      logger->Info("statistics.seed_worker_threads=" + std::to_string(result.statistics.seed_worker_threads), false);
      logger->Info("statistics.seed_task_count=" + std::to_string(result.statistics.seed_task_count), false);
      logger->Info("statistics.seed_tasks_completed=" + std::to_string(result.statistics.seed_tasks_completed), false);
      logger->Info("statistics.seed_parallel_route=" + result.statistics.seed_parallel_route, false);
      logger->Info("statistics.seed_mam_worker_cap=" + std::to_string(result.statistics.seed_mam_worker_cap), false);
      logger->Info("statistics.seed_mam_tile_bases=" + std::to_string(result.statistics.seed_mam_tile_bases), false);
      logger->Info("statistics.seed_oriented_query_count=" + std::to_string(result.statistics.seed_oriented_query_count), false);
      logger->Info("statistics.seed_mam_tile_task_count=" + std::to_string(result.statistics.seed_mam_tile_task_count), false);
      logger->Info("statistics.seed_mam_tile_tasks_completed=" + std::to_string(result.statistics.seed_mam_tile_tasks_completed), false);
      logger->Info("statistics.seed_mam_short_query_task_count=" + std::to_string(result.statistics.seed_mam_short_query_task_count), false);
      logger->Info("statistics.seed_mam_boundary_task_count=" + std::to_string(result.statistics.seed_mam_boundary_task_count), false);
      logger->Info("statistics.seed_mam_boundary_tasks_completed=" + std::to_string(result.statistics.seed_mam_boundary_tasks_completed), false);
      logger->Info("statistics.seed_mam_tile_raw_count=" + std::to_string(result.statistics.seed_mam_tile_raw_count), false);
      logger->Info("statistics.seed_mam_tile_globally_maximal_count=" + std::to_string(result.statistics.seed_mam_tile_globally_maximal_count), false);
      logger->Info("statistics.seed_mam_boundary_mem_raw_count=" + std::to_string(result.statistics.seed_mam_boundary_mem_raw_count), false);
      logger->Info("statistics.seed_mam_boundary_pattern_count=" + std::to_string(result.statistics.seed_mam_boundary_pattern_count), false);
      logger->Info("statistics.seed_mam_boundary_reference_unique_count=" + std::to_string(result.statistics.seed_mam_boundary_reference_unique_count), false);
      logger->Info("statistics.seed_mam_boundary_recovered_count=" + std::to_string(result.statistics.seed_mam_boundary_recovered_count), false);
      logger->Info("statistics.seed_mam_boundary_mem_occurrence_limit=" + std::to_string(result.statistics.seed_mam_boundary_mem_occurrence_limit), false);
      logger->Info("statistics.seed_mam_workspace_baseline_bytes=" + std::to_string(result.statistics.seed_mam_workspace_baseline_bytes), false);
      logger->Info("statistics.seed_mam_workspace_peak_bytes=" + std::to_string(result.statistics.seed_mam_workspace_peak_bytes), false);
      logger->Info("statistics.seed_mam_workspace_limit_bytes=" + std::to_string(result.statistics.seed_mam_workspace_limit_bytes), false);
      logger->Info("statistics.chaining_route=" + result.statistics.chaining_route, false);
      logger->Info("statistics.chaining_candidate_pairs=" + std::to_string(result.statistics.chaining_candidate_pairs), false);
      logger->Info("statistics.chaining_legal_edges=" + std::to_string(result.statistics.chaining_legal_edges), false);
      logger->Info("statistics.chaining_components=" + std::to_string(result.statistics.chaining_components), false);
      logger->Info("statistics.chaining_max_bucket_occupancy=" + std::to_string(result.statistics.chaining_max_bucket_occupancy), false);
      logger->Info("statistics.chaining_max_component_seeds=" + std::to_string(result.statistics.chaining_max_component_seeds), false);
      logger->Info("statistics.chaining_max_component_edges=" + std::to_string(result.statistics.chaining_max_component_edges), false);
      logger->Info("statistics.chaining_dp_passes=" + std::to_string(result.statistics.chaining_dp_passes), false);
      logger->Info("statistics.chaining_edge_relaxations=" + std::to_string(result.statistics.chaining_edge_relaxations), false);
      logger->Info("statistics.chaining_working_set_peak_bytes=" + std::to_string(result.statistics.chaining_working_set_peak_bytes), false);
      logger->Info("statistics.chaining_working_set_limit_bytes=" + std::to_string(result.statistics.chaining_working_set_limit_bytes), false);
      logger->Info("statistics.chaining_candidate_pair_limit=" + std::to_string(result.statistics.chaining_candidate_pair_limit), false);
      logger->Info("statistics.chaining_legal_edge_limit=" + std::to_string(result.statistics.chaining_legal_edge_limit), false);
      logger->Info("statistics.chaining_edge_relaxation_limit=" + std::to_string(result.statistics.chaining_edge_relaxation_limit), false);
      logger->Info("statistics.chaining_requested_threads=" + std::to_string(result.statistics.chaining_requested_threads), false);
      logger->Info("statistics.chaining_worker_threads=" + std::to_string(result.statistics.chaining_worker_threads), false);
      logger->Info("statistics.extension_worker_threads=" + std::to_string(result.statistics.extension_worker_threads), false);
      logger->Info("statistics.seed_seconds=" + std::to_string(result.statistics.seed_seconds), false);
      logger->Info("statistics.merge_seconds=" + std::to_string(result.statistics.merge_seconds), false);
      logger->Info("statistics.chain_and_extension_seconds=" + std::to_string(result.statistics.chain_and_extension_seconds), false);
      logger->Info("statistics.conflict_resolution_seconds=" + std::to_string(result.statistics.conflict_resolution_seconds), false);
      logger->Info("statistics.total_seconds=" + std::to_string(result.statistics.total_seconds), false);
      logger->Info("statistics.actual_seed_route=" + result.statistics.actual_seed_route, false);
      logger->Info("pairwise.global_ksw_calls=" + std::to_string(result.statistics.pairwise.global_ksw_calls), false);
      logger->Info("pairwise.endpoint_ksw_calls=" + std::to_string(result.statistics.pairwise.endpoint_ksw_calls), false);
      logger->Info("pairwise.global_ksw_seconds=" + std::to_string(result.statistics.pairwise.global_ksw_seconds), false);
      logger->Info("pairwise.endpoint_ksw_seconds=" + std::to_string(result.statistics.pairwise.endpoint_ksw_seconds), false);
      logger->Info("pairwise.link_candidate_checks=" + std::to_string(result.statistics.pairwise.link_candidate_checks), false);
      logger->Info("pairwise.link_direct_attempts=" + std::to_string(result.statistics.pairwise.link_direct_attempts), false);
      logger->Info("pairwise.link_fallback_attempts=" + std::to_string(result.statistics.pairwise.link_fallback_attempts), false);
      logger->Info("pairwise.link_long_gap_rejections=" + std::to_string(result.statistics.pairwise.link_long_gap_rejections), false);
      logger->Info("pairwise.link_closure_failures=" + std::to_string(result.statistics.pairwise.link_closure_failures), false);
      logger->Info("pairwise.seed_grouping_seconds=" + std::to_string(result.statistics.pairwise.seed_grouping_seconds), false);
      logger->Info("pairwise.output_conversion_seconds=" + std::to_string(result.statistics.pairwise.output_conversion_seconds), false);
      logger->Info("pairwise.recovery_enabled=" + std::to_string(result.statistics.pairwise.recovery_enabled), false);
      logger->Info("pairwise.recovery_candidates_checked=" + std::to_string(result.statistics.pairwise.recovery_candidates_checked), false);
      logger->Info("pairwise.recovery_proposed_fragments=" + std::to_string(result.statistics.pairwise.recovery_proposed_fragments), false);
      logger->Info("pairwise.recovery_requeued_fragments=" + std::to_string(result.statistics.pairwise.recovery_requeued_fragments), false);
      logger->Info("pairwise.recovery_accepted_fragments=" + std::to_string(result.statistics.pairwise.recovery_accepted_fragments), false);
      logger->Info("pairwise.recovery_reference_bases=" + std::to_string(result.statistics.pairwise.recovery_reference_bases), false);
      logger->Info("pairwise.recovery_query_bases=" + std::to_string(result.statistics.pairwise.recovery_query_bases), false);
      logger->Info("pairwise.recovery_seconds=" + std::to_string(result.statistics.pairwise.recovery_seconds), false);
      logger->Info("pairwise.gap_fill_strategy=" + result.statistics.pairwise.gap_fill_strategy, false);
      logger->Info("pairwise.gap_fill_adjacent_pairs=" + std::to_string(result.statistics.pairwise.gap_fill_adjacent_pairs), false);
      logger->Info("pairwise.gap_fill_geometry_rejected=" + std::to_string(result.statistics.pairwise.gap_fill_geometry_rejected), false);
      logger->Info("pairwise.gap_fill_occupied_rejected=" + std::to_string(result.statistics.pairwise.gap_fill_occupied_rejected), false);
      logger->Info("pairwise.gap_fill_flank_rejected=" + std::to_string(result.statistics.pairwise.gap_fill_flank_rejected), false);
      logger->Info("pairwise.gap_fill_n_rejected=" + std::to_string(result.statistics.pairwise.gap_fill_n_rejected), false);
      logger->Info("pairwise.gap_fill_candidates=" + std::to_string(result.statistics.pairwise.gap_fill_candidates), false);
      logger->Info("pairwise.gap_fill_exact=" + std::to_string(result.statistics.pairwise.gap_fill_exact), false);
      logger->Info("pairwise.gap_fill_ksw_calls=" + std::to_string(result.statistics.pairwise.gap_fill_ksw_calls), false);
      logger->Info("pairwise.gap_fill_quality_rejected=" + std::to_string(result.statistics.pairwise.gap_fill_quality_rejected), false);
      logger->Info("pairwise.gap_fill_nonexact_rejected=" + std::to_string(result.statistics.pairwise.gap_fill_nonexact_rejected), false);
      logger->Info("pairwise.gap_fill_conflict_rejected=" + std::to_string(result.statistics.pairwise.gap_fill_conflict_rejected), false);
      logger->Info("pairwise.gap_fill_accepted=" + std::to_string(result.statistics.pairwise.gap_fill_accepted), false);
      logger->Info("pairwise.gap_fill_reference_bases=" + std::to_string(result.statistics.pairwise.gap_fill_reference_bases), false);
      logger->Info("pairwise.gap_fill_query_bases=" + std::to_string(result.statistics.pairwise.gap_fill_query_bases), false);
      logger->Info("pairwise.gap_fill_paired_columns=" + std::to_string(result.statistics.pairwise.gap_fill_paired_columns), false);
      logger->Info("pairwise.gap_fill_ksw_seconds=" + std::to_string(result.statistics.pairwise.gap_fill_ksw_seconds), false);
      logger->Info("pairwise.gap_fill_seconds=" + std::to_string(result.statistics.pairwise.gap_fill_seconds), false);
      for (const auto& sample : result.statistics.memory_observations) {
        logger->Info("memory stage=" + sample.stage + " rss_bytes=" + std::to_string(sample.rss_bytes) + " index_bytes=" + std::to_string(sample.index_estimated_bytes) + " seed_capacity_bytes=" + std::to_string(sample.seed_capacity_bytes) + " cluster_capacity_bytes=" + std::to_string(sample.cluster_capacity_bytes) + " anchor_capacity_bytes=" + std::to_string(sample.anchor_capacity_bytes) + " cigar_capacity_bytes=" + std::to_string(sample.cigar_capacity_bytes), false);
      }
      for (const auto& [name,value] : adapter_provenance) logger->Info("index/input " + name + "=" + value, false);
      for (const auto& [name,seconds] : timings) logger->Info("timing stage=" + name + " seconds=" + std::to_string(seconds), false);
      for (const auto& artifact : artifacts) logger->Info("artifact format=" + artifact.format + " path=" + artifact.path.string() + " bytes=" + std::to_string(artifact.bytes), false);
      logger->Info("routes input=" + input_route + " parallel=" + input_parallel_route + " index=" + index_route + " seed=" + seeding_route + " requested_input_workers=" + std::to_string(input_requested_workers) + " actual_threads=" + std::to_string(actual_threads), false);
      CheckInterruption("publication");
      logger->Info("status=success exit_code=0 alignment_records=" + std::to_string(result.statistics.alignment_count) + " total_seconds=" + std::to_string(SecondsBetween(total_begin, Clock::now())));
      logger->Flush();
    }
    published_files.Commit();

    return {paths, result.statistics};
  } catch (const std::exception& error) {
    try { if (logger) logger->Error("status=" + std::string(dynamic_cast<const InterruptedError*>(&error) ? "interrupted" : "failed") + " exit_code=" + std::to_string(FailureExitCode(error)) + " stage=" + stage + " retained_reference_index=" + retained_index.string() + " message=" + error.what()); } catch (...) {}
    throw;
  }
}

}  // namespace ramag

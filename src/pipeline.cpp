#include "ramag/pipeline.hpp"
#if RAMAG_USE_PAIRWISE_CORE
#include "ramag/pairwise_core.hpp"
#endif

#include "ramag/alignment.hpp"
#include "ramag/fasta.hpp"
#include "ramag/manifest.hpp"
#include "ramag/reference_index.hpp"
#include "ramag/runtime.hpp"
#include "ramag/sha256.hpp"
#include "ramag/seqpro_adapter.hpp"
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

std::string CurrentUtcTimestamp() {
  const auto now = std::chrono::system_clock::now();
  const auto time = std::chrono::system_clock::to_time_t(now);
  std::tm utc{};
#if defined(_WIN32)
  gmtime_s(&utc, &time);
#else
  gmtime_r(&time, &utc);
#endif
  char buffer[32]{};
  if (std::strftime(buffer, sizeof(buffer), "%Y-%m-%dT%H:%M:%SZ", &utc) == 0) {
    return "unknown";
  }
  return buffer;
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

struct ParallelSeqProInputs {
  std::array<SeqProInputResult, 2> results;
  std::array<std::exception_ptr, 2> failures;
  std::array<int, 2> worker_ids{-1, -1};
  std::uint32_t requested_workers{1};
  std::uint32_t actual_workers{1};
  std::string route{"serial-reference-query-v1"};
};

ParallelSeqProInputs ReadSeqProInputs(const RunSpec& spec,
                                      const std::filesystem::path& work_dir) {
  ParallelSeqProInputs inputs;
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
      inputs.results[role_index] = ReadFastaWithSeqPro(
          is_reference ? spec.reference_path : spec.query_path,
          SeqProInputOptions{work_dir, std::string(role), 0, false});
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

ArtifactReport Artifact(std::string format,
                        const std::filesystem::path& temporary,
                        const std::filesystem::path& final_path) {
  std::error_code error;
  const auto bytes = std::filesystem::file_size(temporary, error);
  if (error) {
    throw WriterError("cannot inspect temporary artifact " + temporary.string() +
                      ": " + error.message());
  }
  return {std::move(format), final_path, bytes, "validated"};
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
  std::vector<std::filesystem::path> targets{paths.manifest, paths.complete};
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

void ValidateManifestFile(const std::filesystem::path& path,
                          std::string_view run_id,
                          std::string_view expected_content) {
  std::ifstream input(path, std::ios::binary);
  if (!input) {
    throw WriterError("cannot validate manifest: " + path.string());
  }
  const std::string content{std::istreambuf_iterator<char>(input),
                            std::istreambuf_iterator<char>()};
  if (input.bad() || content != expected_content ||
      content.find("\"schema_version\": 1") == std::string::npos ||
      content.find("\"status\": \"success\"") == std::string::npos ||
      content.find(std::string(run_id)) == std::string::npos) {
    throw WriterError("manifest self-validation failed: " + path.string());
  }
}

int FailureExitCode(const std::exception& error) noexcept {
  if (dynamic_cast<const UnsupportedSeedMode*>(&error) != nullptr ||
      dynamic_cast<const UnsupportedFastaFormat*>(&error) != nullptr) {
    return 4;
  }
  if (const auto* interrupted = dynamic_cast<const InterruptedError*>(&error)) {
    return interrupted->ExitCode();
  }
  if (dynamic_cast<const CliError*>(&error) != nullptr) {
    return 2;
  }
  if (dynamic_cast<const FastaError*>(&error) != nullptr) {
    return 3;
  }
  if (dynamic_cast<const DependencyError*>(&error) != nullptr) {
    return 5;
  }
  if (dynamic_cast<const WriterError*>(&error) != nullptr) {
    return 7;
  }
  if (dynamic_cast<const AlignmentError*>(&error) != nullptr) {
    return 6;
  }
  if (dynamic_cast<const std::bad_alloc*>(&error) != nullptr) {
    return 6;
  }
  return 8;
}

std::string EscapeFailureJson(std::string_view value) {
  std::string escaped;
  for (const char character : value) {
    if (character == '\\' || character == '\"') {
      escaped.push_back('\\');
      escaped.push_back(character);
    } else if (character == '\n' || character == '\r' || character == '\t') {
      escaped.push_back(' ');
    } else if (static_cast<unsigned char>(character) >= 0x20U) {
      escaped.push_back(character);
    }
  }
  return escaped;
}

void WriteFailureDiagnostic(const std::filesystem::path& run_directory,
                            std::string_view run_id,
                            std::string_view stage,
                            const std::exception& error) noexcept {
  try {
    std::error_code directory_error;
    std::filesystem::create_directories(run_directory, directory_error);
    if (directory_error) {
      return;
    }
    std::ofstream output(run_directory / "failure.json",
                         std::ios::binary | std::ios::trunc);
    if (!output) {
      return;
    }
    const auto* interrupted = dynamic_cast<const InterruptedError*>(&error);
    output << "{\n"
           << "  \"schema_version\": 1,\n"
           << "  \"status\": \""
           << (interrupted == nullptr ? "failed" : "interrupted") << "\",\n"
           << "  \"run_id\": \"" << EscapeFailureJson(run_id) << "\",\n"
           << "  \"observed_utc\": \""
           << EscapeFailureJson(CurrentUtcTimestamp()) << "\",\n"
           << "  \"exit_code\": " << FailureExitCode(error) << ",\n"
           << "  \"failed_stage\": \"" << EscapeFailureJson(stage) << "\",\n";
    if (interrupted != nullptr) {
      output << "  \"signal\": " << interrupted->SignalNumber() << ",\n";
    }
    output
           << "  \"message\": \"" << EscapeFailureJson(error.what()) << "\"\n"
           << "}\n";
  } catch (...) {
  }
}

}  // namespace

RunOutcome RunAlignmentPipeline(const RunSpec& spec,
                                std::string invocation,
                                std::filesystem::path binary_path) {
  ValidateRunSpec(spec);
  const auto total_begin = Clock::now();
  const auto paths = MakeOutputPaths(spec);
  const auto run_id = MakeNonce();
  const auto started_utc = CurrentUtcTimestamp();
  const auto run_directory = spec.work_dir / "runs" / run_id;
  std::string stage = "configuration";
  EnsureDirectory(spec.work_dir, "work directory");
  EnsureDirectory(run_directory, "run diagnostic directory");
  ProgressSession progress(spec.progress, run_id, spec.threads);

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
    if (SeqProAdapterAvailable() &&
        reference_compression == FastaCompression::Plain &&
        query_compression == FastaCompression::Plain) {
      const auto seqpro_work = run_directory / "seqpro";
      EnsureDirectory(seqpro_work, "SeqPro input work directory");
      auto inputs = ReadSeqProInputs(spec, seqpro_work);
      input_requested_workers = inputs.requested_workers;
      input_actual_workers = inputs.actual_workers;
      input_parallel_route = std::move(inputs.route);
      auto reference_input = std::move(inputs.results[0]);
      auto query_input = std::move(inputs.results[1]);
      adapter_provenance["seqpro.reference.fai"] =
          AbsolutePathString(reference_input.fasta_index_path);
      adapter_provenance["seqpro.reference.metadata"] =
          AbsolutePathString(reference_input.metadata_path);
      adapter_provenance["seqpro.reference.source_fai"] =
          AbsolutePathString(reference_input.source_fasta_index_path);
      adapter_provenance["seqpro.reference.source_fai_bytes"] =
          std::to_string(reference_input.source_fasta_index_size_bytes);
      adapter_provenance["seqpro.reference.source_fai_sha256"] =
          reference_input.source_fasta_index_sha256;
      adapter_provenance["seqpro.reference.source_fai.copy_status"] =
          reference_input.source_fasta_index_copied ? "copied" : "not-present";
      adapter_provenance["seqpro.reference.external_fai.adoption_status"] =
          reference_input.external_fasta_index_adopted ? "adopted"
                                                       : "not-applicable";
      adapter_provenance["seqpro.reference.build_action"] =
          reference_input.build_action;
      adapter_provenance["seqpro.reference.index_origin"] =
          reference_input.index_origin;
      adapter_provenance["seqpro.reference.verification"] =
          reference_input.verification_status;
      adapter_provenance["seqpro.query.fai"] =
          AbsolutePathString(query_input.fasta_index_path);
      adapter_provenance["seqpro.query.metadata"] =
          AbsolutePathString(query_input.metadata_path);
      adapter_provenance["seqpro.query.source_fai"] =
          AbsolutePathString(query_input.source_fasta_index_path);
      adapter_provenance["seqpro.query.source_fai_bytes"] =
          std::to_string(query_input.source_fasta_index_size_bytes);
      adapter_provenance["seqpro.query.source_fai_sha256"] =
          query_input.source_fasta_index_sha256;
      adapter_provenance["seqpro.query.source_fai.copy_status"] =
          query_input.source_fasta_index_copied ? "copied" : "not-present";
      adapter_provenance["seqpro.query.external_fai.adoption_status"] =
          query_input.external_fasta_index_adopted ? "adopted"
                                                   : "not-applicable";
      adapter_provenance["seqpro.query.build_action"] =
          query_input.build_action;
      adapter_provenance["seqpro.query.index_origin"] =
          query_input.index_origin;
      adapter_provenance["seqpro.query.verification"] =
          query_input.verification_status;
      reference = std::move(reference_input.fasta);
      query = std::move(query_input.fasta);
      if (reference_input.source_fasta_index_copied &&
          query_input.source_fasta_index_copied) {
        input_route = "seqpro-external-fai-mmap+caller-buffer-copy";
      } else if (reference_input.source_fasta_index_copied ||
                 query_input.source_fasta_index_copied) {
        input_route = "seqpro-hybrid-fai-mmap+caller-buffer-copy";
      } else {
        input_route = "seqpro-indexed-fasta-mmap+caller-buffer-copy";
      }
    } else {
      reference = ReadFasta(spec.reference_path, 0);
      query = ReadFasta(spec.query_path, 0);
      input_route = reference_compression == FastaCompression::Gzip ||
                            query_compression == FastaCompression::Gzip
                        ? "bundled-strict-zlib-fasta"
                        : "bundled-strict-fasta-oracle";
    }
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
      const auto index_begin = Clock::now();
      const bool load_persistent_index = !spec.reference_index_path.empty();
      stage = load_persistent_index ? "index-load" : "index-build";
      progress.Stage(stage, 0, 0,
                     load_persistent_index ? spec.reference_index_path.string()
                                           : "ephemeral");
      std::string index_creator_commit{kRequiredSufkitCommit};
      SufkitIndexOptions core_index_options{spec.threads};
#if RAMAG_USE_PAIRWISE_CORE
      core_index_options.mam_worker_cap = spec.threads;
      core_index_options.mam_workspace_limit_bytes = UINT64_MAX;
      core_index_options.boundary_mem_occurrence_limit = UINT64_MAX;
#endif
      auto index = [&]() {
        if (load_persistent_index) {
          index_creator_commit =
              ValidateReferenceIndexBundle(spec.reference_index_path, reference);
          return SufkitSeedIndex::Load(
              spec.reference_index_path, reference.sequences,
              core_index_options);
        }
        return SufkitSeedIndex::Build(
            reference.sequences, core_index_options);
      }();
      const double index_seconds = SecondsBetween(index_begin, Clock::now());
      timings[load_persistent_index ? "index_load" : "index_build"] =
          index_seconds;
      adapter_provenance["sufkit.index.action"] =
          load_persistent_index ? "loaded" : "built";
      adapter_provenance["sufkit.index.created_with_commit"] =
          index_creator_commit;
      adapter_provenance["sufkit.index.loaded_with_commit"] =
          std::string(kRequiredSufkitCommit);
      adapter_provenance["sufkit.index.path"] =
          load_persistent_index
              ? AbsolutePathString(spec.reference_index_path)
              : "ephemeral";
      adapter_provenance["sufkit.index.reference_validation"] =
          load_persistent_index ? "manifest+complete+metadata+fingerprint"
                                : "built-from-current-reference";
      if (load_persistent_index) {
        const auto bundle_paths =
            MakeReferenceIndexPaths(spec.reference_index_path);
        adapter_provenance["sufkit.index.bundle_manifest_sha256"] =
            FileSha256(bundle_paths.manifest);
        adapter_provenance["sufkit.index.bundle_complete_path"] =
            AbsolutePathString(bundle_paths.complete);
        adapter_provenance["sufkit.index.normalized_reference_sha256"] =
            NormalizedReferenceSha256(reference.sequences);
      }
      stage = "seed-enumeration";
      progress.Stage(stage, 0,
                     static_cast<std::uint64_t>(query.sequences.size()));
      auto seed_result = index.Enumerate(reference.sequences, query.sequences,
                                         alignment_options);
      RunStatistics seed_statistics;
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
#if RAMAG_USE_PAIRWISE_CORE
      result = AlignPairwiseFromSeeds(
#else
      result = AlignGenomesFromSeeds(
#endif
          reference.sequences, query.sequences, alignment_options,
          std::move(seed_result.seeds), std::move(seed_statistics));
    } else {
      stage = "seed-enumeration";
      progress.Stage(stage);
#if RAMAG_USE_PAIRWISE_CORE
      RunStatistics oracle_statistics;
      auto oracle_seeds = EnumerateSeeds(reference.sequences, query.sequences,
                                        alignment_options, &oracle_statistics);
      result = AlignPairwiseFromSeeds(reference.sequences, query.sequences,
                                   alignment_options, std::move(oracle_seeds),
                                   std::move(oracle_statistics));
#else
      result = AlignGenomes(reference.sequences, query.sequences,
                            alignment_options);
#endif
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
    const auto temporary_manifest = TemporarySibling(paths.manifest, run_id);
    const auto temporary_complete = TemporarySibling(paths.complete, run_id);
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

    stage = "manifest";
    progress.Stage(stage);
    ManifestData manifest;
    manifest.run_id = run_id;
    manifest.started_utc = started_utc;
    manifest.finished_utc = CurrentUtcTimestamp();
    manifest.run_spec = spec;
    manifest.reference = &reference;
    manifest.query = &query;
    manifest.statistics = result.statistics;
    manifest.artifacts = artifacts;
    manifest.stage_wall_seconds = timings;
    manifest.adapter_provenance = std::move(adapter_provenance);
    manifest.invocation = std::move(invocation);
    manifest.binary_path = std::move(binary_path);
    manifest.input_route = std::move(input_route);
    manifest.index_route = std::move(index_route);
    manifest.seeding_route = std::move(seeding_route);
    manifest.chaining_route = result.statistics.chaining_route;
#if RAMAG_USE_PAIRWISE_CORE
    manifest.extension_route = "pairwise-ksw2-certified-global+endpoint-extension";
#endif
    manifest.input_parallel_route = std::move(input_parallel_route);
    manifest.input_requested_workers = input_requested_workers;
    manifest.input_actual_workers = input_actual_workers;
    manifest.actual_threads = actual_threads;
    manifest.stage_wall_seconds["total_before_publish"] =
        SecondsBetween(total_begin, Clock::now());

    temporary_files.Add(temporary_manifest);
    const std::string manifest_json = ManifestJson(manifest);
    WriteFile(temporary_manifest, [&](std::ostream& output) {
      output << manifest_json;
    });
    ValidateManifestFile(temporary_manifest, run_id, manifest_json);
    publish.emplace_back(temporary_manifest, paths.manifest);

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
    ValidateManifestFile(paths.manifest, run_id, manifest_json);

    stage = "complete-marker";
    progress.Stage(stage);
    temporary_files.Add(temporary_complete);
    WriteFile(temporary_complete, [&](std::ostream& output) {
      output << "schema_version=1\n"
             << "run_id=" << run_id << "\n";
    });
    published_files.Add(temporary_complete, paths.complete);
    PublishNewFile(temporary_complete, paths.complete);
    progress.Finish(std::to_string(result.statistics.alignment_count) +
                    " alignments");
    published_files.Commit();

    return {paths, result.statistics};
  } catch (const std::exception& error) {
    WriteFailureDiagnostic(run_directory, run_id, stage, error);
    throw;
  }
}

}  // namespace ramag

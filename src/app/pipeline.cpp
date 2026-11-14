#include "centromere/refinement/refinement.hpp"
#include "app/pipeline.hpp"
#include "ramag/centromere_direction.hpp"
#include "ramag/centromere_blocks.hpp"
#include "ramag/centromere_directed_anchors.hpp"
#include "ramag/centromere_copy_anchors.hpp"
#include "ramag/centromere_stitch.hpp"
#include "ramag/centromere_stitch_formats.hpp"
#include "ramag/pairwise_core.hpp"

#include "ramag/alignment.hpp"
#include "io/fasta.hpp"
#include "index/reference_index.hpp"
#include "app/runtime.hpp"
#include "app/logging.hpp"
#include <map>
#include "index/sufkit_adapter.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdlib>
#include <cstdio>
#include <cerrno>
#include <ctime>
#include <exception>
#include <fstream>
#include <functional>
#include <iterator>
#include <iomanip>
#include <limits>
#include <new>
#include <optional>
#include <random>
#include <sstream>
#include <string_view>
#include <system_error>
#include <utility>
#include <vector>
#if defined(__GLIBC__)
#include <malloc.h>
#endif

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

CentromereBlockOptions DefaultCentromereDirectionOptions(
    const CentromereAlignmentOptions& alignment) {
  CentromereBlockOptions options;
  static_cast<CentromereAlignmentOptions&>(options) = alignment;
  // The direction profile used with the full alignment regression.
  options.planning_strategy = CentromereBlockStrategy::OrientationEvidenceV9;
  options.direction_island_audit = true;
  options.direction_island_context_review = true;
  options.direction_plan_review = true;
  options.overlap_compatibility = true;
  options.orientation_audit_known_blocks = true;
  options.orientation_bridge_unknown_gaps = true;
  options.orientation_connected_gaps = true;
  options.prefer_stronger_sparse_gap_end = true;
  options.recover_long_gap_ends = true;
  options.recover_mixed_gap_ends = true;
  options.recover_short_gap_ends = true;
  options.recover_sparse_gap_ends = true;
  options.recover_sparse_pivots = true;
  options.repair_drift_consistent_gaps = true;
  options.repair_inherited_directions = true;
  options.review_alternate_original_endpoints = true;
  options.review_collinear_fragment_connections = true;
  options.review_collinear_window_gaps = true;
  options.review_gap_interior_opposition = true;
  options.review_isolated_opposite_anchors = true;
  options.review_known_block_counterchains = true;
  options.review_local_competitive_gaps = true;
  options.review_neutral_window_gaps = true;
  options.review_region_chains = true;
  options.review_unified_gaps = true;
  options.supplement_region_anchors = true;
  options.unify_real_anchor_ends = true;
  options.unified_search_input_budget = 13'335'649;
  return options;
}

void VerifyPublishedArtifact(const ArtifactReport& artifact) {
  std::error_code error;
  const auto bytes = std::filesystem::file_size(artifact.path, error);
  if (error || bytes != artifact.bytes) {
    throw WriterError("published artifact size mismatch: " + artifact.path.string());
  }
}

}  // namespace

static RunOutcome RunAlignmentPipelineImpl(const RunSpec& spec,
                                std::string invocation,
                                std::filesystem::path binary_path, RunLogger* logger,
                                const FastaData* shared_reference,
                                const SufkitSeedIndex* shared_index,
                                std::string batch_label = {}) {
  ValidateRunSpec(spec);
  RequireCentromereBackend(spec.centromere);
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
  ProgressSession progress(spec.progress, run_id + batch_label, spec.threads, logger);

  try {
    progress.Stage("configuration");
    if(logger)logger->Info(spec.centromere.whole_input ? EffectiveConfigText(spec) :
                          PairwiseConfigurationText(spec.alignment),false);
    EnsureDirectory(spec.output_prefix.parent_path(), "output directory");
    if (!(spec.centromere.whole_input && (CentromereCliAnchorsOnly() ||
          CentromereCliFullAlignment())))
      EnsureTargetsAbsent(spec, paths);
    std::map<std::string, double> timings;
    stage = "input";
    progress.Stage(stage);
    const auto input_begin = Clock::now();
    FastaData owned_reference;
    const FastaData& reference = shared_reference ? *shared_reference : owned_reference;
    FastaData query;
    std::string input_route;
    std::string input_parallel_route = "serial-reference-query-v1";
    std::uint32_t input_requested_workers = 1;
    std::uint32_t input_actual_workers = 1;
    std::map<std::string, std::string> adapter_provenance;
    const auto reference_compression = shared_reference ? FastaCompression::Plain :
        DetectFastaCompression(spec.reference_path);
    const auto query_compression = DetectFastaCompression(spec.query_path);
    adapter_provenance["input.reference.compression"] =
        shared_reference ? "shared-normalized-reference" :
        (reference_compression == FastaCompression::Gzip ? "gzip" : "plain");
    adapter_provenance["input.query.compression"] =
        query_compression == FastaCompression::Gzip ? "gzip" : "plain";
    if (shared_reference) {
      query = ReadFasta(spec.query_path);
      input_parallel_route = "shared-reference-serial-query";
    } else {
      auto inputs = ReadFastaInputs(spec);
      input_requested_workers = inputs.requested_workers;
      input_actual_workers = inputs.actual_workers;
      input_parallel_route = std::move(inputs.route);
      owned_reference = std::move(inputs.results[0]);
      query = std::move(inputs.results[1]);
    }
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
    if (spec.centromere.whole_input) {
      if (reference.sequences.size()!=1 || query.sequences.size()!=1)
        throw CentromereInputError("align -c requires exactly one FASTA record on each side");
      CentromereAlignmentOptions cen_options;
      cen_options.worker_threads=spec.threads;
      cen_options.progress_callback=alignment_options.progress_callback;
      cen_options.interruption_callback=alignment_options.interruption_callback;
      if (CentromereCliFullAlignment()) {
        using Export = void(*)(std::ostream&, const CentromereDirectedAnchorResult&);
        const std::vector<std::pair<const char*, Export>> exports{
          {".directed-regions.tsv", WriteCentromereDirectedRegions},
          {".initial-anchors.tsv", WriteCentromereDirectedInitialAnchors},
          {".final-anchors.tsv", WriteCentromereDirectedFinalAnchors},
          {".recursion.tsv", WriteCentromereDirectedTrace},
          {".unanchored.tsv", WriteCentromereDirectedUncovered},
          {".reference-overlaps.tsv", WriteCentromereDirectedReferenceOverlaps},
          {".anchor-summary.tsv", WriteCentromereDirectedSummary},
          {".full-reference-initial-anchors.tsv", WriteCentromereDirectedFullReferenceInitialAnchors}};
        using CopyExport = void(*)(std::ostream&, const CentromereCopyAnchorResult&);
        const std::vector<std::pair<const char*, CopyExport>> copy_exports{
          {".copy-anchors.tsv", WriteCentromereCopyAnchors},
          {".copy-ambiguities.tsv", WriteCentromereCopyAmbiguities},
          {".copy-self-repeats.tsv", WriteCentromereCopySelfRepeats},
          {".copy-tasks.tsv", WriteCentromereCopyTasks},
          {".copy-summary.tsv", WriteCentromereCopySummary}};
        using StitchExport = void(*)(std::ostream&, const CentromereStitchResult&);
        const std::vector<std::pair<const char*, StitchExport>> stitch_exports{
          {".alignment-segments.tsv", WriteCentromereStitchSegments},
          {".alignment-gaps.tsv", WriteCentromereStitchGaps},
          {".alignment-summary.tsv", WriteCentromereStitchSummary}};
        using DirectionExport = void(*)(std::ostream&, const CentromereBlockResult&);
        const std::array<std::pair<const char*, DirectionExport>, 2> direction_exports{{
          {".orientation-local-direction-repair-intervals.tsv", WriteCentromereLocalDirectionRepairIntervals},
          {".direction-initial-anchors.tsv", WriteCentromereInitialAnchors}}};
        std::vector<std::filesystem::path> targets;
        for (const auto& item : exports)
          targets.emplace_back(spec.output_prefix.string() + item.first);
        for (const auto& item : copy_exports)
          targets.emplace_back(spec.output_prefix.string() + item.first);
        const auto stitch_offset = targets.size();
        for (const auto& item : stitch_exports)
          targets.emplace_back(spec.output_prefix.string() + item.first);
        using StandardExport = void(*)(std::ostream&, const SequenceRecord&,
                                        const SequenceRecord&, const CentromereStitchResult&);
        std::vector<std::pair<std::filesystem::path, StandardExport>> standard_exports;
        if (spec.formats.sam) standard_exports.emplace_back(paths.sam, WriteCentromereStitchSam);
        if (spec.formats.paf) standard_exports.emplace_back(paths.paf, WriteCentromereStitchPaf);
        if (spec.formats.maf) standard_exports.emplace_back(paths.maf, WriteCentromereStitchMaf);
        const auto standard_offset = targets.size();
        for (const auto& item : standard_exports) targets.push_back(item.first);
        const auto endpoint_offset = targets.size();
        for (const auto* suffix : {".unrefined-initial-anchors.tsv", ".unrefined-final-anchors.tsv", ".endpoint-reviews.tsv"})
          targets.emplace_back(spec.output_prefix.string() + suffix);
        const auto block_offset = targets.size();
        for (const auto* suffix : {".preblock-initial-anchors.tsv", ".preblock-final-anchors.tsv", ".block-reviews.tsv", ".preblock-alignment-segments.tsv"})
          targets.emplace_back(spec.output_prefix.string() + suffix);
        const auto cross_offset = targets.size();
        for (const auto* suffix : {".precross-initial-anchors.tsv", ".precross-final-anchors.tsv", ".cross-region-reviews.tsv", ".precross-alignment-segments.tsv"})
          targets.emplace_back(spec.output_prefix.string() + suffix);
        const auto joint_offset = targets.size();
        for (const auto* suffix : {".prejoint-initial-anchors.tsv", ".prejoint-final-anchors.tsv", ".joint-path-reviews.tsv", ".prejoint-alignment-segments.tsv", ".prejoint-directed-regions.tsv", ".prejoint-recursion.tsv"})
          targets.emplace_back(spec.output_prefix.string() + suffix);
        const auto unplanned_offset = targets.size();
        targets.emplace_back(spec.output_prefix.string() + ".unplanned-core-reviews.tsv");
        const auto direction_offset = targets.size();
        for (const auto& item : direction_exports)
          targets.emplace_back(spec.output_prefix.string() + item.first);
        targets.emplace_back(spec.output_prefix.string() + ".directed-run.tsv");
        for (const auto& target : targets)
          if (std::filesystem::symlink_status(target).type() != std::filesystem::file_type::not_found)
            throw WriterError("refusing to overwrite centromere output: " + target.string());

        const auto& r = reference.sequences.front();
        const auto& q = query.sequences.front();
        TemporaryFiles temporary;
        const auto direction_plan = TemporarySibling(targets[direction_offset], run_id);
        const auto initial_evidence = TemporarySibling(targets[direction_offset + 1], run_id);
        stage = "centromere-direction-plan";
        progress.Stage(stage);
        const auto direction_begin = Clock::now();
        {
          const auto evidence = FindCentromereDirectionRegions(r, q,
              DefaultCentromereDirectionOptions(cen_options));
          if (!evidence.complete || evidence.wfa_calls != 0 || !evidence.refined.empty() ||
              !evidence.alignment.alignments.empty() || evidence.refinement_searches != 0)
            throw AlignmentError("centromere direction prediction did not complete independently of alignment");
          for (std::size_t i = 0; i < direction_exports.size(); ++i) {
            CheckInterruption("centromere-direction-export");
            const auto tmp = TemporarySibling(targets[direction_offset + i], run_id);
            temporary.Add(tmp);
            WriteFile(tmp, [&](std::ostream& out) { direction_exports[i].second(out, evidence); });
          }
        }
        const double direction_seconds = SecondsBetween(direction_begin, Clock::now());
        const auto directions = ReadCentromereDirectedRegions(direction_plan);
        auto anchors = FindCentromereAnchorsInDirections(r, q, directions, cen_options);
        if (!anchors.complete || anchors.wfa_calls != 0)
          throw AlignmentError("fixed-direction anchor stage did not complete without WFA");

        stage = "copy-anchors";
        progress.Stage(stage);
        const auto copies = FindCentromereCopyAnchors(r, q, anchors, {}, cen_options);
        if (!copies.complete || copies.wfa_calls != 0)
          throw AlignmentError("copy anchor stage did not complete without WFA");

        CentromereEndpointOptions review_options;
        review_options.kernel.interruption_callback = alignment_options.interruption_callback;
        stage = "endpoint-review";
        progress.Stage(stage);
        const auto raw_anchors = anchors;
        auto endpoint_review = RefineCentromereRootEndpoints(r, q, anchors, &copies, review_options);
        anchors = std::move(endpoint_review.anchors);

        stage = "anchor-stitch";
        progress.Stage(stage);
        CentromereStitchOptions stitch_options;
        stitch_options.progress_callback = alignment_options.progress_callback;
        stitch_options.interruption_callback = alignment_options.interruption_callback;
        // Retain diagnostic checkpoints even if final publication is aborted.
        for (const auto& item : exports)
          WriteFile(run_directory / (std::string("stitch-frozen") + item.first),
            [&](std::ostream& out) { item.second(out, anchors); });
        for (const auto& item : copy_exports)
          WriteFile(run_directory / (std::string("stitch-frozen") + item.first),
            [&](std::ostream& out) { item.second(out, copies); });
        std::size_t stitch_observations = 0;
        stitch_options.observation_callback = [&](const CentromereStitchResult& observed) {
          if ((stitch_observations++ % 64) != 0 && !observed.complete && !observed.aborted) return;
          for (const auto& item : stitch_exports) {
            const auto target = run_directory / (std::string("stitch-snapshot") + item.first);
            const auto temporary_snapshot = target.string() + ".pending";
            WriteFile(temporary_snapshot, [&](std::ostream& out) { item.second(out, observed); });
            std::filesystem::rename(temporary_snapshot, target);
          }
        };
        auto stitched = AlignCentromereAnchoredPaths(r, q, anchors, &copies, stitch_options);

        stage = "block-review";
        progress.Stage(stage);
        const auto preblock_anchors = anchors;
        const auto preblock_stitched = stitched;
        auto block_review = RefineCentromereAnchorBlocks(r, q, anchors, &copies, stitched, review_options);
        if (block_review.accepted) {
          anchors = std::move(block_review.anchors);
          stage = "block-anchor-stitch";
          progress.Stage(stage);
          for (const auto& item : exports)
            WriteFile(run_directory / (std::string("block-stitch-frozen") + item.first),
              [&](std::ostream& out) { item.second(out, anchors); });
          stitched = AlignCentromereAnchoredPaths(r, q, anchors, &copies, stitch_options);
        }

        stage = "cross-region-review";
        progress.Stage(stage);
        const auto precross_anchors = anchors;
        const auto precross_stitched = stitched;
        auto cross_review = RefineCentromereCrossRegionBlocks(r, q, anchors, &copies, stitched, review_options);
        if (cross_review.accepted) {
          anchors = std::move(cross_review.anchors);
          stage = "cross-region-anchor-stitch";
          progress.Stage(stage);
          if (!cross_review.stitched)
            throw AlignmentError("accepted cross-region review lacks validated stitch");
          stitched = std::move(*cross_review.stitched);
        }

        stage = "joint-path-review";
        progress.Stage(stage);
        const auto prejoint_anchors = anchors;
        const auto prejoint_stitched = stitched;
        auto joint_review = RefineCentromereJointPaths(r, q, anchors, &copies, stitched, review_options);
        if (joint_review.accepted) {
          anchors = std::move(joint_review.anchors);
          if (!joint_review.stitched) throw AlignmentError("accepted joint review lacks validated stitch");
          stitched = std::move(*joint_review.stitched);
        }

        stage = "unplanned-core-review";
        progress.Stage(stage);
        CentromereEndpointResult unplanned_review;
        {
          auto options = review_options;
          std::vector<SequenceRecord> search_reference;
          std::optional<SufkitSeedIndex> search_index;
          std::function<void()> check_search_build;
          if (SufkitAdapterAvailable()) {
            options.unplanned_candidates = [&](Position begin, Position end,
                                                const std::function<void()> &check)
                -> std::optional<std::vector<Seed>> {
              const auto bytes = options.kernel.gap_memory_limit_bytes;
              // Reserve separate shares for index construction, enumeration
              // and the protected assembly. Large inputs retain the old plan.
              if (r.size() > bytes / 1024 || end - begin > 65536) return std::nullopt;
              try {
                if (!search_index) {
                  SufkitIndexOptions index_options;
                  index_options.boundary_mem_occurrence_limit = 100000;
                  index_options.mam_workspace_limit_bytes = bytes / 8;
                  check_search_build = check;
                  index_options.build_stage_context = &check_search_build;
                  index_options.build_stage_callback = [](const char*, void *context) {
                    (*static_cast<std::function<void()>*>(context))();
                  };
                  search_reference = {r};
                  check();
                  search_index.emplace(SufkitSeedIndex::Build(search_reference, index_options));
                  check_search_build = {};
                  if (search_index->BuildStatistics().resident_core_bytes > bytes / 8) {
                    search_index.reset();
                    return std::nullopt;
                  }
                  check();
                }
                AlignmentOptions search_options;
                search_options.seed_mode = SeedMode::MumReference;
                search_options.min_match = 128;
                search_options.interruption_callback = [&](std::string_view) { check(); };
                std::vector<SequenceRecord> query_window{
                    {q.numeric_id, q.name, q.header, q.bases.substr(begin, end - begin)}};
                auto found = search_index->Enumerate(search_reference, query_window, search_options);
                check();
                if (found.seeds.size() > 3000) return std::nullopt;
                return std::move(found.seeds);
              } catch (const std::bad_alloc&) {
                check_search_build = {};
                return std::nullopt;
              } catch (const AlignmentError &error) {
                check_search_build = {};
                const std::string_view reason = error.what();
                if (reason.starts_with("reference-MAM workspace ") ||
                    reason.starts_with("reference-MAM boundary MEM occurrences exceed limit "))
                  return std::nullopt;
                throw;
              }
            };
          }
          unplanned_review = RecoverCentromereUnplannedCores(r, q, anchors, &copies, stitched,
              initial_evidence, options);
        }
        if (unplanned_review.accepted) {
          anchors = std::move(unplanned_review.anchors);
          if (!unplanned_review.stitched)
            throw AlignmentError("accepted unplanned core review lacks validated stitch");
          stitched = std::move(*unplanned_review.stitched);
        }
        PublishedFiles published;
        for (std::size_t i = 0; i < targets.size(); ++i) {
          stage = "directed-anchor-export";
          CheckInterruption(stage);
          const auto tmp = TemporarySibling(targets[i], run_id);
          if (i >= direction_offset && i < direction_offset + direction_exports.size()) {
            MaybeInjectWriterFailure("direction-evidence");
            published.Add(tmp, targets[i]);
            PublishNewFile(tmp, targets[i]);
            continue;
          }
          temporary.Add(tmp);
          WriteFile(tmp, [&](std::ostream& out) {
            if (i < exports.size()) exports[i].second(out, anchors);
            else if (i < exports.size() + copy_exports.size())
              copy_exports[i - exports.size()].second(out, copies);
            else if (i >= stitch_offset && i < stitch_offset + stitch_exports.size())
              stitch_exports[i - stitch_offset].second(out, stitched);
            else if (i >= standard_offset && i < standard_offset + standard_exports.size())
              standard_exports[i - standard_offset].second(out, r, q, stitched);
            else if (i >= endpoint_offset && i < endpoint_offset + 3) {
              if (i == endpoint_offset) WriteCentromereDirectedInitialAnchors(out, raw_anchors);
              else if (i == endpoint_offset + 1) WriteCentromereDirectedFinalAnchors(out, raw_anchors);
              else WriteCentromereEndpointReviews(out, endpoint_review);
            }
            else if (i >= block_offset && i < block_offset + 4) {
              if (i == block_offset) WriteCentromereDirectedInitialAnchors(out, preblock_anchors);
              else if (i == block_offset + 1) WriteCentromereDirectedFinalAnchors(out, preblock_anchors);
              else if (i == block_offset + 2) WriteCentromereEndpointReviews(out, block_review);
              else WriteCentromereStitchSegments(out, preblock_stitched);
            }
            else if (i >= cross_offset && i < cross_offset + 4) {
              if (i == cross_offset) WriteCentromereDirectedInitialAnchors(out, precross_anchors);
              else if (i == cross_offset + 1) WriteCentromereDirectedFinalAnchors(out, precross_anchors);
              else if (i == cross_offset + 2) WriteCentromereCrossRegionReviews(out, cross_review);
              else WriteCentromereStitchSegments(out, precross_stitched);
            }
            else if (i >= joint_offset && i < joint_offset + 6) {
              if (i == joint_offset) WriteCentromereDirectedInitialAnchors(out, prejoint_anchors);
              else if (i == joint_offset + 1) WriteCentromereDirectedFinalAnchors(out, prejoint_anchors);
              else if (i == joint_offset + 2) WriteCentromereCrossRegionReviews(out, joint_review);
              else if (i == joint_offset + 3) WriteCentromereStitchSegments(out, prejoint_stitched);
              else if (i == joint_offset + 4) WriteCentromereDirectedRegions(out, prejoint_anchors);
              else WriteCentromereDirectedTrace(out, prejoint_anchors);
            }
            else if (i == unplanned_offset)
              WriteCentromereCrossRegionReviews(out, unplanned_review);
            else {
              out << std::setprecision(std::numeric_limits<double>::max_digits10)
                  << "field\tvalue\nstatus\t"
                  << (!stitched.complete ? "partial" : "complete")
                  << "\nwfa_calls\t0\n"
                  << "source_commit\t" << RAMAG_GIT_COMMIT << '\n'
                  << "version\t" << RAMAG_VERSION << '\n'
                  << "orientation_model\torientation-evidence-v9\n"
                  << "direction_plan\t" << std::filesystem::absolute(targets[direction_offset]).string() << '\n'
                  << "direction_initial_anchors\t" << std::filesystem::absolute(targets[direction_offset + 1]).string() << '\n'
                  << "direction_seconds\t" << direction_seconds << '\n'
                  << "reference\t" << std::filesystem::absolute(spec.reference_path).string() << '\n'
                  << "query\t" << std::filesystem::absolute(spec.query_path).string() << '\n'
                  << "reference_length\t" << r.size() << "\nquery_length\t" << q.size()
                  << "\nregions\t" << anchors.regions.size()
                  << "\ninitial_anchors\t" << anchors.initial.size()
                  << "\nfinal_anchors\t" << anchors.final.size()
                  << "\nthreads\t" << spec.threads
                  << "\nanchor_seconds\t" << anchors.total_seconds
                  << "\nreference_scope\tforward-flanked-reverse-window-v3\nrecursion_boundaries\tchecked-half-open-v1\n";
              out << "endpoint_model\tcompensating-one-sided-root-v1\nendpoint_seconds\t" << endpoint_review.seconds
                    << "\nendpoint_kernel_calls\t" << endpoint_review.kernel_calls
                    << "\nendpoint_accepted\t" << endpoint_review.accepted
                    << "\ninitial_anchor_stage\treviewed\nrecursion_trace_stage\tunrefined-search\n";
              out << "block_model\tbounded-context-block-v1\nblock_kernel\tlinear-7state-dp-9bit-two-strip-v1\nblock_seconds\t" << block_review.seconds
                    << "\nblock_kernel_calls\t" << block_review.kernel_calls
                    << "\nblock_accepted\t" << block_review.accepted
                    << "\nblock_pre_stitch_seconds\t" << preblock_stitched.total_seconds << '\n';
              out << "cross_region_model\tbounded-positive-context-v1\ncross_region_directions\tfrozen\ncross_region_score_baseline\tsum-existing-forward-fragments\ncross_region_seconds\t" << cross_review.seconds
                    << "\ncross_region_kernel_calls\t" << cross_review.kernel_calls
                    << "\ncross_region_accepted\t" << cross_review.accepted << '\n';
              out << "copy_model\tcopy-aware-v1\ncopy_seconds\t" << copies.total_seconds << '\n';
              out << "wfa_calls_scope\tanchor_stage\nanchor_status\tcomplete\n"
                    << "stitch_status\t" << (stitched.complete ? "complete" : "partial")
                    << "\nstitch_complete_scope\teligible_internal_gaps\nstitch_score_model\t"
                    << CentromereStitchScoreModelName(stitched.score_model)
                    << "\nstitch_dp_calls\t" << stitched.dp_calls
                    << "\nstitch_wavefront_calls\t" << stitched.wavefront_calls << '\n'
                    << "stitch_alignment_formats\t" << spec.formats.ToString() << '\n'
                    << "stitch_format_contract\tpath-segments-v1\n";
            }
          });
          if (i >= standard_offset && i < standard_offset + standard_exports.size())
            MaybeInjectWriterFailure("stitch-" + targets[i].extension().string().substr(1));
          if (i == endpoint_offset + 2)
            MaybeInjectWriterFailure("endpoint-review");
          if (i == block_offset + 2)
            MaybeInjectWriterFailure("block-review");
          if (i == cross_offset + 2)
            MaybeInjectWriterFailure("cross-region-review");
          MaybeInjectWriterFailure(i == stitch_offset
            ? "stitch-segments" : (i == 2 ? "directed-final-anchors" : "directed-anchors"));
          published.Add(tmp, targets[i]);
          PublishNewFile(tmp, targets[i]);
        }
        CheckInterruption("publication");
        RunStatistics stats;
        stats.reference_contigs = stats.query_contigs = 1;
        stats.reference_bases = r.size(); stats.query_bases = q.size();
        stats.selected_seed_count = stats.merged_seed_count = anchors.final.size();
        stats.chain_count = static_cast<std::uint64_t>(std::count_if(anchors.tasks.begin(), anchors.tasks.end(),
          [](const auto& task) { return task.final_anchor_count != 0; }));
        stats.seed_seconds = anchors.total_seconds;
        stats.total_seconds = SecondsBetween(total_begin, Clock::now());
        stats.actual_seed_route = "rama-fixed-direction-recursive-v1";
        stats.chaining_route = "rama-source-collinear-v1";
        stats.alignment_count = stitched.segments.size();
        stats.chain_and_extension_seconds = copies.total_seconds + stitched.total_seconds;
        progress.Finish(stitched.complete ? "supported internal gaps completed"
          : "partial internal gaps; unresolved intervals retained");
        if (logger) {
          const auto status = !stitched.complete ? "partial" : "success";
          logger->Info(std::string("status=") + status + " output_kind=branch-segments" +
            " wfa_calls=" + std::to_string(stitched.wavefront_calls) + " regions=" +
            std::to_string(anchors.regions.size()) + " initial_anchors=" + std::to_string(anchors.initial.size()) +
            " final_anchors=" + std::to_string(anchors.final.size()) + " total_seconds=" + std::to_string(stats.total_seconds));
          for (const auto& target : targets) logger->Info("anchor_output=" + target.string());
          logger->Flush();
        }
        OutputPaths exported_paths;
        if (spec.formats.sam) exported_paths.sam = paths.sam;
        if (spec.formats.paf) exported_paths.paf = paths.paf;
        if (spec.formats.maf) exported_paths.maf = paths.maf;
        RunOutcome outcome{std::move(exported_paths), std::move(stats)};
        published.Commit();
        return outcome;
      }
      if(CentromereCliBlockOptions().planning_strategy==CentromereBlockStrategy::CountClusterDirectionsV8 ||
         CentromereCliBlockOptions().planning_strategy==CentromereBlockStrategy::OrientationEvidenceV9) {
        auto options=CentromereCliBlockOptions();
        static_cast<CentromereAlignmentOptions&>(options)=cen_options;
        using Export=void(*)(std::ostream&,const CentromereBlockResult&);
        std::vector<std::pair<const char*,Export>> exports{
          {".initial-anchors.tsv",WriteCentromereInitialAnchors},{".root-candidates.tsv",WriteCentromereRootCandidates},
          {".blocks.tsv",WriteCentromereBlocks},{".planning-anchors.tsv",WriteCentromerePlanningAnchors}};
        if(options.planning_strategy==CentromereBlockStrategy::OrientationEvidenceV9)exports.insert(exports.end(),{
          {".orientation-windows.tsv",WriteCentromereOrientationWindows},{".orientation-evidence.tsv",WriteCentromereOrientationIntervals},
          {".orientation-sources.tsv",WriteCentromereOrientationSources},{".orientation-gaps.tsv",WriteCentromereOrientationGaps},
          {".orientation-unknown.tsv",WriteCentromereOrientationUnknown}});
        else exports.insert(exports.end(),{
          {".count-windows.tsv",WriteCentromereCountWindows},{".count-nodes.tsv",WriteCentromereCountNodes},
          {".count-candidates.tsv",WriteCentromereCountCandidates},{".count-gaps.tsv",WriteCentromereCountGaps},
          {".supported-intervals.tsv",WriteCentromereSupportedIntervals},{".direction-regions.tsv",WriteCentromereDirectionRegions}});
        if(options.direction_plan_review)exports.push_back({".orientation-plan-reviews.tsv",WriteCentromereDirectionPlanReviews});
        if(options.direction_local_protection)exports.push_back({".orientation-local-protection.tsv",WriteCentromereDirectionLocalProtectionReviews});
        if(options.direction_island_context_review)exports.push_back({".orientation-island-context-reviews.tsv",WriteCentromereDirectionIslandContextReviews});
        if(options.direction_island_audit)exports.push_back({".orientation-island-reviews.tsv",WriteCentromereDirectionIslandAudits});
        if(options.overlap_compatibility)exports.insert(exports.end(),{
          {".overlap-variants.tsv",WriteCentromereOverlapVariants},
          {".overlap-summary.tsv",WriteCentromereOverlapSummary}});
        if(options.orientation_audit_known_blocks)exports.insert(exports.end(),{
          {".orientation-block-audit.tsv",WriteCentromereOrientationBlockAudit},
          {".orientation-known-blocks.tsv",WriteCentromereOrientationKnownIntervals}});
        if(options.repair_inherited_directions)exports.insert(exports.end(),{
          {".orientation-local-direction-repair-reviews.tsv",WriteCentromereLocalDirectionRepairReviews},
          {".orientation-local-direction-repair-witnesses.tsv",WriteCentromereLocalDirectionRepairWitnesses},
          {".orientation-local-direction-repair-intervals.tsv",WriteCentromereLocalDirectionRepairIntervals},
          {".orientation-local-direction-repair-changes.tsv",WriteCentromereLocalDirectionRepairChanges}});
        if(options.review_known_block_counterchains)exports.insert(exports.end(),{
          {".orientation-known-counterchain-reviews.tsv",WriteCentromereKnownCounterchainReviews},
          {".orientation-known-counterchain-witnesses.tsv",WriteCentromereKnownCounterchainWitnesses},
          {".orientation-known-counterchain-intervals.tsv",WriteCentromereKnownCounterchainIntervals}});
        if(options.recover_long_gap_ends)exports.insert(exports.end(),{
          {".orientation-long-gap-ends.tsv",WriteCentromereOrientationLongGapEnds},
          {".orientation-long-gap-gaps.tsv",WriteCentromereOrientationLongGapGaps}});
        if(options.recover_sparse_pivots)exports.insert(exports.end(),{
          {".orientation-sparse-pivots.tsv",WriteCentromereOrientationSparsePivots},
          {".orientation-sparse-gaps.tsv",WriteCentromereOrientationSparseGaps},
          {".orientation-sparse-sources.tsv",WriteCentromereOrientationSparseSources}});
        if(options.recover_short_gap_ends)exports.insert(exports.end(),{
          {".orientation-short-gap-ends.tsv",WriteCentromereOrientationShortGapEnds},
          {".orientation-short-gap-gaps.tsv",WriteCentromereOrientationShortGapGaps}});
        if(options.recover_mixed_gap_ends)exports.insert(exports.end(),{
          {".orientation-mixed-gap-ends.tsv",WriteCentromereOrientationMixedGapEnds},
          {".orientation-mixed-gap-gaps.tsv",WriteCentromereOrientationMixedGapGaps}});
        if(options.recover_sparse_gap_ends)exports.insert(exports.end(),{
          {".orientation-sparse-gap-ends.tsv",WriteCentromereOrientationSparseGapEnds},
          {".orientation-sparse-gap-gaps.tsv",WriteCentromereOrientationSparseGapGaps}});
        if(options.unify_real_anchor_ends)exports.insert(exports.end(),{
          {".orientation-unified-pivots.tsv",WriteCentromereOrientationUnifiedPivots},
          {".orientation-unified-proposals.tsv",WriteCentromereOrientationUnifiedScans},
          {".orientation-unified-ends.tsv",WriteCentromereOrientationUnifiedEnds}});
        if(options.review_local_competitive_gaps)exports.insert(exports.end(),{
          {".orientation-local-reviews.tsv",WriteCentromereOrientationLocalReviews},
          {".orientation-local-gaps.tsv",WriteCentromereOrientationLocalGaps}});
        if(options.review_collinear_fragment_connections)exports.insert(exports.end(),{
          {".orientation-collinear-reviews.tsv",WriteCentromereOrientationCollinearReviews},
          {".orientation-collinear-counters.tsv",WriteCentromereOrientationCollinearCounters},
          {".orientation-collinear-witnesses.tsv",WriteCentromereOrientationCollinearWitnesses},
          {".orientation-collinear-gaps.tsv",WriteCentromereOrientationCollinearGaps}});
        if(options.review_alternate_original_endpoints)exports.insert(exports.end(),{
          {".orientation-alternate-reviews.tsv",WriteCentromereOrientationAlternateReviews},
          {".orientation-alternate-counters.tsv",WriteCentromereOrientationAlternateCounters},
          {".orientation-alternate-witnesses.tsv",WriteCentromereOrientationAlternateWitnesses},
          {".orientation-alternate-gaps.tsv",WriteCentromereOrientationAlternateGaps}});
        if(options.review_isolated_opposite_anchors)exports.insert(exports.end(),{
          {".orientation-isolated-reviews.tsv",WriteCentromereOrientationIsolatedReviews},
          {".orientation-opposite-reviews.tsv",WriteCentromereOrientationOppositeReviews},
          {".orientation-isolated-gaps.tsv",WriteCentromereOrientationIsolatedGaps}});
        for(const auto& [suffix,writer]:exports) {
          (void)writer;const auto target=spec.output_prefix.string()+suffix;
          if(std::filesystem::symlink_status(target).type()!=std::filesystem::file_type::not_found)
            throw WriterError("refusing to overwrite direction evidence: "+target);
        }
        const auto evidence=FindCentromereDirectionRegions(reference.sequences.front(),query.sequences.front(),options);
        if(evidence.wfa_calls || !evidence.refined.empty() || !evidence.alignment.alignments.empty() || evidence.refinement_searches)
          throw AlignmentError("direction-only strategy performed alignment or local search");
        TemporaryFiles temporary;PublishedFiles published;
        std::vector<std::pair<std::filesystem::path,std::filesystem::path>> publication;
        for(const auto& [suffix,writer]:exports) {
          const std::filesystem::path target=spec.output_prefix.string()+suffix;
          const auto tmp=TemporarySibling(target,run_id);temporary.Add(tmp);
          CheckInterruption("direction-export");WriteFile(tmp,[&](std::ostream& out){writer(out,evidence);});
          MaybeInjectWriterFailure("direction-evidence");publication.emplace_back(tmp,target);
        }
        progress.Stage("publication",0,publication.size(),"direction evidence");
        for(const auto& [tmp,target]:publication) {
          CheckInterruption("publication");published.Add(tmp,target);
          MaybeInjectPublicationRace(target);PublishNewFile(tmp,target);
        }
        CheckInterruption("publication");
        auto stats=evidence.alignment.statistics;
        stats.total_seconds=SecondsBetween(total_begin,Clock::now());
        if(logger) {
          if(options.repair_inherited_directions)logger->Info("repair_inherited_directions=1 repair_drift_consistent_gaps="+
            std::to_string(evidence.repair_drift_consistent_gaps)+" local_direction_repair_seconds="+
            std::to_string(evidence.local_direction_repair_seconds)+" local_direction_repair_connection_checks="+
            std::to_string(evidence.local_direction_repair_connection_checks)+" local_direction_repair_capacity_bytes="+
            std::to_string(evidence.local_direction_repair_capacity_bytes)+" local_direction_repair_reviews="+
            std::to_string(evidence.local_direction_repair_reviews.size())+" local_direction_repair_changes="+
            std::to_string(evidence.local_direction_repair_changes.size())+" local_direction_repair_intervals="+
            std::to_string(evidence.local_direction_repair_intervals.size())+
            " orientation_known_view=orientation-local-direction-repair-intervals.tsv");
          if(options.review_known_block_counterchains)logger->Info("review_known_block_counterchains=1 known_counterchain_seconds="+
            std::to_string(evidence.known_counterchain_seconds)+" known_counterchain_connection_checks="+
            std::to_string(evidence.known_counterchain_connection_checks)+" known_counterchain_capacity_bytes="+
            std::to_string(evidence.known_counterchain_capacity_bytes)+" known_counterchain_reviews="+
            std::to_string(evidence.known_counterchain_reviews.size())+" known_counterchain_intervals="+
            std::to_string(evidence.known_counterchain_intervals.size())+
            (options.repair_inherited_directions?" orientation_provenance_view=orientation-known-counterchain-intervals.tsv":
              " orientation_known_view=orientation-known-counterchain-intervals.tsv"));
          if(options.unify_real_anchor_ends)logger->Info("unify_real_anchor_ends=1 orientation_unified_seconds="+
            std::to_string(evidence.orientation_unified_seconds)+" orientation_unified_pivots="+
            std::to_string(evidence.orientation_unified_pivots.size())+" orientation_unified_intervals="+
            std::to_string(evidence.orientation_unified_intervals.size()));
          if(options.review_local_competitive_gaps)logger->Info("review_local_competitive_gaps=1 orientation_local_review_seconds="+
            std::to_string(evidence.orientation_local_review_seconds)+" orientation_local_reviews="+
            std::to_string(evidence.orientation_local_reviews.size()));
          if(options.review_collinear_window_gaps)logger->Info("review_collinear_window_gaps=1 policy=positive-window-break-full-audit-v1");
          if(options.review_neutral_window_gaps)logger->Info("review_neutral_window_gaps=1 policy=bounded-uncontested-window-full-review-v1");
          if(options.review_collinear_fragment_connections)logger->Info("review_collinear_fragment_connections=1 orientation_collinear_review_seconds="+
            std::to_string(evidence.orientation_collinear_review_seconds)+" orientation_collinear_reviews="+
            std::to_string(evidence.orientation_collinear_reviews.size()));
          if(options.review_isolated_opposite_anchors)logger->Info("review_isolated_opposite_anchors=1 orientation_isolated_review_seconds="+
            std::to_string(evidence.orientation_isolated_review_seconds)+" orientation_isolated_reviews="+
            std::to_string(evidence.orientation_isolated_reviews.size()));
          logger->Info("status=success output_kind=direction-evidence strategy="+std::string(CentromereBlockStrategyName(options.planning_strategy))+
            " wfa_calls=0 local_searches=0 records=0");
          logger->Info("count_window="+std::to_string(options.count_window)+" count_step="+std::to_string(options.count_step)+
            " count_min_anchors="+std::to_string(options.count_min_anchors)+" count_dominance_ratio="+std::to_string(options.count_dominance_ratio)+
            " count_cluster_gap="+std::to_string(options.count_cluster_gap)+" count_inferred_gap="+std::to_string(options.count_inferred_gap)+
            " orientation_connected_gaps="+std::to_string(options.orientation_connected_gaps)+
            " orientation_bridge_unknown_gaps="+std::to_string(options.orientation_bridge_unknown_gaps)+
            " orientation_audit_known_blocks="+std::to_string(options.orientation_audit_known_blocks)+
            " recover_long_gap_ends="+std::to_string(options.recover_long_gap_ends)+
            " orientation_long_gap_seconds="+std::to_string(evidence.orientation_long_gap_seconds)+
            " orientation_long_gap_intervals="+std::to_string(evidence.orientation_long_gap_intervals.size())+
            " recover_sparse_pivots="+std::to_string(options.recover_sparse_pivots)+
            " orientation_sparse_seconds="+std::to_string(evidence.orientation_sparse_seconds)+
            " orientation_sparse_intervals="+std::to_string(evidence.orientation_sparse_intervals.size())+
            " recover_short_gap_ends="+std::to_string(options.recover_short_gap_ends)+
            " orientation_short_gap_seconds="+std::to_string(evidence.orientation_short_gap_seconds)+
            " orientation_short_gap_intervals="+std::to_string(evidence.orientation_short_gap_intervals.size())+
            " recover_mixed_gap_ends="+std::to_string(options.recover_mixed_gap_ends)+
            " orientation_mixed_gap_seconds="+std::to_string(evidence.orientation_mixed_gap_seconds)+
            " orientation_mixed_gap_intervals="+std::to_string(evidence.orientation_mixed_gap_intervals.size())+
            " recover_sparse_gap_ends="+std::to_string(options.recover_sparse_gap_ends)+
            (options.prefer_stronger_sparse_gap_end?" prefer_stronger_sparse_gap_end=1":"")+
            " orientation_sparse_gap_seconds="+std::to_string(evidence.orientation_sparse_gap_seconds)+
            " orientation_sparse_gap_intervals="+std::to_string(evidence.orientation_sparse_gap_intervals.size())+
            " orientation_audit_seconds="+std::to_string(evidence.orientation_audit_seconds)+
            " orientation_retractions="+std::to_string(evidence.orientation_block_audits.size())+
            " initial_seconds="+std::to_string(evidence.initial_seconds)+" planning_seconds="+std::to_string(evidence.planning_seconds)+
            " orientation_seconds="+std::to_string(evidence.orientation_seconds)+
            " total_seconds="+std::to_string(stats.total_seconds));logger->Flush();
        }
        progress.Finish("direction evidence only; no recursion or WFA");
        published.Commit();return {{},stats};
      }
      if (CentromereCliBlocks()) {
        const auto& r=reference.sequences.front();
        const auto& q=query.sequences.front();
        auto options=CentromereCliBlockOptions();
        std::vector<std::string> suffixes{".initial-anchors.tsv",".blocks.tsv",".refined-anchors.tsv",".unresolved.tsv"};
        using Export=void(*)(std::ostream&,const CentromereBlockResult&);
        std::vector<Export> exports{WriteCentromereInitialAnchors,WriteCentromereBlocks,
          WriteCentromereRefinedAnchors,WriteCentromereUnresolved};
        if(options.planning_strategy!=CentromereBlockStrategy::InitialDirectionBlocksV1) {
          suffixes.insert(suffixes.end(),{".root-candidates.tsv",".planning-anchors.tsv"});
          exports.insert(exports.end(),{WriteCentromereRootCandidates,WriteCentromerePlanningAnchors});
        }
        if(options.planning_strategy==CentromereBlockStrategy::HierarchicalDirectionRegionsV4 ||
           options.planning_strategy==CentromereBlockStrategy::DirectionRegionsOnlyV5) {
          suffixes.insert(suffixes.end(),{".direction-regions.tsv",".refinement-regions.tsv",".local-candidates.tsv"});
          exports.insert(exports.end(),{WriteCentromereDirectionRegions,WriteCentromereRefinementRegions,WriteCentromereLocalCandidates});
        }
        std::vector<std::filesystem::path> evidence_paths(suffixes.size());
        for (std::size_t i=0;i<suffixes.size();++i) {
          evidence_paths[i]=spec.output_prefix.string()+suffixes[i];
          if (std::filesystem::symlink_status(evidence_paths[i]).type()!=std::filesystem::file_type::not_found)
            throw WriterError("refusing to overwrite block evidence: "+evidence_paths[i].string());
        }
        const auto partial_path=[](const std::filesystem::path& path) {
          return path.parent_path()/(path.stem().string()+".partial"+path.extension().string());
        };
        const std::array<std::pair<bool,std::filesystem::path>,5> selected{{
          {spec.formats.sam,paths.sam},{spec.formats.paf,paths.paf},{spec.formats.delta,paths.delta},
          {spec.formats.maf,paths.maf},{spec.formats.chain,paths.chain}}};
        for (const auto& [enabled,path]:selected) if (enabled &&
          std::filesystem::symlink_status(partial_path(path)).type()!=std::filesystem::file_type::not_found)
            throw WriterError("refusing to overwrite partial output: "+partial_path(path).string());
        static_cast<CentromereAlignmentOptions&>(options)=cen_options;
        std::size_t observations=0;
        bool plan_observed=false;
        options.observation_callback=[&](const CentromereBlockResult& observed) {
          // Diagnostic snapshots are not formal published results. Preserve the
          // frozen plan even when a later correctness error aborts publication.
          const bool first_plan=!plan_observed && !observed.blocks.empty();
          if(!observed.blocks.empty())plan_observed=true;
          if ((observations++%64)!=0 && !first_plan) return;
          for (std::size_t i=0;i<exports.size();++i)
            WriteFile(run_directory/("block-snapshot"+suffixes[i]),[&](std::ostream& out){exports[i](out,observed);});
        };
        auto evidence=AlignCentromeresWithDirections(r,q,options);
        TemporaryFiles temporary;
        PublishedFiles published;
        std::vector<std::pair<std::filesystem::path,std::filesystem::path>> publish;
        const auto stage_file=[&](const std::filesystem::path& target,const auto& writer,std::string_view kind) {
          const auto tmp=TemporarySibling(target,run_id);temporary.Add(tmp);
          CheckInterruption("block-export");WriteFile(tmp,writer);MaybeInjectWriterFailure(kind);
          publish.emplace_back(tmp,target);return tmp;
        };
        for (std::size_t i=0;i<exports.size();++i)
          stage_file(evidence_paths[i],[&](std::ostream& out){exports[i](out,evidence);},"block-evidence");
        const std::array<std::string_view,5> names{"sam","paf","delta","maf","chain"};
        const std::array<void(*)(const std::filesystem::path&),5> validators{
          ValidateSamFile,ValidatePafFile,ValidateDeltaFile,ValidateMafFile,ValidateChainFile};
        OutputPaths final_paths=paths;
        std::array<std::filesystem::path*,5> returned{&final_paths.sam,&final_paths.paf,&final_paths.delta,&final_paths.maf,&final_paths.chain};
        const auto output_start=Clock::now();
        for (std::size_t i=0;i<selected.size();++i) if (selected[i].first) {
          const auto target=evidence.complete ? selected[i].second : partial_path(selected[i].second);
          *returned[i]=target;
          const auto tmp=stage_file(target,[&](std::ostream& out) {
            const auto& a=evidence.alignment.alignments;
            switch(i) {
              case 0:WriteSam(out,reference,query,a,invocation);break;
              case 1:WritePaf(out,reference,query,a);break;
              case 2:WriteDelta(out,reference,query,a);break;
              case 3:WriteMaf(out,reference,query,a);break;
              case 4:WriteChain(out,reference,query,a);break;
            }
          },names[i]);
          MaybeInjectValidatorFailure(names[i]);validators[i](tmp);
        }
        for (const auto& [tmp,target]:publish) {
          CheckInterruption("publication");published.Add(tmp,target);
          MaybeInjectPublicationRace(target);PublishNewFile(tmp,target);
        }
        CheckInterruption("publication");
        auto& stats=evidence.alignment.statistics;stats.total_seconds=SecondsBetween(total_begin,Clock::now());
        if(logger) {
          logger->Info(std::string("status=")+(evidence.complete ? "success exit_code=0" : "partial exit_code=9")+
            " strategy="+std::string(CentromereBlockStrategyName(evidence.strategy))+" initial_forward="+std::to_string(evidence.initial_forward.size())+
            " initial_reverse="+std::to_string(evidence.initial_reverse.size())+" blocks="+std::to_string(evidence.blocks.size())+
            " refined_anchors="+std::to_string(evidence.refined.size())+" records="+std::to_string(stats.alignment_count)+
            " unresolved_or_unpaired="+std::to_string(evidence.unresolved.size())+" wfa_calls="+std::to_string(evidence.wfa_calls)+
            " root_candidates="+std::to_string(evidence.root_candidates.size())+
            " planner_reserved_bytes="+std::to_string(evidence.planner_reserved_bytes)+
            " planner_capacity_bytes="+std::to_string(evidence.planner_capacity_bytes)+
            " local_chain_states="+std::to_string(evidence.local_chain_states)+
            " local_chain_transitions="+std::to_string(evidence.local_chain_transitions)+
            " local_chain_dominated="+std::to_string(evidence.local_chain_dominated)+
            " direction_regions="+std::to_string(evidence.regions.size())+
            " refinement_rounds="+std::to_string(evidence.refinement_rounds)+
            " refinement_searches="+std::to_string(evidence.refinement_searches)+
            " refinement_new_mappings="+std::to_string(evidence.refinement_new_mappings)+
            " refinement_partial_reason="+evidence.refinement_partial_reason);
          logger->Info("initial_seconds="+std::to_string(evidence.initial_seconds)+" planning_seconds="+std::to_string(evidence.planning_seconds)+
            " recursion_seconds="+std::to_string(evidence.recursion_seconds)+" wfa_seconds="+std::to_string(evidence.wfa_seconds)+
            " output_seconds="+std::to_string(SecondsBetween(output_start,Clock::now()))+" total_seconds="+std::to_string(stats.total_seconds));
          logger->Flush();
        }
        progress.Finish(evidence.complete ? "fixed-direction block alignment complete" : "partial alignment; see unresolved intervals");
        published.Commit();return {final_paths,stats,evidence.complete ? 0 : 9};
      }
      if (CentromereCliDirections()) {
        const auto& r=reference.sequences.front();
        const auto& q=query.sequences.front();
        const std::array<std::string,5> suffixes{".direction-anchors.tsv",".direction-regions.tsv",
          ".direction-chains.tsv",".direction-tasks.tsv",".inversion-candidates.tsv"};
        std::array<std::filesystem::path,5> targets;
        for (std::size_t i=0;i<targets.size();++i) {
          targets[i]=spec.output_prefix.string()+suffixes[i];
          if (std::filesystem::symlink_status(targets[i]).type()!=std::filesystem::file_type::not_found)
            throw WriterError("refusing to overwrite direction evidence: "+targets[i].string());
        }
        CentromereDirectionOptions options;
        options.search_strategy=CentromereCliDirectionStrategy();
        static_cast<CentromereAlignmentOptions&>(options)=cen_options;
        if(logger)options.memory_diagnostic_callback=[&](std::string_view detail) {
          logger->Info("direction-memory "+std::string(detail));
        };
        auto evidence=AnalyzeCentromereDirections(r,q,options);
        using Export=void(*)(std::ostream&,const CentromereDirectionResult&);
        const std::array<Export,5> writers{WriteDirectionAnchors,WriteDirectionRegions,
          WriteDirectionChains,WriteDirectionTasks,WriteInversionEvidence};
        TemporaryFiles temporary;
        PublishedFiles published;
        for (std::size_t i=0;i<targets.size();++i) {
          stage="direction-export";CheckInterruption(stage);
          const auto tmp=TemporarySibling(targets[i],run_id);temporary.Add(tmp);
          WriteFile(tmp,[&](std::ostream& out){writers[i](out,evidence);});
          MaybeInjectWriterFailure(i==0 ? "direction-anchors" : "direction-evidence");
          published.Add(tmp,targets[i]);PublishNewFile(tmp,targets[i]);
        }
        CheckInterruption("publication");
        evidence.statistics.total_seconds=SecondsBetween(total_begin,Clock::now());
        progress.Finish("bidirectional evidence only; WFA not executed");
        if (logger) {
          logger->Info("status=success output_kind=direction-evidence wfa_calls=0 strategy="+
            std::string(DirectionStrategyName(options.search_strategy))+" objective=effective-collinear-length-v1 actual_workers=1 requested_workers="+
            std::to_string(spec.threads)+" anchors="+std::to_string(evidence.anchors.size())+
            " regions="+std::to_string(evidence.regions.size())+" candidates="+std::to_string(evidence.candidates.size())+
            " reference="+r.name+" reference_length="+std::to_string(r.size())+" query="+q.name+
            " query_length="+std::to_string(q.size())+" total_seconds="+std::to_string(evidence.statistics.total_seconds));
          for (const auto& target:targets) logger->Info("direction_output="+target.string());
          logger->Flush();
        }
        published.Commit();return {{},evidence.statistics};
      }
      {  // Anchors-only diagnostic builds.
        const auto& r=reference.sequences.front();
        const auto& q=query.sequences.front();
        const std::array<std::filesystem::path,2> targets{
          spec.output_prefix.string()+".initial-anchors.tsv",
          spec.output_prefix.string()+".final-anchors.tsv"};
        for (const auto& target:targets)
          if (std::filesystem::symlink_status(target).type()!=std::filesystem::file_type::not_found)
            throw WriterError("refusing to overwrite anchor output: "+target.string());
        auto anchors=FindCentromereAnchors(r,q,cen_options);
        TemporaryFiles temporary;
        PublishedFiles published;
        const std::array<const std::vector<CentromereAnchor>*,2> sets{&anchors.initial,&anchors.final};
        for (std::size_t i=0;i<targets.size();++i) {
          stage="anchor-export";
          CheckInterruption(stage);
          const auto tmp=TemporarySibling(targets[i],run_id);
          temporary.Add(tmp);
          WriteFile(tmp,[&](std::ostream& out) {
            out<<std::setprecision(std::numeric_limits<double>::max_digits10);
            out<<"anchor_id\treference_id\treference_length\treference_begin\treference_end\tquery_id\tquery_length\tquery_begin\tquery_end\tstrand\tmatch_length\tweight\n";
            std::size_t id=0;
            for (const auto& a:*sets[i]) {
              if ((id%65536)==0) CheckInterruption(stage);
              out<<++id<<'\t'<<r.name<<'\t'<<r.size()<<'\t'<<a.reference_begin<<'\t'
                 <<a.reference_begin+a.length<<'\t'<<q.name<<'\t'<<q.size()<<'\t'
                 <<a.query_begin<<'\t'<<a.query_begin+a.length<<"\t+\t"<<a.length<<'\t'<<a.weight<<'\n';
            }
          });
          MaybeInjectWriterFailure(i==0 ? "initial-anchors" : "final-anchors");
          published.Add(tmp,targets[i]);
          PublishNewFile(tmp,targets[i]);
        }
        CheckInterruption("publication");
        progress.Finish("anchors only; WFA not executed");
        anchors.statistics.total_seconds=SecondsBetween(total_begin,Clock::now());
        if (logger) {
          logger->Info("status=success output_kind=anchors-only wfa_calls=0 initial_anchors="+
            std::to_string(anchors.initial.size())+" final_anchors="+std::to_string(anchors.final.size())+
            " total_seconds="+std::to_string(anchors.statistics.total_seconds));
          for (const auto& target:targets) logger->Info("anchor_output="+target.string());
          logger->Flush();
        }
        published.Commit();
        return {{},anchors.statistics};
      }
    } else if (SufkitAdapterAvailable()) {
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
      stage = shared_index ? "index-reuse" : (load_persistent_index ? "index-load" : "index-build");
      progress.Stage(stage, 0, 0,
                     load_persistent_index ? spec.reference_index_path.string()
                                           : "ephemeral");
      SufkitIndexOptions core_index_options{spec.threads};
      core_index_options.mam_worker_cap = spec.threads;
      core_index_options.mam_workspace_limit_bytes = UINT64_MAX;
      core_index_options.boundary_mem_occurrence_limit = UINT64_MAX;
      std::string last_load_stage;
      core_index_options.load_stage_callback = [&](std::string_view phase) {
        CheckInterruption("index-load");
        if (phase != last_load_stage) {
          last_load_stage = phase;
          progress.Stage("index-load:" + last_load_stage);
        }
      };
      std::optional<SufkitSeedIndex> owned_index;
      if (!shared_index) owned_index.emplace([&]() {
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
      }());
      const auto& index = shared_index ? *shared_index : *owned_index;
      const double index_seconds = SecondsBetween(index_begin, Clock::now());
      timings[shared_index ? "index_reuse" : (load_persistent_index ? "index_load" : "index_build")] =
          index_seconds;
      adapter_provenance["sufkit.index.save_path"] = spec.save_index_path.string();
      if (load_persistent_index && !shared_index) {
        const auto& loaded = index.BuildStatistics();
        for (const auto& [phase, seconds] : loaded.load_stage_seconds)
          timings["index_load_" + phase] = seconds;
        timings["sufkit_load_total"] = loaded.sufkit_load_seconds;
        timings["index_reference_validation"] = loaded.reference_validation_seconds;
        if (logger) logger->Info("index_action=loaded requested_threads=" +
            std::to_string(spec.threads) + " lcp_encoding=" + loaded.lcp_encoding +
            " logical_read_bytes=" + std::to_string(loaded.load_logical_read_bytes) +
            " crc_cpu_seconds=" + std::to_string(loaded.load_crc_seconds), false);
      }
      adapter_provenance["sufkit.index.persisted"] = "false";
      if (!spec.save_index_path.empty() && !shared_index) {
        stage = "index-save";
        IndexSpec save_spec{spec.reference_path, spec.save_index_path,
                            spec.work_dir, spec.threads, spec.progress};
        const auto saved = SaveReferenceIndex(save_spec, reference, index,
            invocation, binary_path, progress, index_seconds, timings, logger);
        retained_index = saved.index;
        adapter_provenance["sufkit.index.persisted"] = "true";
      }
      adapter_provenance["sufkit.index.action"] =
          shared_index ? "reused" : (load_persistent_index ? "loaded" : "built");
      if (shared_index && logger) logger->Info("index_action=reused index_load_calls=0 index_build_calls=0", false);
      adapter_provenance["sufkit.index.loaded_with_commit"] =
          std::string(kRequiredSufkitCommit);
      adapter_provenance["sufkit.index.path"] =
          load_persistent_index
              ? AbsolutePathString(spec.reference_index_path)
              : (spec.save_index_path.empty() ? "ephemeral" : AbsolutePathString(spec.save_index_path));
      adapter_provenance["sufkit.index.reference_validation"] =
          shared_index ? "validated-once-at-batch-initialization" :
          load_persistent_index ? "sufkit-crc+structure+metadata+fingerprint"
                                : "built-from-current-reference";
      stage = "seed-enumeration";
      progress.Stage(stage, 0,
                     static_cast<std::uint64_t>(query.sequences.size()));
      observe_index("index-ready", index.BuildStatistics().resident_core_bytes);
      seed_result = shared_index
          ? index.EnumerateOwned(query.sequences, alignment_options)
          : index.Enumerate(reference.sequences, query.sequences, alignment_options);
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
      observe_index(shared_index ? "index-retained-for-batch" : "index-released",
                    shared_index ? shared_index->BuildStatistics().resident_core_bytes : 0);
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
                 result.statistics.conflict_resolution_seconds -
                 result.statistics.mainpath_link.seconds);

    // Alignment-core conflict resolution is the sole authority for
    // de-duplication, primary assignment, and deterministic record order.
    // Writers consume this immutable result and never select independently.
    timings["conflict_resolution"] =
        result.statistics.conflict_resolution_seconds;

    timings["mainpath_link"] = result.statistics.mainpath_link.seconds;
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
      logger->Info("mainpath_link.input_records=" + std::to_string(result.statistics.mainpath_link.input_records), false);
      logger->Info("mainpath_link.output_records=" + std::to_string(result.statistics.mainpath_link.output_records), false);
      logger->Info("mainpath_link.adjacent_pairs=" + std::to_string(result.statistics.mainpath_link.adjacent_pairs), false);
      logger->Info("mainpath_link.geometry_rejected=" + std::to_string(result.statistics.mainpath_link.geometry_rejected), false);
      logger->Info("mainpath_link.flank_rejected=" + std::to_string(result.statistics.mainpath_link.flank_rejected), false);
      logger->Info("mainpath_link.n_rejected=" + std::to_string(result.statistics.mainpath_link.n_rejected), false);
      logger->Info("mainpath_link.length_rejected=" + std::to_string(result.statistics.mainpath_link.length_rejected), false);
      logger->Info("mainpath_link.budget_skipped=" + std::to_string(result.statistics.mainpath_link.budget_skipped), false);
      logger->Info("mainpath_link.quality_rejected=" + std::to_string(result.statistics.mainpath_link.quality_rejected), false);
      logger->Info("mainpath_link.closure_rejected=" + std::to_string(result.statistics.mainpath_link.closure_rejected), false);
      logger->Info("mainpath_link.memory_rejected=" + std::to_string(result.statistics.mainpath_link.memory_rejected), false);
      logger->Info("mainpath_link.zero_gap_links=" + std::to_string(result.statistics.mainpath_link.zero_gap_links), false);
      logger->Info("mainpath_link.exact_gap_links=" + std::to_string(result.statistics.mainpath_link.exact_gap_links), false);
      logger->Info("mainpath_link.single_gap_links=" + std::to_string(result.statistics.mainpath_link.single_gap_links), false);
      logger->Info("mainpath_link.dp_gap_links=" + std::to_string(result.statistics.mainpath_link.dp_gap_links), false);
      logger->Info("mainpath_link.added_reference_bases=" + std::to_string(result.statistics.mainpath_link.added_reference_bases), false);
      logger->Info("mainpath_link.added_query_bases=" + std::to_string(result.statistics.mainpath_link.added_query_bases), false);
      logger->Info("mainpath_link.added_paired_columns=" + std::to_string(result.statistics.mainpath_link.added_paired_columns), false);
      logger->Info("mainpath_link.work_budget=" + std::to_string(result.statistics.mainpath_link.work_budget), false);
      logger->Info("mainpath_link.work_reserved=" + std::to_string(result.statistics.mainpath_link.work_reserved), false);
      logger->Info("mainpath_link.ksw_calls=" + std::to_string(result.statistics.mainpath_link.ksw_calls), false);
      logger->Info("mainpath_link.ksw_seconds=" + std::to_string(result.statistics.mainpath_link.ksw_seconds), false);
      logger->Info("mainpath_link.planning_seconds=" + std::to_string(result.statistics.mainpath_link.planning_seconds), false);
      logger->Info("mainpath_link.assembly_seconds=" + std::to_string(result.statistics.mainpath_link.assembly_seconds), false);
      logger->Info("mainpath_link.seconds=" + std::to_string(result.statistics.mainpath_link.seconds), false);
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
      logger->Info("pairwise.link_input_candidates=" + std::to_string(result.statistics.pairwise.link_input_candidates), false);
      logger->Info("pairwise.link_consumed_candidates=" + std::to_string(result.statistics.pairwise.link_consumed_candidates), false);
      logger->Info("pairwise.link_emitted_candidates=" + std::to_string(result.statistics.pairwise.link_emitted_candidates), false);
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

RunOutcome RunAlignmentPipeline(const RunSpec& spec, std::string invocation,
    std::filesystem::path binary_path, RunLogger* logger) {
  return RunAlignmentPipelineImpl(spec, std::move(invocation),
      std::move(binary_path), logger, nullptr, nullptr);
}

BatchOutcome RunBatchPipeline(const BatchSpec& spec, std::string invocation,
    std::filesystem::path binary_path, RunLogger* logger) {
  ValidateBatchSpec(spec);
  if (!SufkitAdapterAvailable()) throw DependencyError("batch requires the pinned Sufkit adapter");
  EnsureDirectory(spec.output_dir, "batch output directory");
  // Exclusive creation prevents concurrent batches from claiming the same table.
  const auto close_table = [](std::FILE* file) noexcept { (void)std::fclose(file); };
  std::unique_ptr<std::FILE, decltype(close_table)> table(
      std::fopen((spec.output_dir / "batch.tsv").c_str(), "wx"), close_table);
  if (!table) throw WriterError("cannot exclusively create batch.tsv");
  const auto write = [&](const std::string& text) {
    if (std::fwrite(text.data(), 1, text.size(), table.get()) != text.size() ||
        std::fflush(table.get()) != 0) throw WriterError("cannot reliably write batch.tsv");
  };
  const auto escape = [](std::string_view text) {
    std::string result;
    for (const char c : text) {
      if (c == '\\') result += "\\\\";
      else if (c == '\t') result += "\\t";
      else if (c == '\n') result += "\\n";
      else if (c == '\r') result += "\\r";
      else result += c;
    }
    return result;
  };
  write("order\tname\tquery_path\tstatus\texit_code\telapsed_seconds\talignment_records\toutput_dir\tlog_path\tmessage\n");
  BatchOutcome outcome;
  for (const auto& query : spec.queries) {
    BatchQueryOutcome entry;
    entry.name = query.name;
    outcome.queries.push_back(std::move(entry));
  }
  std::size_t recorded = 0;
  const auto record = [&](std::size_t i) {
    const auto& item = outcome.queries[i];
    std::ostringstream row;
    row << i + 1 << '\t' << item.name << '\t' << escape(spec.queries[i].path.string())
        << '\t' << item.status << '\t';
    if (item.status != "not-run") row << item.exit_code;
    row << '\t' << item.elapsed_seconds << '\t' << item.alignment_count << '\t'
        << escape((spec.output_dir / item.name).string()) << '\t'
        << escape(item.log_path.string()) << '\t' << escape(item.message) << '\n';
    write(row.str());
  };
  const auto batch_begin = Clock::now();
  try {
    // No writable alias survives reference initialization. The index aliases
    // the sequence vector in this owner, so no reference bases are copied.
    std::shared_ptr<const FastaData> reference_owner;
    std::optional<SufkitSeedIndex> index;
    {
      ProgressSession progress(spec.common.progress, logger ? logger->Id() : MakeNonce(), spec.common.threads, logger);
      progress.Stage("batch-reference-input", 0, spec.queries.size());
      const auto begin = Clock::now();
      reference_owner = std::make_shared<const FastaData>(ReadFasta(spec.common.reference_path));
      const auto& reference = *reference_owner;
      const std::shared_ptr<const std::vector<SequenceRecord>> records(
          reference_owner, &reference.sequences);
      if (logger) logger->Info("batch shared_reference_input_seconds=" + std::to_string(SecondsBetween(begin, Clock::now())));
      const auto index_begin = Clock::now();
      SufkitIndexOptions options{spec.common.threads};
      options.mam_worker_cap = spec.common.threads;
      options.mam_workspace_limit_bytes = UINT64_MAX;
      options.boundary_mem_occurrence_limit = UINT64_MAX;
      // Captures outlive the index, including any callbacks retained by adapter.
      options.load_stage_callback = [](std::string_view phase) { CheckInterruption(phase); };
      if (spec.common.reference_index_path.empty()) {
        progress.Stage("batch-index-build");
        ValidateCpuAffinityBudget(CurrentCpuAffinity(), spec.common.threads, "batch index build");
        ++outcome.index_build_calls;
        index.emplace(SufkitSeedIndex::BuildOwned(records, options));
      } else {
        progress.Stage("batch-index-load");
        ++outcome.index_load_calls;
        index.emplace(SufkitSeedIndex::LoadOwned(spec.common.reference_index_path, records, options));
      }
      const auto seconds = SecondsBetween(index_begin, Clock::now());
      if (logger) {
        const auto stats = index->BuildStatistics();
        logger->Info("batch index_action=" + std::string(outcome.index_load_calls ? "loaded" : "built") +
            " index_load_calls=" + std::to_string(outcome.index_load_calls) +
            " index_build_calls=" + std::to_string(outcome.index_build_calls) +
            " shared_index_seconds=" + std::to_string(seconds) + " lcp_encoding=" + stats.lcp_encoding);
        for (const auto& [phase, value] : stats.load_stage_seconds)
          logger->Info("batch timing stage=" + phase + " seconds=" + std::to_string(value), false);
        logger->Info("batch index retained across pairwise stages; estimated_resident_bytes=" + std::to_string(stats.resident_core_bytes), false);
        logger->Info("batch timing stage=reference-validation seconds=" + std::to_string(stats.reference_validation_seconds) +
            " sufkit_load_seconds=" + std::to_string(stats.sufkit_load_seconds), false);
      }
      if (!spec.common.save_index_path.empty()) {
        IndexSpec save{spec.common.reference_path, spec.common.save_index_path, spec.common.work_dir,
                       spec.common.threads, spec.common.progress};
        std::map<std::string, double> timings;
        (void)SaveReferenceIndex(save, reference, *index, invocation, binary_path, progress, seconds, timings, logger);
      }
    }
    if (logger) logger->Info("batch scheduling=serial active_query_limit=1 reference_validation=once; threads apply within each query");
    for (std::size_t i = 0; i < spec.queries.size(); ++i) {
      // Previous query objects have been destroyed before admitting the next.
      // Return freed allocator pages on glibc; never touch the live index.
#if defined(__GLIBC__)
      (void)malloc_trim(0);
#endif
      CheckInterruption("batch-query-boundary");
      auto& item = outcome.queries[i];
      const auto run = BatchQuerySpec(spec, spec.queries[i]);
      const auto begin = Clock::now();
      std::unique_ptr<RunLogger> child_logger;
      try {
        if (logger) {
          logger->Info("batch query=" + item.name + " position=" + std::to_string(i + 1) + "/" + std::to_string(spec.queries.size()));
          child_logger = std::make_unique<RunLogger>(run.work_dir);
          item.log_path = child_logger->Path();
          child_logger->Info("batch query=" + item.name + " invocation=" + invocation, false);
        }
        if (!std::filesystem::is_regular_file(run.query_path)) throw FastaError("query is not a readable regular file: " + run.query_path.string());
        const auto result = RunAlignmentPipelineImpl(run, invocation, binary_path,
            child_logger.get(), reference_owner.get(), &*index,
            ":query=" + item.name + ":" + std::to_string(i + 1) + "/" + std::to_string(spec.queries.size()));
        item.alignment_count = result.statistics.alignment_count;
        item.status = "success";
      } catch (const std::exception& error) {
        item.exit_code = FailureExitCode(error);
        item.status = dynamic_cast<const InterruptedError*>(&error) ? "interrupted" : "failed";
        item.message = error.what();
        item.elapsed_seconds = SecondsBetween(begin, Clock::now());
        record(i); ++recorded;
        if (logger) logger->Error("batch query=" + item.name + " status=" + item.status + " message=" + item.message);
        if (dynamic_cast<const InterruptedError*>(&error) || dynamic_cast<const std::bad_alloc*>(&error)) throw;
        if (outcome.exit_code == 0) outcome.exit_code = item.exit_code;
        continue;
      }
      item.elapsed_seconds = SecondsBetween(begin, Clock::now());
      record(i); ++recorded;
    }
    if (logger) {
      logger->Info("batch status=" + std::string(outcome.exit_code ? "failed" : "success") +
          " exit_code=" + std::to_string(outcome.exit_code) + " total_seconds=" +
          std::to_string(SecondsBetween(batch_begin, Clock::now())));
      logger->Flush();
    }
  } catch (const std::exception& error) {
    for (std::size_t i = recorded; i < outcome.queries.size(); ++i) {
      outcome.queries[i].message = "batch stopped: " + std::string(error.what());
      record(i);
    }
    throw;
  }
  return outcome;
}

}  // namespace ramag

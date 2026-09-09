#include "ramag/reference_index.hpp"

#include "ramag/logging.hpp"
#include "ramag/sufkit_adapter.hpp"
#include "ramag/writers.hpp"
#include "ramag/runtime.hpp"

#include <algorithm>
#include <memory>
#include <optional>
#include <cstdlib>
#include <chrono>
#include <exception>
#include <fstream>
#include <iomanip>
#include <random>
#include <sstream>
#include <system_error>
#include <utility>
#include <vector>

#if defined(__unix__) || defined(__APPLE__)
#include <unistd.h>
#endif

namespace ramag {
namespace {

using Clock = std::chrono::steady_clock;

void IndexFailureProbe(std::string_view stage) {
  const char* requested = std::getenv("RAMAG_TEST_FAIL_INDEX_STAGE");
  if (requested && stage == requested) throw WriterError("injected index failure: " + std::string(stage));
}

std::string Nonce() {
  std::random_device device;
  std::mt19937_64 generator(device());
  std::ostringstream output;
#if defined(__unix__) || defined(__APPLE__)
  output << static_cast<unsigned long long>(getpid()) << '-';
#endif
  output << std::hex << generator();
  return output.str();
}

void PublishNoReplace(const std::filesystem::path& temporary,
                      const std::filesystem::path& target) {
  std::error_code error;
  std::filesystem::create_hard_link(temporary, target, error);
  if (error) {
    throw WriterError("cannot publish reference index artifact without overwrite: " +
                      target.string() + ": " + error.message());
  }
}

std::filesystem::path Temporary(const std::filesystem::path& target,
                                std::string_view nonce) {
  return std::filesystem::path(target.string() + ".partial." + std::string(nonce));
}

}  // namespace

ReferenceIndexPaths MakeReferenceIndexPaths(const std::filesystem::path& index) {
  return {index};
}

static ReferenceIndexPaths RunReferenceIndexPipelineImpl(
    const IndexSpec& spec, std::string invocation,
    const std::filesystem::path& binary_path,
    const CpuAffinityInfo& launch_affinity,
    const FastaData* supplied_reference, const SufkitSeedIndex* supplied_index,
    ProgressSession* supplied_progress, double supplied_build_seconds,
    std::map<std::string, double>* timings, RunLogger* logger) {
  ValidateIndexSpec(spec, supplied_index == nullptr);
  const auto paths = MakeReferenceIndexPaths(spec.output_path);
  for (const auto& path : {paths.index}) {
    std::error_code error;
    if (std::filesystem::exists(path, error)) {
      throw WriterError("refusing to overwrite existing index artifact: " + path.string());
    }
    if (error) throw WriterError("cannot inspect index artifact: " + error.message());
  }
  std::error_code error;
  if (!paths.index.parent_path().empty()) {
    std::filesystem::create_directories(paths.index.parent_path(), error);
    if (error) throw WriterError("cannot create index output directory: " + error.message());
  }
  std::filesystem::create_directories(spec.work_dir, error);
  if (error) throw WriterError("cannot create index work directory: " + error.message());

  const auto nonce = Nonce();
  std::unique_ptr<ProgressSession> owned_progress;
  if (!supplied_progress) owned_progress = std::make_unique<ProgressSession>(spec.progress, nonce, spec.threads, logger);
  auto& progress = supplied_progress ? *supplied_progress : *owned_progress;
  const auto temporary_index = Temporary(paths.index, nonce);
  std::vector<std::filesystem::path> temporaries{
      temporary_index};
  std::vector<std::pair<std::filesystem::path, std::filesystem::path>> published;
  bool bundle_committed = false;
  std::string stage{"configuration"};
  try {
    const auto started = Clock::now();
    stage = "input";
    if (!supplied_reference) progress.Stage(stage);
    std::optional<FastaData> owned_reference;
    if (!supplied_reference) owned_reference.emplace(ReadFasta(spec.reference_path, 0));
    const auto& reference = supplied_reference ? *supplied_reference : *owned_reference;
    const auto pre_sufkit_affinity = CurrentCpuAffinity();
    ValidateCpuAffinityBudget(pre_sufkit_affinity, spec.threads,
                              "Sufkit index build");
    if (launch_affinity.supported &&
        pre_sufkit_affinity.logical_cpus != launch_affinity.logical_cpus) {
      throw CliError(
          "CPU affinity changed before Sufkit index build "
          "(launch=" + launch_affinity.CpuList() +
          ", pre_sufkit=" + pre_sufkit_affinity.CpuList() + ")");
    }
    SufkitIndexOptions index_options{spec.threads};
    index_options.build_stage_context = &progress;
    index_options.build_stage_callback = [](const char* phase, void* context) {
      // Observability is best effort; cancellation is checked by the pipeline.
      try {
        static_cast<ProgressSession*>(context)->Stage(phase);
      } catch (...) {
      }
    };
    const bool expect_caps =
        spec.threads > 1 &&
        reference.total_bases >=
            index_options.parallel_caps_min_reference_bases;
    std::ostringstream build_detail;
    build_detail << "backend=" << (expect_caps ? "caps" : "divsufsort")
                 << " requested_threads=" << spec.threads
                 << " allowed_cpu_count="
                 << pre_sufkit_affinity.logical_cpus.size()
                 << " allowed=" << pre_sufkit_affinity.CpuList();
    const auto build_begin = Clock::now();
    stage = "index-build";
    if (!supplied_index) progress.Stage(stage, 0, 0, build_detail.str());
    std::optional<SufkitSeedIndex> owned_index;
    if (!supplied_index) owned_index.emplace(SufkitSeedIndex::Build(reference.sequences, index_options));
    const auto& index = supplied_index ? *supplied_index : *owned_index;
    const double build_seconds = supplied_index ? supplied_build_seconds :
        std::chrono::duration<double>(Clock::now() - build_begin).count();
    const auto stats = index.BuildStatistics();
    if (expect_caps && !stats.backend.starts_with("caps")) {
      throw DependencyError(
          "large-reference parallel index policy selected CaPS, but Sufkit "
          "reported backend '" + stats.backend + "'");
    }
    stage = "index-save";
    progress.Stage(stage);
    const auto save_begin = Clock::now();
    index.Save(temporary_index);
    IndexFailureProbe("save");
    const double save_seconds =
        std::chrono::duration<double>(Clock::now() - save_begin).count();
    stage = "index-validation";
    progress.Stage(stage);
    const auto validation_begin = Clock::now();
    static_cast<void>(SufkitSeedIndex::Load(
        temporary_index, reference.sequences, index_options));
    IndexFailureProbe("validation");
    const double validation_seconds =
        std::chrono::duration<double>(Clock::now() - validation_begin).count();
    const auto serialized_bytes = std::filesystem::file_size(temporary_index);
    stage = "publication";
    progress.Stage(stage, 0, 1);
    const auto publication_begin = Clock::now();
    published.emplace_back(temporary_index, paths.index);
    PublishNoReplace(temporary_index, paths.index);
    IndexFailureProbe("publication");
    const double index_publication_seconds =
        std::chrono::duration<double>(Clock::now() - publication_begin).count();
    progress.Update(1, 1);
    if (!supplied_progress) progress.Finish(std::to_string(serialized_bytes) + " bytes");
    if (logger) {
      logger->Info("index_action=built persisted=true index=" + paths.index.string() + " backend=" + stats.backend + " lcp_encoding=" + stats.lcp_encoding + " index_bytes=" + std::to_string(serialized_bytes));
      logger->Info("index_build_seconds=" + std::to_string(build_seconds) + " sufkit_build_seconds=" + std::to_string(stats.sufkit_build_seconds) + " index_save_seconds=" + std::to_string(save_seconds) + " index_validation_seconds=" + std::to_string(validation_seconds) + " index_publication_seconds=" + std::to_string(index_publication_seconds) + " total_seconds=" + std::to_string(std::chrono::duration<double>(Clock::now() - started).count()), false);
      logger->Info("binary=" + binary_path.string() + " invocation=" + invocation, false);
      logger->Flush();
    }
    CheckInterruption("index-publication");
    bundle_committed = true;
    if (timings) {
      (*timings)["index_save"] = save_seconds;
      (*timings)["index_validation"] = validation_seconds;
      (*timings)["index_publication"] = std::chrono::duration<double>(Clock::now() - publication_begin).count();
    }
    for (const auto& temporary : temporaries) {
      std::error_code ignored;
      std::filesystem::remove(temporary, ignored);
    }
    return paths;
  } catch (...) {
    try {
      if (logger) {
        try { throw; } catch (const std::exception& failure) {
          logger->Error("status=failed stage=" + stage + " exit_code=" + std::to_string(FailureExitCode(failure)) + " message=" + failure.what());
        }
      }
    } catch (...) {}
    for (auto it = published.rbegin(); !bundle_committed && it != published.rend(); ++it) {
      std::error_code ignored;
      if (std::filesystem::equivalent(it->first, it->second, ignored) && !ignored) {
        std::filesystem::remove(it->second, ignored);
      }
    }
    for (const auto& path : temporaries) {
      std::error_code ignored;
      std::filesystem::remove(path, ignored);
    }
    throw;
  }
}

ReferenceIndexPaths RunReferenceIndexPipeline(
    const IndexSpec& spec, std::string invocation,
    const std::filesystem::path& binary_path,
    const CpuAffinityInfo& launch_affinity, RunLogger* logger) {
  return RunReferenceIndexPipelineImpl(spec, std::move(invocation), binary_path,
      launch_affinity, nullptr, nullptr, nullptr, 0.0, nullptr, logger);
}

ReferenceIndexPaths SaveReferenceIndex(
    const IndexSpec& spec, const FastaData& reference,
    const SufkitSeedIndex& index, std::string invocation,
    const std::filesystem::path& binary_path, ProgressSession& progress,
    double build_seconds, std::map<std::string, double>& timings, RunLogger* logger) {
  return RunReferenceIndexPipelineImpl(spec, std::move(invocation), binary_path,
      CurrentCpuAffinity(), &reference, &index, &progress, build_seconds, &timings, logger);
}

}  // namespace ramag

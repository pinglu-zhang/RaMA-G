#include "ramag/reference_index.hpp"

#include "ramag/sha256.hpp"
#include "ramag/sufkit_adapter.hpp"
#include "ramag/writers.hpp"
#include "ramag/runtime.hpp"

#include <algorithm>
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

std::string UtcTimestamp() {
  const auto now = std::chrono::system_clock::now();
  const auto time = std::chrono::system_clock::to_time_t(now);
  std::tm value{};
#if defined(_WIN32)
  gmtime_s(&value, &time);
#else
  gmtime_r(&time, &value);
#endif
  std::ostringstream output;
  output << std::put_time(&value, "%Y-%m-%dT%H:%M:%SZ");
  return output.str();
}

std::string Escape(std::string_view text) {
  std::ostringstream output;
  for (const unsigned char value : text) {
    switch (value) {
      case '\\': output << "\\\\"; break;
      case '"': output << "\\\""; break;
      case '\n': output << "\\n"; break;
      case '\r': output << "\\r"; break;
      case '\t': output << "\\t"; break;
      default:
        if (value < 0x20U) {
          output << "\\u" << std::hex << std::setw(4) << std::setfill('0')
                 << static_cast<unsigned int>(value) << std::dec;
        } else {
          output << static_cast<char>(value);
        }
    }
  }
  return output.str();
}

void WriteText(const std::filesystem::path& path, std::string_view content) {
  std::ofstream output(path, std::ios::binary | std::ios::trunc);
  if (!output) throw WriterError("cannot create index companion: " + path.string());
  output << content;
  output.close();
  if (!output) throw WriterError("cannot finish index companion: " + path.string());
}

std::string ReadText(const std::filesystem::path& path) {
  std::ifstream input(path, std::ios::binary);
  if (!input) throw DependencyError("missing reference index companion: " + path.string());
  std::ostringstream output;
  output << input.rdbuf();
  if (input.bad()) throw DependencyError("cannot read reference index companion: " + path.string());
  return output.str();
}

void RequireContains(std::string_view content, std::string_view expected,
                     const std::filesystem::path& path) {
  if (content.find(expected) == std::string_view::npos) {
    throw DependencyError("reference index companion does not contain required identity '" +
                          std::string(expected) + "': " + path.string());
  }
}

std::string MarkerValue(std::string_view content, std::string_view key,
                        const std::filesystem::path& path) {
  const std::string prefix = std::string(key) + "=";
  std::size_t begin = 0;
  std::string value;
  bool found = false;
  while (begin <= content.size()) {
    const auto end = content.find('\n', begin);
    const auto line = content.substr(
        begin, end == std::string_view::npos ? content.size() - begin : end - begin);
    if (line.starts_with(prefix)) {
      if (found) {
        throw DependencyError("duplicate reference index marker key '" +
                              std::string(key) + "': " + path.string());
      }
      value = std::string(line.substr(prefix.size()));
      found = true;
    }
    if (end == std::string_view::npos) break;
    begin = end + 1;
  }
  if (!found || value.empty()) {
    throw DependencyError("missing reference index marker key '" +
                          std::string(key) + "': " + path.string());
  }
  return value;
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
  return {index, std::filesystem::path(index.string() + ".manifest.json"),
          std::filesystem::path(index.string() + ".complete")};
}

std::string ValidateReferenceIndexBundle(const std::filesystem::path& index,
                                  const FastaData& reference) {
  const auto paths = MakeReferenceIndexPaths(index);
  const auto manifest = ReadText(paths.manifest);
  const auto complete = ReadText(paths.complete);
  const auto digest = NormalizedReferenceSha256(reference.sequences);
  const auto creator = MarkerValue(complete, "sufkit_commit", paths.complete);
  // Format 1.4 / Fast full-SA compatibility is tested in both directions.
  constexpr std::string_view previous =
      "bdb67c6de5daddd8a005640de73d96549d2575f4";
  if (creator != kRequiredSufkitCommit && creator != previous) {
    throw DependencyError("reference index was created by an incompatible Sufkit commit");
  }
  if (MarkerValue(complete, "schema_version", paths.complete) != "1" ||
      MarkerValue(complete, "normalized_reference_sha256", paths.complete) !=
          digest ||
      MarkerValue(complete, "manifest_sha256", paths.complete) !=
          Sha256Hex(manifest)) {
    throw DependencyError(
        "reference index complete marker does not match its manifest or reference");
  }
  RequireContains(manifest, "\"schema_version\": 1", paths.manifest);
  RequireContains(manifest, "\"status\": \"success\"", paths.manifest);
  RequireContains(manifest,
                  "\"index_kind\": \"standalone_suffix_array\"",
                  paths.manifest);
  RequireContains(manifest, "\"resource_profile\": \"fast\"",
                  paths.manifest);
  RequireContains(manifest, "\"sampling_rate\": 1", paths.manifest);
  RequireContains(manifest, "\"acceleration\": \"suffix-link\"",
                  paths.manifest);
  RequireContains(manifest, "\"sufkit_commit\": \"" +
                    creator + "\"", paths.manifest);
  RequireContains(manifest, "\"sufkit_version\": \"0.3.0\"", paths.manifest);
  RequireContains(manifest, "\"normalized_reference_sha256\": \"" + digest + "\"",
                  paths.manifest);
  RequireContains(manifest,
                  "\"reference_contigs\": " +
                      std::to_string(reference.sequences.size()),
                  paths.manifest);
  RequireContains(manifest,
                  "\"reference_bases\": " +
                      std::to_string(reference.total_bases),
                  paths.manifest);
  RequireContains(manifest,
                  "\"reference_ambiguous_bases\": " +
                      std::to_string(reference.ambiguous_bases),
                  paths.manifest);
  for (std::size_t ordinal = 0; ordinal < reference.sequences.size(); ++ordinal) {
    const auto& record = reference.sequences[ordinal];
    const auto ambiguous = static_cast<std::uint64_t>(
        std::count(record.bases.begin(), record.bases.end(), 'N'));
    std::ostringstream catalog_entry;
    catalog_entry << "{\"ordinal\": " << ordinal << ", \"name\": \""
                  << Escape(record.name) << "\", \"header\": \""
                  << Escape(record.header) << "\", \"length\": "
                  << record.size() << ", \"ambiguous_bases\": "
                  << ambiguous << "}";
    RequireContains(manifest, catalog_entry.str(), paths.manifest);
  }
  return creator;
}

ReferenceIndexPaths RunReferenceIndexPipeline(
    const IndexSpec& spec, std::string invocation,
    const std::filesystem::path& binary_path,
    const CpuAffinityInfo& launch_affinity) {
  ValidateIndexSpec(spec);
  const auto paths = MakeReferenceIndexPaths(spec.output_path);
  for (const auto& path : {paths.index, paths.manifest, paths.complete}) {
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
  ProgressSession progress(spec.progress, nonce, spec.threads);
  const auto temporary_index = Temporary(paths.index, nonce);
  const auto temporary_manifest = Temporary(paths.manifest, nonce);
  const auto temporary_complete = Temporary(paths.complete, nonce);
  std::vector<std::filesystem::path> temporaries{
      temporary_index, temporary_manifest, temporary_complete};
  std::vector<std::pair<std::filesystem::path, std::filesystem::path>> published;
  std::string stage{"configuration"};
  try {
    const auto started = Clock::now();
    stage = "input";
    progress.Stage(stage);
    const auto reference_compression =
        DetectFastaCompression(spec.reference_path);
    auto reference = ReadFasta(spec.reference_path, 0);
    const auto digest = NormalizedReferenceSha256(reference.sequences);
    const auto source_digest = FileSha256(spec.reference_path);
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
    progress.Stage(stage, 0, 0, build_detail.str());
    auto index = SufkitSeedIndex::Build(
        reference.sequences, index_options);
    const double build_seconds =
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
    const double save_seconds =
        std::chrono::duration<double>(Clock::now() - save_begin).count();
    stage = "index-validation";
    progress.Stage(stage);
    const auto validation_begin = Clock::now();
    static_cast<void>(SufkitSeedIndex::Load(
        temporary_index, reference.sequences, index_options));
    const double validation_seconds =
        std::chrono::duration<double>(Clock::now() - validation_begin).count();
    const auto serialized_bytes = std::filesystem::file_size(temporary_index);
    const auto source_bytes = std::filesystem::file_size(spec.reference_path);
    const auto binary_bytes = std::filesystem::file_size(binary_path);
    const auto binary_digest = FileSha256(binary_path);

    stage = "publication";
    progress.Stage(stage, 0, 3);
    const auto publication_begin = Clock::now();
    published.emplace_back(temporary_index, paths.index);
    PublishNoReplace(temporary_index, paths.index);
    const double index_publication_seconds =
        std::chrono::duration<double>(Clock::now() - publication_begin).count();
    progress.Update(1, 3);

    std::ostringstream manifest;
    manifest << "{\n"
             << "  \"schema_version\": 1,\n"
             << "  \"status\": \"success\",\n"
             << "  \"created_utc\": \"" << UtcTimestamp() << "\",\n"
             << "  \"sufkit_commit\": \"" << kRequiredSufkitCommit << "\",\n"
             << "  \"sufkit_version\": \"" << stats.library_version << "\",\n"
             << "  \"index_kind\": \"standalone_suffix_array\",\n"
             << "  \"index_format\": \"" << stats.format_version << "\",\n"
             << "  \"backend\": \"" << Escape(stats.backend) << "\",\n"
             << "  \"resource_profile\": \"" << Escape(stats.resource_profile) << "\",\n"
             << "  \"sampling_rate\": " << stats.sampling_rate << ",\n"
             << "  \"acceleration\": \"" << Escape(stats.acceleration) << "\",\n"
             << "  \"lookup_acceleration\": \""
             << Escape(stats.lookup_acceleration) << "\",\n"
             << "  \"lcp_encoding\": \"" << Escape(stats.lcp_encoding) << "\",\n"
             << "  \"reference_path\": \"" << Escape(reference.source.string()) << "\",\n"
             << "  \"reference_compression\": \""
             << (reference_compression == FastaCompression::Gzip ? "gzip" : "plain")
             << "\",\n"
             << "  \"reference_source_bytes\": " << source_bytes << ",\n"
             << "  \"reference_source_sha256\": \"" << source_digest << "\",\n"
             << "  \"normalized_reference_sha256\": \"" << digest << "\",\n"
             << "  \"reference_fingerprint\": " << stats.reference_fingerprint << ",\n"
             << "  \"reference_contigs\": " << reference.sequences.size() << ",\n"
             << "  \"reference_bases\": " << reference.total_bases << ",\n"
             << "  \"reference_ambiguous_bases\": " << reference.ambiguous_bases << ",\n"
             << "  \"index_bytes\": " << serialized_bytes << ",\n"
             << "  \"threads\": " << spec.threads << ",\n"
             << "  \"requested_threads\": " << spec.threads << ",\n"
             << "  \"launch_allowed_cpu_count\": "
             << launch_affinity.logical_cpus.size() << ",\n"
             << "  \"launch_allowed_cpu_list\": \""
             << Escape(launch_affinity.CpuList()) << "\",\n"
             << "  \"launch_allowed_cpu_source\": \""
             << Escape(launch_affinity.source) << "\",\n"
             << "  \"pre_sufkit_allowed_cpu_count\": "
             << pre_sufkit_affinity.logical_cpus.size() << ",\n"
             << "  \"pre_sufkit_allowed_cpu_list\": \""
             << Escape(pre_sufkit_affinity.CpuList()) << "\",\n"
             << "  \"index_backend\": \"" << Escape(stats.backend) << "\",\n"
             << "  \"build_seconds\": " << build_seconds << ",\n"
             << "  \"sufkit_build_seconds\": "
             << stats.sufkit_build_seconds << ",\n"
             << "  \"suffix_array_seconds\": "
             << stats.suffix_array_seconds << ",\n"
             << "  \"isa_seconds\": " << stats.isa_seconds << ",\n"
             << "  \"caps_construct_seconds\": " << stats.caps_construct_seconds << ",\n"
             << "  \"caps_output_allocation_seconds\": " << stats.caps_output_allocation_seconds << ",\n"
             << "  \"text_prepare_seconds\": " << stats.text_prepare_seconds << ",\n"
             << "  \"lcp_finalize_seconds\": " << stats.lcp_finalize_seconds << ",\n"
             << "  \"prefix_directory_seconds\": " << stats.prefix_directory_seconds << ",\n"
             << "  \"lcp_seconds\": " << stats.lcp_seconds << ",\n"
             << "  \"index_save_seconds\": " << save_seconds << ",\n"
             << "  \"index_validation_seconds\": "
             << validation_seconds << ",\n"
             << "  \"index_publication_seconds\": "
             << index_publication_seconds << ",\n"
             << "  \"total_seconds\": "
             << std::chrono::duration<double>(Clock::now() - started).count() << ",\n"
             << "  \"ramag_commit\": \"" << RAMAG_GIT_COMMIT << "\",\n"
             << "  \"binary_path\": \"" << Escape(binary_path.string()) << "\",\n"
             << "  \"binary_bytes\": " << binary_bytes << ",\n"
             << "  \"binary_sha256\": \"" << binary_digest << "\",\n"
             << "  \"invocation\": \"" << Escape(invocation) << "\",\n"
             << "  \"contigs\": [\n";
    for (std::size_t ordinal = 0; ordinal < reference.sequences.size(); ++ordinal) {
      const auto& record = reference.sequences[ordinal];
      const auto ambiguous = static_cast<std::uint64_t>(
          std::count(record.bases.begin(), record.bases.end(), 'N'));
      manifest << "    {\"ordinal\": " << ordinal << ", \"name\": \""
               << Escape(record.name) << "\", \"header\": \""
               << Escape(record.header) << "\", \"length\": " << record.size()
               << ", \"ambiguous_bases\": " << ambiguous << "}"
               << (ordinal + 1 == reference.sequences.size() ? "\n" : ",\n");
    }
    manifest << "  ]\n}\n";
    const auto manifest_text = manifest.str();
    WriteText(temporary_manifest, manifest_text);
    std::ostringstream complete;
    complete << "schema_version=1\n"
             << "sufkit_commit=" << kRequiredSufkitCommit << "\n"
             << "normalized_reference_sha256=" << digest << "\n"
             << "manifest_sha256=" << Sha256Hex(manifest_text) << "\n";
    WriteText(temporary_complete, complete.str());

    published.emplace_back(temporary_manifest, paths.manifest);
    PublishNoReplace(temporary_manifest, paths.manifest);
    progress.Update(2, 3);
    published.emplace_back(temporary_complete, paths.complete);
    PublishNoReplace(temporary_complete, paths.complete);
    progress.Update(3, 3);
    progress.Finish(std::to_string(serialized_bytes) + " bytes");
    for (const auto& temporary : temporaries) {
      std::error_code ignored;
      std::filesystem::remove(temporary, ignored);
    }
    return paths;
  } catch (...) {
    try {
      const auto failure = std::current_exception();
      std::filesystem::create_directories(spec.work_dir / "runs" / nonce);
      try {
        std::rethrow_exception(failure);
      } catch (const InterruptedError& interrupted) {
        std::ostringstream diagnostic;
        diagnostic << "{\n  \"schema_version\": 1,\n"
                   << "  \"status\": \"interrupted\",\n"
                   << "  \"run_id\": \"" << nonce << "\",\n"
                   << "  \"observed_utc\": \"" << UtcTimestamp() << "\",\n"
                   << "  \"signal\": " << interrupted.SignalNumber() << ",\n"
                   << "  \"exit_code\": " << interrupted.ExitCode() << ",\n"
                   << "  \"interrupted_stage\": \"" << Escape(stage) << "\",\n"
                   << "  \"message\": \"" << Escape(interrupted.what()) << "\"\n}\n";
        WriteText(spec.work_dir / "runs" / nonce / "failure.json",
                  diagnostic.str());
      } catch (const std::exception& failure_error) {
        std::ostringstream diagnostic;
        diagnostic << "{\n  \"schema_version\": 1,\n"
                   << "  \"status\": \"failed\",\n"
                   << "  \"run_id\": \"" << nonce << "\",\n"
                   << "  \"observed_utc\": \"" << UtcTimestamp() << "\",\n"
                   << "  \"failed_stage\": \"" << Escape(stage) << "\",\n"
                   << "  \"message\": \"" << Escape(failure_error.what()) << "\"\n}\n";
        WriteText(spec.work_dir / "runs" / nonce / "failure.json",
                  diagnostic.str());
      }
    } catch (...) {
    }
    for (auto it = published.rbegin(); it != published.rend(); ++it) {
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

}  // namespace ramag

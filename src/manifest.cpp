#include "ramag/manifest.hpp"

#include <chrono>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace ramag {
namespace {

std::string EscapeJson(std::string_view value) {
  std::ostringstream escaped;
  for (const unsigned char character : value) {
    switch (character) {
      case '\"':
        escaped << "\\\"";
        break;
      case '\\':
        escaped << "\\\\";
        break;
      case '\b':
        escaped << "\\b";
        break;
      case '\f':
        escaped << "\\f";
        break;
      case '\n':
        escaped << "\\n";
        break;
      case '\r':
        escaped << "\\r";
        break;
      case '\t':
        escaped << "\\t";
        break;
      default:
        if (character < 0x20U) {
          escaped << "\\u" << std::hex << std::setw(4) << std::setfill('0')
                  << static_cast<unsigned int>(character) << std::dec;
        } else {
          escaped << static_cast<char>(character);
        }
    }
  }
  return escaped.str();
}

std::string Quote(std::string_view value) {
  return "\"" + EscapeJson(value) + "\"";
}

std::string AbsolutePath(const std::filesystem::path& path) {
  if (path.empty()) {
    return {};
  }
  std::error_code error;
  const auto absolute = std::filesystem::absolute(path, error);
  return (error ? path : absolute.lexically_normal()).string();
}

std::string UtcTimestamp() {
  const auto now = std::chrono::system_clock::now();
  const auto time = std::chrono::system_clock::to_time_t(now);
  std::tm utc{};
#if defined(_WIN32)
  gmtime_s(&utc, &time);
#else
  gmtime_r(&time, &utc);
#endif
  std::ostringstream result;
  result << std::put_time(&utc, "%Y-%m-%dT%H:%M:%SZ");
  return result.str();
}

void WriteInput(std::ostringstream& output,
                std::string_view role,
                const FastaData& data,
                bool trailing_comma) {
  output << "    " << Quote(role) << ": {\n"
         << "      \"path\": " << Quote(AbsolutePath(data.source)) << ",\n"
         << "      \"contigs\": " << data.sequences.size() << ",\n"
         << "      \"bases\": " << data.total_bases << ",\n"
         << "      \"ambiguous_bases\": " << data.ambiguous_bases << ",\n"
         << "      \"catalog\": [";
  for (std::size_t index = 0; index < data.sequences.size(); ++index) {
    const auto& sequence = data.sequences[index];
    if (index != 0) {
      output << ',';
    }
    output << "{\"id\":" << sequence.numeric_id << ",\"name\":"
           << Quote(sequence.name) << ",\"length\":" << sequence.size()
           << '}';
  }
  output << "]\n    }" << (trailing_comma ? "," : "") << '\n';
}

}  // namespace

std::string ManifestJson(const ManifestData& data) {
  if (data.reference == nullptr || data.query == nullptr) {
    throw std::invalid_argument("manifest requires reference and query metadata");
  }
  const auto& stats = data.statistics;
  const auto openmp = CurrentOpenMpRuntimeInfo();
  const auto created = data.finished_utc.empty() ? UtcTimestamp() : data.finished_utc;
  const auto started = data.started_utc.empty() ? created : data.started_utc;
  std::ostringstream output;
  output << std::fixed << std::setprecision(9);
  output << "{\n"
         << "  \"schema_version\": 1,\n"
         << "  \"schema\": \"ramag.run-manifest.v1\",\n"
         << "  \"created_utc\": " << Quote(created) << ",\n"
         << "  \"status\": " << Quote(data.status) << ",\n"
         << "  \"exit_code\": " << data.exit_code << ",\n"
         << "  \"run\": {\"id\":" << Quote(data.run_id)
         << ",\"started_at\":" << Quote(started)
         << ",\"finished_at\":" << Quote(created)
         << ",\"exit_code\":" << data.exit_code << "},\n"
         << "  \"binary\": {\n"
         << "    \"path\": " << Quote(AbsolutePath(data.binary_path)) << ",\n"
         << "    \"ramag_version\": " << Quote(RAMAG_VERSION) << ",\n"
         << "    \"ramag_commit\": " << Quote(RAMAG_GIT_COMMIT) << ",\n"
         << "    \"build_type\": " << Quote(RAMAG_BUILD_TYPE) << ",\n"
         << "    \"compiler\": " << Quote(RAMAG_CXX_COMPILER) << ",\n"
         << "    \"compile_flags\": [" << Quote(RAMAG_COMPILE_FLAGS) << "],\n"
         << "    \"openmp_enabled\": "
         << (openmp.enabled ? "true" : "false") << ",\n"
         << "    \"openmp_runtime\": " << Quote(openmp.runtime) << ",\n"
         << "    \"openmp_specification_date\": "
         << openmp.specification_date << "\n"
         << "  },\n"
         << "  \"program\": {\n"
         << "    \"name\": \"RaMA-G\",\n"
         << "    \"cli\": \"ramag\",\n"
         << "    \"version\": " << Quote(RAMAG_VERSION) << ",\n"
         << "    \"commit\": " << Quote(RAMAG_GIT_COMMIT) << ",\n"
         << "    \"binary_path\": " << Quote(AbsolutePath(data.binary_path))
         << ",\n"
         << "    \"invocation\": " << Quote(data.invocation) << "\n"
         << "  },\n"
         << "  \"build\": {\n"
         << "    \"compiler\": " << Quote(RAMAG_CXX_COMPILER) << ",\n"
         << "    \"build_type\": " << Quote(RAMAG_BUILD_TYPE) << ",\n"
         << "    \"flags\": " << Quote(RAMAG_COMPILE_FLAGS) << ",\n"
         << "    \"openmp_enabled\": "
         << (openmp.enabled ? "true" : "false") << ",\n"
         << "    \"openmp_runtime\": " << Quote(openmp.runtime) << ",\n"
         << "    \"openmp_specification_date\": "
         << openmp.specification_date << ",\n"
         << "    \"sufkit_divsufsort_openmp\": "
         << (RAMAG_SUFKIT_DIVSUFSORT_OPENMP ? "true" : "false") << "\n"
         << "  },\n"
         << "  \"dependencies\": {\n"
         << "    \"sufkit_commit\": " << Quote(RAMAG_SUFKIT_COMMIT) << ",\n"
         << "    \"sufkit_source_path\": " << Quote(RAMAG_SUFKIT_SOURCE_PATH)
         << ",\n"
         << "    \"seqpro_commit\": " << Quote(RAMAG_SEQPRO_COMMIT) << ",\n"
         << "    \"seqpro_source_path\": " << Quote(RAMAG_SEQPRO_SOURCE_PATH)
         << ",\n"
         << "    \"dependency_mode\": " << Quote(RAMAG_DEPENDENCY_MODE) << ",\n"
#if RAMAG_USE_PAIRWISE_CORE
         << "    \"ksw2_source\": \"7d08359e0df7f7e6ffcfe67217c3399761cb2129:src/align/ksw2_extz2_sse.c\"\n"
#else
         << "    \"ksw2_source\": \"not-linked; bounded-scalar-baseline\"\n"
#endif
         << "  },\n"
         << "  \"effective_config\": {\n"
         << "    \"reference\": "
         << Quote(AbsolutePath(data.run_spec.reference_path)) << ",\n"
         << "    \"reference_index\": "
         << Quote(AbsolutePath(data.run_spec.reference_index_path)) << ",\n"
         << "    \"index_action\": "
         << Quote(data.run_spec.reference_index_path.empty() ? "built" : "loaded")
         << ",\n"
         << "    \"query\": " << Quote(AbsolutePath(data.run_spec.query_path))
         << ",\n"
         << "    \"output_prefix\": "
         << Quote(AbsolutePath(data.run_spec.output_prefix)) << ",\n"
         << "    \"work_dir\": "
         << Quote(AbsolutePath(data.run_spec.work_dir)) << ",\n"
         << "    \"threads\": " << data.run_spec.threads << ",\n"
         << "    \"openmp_requested_threads\": "
         << data.run_spec.threads << ",\n"
         << "    \"formats\": " << Quote(data.run_spec.formats.ToString())
         << ",\n"
         << "    \"progress\": "
         << Quote(ProgressModeName(data.run_spec.progress.mode)) << ",\n"
         << "    \"progress_interval_seconds\": "
         << data.run_spec.progress.interval_seconds << ",\n"
         << "    \"seed_mode\": "
         << Quote(SeedModeName(data.run_spec.alignment.seed_mode)) << ",\n"
         << "    \"selection_mode\": "
         << Quote(AlignmentSelectionName(data.run_spec.alignment.selection))
         << ",\n"
         << "    \"min_match\": " << data.run_spec.alignment.min_match << ",\n"
         << "    \"smem_min_occurrences\": "
         << data.run_spec.alignment.smem_min_occurrences << ",\n"
         << "    \"max_gap\": " << data.run_spec.alignment.max_gap << ",\n"
         << "    \"diag_diff\": " << data.run_spec.alignment.diag_diff
         << ",\n"
         << "    \"diag_factor\": " << data.run_spec.alignment.diag_factor
         << ",\n"
         << "    \"min_cluster\": " << data.run_spec.alignment.min_cluster
         << ",\n"
         << "    \"break_length\": " << data.run_spec.alignment.break_length
         << ",\n"
         << "    \"max_dp_cells\": "
         << data.run_spec.alignment.max_dp_cells << "\n"
         << "  },\n"
         << "  \"inputs\": {\n";
  WriteInput(output, "reference", *data.reference, true);
  WriteInput(output, "query", *data.query, false);
  output << "  },\n"
         << "  \"routes\": {\n"
         << "    \"input\": " << Quote(data.input_route) << ",\n"
         << "    \"input_parallel\": "
         << Quote(data.input_parallel_route) << ",\n"
         << "    \"index\": " << Quote(data.index_route) << ",\n"
         << "    \"seed\": " << Quote(data.seeding_route) << ",\n"
         << "    \"chain\": " << Quote(data.chaining_route) << ",\n"
#if RAMAG_USE_PAIRWISE_CORE
         << "    \"extension\": [\"pairwise-certified-global-ksw2\",\"pairwise-endpoint-extension\"],\n"
         << "    \"alignment_core\": \"pairwise\",\n"
         << "    \"source_commit\": \"7d08359e0df7f7e6ffcfe67217c3399761cb2129\",\n"
         << "    \"semantic_exceptions\": [\"signed-coordinate-arithmetic\",\"floating-point-dp-best\"],\n"
         << "    \"scoring_contract\": \"scaled-HOXD70;open=40;extend=3;N=mismatch\",\n"
         << "    \"selection_contract\": \"pairwise-dual-dp-intersection;all=pre-selection-records\"\n"
#else
         << "    \"extension\": [\"exact\",\"ungapped\",\"bounded-scalar-dp\"]\n"
#endif
         << "  },\n"
         << "  \"actual_routes\": {\n"
         << "    \"input\": " << Quote(data.input_route) << ",\n"
         << "    \"input_parallel\": "
         << Quote(data.input_parallel_route) << ",\n"
         << "    \"index\": " << Quote(data.index_route) << ",\n"
         << "    \"seed\": " << Quote(data.seeding_route) << ",\n"
         << "    \"chain\": " << Quote(data.chaining_route) << ",\n"
         << "    \"extension\": " << Quote(data.extension_route) << ",\n"
         << "    \"actual_threads\": " << data.actual_threads << "\n"
         << "  },\n"
         << "  \"threading\": {\n"
         << "    \"requested_budget\": " << data.run_spec.threads << ",\n"
         << "    \"openmp_actual_requested_threads\": "
         << openmp.configured_requested_threads << ",\n"
         << "    \"openmp_runtime_max_threads\": "
         << openmp.max_threads << ",\n"
         << "    \"pipeline_worker_threads\": " << data.actual_threads << ",\n"
         << "    \"input_requested_workers\": "
         << data.input_requested_workers << ",\n"
         << "    \"input_actual_workers\": "
         << data.input_actual_workers << ",\n"
         << "    \"input_parallel_route\": "
         << Quote(data.input_parallel_route) << ",\n"
         << "    \"index_requested_threads\": " << data.run_spec.threads << ",\n"
         << "    \"index_observed_threads\": null,\n"
         << "    \"seed_requested_threads\": "
         << stats.seed_requested_threads << ",\n"
         << "    \"seed_scheduled_threads\": "
         << stats.seed_scheduled_threads << ",\n"
         << "    \"seed_worker_threads\": "
         << stats.seed_worker_threads << ",\n"
         << "    \"seed_task_count\": " << stats.seed_task_count << ",\n"
         << "    \"seed_tasks_completed\": "
         << stats.seed_tasks_completed << ",\n"
         << "    \"seed_parallel_route\": "
         << Quote(stats.seed_parallel_route) << ",\n"
         << "    \"chaining_requested_threads\": "
         << stats.chaining_requested_threads << ",\n"
         << "    \"chaining_worker_threads\": "
         << stats.chaining_worker_threads << ",\n"
         << "    \"extension_worker_threads\": "
         << stats.extension_worker_threads << ",\n"
         << "    \"seed_chain_extension_parallel\": "
         << (stats.seed_worker_threads > 1 &&
                     stats.chaining_worker_threads > 1 &&
                     stats.extension_worker_threads > 1
                 ? "true"
                 : "false")
         << "\n"
         << "  },\n"
         << "  \"mam_tiling\": {\n"
         << "    \"worker_cap\": " << stats.seed_mam_worker_cap << ",\n"
         << "    \"tile_bases\": " << stats.seed_mam_tile_bases << ",\n"
         << "    \"oriented_queries\": "
         << stats.seed_oriented_query_count << ",\n"
         << "    \"tile_tasks\": " << stats.seed_mam_tile_task_count
         << ",\n"
         << "    \"tile_tasks_completed\": "
         << stats.seed_mam_tile_tasks_completed << ",\n"
         << "    \"short_query_tasks\": "
         << stats.seed_mam_short_query_task_count << ",\n"
         << "    \"boundary_tasks\": "
         << stats.seed_mam_boundary_task_count << ",\n"
         << "    \"boundary_tasks_completed\": "
         << stats.seed_mam_boundary_tasks_completed << ",\n"
         << "    \"tile_raw_mams\": "
         << stats.seed_mam_tile_raw_count << ",\n"
         << "    \"tile_globally_maximal_mams\": "
         << stats.seed_mam_tile_globally_maximal_count << ",\n"
         << "    \"boundary_raw_mems\": "
         << stats.seed_mam_boundary_mem_raw_count << ",\n"
         << "    \"boundary_patterns\": "
         << stats.seed_mam_boundary_pattern_count << ",\n"
         << "    \"boundary_reference_unique_patterns\": "
         << stats.seed_mam_boundary_reference_unique_count << ",\n"
         << "    \"boundary_recovered_mams\": "
         << stats.seed_mam_boundary_recovered_count << ",\n"
         << "    \"workspace_baseline_bytes\": "
         << stats.seed_mam_workspace_baseline_bytes << ",\n"
         << "    \"workspace_peak_bytes\": "
         << stats.seed_mam_workspace_peak_bytes << ",\n"
         << "    \"resource_limits\": {\n"
         << "      \"boundary_mem_occurrences\": "
         << stats.seed_mam_boundary_mem_occurrence_limit << ",\n"
         << "      \"workspace_bytes\": "
         << stats.seed_mam_workspace_limit_bytes << "\n"
         << "    }\n"
         << "  },\n"
         << "  \"openmp\": {\n"
         << "    \"enabled\": " << (openmp.enabled ? "true" : "false")
         << ",\n"
         << "    \"runtime\": " << Quote(openmp.runtime) << ",\n"
         << "    \"specification_date\": "
         << openmp.specification_date << ",\n"
         << "    \"runtime_max_threads\": " << openmp.max_threads << ",\n"
         << "    \"thread_limit\": " << openmp.thread_limit << ",\n"
         << "    \"dynamic\": " << (openmp.dynamic != 0 ? "true" : "false")
         << ",\n"
         << "    \"max_active_levels\": " << openmp.max_active_levels << ",\n"
         << "    \"run_requested_threads\": " << data.run_spec.threads
         << ",\n"
         << "    \"configured_requested_threads\": "
         << openmp.configured_requested_threads << ",\n"
         << "    \"sufkit_divsufsort_openmp\": "
         << (RAMAG_SUFKIT_DIVSUFSORT_OPENMP ? "true" : "false") << "\n"
         << "  },\n"
         << "  \"adapter_provenance\": {\n";
  std::size_t provenance_index = 0;
  for (const auto& [name, value] : data.adapter_provenance) {
    output << "    " << Quote(name) << ": " << Quote(value);
    if (++provenance_index != data.adapter_provenance.size()) {
      output << ',';
    }
    output << '\n';
  }
  output << "  },\n"
         << "  \"chaining\": {\n"
         << "    \"route\": " << Quote(data.chaining_route) << ",\n"
         << "    \"candidate_pairs\": "
         << stats.chaining_candidate_pairs << ",\n"
         << "    \"legal_edges\": " << stats.chaining_legal_edges << ",\n"
         << "    \"components\": " << stats.chaining_components << ",\n"
         << "    \"max_bucket_occupancy\": "
         << stats.chaining_max_bucket_occupancy << ",\n"
         << "    \"max_component_seeds\": "
         << stats.chaining_max_component_seeds << ",\n"
         << "    \"max_component_edges\": "
         << stats.chaining_max_component_edges << ",\n"
         << "    \"dp_passes\": " << stats.chaining_dp_passes << ",\n"
         << "    \"edge_relaxations\": "
         << stats.chaining_edge_relaxations << ",\n"
         << "    \"requested_threads\": "
         << stats.chaining_requested_threads << ",\n"
         << "    \"worker_threads\": "
         << stats.chaining_worker_threads << ",\n"
         << "    \"extension_worker_threads\": "
         << stats.extension_worker_threads << ",\n"
         << "    \"working_set_peak_bytes\": "
         << stats.chaining_working_set_peak_bytes << ",\n"
         << "    \"resource_limits\": {\n"
         << "      \"working_set_bytes\": "
         << stats.chaining_working_set_limit_bytes << ",\n"
         << "      \"candidate_pairs\": "
         << stats.chaining_candidate_pair_limit << ",\n"
         << "      \"legal_edges\": "
         << stats.chaining_legal_edge_limit << ",\n"
         << "      \"edge_relaxations\": "
         << stats.chaining_edge_relaxation_limit << "\n"
         << "    }\n"
         << "  },\n"
         << "  \"counts\": {\n"
         << "    \"reference_contigs\": " << stats.reference_contigs
         << ",\n"
         << "    \"reference_bases\": " << stats.reference_bases << ",\n"
         << "    \"reference_ambiguous_bases\": "
         << stats.reference_ambiguous_bases << ",\n"
         << "    \"query_contigs\": " << stats.query_contigs << ",\n"
         << "    \"query_bases\": " << stats.query_bases << ",\n"
         << "    \"query_ambiguous_bases\": "
         << stats.query_ambiguous_bases << ",\n"
         << "    \"mem_seeds\": " << stats.mem_seed_count << ",\n"
         << "    \"mam_seeds\": " << stats.mam_seed_count << ",\n"
         << "    \"mum_seeds\": " << stats.mum_seed_count << ",\n"
         << "    \"smem_intervals\": " << stats.smem_interval_count
         << ",\n"
         << "    \"smem_coordinate_seeds\": "
         << stats.smem_coordinate_seed_count << ",\n"
         << "    \"selected_seeds\": " << stats.selected_seed_count << ",\n"
         << "    \"merged_seeds\": " << stats.merged_seed_count << ",\n"
         << "    \"chains\": " << stats.chain_count << ",\n"
         << "    \"candidate_alignments\": "
         << stats.candidate_alignment_count << ",\n"
         << "    \"conflict_rejected_alignments\": "
         << stats.conflict_rejected_alignment_count << ",\n"
         << "    \"alignments\": " << stats.alignment_count << ",\n"
         << "    \"exact_gaps\": " << stats.exact_gap_count << ",\n"
         << "    \"ungapped_extensions\": " << stats.ungapped_gap_count
         << ",\n"
         << "    \"dp_extensions\": " << stats.dp_gap_count << "\n"
         << "  },\n"
         << "  \"stage_wall_seconds\": {\n";
  std::size_t timing_index = 0;
  for (const auto& [name, seconds] : data.stage_wall_seconds) {
    output << "    " << Quote(name) << ": " << seconds;
    if (++timing_index != data.stage_wall_seconds.size()) {
      output << ',';
    }
    output << '\n';
  }
  output << "  },\n"
         << "  \"stages\": [\n";
  timing_index = 0;
  for (const auto& [name, seconds] : data.stage_wall_seconds) {
    output << "    {\"name\":" << Quote(name)
           << ",\"wall_seconds\":" << seconds
           << ",\"cpu_seconds\":null}";
    if (++timing_index != data.stage_wall_seconds.size()) {
      output << ',';
    }
    output << '\n';
  }
  output << "  ],\n"
         << "  \"artifacts\": [\n";
  for (std::size_t index = 0; index < data.artifacts.size(); ++index) {
    const auto& artifact = data.artifacts[index];
    output << "    {\"format\":" << Quote(artifact.format)
           << ",\"path\":" << Quote(AbsolutePath(artifact.path))
           << ",\"bytes\":" << artifact.bytes << ",\"validated\":true"
           << ",\"state\":" << Quote(artifact.state) << '}';
    if (index + 1 != data.artifacts.size()) {
      output << ',';
    }
    output << '\n';
  }
  output << "  ]\n}\n";
  return output.str();
}

}  // namespace ramag

#pragma once

#include "ramag/alignment.hpp"
#include "ramag/model.hpp"

#include <cstdint>
#include <filesystem>
#include <map>
#include <functional>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

namespace ramag {

// The production adapter is written against this exact public sufkit API.
// CMake is responsible for rejecting a source-directory override whose clean
// HEAD differs from this value.
inline constexpr std::string_view kRequiredSufkitCommit =
    "028075e6f2d622fcbdcf76b153bbde0f069ca64f";

struct SufkitIndexOptions {
    // One shared RaMA-G thread budget is reused by sequential stages: first
    // the complete-SA builder, then reference-MAM query/strand tasks.  The two
    // stages never own separate teams at the same time.
    std::uint32_t threads{1};

    // RaMA-G's large-reference policy prefers CaPS once there is enough work
    // to amortize its shared-memory setup.  The 64 MiB threshold is below the
    // 185 Mb Human benchmark where an isolated A/B showed a large wall-time
    // win, while keeping tiny fixtures and small genomes on divsufsort.
    // This selects only the constructor; the complete SA/ISA/LCP query
    // representation and all MAM semantics remain identical.
    std::uint64_t parallel_caps_min_reference_bases{64ULL * 1024ULL * 1024ULL};

    // Reference-MAM queries use non-overlapping query tiles and exact
    // boundary recovery.  The public setting exists primarily so exhaustive
    // tests can force tiny tiles; production keeps the 4 MiB default.
    Length mam_tile_bases{4ULL * 1024ULL * 1024ULL};

    // The seed-query team is capped independently from the requested shared
    // budget.  This bounds simultaneously encoded tile workspaces without
    // changing the suffix-array build budget.
    std::uint32_t mam_worker_cap{16};

    // These are hard failure limits, never permission to truncate or sample.
    // Boundary MEM callbacks are counted before reference-uniqueness
    // filtering.  Workspace covers adapter-owned/query-proportional temporary
    // storage and returned seed vectors, but not the immutable suffix array.
    std::uint64_t boundary_mem_occurrence_limit{100'000'000};
    std::uint64_t mam_workspace_limit_bytes{
        6ULL * 1024ULL * 1024ULL * 1024ULL};
    void (*build_stage_callback)(const char*, void*) = nullptr;
    void* build_stage_context = nullptr;
    std::function<void(std::string_view)> load_stage_callback{};
};

struct SufkitIndexStatistics {
    std::map<std::string, double> load_stage_seconds;
    double sufkit_load_seconds{};
    double reference_validation_seconds{};
    double load_crc_seconds{};
    std::uint64_t load_logical_read_bytes{};
    std::uint32_t requested_threads{1};
    std::uint64_t parallel_caps_min_reference_bases{};
    std::uint32_t sampling_rate{1};
    std::uint8_t coordinate_width{};
    std::uint8_t stored_coordinate_width{};

    std::uint64_t reference_contigs{};
    std::uint64_t reference_bases{};
    std::uint64_t reference_ambiguous_bases{};
    std::uint64_t text_symbols{};
    std::uint64_t suffix_count{};
    std::uint64_t reference_fingerprint{};

    std::uint64_t serialized_bytes{};
    std::uint64_t text_bytes{};
    std::uint64_t suffix_array_bytes{};
    std::uint64_t inverse_suffix_array_bytes{};
    std::uint64_t lcp_bytes{};
    std::uint64_t lcp_primary_bytes{};
    std::uint64_t lcp_overflow_anchors{};
    std::uint64_t lcp_overflow_bytes{};
    std::uint64_t lcp_guide_bytes{};
    std::uint64_t child_bytes{};
    std::uint64_t auxiliary_bytes{};
    std::uint64_t resident_core_bytes{};
    std::uint64_t learned_index_bytes{};
    std::uint64_t prefix_directory_bytes{};
    bool prefix_directory_enabled{};

    double total_build_seconds{};
    double sufkit_build_seconds{};
    double caps_construct_seconds{};
    double caps_output_allocation_seconds{};
    double text_prepare_seconds{};
    double lcp_finalize_seconds{};
    double prefix_directory_seconds{};
    double suffix_array_seconds{};
    double storage_compaction_seconds{};
    double isa_seconds{};
    double lcp_seconds{};
    double child_seconds{};
    double learned_index_seconds{};

    std::string library_version;
    std::string format_version;
    std::string backend;
    std::string backend_signature;
    std::string resource_profile;
    std::string lcp_encoding;
    std::string acceleration;
    std::string lookup_acceleration;
};

struct SufkitSeedStatistics {
    std::vector<MemoryObservation> memory_observations;
    SufkitIndexStatistics index;

    std::uint64_t query_contigs{};
    std::uint64_t query_bases{};
    // Callback-level counts before final deterministic de-duplication.
    std::uint64_t mam_occurrence_count{};
    std::uint64_t mem_occurrence_count{};
    std::uint64_t mum_occurrence_count{};
    std::uint64_t smem_interval_count{};
    std::uint64_t smem_occurrence_count{};
    std::uint64_t mem_covered_by_mam_count{};
    std::uint64_t raw_selected_seed_count{};
    std::uint64_t selected_seed_count{};

    // Reference-MAM enumeration uses stable task slots for every oriented
    // query tile and internal-boundary recovery window.  Other seed modes
    // leave task counts at zero and report one actual query thread.
    std::uint32_t query_requested_threads{1};
    std::uint32_t query_scheduled_threads{1};
    std::uint32_t query_actual_threads{1};
    std::uint64_t query_strand_task_count{};
    std::uint64_t query_strand_tasks_completed{};

    std::uint32_t mam_worker_cap{16};
    std::uint64_t mam_tile_bases{};
    std::uint64_t oriented_query_count{};
    std::uint64_t mam_tile_task_count{};
    std::uint64_t mam_tile_tasks_completed{};
    std::uint64_t mam_short_query_task_count{};
    std::uint64_t mam_boundary_task_count{};
    std::uint64_t mam_boundary_tasks_completed{};
    std::uint64_t mam_tile_raw_count{};
    std::uint64_t mam_tile_globally_maximal_count{};
    std::uint64_t mam_boundary_mem_raw_count{};
    std::uint64_t mam_boundary_pattern_count{};
    std::uint64_t mam_boundary_reference_unique_count{};
    std::uint64_t mam_boundary_recovered_count{};
    std::uint64_t mam_boundary_mem_occurrence_limit{};
    std::uint64_t mam_workspace_baseline_bytes{};
    std::uint64_t mam_workspace_peak_bytes{};
    std::uint64_t mam_workspace_limit_bytes{};

    double enumeration_seconds{};
    std::string actual_route;
    std::string query_parallel_route{"serial"};
};

struct SufkitSeedResult {
    std::vector<Seed> seeds;
    SufkitSeedStatistics statistics;
};

// True only when this translation unit was compiled and linked with the
// pinned sufkit-compatible API.  A build without sufkit remains possible for
// the small exhaustive oracle tests, but cannot use this production adapter.
[[nodiscard]] bool SufkitAdapterAvailable() noexcept;

// Owns the normalized sufkit reference and one immutable, complete suffix
// array.  Forward and reverse-complement query enumeration share this index;
// no reverse-reference index is constructed.
class SufkitSeedIndex {
public:
    [[nodiscard]] static SufkitSeedIndex Build(
        const std::vector<SequenceRecord>& references,
        const SufkitIndexOptions& options = {});
    [[nodiscard]] static SufkitSeedIndex Load(
        const std::filesystem::path& path,
        const std::vector<SequenceRecord>& references,
        const SufkitIndexOptions& options = {});

    void Save(const std::filesystem::path& path) const;

    SufkitSeedIndex(SufkitSeedIndex&&) noexcept;
    SufkitSeedIndex& operator=(SufkitSeedIndex&&) noexcept;
    ~SufkitSeedIndex();

    SufkitSeedIndex(const SufkitSeedIndex&) = delete;
    SufkitSeedIndex& operator=(const SufkitSeedIndex&) = delete;

    // MumReference shares one immutable index across stable oriented tile and
    // boundary tasks.  Tile-local MAMs are globally maximality-filtered;
    // crossing matches are recovered with occurrence-complete boundary MEMs,
    // full-query/reference extension, and full-pattern reference Count().
    // Strict MUM and generalized SMEM are enumerated on complete query
    // records so record-local uniqueness and containment are not weakened by
    // tiling. `references` must be the same ordered records used by Build();
    // it is accepted by const reference so the index does not retain a second
    // copy of the reference text.  Other supported modes retain their kBoth
    // enumeration.
    [[nodiscard]] SufkitSeedResult Enumerate(
        const std::vector<SequenceRecord>& references,
        const std::vector<SequenceRecord>& queries,
        const AlignmentOptions& options) const;

    [[nodiscard]] SufkitIndexStatistics BuildStatistics() const;

private:
    struct Impl;
    explicit SufkitSeedIndex(std::unique_ptr<Impl> implementation) noexcept;
    std::unique_ptr<Impl> implementation_;
};

}  // namespace ramag

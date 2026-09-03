#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace ramag {

// All coordinates inside RaMA-G are zero-based, half-open, and 64 bit.  The
// numeric sequence id is deliberately independent from a FASTA record name.
using Position = std::uint64_t;
using Length = std::uint64_t;
using SequenceId = std::uint32_t;

enum class Strand : std::uint8_t {
    Forward,
    Reverse,
};

enum class SeedMode : std::uint8_t {
    FastHierarchical,
    MumReference,
    Mum,
    Smem,
    MaxMatch,
};

struct SequenceRecord {
    SequenceId numeric_id{};
    std::string name;    // First whitespace-delimited token in the FASTA header.
    std::string header;  // Complete header text, excluding the leading '>'.
    std::string bases;   // Upper-case A/C/G/T/N; N is a seed hard break.

    [[nodiscard]] Length size() const noexcept {
        return static_cast<Length>(bases.size());
    }

    friend bool operator==(const SequenceRecord&, const SequenceRecord&) = default;
};

struct Seed {
    SequenceId reference_id{};
    Position reference_begin{};
    SequenceId query_id{};
    // This is always a coordinate on the original forward query, including
    // when strand == Reverse.
    Position query_begin{};
    Length length{};
    Strand strand{Strand::Forward};

    friend bool operator==(const Seed&, const Seed&) = default;
};

// RaMA-G's canonical internal CIGAR uses extended operations (=, X, I, D).
// Ambiguous M is deliberately excluded so NM/MD can be derived unambiguously.
struct CigarOp {
    char operation{'='};
    Length length{};

    friend bool operator==(const CigarOp&, const CigarOp&) = default;
};

struct AlignmentRecord {
    SequenceId reference_id{};
    Position reference_begin{};
    Position reference_end{};
    SequenceId query_id{};
    // Query coordinates are always on the original forward query and obey
    // query_begin <= query_end.  CIGAR traversal follows `strand`.
    Position query_begin{};
    Position query_end{};
    Strand strand{Strand::Forward};
    std::vector<CigarOp> cigar;
    std::int64_t score{};
    std::uint64_t edit_distance{};
    bool primary{false};

    friend bool operator==(const AlignmentRecord&, const AlignmentRecord&) = default;
};

// Counts and timings that are intrinsic to the alignment core.  CLI/writer
// layers may add their own provenance and output statistics to the manifest.
struct RunStatistics {
    std::uint64_t reference_contigs{};
    std::uint64_t query_contigs{};
    Length reference_bases{};
    Length query_bases{};
    Length reference_ambiguous_bases{};
    Length query_ambiguous_bases{};

    std::uint64_t mem_seed_count{};
    std::uint64_t mam_seed_count{};
    std::uint64_t mum_seed_count{};
    std::uint64_t smem_interval_count{};
    std::uint64_t smem_coordinate_seed_count{};
    std::uint64_t selected_seed_count{};
    std::uint64_t merged_seed_count{};
    std::uint64_t chain_count{};
    std::uint64_t candidate_alignment_count{};
    std::uint64_t conflict_rejected_alignment_count{};
    std::uint64_t alignment_count{};

    std::uint64_t exact_gap_count{};
    std::uint64_t ungapped_gap_count{};
    std::uint64_t dp_gap_count{};

    std::uint32_t seed_requested_threads{1};
    std::uint32_t seed_scheduled_threads{1};
    std::uint32_t seed_worker_threads{1};
    std::uint64_t seed_task_count{};
    std::uint64_t seed_tasks_completed{};
    std::string seed_parallel_route{"serial"};

    std::uint32_t seed_mam_worker_cap{};
    Length seed_mam_tile_bases{};
    std::uint64_t seed_oriented_query_count{};
    std::uint64_t seed_mam_tile_task_count{};
    std::uint64_t seed_mam_tile_tasks_completed{};
    std::uint64_t seed_mam_short_query_task_count{};
    std::uint64_t seed_mam_boundary_task_count{};
    std::uint64_t seed_mam_boundary_tasks_completed{};
    std::uint64_t seed_mam_tile_raw_count{};
    std::uint64_t seed_mam_tile_globally_maximal_count{};
    std::uint64_t seed_mam_boundary_mem_raw_count{};
    std::uint64_t seed_mam_boundary_pattern_count{};
    std::uint64_t seed_mam_boundary_reference_unique_count{};
    std::uint64_t seed_mam_boundary_recovered_count{};
    std::uint64_t seed_mam_boundary_mem_occurrence_limit{};
    std::uint64_t seed_mam_workspace_baseline_bytes{};
    std::uint64_t seed_mam_workspace_peak_bytes{};
    std::uint64_t seed_mam_workspace_limit_bytes{};

    std::string chaining_route;
    std::uint64_t chaining_candidate_pairs{};
    std::uint64_t chaining_legal_edges{};
    std::uint64_t chaining_components{};
    std::uint64_t chaining_max_bucket_occupancy{};
    std::uint64_t chaining_max_component_seeds{};
    std::uint64_t chaining_max_component_edges{};
    std::uint64_t chaining_dp_passes{};
    std::uint64_t chaining_edge_relaxations{};
    std::uint64_t chaining_working_set_peak_bytes{};
    std::uint64_t chaining_working_set_limit_bytes{};
    std::uint64_t chaining_candidate_pair_limit{};
    std::uint64_t chaining_legal_edge_limit{};
    std::uint64_t chaining_edge_relaxation_limit{};
    std::uint32_t chaining_requested_threads{1};
    std::uint32_t chaining_worker_threads{1};
    std::uint32_t extension_worker_threads{1};

    double seed_seconds{};
    double merge_seconds{};
    double chain_and_extension_seconds{};
    double conflict_resolution_seconds{};
    double total_seconds{};
    std::string actual_seed_route;
};

}  // namespace ramag

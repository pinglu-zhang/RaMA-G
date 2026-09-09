#pragma once

#include "ramag/model.hpp"

#include <cstdint>
#include <functional>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace ramag {

class AlignmentError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class DependencyError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class UnsupportedSeedMode : public AlignmentError {
public:
    using AlignmentError::AlignmentError;
};

struct ChainingResourceLimits {
    std::uint64_t working_set_bytes{6ULL * 1024ULL * 1024ULL * 1024ULL};
    std::uint64_t candidate_pairs{500'000'000};
    std::uint64_t legal_edges{100'000'000};
    std::uint64_t edge_relaxations{1'000'000'000};
};

enum class AlignmentSelection : std::uint8_t {
    All,
    ReciprocalOneToOne,
};

struct AlignmentOptions {
    SeedMode seed_mode{SeedMode::FastHierarchical};
    AlignmentSelection selection{AlignmentSelection::All};
    Length min_match{20};
    std::uint64_t smem_min_occurrences{1};
    Length max_gap{90};
    Length diag_diff{5};
    double diag_factor{0.12};
    Length min_cluster{65};
    Length break_length{200};

    // Maximum number of cells in one scalar gap-DP matrix.  A seed edge that
    // would exceed this bound is not chained, so no unbounded matrix is built.
    std::uint64_t max_dp_cells{4'000'000};

    std::int32_t match_score{2};
    std::int32_t mismatch_penalty{4};
    std::int32_t gap_open_penalty{4};
    std::int32_t gap_extend_penalty{2};

    // One shared worker budget is used first for independent chaining groups
    // and then for independent chain-extension tasks.
    std::uint32_t worker_threads{1};

    // Hard failure limits for exact sparse chaining.  Exceeding one of these
    // bounds is never treated as permission to truncate or sample anchors.
    ChainingResourceLimits chaining_limits{};

    // Optional orchestration hooks. They do not participate in alignment
    // decisions and are called only at deterministic phase/task boundaries.
    std::function<void(std::string_view, std::uint64_t, std::uint64_t)>
        progress_callback;
    std::function<void(std::string_view)> interruption_callback;
    // Optional host-supplied RSS sampler. The library does not inspect /proc.
    std::function<std::uint64_t()> resident_bytes_callback;
};

struct AlignmentResult {
    std::vector<AlignmentRecord> alignments;
    RunStatistics statistics;
};

// Apply the current pairwise two-sided DP selection to extended-CIGAR records.
// All mode retains valid records; first retained record per query is primary.
void ResolveAlignmentRecords(
    std::vector<AlignmentRecord>& alignments,
    AlignmentSelection selection = AlignmentSelection::All);

// Enumerate formal occurrence-level MEMs, reference-unique MAMs, strict MUMs,
// or generalized SMEMs. The fast mode returns a MAM skeleton plus MEMs whose
// forward-query interval is not already covered by that skeleton. MUM and
// SMEM uniqueness/containment are scoped to one complete query record.
[[nodiscard]] std::vector<Seed> EnumerateSeeds(
    const std::vector<SequenceRecord>& references,
    const std::vector<SequenceRecord>& queries,
    const AlignmentOptions& options,
    RunStatistics* statistics = nullptr);

// Merge overlapping or adjacent seeds on the same contig pair, strand, and
// exact diagonal.  Reverse seeds are ordered in reverse-complement query space
// but remain represented with original-forward query coordinates.
[[nodiscard]] std::vector<Seed> MergeDiagonalSeeds(
    std::vector<Seed> seeds,
    const std::vector<SequenceRecord>& queries);

// High-level alignment using the same pairwise + KSW2 core as the CLI.
[[nodiscard]] AlignmentResult AlignGenomes(
    const std::vector<SequenceRecord>& references,
    const std::vector<SequenceRecord>& queries,
    const AlignmentOptions& options = {});

// Adapter entry for an external exact-seed provider (for example, a pinned
// sufkit build).  Seeds must already use RaMA-G's contig-local coordinates and
// original-forward query coordinates.  This function deliberately has no
// dependency on the provider's headers and reuses the same merge/chaining/gap
// extension path as AlignGenomes().
[[nodiscard]] AlignmentResult AlignGenomesFromSeeds(
    const std::vector<SequenceRecord>& references,
    const std::vector<SequenceRecord>& queries,
    const AlignmentOptions& options,
    std::vector<Seed> seeds,
    RunStatistics seed_statistics = {});

[[nodiscard]] std::string CigarToString(std::span<const CigarOp> cigar);
[[nodiscard]] Length CigarReferenceLength(std::span<const CigarOp> cigar);
[[nodiscard]] Length CigarQueryLength(std::span<const CigarOp> cigar);
[[nodiscard]] Length CigarMatchLength(std::span<const CigarOp> cigar);
[[nodiscard]] Length CigarEditDistance(std::span<const CigarOp> cigar);



}  // namespace ramag

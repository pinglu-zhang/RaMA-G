#include "ramag/sufkit_adapter.hpp"

#include "ramag/fasta.hpp"
#include "ramag/runtime.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cctype>
#include <cstddef>
#include <exception>
#include <iterator>
#include <limits>
#include <map>
#include <new>
#include <numeric>
#include <string>
#include <string_view>
#include <tuple>
#include <unordered_set>
#include <utility>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef RAMAG_HAVE_SUFKIT
#define RAMAG_HAVE_SUFKIT 0
#endif

#if RAMAG_HAVE_SUFKIT
#include <sufkit/suffix_array.hpp>
#endif

namespace ramag {
namespace {

#if RAMAG_HAVE_SUFKIT

using Clock = std::chrono::steady_clock;
using Interval = std::pair<Position, Position>;

[[nodiscard]] double ElapsedSeconds(Clock::time_point begin) {
    return std::chrono::duration<double>(Clock::now() - begin).count();
}

[[nodiscard]] Position CheckedEnd(Position begin,
                                  Length length,
                                  std::string_view description) {
    if (length > std::numeric_limits<Position>::max() - begin) {
        throw AlignmentError(std::string{description} +
                             " exceeds the 64-bit coordinate range");
    }
    return begin + length;
}

void ValidateRecords(const std::vector<SequenceRecord>& records,
                     std::string_view role) {
    if (records.empty()) {
        throw AlignmentError(std::string{role} + " sequence collection is empty");
    }

    std::unordered_set<SequenceId> ids;
    std::unordered_set<std::string> names;
    for (const SequenceRecord& record : records) {
        if (record.name.empty()) {
            throw AlignmentError(std::string{role} +
                                 " contains an empty sequence name");
        }
        if (record.bases.empty()) {
            throw AlignmentError(std::string{role} + " sequence '" + record.name +
                                 "' is empty");
        }
        if (!ids.insert(record.numeric_id).second) {
            throw AlignmentError(std::string{role} +
                                 " contains duplicate numeric id " +
                                 std::to_string(record.numeric_id));
        }
        if (!names.insert(record.name).second) {
            throw AlignmentError(std::string{role} +
                                 " contains duplicate sequence name '" +
                                 record.name + "'");
        }
        for (const char base : record.bases) {
            if (base != 'A' && base != 'C' && base != 'G' && base != 'T' &&
                base != 'N') {
                throw AlignmentError(std::string{role} + " sequence '" +
                                     record.name +
                                     "' is not normalized to upper-case A/C/G/T/N");
            }
        }
    }
}

[[nodiscard]] std::uint64_t SumBases(
    const std::vector<SequenceRecord>& records,
    std::string_view role) {
    std::uint64_t total = 0;
    for (const SequenceRecord& record : records) {
        if (record.size() > std::numeric_limits<std::uint64_t>::max() - total) {
            throw AlignmentError(std::string{role} +
                                 " base count exceeds the 64-bit range");
        }
        total += record.size();
    }
    return total;
}

[[nodiscard]] bool SeedLess(const Seed& left, const Seed& right) noexcept {
    return std::tie(left.reference_id,
                    left.query_id,
                    left.strand,
                    left.query_begin,
                    left.reference_begin,
                    left.length) <
           std::tie(right.reference_id,
                    right.query_id,
                    right.strand,
                    right.query_begin,
                    right.reference_begin,
                    right.length);
}

void MergeIntervals(std::vector<Interval>& intervals) {
    std::sort(intervals.begin(), intervals.end());
    std::vector<Interval> merged;
    merged.reserve(intervals.size());
    for (const Interval& interval : intervals) {
        if (merged.empty() || interval.first > merged.back().second) {
            merged.push_back(interval);
        } else {
            merged.back().second = std::max(merged.back().second, interval.second);
        }
    }
    intervals = std::move(merged);
}

[[nodiscard]] bool IsCovered(const std::vector<Interval>& intervals,
                             Position begin,
                             Position end) {
    const auto after_begin = std::upper_bound(
        intervals.begin(), intervals.end(), begin,
        [](Position value, const Interval& interval) {
            return value < interval.first;
        });
    if (after_begin == intervals.begin()) {
        return false;
    }
    const Interval& candidate = *std::prev(after_begin);
    return candidate.first <= begin && candidate.second >= end;
}

[[nodiscard]] std::size_t StrandIndex(Strand strand) {
    return strand == Strand::Forward ? 0U : 1U;
}

struct alignas(64) MamTaskWorkspace {
    std::vector<Seed> seeds;
    std::uint64_t mam_occurrence_count{};
    std::uint64_t mem_occurrence_count{};
    std::uint64_t tile_raw_count{};
    std::uint64_t tile_globally_maximal_count{};
    std::uint64_t boundary_mem_raw_count{};
    std::uint64_t boundary_pattern_count{};
    std::uint64_t boundary_reference_unique_count{};
    std::uint64_t boundary_recovered_count{};
    std::exception_ptr failure;
    bool completed{};
};

enum class MamTaskKind : std::uint8_t {
    Tile,
    Boundary,
};

struct OrientedQueryView {
    std::size_t query_index{};
    Strand strand{Strand::Forward};
    std::string_view bases;
};

struct MamTask {
    MamTaskKind kind{MamTaskKind::Tile};
    std::size_t oriented_query_index{};
    Position begin{};
    Position end{};
    Position boundary{};
};

struct ExtendedMatch {
    std::size_t reference_index{};
    Position reference_begin{};
    Position oriented_query_begin{};
    Length length{};
};

// Internal control flow only.  It is deliberately not derived from
// std::exception so cancellation can never be mistaken for, or overwrite,
// the first real task failure retained in a stable task slot.
struct MamTaskCancelled {};

[[nodiscard]] bool CanonicalEqual(char left, char right) noexcept {
    return IsCanonicalBase(left) && left == right;
}

[[nodiscard]] std::uint64_t StableSequenceHash(std::string_view bases) noexcept {
    // FNV-1a is deliberately specified here instead of std::hash so the
    // ordered-reference fingerprint is stable across processes and standard
    // library versions.  IDs, names, and lengths are compared independently.
    std::uint64_t hash = 14695981039346656037ULL;
    for (const unsigned char base : bases) {
        hash ^= base;
        hash *= 1099511628211ULL;
    }
    return hash;
}

[[nodiscard]] std::uint64_t CheckedMultiply(std::uint64_t left,
                                            std::uint64_t right,
                                            std::string_view description) {
    if (left != 0 &&
        right > std::numeric_limits<std::uint64_t>::max() / left) {
        throw AlignmentError(std::string{description} +
                             " exceeds the 64-bit size range");
    }
    return left * right;
}

[[nodiscard]] std::uint64_t CheckedAdd(std::uint64_t left,
                                       std::uint64_t right,
                                       std::string_view description) {
    if (right > std::numeric_limits<std::uint64_t>::max() - left) {
        throw AlignmentError(std::string{description} +
                             " exceeds the 64-bit size range");
    }
    return left + right;
}

class MamWorkspaceBudget {
public:
    MamWorkspaceBudget(std::uint64_t baseline, std::uint64_t limit)
        : current_(baseline), peak_(baseline), limit_(limit) {
        if (baseline > limit) {
            throw AlignmentError(
                "reference-MAM workspace baseline " +
                std::to_string(baseline) + " bytes exceeds limit " +
                std::to_string(limit) + " bytes");
        }
    }

    void Acquire(std::uint64_t bytes, std::string_view purpose) {
        std::uint64_t observed = current_.load(std::memory_order_relaxed);
        for (;;) {
            const std::uint64_t required =
                CheckedAdd(observed, bytes, "reference-MAM workspace");
            if (required > limit_) {
                throw AlignmentError(
                    "reference-MAM workspace for " + std::string{purpose} +
                    " would require " + std::to_string(required) +
                    " bytes, exceeding limit " + std::to_string(limit_) +
                    " bytes");
            }
            if (current_.compare_exchange_weak(
                    observed, required, std::memory_order_relaxed,
                    std::memory_order_relaxed)) {
                std::uint64_t previous_peak =
                    peak_.load(std::memory_order_relaxed);
                while (previous_peak < required &&
                       !peak_.compare_exchange_weak(
                           previous_peak, required, std::memory_order_relaxed,
                           std::memory_order_relaxed)) {
                }
                return;
            }
        }
    }

    void Release(std::uint64_t bytes) noexcept {
        current_.fetch_sub(bytes, std::memory_order_relaxed);
    }

    [[nodiscard]] std::uint64_t Peak() const noexcept {
        return peak_.load(std::memory_order_relaxed);
    }

private:
    std::atomic<std::uint64_t> current_;
    std::atomic<std::uint64_t> peak_;
    std::uint64_t limit_{};
};

void ReserveOneSeed(std::vector<Seed>& seeds,
                    MamWorkspaceBudget& budget) {
    if (seeds.size() != seeds.capacity()) {
        return;
    }
    const std::size_t old_capacity = seeds.capacity();
    const std::size_t new_capacity = old_capacity == 0 ? 8U :
        (old_capacity > std::numeric_limits<std::size_t>::max() / 2U
             ? std::numeric_limits<std::size_t>::max()
             : old_capacity * 2U);
    if (new_capacity == old_capacity) {
        throw AlignmentError("reference-MAM seed vector exceeds the size range");
    }
    const std::uint64_t added = CheckedMultiply(
        static_cast<std::uint64_t>(new_capacity - old_capacity), sizeof(Seed),
        "reference-MAM seed capacity");
    budget.Acquire(added, "task-local seed storage");
    seeds.reserve(new_capacity);
    if (seeds.capacity() > new_capacity) {
        budget.Acquire(
            CheckedMultiply(seeds.capacity() - new_capacity, sizeof(Seed),
                            "reference-MAM allocator capacity growth"),
            "allocator-expanded task-local seed storage");
    }
}

void IncrementBoundaryOccurrence(std::atomic<std::uint64_t>& counter,
                                 std::uint64_t limit) {
    std::uint64_t observed = counter.load(std::memory_order_relaxed);
    for (;;) {
        if (observed >= limit) {
            throw AlignmentError(
                "reference-MAM boundary MEM occurrences exceed limit " +
                std::to_string(limit) + " (observed at least " +
                std::to_string(CheckedAdd(observed, 1,
                                          "boundary MEM occurrence count")) +
                ")");
        }
        if (counter.compare_exchange_weak(
                observed, observed + 1U, std::memory_order_relaxed,
                std::memory_order_relaxed)) {
            return;
        }
    }
}

template <typename Match>
[[nodiscard]] ExtendedMatch ValidateAndMapLocalMatch(
    const Match& match,
    Position window_begin,
    std::string_view window,
    std::string_view oriented_query,
    const std::vector<Length>& reference_lengths,
    std::string_view match_kind) {
    const auto reference_index = static_cast<std::size_t>(match.sequence_id);
    if (reference_index >= reference_lengths.size()) {
        throw DependencyError(
            "sufkit " + std::string{match_kind} +
            " callback returned an unknown reference id " +
            std::to_string(match.sequence_id));
    }
    if (match.strand != sufkit::Strand::kForward) {
        throw DependencyError(
            "sufkit oriented " + std::string{match_kind} +
            " query returned a non-forward strand");
    }
    if (match.length == 0) {
        throw DependencyError("sufkit " + std::string{match_kind} +
                              " callback returned a zero-length match");
    }
    const Position local_end = CheckedEnd(
        match.query_position, match.length, "sufkit local query match");
    if (local_end > window.size()) {
        throw DependencyError("sufkit " + std::string{match_kind} +
                              " callback exceeded its query window");
    }
    const Position oriented_begin = CheckedEnd(
        window_begin, match.query_position, "oriented query match start");
    const Position oriented_end = CheckedEnd(
        oriented_begin, match.length, "oriented query match");
    if (oriented_end > oriented_query.size()) {
        throw DependencyError("sufkit " + std::string{match_kind} +
                              " callback exceeded the oriented query");
    }
    const Position reference_end = CheckedEnd(
        match.reference_position, match.length, "sufkit reference match");
    if (reference_end > reference_lengths[reference_index]) {
        throw DependencyError("sufkit " + std::string{match_kind} +
                              " callback exceeded its reference contig");
    }
    return ExtendedMatch{reference_index, match.reference_position,
                         oriented_begin, match.length};
}

[[nodiscard]] bool IsGloballyTwoSidedMaximal(
    const ExtendedMatch& match,
    std::string_view reference,
    std::string_view oriented_query) {
    if (match.reference_begin > 0 && match.oriented_query_begin > 0 &&
        CanonicalEqual(reference[match.reference_begin - 1],
                       oriented_query[match.oriented_query_begin - 1])) {
        return false;
    }
    const Position reference_end = match.reference_begin + match.length;
    const Position query_end = match.oriented_query_begin + match.length;
    return reference_end >= reference.size() || query_end >= oriented_query.size() ||
           !CanonicalEqual(reference[reference_end], oriented_query[query_end]);
}

[[nodiscard]] ExtendedMatch ExtendGlobally(
    ExtendedMatch match,
    std::string_view reference,
    std::string_view oriented_query) {
    while (match.reference_begin > 0 && match.oriented_query_begin > 0 &&
           CanonicalEqual(reference[match.reference_begin - 1],
                          oriented_query[match.oriented_query_begin - 1])) {
        --match.reference_begin;
        --match.oriented_query_begin;
        ++match.length;
    }
    Position reference_end = match.reference_begin + match.length;
    Position query_end = match.oriented_query_begin + match.length;
    while (reference_end < reference.size() &&
           query_end < oriented_query.size() &&
           CanonicalEqual(reference[reference_end], oriented_query[query_end])) {
        ++reference_end;
        ++query_end;
        ++match.length;
    }
    return match;
}

[[nodiscard]] Seed ToPublicSeed(
    const ExtendedMatch& match,
    const OrientedQueryView& oriented,
    const SequenceRecord& query,
    const std::vector<SequenceId>& reference_ids) {
    const Position oriented_end = CheckedEnd(
        match.oriented_query_begin, match.length, "oriented MAM result");
    const Position query_begin =
        oriented.strand == Strand::Forward
            ? match.oriented_query_begin
            : query.size() - oriented_end;
    return Seed{reference_ids[match.reference_index], match.reference_begin,
                query.numeric_id, query_begin, match.length, oriented.strand};
}

void CheckedAccumulate(std::uint64_t& destination,
                       std::uint64_t value,
                       std::string_view description) {
    if (value > std::numeric_limits<std::uint64_t>::max() - destination) {
        throw AlignmentError(std::string{description} +
                             " exceeds the 64-bit counter range");
    }
    destination += value;
}

[[nodiscard]] const char* RouteName(SeedMode mode) {
    switch (mode) {
        case SeedMode::FastHierarchical:
            return "sufkit-full-sa:mam-skeleton+whole-query-mem-filter";
        case SeedMode::MumReference:
            return "sufkit-full-sa:mumreference";
        case SeedMode::MaxMatch:
            return "sufkit-full-sa:maxmatch";
        case SeedMode::Mum:
            return "sufkit-full-sa:mum+whole-query-record-v1+serial";
        case SeedMode::Smem:
            return "sufkit-full-sa:smem+whole-query-record-v1+serial";
    }
    return "invalid";
}

[[nodiscard]] std::string RecordDescription(const SequenceRecord& record) {
    if (record.header.empty() || record.header == record.name) {
        return {};
    }
    if (record.header.size() > record.name.size() &&
        record.header.compare(0, record.name.size(), record.name) == 0 &&
        std::isspace(
            static_cast<unsigned char>(record.header[record.name.size()])) != 0) {
        std::size_t begin = record.name.size();
        while (begin < record.header.size() &&
               std::isspace(static_cast<unsigned char>(record.header[begin])) != 0) {
            ++begin;
        }
        return record.header.substr(begin);
    }
    return record.header;
}

[[nodiscard]] Strand ConvertStrand(sufkit::Strand strand) {
    switch (strand) {
        case sufkit::Strand::kForward:
            return Strand::Forward;
        case sufkit::Strand::kReverseComplement:
            return Strand::Reverse;
        case sufkit::Strand::kBoth:
            break;
    }
    throw AlignmentError(
        "sufkit maximal-match callback returned non-directional strand 'both'");
}

[[nodiscard]] std::string SufkitFailure(std::string_view operation,
                                        const sufkit::Error& error) {
    return "sufkit " + std::string{operation} + " failed [" +
           sufkit::ToString(error.Code()) + "]: " + error.what();
}

template <typename Match>
[[nodiscard]] Seed ConvertMatch(
    const Match& match,
    const SequenceRecord& query,
    const std::vector<SequenceId>& reference_ids,
    const std::vector<Length>& reference_lengths) {
    const auto reference_index = static_cast<std::size_t>(match.sequence_id);
    if (reference_index >= reference_ids.size()) {
        throw AlignmentError(
            "sufkit maximal-match callback returned an unknown reference id " +
            std::to_string(match.sequence_id));
    }
    if (match.length == 0) {
        throw AlignmentError(
            "sufkit maximal-match callback returned a zero-length seed");
    }
    const Position reference_end = CheckedEnd(
        match.reference_position, match.length, "sufkit reference seed");
    const Position query_end =
        CheckedEnd(match.query_position, match.length, "sufkit query seed");
    if (reference_end > reference_lengths[reference_index]) {
        throw AlignmentError(
            "sufkit maximal-match callback returned a seed outside reference id " +
            std::to_string(reference_ids[reference_index]));
    }
    if (query_end > query.size()) {
        throw AlignmentError(
            "sufkit maximal-match callback returned a seed outside query id " +
            std::to_string(query.numeric_id));
    }

    // sufkit's public MEM/MAM contract already maps reverse-complement hits to
    // the original forward-query coordinate.  Reversing it again here would
    // be an off-by-orientation defect.
    return Seed{reference_ids[reference_index],
                match.reference_position,
                query.numeric_id,
                match.query_position,
                match.length,
                ConvertStrand(match.strand)};
}

SufkitIndexStatistics ConvertIndexStatistics(
    const sufkit::IndexInfo& info,
    const SufkitIndexOptions& options,
    double total_seconds,
    const sufkit::SuffixArrayBuildStatistics* phases = nullptr) {
    SufkitIndexStatistics statistics;
    statistics.requested_threads = options.threads;
    statistics.parallel_caps_min_reference_bases =
        options.parallel_caps_min_reference_bases;
    statistics.sampling_rate = info.sa_sampling_rate;
    statistics.coordinate_width = info.coordinate_width;
    statistics.stored_coordinate_width = info.stored_coordinate_width;
    statistics.reference_contigs = info.sequence_count;
    statistics.reference_bases = info.total_bases;
    statistics.reference_ambiguous_bases = info.ambiguous_bases;
    statistics.text_symbols = info.text_symbols;
    statistics.suffix_count = info.suffix_count;
    statistics.reference_fingerprint = info.fingerprint;
    statistics.serialized_bytes = info.serialized_bytes;
    statistics.text_bytes = info.text_bytes;
    statistics.suffix_array_bytes = info.sa_bytes;
    statistics.inverse_suffix_array_bytes = info.isa_bytes;
    statistics.lcp_bytes = info.lcp_bytes;
    statistics.lcp_primary_bytes = info.lcp_primary_bytes;
    statistics.lcp_overflow_anchors = info.lcp_overflow_anchors;
    statistics.lcp_overflow_bytes = info.lcp_overflow_bytes;
    statistics.lcp_guide_bytes = info.lcp_guide_bytes;
    statistics.child_bytes = info.child_bytes;
    statistics.auxiliary_bytes = info.auxiliary_bytes;
    statistics.resident_core_bytes = info.resident_core_bytes;
    statistics.learned_index_bytes = info.learned_index_bytes;
    const std::uint64_t persistent_auxiliary_bytes =
        info.isa_bytes + info.lcp_bytes + info.child_bytes;
    if (info.auxiliary_bytes < persistent_auxiliary_bytes) {
        throw DependencyError(
            "sufkit index auxiliary-memory accounting is inconsistent");
    }
    statistics.prefix_directory_bytes =
        info.auxiliary_bytes - persistent_auxiliary_bytes;
    statistics.prefix_directory_enabled =
        statistics.prefix_directory_bytes != 0;
    statistics.total_build_seconds = total_seconds;
    if (phases != nullptr) {
        statistics.suffix_array_seconds = phases->sa_seconds;
        statistics.storage_compaction_seconds = phases->storage_compaction_seconds;
        statistics.isa_seconds = phases->isa_seconds;
        statistics.lcp_seconds = phases->lcp_seconds;
        statistics.child_seconds = phases->child_seconds;
        statistics.learned_index_seconds = phases->learned_index_seconds;
    }
    statistics.library_version = info.library_version;
    statistics.format_version = info.format_version;
    statistics.backend = info.backend;
    statistics.backend_signature = info.backend_signature;
    statistics.resource_profile = sufkit::ToString(info.sa_resource_profile);
    statistics.lcp_encoding = sufkit::ToString(info.lcp_encoding);
    statistics.acceleration = sufkit::ToString(info.sa_acceleration);
    statistics.lookup_acceleration =
        sufkit::ToString(info.sa_lookup_acceleration);
    return statistics;
}

#else

[[noreturn]] void ThrowUnavailable() {
    throw DependencyError(
        "the production seed adapter is unavailable because RaMA-G was built "
        "without the pinned sufkit dependency");
}

#endif

}  // namespace

#if RAMAG_HAVE_SUFKIT

struct SufkitSeedIndex::Impl {
    std::vector<SequenceId> reference_ids;
    std::vector<std::string> reference_names;
    std::vector<Length> reference_lengths;
    std::vector<std::uint64_t> reference_content_hashes;
    sufkit::SuffixArray index;
    SufkitIndexStatistics statistics;
    SufkitIndexOptions options;

    Impl(std::vector<SequenceId> ids,
         std::vector<std::string> names,
         std::vector<Length> lengths,
         std::vector<std::uint64_t> content_hashes,
         sufkit::SuffixArray suffix_array,
         SufkitIndexStatistics build_statistics,
         SufkitIndexOptions build_options)
        : reference_ids(std::move(ids)),
          reference_names(std::move(names)),
          reference_lengths(std::move(lengths)),
          reference_content_hashes(std::move(content_hashes)),
          index(std::move(suffix_array)),
          statistics(std::move(build_statistics)),
          options(std::move(build_options)) {}
};

#else

struct SufkitSeedIndex::Impl {};

#endif

bool SufkitAdapterAvailable() noexcept {
#if RAMAG_HAVE_SUFKIT
    return true;
#else
    return false;
#endif
}

SufkitSeedIndex::SufkitSeedIndex(
    std::unique_ptr<Impl> implementation) noexcept
    : implementation_(std::move(implementation)) {}

SufkitSeedIndex::SufkitSeedIndex(SufkitSeedIndex&&) noexcept = default;
SufkitSeedIndex& SufkitSeedIndex::operator=(SufkitSeedIndex&&) noexcept = default;
SufkitSeedIndex::~SufkitSeedIndex() = default;

SufkitSeedIndex SufkitSeedIndex::Build(
    const std::vector<SequenceRecord>& references,
    const SufkitIndexOptions& options) {
#if RAMAG_HAVE_SUFKIT
    ValidateRecords(references, "reference");
    if (options.threads == 0) {
        throw AlignmentError("sufkit index thread count must be greater than zero");
    }
    if (options.mam_tile_bases == 0) {
        throw AlignmentError("reference-MAM tile size must be greater than zero");
    }
    if (options.mam_worker_cap == 0) {
        throw AlignmentError("reference-MAM worker cap must be greater than zero");
    }
    if (options.boundary_mem_occurrence_limit == 0) {
        throw AlignmentError(
            "reference-MAM boundary MEM occurrence limit must be greater than zero");
    }
    if (options.mam_workspace_limit_bytes == 0) {
        throw AlignmentError(
            "reference-MAM workspace limit must be greater than zero");
    }

    try {
        const auto begin = Clock::now();
        std::vector<sufkit::SequenceRecord> sufkit_records;
        std::vector<SequenceId> reference_ids;
        std::vector<std::string> reference_names;
        std::vector<Length> reference_lengths;
        std::vector<std::uint64_t> reference_content_hashes;
        sufkit_records.reserve(references.size());
        reference_ids.reserve(references.size());
        reference_names.reserve(references.size());
        reference_lengths.reserve(references.size());
        reference_content_hashes.reserve(references.size());
        for (const SequenceRecord& record : references) {
            sufkit_records.push_back(sufkit::SequenceRecord{
                record.name, RecordDescription(record), record.bases});
            reference_ids.push_back(record.numeric_id);
            reference_names.push_back(record.name);
            reference_lengths.push_back(record.size());
            reference_content_hashes.push_back(
                StableSequenceHash(record.bases));
        }

        auto genome =
            sufkit::GenomeReference::FromRecords(std::move(sufkit_records));
        sufkit::SuffixArrayBuildStatistics phase_statistics;
        auto build_options = sufkit::FastSuffixArrayBuildOptions();
        build_options.threads = options.threads;
        const std::uint64_t reference_bases =
            SumBases(references, "reference");
        build_options.backend =
            options.threads > 1 &&
                    reference_bases >=
                        options.parallel_caps_min_reference_bases
                ? sufkit::SaBackend::kCaps
                : sufkit::SaBackend::kDivsufsort;
        build_options.sampling_rate = 1;
        build_options.learned_index.enabled = false;
        build_options.statistics = &phase_statistics;

        auto index = sufkit::SuffixArray::Build(genome, build_options);
        const sufkit::IndexInfo info = index.GetInfo();
        if (index.SamplingRate() != 1 || info.sa_sampling_rate != 1) {
            throw DependencyError(
                "sufkit violated the RaMA-G complete suffix-array contract: "
                "sampling_rate is not 1");
        }
        if (info.sequence_count != references.size()) {
            throw DependencyError(
                "sufkit index sequence count differs from the supplied reference");
        }

        SufkitIndexStatistics statistics;
        statistics.requested_threads = options.threads;
        statistics.parallel_caps_min_reference_bases =
            options.parallel_caps_min_reference_bases;
        statistics.sampling_rate = info.sa_sampling_rate;
        statistics.coordinate_width = info.coordinate_width;
        statistics.stored_coordinate_width = info.stored_coordinate_width;
        statistics.reference_contigs = info.sequence_count;
        statistics.reference_bases = info.total_bases;
        statistics.reference_ambiguous_bases = info.ambiguous_bases;
        statistics.text_symbols = info.text_symbols;
        statistics.suffix_count = info.suffix_count;
        statistics.reference_fingerprint = info.fingerprint;
        statistics.serialized_bytes = info.serialized_bytes;
        statistics.text_bytes = info.text_bytes;
        statistics.suffix_array_bytes = info.sa_bytes;
        statistics.inverse_suffix_array_bytes = info.isa_bytes;
        statistics.lcp_bytes = info.lcp_bytes;
        statistics.lcp_primary_bytes = info.lcp_primary_bytes;
        statistics.lcp_overflow_anchors = info.lcp_overflow_anchors;
        statistics.lcp_overflow_bytes = info.lcp_overflow_bytes;
        statistics.lcp_guide_bytes = info.lcp_guide_bytes;
        statistics.child_bytes = info.child_bytes;
        statistics.auxiliary_bytes = info.auxiliary_bytes;
        statistics.resident_core_bytes = info.resident_core_bytes;
        statistics.learned_index_bytes = info.learned_index_bytes;
        const std::uint64_t persistent_auxiliary_bytes =
            info.isa_bytes + info.lcp_bytes + info.child_bytes;
        if (info.auxiliary_bytes < persistent_auxiliary_bytes) {
            throw DependencyError(
                "sufkit index auxiliary-memory accounting is inconsistent");
        }
        statistics.prefix_directory_bytes =
            info.auxiliary_bytes - persistent_auxiliary_bytes;
        statistics.prefix_directory_enabled =
            statistics.prefix_directory_bytes != 0;
        statistics.total_build_seconds = ElapsedSeconds(begin);
        statistics.suffix_array_seconds = phase_statistics.sa_seconds;
        statistics.storage_compaction_seconds =
            phase_statistics.storage_compaction_seconds;
        statistics.isa_seconds = phase_statistics.isa_seconds;
        statistics.lcp_seconds = phase_statistics.lcp_seconds;
        statistics.child_seconds = phase_statistics.child_seconds;
        statistics.learned_index_seconds =
            phase_statistics.learned_index_seconds;
        statistics.library_version = info.library_version;
        statistics.format_version = info.format_version;
        statistics.backend = info.backend;
        statistics.backend_signature = info.backend_signature;
        statistics.resource_profile =
            sufkit::ToString(info.sa_resource_profile);
        statistics.lcp_encoding = sufkit::ToString(info.lcp_encoding);
        statistics.acceleration = sufkit::ToString(info.sa_acceleration);
        statistics.lookup_acceleration =
            sufkit::ToString(info.sa_lookup_acceleration);

        return SufkitSeedIndex(std::make_unique<Impl>(
            std::move(reference_ids), std::move(reference_names),
            std::move(reference_lengths),
            std::move(reference_content_hashes), std::move(index),
            std::move(statistics), options));
    } catch (const InterruptedError&) {
        throw;
    } catch (const DependencyError&) {
        throw;
    } catch (const AlignmentError&) {
        throw;
    } catch (const std::bad_alloc&) {
        throw;
    } catch (const sufkit::Error& error) {
        throw DependencyError(SufkitFailure("complete suffix-array build", error));
    } catch (const std::exception& error) {
        throw DependencyError(
            "sufkit complete suffix-array build failed: " +
            std::string{error.what()});
    }
#else
    static_cast<void>(references);
    static_cast<void>(options);
    ThrowUnavailable();
#endif
}

SufkitSeedIndex SufkitSeedIndex::Load(
    const std::filesystem::path& path,
    const std::vector<SequenceRecord>& references,
    const SufkitIndexOptions& options) {
#if RAMAG_HAVE_SUFKIT
    ValidateRecords(references, "reference");
    if (options.threads == 0 || options.mam_tile_bases == 0 ||
        options.mam_worker_cap == 0 ||
        options.boundary_mem_occurrence_limit == 0 ||
        options.mam_workspace_limit_bytes == 0) {
        throw AlignmentError("invalid zero-valued sufkit index option");
    }
    try {
        const auto begin = Clock::now();
        auto index = sufkit::SuffixArray::Load(path);
        const auto info = index.GetInfo();
        if (info.kind != sufkit::IndexKind::kSuffixArray ||
            index.SamplingRate() != 1 || info.sa_sampling_rate != 1) {
            throw DependencyError(
                "reference index is not a complete standalone suffix array");
        }
        if (info.sa_resource_profile != sufkit::SaResourceProfile::kFast ||
            info.sa_acceleration != sufkit::SaAcceleration::kLcpSuffixLink) {
            throw DependencyError(
                "reference index lacks the RaMA-G Fast SA+ISA+LCP capability");
        }
        if (info.library_version != "0.3.0") {
            throw DependencyError("reference index was not written by sufkit 0.3.0");
        }
        if (info.sequence_count != references.size()) {
            throw DependencyError(
                "reference index contig count differs from --reference");
        }

        std::vector<sufkit::SequenceRecord> sufkit_records;
        std::vector<SequenceId> reference_ids;
        std::vector<std::string> reference_names;
        std::vector<Length> reference_lengths;
        std::vector<std::uint64_t> reference_content_hashes;
        sufkit_records.reserve(references.size());
        reference_ids.reserve(references.size());
        reference_names.reserve(references.size());
        reference_lengths.reserve(references.size());
        reference_content_hashes.reserve(references.size());
        for (std::size_t ordinal = 0; ordinal < references.size(); ++ordinal) {
            const auto& record = references[ordinal];
            const auto stored = index.GetSequenceInfo(
                static_cast<sufkit::SequenceId>(ordinal));
            const auto description = RecordDescription(record);
            const auto ambiguous = static_cast<std::uint64_t>(
                std::count(record.bases.begin(), record.bases.end(), 'N'));
            if (stored.id != static_cast<sufkit::SequenceId>(ordinal) ||
                stored.name != record.name ||
                stored.description != description ||
                stored.length != record.size() ||
                stored.ambiguous_bases != ambiguous) {
                throw DependencyError(
                    "reference index metadata differs from --reference at contig " +
                    std::to_string(ordinal) + " ('" + record.name + "')");
            }
            sufkit_records.push_back(
                sufkit::SequenceRecord{record.name, description, record.bases});
            reference_ids.push_back(record.numeric_id);
            reference_names.push_back(record.name);
            reference_lengths.push_back(record.size());
            reference_content_hashes.push_back(StableSequenceHash(record.bases));
        }
        auto normalized_reference =
            sufkit::GenomeReference::FromRecords(std::move(sufkit_records));
        if (normalized_reference.Fingerprint() != info.fingerprint) {
            throw DependencyError(
                "reference index normalized-content fingerprint differs from --reference");
        }
        auto statistics = ConvertIndexStatistics(
            info, options, ElapsedSeconds(begin));
        return SufkitSeedIndex(std::make_unique<Impl>(
            std::move(reference_ids), std::move(reference_names),
            std::move(reference_lengths), std::move(reference_content_hashes),
            std::move(index), std::move(statistics), options));
    } catch (const DependencyError&) {
        throw;
    } catch (const AlignmentError&) {
        throw;
    } catch (const std::bad_alloc&) {
        throw;
    } catch (const sufkit::Error& error) {
        throw DependencyError(SufkitFailure("reference index load", error));
    } catch (const std::exception& error) {
        throw DependencyError("sufkit reference index load failed: " +
                              std::string{error.what()});
    }
#else
    static_cast<void>(path);
    static_cast<void>(references);
    static_cast<void>(options);
    ThrowUnavailable();
#endif
}

void SufkitSeedIndex::Save(const std::filesystem::path& path) const {
#if RAMAG_HAVE_SUFKIT
    if (!implementation_) {
        throw AlignmentError("cannot save a moved-from sufkit seed index");
    }
    try {
        implementation_->index.Save(path);
    } catch (const sufkit::Error& error) {
        throw DependencyError(SufkitFailure("reference index save", error));
    } catch (const std::exception& error) {
        throw DependencyError("sufkit reference index save failed: " +
                              std::string{error.what()});
    }
#else
    static_cast<void>(path);
    ThrowUnavailable();
#endif
}

SufkitSeedResult SufkitSeedIndex::Enumerate(
    const std::vector<SequenceRecord>& references,
    const std::vector<SequenceRecord>& queries,
    const AlignmentOptions& options) const {
    if (options.min_match == 0) {
        throw AlignmentError("sufkit seed min_match must be greater than zero");
    }
    if (options.smem_min_occurrences == 0) {
        throw AlignmentError(
            "sufkit SMEM min_occurrences must be greater than zero");
    }
#if RAMAG_HAVE_SUFKIT
    if (!implementation_) {
        throw AlignmentError("cannot query a moved-from sufkit seed index");
    }
    ValidateRecords(references, "reference");
    ValidateRecords(queries, "query");
    if (references.size() != implementation_->reference_ids.size()) {
        throw AlignmentError(
            "reference collection differs from the ordered records used to "
            "build the sufkit index: contig count changed");
    }
    for (std::size_t index = 0; index < references.size(); ++index) {
        const SequenceRecord& reference = references[index];
        if (reference.numeric_id != implementation_->reference_ids[index] ||
            reference.name != implementation_->reference_names[index] ||
            reference.size() != implementation_->reference_lengths[index] ||
            StableSequenceHash(reference.bases) !=
                implementation_->reference_content_hashes[index]) {
            throw AlignmentError(
                "reference collection differs from the ordered records used "
                "to build the sufkit index at contig " +
                std::to_string(index) + " ('" + reference.name + "')");
        }
    }

    try {
        const auto enumeration_begin = Clock::now();
        SufkitSeedResult result;
        result.statistics.index = implementation_->statistics;
        result.statistics.query_contigs = queries.size();
        result.statistics.query_bases = SumBases(queries, "query");
        result.statistics.actual_route = RouteName(options.seed_mode);
        result.statistics.query_requested_threads =
            implementation_->statistics.requested_threads;
        result.statistics.mam_worker_cap = implementation_->options.mam_worker_cap;
        result.statistics.mam_tile_bases =
            implementation_->options.mam_tile_bases;
        result.statistics.mam_boundary_mem_occurrence_limit =
            implementation_->options.boundary_mem_occurrence_limit;
        result.statistics.mam_workspace_limit_bytes =
            implementation_->options.mam_workspace_limit_bytes;

        sufkit::MamOptions mam_options;
        mam_options.min_length = options.min_match;
        mam_options.strands = sufkit::StrandMode::kBoth;
        mam_options.algorithm = sufkit::MemSearchAlgorithm::kAutoSelect;
        mam_options.lookup_algorithm = sufkit::SaSearchAlgorithm::kAutoSelect;

        sufkit::MemOptions mem_options;
        mem_options.min_length = options.min_match;
        mem_options.strands = sufkit::StrandMode::kBoth;
        mem_options.algorithm = sufkit::MemSearchAlgorithm::kAutoSelect;
        mem_options.lookup_algorithm = sufkit::SaSearchAlgorithm::kAutoSelect;
        // An absent skip multiplier asks sufkit for its complete, correctness-
        // preserving algorithmic default.  It is not an occurrence cap.
        mem_options.skip_multiplier.reset();

        sufkit::MumOptions mum_options;
        mum_options.min_length = options.min_match;
        mum_options.strands = sufkit::StrandMode::kBoth;
        mum_options.algorithm = sufkit::MemSearchAlgorithm::kAutoSelect;
        mum_options.lookup_algorithm =
            sufkit::SaSearchAlgorithm::kAutoSelect;

        sufkit::SmemOptions smem_options;
        smem_options.min_length = options.min_match;
        smem_options.min_occurrences = options.smem_min_occurrences;
        smem_options.strands = sufkit::StrandMode::kBoth;
        smem_options.algorithm = sufkit::MemSearchAlgorithm::kAutoSelect;
        smem_options.lookup_algorithm =
            sufkit::SaSearchAlgorithm::kAutoSelect;

        if (options.seed_mode == SeedMode::MumReference) {
            if (queries.size() >
                std::numeric_limits<std::size_t>::max() / 2U) {
                throw AlignmentError(
                    "oriented query count exceeds the size range");
            }

            std::vector<std::string> reverse_queries = ReverseComplements(
                queries, implementation_->statistics.requested_threads);

            std::vector<OrientedQueryView> oriented_queries;
            oriented_queries.reserve(queries.size() * 2U);
            for (std::size_t query_index = 0; query_index < queries.size();
                 ++query_index) {
                oriented_queries.push_back(
                    OrientedQueryView{query_index, Strand::Forward,
                                      queries[query_index].bases});
                oriented_queries.push_back(
                    OrientedQueryView{query_index, Strand::Reverse,
                                      reverse_queries[query_index]});
            }
            result.statistics.oriented_query_count = oriented_queries.size();

            const Length tile_bases = implementation_->options.mam_tile_bases;
            std::uint64_t tile_task_count = 0;
            std::uint64_t boundary_task_count = 0;
            std::uint64_t short_query_count = 0;
            for (const OrientedQueryView& oriented : oriented_queries) {
                const std::uint64_t length = oriented.bases.size();
                const std::uint64_t tiles =
                    length / tile_bases + (length % tile_bases != 0 ? 1U : 0U);
                tile_task_count = CheckedAdd(
                    tile_task_count, tiles, "reference-MAM tile task count");
                boundary_task_count = CheckedAdd(
                    boundary_task_count, tiles > 0 ? tiles - 1U : 0U,
                    "reference-MAM boundary task count");
                if (tiles == 1U) {
                    ++short_query_count;
                }
            }
            const std::uint64_t task_count_u64 = CheckedAdd(
                tile_task_count, boundary_task_count,
                "reference-MAM total task count");
            if (task_count_u64 >
                static_cast<std::uint64_t>(
                    std::numeric_limits<std::size_t>::max())) {
                throw AlignmentError(
                    "reference-MAM task count exceeds the size range");
            }
            const std::size_t task_count =
                static_cast<std::size_t>(task_count_u64);

            std::vector<MamTask> task_descriptors;
            task_descriptors.reserve(task_count);
            const Length recovery_flank =
                std::max<Length>(1U, options.min_match - 1U);
            for (std::size_t oriented_index = 0;
                 oriented_index < oriented_queries.size(); ++oriented_index) {
                const Length length = oriented_queries[oriented_index].bases.size();
                for (Position begin = 0; begin < length;) {
                    const Length remaining = length - begin;
                    const Length width = std::min(tile_bases, remaining);
                    const Position end = CheckedEnd(
                        begin, width, "reference-MAM tile interval");
                    task_descriptors.push_back(MamTask{
                        MamTaskKind::Tile, oriented_index, begin, end, 0});
                    begin = end;
                }
                for (Position boundary = tile_bases; boundary < length;) {
                    const Length left_flank =
                        std::min<Length>(recovery_flank, boundary);
                    const Length right_flank =
                        std::min<Length>(recovery_flank, length - boundary);
                    task_descriptors.push_back(MamTask{
                        MamTaskKind::Boundary, oriented_index,
                        boundary - left_flank,
                        CheckedEnd(boundary, right_flank,
                                   "reference-MAM boundary window"),
                        boundary});
                    if (tile_bases > length - boundary) {
                        break;
                    }
                    boundary += tile_bases;
                }
            }
            if (task_descriptors.size() != task_count) {
                throw DependencyError(
                    "reference-MAM stable task construction count mismatch");
            }

            std::vector<MamTaskWorkspace> tasks(task_count);
            std::uint64_t workspace_baseline = 0;
            workspace_baseline = CheckedAdd(
                workspace_baseline,
                CheckedMultiply(reverse_queries.capacity(), sizeof(std::string),
                                "reference-MAM reverse query vector"),
                "reference-MAM workspace baseline");
            for (const std::string& reverse : reverse_queries) {
                workspace_baseline = CheckedAdd(
                    workspace_baseline, reverse.capacity(),
                    "reference-MAM reverse query storage");
            }
            workspace_baseline = CheckedAdd(
                workspace_baseline,
                CheckedMultiply(oriented_queries.capacity(),
                                sizeof(OrientedQueryView),
                                "reference-MAM oriented query views"),
                "reference-MAM workspace baseline");
            workspace_baseline = CheckedAdd(
                workspace_baseline,
                CheckedMultiply(task_descriptors.capacity(), sizeof(MamTask),
                                "reference-MAM task descriptors"),
                "reference-MAM workspace baseline");
            workspace_baseline = CheckedAdd(
                workspace_baseline,
                CheckedMultiply(tasks.capacity(), sizeof(MamTaskWorkspace),
                                "reference-MAM task workspaces"),
                "reference-MAM workspace baseline");
            MamWorkspaceBudget workspace_budget(
                workspace_baseline,
                implementation_->options.mam_workspace_limit_bytes);

            const std::size_t scheduled_size = std::min<std::size_t>(
                {task_count,
                 implementation_->statistics.requested_threads,
                 implementation_->options.mam_worker_cap,
                 static_cast<std::size_t>(
                     std::numeric_limits<int>::max())});
            const int scheduled_threads =
                static_cast<int>(std::max<std::size_t>(scheduled_size, 1U));
            result.statistics.query_scheduled_threads =
                static_cast<std::uint32_t>(scheduled_threads);
            result.statistics.query_strand_task_count = task_count;
            result.statistics.mam_tile_task_count = tile_task_count;
            result.statistics.mam_short_query_task_count = short_query_count;
            result.statistics.mam_boundary_task_count = boundary_task_count;
            result.statistics.mam_workspace_baseline_bytes = workspace_baseline;

            std::atomic<std::uint64_t> boundary_occurrences{0};
            std::atomic<bool> cancellation_requested{false};

            const auto run_task = [&](std::size_t task_id) {
                MamTaskWorkspace& task = tasks[task_id];
                try {
                    if (cancellation_requested.load(std::memory_order_relaxed)) {
                        return;
                    }
                    const MamTask& descriptor = task_descriptors[task_id];
                    const OrientedQueryView& oriented =
                        oriented_queries[descriptor.oriented_query_index];
                    const SequenceRecord& query = queries[oriented.query_index];
                    const std::string_view window = oriented.bases.substr(
                        static_cast<std::size_t>(descriptor.begin),
                        static_cast<std::size_t>(descriptor.end -
                                                 descriptor.begin));

                    if (descriptor.kind == MamTaskKind::Tile) {
                        sufkit::MamOptions task_options = mam_options;
                        task_options.strands = sufkit::StrandMode::kForward;
                        implementation_->index.ForEachMam(
                            window, task_options,
                            [&](const sufkit::MamMatch& match) {
                                if (options.interruption_callback) {
                                    options.interruption_callback("seed-enumeration");
                                }
                                if (cancellation_requested.load(
                                        std::memory_order_relaxed)) {
                                    throw MamTaskCancelled{};
                                }
                                ++task.tile_raw_count;
                                if (match.length < options.min_match) {
                                    throw DependencyError(
                                        "sufkit tile MAM callback returned a "
                                        "match shorter than min_match");
                                }
                                const ExtendedMatch global =
                                    ValidateAndMapLocalMatch(
                                        match, descriptor.begin, window,
                                        oriented.bases,
                                        implementation_->reference_lengths,
                                        "tile MAM");
                                const std::string_view reference =
                                    references[global.reference_index].bases;
                                if (!IsGloballyTwoSidedMaximal(
                                        global, reference, oriented.bases)) {
                                    return;
                                }
                                ++task.tile_globally_maximal_count;
                                ++task.mam_occurrence_count;
                                ReserveOneSeed(task.seeds, workspace_budget);
                                task.seeds.push_back(ToPublicSeed(
                                    global, oriented, query,
                                    implementation_->reference_ids));
                            });
                    } else {
                        sufkit::MemOptions task_options = mem_options;
                        task_options.strands = sufkit::StrandMode::kForward;
                        implementation_->index.ForEachMem(
                            window, task_options,
                            [&](const sufkit::MemMatch& match) {
                                if (options.interruption_callback) {
                                    options.interruption_callback("seed-enumeration");
                                }
                                if (cancellation_requested.load(
                                        std::memory_order_relaxed)) {
                                    throw MamTaskCancelled{};
                                }
                                IncrementBoundaryOccurrence(
                                    boundary_occurrences,
                                    implementation_->options
                                        .boundary_mem_occurrence_limit);
                                ++task.boundary_mem_raw_count;
                                ++task.mem_occurrence_count;
                                if (match.length < options.min_match) {
                                    throw DependencyError(
                                        "sufkit boundary MEM callback returned "
                                        "a match shorter than min_match");
                                }
                                ExtendedMatch global = ValidateAndMapLocalMatch(
                                    match, descriptor.begin, window,
                                    oriented.bases,
                                    implementation_->reference_lengths,
                                    "boundary MEM");
                                global = ExtendGlobally(
                                    global,
                                    references[global.reference_index].bases,
                                    oriented.bases);
                                const Position query_end = CheckedEnd(
                                    global.oriented_query_begin, global.length,
                                    "globally extended boundary MEM");
                                if (global.oriented_query_begin >=
                                        descriptor.boundary ||
                                    query_end <= descriptor.boundary) {
                                    return;
                                }
                                ++task.boundary_pattern_count;
                                const std::string_view pattern =
                                    oriented.bases.substr(
                                        static_cast<std::size_t>(
                                            global.oriented_query_begin),
                                        static_cast<std::size_t>(global.length));
                                if (implementation_->index.Count(
                                        pattern, sufkit::StrandMode::kForward) !=
                                    1U) {
                                    return;
                                }
                                ++task.boundary_reference_unique_count;
                                ++task.boundary_recovered_count;
                                ++task.mam_occurrence_count;
                                ReserveOneSeed(task.seeds, workspace_budget);
                                task.seeds.push_back(ToPublicSeed(
                                    global, oriented, query,
                                    implementation_->reference_ids));
                            });
                    }
                    task.completed = true;
                } catch (const MamTaskCancelled&) {
                    // A different stable task owns the real failure.  Leave
                    // this slot incomplete and do not manufacture a competing
                    // exception that could obscure the root cause.
                } catch (...) {
                    task.failure = std::current_exception();
                    cancellation_requested.store(true,
                                                 std::memory_order_relaxed);
                }
            };

            std::uint32_t actual_threads = 1;
#ifdef _OPENMP
            const bool enter_parallel_region =
                scheduled_threads > 1 && omp_in_parallel() == 0;
#pragma omp parallel if(enter_parallel_region) num_threads(scheduled_threads)
            {
#pragma omp single
                {
                    actual_threads =
                        static_cast<std::uint32_t>(omp_get_num_threads());
                }
#pragma omp for schedule(dynamic, 1)
                for (std::size_t task_id = 0; task_id < task_count; ++task_id) {
                    run_task(task_id);
                }
            }
#else
            for (std::size_t task_id = 0; task_id < task_count; ++task_id) {
                run_task(task_id);
            }
#endif
            result.statistics.query_actual_threads = actual_threads;
            result.statistics.query_parallel_route =
                actual_threads > 1
                    ? "openmp-dynamic-stable-tile-boundary-tasks"
                    : "serial-stable-tile-boundary-tasks";
            result.statistics.actual_route =
                "sufkit-full-sa:mumreference+tiled-mam-boundary-mem-v1+" +
                result.statistics.query_parallel_route;

            std::exception_ptr stable_failure;
            for (const MamTaskWorkspace& task : tasks) {
                if (task.failure) {
                    stable_failure = task.failure;
                    break;
                }
            }
            if (stable_failure) {
                std::rethrow_exception(stable_failure);
            }
            if (cancellation_requested.load(std::memory_order_relaxed)) {
                throw DependencyError(
                    "reference-MAM task cancellation was requested without a "
                    "retained root-cause exception");
            }

            std::uint64_t result_seed_count = 0;
            for (std::size_t task_id = 0; task_id < task_count; ++task_id) {
                MamTaskWorkspace& task = tasks[task_id];
                if (!task.completed) {
                    throw DependencyError(
                        "reference-MAM query/strand task did not complete");
                }
                ++result.statistics.query_strand_tasks_completed;
                if (task_descriptors[task_id].kind == MamTaskKind::Tile) {
                    ++result.statistics.mam_tile_tasks_completed;
                } else {
                    ++result.statistics.mam_boundary_tasks_completed;
                }
                CheckedAccumulate(result.statistics.mam_occurrence_count,
                                  task.mam_occurrence_count,
                                  "MAM occurrence count");
                CheckedAccumulate(result.statistics.mem_occurrence_count,
                                  task.mem_occurrence_count,
                                  "boundary MEM occurrence count");
                CheckedAccumulate(result.statistics.mam_tile_raw_count,
                                  task.tile_raw_count,
                                  "tile raw MAM count");
                CheckedAccumulate(
                    result.statistics.mam_tile_globally_maximal_count,
                    task.tile_globally_maximal_count,
                    "globally maximal tile MAM count");
                CheckedAccumulate(
                    result.statistics.mam_boundary_mem_raw_count,
                    task.boundary_mem_raw_count,
                    "raw boundary MEM count");
                CheckedAccumulate(result.statistics.mam_boundary_pattern_count,
                                  task.boundary_pattern_count,
                                  "boundary pattern count");
                CheckedAccumulate(
                    result.statistics.mam_boundary_reference_unique_count,
                    task.boundary_reference_unique_count,
                    "reference-unique boundary pattern count");
                CheckedAccumulate(
                    result.statistics.mam_boundary_recovered_count,
                    task.boundary_recovered_count,
                    "boundary recovered MAM count");
                CheckedAccumulate(result_seed_count, task.seeds.size(),
                                  "reference-MAM result seed count");
            }
            result.statistics.raw_selected_seed_count =
                result.statistics.mam_occurrence_count;
            if (result_seed_count >
                static_cast<std::uint64_t>(
                    std::numeric_limits<std::size_t>::max())) {
                throw AlignmentError(
                    "reference-MAM result seed count exceeds the size range");
            }
            const std::uint64_t result_capacity_bytes = CheckedMultiply(
                result_seed_count, sizeof(Seed),
                "reference-MAM merged seed capacity");
            workspace_budget.Acquire(result_capacity_bytes,
                                     "merged seed storage");
            result.seeds.reserve(static_cast<std::size_t>(result_seed_count));
            if (result.seeds.capacity() > result_seed_count) {
                workspace_budget.Acquire(
                    CheckedMultiply(result.seeds.capacity() - result_seed_count,
                                    sizeof(Seed),
                                    "reference-MAM merged allocator growth"),
                    "allocator-expanded merged seed storage");
            }
            for (MamTaskWorkspace& task : tasks) {
                result.seeds.insert(
                    result.seeds.end(),
                    std::make_move_iterator(task.seeds.begin()),
                    std::make_move_iterator(task.seeds.end()));
                const std::uint64_t task_capacity_bytes = CheckedMultiply(
                    task.seeds.capacity(), sizeof(Seed),
                    "reference-MAM task seed capacity");
                std::vector<Seed>().swap(task.seeds);
                workspace_budget.Release(task_capacity_bytes);
            }
            result.statistics.mam_workspace_peak_bytes =
                workspace_budget.Peak();
        } else {
            for (const SequenceRecord& query : queries) {
                std::array<std::vector<Interval>, 2> mam_coverage;

            if (options.seed_mode == SeedMode::Mum) {
                implementation_->index.ForEachMum(
                    query.bases, mum_options,
                    [&](const sufkit::MumMatch& match) {
                        if (options.interruption_callback) {
                            options.interruption_callback("seed-enumeration");
                        }
                        Seed seed = ConvertMatch(
                            match, query, implementation_->reference_ids,
                            implementation_->reference_lengths);
                        if (seed.length < options.min_match) {
                            throw DependencyError(
                                "sufkit MUM callback returned a seed shorter "
                                "than min_match");
                        }
                        ++result.statistics.mum_occurrence_count;
                        ++result.statistics.raw_selected_seed_count;
                        result.seeds.push_back(std::move(seed));
                    });
            }

            if (options.seed_mode == SeedMode::Smem) {
                using SmemInterval =
                    std::tuple<Position, Length, Strand>;
                // expected reference occurrences, observed callback records
                std::map<SmemInterval,
                         std::pair<std::uint64_t, std::uint64_t>>
                    intervals;
                implementation_->index.ForEachSmem(
                    query.bases, smem_options,
                    [&](const sufkit::SmemMatch& match) {
                        if (options.interruption_callback) {
                            options.interruption_callback("seed-enumeration");
                        }
                        if (match.reference_occurrences == 0 ||
                            match.reference_occurrences <
                                options.smem_min_occurrences) {
                            throw DependencyError(
                                "sufkit SMEM callback returned an invalid "
                                "reference occurrence count");
                        }
                        Seed seed = ConvertMatch(
                            match, query, implementation_->reference_ids,
                            implementation_->reference_lengths);
                        if (seed.length < options.min_match) {
                            throw DependencyError(
                                "sufkit SMEM callback returned a seed shorter "
                                "than min_match");
                        }
                        const SmemInterval interval{seed.query_begin,
                                                    seed.length,
                                                    seed.strand};
                        auto [position, inserted] = intervals.try_emplace(
                            interval,
                            std::make_pair(match.reference_occurrences,
                                           std::uint64_t{0}));
                        if (!inserted &&
                            position->second.first !=
                                match.reference_occurrences) {
                            throw DependencyError(
                                "sufkit SMEM callback changed an interval's "
                                "reference occurrence count");
                        }
                        CheckedAccumulate(position->second.second, 1U,
                                          "SMEM interval callback count");
                        ++result.statistics.smem_occurrence_count;
                        ++result.statistics.raw_selected_seed_count;
                        result.seeds.push_back(std::move(seed));
                    });
                for (const auto& [interval, counts] : intervals) {
                    static_cast<void>(interval);
                    if (counts.first != counts.second) {
                        throw DependencyError(
                            "sufkit SMEM callback did not expand an interval "
                            "to every reference coordinate");
                    }
                }
                CheckedAccumulate(
                    result.statistics.smem_interval_count,
                    static_cast<std::uint64_t>(intervals.size()),
                    "SMEM interval count");
            }

            if (options.seed_mode == SeedMode::FastHierarchical) {
                implementation_->index.ForEachMam(
                    query.bases, mam_options,
                    [&](const sufkit::MamMatch& match) {
                        if (options.interruption_callback) {
                            options.interruption_callback("seed-enumeration");
                        }
                        Seed seed = ConvertMatch(
                            match, query, implementation_->reference_ids,
                            implementation_->reference_lengths);
                        if (seed.length < options.min_match) {
                            throw AlignmentError(
                                "sufkit MAM callback returned a seed shorter "
                                "than min_match");
                        }
                        ++result.statistics.mam_occurrence_count;
                        ++result.statistics.raw_selected_seed_count;
                        if (options.seed_mode == SeedMode::FastHierarchical) {
                            mam_coverage[StrandIndex(seed.strand)].push_back(
                                {seed.query_begin,
                                 CheckedEnd(seed.query_begin, seed.length,
                                            "MAM query interval")});
                        }
                        result.seeds.push_back(std::move(seed));
                    });
            }

            if (options.seed_mode == SeedMode::FastHierarchical) {
                MergeIntervals(mam_coverage[0]);
                MergeIntervals(mam_coverage[1]);
            }

            if (options.seed_mode == SeedMode::FastHierarchical ||
                options.seed_mode == SeedMode::MaxMatch) {
                implementation_->index.ForEachMem(
                    query.bases, mem_options,
                    [&](const sufkit::MemMatch& match) {
                        if (options.interruption_callback) {
                            options.interruption_callback("seed-enumeration");
                        }
                        Seed seed = ConvertMatch(
                            match, query, implementation_->reference_ids,
                            implementation_->reference_lengths);
                        if (seed.length < options.min_match) {
                            throw AlignmentError(
                                "sufkit MEM callback returned a seed shorter "
                                "than min_match");
                        }
                        ++result.statistics.mem_occurrence_count;
                        if (options.seed_mode == SeedMode::FastHierarchical) {
                            const Position end = CheckedEnd(
                                seed.query_begin, seed.length,
                                "MEM query interval");
                            if (IsCovered(mam_coverage[StrandIndex(seed.strand)],
                                          seed.query_begin, end)) {
                                ++result.statistics.mem_covered_by_mam_count;
                                return;
                            }
                        }
                        ++result.statistics.raw_selected_seed_count;
                        result.seeds.push_back(std::move(seed));
                    });
                }
            }
        }

        std::sort(result.seeds.begin(), result.seeds.end(), SeedLess);
        result.seeds.erase(
            std::unique(result.seeds.begin(), result.seeds.end()),
            result.seeds.end());
        result.statistics.selected_seed_count = result.seeds.size();
        result.statistics.enumeration_seconds =
            ElapsedSeconds(enumeration_begin);
        return result;
    } catch (const InterruptedError&) {
        throw;
    } catch (const DependencyError&) {
        throw;
    } catch (const AlignmentError&) {
        throw;
    } catch (const std::bad_alloc&) {
        throw;
    } catch (const sufkit::Error& error) {
        throw DependencyError(
            SufkitFailure("MEM/MAM/MUM/SMEM enumeration", error));
    } catch (const std::exception& error) {
        throw DependencyError("sufkit MEM/MAM enumeration failed: " +
                              std::string{error.what()});
    }
#else
    static_cast<void>(references);
    static_cast<void>(queries);
    static_cast<void>(options);
    ThrowUnavailable();
#endif
}

SufkitIndexStatistics SufkitSeedIndex::BuildStatistics() const {
#if RAMAG_HAVE_SUFKIT
    if (!implementation_) {
        throw AlignmentError(
            "cannot inspect a moved-from sufkit seed index");
    }
    return implementation_->statistics;
#else
    ThrowUnavailable();
#endif
}

}  // namespace ramag

#include "ramag/alignment.hpp"

#include "extension_backend.hpp"

#include "ramag/fasta.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <exception>
#include <functional>
#include <limits>
#include <memory>
#include <queue>
#include <sstream>
#include <string_view>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#if defined(_OPENMP)
#include <omp.h>
#endif

namespace ramag {
namespace {

using Clock = std::chrono::steady_clock;

[[nodiscard]] double SecondsBetween(Clock::time_point begin, Clock::time_point end) {
    return std::chrono::duration<double>(end - begin).count();
}

[[nodiscard]] Position CheckedEnd(Position begin, Length length, std::string_view what) {
    if (length > std::numeric_limits<Position>::max() - begin) {
        throw AlignmentError(std::string{what} + " exceeds the 64-bit coordinate range");
    }
    return begin + length;
}

[[nodiscard]] Length CheckedLengthAdd(Length left, Length right, std::string_view what) {
    if (right > std::numeric_limits<Length>::max() - left) {
        throw AlignmentError(std::string{what} + " exceeds the 64-bit coordinate range");
    }
    return left + right;
}

[[nodiscard]] std::int64_t CheckedScoreProduct(Length length,
                                               std::int32_t value,
                                               std::string_view what) {
    const std::int64_t absolute_value =
        value < 0 ? -static_cast<std::int64_t>(value) : value;
    if (absolute_value != 0 &&
        length > static_cast<Length>(std::numeric_limits<std::int64_t>::max() /
                                     absolute_value)) {
        throw AlignmentError(std::string{what} + " exceeds the signed 64-bit score range");
    }
    const auto score = static_cast<std::int64_t>(length) * absolute_value;
    return value < 0 ? -score : score;
}

[[nodiscard]] std::int64_t CheckedScoreAdd(std::int64_t left,
                                           std::int64_t right,
                                           std::string_view what) {
    if ((right > 0 && left > std::numeric_limits<std::int64_t>::max() - right) ||
        (right < 0 && left < std::numeric_limits<std::int64_t>::min() - right)) {
        throw AlignmentError(std::string{what} + " exceeds the signed 64-bit score range");
    }
    return left + right;
}

void ValidateOptions(const AlignmentOptions& options) {
    if (options.min_match == 0) {
        throw AlignmentError("min_match must be greater than zero");
    }
    if (options.smem_min_occurrences == 0) {
        throw AlignmentError("smem_min_occurrences must be greater than zero");
    }
    if (options.min_cluster == 0) {
        throw AlignmentError("min_cluster must be greater than zero");
    }
    if (options.break_length == 0) {
        throw AlignmentError("break_length must be greater than zero");
    }
    if (!std::isfinite(options.diag_factor) || options.diag_factor < 0.0) {
        throw AlignmentError("diag_factor must be finite and non-negative");
    }
    if (options.max_dp_cells == 0) {
        throw AlignmentError("max_dp_cells must be greater than zero");
    }
    if (options.match_score <= 0 || options.mismatch_penalty < 0 ||
        options.gap_open_penalty < 0 || options.gap_extend_penalty < 0) {
        throw AlignmentError(
            "match_score must be positive and alignment penalties must be non-negative");
    }
    if (options.worker_threads == 0) {
        throw AlignmentError("worker_threads must be greater than zero");
    }
}

void ValidateSequenceCollection(const std::vector<SequenceRecord>& sequences,
                                std::string_view role) {
    if (sequences.empty()) {
        throw AlignmentError(std::string{role} + " sequence collection is empty");
    }
    std::unordered_set<SequenceId> numeric_ids;
    std::unordered_set<std::string> names;
    for (const SequenceRecord& sequence : sequences) {
        if (sequence.name.empty()) {
            throw AlignmentError(std::string{role} + " contains an empty sequence name");
        }
        if (sequence.bases.empty()) {
            throw AlignmentError(std::string{role} + " sequence '" + sequence.name +
                                 "' is empty");
        }
        if (!numeric_ids.insert(sequence.numeric_id).second) {
            throw AlignmentError(std::string{role} + " contains duplicate numeric id " +
                                 std::to_string(sequence.numeric_id));
        }
        if (!names.insert(sequence.name).second) {
            throw AlignmentError(std::string{role} + " contains duplicate name '" +
                                 sequence.name + "'");
        }
        for (char base : sequence.bases) {
            if (base != 'A' && base != 'C' && base != 'G' && base != 'T' && base != 'N') {
                throw AlignmentError(std::string{role} + " sequence '" + sequence.name +
                                     "' is not normalized to upper-case A/C/G/T/N");
            }
        }
    }
}

[[nodiscard]] Length CountAmbiguous(const std::vector<SequenceRecord>& sequences) {
    Length count = 0;
    for (const SequenceRecord& sequence : sequences) {
        const auto local = static_cast<Length>(
            std::count(sequence.bases.begin(), sequence.bases.end(), 'N'));
        count = CheckedLengthAdd(count, local, "ambiguous-base count");
    }
    return count;
}

[[nodiscard]] Length CountBases(const std::vector<SequenceRecord>& sequences) {
    Length count = 0;
    for (const SequenceRecord& sequence : sequences) {
        count = CheckedLengthAdd(count, sequence.size(), "sequence base count");
    }
    return count;
}

[[nodiscard]] const char* SeedRouteName(SeedMode mode) {
    switch (mode) {
        case SeedMode::FastHierarchical:
            return "mam-skeleton+whole-query-mem-filter";
        case SeedMode::MumReference:
            return "mumreference";
        case SeedMode::Mum:
            return "mum";
        case SeedMode::Smem:
            return "smem";
        case SeedMode::MaxMatch:
            return "maxmatch";
    }
    return "unknown";
}

[[nodiscard]] bool EqualCanonical(char left, char right) noexcept {
    return left == right && IsCanonicalBase(left);
}

struct RefLocation {
    std::size_t record_index{};
    std::size_t begin{};
};

struct StringViewHash {
    [[nodiscard]] std::size_t operator()(std::string_view value) const noexcept {
        return std::hash<std::string_view>{}(value);
    }
};

using KmerBuckets =
    std::unordered_map<std::string_view, std::vector<RefLocation>, StringViewHash>;

[[nodiscard]] std::vector<std::size_t> CanonicalRuns(std::string_view bases) {
    std::vector<std::size_t> runs(bases.size() + 1, 0);
    for (std::size_t index = bases.size(); index > 0; --index) {
        const std::size_t position = index - 1;
        runs[position] = IsCanonicalBase(bases[position]) ? runs[position + 1] + 1 : 0;
    }
    return runs;
}

[[nodiscard]] KmerBuckets BuildReferenceKmers(
    const std::vector<SequenceRecord>& references,
    std::size_t kmer_length) {
    KmerBuckets buckets;
    for (std::size_t reference_index = 0; reference_index < references.size();
         ++reference_index) {
        const std::string& bases = references[reference_index].bases;
        if (bases.size() < kmer_length) {
            continue;
        }
        const std::vector<std::size_t> runs = CanonicalRuns(bases);
        for (std::size_t begin = 0; begin + kmer_length <= bases.size(); ++begin) {
            if (runs[begin] < kmer_length) {
                continue;
            }
            buckets[std::string_view{bases}.substr(begin, kmer_length)].push_back(
                RefLocation{reference_index, begin});
        }
    }
    return buckets;
}

[[nodiscard]] std::vector<RefLocation> FindReferenceOccurrences(
    const std::vector<SequenceRecord>& references,
    std::string_view pattern) {
    std::vector<RefLocation> occurrences;
    for (std::size_t reference_index = 0;
         reference_index < references.size(); ++reference_index) {
        const std::string& bases = references[reference_index].bases;
        if (pattern.size() > bases.size()) {
            continue;
        }
        for (std::size_t begin = 0;
             begin + pattern.size() <= bases.size(); ++begin) {
            if (std::equal(pattern.begin(), pattern.end(),
                           bases.begin() +
                               static_cast<std::ptrdiff_t>(begin))) {
                occurrences.push_back({reference_index, begin});
            }
        }
    }
    return occurrences;
}

struct RawSeed {
    Seed seed;
    std::size_t reference_index{};
    std::size_t query_index{};
    Position oriented_query_begin{};
};

[[nodiscard]] Position ForwardQueryBegin(Length query_length,
                                         Position oriented_begin,
                                         Length seed_length,
                                         Strand strand) {
    const Position oriented_end =
        CheckedEnd(oriented_begin, seed_length, "oriented query seed");
    if (oriented_end > query_length) {
        throw AlignmentError("seed lies outside its query sequence");
    }
    return strand == Strand::Forward ? oriented_begin : query_length - oriented_end;
}

[[nodiscard]] std::uint64_t CountQueryOccurrences(
    std::string_view query,
    std::string_view pattern) {
    if (pattern.size() > query.size()) {
        return 0;
    }
    std::uint64_t count = 0;
    for (std::size_t begin = 0;
         begin + pattern.size() <= query.size(); ++begin) {
        if (std::equal(pattern.begin(), pattern.end(),
                       query.begin() +
                           static_cast<std::ptrdiff_t>(begin))) {
            count = CheckedLengthAdd(count, 1U,
                                     "query pattern occurrence count");
        }
    }
    return count;
}

struct QualifyingSmem {
    std::size_t begin{};
    std::size_t end{};
    std::vector<RefLocation> occurrences;
};

struct SmemOracleResult {
    std::vector<Seed> seeds;
    std::uint64_t interval_count{};
};

[[nodiscard]] SmemOracleResult EnumerateSmemOracle(
    const std::vector<SequenceRecord>& references,
    const std::vector<SequenceRecord>& queries,
    std::size_t min_length,
    std::uint64_t min_occurrences) {
    SmemOracleResult result;
    for (const SequenceRecord& query : queries) {
        const std::string reverse = ReverseComplement(query.bases);
        for (const Strand strand : {Strand::Forward, Strand::Reverse}) {
            const std::string_view oriented =
                strand == Strand::Forward ? std::string_view{query.bases}
                                          : std::string_view{reverse};
            const std::vector<std::size_t> runs = CanonicalRuns(oriented);
            std::vector<QualifyingSmem> qualifying;
            for (std::size_t begin = 0; begin < oriented.size(); ++begin) {
                if (runs[begin] < min_length) {
                    continue;
                }
                const std::size_t run_end = begin + runs[begin];
                for (std::size_t end = begin + min_length;
                     end <= run_end; ++end) {
                    auto occurrences = FindReferenceOccurrences(
                        references, oriented.substr(begin, end - begin));
                    if (occurrences.size() >= min_occurrences) {
                        qualifying.push_back(
                            {begin, end, std::move(occurrences)});
                    }
                }
            }
            for (std::size_t index = 0; index < qualifying.size(); ++index) {
                const QualifyingSmem& candidate = qualifying[index];
                bool contained = false;
                for (std::size_t other_index = 0;
                     other_index < qualifying.size(); ++other_index) {
                    if (index == other_index) {
                        continue;
                    }
                    const QualifyingSmem& other = qualifying[other_index];
                    if (other.begin <= candidate.begin &&
                        other.end >= candidate.end &&
                        (other.begin < candidate.begin ||
                         other.end > candidate.end)) {
                        contained = true;
                        break;
                    }
                }
                if (contained) {
                    continue;
                }
                result.interval_count = CheckedLengthAdd(
                    result.interval_count, 1U, "SMEM interval count");
                const Length length =
                    static_cast<Length>(candidate.end - candidate.begin);
                const Position query_begin = ForwardQueryBegin(
                    query.size(), static_cast<Position>(candidate.begin),
                    length, strand);
                for (const RefLocation& occurrence :
                     candidate.occurrences) {
                    result.seeds.push_back(
                        Seed{references[occurrence.record_index].numeric_id,
                             static_cast<Position>(occurrence.begin),
                             query.numeric_id, query_begin, length, strand});
                }
            }
        }
    }
    std::sort(result.seeds.begin(), result.seeds.end(),
              [](const Seed& left, const Seed& right) {
                  return std::tie(left.reference_id, left.query_id,
                                  left.strand, left.query_begin,
                                  left.reference_begin, left.length) <
                         std::tie(right.reference_id, right.query_id,
                                  right.strand, right.query_begin,
                                  right.reference_begin, right.length);
              });
    result.seeds.erase(
        std::unique(result.seeds.begin(), result.seeds.end()),
        result.seeds.end());
    return result;
}

void EnumerateOrientedMems(const std::vector<SequenceRecord>& references,
                           const SequenceRecord& query,
                           std::size_t query_index,
                           std::string_view oriented_query,
                           Strand strand,
                           std::size_t kmer_length,
                           const KmerBuckets& reference_kmers,
                           std::vector<RawSeed>& output) {
    if (oriented_query.size() < kmer_length) {
        return;
    }
    const std::vector<std::size_t> query_runs = CanonicalRuns(oriented_query);
    for (std::size_t query_begin = 0;
         query_begin + kmer_length <= oriented_query.size();
         ++query_begin) {
        if (query_runs[query_begin] < kmer_length) {
            continue;
        }
        const std::string_view key = oriented_query.substr(query_begin, kmer_length);
        const auto bucket = reference_kmers.find(key);
        if (bucket == reference_kmers.end()) {
            continue;
        }
        for (const RefLocation& location : bucket->second) {
            const SequenceRecord& reference = references[location.record_index];
            if (location.begin > 0 && query_begin > 0 &&
                EqualCanonical(reference.bases[location.begin - 1],
                               oriented_query[query_begin - 1])) {
                // This k-mer is inside a longer match.  Only its left-maximal
                // start is allowed to produce the occurrence-level MEM.
                continue;
            }

            std::size_t length = kmer_length;
            while (location.begin + length < reference.bases.size() &&
                   query_begin + length < oriented_query.size() &&
                   EqualCanonical(reference.bases[location.begin + length],
                                  oriented_query[query_begin + length])) {
                ++length;
            }

            const Length public_length = static_cast<Length>(length);
            const Position oriented_position = static_cast<Position>(query_begin);
            output.push_back(RawSeed{
                Seed{reference.numeric_id,
                     static_cast<Position>(location.begin),
                     query.numeric_id,
                     ForwardQueryBegin(query.size(), oriented_position, public_length, strand),
                     public_length,
                     strand},
                location.record_index,
                query_index,
                oriented_position});
        }
    }
}

[[nodiscard]] bool IsReferenceUnique(const RawSeed& raw_seed,
                                     const std::vector<SequenceRecord>& references,
                                     std::size_t kmer_length,
                                     const KmerBuckets& reference_kmers) {
    const SequenceRecord& source = references[raw_seed.reference_index];
    const std::size_t begin = static_cast<std::size_t>(raw_seed.seed.reference_begin);
    const std::size_t length = static_cast<std::size_t>(raw_seed.seed.length);
    const std::string_view pattern{source.bases.data() + begin, length};
    const auto bucket = reference_kmers.find(pattern.substr(0, kmer_length));
    if (bucket == reference_kmers.end()) {
        throw AlignmentError("internal error: seed prefix is missing from reference index");
    }

    std::uint32_t occurrences = 0;
    for (const RefLocation& candidate : bucket->second) {
        const std::string& candidate_bases = references[candidate.record_index].bases;
        if (length > candidate_bases.size() - candidate.begin) {
            continue;
        }
        if (std::equal(pattern.begin(), pattern.end(), candidate_bases.begin() +
                                                        static_cast<std::ptrdiff_t>(candidate.begin))) {
            ++occurrences;
            if (occurrences > 1) {
                return false;
            }
        }
    }
    return occurrences == 1;
}

[[nodiscard]] bool IsQueryUnique(
    const RawSeed& raw_seed,
    const std::vector<SequenceRecord>& references,
    const std::vector<SequenceRecord>& queries) {
    const SequenceRecord& reference = references[raw_seed.reference_index];
    const SequenceRecord& query = queries[raw_seed.query_index];
    const std::size_t reference_begin =
        static_cast<std::size_t>(raw_seed.seed.reference_begin);
    const std::size_t length =
        static_cast<std::size_t>(raw_seed.seed.length);
    const std::string_view pattern{reference.bases.data() + reference_begin,
                                   length};
    if (raw_seed.seed.strand == Strand::Forward) {
        return CountQueryOccurrences(query.bases, pattern) == 1;
    }
    const std::string reverse = ReverseComplement(query.bases);
    return CountQueryOccurrences(reverse, pattern) == 1;
}

[[nodiscard]] bool PublicSeedLess(const Seed& left, const Seed& right) {
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

[[nodiscard]] std::uint64_t QueryStrandKey(const Seed& seed) noexcept {
    return (static_cast<std::uint64_t>(seed.query_id) << 1U) |
           static_cast<std::uint64_t>(seed.strand == Strand::Reverse);
}

using QueryLengths = std::unordered_map<SequenceId, Length>;

[[nodiscard]] QueryLengths MakeQueryLengths(
    const std::vector<SequenceRecord>& queries) {
    QueryLengths lengths;
    for (const SequenceRecord& query : queries) {
        if (!lengths.emplace(query.numeric_id, query.size()).second) {
            throw AlignmentError("duplicate query numeric id " +
                                 std::to_string(query.numeric_id));
        }
    }
    return lengths;
}

[[nodiscard]] Position OrientedQueryBegin(const Seed& seed, Length query_length) {
    const Position query_end = CheckedEnd(seed.query_begin, seed.length, "query seed");
    if (query_end > query_length) {
        throw AlignmentError("seed lies outside query id " + std::to_string(seed.query_id));
    }
    return seed.strand == Strand::Forward ? seed.query_begin : query_length - query_end;
}

struct OrientedSeed {
    Seed seed;
    Position query_begin{};
};

[[nodiscard]] bool SameSeedGroup(const OrientedSeed& left,
                                 const OrientedSeed& right) noexcept {
    return left.seed.reference_id == right.seed.reference_id &&
           left.seed.query_id == right.seed.query_id &&
           left.seed.strand == right.seed.strand;
}

struct SignedMagnitude {
    bool negative{};
    std::uint64_t magnitude{};
};

[[nodiscard]] SignedMagnitude DiagonalValue(const OrientedSeed& seed) noexcept {
    if (seed.seed.reference_begin >= seed.query_begin) {
        return SignedMagnitude{false, seed.seed.reference_begin - seed.query_begin};
    }
    return SignedMagnitude{true, seed.query_begin - seed.seed.reference_begin};
}

[[nodiscard]] bool SameDiagonal(const OrientedSeed& left,
                                const OrientedSeed& right) noexcept {
    const SignedMagnitude left_value = DiagonalValue(left);
    const SignedMagnitude right_value = DiagonalValue(right);
    return left_value.negative == right_value.negative &&
           left_value.magnitude == right_value.magnitude;
}

[[nodiscard]] bool DiagonalLess(const OrientedSeed& left,
                                const OrientedSeed& right) noexcept {
    const SignedMagnitude left_value = DiagonalValue(left);
    const SignedMagnitude right_value = DiagonalValue(right);
    if (left_value.negative != right_value.negative) {
        return left_value.negative;
    }
    return left_value.negative ? left_value.magnitude > right_value.magnitude
                               : left_value.magnitude < right_value.magnitude;
}

[[nodiscard]] Seed ToPublicSeed(const OrientedSeed& oriented, Length query_length) {
    Seed result = oriented.seed;
    result.query_begin = ForwardQueryBegin(
        query_length, oriented.query_begin, oriented.seed.length, oriented.seed.strand);
    return result;
}

[[nodiscard]] bool FitsDpBound(Length reference_gap,
                               Length query_gap,
                               std::uint64_t max_cells) noexcept {
    if (reference_gap == 0 || query_gap == 0 || reference_gap == query_gap) {
        return true;
    }
    if (reference_gap == std::numeric_limits<Length>::max() ||
        query_gap == std::numeric_limits<Length>::max()) {
        return false;
    }
    const std::uint64_t rows = reference_gap + 1;
    const std::uint64_t columns = query_gap + 1;
    return rows <= max_cells / columns;
}

[[nodiscard]] bool CanChain(const OrientedSeed& previous,
                            const OrientedSeed& next,
                            const AlignmentOptions& options) {
    const Position previous_reference_end =
        CheckedEnd(previous.seed.reference_begin, previous.seed.length, "reference seed");
    const Position previous_query_end =
        CheckedEnd(previous.query_begin, previous.seed.length, "oriented query seed");
    if (previous_reference_end > next.seed.reference_begin ||
        previous_query_end > next.query_begin) {
        return false;
    }
    const Length reference_gap = next.seed.reference_begin - previous_reference_end;
    const Length query_gap = next.query_begin - previous_query_end;
    const Length separation = std::max(reference_gap, query_gap);
    if (separation > options.max_gap || separation > options.break_length) {
        return false;
    }
    const Length diagonal_change = reference_gap > query_gap
                                       ? reference_gap - query_gap
                                       : query_gap - reference_gap;
    const long double allowed = static_cast<long double>(options.diag_diff) +
                                static_cast<long double>(options.diag_factor) *
                                    static_cast<long double>(separation);
    if (static_cast<long double>(diagonal_change) > allowed) {
        return false;
    }
    return FitsDpBound(reference_gap, query_gap, options.max_dp_cells);
}

using Chain = std::vector<OrientedSeed>;

[[nodiscard]] std::vector<OrientedSeed> MakeSortedOrientedSeeds(
    const std::vector<Seed>& seeds,
    const QueryLengths& query_lengths) {
    std::vector<OrientedSeed> oriented;
    oriented.reserve(seeds.size());
    for (const Seed& seed : seeds) {
        const auto length = query_lengths.find(seed.query_id);
        if (length == query_lengths.end()) {
            throw AlignmentError("seed refers to unknown query id " +
                                 std::to_string(seed.query_id));
        }
        oriented.push_back(
            OrientedSeed{seed, OrientedQueryBegin(seed, length->second)});
    }
    std::sort(oriented.begin(), oriented.end(), [](const OrientedSeed& left,
                                                   const OrientedSeed& right) {
        return std::tie(left.seed.reference_id,
                        left.seed.query_id,
                        left.seed.strand,
                        left.seed.reference_begin,
                        left.query_begin,
                        left.seed.length) <
               std::tie(right.seed.reference_id,
                        right.seed.query_id,
                        right.seed.strand,
                        right.seed.reference_begin,
                        right.query_begin,
                        right.seed.length);
    });
    return oriented;
}

// Retained only as a small-fixture semantic oracle.  The production alignment
// path below never selects or falls back to this quadratic implementation.
[[nodiscard]] std::vector<Chain> BuildChainsQuadraticOracle(
    const std::vector<Seed>& seeds,
    const QueryLengths& query_lengths,
    const AlignmentOptions& options) {
    const std::vector<OrientedSeed> oriented =
        MakeSortedOrientedSeeds(seeds, query_lengths);

    std::vector<Chain> chains;
    std::size_t group_begin = 0;
    while (group_begin < oriented.size()) {
        std::size_t group_end = group_begin + 1;
        while (group_end < oriented.size() &&
               SameSeedGroup(oriented[group_begin], oriented[group_end])) {
            ++group_end;
        }

        const std::size_t group_size = group_end - group_begin;
        std::vector<bool> used(group_size, false);
        while (true) {
            std::vector<Length> best_score(group_size, 0);
            std::vector<std::size_t> predecessor(group_size, group_size);
            std::size_t best_end = group_size;

            for (std::size_t local_i = 0; local_i < group_size; ++local_i) {
                if (used[local_i]) {
                    continue;
                }
                const OrientedSeed& current = oriented[group_begin + local_i];
                best_score[local_i] = current.seed.length;
                for (std::size_t local_j = 0; local_j < local_i; ++local_j) {
                    if (used[local_j] || best_score[local_j] == 0) {
                        continue;
                    }
                    const OrientedSeed& previous = oriented[group_begin + local_j];
                    if (!CanChain(previous, current, options)) {
                        continue;
                    }
                    const Length candidate =
                        CheckedLengthAdd(best_score[local_j], current.seed.length,
                                         "chain seed length");
                    if (candidate > best_score[local_i] ||
                        (candidate == best_score[local_i] &&
                         local_j < predecessor[local_i])) {
                        best_score[local_i] = candidate;
                        predecessor[local_i] = local_j;
                    }
                }
                if (best_end == group_size || best_score[local_i] > best_score[best_end] ||
                    (best_score[local_i] == best_score[best_end] && local_i < best_end)) {
                    best_end = local_i;
                }
            }

            if (best_end == group_size || best_score[best_end] < options.min_cluster) {
                break;
            }

            std::vector<std::size_t> chain_indices;
            for (std::size_t cursor = best_end; cursor != group_size;
                 cursor = predecessor[cursor]) {
                chain_indices.push_back(cursor);
            }
            std::reverse(chain_indices.begin(), chain_indices.end());

            Chain chain;
            chain.reserve(chain_indices.size());
            for (std::size_t local_index : chain_indices) {
                used[local_index] = true;
                chain.push_back(oriented[group_begin + local_index]);
            }
            chains.push_back(std::move(chain));
        }
        group_begin = group_end;
    }
    return chains;
}

struct EndpointCell {
    std::uint64_t reference{};
    std::uint64_t query{};

    friend bool operator==(const EndpointCell&, const EndpointCell&) = default;
};

struct EndpointCellHash {
    [[nodiscard]] std::size_t operator()(const EndpointCell& cell) const noexcept {
        // SplitMix64 finalizers make the pair hash insensitive to regular
        // genomic coordinate strides.  Hash traversal order never participates
        // in a chaining tie-break.
        auto mix = [](std::uint64_t value) noexcept {
            value += 0x9e3779b97f4a7c15ULL;
            value = (value ^ (value >> 30U)) * 0xbf58476d1ce4e5b9ULL;
            value = (value ^ (value >> 27U)) * 0x94d049bb133111ebULL;
            return value ^ (value >> 31U);
        };
        const std::uint64_t left = mix(cell.reference);
        const std::uint64_t right = mix(cell.query);
        return static_cast<std::size_t>(left ^ (right + 0x9e3779b97f4a7c15ULL +
                                                (left << 6U) + (left >> 2U)));
    }
};

class DisjointSets {
public:
    explicit DisjointSets(std::size_t size) : parent_(size), rank_(size, 0) {
        for (std::size_t index = 0; index < size; ++index) {
            parent_[index] = index;
        }
    }

    [[nodiscard]] std::size_t Find(std::size_t value) {
        std::size_t root = value;
        while (parent_[root] != root) {
            root = parent_[root];
        }
        while (parent_[value] != value) {
            const std::size_t next = parent_[value];
            parent_[value] = root;
            value = next;
        }
        return root;
    }

    void Unite(std::size_t left, std::size_t right) {
        left = Find(left);
        right = Find(right);
        if (left == right) {
            return;
        }
        if (rank_[left] < rank_[right] ||
            (rank_[left] == rank_[right] && left > right)) {
            std::swap(left, right);
        }
        parent_[right] = left;
        if (rank_[left] == rank_[right]) {
            ++rank_[left];
        }
    }

private:
    std::vector<std::size_t> parent_;
    std::vector<std::uint8_t> rank_;
};

[[nodiscard]] std::string ChainingGroupDescription(const OrientedSeed& seed) {
    std::ostringstream output;
    output << "reference_id=" << seed.seed.reference_id
           << ", query_id=" << seed.seed.query_id << ", strand="
           << (seed.seed.strand == Strand::Forward ? "forward" : "reverse");
    return output.str();
}

[[noreturn]] void ThrowChainingLimit(std::string_view resource,
                                     std::uint64_t observed,
                                     std::uint64_t limit,
                                     std::string_view group) {
    std::ostringstream message;
    message << "chaining resource limit exceeded for " << group << ": "
            << resource << " observed=" << observed << ", limit=" << limit
            << "; exact chaining was aborted without truncation or sampling";
    throw AlignmentError(message.str());
}

void IncrementUnbounded(std::uint64_t& value, std::string_view resource) {
    if (value == std::numeric_limits<std::uint64_t>::max()) {
        throw AlignmentError(std::string{resource} + " exceeds the 64-bit counter range");
    }
    ++value;
}

[[nodiscard]] std::uint64_t CheckedByteProduct(std::size_t count,
                                               std::size_t width) {
    const std::uint64_t public_count = static_cast<std::uint64_t>(count);
    const std::uint64_t public_width = static_cast<std::uint64_t>(width);
    if (public_width != 0 &&
        public_count > std::numeric_limits<std::uint64_t>::max() / public_width) {
        throw AlignmentError("chaining working-set estimate exceeds 64-bit range");
    }
    return public_count * public_width;
}

[[nodiscard]] std::uint64_t CheckedByteAdd(std::uint64_t left,
                                           std::uint64_t right) {
    if (right > std::numeric_limits<std::uint64_t>::max() - left) {
        throw AlignmentError("chaining working-set estimate exceeds 64-bit range");
    }
    return left + right;
}

class alignas(64) SharedWorkingSetBudget {
public:
    explicit SharedWorkingSetBudget(std::uint64_t limit) : limit_(limit) {}

    void Acquire(std::uint64_t bytes, std::string_view group) {
        std::uint64_t current = current_.load(std::memory_order_relaxed);
        while (true) {
            if (bytes > std::numeric_limits<std::uint64_t>::max() - current) {
                ThrowChainingLimit("estimated working-set bytes",
                                   std::numeric_limits<std::uint64_t>::max(),
                                   limit_, group);
            }
            const std::uint64_t observed = current + bytes;
            if (observed > limit_) {
                ThrowChainingLimit("estimated working-set bytes", observed,
                                   limit_, group);
            }
            if (current_.compare_exchange_weak(current, observed,
                                               std::memory_order_acq_rel,
                                               std::memory_order_relaxed)) {
                std::uint64_t peak = peak_.load(std::memory_order_relaxed);
                while (peak < observed &&
                       !peak_.compare_exchange_weak(peak, observed,
                                                    std::memory_order_relaxed,
                                                    std::memory_order_relaxed)) {
                }
                return;
            }
        }
    }

    void Release(std::uint64_t bytes) noexcept {
        if (bytes != 0) {
            static_cast<void>(current_.fetch_sub(bytes, std::memory_order_acq_rel));
        }
    }

    [[nodiscard]] std::uint64_t peak() const noexcept {
        return peak_.load(std::memory_order_relaxed);
    }

private:
    std::uint64_t limit_{};
    alignas(64) std::atomic<std::uint64_t> current_{0};
    alignas(64) std::atomic<std::uint64_t> peak_{0};
};

// Candidate checks, legal edges, and DP relaxations are task-wide hard
// budgets.  A per-group limit followed by end-of-stage aggregation is exact,
// but it can let every parallel group independently consume the full budget
// before the failure is discovered.  Reserve each unit globally at the point
// of work so cancellation remains fail-fast and no completed run can have
// transiently exceeded the advertised resource limit.
class alignas(64) SharedCountBudget {
public:
    SharedCountBudget(std::uint64_t limit, std::string_view resource)
        : limit_(limit), resource_(resource) {}

    void Acquire(std::string_view group) {
        std::uint64_t current = current_.load(std::memory_order_relaxed);
        while (true) {
            if (current >= limit_) {
                const std::uint64_t observed =
                    limit_ == std::numeric_limits<std::uint64_t>::max()
                        ? limit_
                        : limit_ + 1;
                ThrowChainingLimit(resource_, observed, limit_, group);
            }
            if (current_.compare_exchange_weak(current, current + 1,
                                               std::memory_order_acq_rel,
                                               std::memory_order_relaxed)) {
                return;
            }
        }
    }

    [[nodiscard]] std::uint64_t value() const noexcept {
        return current_.load(std::memory_order_relaxed);
    }

private:
    std::uint64_t limit_{};
    std::string resource_;
    std::atomic<std::uint64_t> current_{0};
};

void IncrementGloballyLimited(std::uint64_t& local_value,
                              SharedCountBudget& global_budget,
                              std::string_view resource,
                              std::string_view group) {
    global_budget.Acquire(group);
    IncrementUnbounded(local_value, resource);
}

class GroupWorkingReservation {
public:
    GroupWorkingReservation(SharedWorkingSetBudget& budget, std::string_view group)
        : budget_(budget), group_(group) {}

    GroupWorkingReservation(const GroupWorkingReservation&) = delete;
    GroupWorkingReservation& operator=(const GroupWorkingReservation&) = delete;

    ~GroupWorkingReservation() {
        budget_.Release(temporary_bytes_);
        if (!committed_) {
            budget_.Release(retained_bytes_);
        }
    }

    void AddTemporary(std::uint64_t bytes) {
        budget_.Acquire(bytes, group_);
        temporary_bytes_ = CheckedByteAdd(temporary_bytes_, bytes);
        temporary_peak_bytes_ =
            std::max(temporary_peak_bytes_, temporary_bytes_);
    }

    void ReleaseTemporary(std::uint64_t bytes) {
        if (bytes > temporary_bytes_) {
            throw AlignmentError(
                "internal error: temporary chaining reservation underflow");
        }
        budget_.Release(bytes);
        temporary_bytes_ -= bytes;
    }

    void AddRetained(std::uint64_t bytes) {
        budget_.Acquire(bytes, group_);
        retained_bytes_ = CheckedByteAdd(retained_bytes_, bytes);
    }

    void CommitRetained() {
        budget_.Release(temporary_bytes_);
        temporary_bytes_ = 0;
        committed_ = true;
    }

    [[nodiscard]] std::uint64_t temporary_bytes() const noexcept {
        return temporary_peak_bytes_;
    }

    [[nodiscard]] std::uint64_t retained_bytes() const noexcept {
        return retained_bytes_;
    }

private:
    SharedWorkingSetBudget& budget_;
    std::string group_;
    std::uint64_t temporary_bytes_{};
    std::uint64_t temporary_peak_bytes_{};
    std::uint64_t retained_bytes_{};
    bool committed_{};
};

class ChainingCancelled final : public std::exception {
public:
    [[nodiscard]] const char* what() const noexcept override {
        return "parallel chaining cancelled after another group failed";
    }
};

void CheckChainingCancellation(const std::atomic<bool>* cancellation) {
    if (cancellation != nullptr &&
        cancellation->load(std::memory_order_relaxed)) {
        throw ChainingCancelled{};
    }
}

[[nodiscard]] std::uint64_t CellCoordinate(Position value, Length gap_bound) noexcept {
    if (gap_bound == std::numeric_limits<Length>::max()) {
        return 0;
    }
    return value / (gap_bound + 1);
}

struct CellRange {
    std::array<std::uint64_t, 2> values{};
    std::size_t size{};
};

[[nodiscard]] CellRange EndpointSearchCells(Position begin, Length gap_bound) noexcept {
    const Position lower = begin > gap_bound ? begin - gap_bound : 0;
    const std::uint64_t lower_cell = CellCoordinate(lower, gap_bound);
    const std::uint64_t upper_cell = CellCoordinate(begin, gap_bound);
    CellRange result;
    result.values[0] = lower_cell;
    result.size = 1;
    if (upper_cell != lower_cell) {
        result.values[1] = upper_cell;
        result.size = 2;
    }
    return result;
}

struct ComponentBest {
    Length score{};
    std::size_t end{};
    bool available{};
};

struct ComponentQueueEntry {
    Length score{};
    std::size_t end{};
    std::size_t component{};
    std::uint64_t generation{};
};

struct WorseComponentQueueEntry {
    [[nodiscard]] bool operator()(const ComponentQueueEntry& left,
                                  const ComponentQueueEntry& right) const noexcept {
        if (left.score != right.score) {
            return left.score < right.score;
        }
        if (left.end != right.end) {
            return left.end > right.end;
        }
        return left.component > right.component;
    }
};

[[nodiscard]] bool BetterComponentBest(const ComponentBest& left,
                                       const ComponentBest& right) noexcept {
    if (!left.available) {
        return false;
    }
    if (!right.available) {
        return true;
    }
    return left.score > right.score ||
           (left.score == right.score && left.end < right.end);
}

void AddSparseGroupChains(const std::vector<OrientedSeed>& oriented,
                          std::size_t group_begin,
                          std::size_t group_end,
                          const AlignmentOptions& options,
                          SharedWorkingSetBudget& working_budget,
                          SharedCountBudget& candidate_budget,
                          SharedCountBudget& legal_edge_budget,
                          SharedCountBudget& relaxation_budget,
                          const std::atomic<bool>* cancellation,
                          RunStatistics& statistics,
                          std::uint64_t& temporary_working_bytes,
                          std::uint64_t& retained_chain_bytes,
                          std::vector<Chain>& chains) {
    const std::size_t group_size = group_end - group_begin;
    if (group_size == 0) {
        return;
    }
    const std::string group = ChainingGroupDescription(oriented[group_begin]);
    const Length gap_bound = std::min(options.max_gap, options.break_length);
    GroupWorkingReservation working{working_budget, group};

    // Conservative deterministic payload plus hash-node overhead.  Common
    // oriented/input storage is reserved once by the caller; this reservation
    // covers only the independently executing group.
    constexpr std::size_t kPerSeedFixedBytes =
        sizeof(std::size_t) * 10 + sizeof(Length) + sizeof(bool) +
        sizeof(std::uint8_t) + sizeof(ComponentBest) + sizeof(EndpointCell) +
        sizeof(std::vector<std::size_t>) + sizeof(void*) * 4;
    working.AddTemporary(CheckedByteProduct(group_size, kPerSeedFixedBytes));
    working.AddTemporary(sizeof(std::size_t) * 2);

    std::vector<std::size_t> incoming_offsets(group_size + 1, 0);
    std::vector<std::size_t> incoming_predecessors;
    const auto ensure_incoming_capacity = [&] {
        if (incoming_predecessors.size() < incoming_predecessors.capacity()) {
            return;
        }
        const std::size_t old_capacity = incoming_predecessors.capacity();
        constexpr std::size_t kInitialEdgeCapacity = 4096;
        std::size_t new_capacity = 0;
        if (old_capacity == 0) {
            new_capacity = std::min(group_size, kInitialEdgeCapacity);
            new_capacity = std::max<std::size_t>(new_capacity, 1);
        } else {
            if (old_capacity > std::numeric_limits<std::size_t>::max() / 2) {
                throw AlignmentError(
                    "chaining predecessor capacity exceeds size_t range");
            }
            new_capacity = old_capacity * 2;
        }
        const std::uint64_t old_bytes =
            CheckedByteProduct(old_capacity, sizeof(std::size_t));
        const std::uint64_t new_bytes =
            CheckedByteProduct(new_capacity, sizeof(std::size_t));

        // reserve() may briefly hold the old and new buffers together.  Add
        // the full new allocation first, then release the old accounting so
        // the shared peak remains a conservative bound without one CAS per
        // legal edge.
        working.AddTemporary(new_bytes);
        incoming_predecessors.reserve(new_capacity);
        working.ReleaseTemporary(old_bytes);
    };
    DisjointSets components(group_size);
    std::unordered_map<EndpointCell, std::vector<std::size_t>, EndpointCellHash>
        endpoint_buckets;
    endpoint_buckets.reserve(group_size);

    for (std::size_t local_i = 0; local_i < group_size; ++local_i) {
        CheckChainingCancellation(cancellation);
        const OrientedSeed& current = oriented[group_begin + local_i];
        const CellRange reference_cells =
            EndpointSearchCells(current.seed.reference_begin, gap_bound);
        const CellRange query_cells =
            EndpointSearchCells(current.query_begin, gap_bound);
        for (std::size_t ref_cell_index = 0;
             ref_cell_index < reference_cells.size; ++ref_cell_index) {
            for (std::size_t query_cell_index = 0;
                 query_cell_index < query_cells.size; ++query_cell_index) {
                const EndpointCell key{reference_cells.values[ref_cell_index],
                                       query_cells.values[query_cell_index]};
                const auto bucket = endpoint_buckets.find(key);
                if (bucket == endpoint_buckets.end()) {
                    continue;
                }
                for (const std::size_t local_j : bucket->second) {
                    IncrementGloballyLimited(
                        statistics.chaining_candidate_pairs, candidate_budget,
                        "candidate pair checks", group);
                    const OrientedSeed& previous = oriented[group_begin + local_j];
                    if (!CanChain(previous, current, options)) {
                        continue;
                    }
                    IncrementGloballyLimited(statistics.chaining_legal_edges,
                                             legal_edge_budget,
                                             "legal edges", group);
                    ensure_incoming_capacity();
                    incoming_predecessors.push_back(local_j);
                    components.Unite(local_j, local_i);
                }
            }
        }
        incoming_offsets[local_i + 1] = incoming_predecessors.size();

        const Position reference_end =
            CheckedEnd(current.seed.reference_begin, current.seed.length,
                       "reference seed endpoint for sparse chaining");
        const Position query_end =
            CheckedEnd(current.query_begin, current.seed.length,
                       "query seed endpoint for sparse chaining");
        const EndpointCell endpoint{CellCoordinate(reference_end, gap_bound),
                                    CellCoordinate(query_end, gap_bound)};
        auto [bucket, inserted] = endpoint_buckets.try_emplace(endpoint);
        static_cast<void>(inserted);
        bucket->second.push_back(local_i);
        statistics.chaining_max_bucket_occupancy = std::max(
            statistics.chaining_max_bucket_occupancy,
            static_cast<std::uint64_t>(bucket->second.size()));
    }

    const std::size_t no_component = group_size;
    std::vector<std::size_t> component_of_root(group_size, no_component);
    std::vector<std::size_t> component_of_node(group_size, no_component);
    std::size_t component_count = 0;
    for (std::size_t node = 0; node < group_size; ++node) {
        const std::size_t root = components.Find(node);
        if (component_of_root[root] == no_component) {
            component_of_root[root] = component_count;
            ++component_count;
        }
        component_of_node[node] = component_of_root[root];
    }
    if (static_cast<std::uint64_t>(component_count) >
        std::numeric_limits<std::uint64_t>::max() - statistics.chaining_components) {
        throw AlignmentError("chaining component count exceeds 64-bit range");
    }
    statistics.chaining_components += static_cast<std::uint64_t>(component_count);

    std::vector<std::size_t> component_counts(component_count, 0);
    for (const std::size_t component : component_of_node) {
        ++component_counts[component];
    }
    std::vector<std::size_t> component_offsets(component_count + 1, 0);
    for (std::size_t component = 0; component < component_count; ++component) {
        if (component_counts[component] >
            std::numeric_limits<std::size_t>::max() - component_offsets[component]) {
            throw AlignmentError("chaining component offsets exceed size_t range");
        }
        component_offsets[component + 1] =
            component_offsets[component] + component_counts[component];
        statistics.chaining_max_component_seeds = std::max(
            statistics.chaining_max_component_seeds,
            static_cast<std::uint64_t>(component_counts[component]));
    }
    std::vector<std::size_t> component_nodes(group_size, 0);
    std::vector<std::size_t> component_cursor = component_offsets;
    for (std::size_t node = 0; node < group_size; ++node) {
        const std::size_t component = component_of_node[node];
        component_nodes[component_cursor[component]++] = node;
    }
    std::vector<std::size_t> component_edge_counts(component_count, 0);
    for (std::size_t node = 0; node < group_size; ++node) {
        const std::size_t edges = incoming_offsets[node + 1] - incoming_offsets[node];
        const std::size_t component = component_of_node[node];
        if (edges > std::numeric_limits<std::size_t>::max() -
                        component_edge_counts[component]) {
            throw AlignmentError("component edge count exceeds size_t range");
        }
        component_edge_counts[component] += edges;
    }
    for (const std::size_t edge_count : component_edge_counts) {
        statistics.chaining_max_component_edges = std::max(
            statistics.chaining_max_component_edges,
            static_cast<std::uint64_t>(edge_count));
    }

    std::vector<bool> used(group_size, false);
    std::vector<Length> best_score(group_size, 0);
    std::vector<std::size_t> predecessor(group_size, group_size);
    std::vector<ComponentBest> component_best(component_count);
    std::vector<std::uint64_t> component_generation(component_count, 0);

    const auto recompute_component = [&](std::size_t component) {
        CheckChainingCancellation(cancellation);
        IncrementUnbounded(statistics.chaining_dp_passes, "chaining DP pass count");
        ComponentBest best;
        best.end = group_size;
        for (std::size_t cursor = component_offsets[component];
             cursor < component_offsets[component + 1]; ++cursor) {
            const std::size_t local_i = component_nodes[cursor];
            if (used[local_i]) {
                best_score[local_i] = 0;
                predecessor[local_i] = group_size;
                continue;
            }
            const OrientedSeed& current = oriented[group_begin + local_i];
            best_score[local_i] = current.seed.length;
            predecessor[local_i] = group_size;
            for (std::size_t edge = incoming_offsets[local_i];
                 edge < incoming_offsets[local_i + 1]; ++edge) {
                const std::size_t local_j = incoming_predecessors[edge];
                if (used[local_j] || best_score[local_j] == 0) {
                    continue;
                }
                IncrementGloballyLimited(
                    statistics.chaining_edge_relaxations, relaxation_budget,
                    "edge relaxations", group);
                const Length candidate =
                    CheckedLengthAdd(best_score[local_j], current.seed.length,
                                     "chain seed length");
                if (candidate > best_score[local_i] ||
                    (candidate == best_score[local_i] &&
                     local_j < predecessor[local_i])) {
                    best_score[local_i] = candidate;
                    predecessor[local_i] = local_j;
                }
            }
            const ComponentBest candidate{best_score[local_i], local_i, true};
            if (BetterComponentBest(candidate, best)) {
                best = candidate;
            }
        }
        component_best[component] = best;
    };

    for (std::size_t component = 0; component < component_count; ++component) {
        recompute_component(component);
    }

    working.AddTemporary(
        CheckedByteProduct(component_count, sizeof(ComponentQueueEntry)));
    std::vector<ComponentQueueEntry> initial_queue;
    initial_queue.reserve(component_count);
    for (std::size_t component = 0; component < component_count; ++component) {
        if (component_best[component].available) {
            initial_queue.push_back(ComponentQueueEntry{
                component_best[component].score,
                component_best[component].end,
                component,
                component_generation[component]});
        }
    }
    std::priority_queue<ComponentQueueEntry,
                        std::vector<ComponentQueueEntry>,
                        WorseComponentQueueEntry>
        ready{WorseComponentQueueEntry{}, std::move(initial_queue)};

    while (!ready.empty()) {
        CheckChainingCancellation(cancellation);
        const ComponentQueueEntry entry = ready.top();
        ready.pop();
        if (entry.generation != component_generation[entry.component]) {
            continue;
        }
        if (entry.score < options.min_cluster) {
            break;
        }
        const std::size_t winner = entry.component;
        const ComponentBest winning_best = component_best[winner];

        std::size_t path_length = 0;
        for (std::size_t cursor = winning_best.end; cursor != group_size;
             cursor = predecessor[cursor]) {
            ++path_length;
            if (path_length > component_counts[winner]) {
                throw AlignmentError("internal error: sparse chain predecessor cycle");
            }
        }
        const std::uint64_t retained_seed_bytes =
            CheckedByteProduct(path_length, sizeof(OrientedSeed));
        working.AddRetained(
            CheckedByteAdd(retained_seed_bytes, sizeof(Chain)));

        std::vector<std::size_t> chain_indices;
        chain_indices.reserve(path_length);
        for (std::size_t cursor = winning_best.end; cursor != group_size;
             cursor = predecessor[cursor]) {
            chain_indices.push_back(cursor);
        }
        std::reverse(chain_indices.begin(), chain_indices.end());

        Chain chain;
        chain.reserve(chain_indices.size());
        for (const std::size_t local_index : chain_indices) {
            used[local_index] = true;
            chain.push_back(oriented[group_begin + local_index]);
        }
        chains.push_back(std::move(chain));
        recompute_component(winner);
        IncrementUnbounded(component_generation[winner],
                           "chaining component generation");
        if (component_best[winner].available) {
            ready.push(ComponentQueueEntry{
                component_best[winner].score,
                component_best[winner].end,
                winner,
                component_generation[winner]});
        }
    }
    temporary_working_bytes = working.temporary_bytes();
    retained_chain_bytes = working.retained_bytes();
    working.CommitRetained();
}

struct SeedGroup {
    std::size_t begin{};
    std::size_t end{};
};

struct GroupChainingResult {
    std::vector<Chain> chains;
    RunStatistics statistics;
    std::uint64_t temporary_working_bytes{};
    std::uint64_t retained_chain_bytes{};
};

struct alignas(64) WorkerActivity {
    bool processed{};
};

[[nodiscard]] std::uint32_t PlannedWorkerCount(std::uint32_t requested,
                                               std::size_t task_count) {
    if (task_count == 0) {
        return 1;
    }
    const std::uint64_t bounded = std::min<std::uint64_t>(
        requested, static_cast<std::uint64_t>(task_count));
#if defined(_OPENMP)
    return static_cast<std::uint32_t>(std::min<std::uint64_t>(
        bounded, static_cast<std::uint64_t>(std::numeric_limits<int>::max())));
#else
    static_cast<void>(bounded);
    return 1;
#endif
}

void CheckedAggregate(std::uint64_t& total,
                      std::uint64_t value,
                      std::string_view resource) {
    if (value > std::numeric_limits<std::uint64_t>::max() - total) {
        throw AlignmentError(std::string{resource} +
                             " exceeds the 64-bit counter range");
    }
    total += value;
}

void CheckedAggregateLimited(std::uint64_t& total,
                             std::uint64_t value,
                             std::uint64_t limit,
                             std::string_view resource,
                             std::string_view group) {
    if (value > std::numeric_limits<std::uint64_t>::max() - total) {
        ThrowChainingLimit(resource, std::numeric_limits<std::uint64_t>::max(),
                           limit, group);
    }
    const std::uint64_t observed = total + value;
    if (observed > limit) {
        ThrowChainingLimit(resource, observed, limit, group);
    }
    total = observed;
}

[[nodiscard]] std::uint64_t DeterministicWorkingSetBound(
    std::uint64_t common_bytes,
    const std::vector<GroupChainingResult>& groups,
    std::uint32_t worker_count) {
    std::uint64_t result = common_bytes;
    std::vector<std::uint64_t> temporary;
    temporary.reserve(groups.size());
    for (const GroupChainingResult& group : groups) {
        result = CheckedByteAdd(result, group.retained_chain_bytes);
        temporary.push_back(group.temporary_working_bytes);
    }
    std::sort(temporary.begin(), temporary.end(), std::greater<>{});
    const std::size_t simultaneous = std::min<std::size_t>(
        temporary.size(), static_cast<std::size_t>(worker_count));
    for (std::size_t index = 0; index < simultaneous; ++index) {
        result = CheckedByteAdd(result, temporary[index]);
    }
    return result;
}

[[nodiscard]] std::vector<Chain> BuildChainsSparseExactEdgeComponents(
    const std::vector<Seed>& seeds,
    const QueryLengths& query_lengths,
    const AlignmentOptions& options,
    RunStatistics& statistics) {
    statistics.chaining_route = "sparse-exact-edge-components-v1";
    statistics.chaining_candidate_pairs = 0;
    statistics.chaining_legal_edges = 0;
    statistics.chaining_components = 0;
    statistics.chaining_max_bucket_occupancy = 0;
    statistics.chaining_max_component_seeds = 0;
    statistics.chaining_max_component_edges = 0;
    statistics.chaining_dp_passes = 0;
    statistics.chaining_edge_relaxations = 0;
    statistics.chaining_working_set_peak_bytes = 0;
    statistics.chaining_working_set_limit_bytes =
        options.chaining_limits.working_set_bytes;
    statistics.chaining_candidate_pair_limit =
        options.chaining_limits.candidate_pairs;
    statistics.chaining_legal_edge_limit = options.chaining_limits.legal_edges;
    statistics.chaining_edge_relaxation_limit =
        options.chaining_limits.edge_relaxations;
    statistics.chaining_requested_threads = options.worker_threads;
    statistics.chaining_worker_threads = 1;

    const std::vector<OrientedSeed> oriented =
        MakeSortedOrientedSeeds(seeds, query_lengths);
    std::vector<SeedGroup> groups;
    std::size_t group_begin = 0;
    while (group_begin < oriented.size()) {
        std::size_t group_end = group_begin + 1;
        while (group_end < oriented.size() &&
               SameSeedGroup(oriented[group_begin], oriented[group_end])) {
            ++group_end;
        }
        groups.push_back(SeedGroup{group_begin, group_end});
        group_begin = group_end;
    }
    if (groups.size() >
        static_cast<std::size_t>(std::numeric_limits<std::ptrdiff_t>::max())) {
        throw AlignmentError("chaining group count exceeds the OpenMP loop range");
    }

    std::uint64_t common_working_bytes =
        CheckedByteProduct(oriented.size(), sizeof(OrientedSeed));
    common_working_bytes = CheckedByteAdd(
        common_working_bytes, CheckedByteProduct(seeds.size(), sizeof(Seed)));
    common_working_bytes = CheckedByteAdd(
        common_working_bytes, CheckedByteProduct(groups.size(), sizeof(SeedGroup)));
    common_working_bytes = CheckedByteAdd(
        common_working_bytes,
        CheckedByteProduct(groups.size(), sizeof(GroupChainingResult)));
    common_working_bytes = CheckedByteAdd(
        common_working_bytes,
        CheckedByteProduct(groups.size(), sizeof(std::exception_ptr)));

    auto working_budget = std::make_unique<SharedWorkingSetBudget>(
        options.chaining_limits.working_set_bytes);
    auto candidate_budget = std::make_unique<SharedCountBudget>(
        options.chaining_limits.candidate_pairs, "candidate pair checks");
    auto legal_edge_budget = std::make_unique<SharedCountBudget>(
        options.chaining_limits.legal_edges, "legal edges");
    auto relaxation_budget = std::make_unique<SharedCountBudget>(
        options.chaining_limits.edge_relaxations, "edge relaxations");
    working_budget->Acquire(common_working_bytes, "all chaining groups");

    std::vector<GroupChainingResult> group_results(groups.size());
    std::vector<std::exception_ptr> group_errors(groups.size());
    std::atomic<bool> cancellation_requested{false};
    const std::uint32_t planned_workers =
        PlannedWorkerCount(options.worker_threads, groups.size());
    std::vector<WorkerActivity> chaining_activity(
        static_cast<std::size_t>(planned_workers));
    std::uint32_t actual_workers = 1;

    const auto process_group = [&](std::size_t group_index) {
        if (cancellation_requested.load(std::memory_order_relaxed)) {
            return;
        }
        try {
            if (options.interruption_callback) {
                options.interruption_callback("chaining");
            }
            const SeedGroup group = groups[group_index];
            GroupChainingResult local;
            AddSparseGroupChains(oriented, group.begin, group.end, options,
                                 *working_budget, *candidate_budget,
                                 *legal_edge_budget, *relaxation_budget,
                                 &cancellation_requested,
                                 local.statistics,
                                 local.temporary_working_bytes,
                                 local.retained_chain_bytes, local.chains);
            group_results[group_index] = std::move(local);
        } catch (const ChainingCancelled&) {
            return;
        } catch (...) {
            group_errors[group_index] = std::current_exception();
            cancellation_requested.store(true, std::memory_order_relaxed);
        }
    };

#if defined(_OPENMP)
    if (!groups.empty()) {
        const int omp_workers = static_cast<int>(planned_workers);
#pragma omp parallel num_threads(omp_workers)
        {
#pragma omp for schedule(dynamic, 1)
            for (std::ptrdiff_t signed_index = 0;
                 signed_index < static_cast<std::ptrdiff_t>(groups.size());
                 ++signed_index) {
                chaining_activity[static_cast<std::size_t>(
                    omp_get_thread_num())]
                    .processed = true;
                process_group(static_cast<std::size_t>(signed_index));
            }
        }
    }
#else
    static_cast<void>(planned_workers);
    for (std::size_t group_index = 0; group_index < groups.size(); ++group_index) {
        chaining_activity[0].processed = true;
        process_group(group_index);
    }
#endif
    actual_workers = 0;
    for (const WorkerActivity& activity : chaining_activity) {
        if (activity.processed) {
            ++actual_workers;
        }
    }
    actual_workers = std::max<std::uint32_t>(actual_workers, 1);
    statistics.chaining_worker_threads = actual_workers;

    for (std::size_t group_index = 0; group_index < group_errors.size();
         ++group_index) {
        if (group_errors[group_index] != nullptr) {
            std::rethrow_exception(group_errors[group_index]);
        }
    }

    std::size_t chain_count = 0;
    for (std::size_t group_index = 0; group_index < group_results.size();
         ++group_index) {
        const GroupChainingResult& group = group_results[group_index];
        const std::string group_description =
            ChainingGroupDescription(oriented[groups[group_index].begin]);
        CheckedAggregateLimited(statistics.chaining_candidate_pairs,
                                group.statistics.chaining_candidate_pairs,
                                options.chaining_limits.candidate_pairs,
                                "candidate pair checks", group_description);
        CheckedAggregateLimited(statistics.chaining_legal_edges,
                                group.statistics.chaining_legal_edges,
                                options.chaining_limits.legal_edges,
                                "legal edges", group_description);
        CheckedAggregateLimited(statistics.chaining_edge_relaxations,
                                group.statistics.chaining_edge_relaxations,
                                options.chaining_limits.edge_relaxations,
                                "edge relaxations", group_description);
        CheckedAggregate(statistics.chaining_components,
                         group.statistics.chaining_components,
                         "chaining component count");
        CheckedAggregate(statistics.chaining_dp_passes,
                         group.statistics.chaining_dp_passes,
                         "chaining DP pass count");
        statistics.chaining_max_bucket_occupancy = std::max(
            statistics.chaining_max_bucket_occupancy,
            group.statistics.chaining_max_bucket_occupancy);
        statistics.chaining_max_component_seeds = std::max(
            statistics.chaining_max_component_seeds,
            group.statistics.chaining_max_component_seeds);
        statistics.chaining_max_component_edges = std::max(
            statistics.chaining_max_component_edges,
            group.statistics.chaining_max_component_edges);
        if (group.chains.size() >
            std::numeric_limits<std::size_t>::max() - chain_count) {
            throw AlignmentError("chain count exceeds size_t range");
        }
        chain_count += group.chains.size();
    }

    if (statistics.chaining_candidate_pairs != candidate_budget->value() ||
        statistics.chaining_legal_edges != legal_edge_budget->value() ||
        statistics.chaining_edge_relaxations != relaxation_budget->value()) {
        throw AlignmentError(
            "internal error: global chaining work budgets differ from "
            "deterministic group aggregation");
    }

    const std::uint64_t deterministic_working_set =
        DeterministicWorkingSetBound(common_working_bytes, group_results,
                                     planned_workers);
    if (working_budget->peak() > deterministic_working_set) {
        throw AlignmentError(
            "internal error: concurrent chaining exceeded its deterministic "
            "working-set bound");
    }
    if (deterministic_working_set > options.chaining_limits.working_set_bytes) {
        ThrowChainingLimit("estimated working-set bytes",
                           deterministic_working_set,
                           options.chaining_limits.working_set_bytes,
                           "all chaining groups");
    }
    statistics.chaining_working_set_peak_bytes = deterministic_working_set;

    std::vector<Chain> chains;
    chains.reserve(chain_count);
    for (GroupChainingResult& group : group_results) {
        for (Chain& chain : group.chains) {
            chains.push_back(std::move(chain));
        }
    }
    return chains;
}

void AppendCigar(std::vector<CigarOp>& cigar, char operation, Length length) {
    if (length == 0) {
        return;
    }
    if (operation != '=' && operation != 'X' && operation != 'I' && operation != 'D') {
        throw AlignmentError("internal error: unsupported CIGAR operation");
    }
    if (!cigar.empty() && cigar.back().operation == operation) {
        cigar.back().length =
            CheckedLengthAdd(cigar.back().length, length, "CIGAR operation length");
    } else {
        cigar.push_back(CigarOp{operation, length});
    }
}

using internal::ExtensionCallMetrics;
using internal::GapAlignment;

[[nodiscard]] std::int64_t GapPenalty(Length length,
                                      const AlignmentOptions& options) {
    const std::int64_t extension =
        CheckedScoreProduct(length, options.gap_extend_penalty, "gap extension score");
    const std::int64_t total = CheckedScoreAdd(
        static_cast<std::int64_t>(options.gap_open_penalty), extension,
        "gap penalty");
    return -total;
}

[[nodiscard]] GapAlignment AlignScalarAffineGap(std::string_view reference,
                                                std::string_view query,
                                                const AlignmentOptions& options,
                                                bool force_affine = false) {
    if (reference.empty()) {
        GapAlignment result;
        AppendCigar(result.cigar, 'I', static_cast<Length>(query.size()));
        result.score = GapPenalty(static_cast<Length>(query.size()), options);
        return result;
    }
    if (query.empty()) {
        GapAlignment result;
        AppendCigar(result.cigar, 'D', static_cast<Length>(reference.size()));
        result.score = GapPenalty(static_cast<Length>(reference.size()), options);
        return result;
    }

    if (reference.size() == query.size() && !force_affine) {
        GapAlignment result;
        bool all_exact = true;
        for (std::size_t index = 0; index < reference.size(); ++index) {
            if (!EqualCanonical(reference[index], query[index])) {
                all_exact = false;
                break;
            }
        }
        if (all_exact) {
            AppendCigar(result.cigar, '=', static_cast<Length>(reference.size()));
            result.score = CheckedScoreProduct(
                static_cast<Length>(reference.size()), options.match_score,
                "exact gap score");
            return result;
        }
        for (std::size_t index = 0; index < reference.size(); ++index) {
            const bool match = EqualCanonical(reference[index], query[index]);
            AppendCigar(result.cigar, match ? '=' : 'X', 1);
            result.score = CheckedScoreAdd(
                result.score,
                match ? static_cast<std::int64_t>(options.match_score)
                      : -static_cast<std::int64_t>(options.mismatch_penalty),
                "ungapped alignment score");
        }
        return result;
    }

    if (reference.size() == std::numeric_limits<std::size_t>::max() ||
        query.size() == std::numeric_limits<std::size_t>::max()) {
        throw AlignmentError("bounded gap DP dimensions overflow size_t");
    }
    const std::size_t rows = reference.size() + 1;
    const std::size_t columns = query.size() + 1;
    if (rows > std::numeric_limits<std::size_t>::max() / columns ||
        rows > options.max_dp_cells / columns) {
        throw AlignmentError("bounded gap DP limit exceeded; the chain should have been split");
    }
    const std::size_t cells = rows * columns;
    constexpr std::int64_t kNegativeInfinity =
        std::numeric_limits<std::int64_t>::lowest() / 4;
    constexpr std::uint8_t kMatchState = 0;
    constexpr std::uint8_t kInsertionState = 1;
    constexpr std::uint8_t kDeletionState = 2;
    constexpr std::uint8_t kNoState = 3;

    std::vector<std::int64_t> match(cells, kNegativeInfinity);
    std::vector<std::int64_t> insertion(cells, kNegativeInfinity);
    std::vector<std::int64_t> deletion(cells, kNegativeInfinity);
    std::vector<std::uint8_t> trace_match(cells, kNoState);
    std::vector<std::uint8_t> trace_insertion(cells, kNoState);
    std::vector<std::uint8_t> trace_deletion(cells, kNoState);

    const auto index_of = [columns](std::size_t row, std::size_t column) {
        return row * columns + column;
    };
    const std::int64_t open_and_extend =
        static_cast<std::int64_t>(options.gap_open_penalty) +
        static_cast<std::int64_t>(options.gap_extend_penalty);
    const std::int64_t extend = options.gap_extend_penalty;
    const auto subtract = [&](std::int64_t score, std::int64_t penalty) {
        return score == kNegativeInfinity ? kNegativeInfinity : score - penalty;
    };
    const auto choose = [=](std::int64_t first,
                            std::int64_t second,
                            std::int64_t third) -> std::pair<std::int64_t, std::uint8_t> {
        // Stable tie order is match, insertion, deletion.
        if (first >= second && first >= third) {
            return {first, kMatchState};
        }
        if (second >= third) {
            return {second, kInsertionState};
        }
        return {third, kDeletionState};
    };

    match[0] = 0;
    for (std::size_t column = 1; column < columns; ++column) {
        const std::size_t cell = index_of(0, column);
        insertion[cell] = GapPenalty(static_cast<Length>(column), options);
        trace_insertion[cell] = column == 1 ? kMatchState : kInsertionState;
    }
    for (std::size_t row = 1; row < rows; ++row) {
        const std::size_t cell = index_of(row, 0);
        deletion[cell] = GapPenalty(static_cast<Length>(row), options);
        trace_deletion[cell] = row == 1 ? kMatchState : kDeletionState;
    }

    for (std::size_t row = 1; row < rows; ++row) {
        for (std::size_t column = 1; column < columns; ++column) {
            const std::size_t cell = index_of(row, column);
            const std::size_t diagonal = index_of(row - 1, column - 1);
            const auto best_diagonal =
                choose(match[diagonal], insertion[diagonal], deletion[diagonal]);
            const std::int64_t substitution =
                EqualCanonical(reference[row - 1], query[column - 1])
                    ? static_cast<std::int64_t>(options.match_score)
                    : -static_cast<std::int64_t>(options.mismatch_penalty);
            match[cell] = CheckedScoreAdd(best_diagonal.first, substitution,
                                          "affine alignment score");
            trace_match[cell] = best_diagonal.second;

            const std::size_t left = index_of(row, column - 1);
            const auto best_insertion = choose(subtract(match[left], open_and_extend),
                                               subtract(insertion[left], extend),
                                               subtract(deletion[left], open_and_extend));
            insertion[cell] = best_insertion.first;
            trace_insertion[cell] = best_insertion.second;

            const std::size_t above = index_of(row - 1, column);
            const auto best_deletion = choose(subtract(match[above], open_and_extend),
                                              subtract(insertion[above], open_and_extend),
                                              subtract(deletion[above], extend));
            deletion[cell] = best_deletion.first;
            trace_deletion[cell] = best_deletion.second;
        }
    }

    const std::size_t final_cell = index_of(reference.size(), query.size());
    const auto final = choose(match[final_cell], insertion[final_cell], deletion[final_cell]);
    GapAlignment result;
    result.score = final.first;
    std::size_t row = reference.size();
    std::size_t column = query.size();
    std::uint8_t state = final.second;
    std::vector<CigarOp> reverse_cigar;
    while (row != 0 || column != 0) {
        const std::size_t cell = index_of(row, column);
        if (state == kMatchState) {
            if (row == 0 || column == 0) {
                throw AlignmentError("internal error while tracing affine match state");
            }
            const bool is_match = EqualCanonical(reference[row - 1], query[column - 1]);
            AppendCigar(reverse_cigar, is_match ? '=' : 'X', 1);
            state = trace_match[cell];
            --row;
            --column;
        } else if (state == kInsertionState) {
            if (column == 0) {
                throw AlignmentError("internal error while tracing affine insertion state");
            }
            AppendCigar(reverse_cigar, 'I', 1);
            state = trace_insertion[cell];
            --column;
        } else if (state == kDeletionState) {
            if (row == 0) {
                throw AlignmentError("internal error while tracing affine deletion state");
            }
            AppendCigar(reverse_cigar, 'D', 1);
            state = trace_deletion[cell];
            --row;
        } else {
            throw AlignmentError("internal error: affine traceback reached no-state");
        }
    }
    std::reverse(reverse_cigar.begin(), reverse_cigar.end());
    for (const CigarOp& operation : reverse_cigar) {
        AppendCigar(result.cigar, operation.operation, operation.length);
    }
    return result;
}

[[nodiscard]] std::int64_t ValidateAndScoreGapAlignment(
    std::string_view reference,
    std::string_view query,
    const AlignmentOptions& options,
    const GapAlignment& alignment) {
    Position reference_offset = 0;
    Position query_offset = 0;
    std::int64_t score = 0;
    char previous_operation = '\0';
    for (const CigarOp& operation : alignment.cigar) {
        if (operation.length == 0) {
            throw AlignmentError("extension backend returned a zero-length CIGAR operation");
        }
        if (operation.operation == previous_operation) {
            throw AlignmentError("extension backend returned a non-canonical CIGAR");
        }
        previous_operation = operation.operation;
        if (operation.operation == '=' || operation.operation == 'X') {
            const Position reference_end = CheckedEnd(
                reference_offset, operation.length, "extension CIGAR reference");
            const Position query_end = CheckedEnd(
                query_offset, operation.length, "extension CIGAR query");
            if (reference_end > reference.size() || query_end > query.size()) {
                throw AlignmentError("extension backend CIGAR exceeds its input sequences");
            }
            for (Length offset = 0; offset < operation.length; ++offset) {
                const bool equal = EqualCanonical(
                    reference[static_cast<std::size_t>(reference_offset + offset)],
                    query[static_cast<std::size_t>(query_offset + offset)]);
                if ((operation.operation == '=') != equal) {
                    throw AlignmentError(
                        "extension backend CIGAR disagrees with canonical bases");
                }
                score = CheckedScoreAdd(
                    score,
                    equal ? static_cast<std::int64_t>(options.match_score)
                          : -static_cast<std::int64_t>(options.mismatch_penalty),
                    "extension backend score");
            }
            reference_offset = reference_end;
            query_offset = query_end;
        } else if (operation.operation == 'I') {
            query_offset = CheckedEnd(
                query_offset, operation.length, "extension CIGAR query insertion");
            if (query_offset > query.size()) {
                throw AlignmentError("extension backend insertion exceeds the query");
            }
            score = CheckedScoreAdd(
                score, GapPenalty(operation.length, options),
                "extension backend insertion score");
        } else if (operation.operation == 'D') {
            reference_offset = CheckedEnd(
                reference_offset, operation.length, "extension CIGAR reference deletion");
            if (reference_offset > reference.size()) {
                throw AlignmentError("extension backend deletion exceeds the reference");
            }
            score = CheckedScoreAdd(
                score, GapPenalty(operation.length, options),
                "extension backend deletion score");
        } else {
            throw AlignmentError("extension backend returned an unsupported CIGAR operation");
        }
    }
    if (reference_offset != reference.size() || query_offset != query.size()) {
        throw AlignmentError("extension backend CIGAR does not consume both sequences");
    }
    if (score != alignment.score) {
        throw AlignmentError(
            "extension backend score disagrees with its CIGAR: reported=" +
            std::to_string(alignment.score) + ", recomputed=" +
            std::to_string(score) + ", cigar=" + CigarToString(alignment.cigar));
    }
    return score;
}

[[nodiscard]] GapAlignment AlignConfiguredUnequalGap(
    std::string_view reference,
    std::string_view query,
    const AlignmentOptions& options,
    ExtensionCallMetrics* metrics) {
    if (reference.empty() || query.empty() || reference.size() == query.size()) {
        throw AlignmentError(
            "internal error: experimental backend requires non-empty unequal gaps");
    }
#if defined(RAMAG_EXTENSION_KSW2_EXACT)
    GapAlignment result = internal::AlignKsw2Gap(
        reference, query, options, false, metrics);
#elif defined(RAMAG_EXTENSION_KSW2_BAND_AUTO)
    GapAlignment result = internal::AlignKsw2Gap(
        reference, query, options, true, metrics);
#elif defined(RAMAG_EXTENSION_BLOCK_EXACT)
    GapAlignment result = internal::AlignBlockGap(
        reference, query, options, true, metrics);
#elif defined(RAMAG_EXTENSION_BLOCK_ADAPTIVE)
    GapAlignment result = internal::AlignBlockGap(
        reference, query, options, false, metrics);
#else
    GapAlignment result = AlignScalarAffineGap(reference, query, options, true);
    if (metrics != nullptr) {
        const std::uint64_t rows = static_cast<std::uint64_t>(reference.size()) + 1;
        const std::uint64_t columns = static_cast<std::uint64_t>(query.size()) + 1;
        metrics->full_matrix_cells = rows * columns;
        metrics->estimated_cells = metrics->full_matrix_cells;
    }
#endif
    static_cast<void>(ValidateAndScoreGapAlignment(
        reference, query, options, result));
    return result;
}

[[nodiscard]] GapAlignment AlignAffineGap(std::string_view reference,
                                          std::string_view query,
                                          const AlignmentOptions& options,
                                          bool force_affine,
                                          ExtensionCallMetrics* metrics) {
#if defined(RAMAG_EXTENSION_SCALAR) && !RAMAG_EXTENSION_METRICS
    static_cast<void>(metrics);
    return AlignScalarAffineGap(reference, query, options, force_affine);
#else
    if (reference.empty() || query.empty() || reference.size() == query.size()) {
        return AlignScalarAffineGap(reference, query, options, force_affine);
    }
    return AlignConfiguredUnequalGap(reference, query, options, metrics);
#endif
}

[[nodiscard]] bool ShouldRealignEqualGap(std::string_view reference,
                                        std::string_view query,
                                        const AlignmentOptions& options) noexcept {
    // A pair of compensating indels can leave equal-length seed gaps while
    // making the direct diagonal substantially worse than an affine path.
    // Derive a conservative trigger from the configured scores: two
    // compensating one-base gaps pay two open-and-extend penalties.  Count
    // only newly recovered match score as evidence, rather than assuming that
    // every diagonal mismatch penalty will disappear; the latter is too
    // eager in repetitive sequence.  Require a strict improvement before
    // entering the expensive DP.
    constexpr std::size_t kMaxEqualGapDpLength = 128;
    if (reference.size() != query.size() ||
        reference.size() > kMaxEqualGapDpLength) {
        return false;
    }
    const std::uint64_t gap_cost =
        2ULL * (static_cast<std::uint64_t>(options.gap_open_penalty) +
                static_cast<std::uint64_t>(options.gap_extend_penalty));
    const std::uint64_t recovered_match_gain =
        static_cast<std::uint64_t>(options.match_score);
    if (recovered_match_gain == 0) {
        return false;
    }
    const std::uint64_t minimum_mismatches =
        gap_cost / recovered_match_gain + 1ULL;
    std::size_t mismatches = 0;
    for (std::size_t index = 0; index < reference.size(); ++index) {
        if (!EqualCanonical(reference[index], query[index]) &&
            ++mismatches >= minimum_mismatches) {
            return true;
        }
    }
    return false;
}

using SequenceLookup = std::unordered_map<SequenceId, const SequenceRecord*>;

// Gap counters are the only per-extension statistics.  Keeping one padded
// accumulator per worker avoids constructing a full RunStatistics object for
// every chain and prevents adjacent workers from contending on one cache line.
struct alignas(64) GapCountAccumulator {
    std::uint64_t exact{};
    std::uint64_t ungapped{};
    std::uint64_t dp{};
    std::uint64_t backend_calls{};
    std::uint64_t estimated_cells{};
    std::uint64_t full_matrix_cells{};
    std::uint64_t band_calls{};
    std::uint64_t band_sum{};
    std::uint64_t band_min{};
    std::uint64_t band_max{};
    std::uint64_t block_calls{};
    std::uint64_t block_sum{};
    std::uint64_t block_min{};
    std::uint64_t block_max{};
    double backend_seconds{};
    bool processed{};
};

[[nodiscard]] SequenceLookup MakeSequenceLookup(
    const std::vector<SequenceRecord>& sequences) {
    SequenceLookup lookup;
    for (const SequenceRecord& sequence : sequences) {
        lookup.emplace(sequence.numeric_id, &sequence);
    }
    return lookup;
}

[[nodiscard]] AlignmentRecord BuildAlignment(
    const Chain& chain,
    const SequenceLookup& references,
    const SequenceLookup& queries,
    const std::unordered_map<SequenceId, std::string>& reverse_queries,
    const AlignmentOptions& options,
    GapCountAccumulator& statistics) {
    if (chain.empty()) {
        throw AlignmentError("internal error: cannot align an empty seed chain");
    }
    const auto reference_entry = references.find(chain.front().seed.reference_id);
    const auto query_entry = queries.find(chain.front().seed.query_id);
    if (reference_entry == references.end() || query_entry == queries.end()) {
        throw AlignmentError("seed chain refers to an unknown sequence id");
    }
    const SequenceRecord& reference = *reference_entry->second;
    const SequenceRecord& query = *query_entry->second;
    const std::string* oriented_query = &query.bases;
    if (chain.front().seed.strand == Strand::Reverse) {
        const auto reverse_entry = reverse_queries.find(query.numeric_id);
        if (reverse_entry == reverse_queries.end()) {
            throw AlignmentError("internal error: reverse query cache is incomplete");
        }
        oriented_query = &reverse_entry->second;
    }

    AlignmentRecord alignment;
    alignment.reference_id = reference.numeric_id;
    alignment.reference_begin = chain.front().seed.reference_begin;
    alignment.query_id = query.numeric_id;
    alignment.strand = chain.front().seed.strand;

    auto append_exact_seed = [&](const OrientedSeed& seed) {
        const Position reference_end =
            CheckedEnd(seed.seed.reference_begin, seed.seed.length, "reference seed");
        const Position query_end =
            CheckedEnd(seed.query_begin, seed.seed.length, "oriented query seed");
        if (reference_end > reference.size() || query_end > oriented_query->size()) {
            throw AlignmentError("seed chain lies outside its source sequence");
        }
        for (Length offset = 0; offset < seed.seed.length; ++offset) {
            const char reference_base =
                reference.bases[static_cast<std::size_t>(seed.seed.reference_begin + offset)];
            const char query_base =
                (*oriented_query)[static_cast<std::size_t>(seed.query_begin + offset)];
            if (!EqualCanonical(reference_base, query_base)) {
                throw AlignmentError("internal error: a seed is not an exact canonical match");
            }
        }
        AppendCigar(alignment.cigar, '=', seed.seed.length);
        alignment.score = CheckedScoreAdd(
            alignment.score,
            CheckedScoreProduct(seed.seed.length, options.match_score, "seed score"),
            "alignment score");
    };

    append_exact_seed(chain.front());
    for (std::size_t index = 1; index < chain.size(); ++index) {
        const OrientedSeed& previous = chain[index - 1];
        const OrientedSeed& current = chain[index];
        const Position previous_reference_end =
            CheckedEnd(previous.seed.reference_begin, previous.seed.length, "reference seed");
        const Position previous_query_end =
            CheckedEnd(previous.query_begin, previous.seed.length, "oriented query seed");
        if (previous_reference_end > current.seed.reference_begin ||
            previous_query_end > current.query_begin) {
            throw AlignmentError("internal error: chain contains overlapping seeds");
        }
        const Length reference_gap = current.seed.reference_begin - previous_reference_end;
        const Length query_gap = current.query_begin - previous_query_end;
        if (reference_gap != 0 || query_gap != 0) {
            const std::string_view reference_gap_bases{
                reference.bases.data() + static_cast<std::size_t>(previous_reference_end),
                static_cast<std::size_t>(reference_gap)};
            const std::string_view query_gap_bases{
                oriented_query->data() + static_cast<std::size_t>(previous_query_end),
                static_cast<std::size_t>(query_gap)};

            bool exact = reference_gap == query_gap;
            if (exact) {
                for (std::size_t offset = 0; offset < reference_gap_bases.size(); ++offset) {
                    if (!EqualCanonical(reference_gap_bases[offset], query_gap_bases[offset])) {
                        exact = false;
                        break;
                    }
                }
            }
            const bool realign_equal =
                !exact && ShouldRealignEqualGap(reference_gap_bases,
                                                query_gap_bases, options);
            if (exact) {
                ++statistics.exact;
            } else if (reference_gap == query_gap && !realign_equal) {
                ++statistics.ungapped;
            } else {
                ++statistics.dp;
            }

#if RAMAG_EXTENSION_METRICS || !defined(RAMAG_EXTENSION_SCALAR)
            ExtensionCallMetrics call_metrics;
#if RAMAG_EXTENSION_METRICS
            const Clock::time_point backend_begin = Clock::now();
#endif
#endif
            GapAlignment gap = AlignAffineGap(
                reference_gap_bases, query_gap_bases, options,
                realign_equal,
#if RAMAG_EXTENSION_METRICS || !defined(RAMAG_EXTENSION_SCALAR)
                &call_metrics);
#else
                nullptr);
#endif
#if RAMAG_EXTENSION_METRICS || !defined(RAMAG_EXTENSION_SCALAR)
            if (!reference_gap_bases.empty() && !query_gap_bases.empty() &&
                reference_gap_bases.size() != query_gap_bases.size()) {
                ++statistics.backend_calls;
                CheckedAggregate(statistics.estimated_cells,
                                 call_metrics.estimated_cells,
                                 "extension estimated cell count");
                CheckedAggregate(statistics.full_matrix_cells,
                                 call_metrics.full_matrix_cells,
                                 "extension full matrix cell count");
                if (call_metrics.effective_band_width != 0) {
                    ++statistics.band_calls;
                    CheckedAggregate(statistics.band_sum,
                                     call_metrics.effective_band_width,
                                     "extension band width sum");
                    statistics.band_min = statistics.band_min == 0
                                              ? call_metrics.effective_band_width
                                              : std::min(statistics.band_min,
                                                         call_metrics.effective_band_width);
                    statistics.band_max = std::max(
                        statistics.band_max, call_metrics.effective_band_width);
                }
                if (call_metrics.effective_block_size != 0) {
                    ++statistics.block_calls;
                    CheckedAggregate(statistics.block_sum,
                                     call_metrics.effective_block_size,
                                     "extension block size sum");
                    statistics.block_min = statistics.block_min == 0
                                               ? call_metrics.effective_block_size
                                               : std::min(statistics.block_min,
                                                          call_metrics.effective_block_size);
                    statistics.block_max = std::max(
                        statistics.block_max, call_metrics.effective_block_size);
                }
#if RAMAG_EXTENSION_METRICS
                statistics.backend_seconds +=
                    SecondsBetween(backend_begin, Clock::now());
#endif
            }
#endif
            for (const CigarOp& operation : gap.cigar) {
                AppendCigar(alignment.cigar, operation.operation, operation.length);
            }
            alignment.score =
                CheckedScoreAdd(alignment.score, gap.score, "alignment score");
        }
        append_exact_seed(current);
    }

    alignment.reference_end = CheckedEnd(
        chain.back().seed.reference_begin, chain.back().seed.length, "alignment reference");
    const Position oriented_query_begin = chain.front().query_begin;
    const Position oriented_query_end =
        CheckedEnd(chain.back().query_begin, chain.back().seed.length, "alignment query");
    if (alignment.strand == Strand::Forward) {
        alignment.query_begin = oriented_query_begin;
        alignment.query_end = oriented_query_end;
    } else {
        alignment.query_begin = query.size() - oriented_query_end;
        alignment.query_end = query.size() - oriented_query_begin;
    }
    alignment.edit_distance = CigarEditDistance(alignment.cigar);

    if (CigarReferenceLength(alignment.cigar) !=
            alignment.reference_end - alignment.reference_begin ||
        CigarQueryLength(alignment.cigar) != alignment.query_end - alignment.query_begin) {
        throw AlignmentError("internal error: CIGAR and alignment coordinates disagree");
    }
    return alignment;
}

struct AlignmentMetrics {
    Length query_length{};
    Length block_length{};
    Length matches{};
};

[[nodiscard]] AlignmentMetrics Metrics(const AlignmentRecord& alignment) {
    AlignmentMetrics metrics;
    for (const CigarOp& operation : alignment.cigar) {
        metrics.block_length =
            CheckedLengthAdd(metrics.block_length, operation.length, "alignment block length");
        if (operation.operation != 'D') {
            metrics.query_length = CheckedLengthAdd(
                metrics.query_length, operation.length, "aligned query length");
        }
        if (operation.operation == '=') {
            metrics.matches =
                CheckedLengthAdd(metrics.matches, operation.length, "alignment match length");
        }
    }
    return metrics;
}

enum class ConflictSide : std::uint8_t {
    Reference,
    Query,
};

[[nodiscard]] int CompareFractions(std::uint64_t left_numerator,
                                   std::uint64_t left_denominator,
                                   std::uint64_t right_numerator,
                                   std::uint64_t right_denominator) {
    if (left_denominator == 0 || right_denominator == 0) {
        throw AlignmentError("cannot compare an alignment identity with zero length");
    }
    bool inverted = false;
    while (true) {
        const std::uint64_t left_quotient = left_numerator / left_denominator;
        const std::uint64_t right_quotient = right_numerator / right_denominator;
        if (left_quotient != right_quotient) {
            const int ordinary = left_quotient < right_quotient ? -1 : 1;
            return inverted ? -ordinary : ordinary;
        }
        left_numerator %= left_denominator;
        right_numerator %= right_denominator;
        if (left_numerator == 0 || right_numerator == 0) {
            if (left_numerator == 0 && right_numerator == 0) {
                return 0;
            }
            const int ordinary = left_numerator == 0 ? -1 : 1;
            return inverted ? -ordinary : ordinary;
        }
        std::swap(left_numerator, left_denominator);
        std::swap(right_numerator, right_denominator);
        inverted = !inverted;
    }
}

[[nodiscard]] bool BetterPrimary(const AlignmentRecord& candidate,
                                 const AlignmentRecord& incumbent) {
    if (candidate.score != incumbent.score) {
        return candidate.score > incumbent.score;
    }
    const AlignmentMetrics candidate_metrics = Metrics(candidate);
    const AlignmentMetrics incumbent_metrics = Metrics(incumbent);
    if (candidate_metrics.query_length != incumbent_metrics.query_length) {
        return candidate_metrics.query_length > incumbent_metrics.query_length;
    }
    const int identity_comparison = CompareFractions(
        candidate_metrics.matches, candidate_metrics.block_length,
        incumbent_metrics.matches, incumbent_metrics.block_length);
    if (identity_comparison != 0) {
        return identity_comparison > 0;
    }
    return std::tie(candidate.reference_id, candidate.reference_begin,
                    candidate.query_begin, candidate.strand) <
           std::tie(incumbent.reference_id, incumbent.reference_begin,
                    incumbent.query_begin, incumbent.strand);
}

struct ConflictInterval {
    std::size_t alignment_index{};
    Position begin{};
    Position end{};
    Length matches{};
    Length block_length{};
};

[[nodiscard]] std::vector<bool> SelectBestCoverageAlignments(
    const std::vector<AlignmentRecord>& alignments,
    ConflictSide side) {
    std::unordered_map<SequenceId, std::vector<ConflictInterval>> groups;
    groups.reserve(alignments.size());
    for (std::size_t index = 0; index < alignments.size(); ++index) {
        const AlignmentRecord& alignment = alignments[index];
        const bool reference_side = side == ConflictSide::Reference;
        const SequenceId sequence_id = reference_side
                                           ? alignment.reference_id
                                           : alignment.query_id;
        const Position begin = reference_side
                                   ? alignment.reference_begin
                                   : alignment.query_begin;
        const Position end = reference_side
                                 ? alignment.reference_end
                                 : alignment.query_end;
        if (begin >= end) {
            throw AlignmentError(
                "cannot resolve an empty alignment interval");
        }
        const AlignmentMetrics metrics = Metrics(alignment);
        groups[sequence_id].push_back(ConflictInterval{
            index, begin, end, metrics.matches, metrics.block_length});
    }

    std::vector<bool> selected(alignments.size(), false);
    for (auto& [sequence_id, intervals] : groups) {
        static_cast<void>(sequence_id);
        std::sort(intervals.begin(), intervals.end(),
                  [](const ConflictInterval& left,
                     const ConflictInterval& right) {
                      return std::tie(left.begin, left.end,
                                      left.alignment_index) <
                             std::tie(right.begin, right.end,
                                      right.alignment_index);
                  });
        std::vector<Position> coordinates;
        coordinates.reserve(intervals.size() * 2U);
        for (const ConflictInterval& interval : intervals) {
            coordinates.push_back(interval.begin);
            coordinates.push_back(interval.end);
        }
        std::sort(coordinates.begin(), coordinates.end());
        coordinates.erase(std::unique(coordinates.begin(), coordinates.end()),
                          coordinates.end());

        const auto worse = [&](std::size_t left_index,
                               std::size_t right_index) {
            const ConflictInterval& left = intervals[left_index];
            const ConflictInterval& right = intervals[right_index];
            if (left.matches != right.matches) {
                return left.matches < right.matches;
            }
            const int identity = CompareFractions(
                left.matches, left.block_length,
                right.matches, right.block_length);
            if (identity != 0) {
                return identity < 0;
            }
            const Length left_span = left.end - left.begin;
            const Length right_span = right.end - right.begin;
            if (left_span != right_span) {
                return left_span < right_span;
            }
            return left.alignment_index > right.alignment_index;
        };
        std::priority_queue<std::size_t, std::vector<std::size_t>,
                            decltype(worse)>
            active(worse);
        std::size_t next_interval = 0;
        for (std::size_t coordinate_index = 0;
             coordinate_index + 1U < coordinates.size(); ++coordinate_index) {
            const Position coordinate = coordinates[coordinate_index];
            while (next_interval < intervals.size() &&
                   intervals[next_interval].begin == coordinate) {
                active.push(next_interval);
                ++next_interval;
            }
            while (!active.empty() &&
                   intervals[active.top()].end <= coordinate) {
                active.pop();
            }
            if (!active.empty() && coordinates[coordinate_index + 1U] > coordinate) {
                selected[intervals[active.top()].alignment_index] = true;
            }
        }
    }
    return selected;
}

void ResolveAlignmentRecordsImpl(std::vector<AlignmentRecord>& alignments,
                                 AlignmentSelection selection) {
    for (AlignmentRecord& alignment : alignments) {
        alignment.primary = false;
    }
    const auto content_less = [](const AlignmentRecord& left,
                                 const AlignmentRecord& right) {
        const auto left_key = std::tie(left.query_id,
                                       left.query_begin,
                                       left.query_end,
                                       left.reference_id,
                                       left.reference_begin,
                                       left.reference_end,
                                       left.strand,
                                       left.score,
                                       left.edit_distance);
        const auto right_key = std::tie(right.query_id,
                                        right.query_begin,
                                        right.query_end,
                                        right.reference_id,
                                        right.reference_begin,
                                        right.reference_end,
                                        right.strand,
                                        right.score,
                                        right.edit_distance);
        if (left_key != right_key) {
            return left_key < right_key;
        }
        return std::lexicographical_compare(
            left.cigar.begin(), left.cigar.end(), right.cigar.begin(), right.cigar.end(),
            [](const CigarOp& left_op, const CigarOp& right_op) {
                return std::tie(left_op.operation, left_op.length) <
                       std::tie(right_op.operation, right_op.length);
            });
    };
    std::sort(alignments.begin(), alignments.end(), content_less);
    alignments.erase(
        std::unique(alignments.begin(), alignments.end()), alignments.end());

    if (selection == AlignmentSelection::ReciprocalOneToOne) {
        const std::vector<bool> reference_selected =
            SelectBestCoverageAlignments(alignments, ConflictSide::Reference);
        const std::vector<bool> query_selected =
            SelectBestCoverageAlignments(alignments, ConflictSide::Query);
        std::size_t retained = 0;
        for (std::size_t index = 0; index < alignments.size(); ++index) {
            if (reference_selected[index] && query_selected[index]) {
                alignments[index].primary = true;
                ++retained;
            }
        }
        std::vector<AlignmentRecord> resolved;
        resolved.reserve(retained);
        for (AlignmentRecord& alignment : alignments) {
            if (alignment.primary) {
                resolved.push_back(std::move(alignment));
            }
        }
        alignments = std::move(resolved);
    } else {
        std::unordered_map<SequenceId, std::size_t> best_by_query;
        for (std::size_t index = 0; index < alignments.size(); ++index) {
            const auto existing = best_by_query.find(alignments[index].query_id);
            if (existing == best_by_query.end() ||
                BetterPrimary(alignments[index], alignments[existing->second])) {
                best_by_query[alignments[index].query_id] = index;
            }
        }
        for (const auto& [query_id, index] : best_by_query) {
            static_cast<void>(query_id);
            alignments[index].primary = true;
        }
    }

    std::sort(alignments.begin(), alignments.end(), [&](const AlignmentRecord& left,
                                                        const AlignmentRecord& right) {
        if (left.query_id != right.query_id) {
            return left.query_id < right.query_id;
        }
        if (left.primary != right.primary) {
            return left.primary > right.primary;
        }
        const auto left_coordinates =
            std::tie(left.query_begin, left.query_end, left.reference_id,
                     left.reference_begin, left.strand);
        const auto right_coordinates =
            std::tie(right.query_begin, right.query_end, right.reference_id,
                     right.reference_begin, right.strand);
        if (left_coordinates != right_coordinates) {
            return left_coordinates < right_coordinates;
        }
        if (left.score != right.score) {
            return left.score > right.score;
        }
        return content_less(left, right);
    });
}

void PopulateInputStatistics(RunStatistics& statistics,
                             const std::vector<SequenceRecord>& references,
                             const std::vector<SequenceRecord>& queries) {
    statistics.reference_contigs = references.size();
    statistics.query_contigs = queries.size();
    statistics.reference_bases = CountBases(references);
    statistics.query_bases = CountBases(queries);
    statistics.reference_ambiguous_bases = CountAmbiguous(references);
    statistics.query_ambiguous_bases = CountAmbiguous(queries);
}

void ValidateExternalSeeds(const std::vector<Seed>& seeds,
                           const SequenceLookup& references,
                           const SequenceLookup& queries,
                           const AlignmentOptions& options) {
    for (const Seed& seed : seeds) {
        if (seed.strand != Strand::Forward && seed.strand != Strand::Reverse) {
            throw AlignmentError("external seed has an invalid strand value");
        }
        if (seed.length < options.min_match) {
            throw AlignmentError("external seed is shorter than min_match");
        }
        const auto reference = references.find(seed.reference_id);
        if (reference == references.end()) {
            throw AlignmentError("external seed refers to unknown reference id " +
                                 std::to_string(seed.reference_id));
        }
        const auto query = queries.find(seed.query_id);
        if (query == queries.end()) {
            throw AlignmentError("external seed refers to unknown query id " +
                                 std::to_string(seed.query_id));
        }
        if (CheckedEnd(seed.reference_begin, seed.length, "external reference seed") >
            reference->second->size()) {
            throw AlignmentError("external seed lies outside reference id " +
                                 std::to_string(seed.reference_id));
        }
        if (CheckedEnd(seed.query_begin, seed.length, "external query seed") >
            query->second->size()) {
            throw AlignmentError("external seed lies outside query id " +
                                 std::to_string(seed.query_id));
        }
    }
}

[[nodiscard]] AlignmentResult FinishAlignmentFromSeeds(
    const std::vector<SequenceRecord>& references,
    const std::vector<SequenceRecord>& queries,
    const AlignmentOptions& options,
    std::vector<Seed> seeds,
    RunStatistics statistics,
    Clock::time_point total_begin,
    double external_seed_seconds) {
    AlignmentResult result;
    result.statistics = std::move(statistics);
    PopulateInputStatistics(result.statistics, references, queries);
    result.statistics.selected_seed_count = static_cast<std::uint64_t>(seeds.size());

    if (options.progress_callback) {
        options.progress_callback("seed-merge", 0,
                                  static_cast<std::uint64_t>(seeds.size()));
    }
    const Clock::time_point merge_begin = Clock::now();
    seeds = MergeDiagonalSeeds(std::move(seeds), queries);
    result.statistics.merged_seed_count = static_cast<std::uint64_t>(seeds.size());
    const Clock::time_point merge_end = Clock::now();
    result.statistics.merge_seconds = SecondsBetween(merge_begin, merge_end);

    if (options.progress_callback) {
        options.progress_callback("chaining", 0,
                                  result.statistics.merged_seed_count);
    }
    const Clock::time_point chain_begin = Clock::now();
    const QueryLengths query_lengths = MakeQueryLengths(queries);
    std::vector<Chain> chains = BuildChainsSparseExactEdgeComponents(
        seeds, query_lengths, options, result.statistics);
    result.statistics.chain_count = static_cast<std::uint64_t>(chains.size());

    if (options.progress_callback) {
        options.progress_callback("extension", 0,
                                  result.statistics.chain_count);
    }

    const SequenceLookup reference_lookup = MakeSequenceLookup(references);
    const SequenceLookup query_lookup = MakeSequenceLookup(queries);
    std::unordered_map<SequenceId, std::string> reverse_queries;
    std::vector<std::string> reverse_query_bases =
        ReverseComplements(queries, options.worker_threads);
    for (std::size_t query_index = 0; query_index < queries.size(); ++query_index) {
        reverse_queries.emplace(queries[query_index].numeric_id,
                                std::move(reverse_query_bases[query_index]));
    }
    if (chains.size() >
        static_cast<std::size_t>(std::numeric_limits<std::ptrdiff_t>::max())) {
        throw AlignmentError("chain count exceeds the OpenMP loop range");
    }
    result.alignments.resize(chains.size());
    std::vector<std::exception_ptr> extension_errors(chains.size());
    const std::uint32_t planned_extension_workers =
        PlannedWorkerCount(options.worker_threads, chains.size());
    std::vector<GapCountAccumulator> extension_statistics(
        static_cast<std::size_t>(planned_extension_workers));
    std::uint32_t actual_extension_workers = 1;
    const auto extend_chain = [&](std::size_t chain_index,
                                  std::size_t worker_index) {
        extension_statistics[worker_index].processed = true;
        try {
            if (options.interruption_callback) {
                options.interruption_callback("extension");
            }
            result.alignments[chain_index] = BuildAlignment(
                chains[chain_index], reference_lookup, query_lookup,
                reverse_queries, options, extension_statistics[worker_index]);
        } catch (...) {
            extension_errors[chain_index] = std::current_exception();
        }
    };
#if defined(_OPENMP)
    if (!chains.empty()) {
        const int omp_workers = static_cast<int>(planned_extension_workers);
#pragma omp parallel num_threads(omp_workers)
        {
#pragma omp for schedule(dynamic, 1)
            for (std::ptrdiff_t signed_index = 0;
                 signed_index < static_cast<std::ptrdiff_t>(chains.size());
                 ++signed_index) {
                extend_chain(static_cast<std::size_t>(signed_index),
                             static_cast<std::size_t>(omp_get_thread_num()));
            }
        }
    }
#else
    static_cast<void>(planned_extension_workers);
    for (std::size_t chain_index = 0; chain_index < chains.size(); ++chain_index) {
        extend_chain(chain_index, 0);
    }
#endif
    actual_extension_workers = 0;
    for (const GapCountAccumulator& local : extension_statistics) {
        if (local.processed) {
            ++actual_extension_workers;
        }
    }
    actual_extension_workers = std::max<std::uint32_t>(
        actual_extension_workers, 1);
    result.statistics.extension_worker_threads = actual_extension_workers;
    for (std::size_t chain_index = 0; chain_index < extension_errors.size();
         ++chain_index) {
        if (extension_errors[chain_index] != nullptr) {
            std::rethrow_exception(extension_errors[chain_index]);
        }
    }
    for (const GapCountAccumulator& local : extension_statistics) {
        CheckedAggregate(result.statistics.exact_gap_count,
                         local.exact, "exact gap count");
        CheckedAggregate(result.statistics.ungapped_gap_count,
                         local.ungapped, "ungapped gap count");
        CheckedAggregate(result.statistics.dp_gap_count,
                         local.dp, "DP gap count");
        CheckedAggregate(result.statistics.extension_backend_gap_count,
                         local.backend_calls, "extension backend gap count");
        CheckedAggregate(result.statistics.extension_estimated_cells,
                         local.estimated_cells, "extension estimated cells");
        CheckedAggregate(result.statistics.extension_full_matrix_cells,
                         local.full_matrix_cells,
                         "extension full matrix cells");
        CheckedAggregate(result.statistics.extension_band_call_count,
                         local.band_calls, "extension band call count");
        CheckedAggregate(result.statistics.extension_band_width_sum,
                         local.band_sum, "extension band width sum");
        CheckedAggregate(result.statistics.extension_block_call_count,
                         local.block_calls, "extension block call count");
        CheckedAggregate(result.statistics.extension_block_size_sum,
                         local.block_sum, "extension block size sum");
        if (local.band_min != 0) {
            result.statistics.extension_band_width_min =
                result.statistics.extension_band_width_min == 0
                    ? local.band_min
                    : std::min(result.statistics.extension_band_width_min,
                               local.band_min);
            result.statistics.extension_band_width_max = std::max(
                result.statistics.extension_band_width_max, local.band_max);
        }
        if (local.block_min != 0) {
            result.statistics.extension_block_size_min =
                result.statistics.extension_block_size_min == 0
                    ? local.block_min
                    : std::min(result.statistics.extension_block_size_min,
                               local.block_min);
            result.statistics.extension_block_size_max = std::max(
                result.statistics.extension_block_size_max, local.block_max);
        }
        result.statistics.extension_backend_call_seconds +=
            local.backend_seconds;
    }
    const Clock::time_point chain_end = Clock::now();
    result.statistics.chain_and_extension_seconds =
        SecondsBetween(chain_begin, chain_end);
    result.statistics.candidate_alignment_count =
        static_cast<std::uint64_t>(result.alignments.size());
    if (options.progress_callback) {
        options.progress_callback("conflict-resolution", 0,
                                  result.statistics.candidate_alignment_count);
    }
    const Clock::time_point conflict_begin = Clock::now();
    ResolveAlignmentRecords(result.alignments, options.selection);
    const Clock::time_point conflict_end = Clock::now();
    result.statistics.alignment_count =
        static_cast<std::uint64_t>(result.alignments.size());
    result.statistics.conflict_rejected_alignment_count =
        result.statistics.candidate_alignment_count -
        result.statistics.alignment_count;
    result.statistics.conflict_resolution_seconds =
        SecondsBetween(conflict_begin, conflict_end);
    result.statistics.total_seconds =
        external_seed_seconds + SecondsBetween(total_begin, Clock::now());
    return result;
}

[[nodiscard]] Length CigarLengthFor(std::span<const CigarOp> cigar,
                                    std::string_view consuming_operations,
                                    std::string_view what) {
    Length result = 0;
    for (const CigarOp& operation : cigar) {
        if (operation.length == 0) {
            throw AlignmentError("zero-length CIGAR operation");
        }
        if (operation.operation != 'M' && operation.operation != '=' &&
            operation.operation != 'X' && operation.operation != 'I' &&
            operation.operation != 'D') {
            throw AlignmentError("unsupported CIGAR operation '" +
                                 std::string(1, operation.operation) + "'");
        }
        if (consuming_operations.find(operation.operation) != std::string_view::npos) {
            result = CheckedLengthAdd(result, operation.length, what);
        }
    }
    return result;
}

}  // namespace

void ResolveAlignmentRecords(std::vector<AlignmentRecord>& alignments,
                             AlignmentSelection selection) {
    ResolveAlignmentRecordsImpl(alignments, selection);
}

std::vector<Seed> EnumerateSeeds(const std::vector<SequenceRecord>& references,
                                 const std::vector<SequenceRecord>& queries,
                                 const AlignmentOptions& options,
                                 RunStatistics* statistics) {
    ValidateOptions(options);
    ValidateSequenceCollection(references, "reference");
    ValidateSequenceCollection(queries, "query");
    if (options.min_match > std::numeric_limits<std::size_t>::max()) {
        if (statistics != nullptr) {
            statistics->mem_seed_count = 0;
            statistics->mam_seed_count = 0;
            statistics->selected_seed_count = 0;
            statistics->actual_seed_route = SeedRouteName(options.seed_mode);
        }
        return {};
    }
    const std::size_t kmer_length = static_cast<std::size_t>(options.min_match);
    if (options.seed_mode == SeedMode::Smem) {
        SmemOracleResult smems = EnumerateSmemOracle(
            references, queries, kmer_length,
            options.smem_min_occurrences);
        if (statistics != nullptr) {
            statistics->smem_interval_count = smems.interval_count;
            statistics->smem_coordinate_seed_count =
                static_cast<std::uint64_t>(smems.seeds.size());
            statistics->selected_seed_count =
                static_cast<std::uint64_t>(smems.seeds.size());
            statistics->actual_seed_route = SeedRouteName(options.seed_mode);
        }
        return smems.seeds;
    }
    const KmerBuckets reference_kmers = BuildReferenceKmers(references, kmer_length);

    std::vector<RawSeed> mems;
    for (std::size_t query_index = 0; query_index < queries.size(); ++query_index) {
        const SequenceRecord& query = queries[query_index];
        EnumerateOrientedMems(references, query, query_index, query.bases,
                              Strand::Forward, kmer_length, reference_kmers, mems);
        const std::string reverse = ReverseComplement(query.bases);
        EnumerateOrientedMems(references, query, query_index, reverse,
                              Strand::Reverse, kmer_length, reference_kmers, mems);
    }

    std::vector<Seed> mams;
    std::vector<RawSeed> raw_mams;
    mams.reserve(mems.size());
    raw_mams.reserve(mems.size());
    for (const RawSeed& mem : mems) {
        if (IsReferenceUnique(mem, references, kmer_length, reference_kmers)) {
            mams.push_back(mem.seed);
            raw_mams.push_back(mem);
        }
    }

    std::vector<Seed> mums;
    if (options.seed_mode == SeedMode::Mum) {
        mums.reserve(raw_mams.size());
        for (const RawSeed& mam : raw_mams) {
            if (IsQueryUnique(mam, references, queries)) {
                mums.push_back(mam.seed);
            }
        }
    }

    std::vector<Seed> selected;
    if (options.seed_mode == SeedMode::MaxMatch) {
        selected.reserve(mems.size());
        for (const RawSeed& mem : mems) {
            selected.push_back(mem.seed);
        }
    } else if (options.seed_mode == SeedMode::MumReference) {
        selected = mams;
    } else if (options.seed_mode == SeedMode::Mum) {
        selected = mums;
    } else {
        selected = mams;
        std::unordered_map<std::uint64_t, std::vector<std::pair<Position, Position>>> coverage;
        for (const Seed& mam : mams) {
            coverage[QueryStrandKey(mam)].push_back(
                {mam.query_begin, CheckedEnd(mam.query_begin, mam.length, "MAM query interval")});
        }
        for (auto& [key, intervals] : coverage) {
            static_cast<void>(key);
            std::sort(intervals.begin(), intervals.end());
            std::vector<std::pair<Position, Position>> merged;
            for (const auto& interval : intervals) {
                if (merged.empty() || interval.first > merged.back().second) {
                    merged.push_back(interval);
                } else {
                    merged.back().second = std::max(merged.back().second, interval.second);
                }
            }
            intervals = std::move(merged);
        }

        for (const RawSeed& raw_mem : mems) {
            const Seed& mem = raw_mem.seed;
            const Position mem_end = CheckedEnd(mem.query_begin, mem.length, "MEM query interval");
            bool covered = false;
            const auto intervals = coverage.find(QueryStrandKey(mem));
            if (intervals != coverage.end()) {
                for (const auto& interval : intervals->second) {
                    if (interval.first > mem.query_begin) {
                        break;
                    }
                    if (interval.first <= mem.query_begin && interval.second >= mem_end) {
                        covered = true;
                        break;
                    }
                }
            }
            if (!covered) {
                selected.push_back(mem);
            }
        }
    }

    std::sort(selected.begin(), selected.end(), PublicSeedLess);
    selected.erase(std::unique(selected.begin(), selected.end()), selected.end());
    if (statistics != nullptr) {
        statistics->mem_seed_count = static_cast<std::uint64_t>(mems.size());
        statistics->mam_seed_count = static_cast<std::uint64_t>(mams.size());
        statistics->mum_seed_count = static_cast<std::uint64_t>(mums.size());
        statistics->selected_seed_count = static_cast<std::uint64_t>(selected.size());
        statistics->actual_seed_route = SeedRouteName(options.seed_mode);
    }
    return selected;
}

std::vector<Seed> MergeDiagonalSeeds(std::vector<Seed> seeds,
                                     const std::vector<SequenceRecord>& queries) {
    if (seeds.empty()) {
        return {};
    }
    const QueryLengths query_lengths = MakeQueryLengths(queries);
    std::vector<OrientedSeed> oriented;
    oriented.reserve(seeds.size());
    for (Seed& seed : seeds) {
        if (seed.length == 0) {
            throw AlignmentError("cannot merge a zero-length seed");
        }
        const auto query_length = query_lengths.find(seed.query_id);
        if (query_length == query_lengths.end()) {
            throw AlignmentError("seed refers to unknown query id " +
                                 std::to_string(seed.query_id));
        }
        const Position oriented_begin =
            OrientedQueryBegin(seed, query_length->second);
        oriented.push_back(OrientedSeed{std::move(seed), oriented_begin});
    }

    std::sort(oriented.begin(), oriented.end(), [](const OrientedSeed& left,
                                                   const OrientedSeed& right) {
        const auto left_group =
            std::tie(left.seed.reference_id, left.seed.query_id, left.seed.strand);
        const auto right_group =
            std::tie(right.seed.reference_id, right.seed.query_id, right.seed.strand);
        if (left_group != right_group) {
            return left_group < right_group;
        }
        if (!SameDiagonal(left, right)) {
            return DiagonalLess(left, right);
        }
        return std::tie(left.query_begin, left.seed.reference_begin, left.seed.length) <
               std::tie(right.query_begin, right.seed.reference_begin, right.seed.length);
    });

    std::vector<OrientedSeed> merged;
    for (OrientedSeed& seed : oriented) {
        if (!merged.empty() && SameSeedGroup(merged.back(), seed) &&
            SameDiagonal(merged.back(), seed)) {
            const Position merged_end =
                CheckedEnd(merged.back().query_begin, merged.back().seed.length,
                           "merged query seed");
            const Position seed_end =
                CheckedEnd(seed.query_begin, seed.seed.length, "query seed");
            if (seed.query_begin <= merged_end) {
                const Position new_end = std::max(merged_end, seed_end);
                merged.back().seed.length = new_end - merged.back().query_begin;
                continue;
            }
        }
        merged.push_back(std::move(seed));
    }

    std::vector<Seed> result;
    result.reserve(merged.size());
    for (const OrientedSeed& seed : merged) {
        const auto query_length = query_lengths.find(seed.seed.query_id);
        result.push_back(ToPublicSeed(seed, query_length->second));
    }
    return result;
}

AlignmentResult AlignGenomes(const std::vector<SequenceRecord>& references,
                             const std::vector<SequenceRecord>& queries,
                             const AlignmentOptions& options) {
    const Clock::time_point total_begin = Clock::now();
    ValidateOptions(options);
    ValidateSequenceCollection(references, "reference");
    ValidateSequenceCollection(queries, "query");

    RunStatistics statistics;
    PopulateInputStatistics(statistics, references, queries);
    const Clock::time_point seed_begin = Clock::now();
    std::vector<Seed> seeds =
        EnumerateSeeds(references, queries, options, &statistics);
    const Clock::time_point seed_end = Clock::now();
    statistics.seed_seconds = SecondsBetween(seed_begin, seed_end);
    return FinishAlignmentFromSeeds(references, queries, options, std::move(seeds),
                                    std::move(statistics), total_begin, 0.0);
}

AlignmentResult AlignGenomesFromSeeds(
    const std::vector<SequenceRecord>& references,
    const std::vector<SequenceRecord>& queries,
    const AlignmentOptions& options,
    std::vector<Seed> seeds,
    RunStatistics seed_statistics) {
    const Clock::time_point total_begin = Clock::now();
    ValidateOptions(options);
    ValidateSequenceCollection(references, "reference");
    ValidateSequenceCollection(queries, "query");
    const SequenceLookup reference_lookup = MakeSequenceLookup(references);
    const SequenceLookup query_lookup = MakeSequenceLookup(queries);
    ValidateExternalSeeds(seeds, reference_lookup, query_lookup, options);
    seed_statistics.selected_seed_count = static_cast<std::uint64_t>(seeds.size());
    if (seed_statistics.actual_seed_route.empty()) {
        seed_statistics.actual_seed_route = "external-exact-seeds";
    }
    const double external_seed_seconds = seed_statistics.seed_seconds;
    return FinishAlignmentFromSeeds(references, queries, options, std::move(seeds),
                                    std::move(seed_statistics), total_begin,
                                    external_seed_seconds);
}

std::string CigarToString(std::span<const CigarOp> cigar) {
    if (cigar.empty()) {
        return "*";
    }
    std::string result;
    for (const CigarOp& operation : cigar) {
        if (operation.length == 0) {
            throw AlignmentError("zero-length CIGAR operation");
        }
        if (operation.operation != 'M' && operation.operation != '=' &&
            operation.operation != 'X' && operation.operation != 'I' &&
            operation.operation != 'D') {
            throw AlignmentError("unsupported CIGAR operation '" +
                                 std::string(1, operation.operation) + "'");
        }
        result += std::to_string(operation.length);
        result.push_back(operation.operation);
    }
    return result;
}

Length CigarReferenceLength(std::span<const CigarOp> cigar) {
    return CigarLengthFor(cigar, "M=XD", "CIGAR reference length");
}

Length CigarQueryLength(std::span<const CigarOp> cigar) {
    return CigarLengthFor(cigar, "M=XI", "CIGAR query length");
}

Length CigarMatchLength(std::span<const CigarOp> cigar) {
    return CigarLengthFor(cigar, "=", "CIGAR exact-match length");
}

Length CigarEditDistance(std::span<const CigarOp> cigar) {
    return CigarLengthFor(cigar, "XID", "CIGAR edit distance");
}

namespace testing {
namespace {

[[nodiscard]] QueryLengths MakeTestingQueryLengths(
    std::span<const std::pair<SequenceId, Length>> query_lengths) {
    QueryLengths result;
    for (const auto& [id, length] : query_lengths) {
        if (!result.emplace(id, length).second) {
            throw AlignmentError("duplicate testing query id " + std::to_string(id));
        }
    }
    return result;
}

[[nodiscard]] std::vector<std::vector<Seed>> PublicChains(
    const std::vector<Chain>& chains) {
    std::vector<std::vector<Seed>> result;
    result.reserve(chains.size());
    for (const Chain& chain : chains) {
        std::vector<Seed> public_chain;
        public_chain.reserve(chain.size());
        for (const OrientedSeed& seed : chain) {
            public_chain.push_back(seed.seed);
        }
        result.push_back(std::move(public_chain));
    }
    return result;
}

}  // namespace

ExtensionGapTestResult AlignExtensionGapForTesting(
    std::string_view reference,
    std::string_view query,
    const AlignmentOptions& options) {
    ValidateOptions(options);
    if (reference.empty() || query.empty() || reference.size() == query.size()) {
        throw AlignmentError(
            "extension differential testing requires non-empty unequal gaps");
    }
    ExtensionCallMetrics metrics;
    const GapAlignment configured = AlignConfiguredUnequalGap(
        reference, query, options, &metrics);
    const GapAlignment scalar = AlignScalarAffineGap(
        reference, query, options, true);
    static_cast<void>(ValidateAndScoreGapAlignment(
        reference, query, options, scalar));
    return ExtensionGapTestResult{
        RAMAG_EXTENSION_BACKEND,
        configured.cigar,
        scalar.cigar,
        configured.score,
        scalar.score,
        metrics.estimated_cells,
        metrics.full_matrix_cells,
        metrics.effective_band_width,
        metrics.effective_block_size};
}

ConfiguredExtensionGapResult AlignConfiguredExtensionGapForTesting(
    std::string_view reference,
    std::string_view query,
    const AlignmentOptions& options) {
    ValidateOptions(options);
    if (reference.empty() || query.empty() || reference.size() == query.size()) {
        throw AlignmentError(
            "extension benchmarking requires non-empty unequal gaps");
    }
    ExtensionCallMetrics metrics;
    GapAlignment configured = AlignConfiguredUnequalGap(
        reference, query, options, &metrics);
    return ConfiguredExtensionGapResult{
        std::move(configured.cigar), configured.score,
        metrics.estimated_cells, metrics.full_matrix_cells,
        metrics.effective_band_width,
        metrics.effective_block_size};
}

Length Ksw2AutomaticBandWidthForTesting(
    Length reference_length,
    Length query_length,
    const AlignmentOptions& options) {
    ValidateOptions(options);
    return internal::Ksw2AutomaticBandWidth(
        reference_length, query_length, options);
}

std::string ConfiguredExtensionBackendName() {
    return std::string{internal::ConfiguredExtensionBackend()};
}

ChainingTestResult BuildSparseChainsForTesting(
    std::span<const Seed> seeds,
    std::span<const std::pair<SequenceId, Length>> query_lengths,
    const AlignmentOptions& options) {
    ValidateOptions(options);
    const QueryLengths lengths = MakeTestingQueryLengths(query_lengths);
    std::vector<Seed> owned_seeds(seeds.begin(), seeds.end());
    ChainingTestResult result;
    result.chains = PublicChains(BuildChainsSparseExactEdgeComponents(
        owned_seeds, lengths, options, result.statistics));
    return result;
}

void VerifySparseChainingAgainstQuadraticOracle(
    std::span<const Seed> seeds,
    std::span<const std::pair<SequenceId, Length>> query_lengths,
    const AlignmentOptions& options) {
    ValidateOptions(options);
    const QueryLengths lengths = MakeTestingQueryLengths(query_lengths);
    std::vector<Seed> owned_seeds(seeds.begin(), seeds.end());
    RunStatistics sparse_statistics;
    const std::vector<Chain> sparse = BuildChainsSparseExactEdgeComponents(
        owned_seeds, lengths, options, sparse_statistics);
    const std::vector<Chain> oracle =
        BuildChainsQuadraticOracle(owned_seeds, lengths, options);
    if (PublicChains(sparse) != PublicChains(oracle)) {
        throw AlignmentError(
            "sparse-exact-edge-components-v1 differs from the quadratic chaining oracle");
    }
}

}  // namespace testing

}  // namespace ramag

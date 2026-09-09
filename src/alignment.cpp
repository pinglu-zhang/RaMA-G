#include "ramag/alignment.hpp"

#include "ramag/fasta.hpp"
#include "ramag/pairwise_core.hpp"

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
using SequenceLookup = std::unordered_map<SequenceId, const SequenceRecord*>;

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

[[nodiscard]] SequenceLookup MakeSequenceLookup(
    const std::vector<SequenceRecord>& sequences) {
    SequenceLookup lookup;
    for (const SequenceRecord& sequence : sequences) {
        lookup.emplace(sequence.numeric_id, &sequence);
    }
    return lookup;
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
    if (options.break_length != 200 || options.max_dp_cells != 4000000 ||
        options.match_score != 2 || options.mismatch_penalty != 4 ||
        options.gap_open_penalty != 4 || options.gap_extend_penalty != 2) {
        throw AlignmentError("legacy extension/scoring overrides are not applicable to the pairwise core");
    }
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
    return AlignPairwiseFromSeeds(references, queries, options, std::move(seeds), std::move(statistics));
}

AlignmentResult AlignGenomesFromSeeds(
    const std::vector<SequenceRecord>& references,
    const std::vector<SequenceRecord>& queries,
    const AlignmentOptions& options,
    std::vector<Seed> seeds,
    RunStatistics seed_statistics) {
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
    return AlignPairwiseFromSeeds(references, queries, options, std::move(seeds), std::move(seed_statistics));
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

} // namespace ramag

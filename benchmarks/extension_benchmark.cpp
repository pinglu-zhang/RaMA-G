#include "ramag/alignment.hpp"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;

struct Case {
    std::string reference;
    std::string query;
};

std::uint64_t ParsePositive(std::string_view text, std::string_view option) {
    std::size_t consumed = 0;
    const unsigned long long value = std::stoull(std::string{text}, &consumed);
    if (consumed != text.size() || value == 0) {
        throw std::runtime_error(std::string{option} + " requires a positive integer");
    }
    return static_cast<std::uint64_t>(value);
}

void GenerateWords(std::size_t length,
                   std::string& current,
                   std::vector<std::string>& words) {
    constexpr std::string_view alphabet = "ACGT";
    if (current.size() == length) {
        words.push_back(current);
        return;
    }
    for (const char base : alphabet) {
        current.push_back(base);
        GenerateWords(length, current, words);
        current.pop_back();
    }
}

std::vector<Case> MakeCases(std::uint64_t random_count, bool exhaustive) {
    std::vector<Case> cases{
        {"A", "AA"},
        {"AA", "A"},
        {"ACGT", "ACGTT"},
        {"ACGTT", "ACGT"},
        {"AAAAAAA", "AAAA"},
        {"ACGTACGT", "ACGACGT"},
        {"ACNGT", "ACNNGT"},
        {"NNNN", "NNNNN"},
        {"ACACACACAC", "ACACACAC"},
        {"GATTACA", "GATTTACA"}};
    cases.push_back(Case{
        std::string(30, 'A') + std::string(30, 'C') + std::string(30, 'G'),
        std::string(30, 'C') + std::string(30, 'G') + std::string(29, 'T')});
    // Near-pure long deletion/insertion cases keep both inputs non-empty, as
    // required by the experimental backend contract, while exercising the
    // largest default gap dimension and a band endpoint set by |r-q|.
    cases.push_back(Case{std::string(90, 'A'), "A"});
    cases.push_back(Case{"A", std::string(90, 'A')});

    if (exhaustive) {
        std::vector<std::vector<std::string>> by_length(5);
        for (std::size_t length = 1; length <= 4; ++length) {
            std::string current;
            GenerateWords(length, current, by_length[length]);
        }
        for (std::size_t reference_length = 1; reference_length <= 4;
             ++reference_length) {
            for (std::size_t query_length = 1; query_length <= 4; ++query_length) {
                if (reference_length == query_length) continue;
                for (const std::string& reference : by_length[reference_length]) {
                    for (const std::string& query : by_length[query_length]) {
                        cases.push_back(Case{reference, query});
                    }
                }
            }
        }
    }

    std::mt19937_64 random(0x52414d41475f4450ULL);
    constexpr std::string_view bases = "ACGT";
    for (std::uint64_t index = 0; index < random_count; ++index) {
        const std::size_t reference_length = 1 + static_cast<std::size_t>(random() % 90);
        const int delta = static_cast<int>(random() % 33) - 16;
        const std::size_t query_length = static_cast<std::size_t>(std::clamp<long long>(
            static_cast<long long>(reference_length) + delta, 1, 90));
        if (query_length == reference_length) {
            --index;
            continue;
        }
        Case item;
        item.reference.resize(reference_length);
        item.query.resize(query_length);
        for (std::size_t position = 0; position < reference_length; ++position) {
            item.reference[position] = bases[static_cast<std::size_t>(random() % 4)];
        }
        const std::size_t shared = std::min(reference_length, query_length);
        for (std::size_t position = 0; position < shared; ++position) {
            const bool mutate = random() % 10 == 0;
            item.query[position] = mutate
                                       ? bases[static_cast<std::size_t>(random() % 4)]
                                       : item.reference[position];
        }
        for (std::size_t position = shared; position < query_length; ++position) {
            item.query[position] = bases[static_cast<std::size_t>(random() % 4)];
        }
        if (index % 17 == 0) {
            std::fill(item.reference.begin(), item.reference.end(), 'A');
            std::fill(item.query.begin(), item.query.end(), 'A');
        }
        if (index % 29 == 0) {
            item.reference[reference_length / 2] = 'N';
            item.query[query_length / 2] = 'N';
        }
        cases.push_back(std::move(item));
    }
    return cases;
}

template <typename Value>
Value Percentile(std::vector<Value> values, double fraction) {
    if (values.empty()) return {};
    std::sort(values.begin(), values.end());
    const std::size_t index = static_cast<std::size_t>(
        fraction * static_cast<double>(values.size() - 1));
    return values[index];
}

template <typename Value>
double Mean(const std::vector<Value>& values) {
    if (values.empty()) return 0.0;
    long double total = 0.0L;
    for (const Value value : values) total += static_cast<long double>(value);
    return static_cast<double>(total / static_cast<long double>(values.size()));
}

void CheckedAdd(std::uint64_t& total,
                std::uint64_t value,
                std::string_view label) {
    if (value > std::numeric_limits<std::uint64_t>::max() - total) {
        throw std::runtime_error(std::string{label} + " exceeds uint64_t");
    }
    total += value;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        std::uint64_t random_cases = 5'000;
        std::uint64_t repetitions = 5;
        bool exhaustive = false;
        for (int index = 1; index < argc; ++index) {
            const std::string_view option = argv[index];
            if (option == "--cases" && index + 1 < argc) {
                random_cases = ParsePositive(argv[++index], option);
            } else if (option == "--repetitions" && index + 1 < argc) {
                repetitions = ParsePositive(argv[++index], option);
            } else if (option == "--exhaustive") {
                exhaustive = true;
            } else {
                throw std::runtime_error("unknown or incomplete option: " +
                                         std::string{option});
            }
        }

        const ramag::AlignmentOptions default_options;
        const auto require_band = [&](ramag::Length reference_length,
                                      ramag::Length query_length,
                                      ramag::Length expected) {
            const ramag::Length actual =
                ramag::testing::Ksw2AutomaticBandWidthForTesting(
                    reference_length, query_length, default_options);
            if (actual != expected) {
                throw std::runtime_error(
                    "ksw2-auto-band-v1 self-check failed for " +
                    std::to_string(reference_length) + "x" +
                    std::to_string(query_length) + ": expected " +
                    std::to_string(expected) + ", got " +
                    std::to_string(actual));
            }
        };
        require_band(0, 0, 0);
        require_band(1, 2, 2);
        require_band(10, 0, 10);
        require_band(90, 90, 16);
        require_band(90, 80, 16);

        const std::vector<Case> cases = MakeCases(random_cases, exhaustive);
        std::uint64_t score_matches = 0;
        std::uint64_t cigar_matches = 0;
        std::vector<std::uint64_t> deficits;
        std::vector<std::uint64_t> bands;
        std::vector<std::uint64_t> blocks;
        std::uint64_t estimated_cells = 0;
        std::uint64_t full_matrix_cells = 0;
        for (const Case& item : cases) {
            ramag::testing::ExtensionGapTestResult result;
            try {
                result = ramag::testing::AlignExtensionGapForTesting(
                    item.reference, item.query);
            } catch (const std::exception& error) {
                throw std::runtime_error(
                    "case reference='" + item.reference + "' query='" +
                    item.query + "': " + error.what());
            }
            if (result.configured_score == result.scalar_score) ++score_matches;
            if (result.configured_cigar == result.scalar_cigar) ++cigar_matches;
            deficits.push_back(result.configured_score < result.scalar_score
                                   ? static_cast<std::uint64_t>(
                                         result.scalar_score - result.configured_score)
                                   : 0);
            if (result.effective_band_width != 0) {
                bands.push_back(result.effective_band_width);
            }
            if (result.effective_block_size != 0) {
                blocks.push_back(result.effective_block_size);
            }
            CheckedAdd(estimated_cells, result.estimated_cells,
                       "estimated cell total");
            CheckedAdd(full_matrix_cells, result.full_matrix_cells,
                       "full matrix cell total");
        }

        std::uint64_t checksum = 0;
        for (const Case& item : cases) {
            const auto result =
                ramag::testing::AlignConfiguredExtensionGapForTesting(
                    item.reference, item.query);
            checksum ^= static_cast<std::uint64_t>(result.score) +
                        result.cigar.size();
        }
        std::vector<std::uint64_t> latency_nanoseconds;
        if (repetitions > std::numeric_limits<std::size_t>::max() /
                              cases.size()) {
            throw std::runtime_error("timed call count exceeds size_t");
        }
        latency_nanoseconds.reserve(
            static_cast<std::size_t>(repetitions) * cases.size());
        const Clock::time_point begin = Clock::now();
        for (std::uint64_t repetition = 0; repetition < repetitions; ++repetition) {
            for (const Case& item : cases) {
                const Clock::time_point call_begin = Clock::now();
                const auto result =
                    ramag::testing::AlignConfiguredExtensionGapForTesting(
                        item.reference, item.query);
                const auto call_nanoseconds =
                    std::chrono::duration_cast<std::chrono::nanoseconds>(
                        Clock::now() - call_begin).count();
                latency_nanoseconds.push_back(static_cast<std::uint64_t>(
                    std::max<std::int64_t>(0, call_nanoseconds)));
                checksum ^= static_cast<std::uint64_t>(result.score) +
                            result.cigar.size() + repetition;
            }
        }
        const double seconds =
            std::chrono::duration<double>(Clock::now() - begin).count();
        const std::uint64_t case_count =
            static_cast<std::uint64_t>(cases.size());
        const std::uint64_t calls = repetitions * case_count;
        const auto maximum_deficit =
            *std::max_element(deficits.begin(), deficits.end());

        std::cout << std::fixed << std::setprecision(9)
                  << "{\n"
                  << "  \"backend\": \""
                  << ramag::testing::ConfiguredExtensionBackendName() << "\",\n"
                  << "  \"cases\": " << cases.size() << ",\n"
                  << "  \"repetitions\": " << repetitions << ",\n"
                  << "  \"timed_calls\": " << calls << ",\n"
                  << "  \"elapsed_seconds\": " << seconds << ",\n"
                  << "  \"calls_per_second\": "
                  << static_cast<double>(calls) / seconds << ",\n"
                  << "  \"latency_ns_p50\": "
                  << Percentile(latency_nanoseconds, 0.50) << ",\n"
                  << "  \"latency_ns_p95\": "
                  << Percentile(latency_nanoseconds, 0.95) << ",\n"
                  << "  \"latency_ns_p99\": "
                  << Percentile(latency_nanoseconds, 0.99) << ",\n"
                  << "  \"exact_score_matches\": " << score_matches << ",\n"
                  << "  \"exact_score_rate\": "
                  << static_cast<double>(score_matches) /
                         static_cast<double>(case_count) << ",\n"
                  << "  \"canonical_cigar_matches\": " << cigar_matches << ",\n"
                  << "  \"canonical_cigar_rate\": "
                  << static_cast<double>(cigar_matches) /
                         static_cast<double>(case_count) << ",\n"
                  << "  \"score_deficit_mean\": " << Mean(deficits) << ",\n"
                  << "  \"score_deficit_p50\": " << Percentile(deficits, 0.50) << ",\n"
                  << "  \"score_deficit_p95\": " << Percentile(deficits, 0.95) << ",\n"
                  << "  \"score_deficit_max\": " << maximum_deficit << ",\n"
                  << "  \"band_width_min\": "
                  << (bands.empty() ? 0 : *std::min_element(bands.begin(), bands.end())) << ",\n"
                  << "  \"band_width_p50\": " << Percentile(bands, 0.50) << ",\n"
                  << "  \"band_width_mean\": " << Mean(bands) << ",\n"
                  << "  \"band_width_p95\": " << Percentile(bands, 0.95) << ",\n"
                  << "  \"band_width_max\": "
                  << (bands.empty() ? 0 : *std::max_element(bands.begin(), bands.end())) << ",\n"
                  << "  \"block_size_min\": "
                  << (blocks.empty() ? 0 : *std::min_element(blocks.begin(), blocks.end())) << ",\n"
                  << "  \"block_size_p50\": " << Percentile(blocks, 0.50) << ",\n"
                  << "  \"block_size_mean\": " << Mean(blocks) << ",\n"
                  << "  \"block_size_p95\": " << Percentile(blocks, 0.95) << ",\n"
                  << "  \"block_size_max\": "
                  << (blocks.empty() ? 0 : *std::max_element(blocks.begin(), blocks.end())) << ",\n"
                  << "  \"estimated_cells\": " << estimated_cells << ",\n"
                  << "  \"full_matrix_cells\": " << full_matrix_cells << ",\n"
                  << "  \"estimated_to_full_cell_ratio\": ";
        if (estimated_cells == 0 || full_matrix_cells == 0) {
            std::cout << "null";
        } else {
            std::cout << static_cast<double>(estimated_cells) /
                             static_cast<double>(full_matrix_cells);
        }
        std::cout << ",\n"
                  << "  \"allocation_events\": null,\n"
                  << "  \"checksum\": " << checksum << "\n"
                  << "}\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "extension benchmark failed: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}

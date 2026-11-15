// Pairwise alignment kernel. Source attribution and MIT notice are retained
// in THIRD_PARTY_NOTICES.md. Do not substitute unrelated gap backends.
#include "ramag/pairwise_core.hpp"
#include "ksw2.h"
#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <chrono>
#include <cstdlib>
#include <exception>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <optional>
#include <queue>
#include <set>
#include <type_traits>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <stdexcept>
#include <tuple>
#include <unordered_map>
#include <variant>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "ramag/internal/genome.hpp"
namespace ramag::pairwise_detail {
uint64_t cigarToInt(char operation, uint64_t len) {
	uint64_t opCode;

	// 将字符型操作转换为对应的数字编码（4 位）
	switch (operation) {
	case 'M': opCode = 0x0; break; // 对齐匹配（可能包含错配）
	case 'I': opCode = 0x1; break; // 插入
	case 'D': opCode = 0x2; break; // 删除
	case '=': opCode = 0x7; break; // 精确匹配（无错配）
	case 'X': opCode = 0x8; break; // 错配
		// TODO: 可以根据需要添加更多操作类型（如 soft clip, hard clip 等）
	default: opCode = 0xF; break; // 未知操作，使用保留码 0xF
	}

	// 将长度左移 4 位，然后用 bitwise OR 合并 opCode
	return (len << 4) | opCode;
}

void intToCigar(CigarUnit cigar, char& operation, uint64_t& len)
{
	uint64_t opCode = cigar & 0xF;  // 提取低 4 位为操作编码
	len = cigar >> 4;               // 高 28 位表示操作长度

	// 将操作编码转换回对应字符
	switch (opCode) {
	case 0x0: operation = 'M'; break;
	case 0x1: operation = 'I'; break;
	case 0x2: operation = 'D'; break;
	case 0x7: operation = '='; break;
	case 0x8: operation = 'X'; break;
	default:  operation = '?'; break; // 未知操作标记为 '?'
	}
}

void appendCigar(Cigar_t& dst, const Cigar_t& src)
{
	if (src.empty()) return;                  // nothing to do
	size_t idx = 0;

	//// 1) 检查“拼接点”——dst 的最后一个元素 vs. src 的首元素
	if (!dst.empty()) {
		char op_dst; uint64_t len_dst;
		intToCigar(dst.back(), op_dst, len_dst);

		char op_src; uint64_t len_src;
		intToCigar(src.front(), op_src, len_src);

		// 若操作码相同，则二者合并成一个
		if (op_dst == op_src) {
			dst.back() = cigarToInt(op_dst, len_dst + len_src);
			idx = 1;                         // src 第 0 元素已被合并，后续从 1 开始
		}
	}
	// 2) 将 src 剩余元素逐个 push
	for (; idx < src.size(); ++idx)
		dst.push_back(src[idx]);
}

void prependCigar(Cigar_t& dst, const Cigar_t& src)
{
    if (src.empty()) return; // nothing to do
    size_t end = src.size();

    if (!dst.empty()) {
        char op_dst; uint64_t len_dst;
        intToCigar(dst.front(), op_dst, len_dst);

        char op_src; uint64_t len_src;
        intToCigar(src.back(), op_src, len_src);

        if (op_dst == op_src) {
            // 合并 src.back() 和 dst.front()
            dst.front() = cigarToInt(op_dst, len_dst + len_src);
            end = src.size() - 1; // 最后一个元素已被合并，不再插入
        }
    }

    // 把 src[0 ... end-1] 插到 dst 前面
    dst.insert(dst.begin(), src.begin(), src.begin() + end);
}

void appendCigarOp(Cigar_t& dst, char op, uint64_t len)
{
	if (len == 0) return;
	if (!dst.empty()) {
		char op_last; uint64_t len_last;
		intToCigar(dst.back(), op_last, len_last);
		if (op_last == op) {                        // 合并
			dst.back() = cigarToInt(op, len_last + len);
			return;
		}
	}
	dst.push_back(cigarToInt(op, len));
}
namespace {
constexpr size_t kMaximumRetainedKswEncoding = 64 * 1024;

struct KswEncodingScratch {
    std::vector<uint8_t> reference;
    std::vector<uint8_t> query;
};

thread_local KswEncodingScratch ksw_encoding_scratch;

char orientedBase(const KswSequenceView sequence, size_t index) {
    if (!sequence.reverse_complement) {
        return sequence.bases[index];
    }
    const size_t reversed_index = sequence.bases.size() - index - 1;
    return BASE_COMPLEMENT[static_cast<uint8_t>(
        sequence.bases[reversed_index])];
}

void encodeSequence(const KswSequenceView sequence,
                    std::vector<uint8_t>& encoded) {
    encoded.resize(sequence.bases.size());
    for (size_t index = 0; index < sequence.bases.size(); ++index) {
        encoded[index] = ScoreChar2Idx[static_cast<uint8_t>(
            orientedBase(sequence, index))];
    }
}

bool isCanonicalExactMatch(const KswSequenceView reference,
                           const KswSequenceView query) {
    if (reference.bases.empty() ||
        reference.bases.size() != query.bases.size()) {
        return false;
    }
    for (size_t index = 0; index < reference.bases.size(); ++index) {
        const char ref_base = orientedBase(reference, index);
        const char query_base = orientedBase(query, index);
        if (ref_base != query_base ||
            (ref_base != 'A' && ref_base != 'C' &&
             ref_base != 'G' && ref_base != 'T')) {
            return false;
        }
    }
    return true;
}

std::vector<uint8_t>& selectEncodingBuffer(
    size_t size, std::vector<uint8_t>& retained,
    std::vector<uint8_t>& temporary) {
    return size <= kMaximumRetainedKswEncoding ? retained : temporary;
}

AlignmentResult copyKswResultAndSummarize(const ksw_extz_t& result) {
    AlignmentResult output;
    output.cigar.reserve(result.n_cigar);
    for (int index = 0; index < result.n_cigar; ++index) {
        const CigarUnit unit = result.cigar[index];
        output.cigar.push_back(unit);
        const uint64_t length = unit >> 4;
        const uint32_t operation = unit & 0xf;
        output.summary.alignment_length += length;
        if (operation != 1) output.summary.reference_length += length;
        if (operation != 2) output.summary.query_length += length;
        if (operation == 0) output.summary.match_length += length;
    }
    return output;
}

struct GlobalKswRun {
    Cigar_t cigar;
    CigarSummary summary;
    int score = std::numeric_limits<int>::min();
};

int maximumSubstitutionScore(const int8_t* matrix, int alphabet_size) {
    int maximum = std::numeric_limits<int>::min();
    for (int index = 0; index < alphabet_size * alphabet_size; ++index) {
        maximum = std::max(maximum, static_cast<int>(matrix[index]));
    }
    return maximum;
}

void appendSimpleGlobalCigar(
    Cigar_t& cigar,
    size_t prefix_length,
    size_t gap_length,
    char gap_operation,
    size_t suffix_length) {
    cigar.clear();
    cigar.reserve(3);
    if (prefix_length > 0) {
        appendCigarOp(
            cigar, 'M', static_cast<uint64_t>(prefix_length));
    }
    if (gap_length > 0) {
        appendCigarOp(
            cigar, gap_operation, static_cast<uint64_t>(gap_length));
    }
    if (suffix_length > 0) {
        appendCigarOp(
            cigar, 'M', static_cast<uint64_t>(suffix_length));
    }
}

// Accept a direct/no-extra-gap path only when it strictly beats an upper
// bound for every path with another gap run. Strict comparison plus a unique
// one-gap placement preserves KSW2's tie-breaking contract.
bool tryProvablyOptimalSimpleGlobalAlignment(
    KswSequenceView reference,
    KswSequenceView query,
    const int8_t* matrix,
    int maximum_substitution_score,
    int gap_open,
    int gap_extend,
    Cigar_t& cigar) {
    const size_t reference_length = reference.bases.size();
    const size_t query_length = query.bases.size();
    const size_t shorter_length =
        std::min(reference_length, query_length);
    const size_t length_difference =
        reference_length > query_length
            ? reference_length - query_length
            : query_length - reference_length;

    if (shorter_length == 0) {
        appendSimpleGlobalCigar(
            cigar, 0, length_difference,
            reference_length > query_length ? 'D' : 'I', 0);
        return true;
    }
    if (gap_open < 0 || gap_extend < 0 ||
        maximum_substitution_score < 0) {
        return false;
    }

    const auto score_pair = [matrix](char reference_base, char query_base) {
        return static_cast<int64_t>(
            matrix[
                ScoreChar2Idx[static_cast<uint8_t>(reference_base)] * 5 +
                ScoreChar2Idx[static_cast<uint8_t>(query_base)]]);
    };

    size_t best_gap_position = 0;
    int64_t best_score = 0;
    bool best_gap_position_is_unique = true;
    if (reference_length == query_length) {
        for (size_t index = 0; index < shorter_length; ++index) {
            best_score += score_pair(
                orientedBase(reference, index),
                orientedBase(query, index));
        }
    } else if (reference_length > query_length) {
        int64_t current_score = 0;
        for (size_t index = 0; index < shorter_length; ++index) {
            current_score += score_pair(
                orientedBase(reference, index + length_difference),
                orientedBase(query, index));
        }
        best_score = current_score;
        for (size_t gap_position = 0;
             gap_position < shorter_length;
             ++gap_position) {
            current_score -= score_pair(
                orientedBase(reference, gap_position + length_difference),
                orientedBase(query, gap_position));
            current_score +=
                score_pair(orientedBase(reference, gap_position),
                           orientedBase(query, gap_position));
            if (current_score > best_score) {
                best_score = current_score;
                best_gap_position = gap_position + 1;
                best_gap_position_is_unique = true;
            } else if (current_score == best_score) {
                best_gap_position_is_unique = false;
            }
        }
    } else {
        int64_t current_score = 0;
        for (size_t index = 0; index < shorter_length; ++index) {
            current_score += score_pair(
                orientedBase(reference, index),
                orientedBase(query, index + length_difference));
        }
        best_score = current_score;
        for (size_t gap_position = 0;
             gap_position < shorter_length;
             ++gap_position) {
            current_score -= score_pair(
                orientedBase(reference, gap_position),
                orientedBase(query, gap_position + length_difference));
            current_score +=
                score_pair(orientedBase(reference, gap_position),
                           orientedBase(query, gap_position));
            if (current_score > best_score) {
                best_score = current_score;
                best_gap_position = gap_position + 1;
                best_gap_position_is_unique = true;
            } else if (current_score == best_score) {
                best_gap_position_is_unique = false;
            }
        }
    }

    if (length_difference > 0) {
        best_score -=
            static_cast<int64_t>(gap_open) +
            static_cast<int64_t>(gap_extend) * length_difference;
    }

    if (length_difference > 0 && !best_gap_position_is_unique) {
        return false;
    }

    int64_t excluded_path_upper_bound;
    if (length_difference == 0) {
        excluded_path_upper_bound =
            static_cast<int64_t>(maximum_substitution_score) *
                (shorter_length - 1) -
            2LL * gap_open - 2LL * gap_extend;
    } else {
        excluded_path_upper_bound =
            static_cast<int64_t>(maximum_substitution_score) *
                shorter_length -
            static_cast<int64_t>(gap_extend) * length_difference -
            2LL * gap_open;
    }
    if (best_score <= excluded_path_upper_bound) {
        return false;
    }

    appendSimpleGlobalCigar(
        cigar,
        best_gap_position,
        length_difference,
        reference_length > query_length ? 'D' : 'I',
        shorter_length - best_gap_position);
    return true;
}

GlobalKswRun runGlobalKsw(
    const std::vector<uint8_t>& reference,
    const std::vector<uint8_t>& query,
    const KSW2AlignConfig& config,
    int band_width) {
    ksw_extz_t result{};
    {
    KswCallTimer timer(true);
    ksw_extz2_sse(
        nullptr,
        static_cast<int>(query.size()), query.data(),
        static_cast<int>(reference.size()), reference.data(),
        config.alphabet_size, config.mat,
        config.gap_open, config.gap_extend,
        band_width, config.zdrop, config.end_bonus,
        config.flag, &result);
    }

    const std::unique_ptr<std::uint32_t, decltype(&std::free)> cigar_owner(result.cigar, &std::free);
    GlobalKswRun run;
    run.score = result.score;
    run.cigar.reserve(result.n_cigar);
    for (int index = 0; index < result.n_cigar; ++index) {
        const CigarUnit unit = result.cigar[index];
        run.cigar.push_back(unit);
        const uint64_t length = unit >> 4;
        const uint32_t operation = unit & 0xf;
        run.summary.alignment_length += length;
        if (operation != 1) run.summary.reference_length += length;
        if (operation != 2) run.summary.query_length += length;
        if (operation == 0) run.summary.match_length += length;
    }
    // cigar_owner releases KSW2 storage on success and allocation/validation failure.
    return run;
}

// A path that leaves a width-w diagonal band and finishes at length
// difference d <= w needs at least 2(w + 1) - d gap bases and two gap opens.
// Give every remaining aligned pair the best matrix score; a strict win over
// that upper bound proves every global optimum stays inside the band.
bool bandExcludesEveryGlobalOptimum(
    size_t reference_length,
    size_t query_length,
    int band_width,
    int score,
    int maximum_substitution_score,
    int gap_open,
    int gap_extend) {
    if (band_width < 0) {
        return true;
    }
    const int64_t length_difference = std::abs(
        static_cast<int64_t>(query_length) -
        static_cast<int64_t>(reference_length));
    if (band_width < length_difference) {
        return false;
    }

    const int64_t minimum_gap_bases =
        2LL * (static_cast<int64_t>(band_width) + 1) -
        length_difference;
    const int64_t total_length =
        static_cast<int64_t>(reference_length) +
        static_cast<int64_t>(query_length);
    if (minimum_gap_bases > total_length) {
        return true;
    }

    const int64_t outside_score_upper_bound_times_two =
        static_cast<int64_t>(maximum_substitution_score) *
            (total_length - minimum_gap_bases) -
        2LL * gap_extend * minimum_gap_bases -
        4LL * gap_open;
    return 2LL * score > outside_score_upper_bound_times_two;
}

int bandRequiredToCertifyScore(
    size_t reference_length,
    size_t query_length,
    int score,
    int maximum_substitution_score,
    int gap_open,
    int gap_extend) {
    const int64_t total_length =
        static_cast<int64_t>(reference_length) +
        static_cast<int64_t>(query_length);
    const int64_t length_difference = std::abs(
        static_cast<int64_t>(query_length) -
        static_cast<int64_t>(reference_length));
    const int64_t coefficient =
        static_cast<int64_t>(maximum_substitution_score) +
        2LL * gap_extend;
    const int64_t numerator =
        static_cast<int64_t>(maximum_substitution_score) * total_length -
        4LL * gap_open - 2LL * score;
    const int64_t required_gap_bases =
        numerator < 0 ? 0 : numerator / coefficient + 1;
    const int64_t required_band =
        (required_gap_bases + length_difference + 1) / 2 - 1;
    return static_cast<int>(
        std::max(length_difference, required_band));
}

}  // namespace

Cigar_t globalAlignKSW2(
    const std::string& reference,
    const std::string& query) {
    return globalAlignKSW2Result(
        {reference, false}, {query, false}).cigar;
}

AlignmentResult globalAlignKSW2Result(
    KswSequenceView reference,
    KswSequenceView query) {
    checkKswLengths(reference,query);
    init_simd_mat();

    KSW2AlignConfig config;
    config.mat = dna5_simd_mat;
    config.alphabet_size = 5;
    config.gap_open = ksw2GapOpenPenalty();
    config.gap_extend = ksw2GapExtendPenalty();
    config.end_bonus = 0;
    config.zdrop = -1;
    config.band_width = -1;
    config.flag = KSW_EZ_RIGHT | KSW_EZ_GENERIC_SC;

    Cigar_t simple_cigar;
    const int maximum_substitution_score =
        maximumSubstitutionScore(
            config.mat, config.alphabet_size);
    if (tryProvablyOptimalSimpleGlobalAlignment(
            reference, query, config.mat,
            maximum_substitution_score,
            config.gap_open, config.gap_extend,
            simple_cigar)) {
        AlignmentResult result;
        result.cigar = std::move(simple_cigar);
        result.summary = summarizeCigar(result.cigar);
        return result;
    }

    std::vector<uint8_t> temporary_reference;
    std::vector<uint8_t> temporary_query;
    auto& encoded_reference = selectEncodingBuffer(
        reference.bases.size(), ksw_encoding_scratch.reference,
        temporary_reference);
    auto& encoded_query = selectEncodingBuffer(
        query.bases.size(), ksw_encoding_scratch.query,
        temporary_query);
    encodeSequence(reference, encoded_reference);
    encodeSequence(query, encoded_query);

    const size_t maximum_length =
        std::max(reference.bases.size(), query.bases.size());
    const int64_t length_difference = std::abs(
        static_cast<int64_t>(query.bases.size()) -
        static_cast<int64_t>(reference.bases.size()));
    const bool can_certify_band =
        maximum_substitution_score >= 0 &&
        config.gap_open >= 0 &&
        config.gap_extend >= 0 &&
        maximum_substitution_score + 2 * config.gap_extend > 0;
    if (!can_certify_band || maximum_length <= 64) {
        GlobalKswRun run = runGlobalKsw(
            encoded_reference, encoded_query, config, -1);
        return {std::move(run.cigar), run.summary};
    }

    int band_width = static_cast<int>(
        std::max<int64_t>(64, length_difference));
    if (static_cast<size_t>(band_width) >= maximum_length) {
        GlobalKswRun run = runGlobalKsw(
            encoded_reference, encoded_query, config, -1);
        return {std::move(run.cigar), run.summary};
    }

    GlobalKswRun banded = runGlobalKsw(
        encoded_reference, encoded_query, config, band_width);
    const auto consumes_complete_input = [&reference, &query](
                                             const GlobalKswRun& run) {
        return run.summary.reference_length == reference.bases.size() &&
               run.summary.query_length == query.bases.size();
    };
    if (consumes_complete_input(banded) &&
        bandExcludesEveryGlobalOptimum(
            reference.bases.size(), query.bases.size(), band_width,
            banded.score, maximum_substitution_score,
            config.gap_open, config.gap_extend)) {
        return {std::move(banded.cigar), banded.summary};
    }

    const int required_band = bandRequiredToCertifyScore(
        reference.bases.size(), query.bases.size(), banded.score,
        maximum_substitution_score,
        config.gap_open, config.gap_extend);
    if (required_band > band_width &&
        static_cast<size_t>(required_band) < maximum_length) {
        banded = runGlobalKsw(
            encoded_reference, encoded_query, config, required_band);
        if (consumes_complete_input(banded) &&
            bandExcludesEveryGlobalOptimum(
                reference.bases.size(), query.bases.size(), required_band,
                banded.score, maximum_substitution_score,
                config.gap_open, config.gap_extend)) {
            return {std::move(banded.cigar), banded.summary};
        }
    }

    GlobalKswRun run = runGlobalKsw(
        encoded_reference, encoded_query, config, -1);
    return {std::move(run.cigar), run.summary};
}

Cigar_t extendAlignKSW2(const std::string& ref,
    const std::string& query,
    int zdrop)
{
    return extendAlignKSW2Result(
        {ref, false}, {query, false}, zdrop).cigar;
}

CigarSummary summarizeCigar(const Cigar_t& cigar) {
    CigarSummary summary;
    for (const CigarUnit unit : cigar) {
        const uint64_t length = unit >> 4;
        const uint32_t operation = unit & 0xf;
        summary.alignment_length += length;
        if (operation != 1) summary.reference_length += length;
        if (operation != 2) summary.query_length += length;
        if (operation == 0) summary.match_length += length;
    }
    return summary;
}

static AlignmentResult extendAlignKSW2Impl(
    KswSequenceView ref, KswSequenceView query, int zdrop,
    bool allow_shortcuts)
{
    checkKswLengths(ref,query);
    if (allow_shortcuts && ref.bases.empty() && query.bases.empty()) {
        return {};
    }
    if (allow_shortcuts && isCanonicalExactMatch(ref, query)) {
        AlignmentResult result;
        const auto length = static_cast<uint64_t>(ref.bases.size());
        result.cigar.push_back(cigarToInt('M', length));
        result.summary = {length, length, length, length};
        return result;
    }

    /* ---------- 1. 序列编码 ---------- */
    std::vector<uint8_t> temporary_reference;
    std::vector<uint8_t> temporary_query;
    auto& ref_enc = selectEncodingBuffer(
        ref.bases.size(), ksw_encoding_scratch.reference,
        temporary_reference);
    auto& qry_enc = selectEncodingBuffer(
        query.bases.size(), ksw_encoding_scratch.query,
        temporary_query);
    encodeSequence(ref, ref_enc);
    encodeSequence(query, qry_enc);

    ///* ---------- 2. 配置 ---------- */
    //KSW2AlignConfig cfg = makeTurboKSW2Config(query.size(), ref.size());
    ////KSW2AlignConfig cfg;
    //cfg.zdrop = zdrop;       // 用于提前终止
    //cfg.flag = KSW_EZ_EXTZ_ONLY     // ends-free extension
    //    | KSW_EZ_APPROX_MAX    // 跟踪 ez.max_q/max_t
    //    | KSW_EZ_APPROX_DROP   // 在 approximate 模式下触发 z-drop 就中断
    //    | KSW_EZ_RIGHT;        // （可选）gap 右对齐     // **关键**：启用 extension/ends-free
    //// 若需要右对齐 gaps 建议保留 KSW_EZ_RIGHT
    //cfg.end_bonus = 100;
    //cfg.band_width = -1;
    init_simd_mat();
    KSW2AlignConfig cfg;
	cfg.mat = dna5_simd_mat;
    cfg.zdrop = zdrop;
    cfg.flag = KSW_EZ_EXTZ_ONLY | KSW_EZ_RIGHT |
        KSW_EZ_APPROX_DROP | KSW_EZ_GENERIC_SC;
    cfg.end_bonus = 50;
    cfg.alphabet_size = 5;
    cfg.gap_open = ksw2GapOpenPenalty();
    cfg.gap_extend = ksw2GapExtendPenalty();
    cfg.band_width = auto_band(ref.bases.size(), query.bases.size());


    /* ---------- 3. 调用 KSW2 ---------- */
    ksw_extz_t ez{};
    {
    KswCallTimer timer(false);
    ksw_extz2_sse(nullptr,
        static_cast<int>(qry_enc.size()), qry_enc.data(),
        static_cast<int>(ref_enc.size()), ref_enc.data(),
        cfg.alphabet_size, cfg.mat,
        cfg.gap_open, cfg.gap_extend,
        cfg.band_width, cfg.zdrop, cfg.end_bonus,
        cfg.flag, &ez);
    }

    // 赋值bool& if_zdrop,int& ref_end,int& qry_end
    /* ---------- 4. 拷贝 & 释放 ---------- */
    const std::unique_ptr<std::uint32_t, decltype(&std::free)> cigar_owner(ez.cigar, &std::free);
    AlignmentResult result = copyKswResultAndSummarize(ez);
    return result;
}

AlignmentResult extendAlignKSW2Result(
    KswSequenceView ref, KswSequenceView query, int zdrop) {
    return extendAlignKSW2Impl(ref, query, zdrop, true);
}


}

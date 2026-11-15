#include "ramag/internal/centromere/support.hpp"
#include "ramag/internal/centromere/gap_alignment.hpp"
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <limits>
#include <new>
#include <vector>
namespace ramag::detail {
namespace {
using Clock = std::chrono::steady_clock;
using Score = std::int32_t;
using Cell = std::array<Score, 7>;
constexpr Score Neg = std::numeric_limits<Score>::min() / 4;
constexpr std::array<Score, 3> Opens{60, 800, 2500};
constexpr std::array<Score, 3> Extends{30, 5, 1};
using ::ramag::detail::IsCanonicalUpper;
struct PackedTrace {
  std::vector<std::uint64_t> words;
  explicit PackedTrace(std::size_t n) : words(n, 0) {}
  void set(std::uint64_t cell, std::uint16_t code) {
    const auto bit = cell * 9;
    const auto word = static_cast<std::size_t>(bit >> 6U);
    const auto shift = static_cast<unsigned>(bit & 63U);
    words[word] |= static_cast<std::uint64_t>(code) << shift;
    if (shift > 55U)
      words[word + 1] |= static_cast<std::uint64_t>(code) >> (64U - shift);
  }
  std::uint16_t get(std::uint64_t cell) const {
    const auto bit = cell * 9;
    const auto word = static_cast<std::size_t>(bit >> 6U);
    const auto shift = static_cast<unsigned>(bit & 63U);
    auto value = words[word] >> shift;
    if (shift > 55U)
      value |= words[word + 1] << (64U - shift);
    return static_cast<std::uint16_t>(value & 511U);
  }
};
void AppendReverse(std::vector<CigarOp> &out, char op) {
  if (!out.empty() && out.back().operation == op)
    ++out.back().length;
  else
    out.push_back({op, 1});
}
template <bool ExactGaps>
StitchGapResult AlignTwoStripImpl(std::string_view, std::string_view,
                                  const CentromereStitchOptions &,
                                  Clock::time_point, std::uint64_t);
} // namespace
namespace {
template <bool ExactGaps>
StitchGapResult AlignBlockDpImpl(std::string_view reference,
                                std::string_view oriented_query,
                                const CentromereStitchOptions &options) {
  if (options.score_model != CentromereStitchScoreModel::Centrolign3Piece ||
      !options.gap_memory_limit_bytes ||
      !std::isfinite(options.gap_seconds_limit) ||
      options.gap_seconds_limit <= 0)
    throw AlignmentError("invalid linear block DP scoring/resources");
  const auto started = Clock::now();
  StitchGapResult out;
  out.engine = ExactGaps ? "linear_7state_dp_exact_gaps" : "linear_7state_dp";
  auto seconds = [&] {
    return std::chrono::duration<double>(Clock::now() - started).count();
  };
  auto limit = [&](const char *reason) -> void {
    StitchGapResourceLimit error(reason);
    error.engine = out.engine;
    error.dp_cells = out.dp_cells;
    error.peak_kernel_bytes = out.peak_kernel_bytes;
    error.elapsed_seconds = seconds();
    throw error;
  };
  auto check = [&] {
    if (options.interruption_callback)
      options.interruption_callback("centromere-block-linear-dp");
    if (seconds() > options.gap_seconds_limit)
      limit("gap_time_budget");
  };
  check();
  // Orient the shorter sequence on score-row columns. Operation labels are
  // restored at traceback; original forward coordinates never change.
  auto r = reference, q = oriented_query;
  const bool swapped = q.size() > r.size();
  if (swapped)
    std::swap(r, q);
  const auto m = static_cast<std::uint64_t>(r.size());
  const auto n = static_cast<std::uint64_t>(q.size());
  const auto max = std::numeric_limits<std::uint64_t>::max();
  if (m == max || n == max || m + 1 > max / (n + 1))
    limit("kernel_numerical_size_limit");
  const auto cells = (m + 1) * (n + 1);
  out.dp_cells = cells;
  if (m > max - n || m + n > static_cast<std::uint64_t>((-Neg - 5000) / 80))
    limit("kernel_numerical_size_limit");
  if (cells > (max - 63) / 9)
    limit("gap_memory_budget");
  const auto words = (cells * 9 + 63) / 64;
  if (words > max / sizeof(std::uint64_t) || n + 1 > max / (2 * sizeof(Cell)) ||
      m + n > max / sizeof(CigarOp))
    limit("gap_memory_budget");
  const auto trace_bytes = words * sizeof(std::uint64_t);
  const auto row_bytes = (n + 1) * 2 * sizeof(Cell);
  const auto output_bytes = (m + n) * sizeof(CigarOp);
  const auto budget = options.gap_memory_limit_bytes;
  if (trace_bytes > budget || row_bytes > budget - trace_bytes ||
      output_bytes > budget - trace_bytes - row_bytes ||
      words > std::numeric_limits<std::size_t>::max())
    return AlignTwoStripImpl<ExactGaps>(reference, oriented_query, options, started, 0);
  try {
    // Reserve all overlapping kernel-owned storage before starting the DP.
    PackedTrace trace(static_cast<std::size_t>(words));
    std::vector<Cell> prev(static_cast<std::size_t>(n + 1)), curr(prev.size());
    out.cigar.reserve(static_cast<std::size_t>(m + n));
    out.peak_kernel_bytes = trace.words.capacity() * sizeof(std::uint64_t) +
                            (prev.capacity() + curr.capacity()) * sizeof(Cell) +
                            out.cigar.capacity() * sizeof(CigarOp);
    if (out.peak_kernel_bytes > budget)
      limit("gap_memory_budget");
    for (const auto text : {r, q})
      for (std::size_t i = 0; i < text.size(); ++i) {
        if ((i & 32767U) == 0)
          check();
        if (!IsCanonicalUpper(text[i]) && text[i] != 'N')
          throw AlignmentError("linear block DP requires uppercase ACGTN");
      }
    prev[0].fill(Neg);
    prev[0][0] = 0;
    for (std::uint64_t j = 1; j <= n; ++j) {
      auto &here = prev[static_cast<std::size_t>(j)];
      const auto &left = prev[static_cast<std::size_t>(j - 1)];
      here.fill(Neg);
      std::uint16_t code = 0;
      for (std::size_t k = 0; k < 3; ++k) {
        const Score open = left[0] - Opens[k] - Extends[k];
        const Score extend = left[k + 4] - Extends[k];
        here[k + 4] = std::max(open, extend);
        if (extend >= open)
          code |= static_cast<std::uint16_t>(1U << (k + 6));
        if (here[k + 4] > here[0]) {
          here[0] = here[k + 4];
          code = static_cast<std::uint16_t>((code & ~7U) | (k + 4));
        }
      }
      trace.set(j, code);
    }
    std::uint64_t tick = 0;
    for (std::uint64_t i = 1; i <= m; ++i) {
      const auto base = i * (n + 1);
      auto &first = curr[0];
      first.fill(Neg);
      std::uint16_t code = 0;
      for (std::size_t k = 0; k < 3; ++k) {
        const Score open = prev[0][0] - Opens[k] - Extends[k];
        const Score extend = prev[0][k + 1] - Extends[k];
        first[k + 1] = std::max(open, extend);
        if (extend >= open)
          code |= static_cast<std::uint16_t>(1U << (k + 3));
        if (first[k + 1] > first[0]) {
          first[0] = first[k + 1];
          code = static_cast<std::uint16_t>((code & ~7U) | (k + 1));
        }
      }
      trace.set(base, code);
      const char rb = r[static_cast<std::size_t>(i - 1)];
      for (std::uint64_t j = 1; j <= n; ++j) {
        if ((++tick & 32767U) == 0)
          check();
        auto &here = curr[static_cast<std::size_t>(j)];
        const auto &up = prev[static_cast<std::size_t>(j)];
        const auto &left = curr[static_cast<std::size_t>(j - 1)];
        const char qb = q[static_cast<std::size_t>(j - 1)];
        const bool exact = rb == qb && IsCanonicalUpper(rb);
        // The explicit exact-gap candidate searches only =/I/D paths.
        // N remains a hard exact-match break, including N against N.
        here[0] = ExactGaps && !exact
                      ? Neg
                      : prev[static_cast<std::size_t>(j - 1)][0] +
                            (exact ? 20 : -80);
        code = 0;
        for (std::size_t k = 0; k < 3; ++k) {
          const Score open = up[0] - Opens[k] - Extends[k];
          const Score extend = up[k + 1] - Extends[k];
          here[k + 1] = std::max(open, extend);
          if (extend >= open)
            code |= static_cast<std::uint16_t>(1U << (k + 3));
          if (here[k + 1] > here[0]) {
            here[0] = here[k + 1];
            code = static_cast<std::uint16_t>((code & ~7U) | (k + 1));
          }
        }
        for (std::size_t k = 0; k < 3; ++k) {
          const Score open = left[0] - Opens[k] - Extends[k];
          const Score extend = left[k + 4] - Extends[k];
          here[k + 4] = std::max(open, extend);
          if (extend >= open)
            code |= static_cast<std::uint16_t>(1U << (k + 6));
          if (here[k + 4] > here[0]) {
            here[0] = here[k + 4];
            code = static_cast<std::uint16_t>((code & ~7U) | (k + 4));
          }
        }
        trace.set(base + j, code);
      }
      std::swap(prev, curr);
    }
    check();
    out.score = prev[static_cast<std::size_t>(n)][0];
    auto i = m, j = n;
    unsigned state = 0;
    std::uint64_t moves = 0;
    while (i || j) {
      if ((++moves & 32767U) == 0)
        check();
      const auto code = trace.get(i * (n + 1) + j);
      if (state == 0) {
        state = code & 7U;
        if (state)
          continue;
        if (!i || !j)
          throw AlignmentError("linear block DP diagonal traceback bounds");
        const char a = r[static_cast<std::size_t>(i - 1)],
                   b = q[static_cast<std::size_t>(j - 1)];
        AppendReverse(out.cigar, a == b && IsCanonicalUpper(a) ? '=' : 'X');
        --i;
        --j;
      } else {
        if (state > 6)
          throw AlignmentError("linear block DP traceback state");
        const bool extend = (code & (1U << (state + 2))) != 0;
        if (state <= 3) {
          if (!i)
            throw AlignmentError("linear block DP deletion traceback bounds");
          AppendReverse(out.cigar, swapped ? 'I' : 'D');
          --i;
        } else {
          if (!j)
            throw AlignmentError("linear block DP insertion traceback bounds");
          AppendReverse(out.cigar, swapped ? 'D' : 'I');
          --j;
        }
        if (!extend)
          state = 0;
      }
    }
    std::reverse(out.cigar.begin(), out.cigar.end());
    if (ScoreCentromereStitchCigar(reference, oriented_query, out.cigar,
                                   options.score_model) != out.score)
      throw AlignmentError("linear block DP score disagrees with CIGAR replay");
    check();
    out.elapsed_seconds = seconds();
    return out;
  } catch (const std::bad_alloc &) {
    limit("allocation_failure");
  }
  throw AlignmentError("unreachable linear block DP control flow");
}

} // namespace
StitchGapResult AlignCentromereBlockDp(std::string_view reference,
                                      std::string_view oriented_query,
                                      const CentromereStitchOptions &options) {
  return AlignBlockDpImpl<false>(reference, oriented_query, options);
}
StitchGapResult AlignCentromereBlockDpExactGaps(
    std::string_view reference, std::string_view oriented_query,
    const CentromereStitchOptions &options) {
  return AlignBlockDpImpl<true>(reference, oriented_query, options);
}
namespace {
template <bool ExactGaps>
StitchGapResult AlignTwoStripImpl(std::string_view reference,
                                  std::string_view oriented_query,
                                  const CentromereStitchOptions &options,
                                  Clock::time_point started,
                                  std::uint64_t forced_suffix_rows) {
  if (options.score_model != CentromereStitchScoreModel::Centrolign3Piece ||
      !options.gap_memory_limit_bytes ||
      !std::isfinite(options.gap_seconds_limit) ||
      options.gap_seconds_limit <= 0)
    throw AlignmentError("invalid two-strip block DP scoring/resources");
  StitchGapResult out;
  out.engine = ExactGaps ? "linear_7state_dp_2stripe_exact_gaps"
                         : "linear_7state_dp_2stripe";
  auto seconds = [&] {
    return std::chrono::duration<double>(Clock::now() - started).count();
  };
  auto limit = [&](const char *reason) -> void {
    StitchGapResourceLimit error(reason);
    error.engine = out.engine;
    error.dp_cells = out.dp_cells;
    error.peak_kernel_bytes = out.peak_kernel_bytes;
    error.elapsed_seconds = seconds();
    throw error;
  };
  auto check = [&] {
    if (options.interruption_callback)
      options.interruption_callback("centromere-block-two-strip-dp");
    if (seconds() > options.gap_seconds_limit)
      limit("gap_time_budget");
  };
  check();
  auto r = reference, q = oriented_query;
  const bool swapped = q.size() > r.size();
  if (swapped)
    std::swap(r, q);
  const auto m = static_cast<std::uint64_t>(r.size());
  const auto n = static_cast<std::uint64_t>(q.size());
  const auto max = std::numeric_limits<std::uint64_t>::max();
  if (m == max || n == max || m + 1 > max / (n + 1))
    limit("kernel_numerical_size_limit");
  const auto cells = (m + 1) * (n + 1);
  out.dp_cells = cells;
  if (m > max - n || m + n > static_cast<std::uint64_t>((-Neg - 5000) / 80))
    limit("kernel_numerical_size_limit");
  if (cells > (max - 63) / 9 || n + 1 > max / (3 * sizeof(Cell)) ||
      m + n > max / sizeof(CigarOp))
    limit("gap_memory_budget");
  constexpr std::uint64_t Reserve = 16ULL * 1024ULL * 1024ULL;
  const auto budget = options.gap_memory_limit_bytes;
  const auto column_count = n + 1;
  const auto score_bytes = column_count * 3 * sizeof(Cell);
  const auto output_bytes = (m + n) * sizeof(CigarOp);
  if (m < 2 || budget < Reserve || score_bytes > budget - Reserve ||
      output_bytes > budget - Reserve - score_bytes)
    limit("gap_memory_budget");
  const auto trace_budget = budget - Reserve - score_bytes - output_bytes;
  auto words_for = [&](std::uint64_t height) {
    return (((height + 1) * column_count) * 9 + 63) / 64;
  };
  auto fits = [&](std::uint64_t height) {
    const auto words = words_for(height);
    return words <= trace_budget / sizeof(std::uint64_t) &&
           words <= std::numeric_limits<std::size_t>::max();
  };
  std::uint64_t height = forced_suffix_rows;
  if (height) {
    if (height >= m || m - height > height)
      throw AlignmentError("invalid forced two-strip height");
    if (!fits(height))
      limit("gap_memory_budget");
  } else {
    std::uint64_t low = 0, high = m - 1;
    while (low < high) {
      const auto mid = low + (high - low + 1) / 2;
      if (fits(mid))
        low = mid;
      else
        high = mid - 1;
    }
    height = low;
    if (!height || m - height > height)
      limit("gap_memory_budget");
  }
  const auto boundary = m - height;
  const auto trace_words = words_for(height);
  try {
    // One trace allocation survives both phases. The shorter prefix is
    // recomputed into it only after suffix traceback has finished.
    PackedTrace trace(static_cast<std::size_t>(trace_words));
    std::vector<Cell> prev(static_cast<std::size_t>(column_count)),
        curr(prev.size()), checkpoint(prev.size());
    out.cigar.reserve(static_cast<std::size_t>(m + n));
    out.peak_kernel_bytes =
        trace.words.capacity() * sizeof(std::uint64_t) +
        (prev.capacity() + curr.capacity() + checkpoint.capacity()) *
            sizeof(Cell) +
        out.cigar.capacity() * sizeof(CigarOp);
    if (out.peak_kernel_bytes > budget - Reserve)
      limit("gap_memory_budget");
    check();
    for (const auto text : {r, q})
      for (std::size_t i = 0; i < text.size(); ++i) {
        if ((i & 32767U) == 0)
          check();
        if (!IsCanonicalUpper(text[i]) && text[i] != 'N')
          throw AlignmentError("two-strip block DP requires uppercase ACGTN");
      }
    auto initialize = [&]<bool Store>(std::uint64_t last_column) {
      prev[0].fill(Neg);
      prev[0][0] = 0;
      for (std::uint64_t j = 1; j <= last_column; ++j) {
        auto &here = prev[static_cast<std::size_t>(j)];
        const auto &left = prev[static_cast<std::size_t>(j - 1)];
        here.fill(Neg);
        std::uint16_t code = 0;
        for (std::size_t k = 0; k < 3; ++k) {
          const Score open = left[0] - Opens[k] - Extends[k];
          const Score extend = left[k + 4] - Extends[k];
          here[k + 4] = std::max(open, extend);
          if (extend >= open)
            code |= static_cast<std::uint16_t>(1U << (k + 6));
          if (here[k + 4] > here[0]) {
            here[0] = here[k + 4];
            code = static_cast<std::uint16_t>((code & ~7U) | (k + 4));
          }
        }
        if constexpr (Store)
          trace.set(j, code);
      }
    };
    std::uint64_t tick = 0;
    auto row = [&]<bool Store>(std::uint64_t i, std::uint64_t origin,
                               std::uint64_t last_column,
                               std::uint64_t stride) {
      const auto base = (i - origin) * stride;
      auto &first = curr[0];
      first.fill(Neg);
      std::uint16_t code = 0;
      for (std::size_t k = 0; k < 3; ++k) {
        const Score open = prev[0][0] - Opens[k] - Extends[k];
        const Score extend = prev[0][k + 1] - Extends[k];
        first[k + 1] = std::max(open, extend);
        if (extend >= open)
          code |= static_cast<std::uint16_t>(1U << (k + 3));
        if (first[k + 1] > first[0]) {
          first[0] = first[k + 1];
          code = static_cast<std::uint16_t>((code & ~7U) | (k + 1));
        }
      }
      if constexpr (Store)
        trace.set(base, code);
      const char rb = r[static_cast<std::size_t>(i - 1)];
      for (std::uint64_t j = 1; j <= last_column; ++j) {
        if ((++tick & 32767U) == 0)
          check();
        auto &here = curr[static_cast<std::size_t>(j)];
        const auto &up = prev[static_cast<std::size_t>(j)];
        const auto &left = curr[static_cast<std::size_t>(j - 1)];
        const char qb = q[static_cast<std::size_t>(j - 1)];
        const bool exact = rb == qb && IsCanonicalUpper(rb);
        // The explicit exact-gap candidate searches only =/I/D paths.
        // N remains a hard exact-match break, including N against N.
        here[0] = ExactGaps && !exact
                      ? Neg
                      : prev[static_cast<std::size_t>(j - 1)][0] +
                            (exact ? 20 : -80);
        code = 0;
        for (std::size_t k = 0; k < 3; ++k) {
          const Score open = up[0] - Opens[k] - Extends[k];
          const Score extend = up[k + 1] - Extends[k];
          here[k + 1] = std::max(open, extend);
          if (extend >= open)
            code |= static_cast<std::uint16_t>(1U << (k + 3));
          if (here[k + 1] > here[0]) {
            here[0] = here[k + 1];
            code = static_cast<std::uint16_t>((code & ~7U) | (k + 1));
          }
        }
        for (std::size_t k = 0; k < 3; ++k) {
          const Score open = left[0] - Opens[k] - Extends[k];
          const Score extend = left[k + 4] - Extends[k];
          here[k + 4] = std::max(open, extend);
          if (extend >= open)
            code |= static_cast<std::uint16_t>(1U << (k + 6));
          if (here[k + 4] > here[0]) {
            here[0] = here[k + 4];
            code = static_cast<std::uint16_t>((code & ~7U) | (k + 4));
          }
        }
        if constexpr (Store)
          trace.set(base + j, code);
      }
      std::swap(prev, curr);
    };
    initialize.template operator()<false>(n);
    for (std::uint64_t i = 1; i <= boundary; ++i)
      row.template operator()<false>(i, 0, n, column_count);
    std::copy(prev.begin(), prev.end(), checkpoint.begin());
    for (std::uint64_t i = boundary + 1; i <= m; ++i)
      row.template operator()<true>(i, boundary, n, column_count);
    check();
    out.score = prev[static_cast<std::size_t>(n)][0];
    auto i = m, j = n;
    unsigned state = 0;
    std::uint64_t moves = 0;
    auto traceback = [&](std::uint64_t origin, std::uint64_t stride) {
      // Do not consume the boundary row. Its horizontal gaps and any
      // carried vertical component belong to the recomputed prefix trace.
      while (i > origin || (origin == 0 && j)) {
        if ((++moves & 32767U) == 0)
          check();
        const auto code = trace.get((i - origin) * stride + j);
        if (state == 0) {
          state = code & 7U;
          if (state)
            continue;
          if (!i || !j)
            throw AlignmentError("two-strip diagonal traceback bounds");
          const char a = r[static_cast<std::size_t>(i - 1)],
                     b = q[static_cast<std::size_t>(j - 1)];
          AppendReverse(out.cigar, a == b && IsCanonicalUpper(a) ? '=' : 'X');
          --i;
          --j;
        } else {
          if (state > 6)
            throw AlignmentError("two-strip traceback state");
          const bool extend = (code & (1U << (state + 2))) != 0;
          if (state <= 3) {
            if (!i)
              throw AlignmentError("two-strip deletion traceback bounds");
            AppendReverse(out.cigar, swapped ? 'I' : 'D');
            --i;
          } else {
            if (!j)
              throw AlignmentError("two-strip insertion traceback bounds");
            AppendReverse(out.cigar, swapped ? 'D' : 'I');
            --j;
          }
          if (!extend)
            state = 0;
        }
      }
    };
    traceback(boundary, column_count);
    check();
    // The carried state at (boundary,j) depends only on the upper-left
    // rectangle. Recompute no columns beyond this known traceback endpoint;
    // their scores cannot influence any predecessor of the carried state.
    const auto prefix_last_column = j;
    const auto prefix_columns = prefix_last_column + 1;
    const auto prefix_words = ((boundary + 1) * prefix_columns * 9 + 63) / 64;
    std::fill(trace.words.begin(),
              trace.words.begin() + static_cast<std::ptrdiff_t>(prefix_words),
              0);
    check();
    initialize.template operator()<true>(prefix_last_column);
    for (std::uint64_t k = 1; k <= boundary; ++k)
      row.template operator()<true>(k, 0, prefix_last_column, prefix_columns);
    if (!std::equal(prev.begin(),
                    prev.begin() + static_cast<std::ptrdiff_t>(prefix_columns),
                    checkpoint.begin()))
      throw AlignmentError("two-strip checkpoint scores differ on replay");
    traceback(0, prefix_columns);
    std::reverse(out.cigar.begin(), out.cigar.end());
    if (ScoreCentromereStitchCigar(reference, oriented_query, out.cigar,
                                   options.score_model) != out.score)
      throw AlignmentError("two-strip DP score disagrees with CIGAR replay");
    check();
    out.elapsed_seconds = seconds();
    return out;
  } catch (const std::bad_alloc &) {
    limit("allocation_failure");
  }
  throw AlignmentError("unreachable two-strip DP control flow");
}
} // namespace
StitchGapResult
AlignCentromereBlockDpTwoStrips(std::string_view reference,
                                std::string_view oriented_query,
                                const CentromereStitchOptions &options,
                                std::uint64_t retained_suffix_rows) {
  if (!retained_suffix_rows)
    throw AlignmentError("forced two-strip height must be positive");
  return AlignTwoStripImpl<false>(reference, oriented_query, options, Clock::now(),
                                  retained_suffix_rows);
}
StitchGapResult AlignCentromereBlockDpExactGapsTwoStrips(
    std::string_view reference, std::string_view oriented_query,
    const CentromereStitchOptions &options, std::uint64_t retained_suffix_rows) {
  if (!retained_suffix_rows)
    throw AlignmentError("forced exact-gap two-strip height must be positive");
  return AlignTwoStripImpl<true>(reference, oriented_query, options, Clock::now(),
                                 retained_suffix_rows);
}
} // namespace ramag::detail

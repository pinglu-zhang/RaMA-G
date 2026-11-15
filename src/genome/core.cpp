#include "ramag/internal/genome.hpp"
// Pairwise alignment kernel. Source attribution and MIT notice are retained
// in LICENSE. Do not substitute unrelated gap backends.
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
namespace ramag::pairwise_detail {
thread_local PairwiseStatistics* active_statistics = nullptr;
UnionFind::UnionFind(std::size_t n) {
    reset(n);
}

// 重置并查集为 n 个独立集合
void UnionFind::reset(std::size_t n) {
    parent_.assign(n, -1);                      // -1 表示单元素集合（根节点，大小=1）
    component_cnt_ = static_cast<int_t>(n);     // 初始有 n 个连通分量
}

// ────────────────────────────────────────────
// Query（查询接口）
// ────────────────────────────────────────────

// 查找：返回元素 x 所在集合的根
int_t UnionFind::find(int_t x) {
    // 路径折半（迭代版）：不断把 x 挂到祖父节点上，加速后续查询
    while (parent_[x] >= 0 && parent_[parent_[x]] >= 0) {
        parent_[x] = parent_[parent_[x]];
        x = parent_[x];
    }
    return (parent_[x] < 0) ? x : parent_[x];
}

// 返回 x 所在集合大小
int_t UnionFind::set_size(int_t x) {
    return -parent_[find(x)];
}

// 判断 a 和 b 是否属于同一集合
bool UnionFind::same(int_t a, int_t b) {
    return find(a) == find(b);
}

// ────────────────────────────────────────────
// Modification（修改接口）
// ────────────────────────────────────────────

// 合并两个集合（按大小合并）
// 返回 true 表示确实发生了合并；false 表示原本就在同一集合
bool UnionFind::unite(int_t a, int_t b) {
    a = find(a);
    b = find(b);
    if (a == b) return false;

    // 按大小合并：parent_ 为负表示大小，数值越小(更负)集合越大
    // 若 a 的集合更小，则交换，让 a 做大根
    if (parent_[a] > parent_[b]) std::swap(a, b);

    parent_[a] += parent_[b];   // 更新根 a 的大小（负值累加）
    parent_[b] = a;             // b 挂到 a
    --component_cnt_;
    return true;
}

// ────────────────────────────────────────────
// filterAndMergeMatches：对 match 列表进行过滤与合并
// 假设：matches 已按 query_begin 排序
// 规则：
// 1) 同一 diagonal 且同向同染色体：合并成更长 match，删除后者
// 2) 同 ref 起点：根据重叠比例过滤较短者；同长时用 tentative 机制处理
// 3) 同 qry 起点：同上
// ────────────────────────────────────────────
void filterAndMergeMatches(std::span<Match>& matches) {
    if (matches.empty()) return;

    const size_t N = matches.size();
    std::vector<bool> good(N, true);        // good[i]=true 表示保留
    std::vector<bool> tentative(N, false);  // tentative 标记：用于“长度相等且重叠较大”的歧义处理

    for (size_t i = 0; i < N; ++i) {
        if (!good[i]) continue;

        const Match& mi = matches[i];
        int_t i_diag = static_cast<int_t>(mi.query_begin-mi.reference_begin);
        Coord_t i_end = mi.query_begin + mi.length;  // i 在 query 维度的结束位置（开区间右端）

        // 由于 matches 按 query_begin 排序，只需检查 query_begin <= i_end 的后续元素
        for (size_t j = i + 1; j < N && matches[j].query_begin <= i_end; ++j) {
            if (!good[j]) continue;

            const Match& mj = matches[j];
            int_t j_diag = static_cast<int_t>(mj.query_begin)-static_cast<int_t>(mj.reference_begin);

            // --- Case 1: 同一 diagonal（同一条对角线） ---
            // 条件：diag 相等 + strand 相等 + ref/qry 染色体一致
            if (i_diag == j_diag &&
                mi.strand == mj.strand &&
                mi.reference_id == mj.reference_id &&
                mi.query_id == mj.query_id) {

                // 合并为更长的：计算 mj 相对 mi 的延伸长度
                Coord_t j_extent = mj.length + mj.query_begin - mi.query_begin;
                if (j_extent > matches[i].length) {
                    matches[i].length = j_extent;
                    i_end = mi.query_begin + j_extent;
                }
                good[j] = false; // 删除 mj
            }

            // --- Case 2: 同一 ref 起点 ---
            else if (mi.reference_begin == mj.reference_begin &&
                     mi.reference_id == mj.reference_id) {

                int_t overlap = mi.query_begin + mi.length - mj.query_begin;

                if (mi.length < mj.length) {
                    // i 更短：若重叠超过 i 一半，丢弃 i
                    if (overlap >= static_cast<int_t>(mi.length / 2)) {
                        good[i] = false;
                        break;
                    }
                }
                else if (mj.length < mi.length) {
                    // j 更短：若重叠超过 j 一半，丢弃 j
                    if (overlap >= static_cast<int_t>(mj.length / 2))
                        good[j] = false;
                }
                else {
                    // 长度相等：若重叠超过一半，标记 tentative
                    if (overlap >= static_cast<int_t>(mi.length / 2)) {
                        tentative[j] = true;
                        if (tentative[i]) {
                            // 如果 i 也已被 tentative，则丢弃 i
                            good[i] = false;
                            break;
                        }
                    }
                }
            }

            // --- Case 3: 同一 qry 起点 ---
            else if (mi.query_begin == mj.query_begin &&
                     mi.query_id == mj.query_id) {

                int64_t overlap = static_cast<int64_t>(mi.reference_begin) +
                                  static_cast<int64_t>(mi.length) -
                                  static_cast<int64_t>(mj.reference_begin);

                if (mi.length < mj.length) {
                    if (overlap >= static_cast<int64_t>(mi.length / 2)) {
                        good[i] = false;
                        break;
                    }
                }
                else if (mj.length < mi.length) {
                    if (overlap >= static_cast<int64_t>(mj.length / 2))
                        good[j] = false;
                }
                else {
                    if (overlap >= static_cast<int64_t>(mi.length / 2)) {
                        tentative[j] = true;
                        if (tentative[i]) {
                            good[i] = false;
                            break;
                        }
                    }
                }
            }
        }
    }

    // 收集 good 的 matches，生成新数组
    size_t retained = 0;
    for (size_t i = 0; i < N; ++i) {
        if (good[i]) matches[retained++] = matches[i];
    }
    matches = matches.first(retained);
}

// ────────────────────────────────────────────
// buildClusters：将 unique_match 根据 max_gap / diagdiff / diagfactor 聚类
// 逻辑：
// 1) 对 match 排序（按 start2，再按 start1）
// 2) filterAndMergeMatches 做合并压缩
// 3) 使用并查集：若两个 match 在 query 维度距离 sep <= max_gap
//    且 diagonal 差 <= max(diagdiff, diagfactor * sep)，则归为一类
// 4) 最后根据并查集根构建簇列表
// ────────────────────────────────────────────
MatchClusterVec buildClusters(std::span<Match> unique_match,
                             int_t  max_gap,
                             int_t  diagdiff,
                             double diagfactor) {
    MatchClusterVec clusters;

    // 特殊情况：0 或 1 个元素，直接返回
    if (unique_match.size() < 2) {
        if (unique_match.size() == 1) {
            clusters.emplace_back();
            clusters.back().push_back(std::move(unique_match[0]));
        }
        return clusters;
    }

    // 判断链方向（假设同一批次一致）
    const bool is_forward = (unique_match.front().strand == FORWARD);

    // 先排序：按 start2 再按 start1
    std::sort(unique_match.begin(), unique_match.end(),
        [](const Match& a, const Match& b) {
            if (start2(a) < start2(b)) return true;
            if (start2(a) > start2(b)) return false;
            return start1(a) < start1(b);
        });

    // 再合并压缩（可能会删除元素、缩短 vector）
    filterAndMergeMatches(unique_match);

    // 重新获取 N
    const uint_t N = static_cast<uint_t>(unique_match.size());
    if (N < 2) {
        if (N == 1) {
            clusters.emplace_back();
            clusters.back().push_back(std::move(unique_match[0]));
        }
        return clusters;
    }

    UnionFind uf(N);

    // 利用排序后的局部性进行聚类：对每个 i，只检查后续 sep 不超过 max_gap 的 j
    for (uint_t i = 0; i < N; ++i) {
        uint_t i_end = start2(unique_match[i]) + len2(unique_match[i]);

        int_t i_diag = 0;
        if (is_forward) {
            i_diag = diag(unique_match[i]);
        }
        else {
            i_diag = diag_reverse(unique_match[i]);
        }

        for (uint_t j = i + 1; j < N; ++j) {
            int_t sep = static_cast<int_t>(start2(unique_match[j])) - static_cast<int_t>(i_end);

            // 早停：gap 太大则后面的 j 也不可能满足
            if (sep > static_cast<int_t>(max_gap)) break;

            uint_t diag_diff = 0;
            if (is_forward) {
                diag_diff = std::abs(diag(unique_match[j]) - i_diag);
            }
            else {
                diag_diff = std::abs(diag_reverse(unique_match[j]) - i_diag);
            }

            // 阈值：max(diagdiff, diagfactor * sep)
            int_t th = std::max(diagdiff, static_cast<int_t>(diagfactor * sep));

            // diagonal 差满足阈值：归并到同一簇
            if (th>=0 && diag_diff <= static_cast<uint_t>(th)) {
                uf.unite(i, j);
            }
        }
    }

    // 根据并查集根构建簇：root -> cluster_id
    std::unordered_map<int_t, int_t> root_to_cluster_id;
    root_to_cluster_id.reserve(N / 4);

    for (uint_t idx = 0; idx < unique_match.size(); idx++) {
        int_t root = uf.find(idx);
        auto it = root_to_cluster_id.find(root);

        int_t cid;
        if (it == root_to_cluster_id.end()) {
            cid = static_cast<int_t>(clusters.size());
            clusters.emplace_back();
            clusters.back().reserve(1);
            root_to_cluster_id[root] = cid;
        }
        else {
            cid = it->second;
        }

        clusters[cid].push_back(std::move(unique_match[idx]));
    }

    return clusters;
}

// ────────────────────────────────────────────
MatchVec bestChainDP(MatchVec& cluster, double diagfactor) {
    if (cluster.empty()) return {};
    if (cluster.size() == 1) {return MatchVec{cluster.front()};}

    Strand strand = cluster.front().strand;

    std::sort(cluster.begin(), cluster.end(),
        [](const Match& a, const Match& b) { return start2(a) < start2(b); });

    const uint_t N = static_cast<uint_t>(cluster.size());
    std::vector<int_t> score(N), pred(N, -1);
    uint_t best_idx = 0;

    for (uint_t i = 0; i < N; ++i) {
        // 初始分：自身长度
        score[i] = len2(cluster[i]);

        for (uint_t j = 0; j < i; ++j) {
            // query 维度必须不重叠（i 的 start2 需在 j 的末尾之后）
            if (start2(cluster[i]) <= start2(cluster[j])+len2(cluster[j])) continue;

            int_t d = 0;

            if (strand == FORWARD) {
                // ref 维度也不重叠：i 的 start1 需在 j 的 ref 末尾之后
                int_t prev_endj = start1(cluster[j]) + len1(cluster[j]);
                if (static_cast<int_t>(start1(cluster[i])) <= prev_endj) continue;

                d = std::abs(diag(cluster[i]) - diag(cluster[j]));
            }
            else {
                // 反向链：按原逻辑计算不交叉条件与 sep
                int_t prev_endi = start1(cluster[i]) + len1(cluster[i]);
                if (prev_endi >= static_cast<int_t>(start1(cluster[j]))) continue;

                d = std::abs(diag_reverse(cluster[i]) - diag_reverse(cluster[j]));
            }



            // 候选分数：前链分数 + 当前长度 - diagonal 差惩罚
            int_t cand = score[j] + len2(cluster[i]) - d;
            if (cand > score[i]) {
                score[i] = cand;
                pred[i] = static_cast<int_t>(j);
            }
        }

        if (score[i] > score[best_idx]) best_idx = i;
    }

    // 回溯构建最优链
    MatchVec chain;
    for (int_t k = static_cast<int_t>(best_idx); k != -1; k = pred[k]) {
        chain.emplace_back(cluster[k]);
    }
    std::reverse(chain.begin(), chain.end());

    return chain;
}
MatchClusterVecPtr clusterChrMatch(std::span<Match> unique_match,
                                   uint_t min_cluster_length,
                                   int_t max_gap,
                                   int_t diagdiff,
                                   double diagfactor) {

    auto best_chain_clusters = std::make_shared<MatchClusterVec>();

    if (unique_match.size() == 0) {
        return best_chain_clusters;
    }


    // 1) 聚簇
    MatchClusterVec clusters = buildClusters(unique_match, max_gap, diagdiff, diagfactor);

    // 释放 unique_match 的内存（保持原逻辑）
    // The owning Seed vector is released once all disjoint groups finish.
    unique_match = {};

    best_chain_clusters->reserve(clusters.size());

    // 2) 每个簇选最佳链
    for (auto& cluster : clusters) {
        if (cluster.empty()) continue;

        MatchVec best_chain = bestChainDP(cluster, diagfactor);


        if (best_chain.empty()) {
            releaseCluster(cluster);
            continue;
        }

        uint_t span = 0;
        // 遍历 best_chain 累加 span（保持原逻辑）
        for (auto& m : best_chain) {
            span += m.length;
        }

        // 满足最小簇长度才保留
        if (span >= min_cluster_length) {
            best_chain_clusters->emplace_back(std::move(best_chain));
        }

        // 回收 cluster 剩余元素
        releaseCluster(cluster);
    }

    best_chain_clusters->shrink_to_fit();
    return best_chain_clusters;
}
namespace {
constexpr Coord_t kMaximumLinkedGapAlignmentLength = 10000;

constexpr bool linkedGapCanBeAligned(Coord_t ref_gap, Coord_t query_gap) {
    return ref_gap <= kMaximumLinkedGapAlignmentLength &&
           query_gap <= kMaximumLinkedGapAlignmentLength;
}

static_assert(linkedGapCanBeAligned(9999, 9999));
static_assert(linkedGapCanBeAligned(10000, 10000));
static_assert(!linkedGapCanBeAligned(10001, 10000));
static_assert(!linkedGapCanBeAligned(10000, 10001));

struct ManagedSequenceSlice {
    std::string storage;
    KswSequenceView view;
};

const SequenceViews::SequenceManager& originalSequenceManager(
    const SequenceViews::ManagerVariant& manager) {
    return std::visit([](const auto& pointer) -> const SequenceViews::SequenceManager& {
        using Pointer = std::decay_t<decltype(pointer)>;
        if constexpr (std::is_same_v<
                          Pointer,
                          std::unique_ptr<SequenceViews::SequenceManager>>) {
            return *pointer;
        } else {
            return pointer->getOriginalManager();
        }
    }, manager);
}

void loadSequenceSlice(const SequenceViews::ManagerVariant& manager,
                       ChrIndex chromosome, Coord_t start, Coord_t length,
                       bool reverse_complement,
                       ManagedSequenceSlice& output) {
    output.storage.clear();
    if (length == 0) {
        output.view = {{}, reverse_complement};
        return;
    }
    const auto& original = originalSequenceManager(manager);
    std::string_view contiguous;
    if (original.tryGetContiguousSubSequence(
            chromosome, start, length, contiguous)) {
        output.view = {contiguous, reverse_complement};
        return;
    }
    original.getSubSequenceInto(
        chromosome, start, length, output.storage);
    output.view = {output.storage, reverse_complement};
}

void appendCigarMove(Cigar_t& destination, Cigar_t& source) {
    if (source.empty()) return;
    if (destination.empty()) {
        destination = std::move(source);
        Cigar_t().swap(source);
        return;
    }
    destination.reserve(destination.size() + source.size());
    appendCigar(destination, source);
    Cigar_t().swap(source);
}
}  // namespace

Anchor extendClusterToAnchor(MatchCluster& cluster,
    const SequenceViews::ManagerVariant& ref_mgr,
    const SequenceViews::ManagerVariant& query_mgr) {
    if (cluster.empty()) return Anchor();
    Anchor anchor;
    const Match& first = cluster.front();

    Strand strand = first.strand;
    bool   fwd = (strand == FORWARD);

    if (!fwd) {
        std::sort(cluster.begin(), cluster.end(),
            [](auto& a, auto& b) { return a.reference_begin < b.reference_begin; });
    }

    ChrIndex ref_chr = first.reference_id;
    ChrIndex qry_chr = first.query_id;

    Cigar_t cigar; cigar.reserve(cluster.size() * 2);  // 预估
    Coord_t aln_len = 0;
    Coord_t match_len = 0;
    ManagedSequenceSlice reference_slice;
    ManagedSequenceSlice query_slice;
    for (size_t i = 0;i < cluster.size();++i) {
        const Match& m = cluster[i];

		uint_t len = len1(m);
        uint_t reference_begin = m.reference_begin + len;
        uint_t query_begin = m.query_begin + len;

        appendCigarOp(cigar, 'M', len);
        aln_len += len;
        match_len += len;

        if (i + 1 == cluster.size()) break;

        const Match& nxt = cluster[i + 1];
        uint_t len2 = len1(nxt);
		uint_t ref_end = nxt.reference_begin;
		uint_t qry_end = nxt.query_begin;
        Coord_t query_gap_begin = 0;
        Coord_t query_gap_length = 0;
        if (fwd) {
            query_gap_begin = query_begin;
            query_gap_length = qry_end - query_begin;
        }
        else {
            query_gap_begin = qry_end + len2;
            query_gap_length = m.query_begin - len2 - qry_end;
        }

        loadSequenceSlice(ref_mgr, ref_chr, reference_begin,
            ref_end - reference_begin, false, reference_slice);
        loadSequenceSlice(query_mgr, qry_chr, query_gap_begin,
            query_gap_length, !fwd, query_slice);
        AlignmentResult gap = globalAlignKSW2Result(
            reference_slice.view, query_slice.view);
        match_len += gap.summary.match_length;
        aln_len += gap.summary.alignment_length;
        appendCigarMove(cigar, gap.cigar);

    }
	const Match& last = cluster.back();
    if (fwd) {
        anchor = Anchor(ref_chr, first.reference_begin, last.reference_begin + last.length - first.reference_begin, qry_chr, first.query_begin, last.query_begin + last.length - first.query_begin, strand, aln_len, match_len, std::move(cigar));
    }
    else {
        anchor = Anchor(ref_chr, first.reference_begin, last.reference_begin + last.length - first.reference_begin, qry_chr, last.query_begin, first.query_begin + first.length - last.query_begin, strand, aln_len, match_len, std::move(cigar));
    }

    return anchor;
}
inline bool intervalOverlap(Coord_t a_lo, Coord_t a_hi, Coord_t b_lo, Coord_t b_hi) {
    // 闭开区间 [lo,hi)，不重叠当且仅当 a_hi <= b_lo 或 b_hi <= a_lo
    return !(a_hi <= b_lo || b_hi <= a_lo);
}
namespace AnchorLinkDetail {

struct Statistics {
    uint64_t candidate_checks{0};
    uint64_t sequence_extractions{0};
    uint64_t direct_ksw_calls{0};
    uint64_t fallback_ksw_calls{0};
    uint64_t long_gap_rejections{0};
    uint64_t maximum_seen_gap{0};
    uint64_t estimated_ksw_cells{0};
};

struct ComponentRange {
    size_t begin{0};
    size_t end{0};
    uint64_t estimated_cost{0};
};

AnchorVec materializeClusterAnchors(
    MatchClusterVec& clusters,
    const SequenceViews::ManagerVariant& ref_mgr,
    const SequenceViews::ManagerVariant& qry_mgr);

std::vector<ComponentRange> splitAnchorComponents(
    const AnchorVec& anchors,
    Statistics* statistics = nullptr);

AnchorVec linkAnchorRange(
    AnchorVec& anchors,
    size_t begin,
    size_t end,
    const SequenceViews::ManagerVariant& ref_mgr,
    const SequenceViews::ManagerVariant& qry_mgr,
    Statistics* statistics = nullptr);

}  // namespace AnchorLinkDetail

namespace AnchorLinkDetail {

namespace {

uint64_t saturatingAdd(uint64_t left, uint64_t right) {
    return right > std::numeric_limits<uint64_t>::max() - left
        ? std::numeric_limits<uint64_t>::max()
        : left + right;
}

uint64_t saturatingMultiply(uint64_t left, uint64_t right) {
    if (left != 0 &&
        right > std::numeric_limits<uint64_t>::max() / left) {
        return std::numeric_limits<uint64_t>::max();
    }
    return left * right;
}

uint64_t estimatedGapCells(Coord_t ref_gap, Coord_t query_gap) {
    const uint64_t maximum = std::max<uint64_t>(ref_gap, query_gap);
    const uint64_t band = static_cast<uint64_t>(auto_band(
        static_cast<int>(ref_gap), static_cast<int>(query_gap)));
    return saturatingMultiply(maximum, saturatingAdd(band, band) + 1);
}

void observeGap(Statistics* statistics, int_t ref_gap, int_t query_gap,
                bool rejected) {
    if (!statistics || ref_gap < 0 || query_gap < 0) return;
    statistics->maximum_seen_gap = std::max<uint64_t>(
        statistics->maximum_seen_gap,
        static_cast<uint64_t>(std::max(ref_gap, query_gap)));
    if (rejected) ++statistics->long_gap_rejections;
}

}  // namespace

AnchorVec materializeClusterAnchors(
    MatchClusterVec& clusters,
    const SequenceViews::ManagerVariant& ref_mgr,
    const SequenceViews::ManagerVariant& qry_mgr) {
    std::sort(clusters.begin(), clusters.end(),
        [](const MatchCluster& left, const MatchCluster& right) {
            if (left.empty() || right.empty()) {
                return left.size() < right.size();
            }
            return left.front().reference_begin < right.front().reference_begin;
        });

    MatchClusterVec cleaned;
    cleaned.reserve(clusters.size());
    bool have_previous = false;
    ChrIndex previous_ref_chromosome = 0;
    ChrIndex previous_query_chromosome = 0;
    Strand previous_strand = FORWARD;
    Coord_t previous_ref_end = 0;
    Coord_t previous_query_low = 0;
    Coord_t previous_query_high = 0;

    const auto query_bounds = [](const MatchCluster& cluster) {
        Coord_t low = std::numeric_limits<Coord_t>::max();
        Coord_t high = 0;
        for (const auto& match : cluster) {
            const Coord_t first = match.query_begin;
            const Coord_t second = match.query_begin + len2(match);
            low = std::min(low, std::min(first, second));
            high = std::max(high, std::max(first, second));
        }
        if (low == std::numeric_limits<Coord_t>::max()) low = 0;
        return std::pair<Coord_t, Coord_t>{low, high};
    };

    for (auto& cluster : clusters) {
        if (cluster.empty()) continue;
        if (!have_previous ||
            cluster.front().reference_id != previous_ref_chromosome ||
            cluster.front().query_id != previous_query_chromosome ||
            cluster.front().strand != previous_strand) {
            MatchCluster kept = std::move(cluster);
            previous_ref_chromosome = kept.front().reference_id;
            previous_query_chromosome = kept.front().query_id;
            previous_strand = kept.front().strand;
            previous_ref_end = kept.back().reference_begin + len1(kept.back());
            const auto [low, high] = query_bounds(kept);
            previous_query_low = low;
            previous_query_high = high;
            cleaned.push_back(std::move(kept));
            have_previous = true;
            continue;
        }

        MatchCluster pruned;
        pruned.reserve(cluster.size());
        for (auto& match : cluster) {
            const bool reference_ok = match.reference_begin >= previous_ref_end;
            const Coord_t first = match.query_begin;
            const Coord_t second = match.query_begin + len2(match);
            const Coord_t low = std::min(first, second);
            const Coord_t high = std::max(first, second);
            const bool query_ok = !intervalOverlap(
                low, high, previous_query_low, previous_query_high);
            if (reference_ok && query_ok) {
                pruned.push_back(std::move(match));
            }
        }

        if (pruned.empty()) continue;
        previous_ref_end = pruned.back().reference_begin + len1(pruned.back());
        const auto [low, high] = query_bounds(pruned);
        previous_query_low = low;
        previous_query_high = high;
        cleaned.push_back(std::move(pruned));
    }
    clusters.swap(cleaned);
    MatchClusterVec().swap(cleaned);

    AnchorVec anchors;
    anchors.reserve(clusters.size());
    for (auto& cluster : clusters) {
        if (cluster.empty()) continue;
        anchors.push_back(extendClusterToAnchor(cluster, ref_mgr, qry_mgr));
        anchors.back().is_linked = false;
        MatchCluster().swap(cluster);
    }
    MatchClusterVec().swap(clusters);
    return anchors;
}

std::vector<ComponentRange> splitAnchorComponents(
    const AnchorVec& anchors,
    Statistics* statistics) {
    std::vector<ComponentRange> components;
    if (anchors.empty()) return components;
    size_t begin = 0;
    uint64_t prefix_maximum_end =
        static_cast<uint64_t>(anchors.front().reference_begin) +
        anchors.front().ref_len;
    const auto append_component = [&](size_t first, size_t last) {
        const uint64_t length = last - first;
        uint64_t cost = saturatingAdd(
            length,
            saturatingMultiply(length, std::min<uint64_t>(length, 2000)));
        for (size_t index = first + 1; index < last; ++index) {
            const Anchor& previous = anchors[index - 1];
            const Anchor& current = anchors[index];
            if (previous.reference_id != current.reference_id ||
                previous.query_id != current.query_id ||
                previous.strand != current.strand) {
                continue;
            }
            const int64_t ref_gap = static_cast<int64_t>(current.reference_begin) -
                static_cast<int64_t>(previous.reference_begin + previous.ref_len);
            const int64_t query_gap = current.strand == FORWARD
                ? static_cast<int64_t>(current.query_begin) -
                    static_cast<int64_t>(
                        previous.query_begin + previous.qry_len)
                : static_cast<int64_t>(previous.query_begin) -
                    static_cast<int64_t>(current.query_begin + current.qry_len);
            if (ref_gap < 0 || query_gap < 0 ||
                !linkedGapCanBeAligned(
                    static_cast<Coord_t>(ref_gap),
                    static_cast<Coord_t>(query_gap))) {
                continue;
            }
            cost = saturatingAdd(cost, estimatedGapCells(
                static_cast<Coord_t>(ref_gap),
                static_cast<Coord_t>(query_gap)));
        }
        components.push_back({first, last, cost});
    };

    for (size_t index = 1; index < anchors.size(); ++index) {
        const Anchor& previous = anchors[index - 1];
        const Anchor& current = anchors[index];
        const bool key_changed =
            previous.reference_id != current.reference_id ||
            previous.query_id != current.query_id ||
            previous.strand != current.strand;
        const uint64_t current_start = current.reference_begin;
        const bool separated = current_start > prefix_maximum_end &&
            current_start - prefix_maximum_end >
                kMaximumLinkedGapAlignmentLength;
        if (key_changed || separated) {
            if (!key_changed && separated && statistics) {
                ++statistics->long_gap_rejections;
                statistics->maximum_seen_gap = std::max<uint64_t>(
                    statistics->maximum_seen_gap,
                    current_start - prefix_maximum_end);
            }
            append_component(begin, index);
            begin = index;
            prefix_maximum_end =
                static_cast<uint64_t>(current.reference_begin) + current.ref_len;
        } else {
            prefix_maximum_end = std::max<uint64_t>(
                prefix_maximum_end,
                static_cast<uint64_t>(current.reference_begin) + current.ref_len);
        }
    }
    append_component(begin, anchors.size());
    return components;
}

AnchorVec linkAnchorRange(
    AnchorVec& anchors,
    size_t begin,
    size_t end,
    const SequenceViews::ManagerVariant& ref_mgr,
    const SequenceViews::ManagerVariant& qry_mgr,
    Statistics* statistics) {
    if (begin > end || end > anchors.size()) {
        throw std::out_of_range("Invalid Anchor linking component range");
    }
    AnchorVec output;
    if (begin == end) return output;
    if(active_statistics)active_statistics->link_input_candidates+=end-begin;
    std::vector<size_t> linked;
    linked.reserve(end - begin);
    constexpr size_t kNoIndex = std::numeric_limits<size_t>::max();
    constexpr int kCandidateLimit = 2000;
    constexpr int kLookBack = 2000;
    constexpr int_t kBreakLength = 200;
    ManagedSequenceSlice reference_slice;
    ManagedSequenceSlice query_slice;

    size_t current_index = begin;
    while (current_index < end) {
        Anchor& current = anchors[current_index];
        if (current.is_linked) {
            ++current_index;
            continue;
        }

        size_t best_index = kNoIndex;
        int_t best_score = std::numeric_limits<int_t>::max();
        int looked = 0;
        for (size_t index = current_index + 1;
             index < end && looked < kCandidateLimit; ++index) {
            Anchor& candidate = anchors[index];
            if (candidate.is_linked) continue;
            ++looked;
            if (statistics) ++statistics->candidate_checks;

            const int_t ref_gap = static_cast<int_t>(candidate.reference_begin) -
                static_cast<int_t>(current.reference_begin + current.ref_len);
            const int_t query_gap = current.strand == FORWARD
                ? static_cast<int_t>(candidate.query_begin) -
                    static_cast<int_t>(current.query_begin + current.qry_len)
                : static_cast<int_t>(current.query_begin) -
                    static_cast<int_t>(candidate.query_begin + candidate.qry_len);
            if (ref_gap < 0 || query_gap < 0) continue;

            const long greater = std::max(ref_gap, query_gap);
            const long lesser = std::min(ref_gap, query_gap);
            if (greater < kBreakLength ||
                greater - lesser <= kBreakLength) {
                best_index = index;
                break;
            }
            const int_t score = (greater << 1) - lesser;
            if (best_score > score) {
                best_score = score;
                best_index = index;
            }
        }

        bool reached = false;
        if (best_index != kNoIndex) {
            Anchor& best = anchors[best_index];
            const Coord_t ref_gap_begin = current.reference_begin + current.ref_len;
            const Coord_t ref_gap_length = best.reference_begin - ref_gap_begin;
            Coord_t query_gap_begin = 0;
            Coord_t query_gap_length = 0;
            if (current.strand == FORWARD) {
                query_gap_begin = current.query_begin + current.qry_len;
                query_gap_length = best.query_begin - query_gap_begin;
            } else {
                query_gap_begin = best.query_begin + best.qry_len;
                query_gap_length = current.query_begin - query_gap_begin;
            }

            const bool alignable = linkedGapCanBeAligned(
                ref_gap_length, query_gap_length);
            observeGap(statistics, ref_gap_length, query_gap_length,
                       !alignable);
            if (alignable) {
                if (statistics) {
                    statistics->sequence_extractions += 2;
                    ++statistics->direct_ksw_calls;
                    statistics->estimated_ksw_cells = saturatingAdd(
                        statistics->estimated_ksw_cells,
                        estimatedGapCells(ref_gap_length, query_gap_length));
                }
                loadSequenceSlice(ref_mgr, current.reference_id,
                    ref_gap_begin, ref_gap_length, false, reference_slice);
                loadSequenceSlice(qry_mgr, current.query_id,
                    query_gap_begin, query_gap_length,
                    current.strand == REVERSE, query_slice);
                AlignmentResult gap = extendAlignKSW2Result(
                    reference_slice.view, query_slice.view,
                    2 * kBreakLength);
                if (gap.summary.reference_length == ref_gap_length &&
                    gap.summary.query_length == query_gap_length) {
                    reached = true;
                    current.ref_len = best.reference_begin + best.ref_len -
                        current.reference_begin;
                    if (current.strand == FORWARD) {
                        current.qry_len = best.query_begin + best.qry_len -
                            current.query_begin;
                    } else {
                        current.qry_len = current.query_begin + current.qry_len -
                            best.query_begin;
                        current.query_begin = best.query_begin;
                    }
                    current.cigar.reserve(current.cigar.size() +
                        gap.cigar.size() + best.cigar.size());
                    appendCigarMove(current.cigar, gap.cigar);
                    appendCigarMove(current.cigar, best.cigar);
                    current.alignment_length += best.alignment_length +
                        gap.summary.alignment_length;
                    current.aligned_base += gap.summary.match_length +
                        best.aligned_base;
                    best.is_linked = true;
                    if(active_statistics)++active_statistics->link_consumed_candidates;
                } else if(active_statistics) {
                    ++active_statistics->link_closure_failures;
                }
            }
        }

        if (!reached) {
            if (!linked.empty()) {
                size_t best_linked_position = kNoIndex;
                int_t linked_best_score = std::numeric_limits<int_t>::max();
                int checked = 0;
                for (size_t position = linked.size();
                     position > 0 && checked < kLookBack;
                     --position, ++checked) {
                    const Anchor& previous = anchors[linked[position - 1]];
                    if (statistics) ++statistics->candidate_checks;
                    if (previous.strand != current.strand ||
                        previous.reference_id != current.reference_id ||
                        previous.query_id != current.query_id) {
                        continue;
                    }
                    const int_t ref_gap = static_cast<int_t>(current.reference_begin) -
                        static_cast<int_t>(
                            previous.reference_begin + previous.ref_len);
                    const int_t query_gap = current.strand == FORWARD
                        ? static_cast<int_t>(current.query_begin) -
                            static_cast<int_t>(
                                previous.query_begin + previous.qry_len)
                        : static_cast<int_t>(previous.query_begin) -
                            static_cast<int_t>(
                                current.query_begin + current.qry_len);
                    if (ref_gap < 0 || query_gap < 0) continue;
                    const bool alignable = linkedGapCanBeAligned(
                        static_cast<Coord_t>(ref_gap),
                        static_cast<Coord_t>(query_gap));
                    observeGap(statistics, ref_gap, query_gap, !alignable);
                    if (!alignable) continue;
                    const long greater = std::max(ref_gap, query_gap);
                    const long lesser = std::min(ref_gap, query_gap);
                    const int_t score = (greater << 1) - lesser;
                    if (score < linked_best_score) {
                        linked_best_score = score;
                        best_linked_position = position - 1;
                    }
                }

                if (best_linked_position != kNoIndex) {
                    const Anchor& previous =
                        anchors[linked[best_linked_position]];
                    const Coord_t ref_gap_begin =
                        previous.reference_begin + previous.ref_len;
                    const Coord_t ref_gap_length =
                        current.reference_begin - ref_gap_begin;
                    Coord_t query_gap_begin = 0;
                    Coord_t query_gap_length = 0;
                    if (current.strand == FORWARD) {
                        query_gap_begin = previous.query_begin + previous.qry_len;
                        query_gap_length = current.query_begin - query_gap_begin;
                    } else {
                        query_gap_begin = current.query_begin + current.qry_len;
                        query_gap_length = previous.query_begin - query_gap_begin;
                    }
                    if (statistics) {
                        statistics->sequence_extractions += 2;
                        ++statistics->fallback_ksw_calls;
                        statistics->estimated_ksw_cells = saturatingAdd(
                            statistics->estimated_ksw_cells,
                            estimatedGapCells(
                                ref_gap_length, query_gap_length));
                    }
                    loadSequenceSlice(ref_mgr, current.reference_id,
                        ref_gap_begin, ref_gap_length, false, reference_slice);
                    loadSequenceSlice(qry_mgr, current.query_id,
                        query_gap_begin, query_gap_length,
                        current.strand == REVERSE, query_slice);
                    AlignmentResult gap = extendAlignKSW2Result(
                        reference_slice.view, query_slice.view,
                        2 * kBreakLength);
                    if (gap.summary.reference_length == ref_gap_length &&
                        gap.summary.query_length == query_gap_length) {
                        current.ref_len += gap.summary.reference_length;
                        current.qry_len += gap.summary.query_length;
                        current.alignment_length +=
                            gap.summary.alignment_length;
                        current.aligned_base += gap.summary.match_length;
                        current.reference_begin -= gap.summary.reference_length;
                        if (current.strand == FORWARD) {
                            current.query_begin -= gap.summary.query_length;
                        }
                        prependCigar(current.cigar, gap.cigar);
                    } else if(active_statistics) {
                        ++active_statistics->link_closure_failures;
                    }
                }
            }

            linked.push_back(current_index);
            const Coord_t current_ref_end = current.reference_begin + current.ref_len;
            const Coord_t current_query_end = current.strand == FORWARD
                ? current.query_begin + current.qry_len
                : current.query_begin;
            const size_t stop = best_index == kNoIndex ? end : best_index;
            for (size_t index = current_index + 1; index < stop; ++index) {
                Anchor& candidate = anchors[index];
                if (candidate.is_linked) continue;
                const Coord_t candidate_ref_end =
                    candidate.reference_begin + candidate.ref_len;
                const Coord_t candidate_query_end = current.strand == FORWARD
                    ? candidate.query_begin + candidate.qry_len
                    : candidate.query_begin;
                if (candidate_ref_end <= current_ref_end &&
                    ((current.strand == FORWARD &&
                      candidate_query_end <= current_query_end) ||
                     (current.strand == REVERSE &&
                      candidate_query_end >= current_query_end))) {
                    candidate.is_linked = true;
                    Cigar_t().swap(candidate.cigar);
                }
            }
            while (current_index < end) {
                const Anchor& candidate = anchors[current_index];
                if (!candidate.is_linked &&
                    candidate.reference_begin >= current_ref_end) {
                    break;
                }
                ++current_index;
            }
        }
    }


    if(active_statistics)active_statistics->link_emitted_candidates+=linked.size();
    output.reserve(linked.size());
    for (const size_t index : linked) {
        anchors[index].is_linked = false;
        output.push_back(std::move(anchors[index]));
    }
    return output;
}

}  // namespace AnchorLinkDetail


constexpr size_t kMaximumRetainedDpAnchors = 131072;
std::atomic<uint64_t> dp_treap_fallback_count{0};

struct DpTreapBest {
    double value = -std::numeric_limits<double>::infinity();
    size_t index = std::numeric_limits<size_t>::max();

    bool valid() const {
        return index != std::numeric_limits<size_t>::max();
    }
};

DpTreapBest betterDpBest(DpTreapBest left, DpTreapBest right) {
    if (!left.valid()) return right;
    if (!right.valid()) return left;
    if (right.value > left.value ||
        (right.value == left.value && right.index < left.index)) {
        return right;
    }
    return left;
}

uint64_t splitmix64(uint64_t value) {
    value += 0x9e3779b97f4a7c15ULL;
    value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31);
}

struct DpTreapNode {
    long long interval_end = 0;
    size_t original_index = 0;
    uint64_t priority = 0;
    double dp = 0;
    int left = -1;
    int right = -1;
    DpTreapBest subtree_best;
};

class FixedWindowDpTreap {
public:
    explicit FixedWindowDpTreap(std::vector<DpTreapNode>& storage)
        : nodes_(storage) {}

    void reset(size_t node_count) {
        nodes_.resize(node_count);
        root_ = -1;
    }

    void insert(size_t slot, long long interval_end,
                size_t original_index, double dp) {
        auto& node = nodes_[slot];
        node.interval_end = interval_end;
        node.original_index = original_index;
        node.priority = splitmix64(original_index);
        node.dp = dp;
        node.left = -1;
        node.right = -1;
        node.subtree_best = std::isnan(dp)
            ? DpTreapBest{}
            : DpTreapBest{dp, original_index};
        root_ = insertNode(root_, static_cast<int>(slot));
    }

    void erase(long long interval_end, size_t original_index) {
        root_ = eraseNode(root_, interval_end, original_index);
    }

    DpTreapBest bestAll() const {
        return bestForNode(root_);
    }

    DpTreapBest bestEndingAtOrBefore(long long coordinate) const {
        DpTreapBest best;
        int current = root_;
        while (current >= 0) {
            const auto& node = nodes_[current];
            if (node.interval_end <= coordinate) {
                best = betterDpBest(best, bestForNode(node.left));
                if (!std::isnan(node.dp)) {
                    best = betterDpBest(
                        best, {node.dp, node.original_index});
                }
                current = node.right;
            } else {
                current = node.left;
            }
        }
        return best;
    }

private:
    bool keyLess(int left, int right) const {
        const auto& lhs = nodes_[left];
        const auto& rhs = nodes_[right];
        return lhs.interval_end < rhs.interval_end ||
            (lhs.interval_end == rhs.interval_end &&
             lhs.original_index < rhs.original_index);
    }

    DpTreapBest bestForNode(int node) const {
        return node < 0 ? DpTreapBest{} : nodes_[node].subtree_best;
    }

    void update(int node) {
        if (node < 0) return;
        auto best = std::isnan(nodes_[node].dp)
            ? DpTreapBest{}
            : DpTreapBest{nodes_[node].dp,
                          nodes_[node].original_index};
        best = betterDpBest(best, bestForNode(nodes_[node].left));
        best = betterDpBest(best, bestForNode(nodes_[node].right));
        nodes_[node].subtree_best = best;
    }

    void split(int node, int key, int& left, int& right) {
        if (node < 0) {
            left = right = -1;
            return;
        }
        if (keyLess(node, key)) {
            left = node;
            split(nodes_[node].right, key,
                  nodes_[node].right, right);
            update(left);
        } else {
            right = node;
            split(nodes_[node].left, key,
                  left, nodes_[node].left);
            update(right);
        }
    }

    int merge(int left, int right) {
        if (left < 0) return right;
        if (right < 0) return left;
        if (nodes_[left].priority > nodes_[right].priority) {
            nodes_[left].right = merge(nodes_[left].right, right);
            update(left);
            return left;
        }
        nodes_[right].left = merge(left, nodes_[right].left);
        update(right);
        return right;
    }

    int insertNode(int root, int node) {
        if (root < 0) return node;
        if (nodes_[node].priority > nodes_[root].priority) {
            split(root, node, nodes_[node].left, nodes_[node].right);
            update(node);
            return node;
        }
        if (keyLess(node, root)) {
            nodes_[root].left = insertNode(nodes_[root].left, node);
        } else {
            nodes_[root].right = insertNode(nodes_[root].right, node);
        }
        update(root);
        return root;
    }

    int eraseNode(int root, long long interval_end,
                  size_t original_index) {
        if (root < 0) return -1;
        const auto& node = nodes_[root];
        if (node.interval_end == interval_end &&
            node.original_index == original_index) {
            return merge(node.left, node.right);
        }
        if (interval_end < node.interval_end ||
            (interval_end == node.interval_end &&
             original_index < node.original_index)) {
            nodes_[root].left = eraseNode(
                nodes_[root].left, interval_end, original_index);
        } else {
            nodes_[root].right = eraseNode(
                nodes_[root].right, interval_end, original_index);
        }
        update(root);
        return root;
    }

    std::vector<DpTreapNode>& nodes_;
    int root_ = -1;
};

struct DpWorkspace {
    std::vector<double> dp;
    std::vector<int_t> pre;
    std::vector<DpTreapNode> treap_nodes;
};

thread_local DpWorkspace retained_dp_workspace;

void filterAnchorsByDpTreap(std::span<size_t> result, AnchorVec& anchors, bool filter_ref) {
    if (result.empty()) return;
    std::sort(result.begin(), result.end(),
        [filter_ref, &anchors](size_t left, size_t right) {
            return filter_ref ? anchors[left].reference_begin < anchors[right].reference_begin
                              : anchors[left].query_begin < anchors[right].query_begin;
        });

    DpWorkspace temporary_workspace;
    DpWorkspace& workspace = result.size() <= kMaximumRetainedDpAnchors
        ? retained_dp_workspace : temporary_workspace;
    workspace.dp.assign(result.size(), 0);
    workspace.pre.assign(result.size(), -1);
    FixedWindowDpTreap treap(workspace.treap_nodes);
    treap.reset(std::min(result.size(), kDpWindow + 1));

    const auto interval = [&](size_t index) {
        if (filter_ref) {
            return std::pair<long long, long long>{
                static_cast<long long>(anchors[result[index]].reference_begin),
                static_cast<long long>(anchors[result[index]].ref_len)};
        }
        return std::pair<long long, long long>{
            static_cast<long long>(anchors[result[index]].query_begin),
            static_cast<long long>(anchors[result[index]].qry_len)};
    };
    const auto legacyPredecessor = [&](size_t current, double score) {
        const size_t begin = current > kDpWindow
            ? current - kDpWindow : 0;
        const auto [current_start, current_length] = interval(current);
        const long long current_end = current_start + current_length;
        for (size_t previous = begin; previous < current; ++previous) {
            const auto [previous_start, previous_length] = interval(previous);
            const long long previous_end = previous_start + previous_length;
            const long long overlap = std::max(
                0LL, std::min(previous_end, current_end) -
                         std::max(previous_start, current_start));
            const long long shorter =
                std::min(previous_length, current_length);
            const double overlap_ratio = shorter > 0
                ? static_cast<double>(overlap) /
                    static_cast<double>(shorter)
                : 0.0;
            if (overlap_ratio <= 0.0) {
                const double candidate = workspace.dp[previous] +
                    score - static_cast<double>(overlap);
                if (candidate > workspace.dp[current]) {
                    workspace.dp[current] = candidate;
                    workspace.pre[current] = static_cast<int_t>(previous);
                }
            }
        }
    };

    uint64_t local_fallbacks = 0;
    for (size_t index = 0; index < result.size(); ++index) {
        const double identity = static_cast<float>(
            anchors[result[index]].aligned_base) /
            anchors[result[index]].alignment_length;
        const double score = anchors[result[index]].alignment_length *
            pow(identity, 2);
        workspace.dp[index] = score;

        if (index > kDpWindow) {
            const size_t expired = index - kDpWindow - 1;
            const auto [start, length] = interval(expired);
            treap.erase(start + length, expired);
        }

        const auto [current_start, current_length] = interval(index);
        const DpTreapBest best = current_length == 0
            ? treap.bestAll()
            : treap.bestEndingAtOrBefore(current_start);
        bool fallback = anchors[result[index]].alignment_length == 0 ||
            !std::isfinite(score) ||
            (best.valid() && !std::isfinite(best.value));
        if (best.valid() && !fallback) {
            const double candidate = best.value + score;
            const double previous_representable = std::nextafter(
                best.value, -std::numeric_limits<double>::infinity());
            fallback = previous_representable + score == candidate;
            if (!fallback && candidate > workspace.dp[index]) {
                workspace.dp[index] = candidate;
                workspace.pre[index] = static_cast<int_t>(best.index);
            }
        }
        if (fallback) {
            workspace.dp[index] = score;
            workspace.pre[index] = -1;
            legacyPredecessor(index, score);
            ++local_fallbacks;
        }

        treap.insert(index % (kDpWindow + 1),
            current_start + current_length, index,
            workspace.dp[index]);
    }
    dp_treap_fallback_count.fetch_add(
        local_fallbacks, std::memory_order_relaxed);

    double best = 0.0;
    size_t best_index = 0;
    for (size_t index = 0; index < result.size(); ++index) {
        if (workspace.dp[index] > best) {
            best = workspace.dp[index];
            best_index = index;
        }
    }
    for (int index = static_cast<int>(best_index); index >= 0;
         index = workspace.pre[index]) {
        if (filter_ref) anchors[result[index]].ref_selected = true;
        else anchors[result[index]].qry_selected = true;
        if (workspace.pre[index] == -1) break;
    }
}


} // namespace ramag::pairwise_detail
// RaMA-G ownership, validation and serialization bridge. Algorithm bodies above
// retain pairwise control flow; this layer neither reconnects nor refilters anchors.
#include <chrono>
#include <type_traits>
namespace ramag {
using namespace genome_detail;

static PairwiseCoreResult AlignPairwiseCoreImpl(std::span<const SequenceRecord> references,std::span<const SequenceRecord> queries,std::vector<Seed> seeds,const PairwiseCoreOptions& o,bool selected_only,std::uint64_t& candidate_count) {
    if(o.threads==0 || o.threads>static_cast<unsigned>(INT32_MAX) || o.min_cluster==0 || !std::isfinite(o.diag_factor) || o.diag_factor<0 || o.max_gap>static_cast<Length>(INT64_MAX/256) || o.diag_diff>static_cast<Length>(INT64_MAX/256)) throw AlignmentError("invalid pairwise core configuration");
    auto refs=CheckSequences(references), qrys=CheckSequences(queries);
    for(const auto& q:queries)if(static_cast<long double>(o.diag_factor)*q.size()>static_cast<long double>(INT64_MAX/2))throw AlignmentError("pairwise diagonal threshold cannot be represented");
    const auto memory_begin=KernelClock::now();
    PairwiseCoreResult result;
    const auto observe=[&](std::string stage, std::uint64_t cluster_bytes=0,
            std::uint64_t anchor_bytes=0, std::uint64_t cigar_bytes=0,
            std::uint64_t auxiliary_bytes=0) {
        MemoryObservation sample;
        sample.stage=std::move(stage);sample.elapsed_seconds=Elapsed(memory_begin);
        sample.rss_bytes=o.resident_bytes_callback?o.resident_bytes_callback():0;
        sample.seed_capacity_bytes=seeds.capacity()*sizeof(Seed);
        sample.cluster_capacity_bytes=cluster_bytes;sample.anchor_capacity_bytes=anchor_bytes;
        sample.cigar_capacity_bytes=cigar_bytes;sample.auxiliary_capacity_bytes=auxiliary_bytes;
        result.memory_observations.push_back(std::move(sample));
    };
    observe("core-seeds-owned");
    using Key=std::tuple<SequenceId,SequenceId,pairwise_detail::Strand>;

    KernelProgress(o,"seed-merge",0,seeds.size());
    for(const auto& s:seeds) {
        KernelCheck(o,"seed-merge");
        if(!refs.contains(s.reference_id)||!qrys.contains(s.query_id))throw AlignmentError("pairwise seed has unknown sequence ID");
        const auto& r=refs.at(s.reference_id)->bases;const auto& q=qrys.at(s.query_id)->bases;
        if(s.length==0||s.reference_begin>r.size()||s.length>r.size()-s.reference_begin||s.query_begin>q.size()||s.length>q.size()-s.query_begin)throw AlignmentError("pairwise seed outside sequence");
        if(s.strand!=Strand::Forward&&s.strand!=Strand::Reverse)throw AlignmentError("pairwise invalid seed strand");
        for(Length i=0;i<s.length;++i) {
            const auto x=r[s.reference_begin+i];const auto y=s.strand==Strand::Forward?q[s.query_begin+i]:Complement(q[s.query_begin+s.length-i-1]);
            if(x=='N'||x!=y)throw AlignmentError("pairwise seed is not an exact canonical match");
        }

    }
    const auto group_key=[](const Seed& seed) {
        return Key{seed.reference_id,seed.query_id,seed.strand};
    };
    // Within a group this is the former bucket ordering, including length.
    // Equivalent elements have identical Seed fields and are interchangeable.
    std::sort(seeds.begin(),seeds.end(),[&](const Seed& a,const Seed& b) {
        return std::tuple{group_key(a),a.query_begin,a.reference_begin,a.length}
             < std::tuple{group_key(b),b.query_begin,b.reference_begin,b.length};
    });
    std::vector<std::span<Seed>> groups;
    for(size_t begin=0;begin<seeds.size();) {
        size_t end=begin+1;
        while(end<seeds.size()&&group_key(seeds[end])==group_key(seeds[begin]))++end;
        groups.emplace_back(seeds.data()+begin,end-begin);
        begin=end;
    }
    result.statistics.seed_grouping_seconds=Elapsed(memory_begin);
    observe("core-seeds-grouped",0,0,0,groups.capacity()*sizeof(std::span<Seed>));
    result.workers=std::min(o.threads,static_cast<std::uint32_t>(std::min<std::size_t>(groups.size(),UINT32_MAX)));
    KernelProgress(o,"chaining",0,groups.size());
    auto start=KernelClock::now();
    std::vector<pairwise_detail::MatchClusterVecPtr> clusters(groups.size());
    ParallelKernel(groups.size(),o,[&](std::size_t i){
        KernelCheck(o,"chaining");
        clusters[i]=pairwise_detail::clusterChrMatch(groups[i],o.min_cluster,static_cast<std::int64_t>(o.max_gap),static_cast<std::int64_t>(o.diag_diff),o.diag_factor);
    });
    result.clustering_seconds=Elapsed(start);
    const size_t group_count=groups.size();
    // All span users have joined. Destroy views before releasing their owner.
    decltype(groups)().swap(groups);
    std::vector<Seed>().swap(seeds);
    std::uint64_t cluster_bytes=0;
    for(const auto& v:clusters) {
        result.clusters+=v->size();
        cluster_bytes+=v->capacity()*sizeof(pairwise_detail::MatchCluster);
        for(const auto& c:*v)cluster_bytes+=c.capacity()*sizeof(Seed);
    }
    observe("core-clustered-seeds-released",cluster_bytes);
    KernelProgress(o,"extension",0,group_count);
    start=KernelClock::now();
    pairwise_detail::SequenceViews::ManagerVariant rm{std::make_unique<pairwise_detail::SequenceViews::SequenceManager>(references)};
    pairwise_detail::SequenceViews::ManagerVariant qm{std::make_unique<pairwise_detail::SequenceViews::SequenceManager>(queries)};
    std::vector<pairwise_detail::AnchorVec> extended(group_count);
    std::vector<PairwiseStatistics> group_statistics(group_count);
    std::vector<pairwise_detail::AnchorLinkDetail::Statistics> link_statistics(group_count);
    ParallelKernel(group_count,o,[&](std::size_t i){
        KernelCheck(o,"extension");
        pairwise_detail::StatisticsScope statistics_scope(group_statistics[i]);
        auto anchors=pairwise_detail::AnchorLinkDetail::materializeClusterAnchors(*clusters[i],rm,qm);
        clusters[i].reset();
        const auto components=pairwise_detail::AnchorLinkDetail::splitAnchorComponents(anchors,&link_statistics[i]);
        for(const auto& component:components){KernelCheck(o,"extension");auto output=pairwise_detail::AnchorLinkDetail::linkAnchorRange(anchors,component.begin,component.end,rm,qm,&link_statistics[i]);extended[i].insert(extended[i].end(),std::make_move_iterator(output.begin()),std::make_move_iterator(output.end()));}

    });
    result.extension_seconds=Elapsed(start);
    std::uint64_t extended_bytes=0,packed_bytes=0;
    for(const auto& group:extended) {
        extended_bytes+=group.capacity()*sizeof(pairwise_detail::Anchor);
        for(const auto& anchor:group)packed_bytes+=anchor.cigar.capacity()*sizeof(pairwise_detail::CigarUnit);
    }
    observe("core-extended",0,extended_bytes,packed_bytes);
    for(size_t i=0;i<group_count;++i) {
        const auto& g=group_statistics[i];const auto& l=link_statistics[i];auto& t=result.statistics;
        t.global_ksw_calls+=g.global_ksw_calls;t.endpoint_ksw_calls+=g.endpoint_ksw_calls;
        t.global_ksw_seconds+=g.global_ksw_seconds;t.endpoint_ksw_seconds+=g.endpoint_ksw_seconds;
        t.link_closure_failures+=g.link_closure_failures;
        t.link_input_candidates+=g.link_input_candidates;
        t.link_consumed_candidates+=g.link_consumed_candidates;
        t.link_emitted_candidates+=g.link_emitted_candidates;
        t.link_candidate_checks+=l.candidate_checks;t.link_direct_attempts+=l.direct_ksw_calls;
        t.link_fallback_attempts+=l.fallback_ksw_calls;t.link_long_gap_rejections+=l.long_gap_rejections;
    }
    KernelProgress(o,"conflict-resolution",0,extended.size());
    start=KernelClock::now();
    size_t total_anchors=0;
    for(const auto& group:extended) {
        if(group.size()>SIZE_MAX-total_anchors)throw AlignmentError("pairwise anchor count overflow");
        total_anchors+=group.size();
    }
    pairwise_detail::AnchorVec anchors;
    anchors.reserve(total_anchors);
    for(auto& group:extended) {
        anchors.insert(anchors.end(),std::make_move_iterator(group.begin()),std::make_move_iterator(group.end()));
        pairwise_detail::AnchorVec().swap(group);
    }
    decltype(extended)().swap(extended);
    decltype(clusters)().swap(clusters);
    for(bool reference_side:{true,false}) {
        std::map<SequenceId,std::vector<size_t>> dimension;
        for(size_t index=0;index<anchors.size();++index) {
            const auto& anchor=anchors[index];
            dimension[reference_side?anchor.reference_id:anchor.query_id].push_back(index);
        }
        std::vector<std::vector<size_t>> tasks;
        for(auto& [id,list]:dimension){(void)id;tasks.push_back(std::move(list));}
        std::uint64_t task_bytes=tasks.capacity()*sizeof(std::vector<size_t>);
        for(const auto& task:tasks)task_bytes+=task.capacity()*sizeof(size_t);
        observe(reference_side?"core-reference-selection":"core-query-selection",0,
            anchors.capacity()*sizeof(pairwise_detail::Anchor),packed_bytes,task_bytes);
        ParallelKernel(tasks.size(),o,[&](std::size_t i){KernelCheck(o,"conflict-resolution");pairwise_detail::filterAnchorsByDpTreap(tasks[i],anchors,reference_side);});
    }
    result.selection_seconds=Elapsed(start);
    observe("core-selected",0,anchors.capacity()*sizeof(pairwise_detail::Anchor),packed_bytes);
    start=KernelClock::now();
    candidate_count=anchors.size();
    const auto output_count=static_cast<size_t>(std::count_if(anchors.begin(),anchors.end(),[&](const auto& a) {
        return !selected_only||(a.ref_selected&&a.qry_selected);
    }));
    result.records.reserve(output_count);
    for(auto& a:anchors) {
        KernelCheck(o,"conflict-resolution");
        const bool keep=!selected_only||(a.ref_selected&&a.qry_selected);
        auto record=ConvertAnchor(a,refs,qrys,keep);
        if(keep)result.records.push_back({std::move(record),a.ref_selected,a.qry_selected,a.aligned_base,a.alignment_length});
        pairwise_detail::Cigar_t().swap(a.cigar);
    }
    pairwise_detail::AnchorVec().swap(anchors);
    std::stable_sort(result.records.begin(),result.records.end(),[](const auto& a,const auto& b){const auto& x=a.record;const auto& y=b.record;return std::tie(x.query_id,x.query_begin,x.reference_id,x.reference_begin,x.strand,x.query_end,x.reference_end)<std::tie(y.query_id,y.query_begin,y.reference_id,y.reference_begin,y.strand,y.query_end,y.reference_end);});
    std::uint64_t output_cigar_bytes=0;
    for(const auto& item:result.records)output_cigar_bytes+=item.record.cigar.capacity()*sizeof(CigarOp);
    observe("core-output-converted",0,result.records.capacity()*sizeof(PairwiseAlignment),output_cigar_bytes);
    result.statistics.output_conversion_seconds=Elapsed(start);
    KernelProgress(o,"conflict-resolution",result.records.size(),result.records.size());
    return result;
}



PairwiseCoreResult AlignPairwiseCore(std::span<const SequenceRecord> references,std::span<const SequenceRecord> queries,std::span<const Seed> seeds,const PairwiseCoreOptions& o) {
    std::uint64_t candidate_count{};
    return AlignPairwiseCoreImpl(references,queries,std::vector<Seed>(seeds.begin(),seeds.end()),o,false,candidate_count);
}

void ResolveAlignmentRecords(std::vector<AlignmentRecord>& records, AlignmentSelection selection) {
    if (selection != AlignmentSelection::All && selection != AlignmentSelection::ReciprocalOneToOne)
        throw AlignmentError("record selection without sequences cannot link records");
    pairwise_detail::AnchorVec anchors;
    anchors.reserve(records.size());
    for (const auto& r : records) {
        if (r.reference_end <= r.reference_begin || r.query_end <= r.query_begin ||
            r.reference_end > INT64_MAX || r.query_end > INT64_MAX) {
            throw AlignmentError("pairwise record selection: invalid coordinate span");
        }
        std::uint64_t columns=0, matches=0, reference=0, query=0;
        for (const auto& op : r.cigar) {
            if (op.length==0 || columns>UINT64_MAX-op.length ||
                (op.operation!='=' && op.operation!='X' && op.operation!='I' && op.operation!='D')) {
                throw AlignmentError("pairwise record selection: invalid extended CIGAR");
            }
            columns+=op.length;
            if(op.operation=='=')matches+=op.length;
            if(op.operation!='I')reference+=op.length;
            if(op.operation!='D')query+=op.length;
        }
        if(reference!=r.reference_end-r.reference_begin || query!=r.query_end-r.query_begin)
            throw AlignmentError("pairwise record selection: CIGAR consumption mismatch");
        anchors.emplace_back(r.reference_id,r.reference_begin,reference,r.query_id,r.query_begin,query,
                             r.strand,columns,matches,pairwise_detail::Cigar_t{});
    }
    if(selection==AlignmentSelection::ReciprocalOneToOne) {
        for(bool reference_side : {true,false}) {
            std::map<SequenceId,std::vector<size_t>> groups;
            for(size_t i=0;i<anchors.size();++i)
                groups[reference_side?anchors[i].reference_id:anchors[i].query_id].push_back(i);
            for(auto& [id,group]:groups) {
                (void)id;
                pairwise_detail::filterAnchorsByDpTreap(group,anchors,reference_side);
            }
        }
    }
    std::map<SequenceId,bool> primary;
    std::vector<AlignmentRecord> selected;
    for(size_t i=0;i<records.size();++i) {
        if(selection==AlignmentSelection::ReciprocalOneToOne && !(anchors[i].ref_selected&&anchors[i].qry_selected))continue;
        auto& record=records[i];record.primary=!primary[record.query_id];primary[record.query_id]=true;
        selected.push_back(std::move(record));
    }
    records=std::move(selected);
}

PairwiseCoreOptions ConfiguredPairwiseCoreOptions(const AlignmentOptions& options) {
    PairwiseCoreOptions o{options.max_gap,options.diag_diff,options.diag_factor,options.min_cluster,options.worker_threads,options.interruption_callback,options.progress_callback,options.resident_bytes_callback};
    o.gap_fill_min_match=options.min_match;
    return o;
}

std::string PairwiseConfigurationText(const AlignmentOptions& options) {
    const auto o=ConfiguredPairwiseCoreOptions(options);
    std::ostringstream out;
    out << "pairwise.selection_weight=legacy-aligned-column-proxy-v1\n"
        << "pairwise.mainpath_work_budget_cap=" << o.gap_fill_work_budget << '\n';
    return out.str();
}

AlignmentResult AlignPairwiseFromSeeds(const std::vector<SequenceRecord>& refs,const std::vector<SequenceRecord>& queries,const AlignmentOptions& options,std::vector<Seed> seeds,RunStatistics stats) {
    if(options.break_length!=200 || options.max_dp_cells!=4000000 || options.match_score!=2 || options.mismatch_penalty!=4 || options.gap_open_penalty!=4 || options.gap_extend_penalty!=2)throw AlignmentError("legacy extension/scoring overrides are not applicable to the pairwise core");
    if (options.selection != AlignmentSelection::All &&
        options.selection != AlignmentSelection::ReciprocalOneToOne &&
        options.selection != AlignmentSelection::LinkedOneToOne)
        throw AlignmentError("invalid alignment selection");
    auto o=ConfiguredPairwiseCoreOptions(options);
    const bool selected_only = options.selection != AlignmentSelection::All;
    std::uint64_t candidate_count{};
    const auto seed_count=seeds.size();
    auto result=AlignPairwiseCoreImpl(refs,queries,std::move(seeds),o,selected_only,candidate_count);

    AlignmentResult output;
    output.statistics=std::move(stats);auto& s=output.statistics;
    s.pairwise=result.statistics;
    s.memory_observations.insert(s.memory_observations.end(),
        std::make_move_iterator(result.memory_observations.begin()),std::make_move_iterator(result.memory_observations.end()));
    s.reference_contigs=refs.size();s.query_contigs=queries.size();
    s.reference_bases=0;s.query_bases=0;s.reference_ambiguous_bases=0;s.query_ambiguous_bases=0;
    for(const auto& r:refs){if(r.size()>UINT64_MAX-s.reference_bases)throw AlignmentError("reference statistics overflow");s.reference_bases+=r.size();s.reference_ambiguous_bases+=static_cast<Length>(std::count(r.bases.begin(),r.bases.end(),'N'));}
    for(const auto& q:queries){if(q.size()>UINT64_MAX-s.query_bases)throw AlignmentError("query statistics overflow");s.query_bases+=q.size();s.query_ambiguous_bases+=static_cast<Length>(std::count(q.bases.begin(),q.bases.end(),'N'));}
    s.selected_seed_count=seed_count;s.chain_count=result.clusters;
    s.candidate_alignment_count=candidate_count;
    s.chaining_route="pairwise-cluster-best-chain+component-link+dual-dp-v1";
    s.chaining_requested_threads=options.worker_threads;s.chaining_worker_threads=result.workers;s.extension_worker_threads=result.workers;
    s.chain_and_extension_seconds=result.clustering_seconds+result.extension_seconds;
    s.conflict_resolution_seconds=result.selection_seconds;
    s.total_seconds=s.seed_seconds+s.chain_and_extension_seconds+s.conflict_resolution_seconds;
    std::map<SequenceId,size_t> primary;
    const auto eligible=[&](const auto& item){return !selected_only||(item.reference_selected&&item.query_selected);};
    for(size_t i=0;i<result.records.size();++i){
        const auto& item=result.records[i];
        if(eligible(item))primary.try_emplace(item.record.query_id,i);
    }
    for(size_t i=0;i<result.records.size();++i){
        auto& item=result.records[i];if(!eligible(item))continue;
        primary.try_emplace(item.record.query_id,i);
        item.record.primary=primary.at(item.record.query_id)==i;
        output.alignments.push_back(std::move(item.record));
    }
    // Selection rejection is measured before linking and never includes joins.
    s.conflict_rejected_alignment_count=s.candidate_alignment_count-output.alignments.size();
    if (options.selection == AlignmentSelection::LinkedOneToOne) {
        s.mainpath_link = LinkSelectedAlignmentRecords(refs, queries, output.alignments, options.link_gap_mode, o);
        s.total_seconds += s.mainpath_link.seconds;
    }
    s.alignment_count=output.alignments.size();
    return output;
}
} // namespace ramag

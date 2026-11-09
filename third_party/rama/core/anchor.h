// Modified for RaMA-G: library ownership/safety and deterministic N adaptation.
// Original algorithm rules and detailed changes: docs/centromere-rama-integration.md.
/*
 * Copyright [2024] [MALABZ_UESTC Pinglu Zhang]
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

 // Author: Pinglu Zhang
 // Contact: pingluzhang@outlook.com
 // Created: 2024-02-04
#pragma once

#include "gsacak.h"
#include "logging.h"
#include "utils.h"
#include "gsacak.h"
#include "rare_match.h"
#include "threadpool.h"
#include "RMQ.h"
#include <thread>
#include <mutex>
#include <cstdint>

namespace ramag::rama_source {



#define SAVE_DIR "save"
#define ANCHORFINDER_NAME "anchorfinder.bin"
#define FIRST_ANCHOR_NAME "first_anchor.csv"
#define FINAL_ANCHOR_NAME "final_anchor.csv"



struct Interval
{
	uint_t pos1;
	uint_t len1;

	uint_t pos2;
	uint_t len2;

	Interval(uint_t p1, uint_t l1, uint_t p2, uint_t l2) : pos1(p1), len1(l1), pos2(p2), len2(l2) {}
	Interval() : pos1(0), len1(0), pos2(0), len2(0) {}
};

using Intervals = std::vector<Interval>;

void saveIntervalsToCSV(const Intervals& intervals, const std::string& filename);

struct RecursionTrace {
    std::uint64_t node_id{}, parent_id{};
    uint_t depth{}, rare_level{};
    Interval interval;
    std::uint64_t candidate_count{}, selected_anchor_count{};
    const char* stop_reason{"no_matches"};
    const char* guard_status{"disabled"};
    std::uint64_t query_added{}, reference_added{}, cache_bytes{};
};
using RecursionTraceSink = std::function<void(const RecursionTrace&)>;

struct NodeCoverageGuardOptions {
    bool enabled{false};
    // Retained proposal-vector capacity; excludes ordinary solver workspace.
    std::uint64_t cache_limit_bytes{67108864};
};
// Resolve PRIVATE build definitions in one translation unit.
NodeCoverageGuardOptions DefaultNodeCoverageGuardOptions() noexcept;

struct NodeCoverageGuardSummary {
    double proposal_seconds_sum{}, selection_seconds{}, refill_seconds{};
    std::uint64_t baseline_nodes{}, cached_nodes{}, cache_exhausted_nodes{}, proposal_budget_exhausted_nodes{};
    std::uint64_t accepted_nodes{}, reused_subtrees{}, refill_searches{}, refill_nodes{};
    std::uint64_t refill_nonempty_nodes{}, cache_bytes_retained{};
};
struct NodeCoverageGain { std::uint64_t query_added{}, reference_added{}; };
// Inputs must be individually ordered, nonoverlapping exact-anchor chains.
// Coverage is actual paired positions, never the bounding box between anchors.
NodeCoverageGain MeasureNodeCoverageGain(const RareMatchPairs& baseline,
                                         const RareMatchPairs& proposal);

struct Anchor {
	uint_t depth;
	Anchor* parent;
    // Stable recursion ownership; no task may outlive the root.
	std::vector<Anchor*> children;
	RareMatchPairs rare_match_pairs;
    // Only this node's worker writes these fields; observers run after pool drain.
    Interval search_interval;
    uint_t rare_level{};
    std::uint64_t candidate_count{};
    const char* stop_reason{"no_matches"};
    const char* guard_status{"disabled"};
    std::uint64_t query_added{}, reference_added{}, cache_bytes{};
    std::uint64_t cache_allowance{};
    RareMatchPairs overlap_proposal;
    double proposal_seconds{};
    bool proposal_ready{false}, cache_exhausted{false}, suppress_proposals{false};
    bool proposal_budget_exhausted{false};
    const char* cache_failure_status{"cache_exhausted"};

	// Constructor initializes depth and parent
	Anchor(uint_t d = 0, Anchor* p = nullptr) : depth(d), parent(p) {}

	// Merges RareMatchPairs from this node and all descendants
	RareMatchPairs mergeRareMatchPairs() {
		RareMatchPairs mergedPairs;

		// Reserve memory to improve efficiency
		mergedPairs.reserve(100);
		uint_t rare_pair_index = 0;

		for (Anchor* child : children) {
			// Recursively merge child pairs
			RareMatchPairs childPairs = child->mergeRareMatchPairs();
			mergedPairs.insert(mergedPairs.end(), childPairs.begin(), childPairs.end());
			if (rare_pair_index < rare_match_pairs.size())
				mergedPairs.insert(mergedPairs.end(), rare_match_pairs.begin() + rare_pair_index, rare_match_pairs.begin() + rare_pair_index + 1);
			rare_pair_index++;
		}

		return mergedPairs;
	}

	// Destructor deletes all child Anchors
	~Anchor() {
		for (Anchor* child : children) {
			delete child;
		}
	}
};

// Select a top-down antichain. A failed proposal or exhausted cache freezes its
// complete old subtree; a node with no proposal may pass through to children.
// No sequence search or truth information is used here.
std::vector<Anchor*> SelectNodeCoverageRepairs(Anchor& root,
    const RareMatchPairs& baseline, const std::function<void()>& interrupt = {});

class AnchorFinder : public Serializable {
private:
	uint_t thread_num; // Indicates whether to use parallel processing

	std::string save_file_path; // Path to the directory where the AnchorFinder state is saved

	uint_t max_match_count; // Maximum number of rare matches to find
	bool initial_only{false}; // Separate root-only entry; original recursion defaults unchanged.
    RareCandidateSink initial_candidate_observer;
    bool corrected_recursive_boundaries{false};
    bool collect_trace{false};
    NodeCoverageGuardOptions coverage_guard_options;
    NodeCoverageGuardSummary coverage_guard_summary;
    bool coverage_guard_active{false}, capture_coverage_proposals{false};

	unsigned char* concat_data{}; // Concatenated sequence data

	uint_t concat_data_length; // Total length of the concatenated data
	uint_t first_seq_len; // Length of the first sequence
	uint_t second_seq_len; // Length of the second sequence

	uint_t* SA{}; // Suffix Array

	int_t* LCP{}; // Longest Common Prefix array

	int_da* DA{}; // Document Array indicating the origin sequence of each suffix

	uint_t* ISA{}; // Inverse Suffix Array

	LinearSparseTable rmq;
    std::unique_ptr<unsigned char[]> text_owner;
    std::vector<uint_t> sa_owner, isa_owner;
    std::vector<int_t> lcp_owner;
    std::vector<int_da> da_owner;
    std::function<void()> interrupt;
    bool has_unknown{}; // Range Minimum Query structure for LCP queries

	// Concatenates sequences from provided data
	void concatSequence(std::vector<SequenceInfo>& data);

	// Prints debug information including SA, LCP, and DA
	void printDebugInfo(const uint_t* SA, const int_t* LCP, const int_da* DA, uint_t concat_data_length);

	// Serialization of AnchorFinder state to an output stream
	void serialize(std::ostream& out) const override;

	// Deserialization of AnchorFinder state from an input stream
	void deserialize(std::istream& in) override;

	// Constructs the Inverse Suffix Array (ISA) for a given range
	void constructISA(uint_t start, uint_t end);

	// Constructs the ISA in parallel
	void constructISAParallel(uint_t thread_num);

	// Locates anchors using a given thread pool, recursive depth, and intervals
	void locateAnchor(ThreadPool& pool, uint_t depth, uint_t task_id, Anchor* root, Interval interval);

    void applyNodeCoverageGuard(ThreadPool& pool, Anchor& root,
                                const RareMatchPairs& baseline);
    void cacheNodeCoverageProposal(Anchor& node, const RareMatchPairs& raw,
                                   const RareMatchPairs& original);
	RareMatchPairs verifyAnchors(const RareMatchPairs& rare_match_pairs);


public:
	// Constructor initializes AnchorFinder with sequence data and optional parallel processing
	explicit AnchorFinder(std::vector<SequenceInfo>& data, std::string save_file_path, uint_t thread_num = 0, bool load_from_disk = false, bool save_to_disk = true, uint_t max_match_count = 100, std::function<void()> interrupt = {}, bool corrected_recursive_boundaries = false, NodeCoverageGuardOptions coverage_guard = DefaultNodeCoverageGuardOptions());

    const NodeCoverageGuardSummary& nodeCoverageGuardSummary() const noexcept {
        return coverage_guard_summary;
    }

	// Destructor cleans up allocated resources
	~AnchorFinder();

	// Launches the anchor searching process
	RareMatchPairs lanuchAnchorSearching(std::function<void(const RareMatchPairs&)> initial_observer = {},
                                         const RecursionTraceSink& trace_observer = {});
	RareMatchPairs findInitialAnchors(const RareCandidateSink& observer = {});

    // One region, no chaining and no recursion. Positions remain in the full
    // concatenated input. The caller owns bidirectional scheduling.
    void visitRegionCandidates(Interval region, bool supplemental_levels,
                               const RareCandidateSink& sink);
    std::pair<uint_t,uint_t> countFullOccurrences(uint_t reference_begin, uint_t length) const;

	// Converts RareMatchPairs to Intervals considering specified intervals
	static Intervals RareMatchPairs2Intervals(const RareMatchPairs& rare_match_pairs, Interval interval, uint_t fst_length);

    // Checked half-open partition for the opt-in directional recursion route.
    // Pair second_pos uses concatenated coordinates; Interval.pos2 is query-local.
    static Intervals SplitHalfOpenGaps(const RareMatchPairs& pairs, Interval interval,
                                      uint_t fst_length);

	// Converts a global index to a local index relative to concatenated sequences
	static uint_t indexFromGlogalToLocal(uint_t index, uint_t fst_length);
};


} // namespace ramag::rama_source


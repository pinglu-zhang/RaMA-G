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

#include "anchor.h"
#include <cstring>

namespace ramag::rama_source {


void saveIntervalsToCSV(const Intervals& intervals, const std::string& filename) {
	std::ofstream file(filename); // Opens the file for writing.
	if (!file.is_open()) { // Checks if the file is successfully opened.
		// Replace logger.error() with your logging mechanism or standard error output
		logger.error() << "Failed to open file: " << filename << std::endl; // Logs error if file cannot be opened.
		return; // Exits the function if file cannot be opened.
	}

	// Writes the header row of the CSV file.
	file << "Index,FirstStart,FirstLength,SecondStart,SecondLength\n";

	// Iterates through each pair in the Intervals and writes its data to the CSV file.
	for (size_t i = 0; i < intervals.size(); ++i) {
		const auto& interval = intervals[i]; // References the current interval pair.

		// Writes the index and details of the current interval pair to the file.
		file << i + 1 << ","
			<< interval.pos1 << "," // First interval's start position
			<< interval.len1 << "," // First interval's length
			<< interval.pos2 << "," // Second interval's start position
			<< interval.len2 << "\n"; // Second interval's length
	}

	file.close(); // Closes the file after writing.
	logger.info() << filename << " has been saved" << std::endl;
}

// Constructor for AnchorFinder class
AnchorFinder::AnchorFinder(std::vector<SequenceInfo>& data, std::string save_file_path, uint_t thread_num, bool load_from_disk, bool save_to_disk, uint_t max_match_count, std::function<void()> interruption) :
	save_file_path(save_file_path),
	thread_num(thread_num),
	max_match_count(max_match_count), interrupt(std::move(interruption)) {
	first_seq_len = data[0].seq_len;
	second_seq_len = data[1].seq_len;
	concatSequence(data); // Concatenate sequences from input data
	logger.info() << "The concated data length is " << concat_data_length << std::endl;


    // Value initialization is required by gSACA-K's sparse PLCP workspace.
    sa_owner.resize(concat_data_length); SA = sa_owner.data();
    lcp_owner.resize(concat_data_length); LCP = lcp_owner.data();
    da_owner.resize(concat_data_length); DA = da_owner.data();
    isa_owner.resize(concat_data_length); ISA = isa_owner.data();
    if (load_from_disk || save_to_disk)
        throw std::invalid_argument("RaMA integration has no persistent index route");
    {
		logger.info() << "The suffix array is constructing..." << std::endl;
		try {
			if (interrupt) interrupt();
            if (gsacak(concat_data, SA, LCP, DA, concat_data_length) < 0)
                throw std::runtime_error("gSACA-K construction failed");
            if (has_unknown) {
                // Unknown symbols are hard boundaries even for same-document LCP.
                std::vector<uint_t> run(concat_data_length);
                for (uint_t i = concat_data_length - 1; i > 0; --i)
                    if (canonical(concat_data[i - 1])) run[i - 1] = run[i] + 1;
                for (uint_t i = 1; i < concat_data_length; ++i)
                    LCP[i] = std::min(LCP[i], static_cast<int_t>(
                        std::min(run[SA[i - 1]], run[SA[i]])));
            }
			logger.info() << "The suffix array construction is finished!" << std::endl;
		}
		catch (std::exception& e) {
			logger.error() << "Error: " << e.what() << std::endl;
			logger.error() << "RaMA Exit!" << std::endl;
			throw;
		}

		logger.info() << "The sparse table is constructing..." << std::endl;
		try {
			this->rmq = LinearSparseTable(LCP, concat_data_length, thread_num);
			logger.info() << "The sparse table construction is finished!" << std::endl;
		}
		catch (std::exception& e) {
			logger.error() << "Error: " << e.what() << std::endl;
			logger.error() << "RaMA Exit!" << std::endl;
			throw;
		}


		if (thread_num)
			constructISAParallel(thread_num);
		else
			constructISA(0, concat_data_length - 1);


	}

	if (logger.isDebugEnabled()) {
		printDebugInfo(SA, LCP, DA, concat_data_length);
	}
}

// Destructor for AnchorFinder class
AnchorFinder::~AnchorFinder() = default;

// Concatenates sequences and prepares data for suffix array construction
void AnchorFinder::concatSequence(std::vector<SequenceInfo>& data) {
	// Calculate total length of concatenated string
	uint_t total_length = 0;
	for (uint_t i = 0; i < data.size(); i++) {
		total_length += data[i].seq_len + 1; // Include separator
	}

	total_length++;  // Add 1 for the terminating 0

	// Allocate memory for concatenated string
	text_owner = std::make_unique<unsigned char[]>(total_length);
    concat_data = text_owner.get();

	if (!concat_data) {
		logger.error() << "concat sequence could not allocate enough space" << std::endl;
		throw std::runtime_error("RaMA anchor construction/validation failed");
	}

	// Concatenate all strings with 1 as separator and terminate with 0
	uint_t index = 0;
	for (const auto& s : data) {
		for (uint_t j = 0; j < s.seq_len; ++j) {
            const auto c = static_cast<unsigned char>(s.sequence[j]);
            if (canonical(c)) concat_data[index + j] = c;
            else { concat_data[index + j] = index == 0 ? 2 : 3; has_unknown = true; }
        }
		index += s.seq_len;
		concat_data[index] = 1; // Use 1 as separator
		index++;
	}
	concat_data[total_length - 1] = 0; // Set terminating 0

	concat_data_length = total_length;
}


// Serializes the state of the AnchorFinder object to an output stream.
// This includes all fundamental configurations, concatenated sequence data,
// various arrays such as Suffix Array (SA), Longest Common Prefix (LCP), 
// Document Array (DA), Inverse Suffix Array (ISA), and the associated RMQ structure.
void AnchorFinder::serialize(std::ostream& out) const {
	// Save basic configuration numbers
	saveNumber(out, concat_data_length);
	saveNumber(out, first_seq_len);
	saveNumber(out, second_seq_len);

	// Save concatenated sequence data and associated arrays
	saveArray(out, concat_data, concat_data_length);
	saveArray(out, SA, concat_data_length);
	saveArray(out, LCP, concat_data_length);
	saveArray(out, DA, concat_data_length);
	saveArray(out, ISA, concat_data_length);

	// Serialize the RMQ structure
	this->rmq.serialize(out);
}


// Deserializes the state of the AnchorFinder object from an input stream.
// This method loads the AnchorFinder configuration, concatenated sequence data,
// various arrays (SA, LCP, DA, ISA), and reconstructs the RMQ state for future queries.
void AnchorFinder::deserialize(std::istream& in) {
	// Load basic configuration numbers
	loadNumber(in, concat_data_length);
	loadNumber(in, first_seq_len);
	loadNumber(in, second_seq_len);

	// Load concatenated sequence data and associated arrays
	loadArray(in, concat_data, concat_data_length);
	loadArray(in, SA, concat_data_length);
	loadArray(in, LCP, concat_data_length);
	loadArray(in, DA, concat_data_length);
	loadArray(in, ISA, concat_data_length);

	// Deserialize the RMQ structure and set the LCP array for RMQ queries
	this->rmq.deserialize(in);
	this->rmq.setLCP(LCP);
}


// Prints the debug information for SA, LCP, DA, and ISA arrays.
// This includes indexing and the values within each array for debugging purposes.
void AnchorFinder::printDebugInfo(const uint_t* SA, const int_t* LCP, const int_da* DA, uint_t concat_data_length) {
	std::stringstream s;

	// Print indices for reference
	s << " index: ";
	for (uint_t i = 0; i < concat_data_length; ++i) {
		s << std::setw(6) << std::left << i << " ";
	}
	logger.debug() << s.str() << std::endl;
	s.str("");

	// Print Suffix Array (SA) values
	s << "    SA: ";
	for (uint_t i = 0; i < concat_data_length; ++i) {
		s << std::setw(6) << std::left << SA[i] << " ";
	}
	logger.debug() << s.str() << std::endl;
	s.str("");

	// Print Longest Common Prefix (LCP) values
	s << "   LCP: ";
	for (uint_t i = 0; i < concat_data_length; ++i) {
		s << std::setw(6) << std::left << LCP[i] << " ";
	}
	logger.debug() << s.str() << std::endl;
	s.str("");

	// Print Document Array (DA) values
	s << "    DA: ";
	for (uint_t i = 0; i < concat_data_length; ++i) {
		s << std::setw(6) << std::left << DA[i] << " ";
	}
	logger.debug() << s.str() << std::endl;
	s.str("");

	// Print Inverse Suffix Array (ISA) values
	s << "   ISA: ";
	for (uint_t i = 0; i < concat_data_length; ++i) {
		s << std::setw(6) << std::left << ISA[i] << " ";
	}
	logger.debug() << s.str() << std::endl;
	s.str("");
}

// Constructs the Inverse Suffix Array (ISA) for the given range.
// ISA array is built by marking the index of each element of SA in the ISA array,
// which is useful for various string processing algorithms.
void AnchorFinder::constructISA(uint_t start, uint_t end) {
	for (uint_t i = start; i <= getMinValue(end, concat_data_length - 1); i++) {
		ISA[SA[i]] = i;
	}
}


// Constructs the Inverse Suffix Array (ISA) in parallel.
// This method divides the task of building the ISA array into smaller chunks
// and processes each chunk in parallel using a thread pool, improving performance on multicore systems.
void AnchorFinder::constructISAParallel(uint_t thread_num) {
	ThreadPool pool(thread_num); // Create a thread pool with a thread for each core
	const uint_t part_size = getMaxValue(std::ceil(static_cast<double>(concat_data_length) / std::max<uint_t>(thread_num, 1)), 1.0); // Calculate chunk size
	uint_t start = 0;
	uint_t end = part_size - 1;
	for (uint_t i = 0; i < concat_data_length; i += part_size) { // Divide work into chunks
		start = i;
		end = getMinValue(start + part_size - 1, concat_data_length - 1); // Ensure end does not exceed array bounds
		pool.enqueue([this, start, end]() { // Enqueue each chunk as a task
			this->constructISA(start, end);
			});
	}
	pool.waitAllTasksDone(); // Wait for all tasks to complete
}


// Initiates the process of searching for anchors in the sequences.
// Depending on the configuration, it either launches a parallel search using a thread pool
// or executes a single-threaded search.
RareMatchPairs AnchorFinder::lanuchAnchorSearching(std::function<void(const RareMatchPairs&)> initial_observer) {
	logger.info() << "Begin to search anchors" << std::endl;
	
	ThreadPool pool(thread_num); // Use thread pool for potential parallel execution
	uint_t depth = 0;
	auto owner = std::make_unique<Anchor>(depth);
    Anchor* root = owner.get(); // Create root anchor node
	Interval interval(0, first_seq_len, 0, second_seq_len); // Define interval
	uint_t task_id = 0;
	if (thread_num) {
		pool.enqueue([this, &pool, depth, task_id, root, interval]() {
			this->locateAnchor(pool, depth, task_id, root, interval);
			});
		pool.waitAllTasksDone();
	}
	else {
		locateAnchor(pool, depth, task_id, root, interval); // Fallback to sequential search
	}
	if (initial_observer) initial_observer(root->rare_match_pairs);


	// RareMatchPairs final_anchors = root->mergeRareMatchPairs(); // Merge rare match pairs from the root anchor
	RareMatchPairs final_anchors = verifyAnchors(root->mergeRareMatchPairs()); // Merge rare match pairs from the root anchor



	logger.info() << "Finish searching anchors" << std::endl;

	return final_anchors;
}

RareMatchPairs AnchorFinder::findInitialAnchors(const RareCandidateSink& observer) {
    ThreadPool pool(0);
    Anchor root;
    initial_only=true;
    initial_candidate_observer=observer;
    try {
        locateAnchor(pool,0,0,&root,Interval(0,first_seq_len,0,second_seq_len));
    } catch (...) {initial_only=false;initial_candidate_observer={};throw;}
    initial_only=false;
    initial_candidate_observer={};
    return std::move(root.rare_match_pairs);
}

std::pair<uint_t,uint_t> AnchorFinder::countFullOccurrences(uint_t begin, uint_t length) const {
    if (!length || begin > first_seq_len || length > first_seq_len-begin)
        throw std::out_of_range("invalid exact occurrence query");
    for (uint_t i=0;i<length;++i)
        if (!canonical(concat_data[begin+i])) throw std::invalid_argument("unknown in occurrence query");
    auto compare = [&](uint_t rank) {
        const uint_t pos=SA[rank];
        const uint_t available=concat_data_length-pos;
        const int c=std::memcmp(concat_data+pos,concat_data+begin,std::min(length,available));
        return c ? c : (available<length ? -1 : 0);
    };
    uint_t lo=0,hi=concat_data_length;
    while(lo<hi) { const auto m=lo+(hi-lo)/2; if(compare(m)<0)lo=m+1;else hi=m; }
    const uint_t first=lo; hi=concat_data_length;
    while(lo<hi) { const auto m=lo+(hi-lo)/2; if(compare(m)<=0)lo=m+1;else hi=m; }
    uint_t r=0,q=0;
    for (uint_t i=first;i<lo;++i) {
        const auto pos=SA[i];
        if(pos<=first_seq_len && length<=first_seq_len-pos) ++r;
        if(pos>first_seq_len && pos-first_seq_len-1<=second_seq_len &&
           length<=second_seq_len-(pos-first_seq_len-1)) ++q;
    }
    return {r,q};
}

void AnchorFinder::visitRegionCandidates(Interval region, bool supplemental_levels,
    const RareCandidateSink& sink) {
    if (region.pos1 > first_seq_len || region.len1 > first_seq_len-region.pos1 ||
        region.pos2 > second_seq_len || region.len2 > second_seq_len-region.pos2)
        throw std::out_of_range("RaMA candidate search region outside input");
    if (!region.len1 || !region.len2) return;
    if (interrupt) interrupt();
    const uint_t second_start = region.pos2 + first_seq_len + 1;
    std::vector<uint_t> indices;
    indices.reserve(static_cast<std::size_t>(region.len1)+region.len2);
    for (uint_t i=0;i<region.len1;++i) indices.push_back(ISA[region.pos1+i]);
    for (uint_t i=0;i<region.len2;++i) indices.push_back(ISA[second_start+i]);
    std::sort(indices.begin(),indices.end());
    std::vector<uint_t> sa(indices.size());
    std::vector<int_t> lcp(indices.size(),0);
    std::vector<int_da> da(indices.size());
    for (std::size_t i=0;i<indices.size();++i) {
        sa[i]=SA[indices[i]]; da[i]=DA[indices[i]];
        if (i) lcp[i]=rmq.queryMin(indices[i-1]+1,indices[i]);
    }
    RareMatchFinder finder(concat_data,sa,lcp,da,region.pos1,region.len1,
                           second_start,region.len2);
    finder.interruption=interrupt;
    finder.visitCandidates(max_match_count,supplemental_levels,sink);
}

// Launches the process of locating anchors within given intervals of two sequences.
// The method explores the given intervals, constructs new arrays based on the ISA,
// sorts them, and finds rare matches to determine new intervals for further exploration.
void AnchorFinder::locateAnchor(ThreadPool& pool, uint_t depth, uint_t task_id, Anchor* root, Interval interval) {
	// Log the start of a new task with its depth and task ID for debugging.
	logger.debug() << "Task " << task_id << " of depth " << depth << " begins" << std::endl;

	// Calculate new depth for recursive calls.
	uint_t new_depth = depth + 1;

	// Extract starting points and lengths from the intervals.
	uint_t first_seq_start = interval.pos1;
	uint_t fst_len = interval.len1;
	uint_t second_seq_start = interval.pos2 + first_seq_len + 1;
	uint_t scd_len = interval.len2;

	// Return early if either sequence segment is empty.
	if (fst_len == 0 || scd_len == 0) {
		return;
	}

	uint_t new_array_len = fst_len + scd_len;
	if (interrupt) interrupt();

	// Prepare arrays to hold new SA, LCP, and DA values.
	std::vector<uint_t> new_index_of_SA;
	new_index_of_SA.reserve(new_array_len);

	for (uint_t i = first_seq_start; i < first_seq_start + fst_len; i++) {
		new_index_of_SA.emplace_back(ISA[i]);
	}
	for (uint_t i = second_seq_start; i < second_seq_start + scd_len; i++) {
		new_index_of_SA.emplace_back(ISA[i]);
	}

	// Sort the new SA indices to maintain the order.
	std::sort(new_index_of_SA.begin(), new_index_of_SA.end());

	// Create and populate new SA, LCP, and DA arrays based on the sorted indices.
	std::vector<uint_t> new_SA;
	std::vector<int_t> new_LCP;
	std::vector<int_da> new_DA;
	/*new_SA.reserve(new_array_len);
	new_LCP.reserve(new_array_len);
	new_DA.reserve(new_array_len);

	if (!new_index_of_SA.empty()) {
		uint_t last_index = new_index_of_SA[0];
		new_SA.emplace_back(SA[last_index]);
		new_DA.emplace_back(DA[last_index]);
		new_LCP.emplace_back(0);

		for (size_t i = 1; i < new_index_of_SA.size(); ++i) {
			auto index = new_index_of_SA[i];
			new_SA.emplace_back(SA[index]);
			new_DA.emplace_back(DA[index]);
			new_LCP.emplace_back(rmq.queryMin(last_index + 1, index));
			last_index = index;
		}
	}*/
	new_SA.resize(new_array_len);
	new_LCP.resize(new_array_len);
	new_DA.resize(new_array_len);

	if (!new_index_of_SA.empty()) {
		uint_t last_index = new_index_of_SA[0];
		new_SA[0] = SA[last_index];
		new_DA[0] = DA[last_index];
		new_LCP[0] = 0;
// The recursion pool already owns the thread budget; no nested team.
		for (size_t i = 1; i < new_index_of_SA.size(); ++i) {
			auto index = new_index_of_SA[i];
			new_SA[i] = SA[index];
			new_DA[i] = DA[index];
			new_LCP[i] = rmq.queryMin(new_index_of_SA[i - 1] + 1, index);
		}

	}

	// Initialize RareMatchFinder and find optimal rare match pairs.
	RareMatchFinder rare_match_finder(concat_data, new_SA, new_LCP, new_DA, first_seq_start, fst_len, second_seq_start, scd_len);
	rare_match_finder.interruption = interrupt;
	RareMatchPairs optimal_pairs = rare_match_finder.findRareMatch(max_match_count,
        initial_only ? initial_candidate_observer : RareCandidateSink{});

	if (optimal_pairs.empty())
		return;

	// Convert rare match pairs to intervals for further exploration.
	Intervals rare_match_intervals = RareMatchPairs2Intervals(optimal_pairs, interval, this->first_seq_len);

	// Update the anchor's rare match pairs with the optimal ones found.
	root->rare_match_pairs = optimal_pairs;
	if (initial_only) return;

	// Recursively explore further intervals with new anchors.
	uint_t new_task_id = 0;

	for (const auto& new_interval : rare_match_intervals) {
		auto child = std::make_unique<Anchor>(new_depth, root);
		root->children.emplace_back(child.get());
		Anchor* new_anchor = child.release();
		// Parallel or sequential execution based on configuration.
		if (thread_num) {
			pool.enqueue([this, &pool, new_depth, new_task_id, new_anchor, new_interval]() {
				this->locateAnchor(pool, new_depth, new_task_id, new_anchor, new_interval);
				});
		}
		else {
			locateAnchor(pool, new_depth, new_task_id, new_anchor, new_interval);
		}
		new_task_id++;
	}
	// Log the end of the current task.
	logger.debug() << "Task " << task_id << " of depth " << depth << " ends" << std::endl;
	return;
}

// Converts rare match pairs to intervals for anchor finding. The function determines
// intervals between rare matches for further analysis.
Intervals AnchorFinder::RareMatchPairs2Intervals(const RareMatchPairs& rare_match_pairs, Interval interval, uint_t fst_length) {
	Intervals intervals; // Store the resulting intervals.

	// Return early if there are no rare match pairs.
	if (rare_match_pairs.empty()) {
		intervals.emplace_back(interval);
		return intervals;
	}


	// Initialize starting points and lengths for the first and second sequence segments.
	uint_t start1 = interval.pos1;
	uint_t seq1_length = interval.len1;
	uint_t start2 = interval.pos2 + fst_length + 1;
	uint_t seq2_length = interval.len2;
	uint_t i = 0;
	// Iterate over each rare match pair to calculate intervals.
	for (const auto& pair : rare_match_pairs) {
		// Calculate the start and end positions for each match in both sequences.
		uint_t match_start1 = pair.first_pos;
		uint_t match_start2 = pair.second_pos;
		uint_t match_end1 = match_start1 + pair.match_length - 1;
		uint_t match_end2 = match_start2 + pair.match_length - 1;

		// If current start points are less than or equal to match start points,
		// add the interval to the list and update the start points.
		if (start1 <= match_start1 && start2 <= match_start2) {
			intervals.emplace_back(Interval(start1, match_start1 - start1, indexFromGlogalToLocal(start2, fst_length), match_start2 - start2));
		}
		else {
			logger.error() << "There is conflict in final anchors" << std::endl;
		}

		start1 = match_end1 + 1;
		start2 = match_end2 + 1;
	}

	// Calculate and add the final interval after the last match.
	Interval end;
	if (start1 >= seq1_length) {
		start1--;
		end.pos1 = start1;
		end.len1 = 0;
	}
	else {
		end.pos1 = start1;
		end.len1 = seq1_length - start1;
	}

	if (start2 >= fst_length + 1 + seq2_length) {
		start2--;
		end.pos2 = indexFromGlogalToLocal(start2, fst_length);
		end.len2 = 0;
	}
	else {
		end.pos2 = indexFromGlogalToLocal(start2, fst_length);
		end.len2 = fst_length + 1 + seq2_length - start2;
	}

	intervals.emplace_back(end);

	return intervals;
}

// Converts global index to local index, adjusting for the concatenation offset.
uint_t AnchorFinder::indexFromGlogalToLocal(uint_t index, uint_t fst_length) {
	// Adjust index if it belongs to the second sequence.
	return index > fst_length ? index - fst_length - 1 : index;
}

RareMatchPairs AnchorFinder::verifyAnchors(const RareMatchPairs& rare_match_pairs) {
	if (rare_match_pairs.empty()) return rare_match_pairs;

	// Sort the input vector (if not already sorted)
	RareMatchPairs sorted_pairs(rare_match_pairs);
	std::sort(sorted_pairs.begin(), sorted_pairs.end());

	RareMatchPairs verified;
	RareMatchPair current = sorted_pairs[0];

	for (size_t i = 1; i < sorted_pairs.size(); ++i) {
		if (current.hasOverlap(sorted_pairs[i])) {
			logger.error() << "Error: Overlapping RareMatchPairs detected.\n";
			throw std::runtime_error("RaMA anchor construction/validation failed");
		}
		else if (current.isAdjacent(sorted_pairs[i])) {
			current.mergeWith(sorted_pairs[i]);
		}
		else {
			verified.push_back(current);
			current = sorted_pairs[i];
		}
	}

	// add the last processed RareMatchPair
	verified.push_back(current);

	return verified;
}

} // namespace ramag::rama_source


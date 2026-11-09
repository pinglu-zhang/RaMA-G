#include "centrolign_stitch/budget.hpp"
#ifndef centrolign_target_reachability_hpp
#define centrolign_target_reachability_hpp

#include "centrolign_stitch/topological_order.hpp"
#include "centrolign_stitch/utility.hpp"

#include <vector>
#include <utility>
#include <limits>

namespace centrolign_stitch_vendor {


// return whether each node can reach at least one of the targets
template<class Graph>
cg_budget::vector<bool> target_reachability(const Graph& graph, const cg_budget::vector<uint64_t>& targets) {
    
    cg_budget::vector<bool> reachable(graph.node_size(), false);
    
    for (auto node_id : targets) {
        reachable[node_id] = true;
    }
    
    auto order = topological_order(graph);
    for (uint64_t node_id : ReverseForEachAdapter<cg_budget::vector<uint64_t>>(order)) {
        for (auto next_id : graph.next(node_id)) {
            reachable[node_id] = reachable[node_id] || reachable[next_id];
        }
    }
    
    return reachable;
}

}

#endif /* centrolign_target_reachability_hpp */

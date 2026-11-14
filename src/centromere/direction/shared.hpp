#pragma once
#include "ramag/centromere_blocks.hpp"
#include <algorithm>
#include <memory_resource>
#include <span>
#include <vector>

// Private coordinate cells and exact-witness utilities shared by direction reviews.
namespace ramag::detail {
struct OrientationEndpoint {
  Position begin{},end{};
  Strand strand{Strand::Forward};
  std::uint64_t block_id{};
};
struct OrientationCoverage {Position begin{},end{};unsigned mask{};};

} // namespace ramag::detail

namespace ramag::detail {
// This is an existence index, not a vote count. Prefix maxima retain the fact
// that one actual exact match spans a window; merging intervals first would
// incorrectly turn several unrelated reference positions into one witness.
class DirectionWindowWitnessIndex {
 public:
  explicit DirectionWindowWitnessIndex(std::pmr::memory_resource* memory) : rows_(memory) {}
  void Add(Position begin,Position end) {
    if(begin>=end)throw AlignmentError("invalid exact orientation witness");
    rows_.push_back({begin,end});
  }
  void Finalize() {
    std::sort(rows_.begin(),rows_.end(),[](const auto& a,const auto& b) {
      return std::pair{a.begin,a.end}<std::pair{b.begin,b.end};
    });
    Position end=0;
    for(auto& row:rows_){end=std::max(end,row.end);row.end=end;}
  }
  bool Covers(Position begin,Position end) const {
    if(begin>=end)return false;
    auto it=std::upper_bound(rows_.begin(),rows_.end(),begin,
      [](auto position,const auto& row){return position<row.begin;});
    return it!=rows_.begin() && (--it)->end>=end;
  }
 private:
  std::pmr::vector<CentromereQueryInterval> rows_;
};

inline bool CoherentDirectionWindowDominates(bool enabled,bool single_exact_witness,
    Length own_length,Length opposite_length,std::uint64_t dominance_ratio) {
  // Retain both-direction ambiguity: an exact match in the other direction
  // that also spans this window prevents either side from passing this test.
  return enabled && single_exact_witness && dominance_ratio>=2 &&
    own_length && opposite_length<=own_length/dominance_ratio;
}
} // namespace ramag::detail

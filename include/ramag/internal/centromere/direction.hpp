#pragma once
#include "ramag/centromere_blocks.hpp"
#include <algorithm>
#include <memory_resource>
#include <span>
#include <vector>
#include "ramag/centromere_direction.hpp"
#include <unordered_map>
#include <tuple>
#include "ramag/internal/centromere/centromere_overlap_compatibility.hpp"
#include "ramag/internal/common/centromere.hpp"
#include <limits>


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


namespace ramag {
class DirectionOriginInterner {
  struct Key {
    const DirectionOriginNode* previous;
    DirectionOrigin value;
    bool operator==(const Key& b) const {
      const auto& a=value;const auto& v=b.value;
      return previous==b.previous && a.task_id==v.task_id && a.depth==v.depth &&
        a.rare_level==v.rare_level && a.reference_interval_count==v.reference_interval_count &&
        a.query_interval_count==v.query_interval_count && a.first_hit_level==v.first_hit_level;
    }
  };
  struct Hash {
    std::size_t operator()(const Key& k) const {
      auto h=std::hash<const void*>{}(k.previous);
      for(auto v:{k.value.task_id,static_cast<std::uint64_t>(k.value.depth),
          static_cast<std::uint64_t>(k.value.rare_level),static_cast<std::uint64_t>(k.value.reference_interval_count),
          static_cast<std::uint64_t>(k.value.query_interval_count),static_cast<std::uint64_t>(k.value.first_hit_level)})
        h^=std::hash<std::uint64_t>{}(v)+0x9e3779b9U+(h<<6)+(h>>2);
      return h;
    }
  };
  std::unordered_map<Key,std::weak_ptr<const DirectionOriginNode>,Hash> nodes_;
  std::size_t created_{};
 public:
  std::size_t size() const {return created_;}
  std::size_t ClearCache() {
    const auto count=nodes_.size();
    decltype(nodes_) empty;nodes_.swap(empty);return count;
  }
  bool Append(DirectionOrigins& history,const DirectionOrigin& value) {
    const Key key{history.tail_.get(),value};
    auto it=nodes_.find(key);
    if(it!=nodes_.end()) if(auto p=it->second.lock()) {history.tail_=std::move(p);return false;}
    auto p=std::make_shared<DirectionOriginNode>(DirectionOriginNode{value,history.tail_,history.size()+1});
    nodes_.insert_or_assign(key,p);history.tail_=std::move(p);++created_;return true;
  }
};

namespace direction_detail {
inline std::uint64_t Mix(std::uint64_t x) {
  x^=x>>30;x*=0xbf58476d1ce4e5b9ULL;x^=x>>27;x*=0x94d049bb133111ebULL;return x^(x>>31);
}
inline std::uint64_t AnchorHash(Position r,Position q,Length n,Strand d) {
  return Mix(r)^Mix(q+0x9e3779b97f4a7c15ULL)^Mix(n+static_cast<std::uint64_t>(d)*0x517cc1b727220a95ULL);
}
// Only stable anchor IDs are stored in the table. Keys are read from the
// anchor array, so there is no duplicate coordinate tuple or tree node.
class AnchorLookup {
  const std::vector<DirectionAnchor>& anchors_;
  std::vector<std::uint64_t> slots_;
  std::size_t size_{};
 public:
  explicit AnchorLookup(const std::vector<DirectionAnchor>& anchors):anchors_(anchors) {}
  std::uint64_t bytes() const {return slots_.capacity()*sizeof(std::uint64_t);}
  std::size_t Find(Position r,Position q,Length n,Strand d) const {
    if(slots_.empty())return kNoDirectionNode;
    auto slot=static_cast<std::size_t>(AnchorHash(r,q,n,d))&(slots_.size()-1);
    while(slots_[slot]) {
      const auto index=static_cast<std::size_t>(slots_[slot]-1);const auto& a=anchors_[index];
      if(a.reference_begin==r && a.query_begin==q && a.length==n && a.strand==d)return index;
      slot=(slot+1)&(slots_.size()-1);
    }
    return kNoDirectionNode;
  }
  template<class Guard> void Insert(std::size_t index,Guard& guard) {
    if(slots_.empty() || size_>=slots_.size()/10*7) {
      const auto count=slots_.empty() ? 1024 : slots_.size()*2;
      guard.Add(count*sizeof(std::uint64_t));
      std::vector<std::uint64_t> grown(count,0);
      for(auto id:slots_) if(id) {
        const auto& a=anchors_[static_cast<std::size_t>(id-1)];
        auto slot=static_cast<std::size_t>(AnchorHash(a.reference_begin,a.query_begin,a.length,a.strand))&(count-1);
        while(grown[slot])slot=(slot+1)&(count-1);
        grown[slot]=id;
      }
      guard.Release(bytes());slots_.swap(grown);
    }
    const auto& a=anchors_[index];
    auto slot=static_cast<std::size_t>(AnchorHash(a.reference_begin,a.query_begin,a.length,a.strand))&(slots_.size()-1);
    while(slots_[slot])slot=(slot+1)&(slots_.size()-1);
    slots_[slot]=static_cast<std::uint64_t>(index)+1;++size_;
  }
};
struct CompactBox {
  std::uint32_t reference_begin{},reference_end{},query_begin{},query_end{};
  bool operator==(const CompactBox&) const = default;
  explicit CompactBox(const DirectionBox& b):reference_begin(static_cast<std::uint32_t>(b.reference_begin)),
    reference_end(static_cast<std::uint32_t>(b.reference_end)),query_begin(static_cast<std::uint32_t>(b.query_begin)),
    query_end(static_cast<std::uint32_t>(b.query_end)) {}
  CompactBox()=default;
  DirectionBox Box() const {return {reference_begin,reference_end,query_begin,query_end};}
};
// Coordinates fit the already checked source range. Compact slots avoid one
// red-black tree allocation and four duplicated 64-bit bounds per window.
class WindowSet {
  std::vector<CompactBox> slots_;
  std::size_t size_{};
  static std::uint64_t Hash(const CompactBox& b) {
    return Mix((static_cast<std::uint64_t>(b.reference_begin)<<32)|b.reference_end)^
      Mix(((static_cast<std::uint64_t>(b.query_begin)<<32)|b.query_end)+0x9e3779b97f4a7c15ULL);
  }
 public:
  std::size_t size() const {return size_;}
  std::uint64_t bytes() const {return slots_.capacity()*sizeof(CompactBox);}
  template<class Guard> void Insert(const DirectionBox& box,Guard& guard) {
    const CompactBox value(box);
    if(!slots_.empty()) {
      auto slot=static_cast<std::size_t>(Hash(value))&(slots_.size()-1);
      while(slots_[slot].reference_end) {
        if(slots_[slot]==value)return;
        slot=(slot+1)&(slots_.size()-1);
      }
    }
    if(slots_.empty() || size_>=slots_.size()/10*7) {
      const auto count=slots_.empty() ? 1024 : slots_.size()*2;
      guard.Add(count*sizeof(CompactBox));std::vector<CompactBox> grown(count);
      for(const auto& b:slots_)if(b.reference_end) {
        auto slot=static_cast<std::size_t>(Hash(b))&(count-1);
        while(grown[slot].reference_end)slot=(slot+1)&(count-1);
        grown[slot]=b;
      }
      guard.Release(slots_.capacity()*sizeof(CompactBox));slots_.swap(grown);
    }
    auto slot=static_cast<std::size_t>(Hash(value))&(slots_.size()-1);
    while(slots_[slot].reference_end) {
      if(slots_[slot]==value)return;
      slot=(slot+1)&(slots_.size()-1);
    }
    slots_[slot]=value;++size_;
  }
  std::vector<CompactBox> Finish() {
    std::sort(slots_.begin(),slots_.end(),[](const auto& a,const auto& b) {
      if(!a.reference_end || !b.reference_end)return a.reference_end && !b.reference_end;
      return std::tie(a.reference_begin,a.reference_end,a.query_begin,a.query_end)<
        std::tie(b.reference_begin,b.reference_end,b.query_begin,b.query_end);
    });
    slots_.resize(size_);return std::move(slots_);
  }
};
struct CoveredInterval {
  std::uint64_t diagonal;
  Position begin,end;
};
inline std::uint64_t Diagonal(const DirectionAnchor& a,Length qlength) {
  const auto u=a.strand==Strand::Forward ? a.query_begin : qlength-a.query_begin-a.length;
  const auto delta=a.reference_begin>=u ? a.reference_begin-u : u-a.reference_begin;
  // Sequence analysis checks the pinned source range (< 2^30) before entry.
  return delta | (static_cast<std::uint64_t>(a.reference_begin>=u)<<62) |
    (static_cast<std::uint64_t>(a.strand==Strand::Reverse)<<63);
}
template<class Guard> std::vector<CoveredInterval> Coverage(
    std::span<const DirectionAnchor> anchors,Length qlength,Guard& guard) {
  guard.Add(anchors.size()*sizeof(CoveredInterval));
  std::vector<CoveredInterval> out;out.reserve(anchors.size());
  for(const auto& a:anchors)out.push_back({Diagonal(a,qlength),a.reference_begin,a.reference_begin+a.length});
  std::sort(out.begin(),out.end(),[](const auto& a,const auto& b){return std::tie(a.diagonal,a.begin,a.end)<std::tie(b.diagonal,b.begin,b.end);});
  std::size_t n=0;
  for(const auto& interval:out) {
    if(n && out[n-1].diagonal==interval.diagonal && out[n-1].end>=interval.begin)
      out[n-1].end=std::max(out[n-1].end,interval.end);
    else out[n++]=interval;
  }
  out.resize(n);
  if(n<out.capacity()/2) {
    guard.Add(n*sizeof(CoveredInterval));
    std::vector<CoveredInterval> compact(out.begin(),out.end());
    guard.Release(out.capacity()*sizeof(CoveredInterval));out.swap(compact);
  }
  return out;
}
inline DirectionAnchor Mapping(const CoveredInterval& interval,Length qlength) {
  DirectionAnchor a;
  a.reference_begin=interval.begin;a.length=interval.end-interval.begin;
  a.strand=(interval.diagonal>>63) ? Strand::Reverse : Strand::Forward;
  const auto delta=interval.diagonal&((1ULL<<62)-1);
  const auto u=(interval.diagonal&(1ULL<<62)) ? interval.begin-delta : interval.begin+delta;
  a.query_begin=a.strand==Strand::Forward ? u : qlength-u-a.length;return a;
}
inline bool Contains(std::span<const CoveredInterval> intervals,const DirectionAnchor& a,Length qlength) {
  const auto diagonal=Diagonal(a,qlength);
  auto it=std::upper_bound(intervals.begin(),intervals.end(),std::pair{diagonal,a.reference_begin},
    [](const auto& key,const auto& interval){return key < std::pair{interval.diagonal,interval.begin};});
  if(it==intervals.begin())return false;
  --it;return it->diagonal==diagonal && it->end>=a.reference_begin+a.length;
}
} // namespace direction_detail
} // namespace ramag


namespace ramag::centromere_planning_detail {
constexpr auto none=std::numeric_limits<std::size_t>::max();
struct Node {
  CentromereAnchor anchor;
  Strand strand;
  bool overlap_variant{};
  bool uniform_density{true};
  detail::AnchorClipWitnesses clipping;
  std::size_t predecessor{none},first{};
  Length support{};
  double weight{};
  Length dp_support{};
  double dp_score{};
  std::size_t dp_predecessor{none};
};
class Budget final : public std::pmr::memory_resource {
  std::uint64_t limit_,live_;
  void* do_allocate(std::size_t bytes,std::size_t alignment) override {
    Account(bytes);
    try{return std::pmr::new_delete_resource()->allocate(bytes,alignment);}
    catch(...){live_-=bytes;throw;}
  }
  void do_deallocate(void* p,std::size_t bytes,std::size_t alignment) override {
    std::pmr::new_delete_resource()->deallocate(p,bytes,alignment);live_-=bytes;
  }
  bool do_is_equal(const std::pmr::memory_resource& r) const noexcept override {return this==&r;}
 public:
  std::uint64_t peak{};
  Budget(std::uint64_t limit,std::uint64_t retained):limit_(limit),live_(retained),peak(retained) {
    if(retained>limit)throw detail::CentromereBudgetExceeded("joint root evidence exceeds planner memory budget");
  }
  void Account(std::uint64_t bytes) {
    if(bytes>limit_ || live_>limit_-bytes)
      throw detail::CentromereBudgetExceeded("joint planner live capacity budget exceeded: retained="+std::to_string(live_)+
        " requested="+std::to_string(bytes)+" limit="+std::to_string(limit_)+"; no candidates were truncated");
    live_+=bytes;peak=std::max(peak,live_);
  }
  void Release(std::uint64_t bytes){live_-=bytes;}
  std::uint64_t Live() const {return live_;}
  void Restore(std::uint64_t checkpoint){live_=checkpoint;}
};
template<class T>using Vec=std::pmr::vector<T>;
template<class T>void Reserve(std::vector<T>& v,std::size_t size,Budget& budget) {
  const auto old=v.capacity()*sizeof(T);budget.Account(size*sizeof(T));
  v.reserve(size);budget.Release(old);
}

void Check(const CentromereBlockOptions& options);
void OrientationEvidence(const Vec<Node>& nodes,Length qlength,std::string_view query,
    CentromereBlockResult& out,Budget& budget,const CentromereBlockOptions& options);
}

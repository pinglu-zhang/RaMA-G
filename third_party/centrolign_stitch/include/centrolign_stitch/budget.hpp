#pragma once
// Local resource instrumentation only. Container ordering and numerical
// recurrences are upstream's. Counts are requested allocation payload bytes,
// not allocator metadata/RSS or caller-owned inputs and callback allocations.
#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <deque>
#include <limits>
#include <map>
#include <memory>
#include <new>
#include <queue>
#include <set>
#include <stack>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <vector>
namespace cg_budget {
struct Limit { const char* reason; };
struct Context {
  std::uint64_t limit{}, live{}, peak{}, ticks{}, states{};
  double seconds{};
  const char* phase{"stitch_gap_begin"};
  std::chrono::steady_clock::time_point start{std::chrono::steady_clock::now()};
  void (*callback)(void*){};
  void* callback_data{};
  void check() {
    if (callback) callback(callback_data);
    if (std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count() > seconds)
      throw Limit{"gap_time_budget"};
  }
  void acquire(std::size_t bytes) {
    check();
    if (bytes > limit || live > limit-bytes) throw Limit{"gap_memory_budget"};
    live += bytes; peak=std::max(peak,live);
  }
  void release(std::size_t bytes) noexcept { assert(bytes<=live); live-=bytes; }
};
inline thread_local Context* active=nullptr;
struct Scope {
  Context* previous;
  explicit Scope(Context& c) noexcept : previous(active) { active=&c; }
  ~Scope() { active=previous; }
};
inline void tick() {
  if (active && ((++active->ticks & 255U)==0)) active->check();
}
inline void visited_state() { if (active) ++active->states; }
template<class T> struct Allocator {
  using value_type=T;
  using propagate_on_container_move_assignment=std::true_type;
  using propagate_on_container_swap=std::true_type;
  Context* owner;
  Allocator() noexcept : owner(active) {}
  template<class U> Allocator(const Allocator<U>& o) noexcept : owner(o.owner) {}
  [[nodiscard]] T* allocate(std::size_t n) {
    if (n > std::numeric_limits<std::size_t>::max()/sizeof(T)) throw Limit{"gap_memory_budget"};
    const auto bytes=n*sizeof(T);
    if(owner)owner->acquire(bytes);
    try { return std::allocator<T>{}.allocate(n); }
    catch (const std::bad_alloc&) { if(owner)owner->release(bytes);throw Limit{"allocation_failure"}; }
  }
  void deallocate(T* p,std::size_t n) noexcept {
    std::allocator<T>{}.deallocate(p,n); if(owner)owner->release(n*sizeof(T));
  }
  template<class U> bool operator==(const Allocator<U>& o)const noexcept{return owner==o.owner;}
  template<class U> bool operator!=(const Allocator<U>& o)const noexcept{return owner!=o.owner;}
};
template<class T> using vector=std::vector<T,Allocator<T>>;
template<class T> using deque=std::deque<T,Allocator<T>>;
template<class T,class C=deque<T>> using queue=std::queue<T,C>;
template<class T,class C=deque<T>> using stack=std::stack<T,C>;
template<class K,class V,class H=std::hash<K>,class E=std::equal_to<K>>
using unordered_map=std::unordered_map<K,V,H,E,Allocator<std::pair<const K,V>>>;
template<class K,class H=std::hash<K>,class E=std::equal_to<K>>
using unordered_set=std::unordered_set<K,H,E,Allocator<K>>;
template<class K,class C=std::less<K>> using set=std::set<K,C,Allocator<K>>;
template<class K,class V,class C=std::less<K>> using map=std::map<K,V,C,Allocator<std::pair<const K,V>>>;
} // namespace cg_budget

#include "rama_wfa_memory.hpp"
#include "centromere_detail.hpp"
#include <cstdlib>
#include <new>
#include <stdexcept>
#include <limits>
namespace {
thread_local ramag::detail::WfaMemoryScope* active{};
void CheckAllocation(std::size_t n,std::size_t replaced=0) {
    if (active && active->interruption) active->interruption();
    if (!active) throw std::logic_error("WFA outside owned scope");
    const auto retained=active->live_bytes-replaced;
    if (n>std::numeric_limits<std::uint64_t>::max()-retained ||
        (active->limit_bytes && (retained>active->limit_bytes || n>active->limit_bytes-retained)))
        throw ramag::detail::CentromereBudgetExceeded("wfa_memory_budget");
}
void Remember(void* p,std::size_t n) {
    if (!p) throw std::bad_alloc();
    if (!active) { std::free(p); throw std::logic_error("WFA outside owned scope"); }
    try { active->allocations.emplace(p,n);active->live_bytes+=n; }
    catch (...) { std::free(p); throw; }
}
}
namespace ramag::detail {
WfaMemoryScope::WfaMemoryScope(std::function<void()> callback,std::uint64_t limit):
    previous(active),interruption(std::move(callback)),limit_bytes(limit) { active=this; }
WfaMemoryScope::~WfaMemoryScope() {
    for (const auto& entry:allocations) std::free(entry.first);
    active=previous;
}
}
extern "C" {
void* ramag_wfa_malloc(std::size_t n) {
    n=n ? n : 1;CheckAllocation(n);
    void* p=std::malloc(n); Remember(p,n); return p;
}
void* ramag_wfa_calloc(std::size_t n,std::size_t s) {
    n=n ? n : 1;s=s ? s : 1;
    if(n>std::numeric_limits<std::size_t>::max()/s)throw std::bad_alloc();
    CheckAllocation(n*s);
    void* p=std::calloc(n,s); Remember(p,n*s); return p;
}
void ramag_wfa_free(void* p) {
    if (active) {
        const auto found=active->allocations.find(p);
        if(found!=active->allocations.end()) {active->live_bytes-=found->second;active->allocations.erase(found);}
    }
    std::free(p);
}
void* ramag_wfa_realloc(void* p,std::size_t n) {
    if (!p) return ramag_wfa_malloc(n);
    if (!n) { ramag_wfa_free(p); return nullptr; }
    if(!active || !active->allocations.contains(p))throw std::logic_error("unowned WFA realloc");
    const auto old=active->allocations.at(p);
    CheckAllocation(n,old);
    auto node=active->allocations.extract(p);
    void* q=std::realloc(p,n);
    if (!q) { active->allocations.insert(std::move(node)); throw std::bad_alloc(); }
    node.key()=q;node.mapped()=n;active->allocations.insert(std::move(node));
    active->live_bytes=active->live_bytes-old+n;
    return q;
}
[[noreturn]] void ramag_wfa_exit(int code) {
    throw std::runtime_error("pinned WFA2 reported fatal error "+std::to_string(code));
}
}

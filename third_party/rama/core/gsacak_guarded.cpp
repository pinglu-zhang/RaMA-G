// Compile the original C algorithm as C++ solely for exception-safe allocation.
#include <algorithm>
#include <cstdlib>
#include <new>
#include <stdexcept>
#include <vector>
#include "gsacak.h"
namespace {
struct AllocationScope;
thread_local AllocationScope* allocations = nullptr;
struct AllocationScope {
    AllocationScope* previous = allocations;
    std::vector<void*> pointers;
    AllocationScope() { allocations = this; }
    ~AllocationScope() {
        for (void* p : pointers) std::free(p);
        allocations = previous;
    }
};
void* allocate(std::size_t n) {
    void* p = std::malloc(n);
    if (!p) throw std::bad_alloc();
    try { allocations->pointers.push_back(p); }
    catch (...) { std::free(p); throw; }
    return p;
}
void release(void* p) {
    auto& v = allocations->pointers;
    auto it = std::find(v.begin(), v.end(), p);
    if (it != v.end()) v.erase(it);
    std::free(p);
}
}
#define malloc allocate
#define free release
#define exit(code) throw std::runtime_error("RaMA gSACA-K stack exhausted")
#include "gsacak.c"
#undef malloc
#undef free
#undef exit

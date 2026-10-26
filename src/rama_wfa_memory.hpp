#pragma once
#include <cstddef>
#include <unordered_map>
#include <functional>
#include <cstdint>
namespace ramag::detail {
// All allocations made by the pinned C library belong to one call/worker.
struct WfaMemoryScope {
    WfaMemoryScope* previous{};
    std::unordered_map<void*,std::size_t> allocations;
    std::function<void()> interruption;
    std::uint64_t limit_bytes{}, live_bytes{};
    explicit WfaMemoryScope(std::function<void()> callback = {}, std::uint64_t limit = 0);
    ~WfaMemoryScope();
    WfaMemoryScope(const WfaMemoryScope&)=delete;
};
}

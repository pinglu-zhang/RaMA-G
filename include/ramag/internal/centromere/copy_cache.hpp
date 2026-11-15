#pragma once

#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

namespace ramag::detail {

enum class CopyCountTarget : std::uint8_t { Reference, Query };

inline char CopyCountComplement(char c) noexcept {
  switch (c) {
  case 'A': return 'T';
  case 'C': return 'G';
  case 'G': return 'C';
  case 'T': return 'A';
  default: return 'N';
  }
}

inline char CopyCountPatternBase(std::string_view source, std::size_t begin,
                                std::size_t length, bool reverse,
                                std::size_t offset) noexcept {
  return reverse ? CopyCountComplement(source[begin + length - 1 - offset])
                 : source[begin + offset];
}

struct CopyCountPatternHash {
  std::uint64_t operator()(std::string_view source, std::size_t begin,
                           std::size_t length, bool reverse) const noexcept {
    std::uint64_t hash = UINT64_C(14695981039346656037);
    for (std::size_t i = 0; i < length; ++i) {
      hash ^= static_cast<unsigned char>(
          CopyCountPatternBase(source, begin, length, reverse, i));
      hash *= UINT64_C(1099511628211);
    }
    return hash;
  }
};

// Pure memoization for overlapping occurrence counts. The source string and
// both target sequences must remain immutable for this object's lifetime.
// The callback must be the same pure counting operation for a given target,
// and must not reenter this cache. No pattern or occurrence list is retained:
// slots contain stable coordinates into source, not views into the RC scratch.
// Hash collisions/evictions only repeat the original count operation.
template <class Hasher = CopyCountPatternHash> class CopyCountCache {
  struct Slot {
    std::uint64_t hash{};
    std::uint64_t count{};
    std::size_t begin{};
    std::size_t length{};
    CopyCountTarget target{};
    bool reverse{};
    bool valid{};
  };

public:
  static constexpr std::size_t kDefaultSlots = 65536;

  explicit CopyCountCache(const std::string &source,
                          std::size_t capacity = kDefaultSlots,
                          Hasher hasher = Hasher{})
      : source_(source), capacity_(capacity), hasher_(std::move(hasher)) {
    if (capacity_ > std::numeric_limits<std::size_t>::max() / sizeof(Slot))
      throw std::length_error("copy count cache slot allocation overflow");
    if (capacity_)
      slots_ = std::make_unique<Slot[]>(capacity_);
  }
  CopyCountCache(std::string &&, std::size_t = kDefaultSlots,
                 Hasher = Hasher{}) = delete;
  CopyCountCache(const std::string &&, std::size_t = kDefaultSlots,
                 Hasher = Hasher{}) = delete;
  CopyCountCache(const CopyCountCache &) = delete;
  CopyCountCache &operator=(const CopyCountCache &) = delete;
  CopyCountCache(CopyCountCache &&) = delete;
  CopyCountCache &operator=(CopyCountCache &&) = delete;

  template <class Counter>
  std::uint64_t Get(CopyCountTarget target, std::size_t begin,
                    std::size_t length, bool reverse, Counter &&counter) {
    if (begin > source_.size() || length > source_.size() - begin)
      throw std::out_of_range("copy count cache pattern outside source");
    if (!capacity_)
      return counter(Pattern(begin, length, reverse));
    const auto hash = hasher_(source_, begin, length, reverse) ^
                      (target == CopyCountTarget::Query
                           ? UINT64_C(0x9e3779b97f4a7c15)
                           : UINT64_C(0));
    auto &slot = slots_[static_cast<std::size_t>(hash % capacity_)];
    if (slot.valid && slot.target == target && slot.hash == hash &&
        slot.length == length && Equal(slot, begin, reverse))
      return slot.count;
    // Store only after the unchanged counting callback completed successfully.
    const std::uint64_t count = counter(Pattern(begin, length, reverse));
    slot = {hash, count, begin, length, target, reverse, true};
    return count;
  }

  static constexpr std::size_t SlotBytes() noexcept { return sizeof(Slot); }
  std::size_t SlotCount() const noexcept { return capacity_; }
  std::size_t SlotStorageBytes() const noexcept {
    return capacity_ * sizeof(Slot);
  }
  std::size_t ScratchCapacity() const noexcept { return scratch_capacity_; }

private:
  bool Equal(const Slot &slot, std::size_t begin, bool reverse) const noexcept {
    for (std::size_t i = 0; i < slot.length; ++i)
      if (CopyCountPatternBase(source_, slot.begin, slot.length, slot.reverse, i) !=
          CopyCountPatternBase(source_, begin, slot.length, reverse, i))
        return false;
    return true;
  }

  std::string_view Pattern(std::size_t begin, std::size_t length, bool reverse) {
    if (!reverse)
      return source_.substr(begin, length);
    if (length > scratch_capacity_) {
      // Release the previous buffer before growth so owned scratch payload
      // never exceeds the longest requested reverse-complement pattern.
      scratch_.reset();
      scratch_capacity_ = 0;
      scratch_.reset(new char[length]);
      scratch_capacity_ = length;
    }
    for (std::size_t i = 0; i < length; ++i)
      scratch_[i] = CopyCountPatternBase(source_, begin, length, true, i);
    return length ? std::string_view(scratch_.get(), length) : std::string_view{};
  }

  std::string_view source_;
  std::size_t capacity_{};
  Hasher hasher_;
  std::unique_ptr<Slot[]> slots_;
  std::unique_ptr<char[]> scratch_;
  std::size_t scratch_capacity_{};
};

} // namespace ramag::detail

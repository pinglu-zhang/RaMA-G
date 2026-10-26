// Source logging is intentionally inert: never install a host-global logger.
#pragma once
#include "utils.h"
namespace ramag::rama_source {
struct NullStream {
    template<class T> const NullStream& operator<<(const T&) const { return *this; }
    const NullStream& operator<<(std::ostream& (*)(std::ostream&)) const { return *this; }
};
struct NullLogger {
    NullStream info() const { return {}; }
    NullStream error() const { return {}; }
    NullStream debug() const { return {}; }
    bool isDebugEnabled() const { return false; }
};
inline constexpr NullLogger logger{};
}


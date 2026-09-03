#pragma once

#include "ramag/model.hpp"

#include <filesystem>
#include <span>
#include <string>
#include <string_view>

namespace ramag {

// Strong, boundary-aware identity for the ordered normalized reference
// collection (name, complete header, and A/C/G/T/N bases).
[[nodiscard]] std::string NormalizedReferenceSha256(
    std::span<const SequenceRecord> records);
[[nodiscard]] std::string Sha256Hex(std::string_view content);
[[nodiscard]] std::string FileSha256(const std::filesystem::path& path);

}  // namespace ramag

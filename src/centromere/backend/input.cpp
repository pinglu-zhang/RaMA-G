#include "ramag/centromere.hpp"

#include "ramag/alignment.hpp"

#include <algorithm>
#include <charconv>
#include <fstream>
#include <numeric>
#include <sstream>
#include <string_view>
#include <unordered_map>
#include <utility>

namespace ramag {
namespace {

[[noreturn]] void BedError(const std::filesystem::path& path,
                          std::size_t line, std::string_view message) {
  throw CentromereInputError(path.string() + ":" + std::to_string(line) +
                            ": " + std::string(message));
}

Position ParseCoordinate(std::string_view text, const std::filesystem::path& path,
                         std::size_t line) {
  Position value{};
  if (text.empty() || text.find_first_not_of("0123456789") != std::string_view::npos) {
    BedError(path, line, "BED coordinates must be unsigned decimal integers");
  }
  const auto result = std::from_chars(text.data(), text.data() + text.size(), value);
  if (result.ec != std::errc{} || result.ptr != text.data() + text.size()) {
    BedError(path, line, "BED coordinate exceeds the 64-bit range");
  }
  return value;
}

}  // namespace

CentromereMode CentromereSpec::Mode() const noexcept {
  if (whole_input) return CentromereMode::WholeInput;
  if (!reference_bed.empty() || !query_bed.empty()) return CentromereMode::AnnotatedGenome;
  return CentromereMode::Genome;
}

std::string CentromereModeName(CentromereMode mode) {
  switch (mode) {
    case CentromereMode::Genome: return "genome";
    case CentromereMode::WholeInput: return "centromere-only";
    case CentromereMode::AnnotatedGenome: return "genome-with-centromeres";
  }
  throw CentromereInputError("invalid centromere mode");
}

std::vector<CentromereRegion> ReadCentromereBed(const std::filesystem::path& path) {
  std::error_code error;
  if (!std::filesystem::is_regular_file(path, error) || error) {
    throw CentromereInputError("centromere BED is not a regular file: " + path.string());
  }
  std::ifstream input(path);
  if (!input) throw CentromereInputError("cannot read centromere BED: " + path.string());
  std::vector<CentromereRegion> regions;
  std::string line;
  std::size_t line_number = 0;
  while (std::getline(input, line)) {
    ++line_number;
    if (line_number == 1 && line.starts_with("\xEF\xBB\xBF")) line.erase(0, 3);
    const auto begin = line.find_first_not_of(" \t\r");
    if (begin == std::string::npos || line[begin] == '#') continue;
    std::istringstream fields(line.substr(begin));
    std::string contig, start, end, name;
    fields >> contig;
    if (!(fields >> start)) BedError(path, line_number, "BED requires contig, start and end");
    if (contig == "track" && start.find('=') != std::string::npos) continue;
    if (contig == "browser" && (start == "position" || start == "hide" ||
        start == "dense" || start == "pack" || start == "squish" || start == "full")) continue;
    if (!(fields >> end)) BedError(path, line_number, "BED requires contig, start and end");
    fields >> name;
    const auto first = ParseCoordinate(start, path, line_number);
    const auto last = ParseCoordinate(end, path, line_number);
    if (first >= last) BedError(path, line_number, "BED requires 0 <= start < end (half-open)");
    regions.push_back({std::move(contig), first, last, std::move(name), line_number});
  }
  if (input.bad()) throw CentromereInputError("error reading centromere BED: " + path.string());
  if (regions.empty()) throw CentromereInputError("centromere BED contains no regions: " + path.string());

  // Sort IDs rather than annotations: preserve input order and reject double dispatch.
  std::vector<std::size_t> order(regions.size());
  std::iota(order.begin(), order.end(), std::size_t{0});
  std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
    const auto& left = regions[a];
    const auto& right = regions[b];
    if (left.contig != right.contig) return left.contig < right.contig;
    if (left.begin != right.begin) return left.begin < right.begin;
    if (left.end != right.end) return left.end < right.end;
    return a < b;
  });
  for (std::size_t i = 1; i < order.size(); ++i) {
    const auto& previous = regions[order[i - 1]];
    const auto& current = regions[order[i]];
    if (previous.contig == current.contig && current.begin < previous.end) {
      BedError(path, current.source_line, "overlapping or duplicate centromere regions; provide non-overlapping intervals");
    }
  }
  return regions;
}

CentromereAnnotations LoadCentromereAnnotations(const CentromereSpec& spec) {
  if (spec.whole_input && (!spec.reference_bed.empty() || !spec.query_bed.empty())) {
    throw CentromereInputError("--centromere cannot be combined with --ref-cen or --query-cen");
  }
  if (spec.reference_bed.empty() != spec.query_bed.empty()) {
    throw CentromereInputError("--ref-cen and --query-cen must be supplied together");
  }
  if (spec.Mode() != CentromereMode::AnnotatedGenome) return {};
  return {ReadCentromereBed(spec.reference_bed), ReadCentromereBed(spec.query_bed)};
}

void ValidateCentromereBounds(std::span<const CentromereRegion> regions,
                             std::span<const SequenceRecord> sequences) {
  std::unordered_map<std::string_view, Length> lengths;
  for (const auto& sequence : sequences) {
    if (!lengths.emplace(sequence.name, sequence.size()).second) {
      throw CentromereInputError("duplicate FASTA ID: " + sequence.name);
    }
  }
  for (const auto& region : regions) {
    const auto found = lengths.find(region.contig);
    if (found == lengths.end()) {
      throw CentromereInputError("BED line " + std::to_string(region.source_line) +
                                ": unknown FASTA ID " + region.contig);
    }
    if (region.begin >= region.end || region.end > found->second) {
      throw CentromereInputError("BED line " + std::to_string(region.source_line) +
                                ": interval outside FASTA record " + region.contig);
    }
  }
}

void RequireCentromereBackend(const CentromereSpec& spec) {
  if (spec.Mode() == CentromereMode::WholeInput && CentromereBackendAvailable()) return;
  if (spec.Mode() != CentromereMode::Genome) {
    throw DependencyError(
        "centromere alignment backend is not integrated yet; use --print-effective-config "
        "to check the centromere interface. No genome alignment fallback was run");
  }
}

}  // namespace ramag

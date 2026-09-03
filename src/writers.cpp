#include "ramag/writers.hpp"

#include "ramag/runtime.hpp"

#include "ramag/alignment.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>
#include <tuple>

namespace ramag {
namespace {

const SequenceRecord& FindSequence(const FastaData& data, SequenceId id) {
  const auto found = std::find_if(
      data.sequences.begin(), data.sequences.end(),
      [id](const SequenceRecord& sequence) { return sequence.numeric_id == id; });
  if (found == data.sequences.end()) {
    throw WriterError("alignment references unknown sequence id " +
                      std::to_string(id));
  }
  return *found;
}

std::string OrientedQuerySegment(const AlignmentRecord& alignment,
                                 const SequenceRecord& query) {
  if (alignment.query_begin > alignment.query_end ||
      alignment.query_end > query.bases.size()) {
    throw WriterError("query interval is outside contig " + query.name);
  }
  const auto begin = static_cast<std::size_t>(alignment.query_begin);
  const auto length = static_cast<std::size_t>(alignment.query_end -
                                               alignment.query_begin);
  std::string segment = query.bases.substr(begin, length);
  if (alignment.strand == Strand::Reverse) {
    segment = ReverseComplement(segment);
  }
  return segment;
}

void ValidateCoordinates(const AlignmentRecord& alignment,
                         const SequenceRecord& reference,
                         const SequenceRecord& query) {
  for (const auto& operation : alignment.cigar) {
    if (operation.length == 0 ||
        (operation.operation != '=' && operation.operation != 'X' &&
         operation.operation != 'I' && operation.operation != 'D')) {
      throw WriterError(
          "alignment CIGAR must contain only non-zero =/X/I/D operations for " +
          query.name);
    }
  }
  if (alignment.reference_begin >= alignment.reference_end ||
      alignment.reference_end > reference.bases.size()) {
    throw WriterError("invalid reference interval for " + reference.name);
  }
  if (alignment.query_begin >= alignment.query_end ||
      alignment.query_end > query.bases.size()) {
    throw WriterError("invalid query interval for " + query.name);
  }
  if (CigarReferenceLength(alignment.cigar) !=
      alignment.reference_end - alignment.reference_begin) {
    throw WriterError("CIGAR/reference span mismatch for " + query.name);
  }
  if (CigarQueryLength(alignment.cigar) !=
      alignment.query_end - alignment.query_begin) {
    throw WriterError("CIGAR/query span mismatch for " + query.name);
  }
  if (CigarEditDistance(alignment.cigar) != alignment.edit_distance) {
    throw WriterError("CIGAR/NM mismatch for " + query.name);
  }
}

Length AlignmentBlockLength(std::span<const CigarOp> cigar) {
  Length length = 0;
  for (const auto& operation : cigar) {
    switch (operation.operation) {
      case '=':
      case 'X':
      case 'M':
      case 'I':
      case 'D':
        length += operation.length;
        break;
      default:
        break;
    }
  }
  return length;
}

Length ExactMatchCount(const AlignmentRecord& alignment,
                       const SequenceRecord& reference,
                       const SequenceRecord& query) {
  const auto oriented_query = OrientedQuerySegment(alignment, query);
  std::size_t reference_offset =
      static_cast<std::size_t>(alignment.reference_begin);
  std::size_t query_offset = 0;
  Length matches = 0;
  for (const auto& operation : alignment.cigar) {
    for (Length index = 0; index < operation.length; ++index) {
      switch (operation.operation) {
        case '=':
          ++matches;
          ++reference_offset;
          ++query_offset;
          break;
        case 'M':
          if (reference.bases.at(reference_offset) ==
              oriented_query.at(query_offset)) {
            ++matches;
          }
          ++reference_offset;
          ++query_offset;
          break;
        case 'X':
          ++reference_offset;
          ++query_offset;
          break;
        case 'D':
          ++reference_offset;
          break;
        case 'I':
          ++query_offset;
          break;
        default:
          throw WriterError("unsupported CIGAR operator: " +
                            std::string(1, operation.operation));
      }
    }
  }
  return matches;
}

Length CountNonAlphaColumns(const AlignmentRecord& alignment,
                            const SequenceRecord& reference,
                            const SequenceRecord& query) {
  const auto oriented_query = OrientedQuerySegment(alignment, query);
  std::size_t reference_offset =
      static_cast<std::size_t>(alignment.reference_begin);
  std::size_t query_offset = 0;
  Length count = 0;
  for (const auto& operation : alignment.cigar) {
    for (Length index = 0; index < operation.length; ++index) {
      if (operation.operation == '=' || operation.operation == 'X' ||
          operation.operation == 'M') {
        if (!IsCanonicalBase(reference.bases.at(reference_offset)) ||
            !IsCanonicalBase(oriented_query.at(query_offset))) {
          ++count;
        }
        ++reference_offset;
        ++query_offset;
      } else if (operation.operation == 'D') {
        if (!IsCanonicalBase(reference.bases.at(reference_offset))) {
          ++count;
        }
        ++reference_offset;
      } else if (operation.operation == 'I') {
        if (!IsCanonicalBase(oriented_query.at(query_offset))) {
          ++count;
        }
        ++query_offset;
      }
    }
  }
  return count;
}

void RequireReadable(const std::filesystem::path& path,
                     std::ifstream& stream,
                     std::string_view format) {
  stream.open(path);
  if (!stream) {
    throw WriterError("cannot reopen " + std::string(format) +
                      " output for validation: " + path.string());
  }
}

std::size_t TabFieldCount(std::string_view line) {
  return static_cast<std::size_t>(
             std::count(line.begin(), line.end(), '\t')) +
         1U;
}

std::string EscapeSamHeaderValue(std::string_view value) {
  std::ostringstream escaped;
  escaped << std::hex << std::uppercase << std::setfill('0');
  for (const unsigned char character : value) {
    if (character == '\\') {
      escaped << "\\\\";
    } else if (character == '\t') {
      escaped << "\\t";
    } else if (character == '\n') {
      escaped << "\\n";
    } else if (character == '\r') {
      escaped << "\\r";
    } else if (character < 0x20U || character == 0x7fU) {
      escaped << "\\x" << std::setw(2) << static_cast<unsigned int>(character);
    } else {
      escaped << static_cast<char>(character);
    }
  }
  return escaped.str();
}

void ValidateSamQueryName(std::string_view name) {
  if (name.empty() || name.size() > 254 || name.front() == '@' || name == "*") {
    throw WriterError(
        "query FASTA id is not representable as a mapped SAM QNAME: " +
        std::string(name));
  }
  for (const unsigned char character : name) {
    // SAM 1.6 QNAME is [!-?A-~]{1,254}: printable ASCII with '@'
    // excluded everywhere, not just in the first position.
    if (character < 0x21U || character > 0x7eU || character == '@') {
      throw WriterError(
          "query FASTA id contains a byte not representable in SAM QNAME: " +
          std::string(name));
    }
  }
}

void ValidateSamReferenceName(std::string_view name) {
  // SAM 1.6 reference names use [!-()+-<>-~][!-~]*.  In particular, '*'
  // and '=' are reserved in the first position.  Accepting '*' as an @SQ SN
  // and mapped RNAME silently turns the record into an unmapped record in
  // common consumers, so fail before publishing any SAM artifact.
  if (name.empty() || name.front() == '*' || name.front() == '=') {
    throw WriterError(
        "reference FASTA id is not representable as a mapped SAM RNAME: " +
        std::string(name));
  }
  for (const unsigned char character : name) {
    if (character < 0x21U || character > 0x7eU) {
      throw WriterError(
          "reference FASTA id contains a byte not representable in SAM RNAME: " +
          std::string(name));
    }
  }
}

void ValidateAuxIntegerRanges(const AlignmentRecord& alignment,
                              std::string_view format) {
  constexpr std::int64_t kMinAuxInteger =
      static_cast<std::int64_t>(std::numeric_limits<std::int32_t>::min());
  constexpr std::uint64_t kMaxAuxInteger =
      static_cast<std::uint64_t>(std::numeric_limits<std::uint32_t>::max());
  if (alignment.score < kMinAuxInteger ||
      (alignment.score >= 0 &&
       static_cast<std::uint64_t>(alignment.score) > kMaxAuxInteger)) {
    throw WriterError(std::string(format) +
                      " AS:i value exceeds the SAM-compatible 32-bit range");
  }
  if (alignment.edit_distance > kMaxAuxInteger) {
    throw WriterError(std::string(format) +
                      " NM:i value exceeds the SAM-compatible 32-bit range");
  }
}

void AppendSamCigarOperation(std::ostringstream& output,
                             Length length,
                             char operation) {
  constexpr Length kMaxBamOperationLength = (Length{1} << 28U) - 1U;
  while (length != 0) {
    const Length chunk = std::min(length, kMaxBamOperationLength);
    output << chunk << operation;
    length -= chunk;
  }
}

}  // namespace

OutputPaths MakeOutputPaths(const std::filesystem::path& prefix) {
  const auto base = prefix.string();
  return {base + ".sam", base + ".paf", base + ".delta", base + ".maf",
          base + ".chain", base + ".manifest.json", base + ".complete"};
}

OutputPaths MakeOutputPaths(const RunSpec& spec) {
  auto paths = MakeOutputPaths(spec.output_prefix);
  for (const auto& request : spec.outputs) {
    switch (request.format) {
      case OutputFormat::Sam: paths.sam = request.path; break;
      case OutputFormat::Paf: paths.paf = request.path; break;
      case OutputFormat::Delta: paths.delta = request.path; break;
      case OutputFormat::Maf: paths.maf = request.path; break;
      case OutputFormat::Chain: paths.chain = request.path; break;
    }
  }
  return paths;
}

std::string SamCigar(const AlignmentRecord& alignment, Length query_length) {
  if (alignment.query_end > query_length) {
    throw WriterError("alignment exceeds query length while rendering SAM");
  }
  Length leading = alignment.query_begin;
  Length trailing = query_length - alignment.query_end;
  if (alignment.strand == Strand::Reverse) {
    std::swap(leading, trailing);
  }
  std::ostringstream cigar;
  if (leading != 0) {
    AppendSamCigarOperation(cigar, leading, 'H');
  }
  for (const auto& operation : alignment.cigar) {
    AppendSamCigarOperation(cigar, operation.length, operation.operation);
  }
  if (trailing != 0) {
    AppendSamCigarOperation(cigar, trailing, 'H');
  }
  return cigar.str();
}

std::string BuildMdTag(const AlignmentRecord& alignment,
                       const SequenceRecord& reference,
                       const SequenceRecord& query) {
  ValidateCoordinates(alignment, reference, query);
  const auto oriented_query = OrientedQuerySegment(alignment, query);
  std::size_t reference_offset =
      static_cast<std::size_t>(alignment.reference_begin);
  std::size_t query_offset = 0;
  Length match_run = 0;
  std::ostringstream md;
  const auto flush_matches = [&]() {
    md << match_run;
    match_run = 0;
  };

  for (const auto& operation : alignment.cigar) {
    if (operation.operation == 'I') {
      query_offset += static_cast<std::size_t>(operation.length);
      continue;
    }
    if (operation.operation == 'D') {
      flush_matches();
      md << '^';
      for (Length index = 0; index < operation.length; ++index) {
        md << reference.bases.at(reference_offset++);
      }
      continue;
    }
    if (operation.operation != '=' && operation.operation != 'X' &&
        operation.operation != 'M') {
      throw WriterError("unsupported CIGAR operator while building MD: " +
                        std::string(1, operation.operation));
    }
    for (Length index = 0; index < operation.length; ++index) {
      const char reference_base = reference.bases.at(reference_offset++);
      const char query_base = oriented_query.at(query_offset++);
      if (operation.operation == '=' ||
          (operation.operation == 'M' && reference_base == query_base)) {
        ++match_run;
      } else {
        flush_matches();
        md << reference_base;
      }
    }
  }
  flush_matches();
  return md.str();
}

void WriteSam(std::ostream& output,
              const FastaData& references,
              const FastaData& queries,
              std::span<const AlignmentRecord> alignments,
              std::string_view command_line) {
  output << "@HD\tVN:1.6\tSO:unknown\n";
  for (const auto& reference : references.sequences) {
    ValidateSamReferenceName(reference.name);
    output << "@SQ\tSN:" << reference.name << "\tLN:" << reference.size()
           << '\n';
  }
  output << "@PG\tID:ramag\tPN:ramag\tVN:" << RAMAG_VERSION
         << "\tCL:" << EscapeSamHeaderValue(command_line) << '\n';
  for (const auto& alignment : alignments) {
    CheckInterruption("writer-sam");
    const auto& reference = FindSequence(references, alignment.reference_id);
    const auto& query = FindSequence(queries, alignment.query_id);
    ValidateSamQueryName(query.name);
    ValidateCoordinates(alignment, reference, query);
    ValidateAuxIntegerRanges(alignment, "SAM");
    unsigned int flag = alignment.strand == Strand::Reverse ? 0x10U : 0U;
    if (!alignment.primary) {
      flag |= 0x800U;
    }
    output << query.name << '\t' << flag << '\t' << reference.name << '\t'
           << alignment.reference_begin + 1 << "\t255\t"
           << SamCigar(alignment, query.size())
           << "\t*\t0\t0\t*\t*\tNM:i:" << alignment.edit_distance
           << "\tMD:Z:" << BuildMdTag(alignment, reference, query)
           << "\tAS:i:" << alignment.score << '\n';
  }
  if (!output) {
    throw WriterError("failed while writing SAM output");
  }
}

void WritePaf(std::ostream& output,
              const FastaData& references,
              const FastaData& queries,
              std::span<const AlignmentRecord> alignments) {
  for (const auto& alignment : alignments) {
    CheckInterruption("writer-paf");
    const auto& reference = FindSequence(references, alignment.reference_id);
    const auto& query = FindSequence(queries, alignment.query_id);
    ValidateCoordinates(alignment, reference, query);
    ValidateAuxIntegerRanges(alignment, "PAF");
    const auto matches = ExactMatchCount(alignment, reference, query);
    output << query.name << '\t' << query.size() << '\t'
           << alignment.query_begin << '\t' << alignment.query_end << '\t'
           << (alignment.strand == Strand::Forward ? '+' : '-') << '\t'
           << reference.name << '\t' << reference.size() << '\t'
           << alignment.reference_begin << '\t' << alignment.reference_end
           << '\t' << matches << '\t' << AlignmentBlockLength(alignment.cigar)
           << "\t255\tcg:Z:" << CigarToString(alignment.cigar)
           << "\tNM:i:" << alignment.edit_distance << "\tAS:i:"
           << alignment.score << "\ttp:A:"
           << (alignment.primary ? 'P' : 'S') << '\n';
  }
  if (!output) {
    throw WriterError("failed while writing PAF output");
  }
}

std::vector<std::int64_t> CigarToDelta(std::span<const CigarOp> cigar) {
  std::vector<std::int64_t> deltas;
  std::int64_t distance = 1;
  for (const auto& operation : cigar) {
    switch (operation.operation) {
      case '=':
      case 'X':
      case 'M':
        if (operation.length > static_cast<Length>(
                                   std::numeric_limits<std::int64_t>::max() -
                                   distance)) {
          throw WriterError("delta distance overflows int64");
        }
        distance += static_cast<std::int64_t>(operation.length);
        break;
      case 'D':
      case 'I':
        for (Length index = 0; index < operation.length; ++index) {
          deltas.push_back(operation.operation == 'D' ? distance : -distance);
          distance = 1;
        }
        break;
      default:
        throw WriterError("unsupported CIGAR operator for delta: " +
                          std::string(1, operation.operation));
    }
  }
  return deltas;
}

void WriteDelta(std::ostream& output,
                const FastaData& references,
                const FastaData& queries,
                std::span<const AlignmentRecord> alignments) {
  std::error_code error;
  auto reference_path = std::filesystem::absolute(references.source, error);
  if (error) {
    reference_path = references.source;
  }
  error.clear();
  auto query_path = std::filesystem::absolute(queries.source, error);
  if (error) {
    query_path = queries.source;
  }
  output << reference_path.string() << ' ' << query_path.string()
         << "\nNUCMER\n";
  bool have_pair = false;
  SequenceId previous_reference_id = 0;
  SequenceId previous_query_id = 0;
  for (const auto& alignment : alignments) {
    CheckInterruption("writer-delta");
    const auto& reference = FindSequence(references, alignment.reference_id);
    const auto& query = FindSequence(queries, alignment.query_id);
    ValidateCoordinates(alignment, reference, query);
    const Position reference_start = alignment.reference_begin + 1;
    const Position reference_end = alignment.reference_end;
    const Position query_start = alignment.strand == Strand::Forward
                                     ? alignment.query_begin + 1
                                     : alignment.query_end;
    const Position query_end = alignment.strand == Strand::Forward
                                   ? alignment.query_end
                                   : alignment.query_begin + 1;
    if (!have_pair || alignment.reference_id != previous_reference_id ||
        alignment.query_id != previous_query_id) {
      output << '>' << reference.name << ' ' << query.name << ' '
             << reference.size() << ' ' << query.size() << '\n';
      previous_reference_id = alignment.reference_id;
      previous_query_id = alignment.query_id;
      have_pair = true;
    }
    output << reference_start << ' ' << reference_end << ' ' << query_start
           << ' ' << query_end << ' ' << alignment.edit_distance << ' '
           << alignment.edit_distance << ' '
           << CountNonAlphaColumns(alignment, reference, query) << '\n';
    for (const auto delta : CigarToDelta(alignment.cigar)) {
      output << delta << '\n';
    }
    output << "0\n";
  }
  if (!output) {
    throw WriterError("failed while writing delta output");
  }
}

void WriteMaf(std::ostream& output,
              const FastaData& references,
              const FastaData& queries,
              std::span<const AlignmentRecord> alignments) {
  output << "##maf version=1 scoring=RaMA-G\n\n";
  for (const auto& alignment : alignments) {
    CheckInterruption("writer-maf");
    const auto& reference = FindSequence(references, alignment.reference_id);
    const auto& query = FindSequence(queries, alignment.query_id);
    ValidateCoordinates(alignment, reference, query);
    const auto oriented_query = OrientedQuerySegment(alignment, query);
    std::size_t reference_offset =
        static_cast<std::size_t>(alignment.reference_begin);
    std::size_t query_offset = 0;
    std::string reference_text;
    std::string query_text;
    const auto columns = AlignmentBlockLength(alignment.cigar);
    if (columns > static_cast<Length>(std::numeric_limits<std::size_t>::max())) {
      throw WriterError("MAF block exceeds this process address range");
    }
    reference_text.reserve(static_cast<std::size_t>(columns));
    query_text.reserve(static_cast<std::size_t>(columns));
    for (const auto& operation : alignment.cigar) {
      for (Length index = 0; index < operation.length; ++index) {
        if (operation.operation == '=' || operation.operation == 'X') {
          reference_text.push_back(reference.bases.at(reference_offset++));
          query_text.push_back(oriented_query.at(query_offset++));
        } else if (operation.operation == 'D') {
          reference_text.push_back(reference.bases.at(reference_offset++));
          query_text.push_back('-');
        } else if (operation.operation == 'I') {
          reference_text.push_back('-');
          query_text.push_back(oriented_query.at(query_offset++));
        } else {
          throw WriterError("unsupported CIGAR operator for MAF: " +
                            std::string(1, operation.operation));
        }
      }
    }
    const Position query_start = alignment.strand == Strand::Forward
                                     ? alignment.query_begin
                                     : query.size() - alignment.query_end;
    output << "a score=" << alignment.score << '\n'
           << "s " << reference.name << ' ' << alignment.reference_begin << ' '
           << (alignment.reference_end - alignment.reference_begin) << " + "
           << reference.size() << ' ' << reference_text << '\n'
           << "s " << query.name << ' ' << query_start << ' '
           << (alignment.query_end - alignment.query_begin) << ' '
           << (alignment.strand == Strand::Forward ? '+' : '-') << ' '
           << query.size() << ' ' << query_text << "\n\n";
  }
  if (!output) throw WriterError("failed while writing MAF output");
}

namespace {

struct ChainBlock {
  Length size{};
  Length target_gap{};
  Length query_gap{};
};

struct ChainProjection {
  Position target_begin{};
  Position target_end{};
  Position query_begin{};
  Position query_end{};
  std::vector<ChainBlock> blocks;
};

ChainProjection ProjectChain(const AlignmentRecord& alignment,
                             const SequenceRecord& query) {
  ChainProjection projection;
  Position target = alignment.reference_begin;
  Position oriented_query = alignment.strand == Strand::Forward
                                ? alignment.query_begin
                                : query.size() - alignment.query_end;
  bool have_block = false;
  Length pending_target_gap = 0;
  Length pending_query_gap = 0;
  for (const auto& operation : alignment.cigar) {
    if (operation.operation == 'D') {
      target += operation.length;
      if (have_block) pending_target_gap += operation.length;
      continue;
    }
    if (operation.operation == 'I') {
      oriented_query += operation.length;
      if (have_block) pending_query_gap += operation.length;
      continue;
    }
    if (operation.operation != '=' && operation.operation != 'X') {
      throw WriterError("unsupported CIGAR operator for chain: " +
                        std::string(1, operation.operation));
    }
    if (!have_block) {
      projection.target_begin = target;
      projection.query_begin = oriented_query;
      projection.blocks.push_back(ChainBlock{operation.length, 0, 0});
      have_block = true;
    } else if (pending_target_gap == 0 && pending_query_gap == 0) {
      projection.blocks.back().size += operation.length;
    } else {
      projection.blocks.back().target_gap = pending_target_gap;
      projection.blocks.back().query_gap = pending_query_gap;
      projection.blocks.push_back(ChainBlock{operation.length, 0, 0});
      pending_target_gap = 0;
      pending_query_gap = 0;
    }
    target += operation.length;
    oriented_query += operation.length;
    projection.target_end = target;
    projection.query_end = oriented_query;
  }
  if (!have_block) {
    throw WriterError("alignment has no aligned block representable in chain");
  }
  return projection;
}

}  // namespace

void WriteChain(std::ostream& output,
                const FastaData& references,
                const FastaData& queries,
                std::span<const AlignmentRecord> alignments) {
  std::uint64_t chain_id = 1;
  for (const auto& alignment : alignments) {
    CheckInterruption("writer-chain");
    const auto& reference = FindSequence(references, alignment.reference_id);
    const auto& query = FindSequence(queries, alignment.query_id);
    ValidateCoordinates(alignment, reference, query);
    const auto projection = ProjectChain(alignment, query);
    output << "chain " << std::max<std::int64_t>(0, alignment.score) << ' '
           << reference.name << ' ' << reference.size() << " + "
           << projection.target_begin << ' ' << projection.target_end << ' '
           << query.name << ' ' << query.size() << ' '
           << (alignment.strand == Strand::Forward ? '+' : '-') << ' '
           << projection.query_begin << ' ' << projection.query_end << ' '
           << chain_id++ << '\n';
    for (std::size_t index = 0; index < projection.blocks.size(); ++index) {
      const auto& block = projection.blocks[index];
      output << block.size;
      if (index + 1 != projection.blocks.size()) {
        output << '\t' << block.target_gap << '\t' << block.query_gap;
      }
      output << '\n';
    }
    output << '\n';
  }
  if (!output) throw WriterError("failed while writing chain output");
}

void ValidateSamFile(const std::filesystem::path& path) {
  std::ifstream input;
  RequireReadable(path, input, "SAM");
  std::string line;
  bool saw_header = false;
  while (std::getline(input, line)) {
    if (line.rfind("@HD\tVN:1.6", 0) == 0) {
      saw_header = true;
    } else if (!line.empty() && line.front() == '@') {
      if (line.rfind("@SQ\t", 0) != 0 && line.rfind("@RG\t", 0) != 0 &&
          line.rfind("@PG\t", 0) != 0 && line.rfind("@CO\t", 0) != 0) {
        throw WriterError("SAM contains an unknown or misplaced header line: " +
                          path.string());
      }
    } else if (!line.empty() && TabFieldCount(line) < 11) {
      throw WriterError("SAM record has fewer than 11 fields: " + path.string());
    }
  }
  if (!saw_header || input.bad()) {
    throw WriterError("SAM validation failed: " + path.string());
  }
}

void ValidatePafFile(const std::filesystem::path& path) {
  std::ifstream input;
  RequireReadable(path, input, "PAF");
  std::string line;
  while (std::getline(input, line)) {
    if (!line.empty() && TabFieldCount(line) < 12) {
      throw WriterError("PAF record has fewer than 12 fields: " + path.string());
    }
  }
  if (input.bad()) {
    throw WriterError("PAF validation failed: " + path.string());
  }
}

void ValidateDeltaFile(const std::filesystem::path& path) {
  std::ifstream input;
  RequireReadable(path, input, "delta");
  std::string first;
  std::string second;
  if (!std::getline(input, first) || !std::getline(input, second) ||
      first.empty() || second != "NUCMER") {
    throw WriterError("delta header validation failed: " + path.string());
  }
}

void ValidateMafFile(const std::filesystem::path& path) {
  std::ifstream input;
  RequireReadable(path, input, "MAF");
  std::string line;
  if (!std::getline(input, line) || line.rfind("##maf version=1", 0) != 0) {
    throw WriterError("MAF header validation failed: " + path.string());
  }
  std::size_t sequence_rows = 0;
  std::size_t first_text_size = 0;
  while (std::getline(input, line)) {
    if (line.empty()) {
      if (sequence_rows != 0 && sequence_rows != 2) {
        throw WriterError("MAF block does not contain exactly two s rows: " +
                          path.string());
      }
      sequence_rows = 0;
      first_text_size = 0;
      continue;
    }
    if (line.rfind("a score=", 0) == 0) {
      if (sequence_rows != 0) {
        throw WriterError("MAF block boundary is malformed: " + path.string());
      }
      continue;
    }
    if (line.rfind("s ", 0) != 0) {
      throw WriterError("MAF contains an unsupported line: " + path.string());
    }
    std::istringstream fields(line);
    char row = 0;
    std::string source;
    std::uint64_t start = 0;
    std::uint64_t size = 0;
    char strand = 0;
    std::uint64_t source_size = 0;
    std::string text;
    fields >> row >> source >> start >> size >> strand >> source_size >> text;
    std::string trailing;
    if (!fields || fields >> trailing || row != 's' || source.empty() ||
        (strand != '+' && strand != '-') || start > source_size ||
        size > source_size - start || text.empty()) {
      throw WriterError("MAF s row validation failed: " + path.string());
    }
    const auto nongap = static_cast<std::uint64_t>(
        std::count_if(text.begin(), text.end(), [](char value) { return value != '-'; }));
    if (nongap != size) {
      throw WriterError("MAF s row size disagrees with text: " + path.string());
    }
    if (sequence_rows == 0) first_text_size = text.size();
    else if (text.size() != first_text_size) {
      throw WriterError("MAF pair rows have different column counts: " + path.string());
    }
    ++sequence_rows;
  }
  if (input.bad() || (sequence_rows != 0 && sequence_rows != 2)) {
    throw WriterError("MAF validation failed: " + path.string());
  }
}

void ValidateChainFile(const std::filesystem::path& path) {
  std::ifstream input;
  RequireReadable(path, input, "chain");
  std::string line;
  bool inside_chain = false;
  std::uint64_t target_remaining = 0;
  std::uint64_t query_remaining = 0;
  while (std::getline(input, line)) {
    if (line.empty()) {
      if (inside_chain || target_remaining != 0 || query_remaining != 0) {
        throw WriterError("chain block span does not close: " + path.string());
      }
      continue;
    }
    if (line.rfind("chain ", 0) == 0) {
      if (inside_chain) throw WriterError("nested chain header: " + path.string());
      std::istringstream fields(line);
      std::string keyword;
      std::uint64_t score = 0;
      std::string target_name;
      std::uint64_t target_size = 0;
      char target_strand = 0;
      std::uint64_t target_begin = 0;
      std::uint64_t target_end = 0;
      std::string query_name;
      std::uint64_t query_size = 0;
      char query_strand = 0;
      std::uint64_t query_begin = 0;
      std::uint64_t query_end = 0;
      std::uint64_t id = 0;
      fields >> keyword >> score >> target_name >> target_size >> target_strand >>
          target_begin >> target_end >> query_name >> query_size >> query_strand >>
          query_begin >> query_end >> id;
      std::string trailing;
      if (!fields || fields >> trailing || keyword != "chain" || id == 0 ||
          target_name.empty() || query_name.empty() || target_strand != '+' ||
          (query_strand != '+' && query_strand != '-') ||
          target_begin >= target_end || target_end > target_size ||
          query_begin >= query_end || query_end > query_size) {
        throw WriterError("chain header validation failed: " + path.string());
      }
      target_remaining = target_end - target_begin;
      query_remaining = query_end - query_begin;
      inside_chain = true;
      continue;
    }
    if (!inside_chain) throw WriterError("chain block without header: " + path.string());
    std::istringstream fields(line);
    std::uint64_t block = 0;
    std::uint64_t target_gap = 0;
    std::uint64_t query_gap = 0;
    fields >> block;
    if (!fields || block == 0 || block > target_remaining || block > query_remaining) {
      throw WriterError("invalid chain block size: " + path.string());
    }
    target_remaining -= block;
    query_remaining -= block;
    if (fields >> target_gap) {
      if (!(fields >> query_gap) || target_gap > target_remaining ||
          query_gap > query_remaining) {
        throw WriterError("invalid chain gap: " + path.string());
      }
      std::string trailing;
      if (fields >> trailing) throw WriterError("extra chain fields: " + path.string());
      target_remaining -= target_gap;
      query_remaining -= query_gap;
    } else {
      if (target_remaining != 0 || query_remaining != 0) {
        throw WriterError("final chain block does not close span: " + path.string());
      }
      inside_chain = false;
    }
  }
  if (input.bad() || inside_chain || target_remaining != 0 || query_remaining != 0) {
    throw WriterError("chain validation failed: " + path.string());
  }
}

}  // namespace ramag

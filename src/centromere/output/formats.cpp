#include "common/centromere.hpp"
#include "ramag/centromere_stitch_formats.hpp"

#include <algorithm>
#include <iterator>
#include <limits>
#include <map>
#include <ostream>
#include <set>
#include <string>
#include <vector>

namespace ramag {
namespace detail {
// Private helper is separately testable without allocating a 256 MiB sequence.
void WriteStitchSamCigarRun(std::ostream& out, Length length, char operation) {
  constexpr Length max_run = (Length{1} << 28U) - 1U;
  while (length) {
    const auto chunk = std::min(length, max_run);
    out << chunk << operation;
    length -= chunk;
  }
}
}  // namespace detail
namespace {
struct Counts {
  std::uint64_t matches{}, paired{}, columns{}, edits{};
};
struct ExportSegment {
  const CentromereStitchSegment* segment{};
  Counts counts;
};

using ::ramag::detail::IsCanonicalUpper;
char Complement(char c) {
  switch (c) {
    case 'A': return 'T'; case 'C': return 'G'; case 'G': return 'C'; case 'T': return 'A';
    default: return 'N';
  }
}
char Sign(Strand strand) { return strand == Strand::Forward ? '+' : '-'; }
void CheckName(const std::string& name) {
  if (name.empty() || name == "*") throw AlignmentError("stitch export requires a nonempty sequence name other than '*'");
  for (unsigned char c : name)
    if (c < 33 || c > 126) throw AlignmentError("stitch export sequence names must be printable ASCII without whitespace");
}
void CheckSamNames(const SequenceRecord& reference, const SequenceRecord& query) {
  // SAM QNAME excludes '@'; reference punctuation is restricted by SAM v1.6.
  if (query.name.size() > 254 || query.name.find('@') != std::string::npos)
    throw AlignmentError("stitch SAM query name exceeds 254 characters or contains '@'");
  if (reference.name.front() == '*' || reference.name.front() == '=')
    throw AlignmentError("stitch SAM reference name must not start with '*' or '='");
  constexpr std::string_view punctuation = "!#$%&*+./:;=?@^_|~-";
  for (unsigned char c : reference.name)
    if (!((c >= '0' && c <= '9') || (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z')) &&
        punctuation.find(static_cast<char>(c)) == std::string_view::npos)
      throw AlignmentError("stitch SAM reference name contains reserved punctuation");
  if (reference.size() == 0 || reference.size() > static_cast<Length>(std::numeric_limits<std::int32_t>::max()))
    throw AlignmentError("stitch SAM reference length is outside the SAM @SQ range");
}
void Add(std::uint64_t& destination, std::uint64_t n) {
  if (n > std::numeric_limits<std::uint64_t>::max() - destination)
    throw AlignmentError("stitch export CIGAR count overflow");
  destination += n;
}
std::vector<ExportSegment> Validate(const SequenceRecord& reference, const SequenceRecord& query,
                                  const CentromereStitchResult& result) {
  if (result.aborted) throw AlignmentError("cannot export an aborted stitch snapshot as alignment results");
  CheckName(reference.name); CheckName(query.name);
  for (const auto* sequence : {&reference, &query})
    for (char c : sequence->bases)
      if (!IsCanonicalUpper(c) && c != 'N') throw AlignmentError("stitch export requires normalized uppercase ACGTN sequence");
  if (result.reference_length != reference.size() || result.query_length != query.size())
    throw AlignmentError("stitch export sequence lengths disagree with stitch result");
  if (result.score_model != CentromereStitchScoreModel::Centrolign3Piece &&
      result.score_model != CentromereStitchScoreModel::Rama2PieceControl)
    throw AlignmentError("invalid stitch export score model");
  std::set<std::uint64_t> ids;
  std::vector<ExportSegment> entries;
  entries.reserve(result.segments.size());
  for (const auto& segment : result.segments) {
    const auto& a = segment.alignment;
    if (!segment.segment_id || !segment.path_id || !ids.insert(segment.segment_id).second ||
        (segment.path_role != "primary" && segment.path_role != "copy"))
      throw AlignmentError("invalid stitch export segment identity or role");
    if (a.reference_id != reference.numeric_id || a.query_id != query.numeric_id ||
        a.reference_begin >= a.reference_end || a.reference_end > reference.size() ||
        a.query_begin >= a.query_end || a.query_end > query.size() || a.cigar.empty() ||
        (a.strand != Strand::Forward && a.strand != Strand::Reverse))
      throw AlignmentError("invalid stitch export segment sequence identity, bounds or strand");
    Position r = a.reference_begin, q = 0;
    Counts counts;
    const Length qspan = a.query_end - a.query_begin;
    for (const auto& op : a.cigar) {
      if (!op.length || (op.operation != '=' && op.operation != 'X' && op.operation != 'I' && op.operation != 'D'))
        throw AlignmentError("stitch export requires nonempty =/X/I/D CIGAR runs");
      if (op.operation != 'I' && op.length > a.reference_end - r)
        throw AlignmentError("stitch export CIGAR overruns reference interval");
      if (op.operation != 'D' && op.length > qspan - q)
        throw AlignmentError("stitch export CIGAR overruns query interval");
      Add(counts.columns, op.length);
      if (op.operation == '=' || op.operation == 'X') {
        Add(counts.paired, op.length);
        for (Length k = 0; k < op.length; ++k) {
          const char rb = reference.bases[r + k];
          const char qb = a.strand == Strand::Forward ? query.bases[a.query_begin + q + k]
              : Complement(query.bases[a.query_end - q - k - 1]);
          if ((IsCanonicalUpper(rb) && rb == qb) != (op.operation == '='))
            throw AlignmentError("stitch export =/X disagrees with sequence characters");
        }
      }
      if (op.operation == '=') Add(counts.matches, op.length);
      else Add(counts.edits, op.length);
      if (op.operation != 'I') r += op.length;
      if (op.operation != 'D') q += op.length;
    }
    if (r != a.reference_end || q != qspan || counts.edits != a.edit_distance)
      throw AlignmentError("stitch export CIGAR consumption or edit distance disagrees with segment");
    entries.push_back({&segment, counts});
  }
  std::sort(entries.begin(), entries.end(), [](const auto& a, const auto& b) {
    return a.segment->segment_id < b.segment->segment_id;
  });
  return entries;
}
void Cigar(std::ostream& out, const std::vector<CigarOp>& cigar) {
  for (const auto& op : cigar) out << op.length << op.operation;
}
void IdList(std::ostream& out, const std::vector<std::uint64_t>& ids) {
  for (std::size_t i = 0; i < ids.size(); ++i) { if (i) out << ','; out << ids[i]; }
}
void Tags(std::ostream& out, const ExportSegment& entry, const CentromereStitchResult& result) {
  const auto& s = *entry.segment;
  out << "\tsi:Z:" << s.segment_id << "\tpi:Z:" << s.path_id << "\tpr:Z:" << s.path_role
      << "\tri:Z:" << s.region_id << "\tbi:Z:" << s.branch_id
      << "\tsm:Z:" << CentromereStitchScoreModelName(result.score_model)
      << "\tst:Z:" << (result.complete ? "complete" : "partial") << "\tsc:Z:" << s.alignment.score;
  if (entry.counts.edits <= static_cast<std::uint64_t>(std::numeric_limits<std::int32_t>::max()))
    out << "\tNM:i:" << entry.counts.edits;
  else out << "\tne:Z:" << entry.counts.edits;
  if (s.alignment.score >= std::numeric_limits<std::int32_t>::min() &&
      s.alignment.score <= std::numeric_limits<std::int32_t>::max()) out << "\tAS:i:" << s.alignment.score;
  if (!s.anchor_ids.empty()) { out << "\tai:Z:"; IdList(out, s.anchor_ids); }
  if (!s.gap_ids.empty()) { out << "\tgi:Z:"; IdList(out, s.gap_ids); }
}
void CheckOutput(std::ostream& out) {
  if (!out) throw AlignmentError("failed to write centromere alignment export");
}
void QueryFragment(std::ostream& out, const SequenceRecord& query, const AlignmentRecord& a) {
  if (a.strand == Strand::Forward)
    out.write(query.bases.data() + a.query_begin, static_cast<std::streamsize>(a.query_end - a.query_begin));
  else for (Position q = a.query_end; q > a.query_begin; ) out.put(Complement(query.bases[--q]));
}
// A map of disjoint query intervals permits O(log n) admission per segment.
bool Admit(std::map<Position, Position>& occupied, Position begin, Position end) {
  const auto next = occupied.lower_bound(begin);
  if (next != occupied.end() && next->first < end) return false;
  if (next != occupied.begin() && std::prev(next)->second > begin) return false;
  occupied.emplace(begin, end);
  return true;
}
void MafRow(std::ostream& out, bool reference_row, const SequenceRecord& reference,
            const SequenceRecord& query, const AlignmentRecord& a) {
  const auto& sequence = reference_row ? reference : query;
  const auto begin = reference_row ? a.reference_begin : a.query_begin;
  const auto end = reference_row ? a.reference_end : a.query_end;
  const bool reverse = !reference_row && a.strand == Strand::Reverse;
  out << "s " << (reference_row ? "ref." : "query.") << sequence.name << ' '
      << (reverse ? sequence.size() - end : begin) << ' ' << end - begin << ' '
      << (reverse ? '-' : '+') << ' ' << sequence.size() << ' ';
  Position r = a.reference_begin, q = 0;
  for (const auto& op : a.cigar) {
    const bool gap = reference_row ? op.operation == 'I' : op.operation == 'D';
    for (Length k = 0; k < op.length; ++k) {
      if (gap) out.put('-');
      else if (reference_row) out.put(reference.bases[r + k]);
      else out.put(reverse ? Complement(query.bases[a.query_end - q - k - 1])
                           : query.bases[a.query_begin + q + k]);
    }
    if (op.operation != 'I') r += op.length;
    if (op.operation != 'D') q += op.length;
  }
  out << '\n';
}
}  // namespace

void WriteCentromereStitchPaf(std::ostream& out, const SequenceRecord& reference,
                             const SequenceRecord& query, const CentromereStitchResult& result) {
  const auto entries = Validate(reference, query, result);
  for (const auto& entry : entries) {
    const auto& a = entry.segment->alignment;
    out << query.name << '\t' << query.size() << '\t' << a.query_begin << '\t' << a.query_end << '\t'
        << Sign(a.strand) << '\t' << reference.name << '\t' << reference.size() << '\t'
        << a.reference_begin << '\t' << a.reference_end << '\t' << entry.counts.matches << '\t'
        << entry.counts.columns << "\t255\tcg:Z:";
    Cigar(out, a.cigar); Tags(out, entry, result); out << '\n';
  }
  CheckOutput(out);
}

void WriteCentromereStitchSam(std::ostream& out, const SequenceRecord& reference,
                             const SequenceRecord& query, const CentromereStitchResult& result) {
  auto entries = Validate(reference, query, result);
  CheckSamNames(reference, query);
  std::sort(entries.begin(), entries.end(), [](const auto& a, const auto& b) {
    if (a.segment->path_role != b.segment->path_role) return a.segment->path_role == "primary";
    return a.segment->segment_id < b.segment->segment_id;
  });
  if (!entries.empty()) {
    const auto representative = std::max_element(entries.begin(), entries.end(), [](const auto& a, const auto& b) {
      if (a.segment->path_role != b.segment->path_role) return b.segment->path_role == "primary";
      if (a.counts.paired != b.counts.paired) return a.counts.paired < b.counts.paired;
      return a.segment->segment_id > b.segment->segment_id;
    });
    std::rotate(entries.begin(), representative, std::next(representative));
  }
  out << "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:" << reference.name << "\tLN:" << reference.size()
      << "\n@PG\tID:ramag-centromere-stitch\tPN:RaMA-G\n"
      << "@CO\tSegments only; hard clipping retains full-query coordinates. MAPQ unavailable.\n"
      << "@CO\tpr=primary/copy is algorithm role, not SAM flag. Secondary records retain accepted extra relations; do not filter for full relation evaluation.\n"
      << "@CO\tNo SA tags: path/segment IDs link output without quadratic tag expansion. Status="
      << (result.complete ? "complete" : "partial") << "; complete scope=eligible_internal_gaps.\n";
  if (entries.empty()) {
    out << query.name << "\t4\t*\t0\t0\t*\t*\t0\t0\t" << (query.bases.empty() ? "*" : query.bases)
        << "\t*\tst:Z:" << (result.complete ? "complete" : "partial") << '\n';
  }
  std::map<Position, Position> occupied;
  for (std::size_t i = 0; i < entries.size(); ++i) {
    const auto& entry = entries[i];
    const auto& a = entry.segment->alignment;
    const bool canonical = Admit(occupied, a.query_begin, a.query_end);
    const unsigned flag = (a.strand == Strand::Reverse ? 16U : 0U) | (i ? (canonical ? 2048U : 256U) : 0U);
    const Length left = a.strand == Strand::Forward ? a.query_begin : query.size() - a.query_end;
    const Length right = a.strand == Strand::Forward ? query.size() - a.query_end : a.query_begin;
    out << query.name << '\t' << flag << '\t' << reference.name << '\t' << a.reference_begin + 1 << "\t255\t";
    detail::WriteStitchSamCigarRun(out, left, 'H');
    for (const auto& op : a.cigar) detail::WriteStitchSamCigarRun(out, op.length, op.operation);
    detail::WriteStitchSamCigarRun(out, right, 'H');
    out << "\t*\t0\t0\t"; QueryFragment(out, query, a); out << "\t*";
    Tags(out, entry, result); out << '\n';
  }
  CheckOutput(out);
}

void WriteCentromereStitchMaf(std::ostream& out, const SequenceRecord& reference,
                             const SequenceRecord& query, const CentromereStitchResult& result) {
  const auto entries = Validate(reference, query, result);
  out << "##maf version=1 program=RaMA-G scoring=" << CentromereStitchScoreModelName(result.score_model)
      << "\n# ref. and query. prefix original FASTA names to distinguish sequence roles.\n"
      << "# status=" << (result.complete ? "complete" : "partial")
      << " complete_scope=eligible_internal_gaps; independent blocks may share coordinates.\n\n";
  for (const auto& entry : entries) {
    const auto& s = *entry.segment;
    out << "# ramag si=" << s.segment_id << " pi=" << s.path_id << " pr=" << s.path_role
        << " ri=" << s.region_id << " bi=" << s.branch_id
        << " sm=" << CentromereStitchScoreModelName(result.score_model)
        << " st=" << (result.complete ? "complete" : "partial") << " sc=" << s.alignment.score;
    if (!s.anchor_ids.empty()) { out << " ai="; IdList(out, s.anchor_ids); }
    if (!s.gap_ids.empty()) { out << " gi="; IdList(out, s.gap_ids); }
    out << "\na score=" << s.alignment.score << '\n';
    MafRow(out, true, reference, query, s.alignment);
    MafRow(out, false, reference, query, s.alignment);
    out << '\n';
  }
  CheckOutput(out);
}
}  // namespace ramag

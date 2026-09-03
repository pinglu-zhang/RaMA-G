#pragma once

#include "ramag/fasta.hpp"

#include <cstdint>
#include <filesystem>
#include <string>
#include <string_view>

namespace ramag {

// The adapter is compiled against this exact public SeqPro API.  CMake must
// reject a source override unless its clean HEAD equals this value.
inline constexpr std::string_view kRequiredSeqProCommit =
    "6781cadcf81a0da53d7573444594c1484947017c";

struct SeqProInputOptions {
    // Required.  The active FAI and any SeqPro metadata sidecar are kept below
    // this directory.  A standard <input>.fai, when present, is copied here
    // before validation; the adapter never writes into the input directory.
    std::filesystem::path work_dir;

    // Required role-specific file name, for example "reference" or "query".
    // It must be a basename rather than a path.  ".fai" is appended when it is
    // absent, preventing reference/query indexes in one work directory from
    // colliding when callers use distinct role names.
    std::string index_basename;

    SequenceId first_id{0};

    // False is the safe default: a stale or malformed existing role index is
    // diagnosed instead of silently replaced.  A caller may explicitly opt
    // into SeqPro's atomic rebuild policy for an isolated run directory.
    bool force_rebuild{false};
};

struct SeqProInputResult {
    FastaData fasta;
    std::filesystem::path fasta_index_path;
    std::filesystem::path metadata_path;
    // Populated only when a standard <input>.fai was copied into the isolated
    // run work directory.  The SHA-256 binds manifest provenance to the exact
    // bytes that were copied and validated.
    std::filesystem::path source_fasta_index_path;
    std::uint64_t source_fasta_index_size_bytes{0};
    std::string source_fasta_index_sha256;
    bool source_fasta_index_copied{false};
    bool external_fasta_index_adopted{false};
    std::string build_action;
    std::string index_origin;
    std::string verification_status;
};

// True only when this translation unit was compiled and linked with the
// pinned SeqPro-compatible API.
[[nodiscard]] bool SeqProAdapterAvailable() noexcept;

// Copy and structurally validate a standard <input>.fai when one is present;
// otherwise build an explicitly work-directory-scoped FAI plus SeqPro
// metadata.  The original FASTA is then memory-mapped and materialized into
// normalized RaMA-G records.  FAI order becomes numeric-id order.  The FAI
// name is the first FASTA-header token; the complete header is unavailable
// through the pinned API and therefore degrades to that name.
//
// A/C/G/T are upper-cased.  Every other ASCII alphabetic symbol becomes N and
// is counted as an ambiguous seed hard break.  Non-alphabetic bytes, empty
// records, coordinate overflows, and SeqPro failures are reported as
// FastaError (UnsupportedFastaFormat for SeqPro's unsupported-format code).
[[nodiscard]] SeqProInputResult ReadFastaWithSeqPro(
    const std::filesystem::path& fasta_path,
    const SeqProInputOptions& options);

}  // namespace ramag

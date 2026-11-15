# Troubleshooting

## Build or dependency configuration fails

Check the [Linux toolchain and CPU requirements](installation.md). Use a fresh
build directory if a cache contains retired options or dependencies from another
configuration. Offline overrides must be clean Git checkouts at the documented
exact commits; arbitrary ZIPs or edited dependencies are rejected. A genome-only
binary supports genome commands but cannot run `align -c`.

## Input validation fails

Check that inputs are FASTA, record IDs are unique, every record has bases, and
compression is valid. Sequence whitespace/digits and empty records are rejected;
non-ACGT letters become `N`. For `align -c`, supply exactly one record on each
side. Delta output requires input paths without spaces. A loaded index must
match the supplied normalized reference; an incompatible index is not rebuilt.

## Outputs already exist

Choose a new output prefix and work directory. Do not reuse an existing
`batch.tsv` or per-query output directory. A handled output failure rolls back
files created by that run while preserving existing files. A forced process kill
can leave a subset of files, so check the process exit code and log before using
outputs. No alignment checkpoint/resume interface is provided.

## No matches, partial coverage, or unresolved gaps

A valid zero-hit result can succeed with empty PAF. Centromere output covers
supported paths and internal gaps; unknown directions, ambiguous copies and
unresolved regions remain in diagnostic TSVs. Inspect gap reasons and the
summary, and interpret `complete` at that documented scope. Do not infer full
coverage from success or from the first/last reported coordinates.

Centromere standard outputs can reuse coordinates for accepted copy paths.
Filtering all SAM secondary records can remove these relationships. See
[formats](formats.md) before converting or combining results. Ordinary genome
`align` selects a linked one-to-one main path, while `batch` retains candidate
records, so their counts and repeat coverage can differ.

## Memory use or runtime is high

Reference/query sequence collections and indexes use memory. More threads can
increase concurrent search workspaces. `smem` and `maxmatch` may enumerate many
repeat occurrences. Reuse a saved reference index for repeated genome queries;
batch keeps one resident index but can have a higher alignment-stage memory
peak. There is no CLI-wide memory cap or query-concurrency setting.
A centromere gap that exceeds its resource budget remains unresolved; the
program does not silently expand that budget or invent a completed gap.

## Progress or interruption

`--progress auto` is quiet without a terminal; use `--progress on` to force
reports or inspect the log. On POSIX systems, SIGUSR1 requests a snapshot;
SIGINT/SIGTERM request cooperative interruption and exit 130/143. A dependency
call can delay interruption until it returns. A forced kill is outside handled
rollback. A separately completed saved index survives a later alignment error.

## Report a problem

Check the [exit-code table](user-guide.md#logs-and-exit-status). Include the
version, command, platform, relevant log/error and a small shareable input when
opening a GitHub issue. Remove private paths or sequences before posting.

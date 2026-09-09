# RaMA-G user guide

## Purpose and installation

RaMA-G performs assembly-to-assembly nucleotide alignment. The reference is the
target coordinate system; the query is aligned in both orientations. It is not a
short-read mapper. Linux x86-64, a C++20 compiler, CMake 3.22+, Git, zlib and C/C++
OpenMP are the production build requirements.

```bash
cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release \
  -DRAMAG_WARNINGS_AS_ERRORS=ON
cmake --build build-release -j 8
./build-release/ramag --version
./build-release/ramag --help
ctest --test-dir build-release --output-on-failure
```

Internal checkouts contain tests; public release trees need not. CTest verifies
only tests available in the selected checkout. Absence of tests does not prevent
building the program. Sufkit is fixed to
`f8c4c386ee883e45ad0f973efc4c8e1148b0068a`; a local override must have that exact
HEAD and a clean worktree. Disable pinned dependencies only for internal oracle
development. RaMA-G no longer depends directly on SeqPro.

## FASTA inputs

Reference and query both accept plain and gzip FASTA, including `.fa`, `.fasta`,
`.fna` and gzip variants. Compression is detected from bytes, not just suffixes.
A gzip magic header without a gzip suffix is accepted; a gzip-looking suffix
without valid gzip content fails. BGZF and concatenated gzip members are read
sequentially. No external decompressor, FAI or `.gzi` file is used or created.
Existing adjacent FAI files are ignored and left unchanged.

kseq parses records and zlib reads the stream. The input wrapper rejects FASTQ,
sequence before the first header, duplicate IDs, empty records, nonletter bases,
damaged/truncated gzip and unsupported `.bz2`, `.xz`, `.zst` or `.zip` input.
IDs are the first whitespace-delimited header token; the complete trimmed header
is retained. Empty lines within records are allowed; a leading blank line is
rejected. CRLF and LF are accepted. Sequence spaces and digits are rejected.

Letters become uppercase `A/C/G/T/N`; other letters, including IUPAC ambiguities,
`U` and `X`, become `N`. `N` is a hard break for seed matching. Record boundaries
are hard breaks. Total coordinates remain 64-bit; the frozen kseq parser has an
`int` return length, so a single record longer than `INT_MAX` bases is rejected
explicitly rather than silently truncated. Multiple records can total more than
that limit. The input collection is materialized in memory.

When delta output is requested, reference and query paths cannot contain
whitespace because delta's first line separates the two paths with whitespace.

## Quick start and output selection

Create or supply `reference.fa` and `query.fa` before running these examples:

```bash
# Legacy interface: defaults to SAM, PAF and delta.
ramag align --reference reference.fa --query query.fa \
  --output-prefix result --work-dir work

# Explicit suffix-driven selection.
ramag align --reference reference.fa.gz --query query.fa.gz \
  --output result.paf --output result.maf --output result.chain --work-dir work
```

Explicit outputs must share one directory and base prefix. Each of `.sam`,
`.paf`, `.delta`, `.maf`, `.chain` can occur once. Unknown suffixes, compressed
outputs, duplicate paths/formats and mixed `--output` with `--output-prefix` or
`--formats` fail before alignment. With the legacy interface, `--formats` accepts
a comma-separated subset of these five formats.

Output files are never overwritten. Use a new prefix for another run. All
requested formats are written privately, validated, then published. A handled
writer/publication failure rolls back this run's published files. Existing files
owned by another run are preserved. There is no all-files atomic commit after an
unhandled crash or `SIGKILL`; files alone cannot prove a successful command.

## Persistent reference indexes

```bash
# Build in memory, save that same index, then align.
ramag align --reference reference.fa --query query.fa \
  --save reference.sufidx --output saved-result.paf --work-dir save-work

# Build only an index.
ramag index --reference reference.fa.gz --output other-reference.sufidx \
  --work-dir index-work --threads 8

# Load an existing index without rebuilding it.
ramag align --reference reference.fa --reference-index reference.sufidx \
  --query another-query.fa --output another-result.paf --work-dir another-work
```

The persistent artifact is a single `.sufidx`. No companion manifest or marker
is required. Old raw-LCP and compressed-LCP indexes remain loadable if Sufkit
accepts their format and capabilities; any historical companions are ignored.
New indexes default to byte-coded LCP, with full SA, ISA and suffix links.

The reference FASTA remains required to provide bases for extension and output.
Sufkit checks section CRC, suffix-array structure and resource capabilities.
RaMA-G also compares stored contig metadata and Sufkit's reference fingerprint
with the supplied normalized reference. SHA-256 is no longer calculated; CRC and
the existing Sufkit fingerprint checks remain. A mismatch fails without rebuild.

`--save` and `--reference-index` are mutually exclusive. Saving uses the same
in-memory object, then self-validates and atomically publishes it. An index that
was successfully saved is retained if later alignment fails. Without `--save`,
there is no persistent index Save or save/load self-validation. RaMA-G neither
discovers indexes automatically nor silently converts loaded LCP encodings.

## Seed and selection modes

| Mode | Meaning |
|---|---|
| `fast` (default) | MAM skeleton with whole-query MEM completion/filtering |
| `mumreference` | Reference-unique MAM; query repetition is allowed |
| `mum` | Strict uniqueness in combined reference and current complete query record |
| `smem` | Generalized `(l,c)` intervals; expand every retained reference occurrence |
| `maxmatch` | Occurrence-complete MEM enumeration |

`--min-match` is minimum seed length and SMEM's `l`.
`--smem-min-occurrences` defaults to 1 and sets SMEM's `c`; explicitly using it
outside `smem`, or setting it to zero, is an error. Strict MUM and SMEM enumeration
uses complete query records and is serial at the Sufkit enumeration layer.
MUMREFERENCE uses stable 4 MiB tiles plus boundary recovery and whole-record
maximality/reference-uniqueness verification.

`--selection-mode all` is the default and exposes valid pairwise candidates.
`one-to-one` uses the current pairwise reference-side and query-side DP selection
and retains their intersection. It does not use the former reciprocal elementary
interval selector and does not claim equivalence to MUMmer4 `delta-filter -1`.
The internal residual recovery and guarded gap-fill build settings remain as
before this input/logging migration; neither is silently enabled by this change.

## Alignment options

| Option | Default | Effect |
|---|---:|---|
| `--threads` | 1 | Shared worker budget |
| `--min-match` | 20 | Minimum exact seed length |
| `--max-gap` | 90 | Cluster/chaining gap constraint |
| `--diag-diff` | 5 | Absolute diagonal slack |
| `--diag-factor` | 0.12 | Relative diagonal slack |
| `--min-cluster` | 65 | Minimum matching support |
| `--break-length` | 200 | Compatibility parsing only; nondefault rejected |
| `--max-dp-cells` | 4000000 | Compatibility parsing only; nondefault rejected |

KSW2 uses the current pairwise scoring, certified global bandwidth expansion and
endpoint extension. CLI parameters do not select an alternate DP backend.
Use the same command with `--print-effective-config` to inspect its configuration
without running alignment. This JSON printed to stdout remains supported; it is
not a runtime manifest. No log is created for help, version or effective config.

## Formats and coordinate interpretation

All writers consume the same finalized alignment records. SAM uses one-based
reference POS, CIGAR/clipping and reverse flags. PAF uses zero-based half-open
forward-coordinate intervals and an explicit query strand; `cg` carries CIGAR.
SAM `NM`/`MD` and PAF scoring/edit tags describe the alignment, not confidence.
MAPQ is currently 255 (unavailable), not a calibrated probability.

Delta uses one-based inclusive endpoints, reversed query endpoints for minus
alignments, and indel displacement encoding. Pairwise MAF uses a positive target
row and oriented query row; minus-query start is `query_size - query_end`.
UCSC chain uses reference as positive target and the query's actual orientation.
Chain removes unrepresentable terminal gaps and normalizes internal block gaps.
MAF and chain preserve the emitted record set rather than filtering secondary
records again; SAM-specific primary/supplementary status does not imply filtering
in another format.

## Logs, progress and signals

CLI execution creates `WORK/runs/<run-id>/run.log` and prints its location to
stderr. spdlog records timestamps and levels, invocation, build identity, effective
configuration, input counts, index action/encoding, stage timings, core statistics,
requested artifacts and final status/exit code. Logs are flushed on exit and
retained after failure. Failure before a logger can be created is reported to
stderr. There are no SHA-256 digests, runtime manifests or completion markers.

`--progress auto` shows terminal progress only when stderr is a TTY. `on` forces
periodic progress; `off` suppresses it. The default interval is 10 seconds; an
explicit interval must be at least 1 and cannot accompany `off`. Stage logs and
final/error messages remain available even when periodic progress is off.
Reports include stage, elapsed time, known counts, thread budget and Linux RSS.
Third-party operations without measurable totals report heartbeats.

```bash
ramag align --reference reference.fa --reference-index reference.sufidx \
  --query query.fa --output job.paf --work-dir job-work --threads 8 \
  --progress on >job.stdout 2>job.stderr &
pid=$!
kill -USR1 "$pid"  # request a snapshot; does not terminate
wait "$pid"
status=$?
test "$status" -eq 0 && test -s job.paf
```

SIGINT and SIGTERM request cooperative cancellation and exit 130/143. An internal
Sufkit/CaPS/divsufsort call may delay cancellation until it returns control.
Signal handlers only set flags. No resume/checkpoint interface is provided;
rerun with a new output prefix and reuse a successfully saved index if available.

## Large references and multiple queries

Define `REF`, `QUERY`, `INDEX`, `PREFIX`, and `WORK` for your own data and paths.

```bash
ramag index --reference "$REF" --output "$INDEX" --work-dir "$WORK/index" --threads 32
ramag align --reference "$REF" --query "$QUERY" --reference-index "$INDEX" \
  --output "$PREFIX.paf" --output "$PREFIX.maf" --work-dir "$WORK/align" \
  --threads 32 --seed-mode mumreference --selection-mode one-to-one \
  --print-effective-config
ramag align --reference "$REF" --query "$QUERY" --reference-index "$INDEX" \
  --output "$PREFIX.paf" --output "$PREFIX.maf" --work-dir "$WORK/align" \
  --threads 32 --seed-mode mumreference --selection-mode one-to-one --progress on
```

For multiple queries, repeat alignment with the same reference/index and separate
query, output prefix and work directory. Index loading and reference validation
remain part of each command's elapsed time. CPU affinity must permit the requested
thread count; the executable preserves its launch-authorized CPU set.

## Exit codes and migration

| Code | Meaning |
|---|---|
| 0 | Requested artifacts validated and published |
| 2 | CLI/configuration error |
| 3 | FASTA/input error |
| 4 | Unsupported input/mode |
| 5 | Dependency/index error |
| 6 | Compute/resource error |
| 7 | Output, validation or logging error |
| 8 | Unexpected internal exception |
| 130 / 143 | Cooperative SIGINT / SIGTERM interruption |

Replace automation that waits for `.complete` or reads `.manifest.json` with
process-exit checks, output parsing and `run.log` inspection. Index bundles become
single files; old companion files can be left in place. Library consumers of
removed manifest/SHA-256/SeqPro headers must migrate to `RunStatistics`, `ReadFasta`
or the existing pairwise interfaces. No version-number change is made here.

Current limitations include no automatic comparison/report pipeline, automatic
index discovery, random-access gzip index, compressed output, BAM/CRAM, VCF,
coverage/HTML report, calibrated MAPQ or cross-query-record strict-MUM uniqueness.
Performance conclusions require separately accepted measurements.

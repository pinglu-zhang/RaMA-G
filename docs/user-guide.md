# RaMA-G user guide

## Purpose and installation

RaMA-G performs assembly-to-assembly nucleotide alignment. The reference is the
target coordinate system; the query is aligned in both orientations. It is not a
short-read mapper. Linux x86-64, a C++20 compiler, CMake 3.22+, Git, zlib and C/C++
OpenMP are the production build requirements.

These instructions cover the pending 0.1.1 source version. Run the following
commands from the source root, using a fresh build directory after upgrading:

```bash
cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release \
  -DRAMAG_WARNINGS_AS_ERRORS=ON
cmake --build build-release -j 8
./build-release/ramag --version
./build-release/ramag --help
ctest --test-dir build-release --output-on-failure
```

The project version supplies the CLI version, SAM `@PG VN` and installed CMake
package version. To install the program and library under a local prefix:

```bash
cmake --install build-release --prefix "$PWD/install"
./install/bin/ramag --version
```

Examples below use `ramag`; use the executable's full path or add its installed
`bin` directory to your `PATH`. See the [changelog](../CHANGELOG.md) for changes.

Internal checkouts contain tests; public release trees need not. CTest verifies
only tests available in the selected checkout. Absence of tests does not prevent
building the program. Sufkit is fixed to
`028075e6f2d622fcbdcf76b153bbde0f069ca64f`; a local override must have that exact
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

`--threads` also supplies the index-loading worker budget for SA/ISA consistency
validation and Fast prefix-directory construction. On Linux, the requested
budget must fit the calling thread's allowed CPU set; check your scheduler or
`taskset` allocation if loading reports insufficient allowed CPUs. Not every
loading phase is parallel, and this option does not promise linear speedup.
CRC is accumulated while sections are read, and the normalized reference is
checked through read-only views. Full validation is retained. Upgrading to 0.1.1
does not require rebuilding compatible old raw- or byte-coded-LCP indexes.

`--save` and `--reference-index` are mutually exclusive. Saving uses the same
in-memory object, then self-validates and atomically publishes it. An index that
was successfully saved is retained if later alignment fails. Without `--save`,
there is no persistent index Save or save/load self-validation. RaMA-G neither
discovers indexes automatically nor silently converts loaded LCP encodings.

## Batch queries

`batch` reads one reference and initializes one shared Sufkit index, then processes
query FASTAs in input order. Each query remains an independent alignment, including
its FASTA record boundaries and strict-MUM uniqueness rules.

```bash
ramag batch --reference human.fa.gz --reference-index human.sufidx \
  --query chimp.fa.gz --query gorilla.fa.gz --query orangutan.fa.gz \
  --output-dir results --threads 32

ramag batch --reference human.fa.gz --reference-index human.sufidx \
  --seqfile queries.txt --output-dir named-results --threads 32
```

The seqfile is a plain-text file with one query name and FASTA path per line,
separated by whitespace. `.txt` is the recommended suffix, but parsing does not
depend on the filename suffix.

An editable [seqfile example](../examples/batch/queries.txt) is included; replace
its placeholder paths with your query FASTAs.

```text
# name      FASTA path
chimp       data/chimp.fa.gz
gorilla     data/gorilla.fa.gz
orangutan   data/orangutan.fa.gz
```

Paths are relative to the seqfile directory, not the process working directory.
The first whitespace separator ends the name; the remaining trimmed text is the
path and may contain internal spaces. LF, CRLF, blank lines and whole-line `#`
comments are supported. Shell variables, glob patterns and quoting are not
expanded. This is a query list, not a Cactus tree/seqfile. The reference is supplied
separately with `--reference`. `--seqfile` and repeated `--query` cannot be mixed.

Names may use ASCII letters, digits, `.`, `_` and `-`, and must be unique. Empty
names, `.`, `..`, `batch.tsv` and `.ramag-work` are rejected. With repeated
`--query`, names come from the filename after stripping `.gz` and `.fa`, `.fna`
or `.fasta`. If filenames collide or contain unsupported name characters, use
unique explicit names in a seqfile. Missing/unreadable query files are individual
failures; malformed seqfiles and duplicate names fail before index initialization.

Each item writes `OUTPUT_DIR/NAME/alignment.paf` by default. `--formats` accepts
the same five formats as `align`; seed, selection, scoring and thread defaults
also match `align` (including `fast` and `all`). To select a different policy,
pass `--seed-mode mumreference --selection-mode one-to-one` explicitly. Batch
does not accept `--output` or `--output-prefix`. The output directory is required;
the optional `--work-dir` defaults to `OUTPUT_DIR/.ramag-work`.

Without `--reference-index`, batch builds one ephemeral byte-coded-LCP index.
Use `--save PATH` to save that same object before processing queries; saving also
performs the usual temporary-file Load self-validation. `--save` and
`--reference-index` are mutually exclusive. A saved index survives later query
failure. Compatible old raw/byte-coded indexes remain usable.

Queries run serially, each using the requested thread budget. Only one query's
sequence, seeds and alignment workspace are retained at a time, but the reference
and index stay resident until batch ends. Standalone `align` releases its index
after enumeration; batch can therefore have a higher pairwise-stage memory peak.
There is no automatic query prefetch, parallel-query execution or resume.
The existing bounded per-thread encoding/selection scratch caches may be reused;
batch does not retain prior query records, seeds or alignment results.
The reference is held by an immutable shared owner and checked once during
index initialization. Later queries reuse that validated binding without
rescanning or copying the reference. Query input is still validated independently.
On glibc, batch returns freed allocator pages between queries. This does not
discard live index pages or impose an RSS cap. The extra resident index during
pairwise processing means that a strict standalone-align peak-memory bound is
not guaranteed. No whole-batch thread pool or memory-admission scheduler is used;
OpenMP parallelism remains within the current query. The experimental
`--query-concurrency` and `--memory-budget` options are not supported.

`OUTPUT_DIR/batch.tsv` is created exclusively and flushed after each terminal
query outcome. Columns are `order`, `name`, `query_path`, `status`, `exit_code`,
`elapsed_seconds`, `alignment_records`, `output_dir`, `log_path`, and `message`.
Tabs, newlines, carriage returns and backslashes inside field values are escaped
as `\t`, `\n`, `\r` and `\\`. Per-query elapsed time excludes shared initialization.
It is a status table, not a completion marker; while a query runs its terminal
row is absent. An abrupt process kill may leave only the completed rows.

Statuses are `success`, `failed`, `interrupted` and `not-run`. Valid zero-hit
outputs are successful with zero records. On early termination, remaining rows
are appended as `not-run` with an empty exit-code field. All successful query
outputs are preserved; handled failures roll back only the current item's files.
Existing output files and an existing `batch.tsv` are never overwritten, and
all planned output collisions are checked before initialization.

The batch exits 0 only if all queries succeed. Ordinary failures are recorded
and later queries continue; the final exit code is that of the first failed
query in input order. Shared initialization errors, allocation failure or an
unreliable batch status table stop the batch. SIGINT/SIGTERM stop scheduling,
roll back the interrupted item and exit 130/143. SIGUSR1 requests a snapshot
containing the current query name, position and processing phase.

The batch log records shared reference input and index initialization separately
from query execution. `index_load_calls` and `index_build_calls` count shared
initialization calls, not the Load self-check used when saving. Each query logs
`index_action=reused` with zero initialization calls. Per-query logs live below
`WORK/queries/NAME/runs/RUN_ID/run.log`. A command with
`--print-effective-config` prints the resolved list and paths without reading
FASTA contents or initializing an index.

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
The residual recovery and guarded gap-fill build settings retain their existing
defaults in 0.1.1. The index-loading update changes neither selection nor scoring.

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

For reused indexes, logs include `index_action=loaded`, the requested loading
thread budget, LCP encoding, logical read bytes and CRC CPU time. Timing summaries
separate the Sufkit load, reference validation and detailed loading phases.
CRC computation is included in section-read times; do not add it again to the
load total or interpret CPU seconds as elapsed wall time.

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
or the existing pairwise interfaces. These input/logging migration changes
precede 0.1.1; they are not additional interface removals in this patch release.

Current limitations include no automatic comparison/report pipeline, automatic
index discovery, random-access gzip index, compressed output, BAM/CRAM, VCF,
coverage/HTML report, calibrated MAPQ or cross-query-record strict-MUM uniqueness.
Performance conclusions require separately accepted measurements.

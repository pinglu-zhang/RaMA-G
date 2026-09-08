# RaMA-G User Guide

> Experimental index-build branch: the pinned Sufkit candidate is a
> server-local, unpushed commit. Until review and publication, configure with
> `-DRAMAG_SUFKIT_SOURCE_DIR="$SUFKIT_CHECKOUT"`, where that clean checkout is
> at the exact candidate SHA listed below. Fetching this commit from the
> upstream repository is not yet a reproducible public installation path.
> Existing indexes created by the original Sufkit
> `bdb67c6de5daddd8a005640de73d96549d2575f4` are explicitly accepted after
> the same reference, format, capability, CRC, and companion checks.

This guide describes the behavior implemented by RaMA-G 0.1.0. It is a
runtime manual, not a roadmap: options or workflows that are not implemented
are identified as limitations rather than presented as available features.

## 1. Purpose and supported use case

RaMA-G performs deterministic pairwise nucleotide alignment between one
reference genome or assembly and one query genome or assembly. Both inputs may
be multi-FASTA files. Their roles are explicit and are never exchanged based on
length, N50, contig count, or any other heuristic.

The production baseline is Linux x86-64 with C++20 and OpenMP. RaMA-G is an
assembly-to-assembly and genome-to-genome aligner. It is not a short-read or
paired-end mapper, a variant caller, a pangenome graph constructor, or an
automatic genome-comparison report generator.

One alignment run performs these operations:

1. validate and normalize the reference and query;
2. build or load a complete Sufkit suffix-array index for the reference;
3. enumerate exact seeds on both query strands;
4. merge compatible exact seeds on the same diagonal;
5. build deterministic chains and extend gaps;
6. resolve duplicate records and apply the requested selection mode;
7. write and validate every requested output;
8. publish a provenance manifest and, last, a completion marker.

## 2. Installation and build

### 2.1 Requirements

The normal production build requires:

- CMake 3.22 or newer;
- a C and C++ compiler with C++20 support;
- Git, including network access unless dependency sources are supplied locally;
- zlib development headers and library;
- OpenMP support for both C and C++ on Linux;
- enough storage for the build, run-local work directory, selected outputs,
  and an optional persistent suffix-array index.

On Linux, configuration fails when either the C or C++ OpenMP component is
unavailable. RaMA-G does not accept `--threads` and silently compile its
parallel stages away.

### 2.2 Convenience build

From the repository root:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

The executable is `build/ramag`.

### 2.3 Fresh release-quality build

Use a new build directory and treat warnings as errors:

```bash
cmake -S . -B build-release \
  -DCMAKE_BUILD_TYPE=Release \
  -DRAMAG_WARNINGS_AS_ERRORS=ON
cmake --build build-release -j
ctest --test-dir build-release --output-on-failure
```

Do not reuse a build cache after changing either pinned dependency or its
source override.

### 2.4 Dependency identity

The production build pins exact Git commits:

| Dependency | Required identity |
|---|---|
| Sufkit | `f8c4c386ee883e45ad0f973efc4c8e1148b0068a` (0.3.0) |
| SeqPro | `6781cadcf81a0da53d7573444594c1484947017c` |

The default CMake route fetches these identities. Local sources can be supplied
with `RAMAG_SUFKIT_SOURCE_DIR` and `RAMAG_SEQPRO_SOURCE_DIR`, but each override
must be a clean Git checkout at the complete required 40-character commit.
Branches, moving tags, abbreviated commits, dirty worktrees, and the wrong
commit are rejected.

For example:

```bash
cmake -S . -B build-release \
  -DCMAKE_BUILD_TYPE=Release \
  -DRAMAG_WARNINGS_AS_ERRORS=ON \
  -DRAMAG_SUFKIT_SOURCE_DIR=/path/to/clean/sufkit \
  -DRAMAG_SEQPRO_SOURCE_DIR=/path/to/clean/seqpro
```

`-DRAMAG_ENABLE_PINNED_DEPENDENCIES=OFF` selects the small bundled oracle
baseline. It is useful for isolated core development but is not a production
or benchmark build and reports that fact in its manifest.

### 2.5 Inspect the binary

```bash
./build-release/ramag --version
./build-release/ramag --help
```

The version line records the RaMA-G version and commit, build type, compiler,
OpenMP state and runtime, configured thread information, and whether Sufkit's
bundled divsufsort C objects were built with OpenMP.

## 3. Input contract

### 3.1 Reference and query

Both paths are mandatory for alignment:

```text
--reference PATH
--query PATH
```

The reference is indexed and is always the target in PAF and UCSC chain. Query
coordinates are always reported with respect to the original forward query
record, even when an alignment is on the reverse strand.

Explicit self-alignment is allowed, but the reference and query roles remain
distinct.

### 3.2 Plain and gzip FASTA

RaMA-G accepts plain FASTA and gzip FASTA, including common names such as:

```text
reference.fa
reference.fasta
reference.fna
reference.fa.gz
reference.fasta.gz
reference.fna.gz
```

Compression is detected from gzip magic bytes and then validated as a stream:

- a gzip stream without a `.gz` suffix is accepted;
- a `.gz` path without gzip magic is rejected;
- truncated streams and CRC/trailer errors are rejected;
- concatenated gzip members are accepted;
- BGZF-compatible input is read sequentially as gzip;
- bzip2, xz, zstd, and zip are rejected explicitly.

The gzip route does not invoke an external `gzip` command. It uses the zlib
library linked into RaMA-G. It is sequential and does not use a `.gzi` random
access index.

### 3.3 FASTA records and normalization

For each header, the first whitespace-delimited token is the sequence name;
the complete text after `>` is retained as the header description. Record
names must be unique within each input.

Alphabetic sequence characters are accepted and normalized as follows:

- `A`, `C`, `G`, and `T`, in either case, become uppercase canonical bases;
- every other alphabetic symbol becomes `N`;
- a non-alphabetic sequence character is an input error.

`N` is a hard break for exact seeding. A seed cannot cross an ambiguous base or
a contig boundary. Empty input, sequence text before a header, a missing header
identifier, a duplicate identifier, and an empty record are errors.

### 3.4 Adjacent FAI files

For plain FASTA, an existing standard `<FASTA>.fai` may be used through the
SeqPro path. RaMA-G treats it as read-only: it hashes and copies the file into
the run directory, validates the copy, and records the source and verification
route. It never creates, repairs, or replaces an index beside the input.

When no adjacent FAI exists, any input metadata created by the run remains
under `--work-dir`. Gzip input is handled by the zlib route rather than the
plain-file SeqPro/FAI route.

### 3.5 Path restriction when delta is selected

The MUMmer delta header contains the reference and query paths separated by a
space and has no path-escaping syntax. If delta output is selected, the
normalized absolute reference and query paths must not contain whitespace.
Omit delta when such paths cannot be avoided.

## 4. Quick start

The following small example is portable. Save these files in an empty
directory:

`reference.fa`:

```text
>reference
ACGTTGCACTGATCGTACGATTCGGAACCTAG
```

`query.fa`:

```text
>query
ACGTTGCACTGATCATACGATTCGGAACCTAG
```

The sequences are short, so use small seed and cluster thresholds:

```bash
mkdir -p tutorial-work
./build-release/ramag align \
  --reference reference.fa \
  --query query.fa \
  --output-prefix tutorial \
  --work-dir tutorial-work \
  --threads 2 \
  --min-match 6 \
  --min-cluster 6
```

The legacy prefix interface creates these default outputs:

```text
tutorial.sam
tutorial.paf
tutorial.delta
tutorial.manifest.json
tutorial.complete
```

A successful stdout message reports the alignment count and marker path. Check
the marker rather than inferring completion from another file:

```bash
test -s tutorial.complete
```

RaMA-G refuses to overwrite any of these final paths. Choose a new prefix or
move the existing result set before rerunning.

### 4.1 Explicit output paths

Select formats by repeating `--output`:

```bash
mkdir -p explicit-work results
./build-release/ramag align \
  --reference reference.fa \
  --query query.fa \
  --output results/example.paf \
  --output results/example.maf \
  --output results/example.chain \
  --work-dir explicit-work \
  --threads 2 \
  --min-match 6 \
  --min-cluster 6
```

Every explicit path must have the same directory and base prefix. In this
example the common prefix is `results/example`, so the shared metadata files
are `results/example.manifest.json` and `results/example.complete`.

The explicit interface is mutually exclusive with `--output-prefix` and
`--formats`.

## 5. Reusable reference index

### 5.1 Build an index bundle

```bash
mkdir -p index-work
./build-release/ramag index \
  --reference reference.fa \
  --output reference.sufidx \
  --work-dir index-work \
  --threads 8
```

Success atomically publishes:

```text
reference.sufidx
reference.sufidx.manifest.json
reference.sufidx.complete
```

All three files form one bundle. RaMA-G refuses to overwrite any existing
member. A `.sufidx` without the RaMA-G companion manifest and completion marker
is not accepted as a reusable production index.

On Linux, RaMA-G also verifies the CPU allocation before starting Sufkit.
OpenMP runtimes may bind the initial thread before `main()` when
`OMP_PROC_BIND` is enabled, while the CaPS scheduler creates ordinary worker
threads that inherit the caller's affinity. For index construction, RaMA-G
captures the launch CPU mask in an ELF pre-initialization hook, before libgomp
initializes. It restores that captured set before CaPS starts, respecting the
allocation supplied by `taskset` or the scheduler. If it has fewer logical CPUs than
`--threads`, the command fails before the expensive build instead of silently
running the requested workers on one core.

The index manifest records the launch and pre-Sufkit CPU sets, selected
backend, Sufkit build time, save time, validation time, publication time, and
total command time. These phases are separate: a fast in-memory suffix-array
construction does not imply that serializing and validating a large persistent
index takes no additional time.

### 5.2 Load the index

```bash
mkdir -p indexed-work
./build-release/ramag align \
  --reference reference.fa \
  --reference-index reference.sufidx \
  --query query.fa \
  --output-prefix indexed \
  --work-dir indexed-work \
  --threads 8 \
  --min-match 6 \
  --min-cluster 6
```

`--reference` remains mandatory because extension and output serialization need
the normalized reference bases and catalog. Before seeding, RaMA-G verifies:

- the index, manifest, and marker are all present and mutually bound;
- the index is a standalone Sufkit full suffix array with sampling rate 1;
- the required Sufkit version, complete commit, profile, and capabilities match;
- Sufkit's format, section CRCs, suffix-array permutation, and internal
  semantics are valid;
- reference contig count, order, names, full descriptions, lengths, and
  normalized ambiguous-base statistics match;
- the SHA-256 of the normalized reference catalog and bases matches.

Because the content check occurs after normalization, case changes and
equivalent ambiguity symbols may reuse an index if the ordered normalized
`A/C/G/T/N` records, names, and descriptions are identical. The original
source path, size, compression, and file identity remain separate provenance.

Without `--reference-index`, alignment builds an ephemeral index for that run.
RaMA-G does not automatically search for a cache and does not automatically
save an ephemeral index.

## 6. Seed modes

All modes search both query strands. Coordinates are converted back to the
original forward query record before downstream processing.

| Mode | Exact-match semantics | Current execution notes |
|---|---|---|
| `fast` | Reference-unique MAM skeleton, followed by complete whole-query MEM enumeration and filtering against skeleton-covered query intervals | Default; designed as the current hierarchical route |
| `mumreference` | A two-sided maximal exact match whose sequence occurs once in the combined reference; it may repeat in the query | Stable 4 MiB MAM tiles plus boundary MEM recovery and whole-record revalidation |
| `mum` | A two-sided maximal exact match unique in the combined reference and in the current complete query record | Whole-record, record-local, serial Sufkit enumeration |
| `smem` | A qualifying `(l,c)` query interval not properly contained in another qualifying interval, expanded to every reference occurrence | Whole-record, serial Sufkit enumeration |
| `maxmatch` | Occurrence-complete two-sided maximal exact matches without a uniqueness requirement | Can create many seeds in repetitive input |

`--min-match N` supplies the minimum exact-match length for every mode and must
be at least 1. The default is 20.

For SMEM only, `--smem-min-occurrences N` supplies `c`, the minimum number of
reference occurrences, and must be at least 1:

```bash
./build-release/ramag align \
  --reference reference.fa \
  --query query.fa \
  --output result.paf \
  --work-dir smem-work \
  --seed-mode smem \
  --smem-min-occurrences 2
```

Explicitly using `--smem-min-occurrences` with another seed mode is a
configuration error. No seed mode is universally best; use controlled data and
an explicit evaluation contract when comparing modes. The default remains
`fast`.

## 7. Alignment selection

`--selection-mode` controls record resolution after chaining and extension:

| Value | Behavior |
|---|---|
| `all` | Remove exact duplicate records, retain all remaining records, and mark one deterministic primary alignment per query contig |
| `one-to-one` | Independently select the best-supported record covering each elementary interval on the reference and query, then retain only records selected on both sides |

The default is `all`.

For `all`, primary ranking prefers higher alignment score, longer aligned query
span, higher exact-match identity, lower reference numeric ID, lower reference
start, lower query start, and then strand order. Other retained records are
supplementary in SAM and carry `tp:A:S` in PAF.

For `one-to-one`, interval priority is determined by exact matching bases,
identity, interval span, and a stable alignment-index tie-break. The selected
set is the intersection of the independently selected reference-side and
query-side records. Every retained record is marked primary because competing
records have already been removed.

This is RaMA-G's reciprocal interval-based selection. It is not specified or
advertised as an implementation of MUMmer4 `delta-filter -1`, whose filtering
and LIS behavior are different.

Example:

```bash
./build-release/ramag align \
  --reference reference.fa \
  --query query.fa \
  --output reciprocal.paf \
  --work-dir reciprocal-work \
  --selection-mode one-to-one \
  --min-match 6 \
  --min-cluster 6
```

## 8. Alignment and resource options

| Option | Default | Effect |
|---|---:|---|
| `--threads N` | `1` | Shared worker budget for supported input, index, seed, chaining, extension, and writer stages |
| `--min-match N` | `20` | Minimum exact seed length |
| `--max-gap N` | `90` | Maximum gap allowed between seeds considered for one chain |
| `--diag-diff N` | `5` | Fixed component of the permitted diagonal difference between chained seeds |
| `--diag-factor X` | `0.12` | Relative component of diagonal tolerance; must be non-negative |
| `--min-cluster N` | `65` | Minimum reference/query span required for a seed chain to become an alignment candidate |
| `--break-length N` | `200` | Chain split threshold; together with `--max-gap`, bounds the separation between consecutive seeds |
| `--max-dp-cells N` | `4000000` | Hard maximum cell count for one scalar gap-DP matrix |

`--threads`, `--min-match`, `--min-cluster`, `--break-length`, and
`--max-dp-cells` must be positive. Resource limits fail the run rather than
silently truncating seeds or returning a partial successful result.

Before a long run, validate and print the normalized configuration without
performing alignment:

```bash
./build-release/ramag align \
  --reference reference.fa \
  --query query.fa \
  --output result.paf \
  --work-dir result-work \
  --threads 32 \
  --print-effective-config
```

The same option is available for `ramag index`.

## 9. Output formats

RaMA-G finalizes one ordered `AlignmentRecord` collection before invoking any
writer. Writers do not choose their own records, so all requested formats
describe the same selected alignments.

### 9.1 Output selection

The legacy interface accepts a comma-separated list:

```bash
--output-prefix result --formats paf,maf,chain
```

Without `--formats`, it selects exactly `sam,paf,delta`.

The suffix-driven interface accepts one instance of each format:

```bash
--output result.sam --output result.paf --output result.delta \
--output result.maf --output result.chain
```

Supported final suffixes are `.sam`, `.paf`, `.delta`, `.maf`, and `.chain`.
Compressed output names such as `.sam.gz` and `.paf.gz`, unknown suffixes,
missing suffixes, duplicate formats, duplicate paths, and different base
prefixes are rejected before computation.

### 9.2 Internal coordinate convention

Internally, all positions are unsigned 64-bit, zero-based, and half-open:

```text
[begin, end)
```

Reference coordinates increase. Query `begin < end` is always expressed on the
original forward query record; strand is stored separately.

### 9.3 SAM 1.6

- `@HD` reports `SO:unknown`.
- `@SQ` contains every reference record.
- `@PG` records RaMA-G, its version, and the command line.
- POS is one-based.
- reverse alignments use flag `0x10`.
- non-primary records in `all` mode use supplementary flag `0x800`.
- CIGAR uses `=`, `X`, `I`, and `D`, with hard clipping outside the aligned
  query span.
- `NM:i`, `MD:Z`, and `AS:i` are emitted.
- RNEXT, PNEXT, TLEN, SEQ, and QUAL use unmapped placeholders because this is
  an assembly alignment record, not a paired-read record.
- MAPQ is `255`, meaning unavailable/unknown; it is not a calibrated confidence
  value.

SAM imposes additional lexical restrictions on query and reference names. A
name that is valid FASTA may still be unrepresentable as SAM QNAME or RNAME;
in that case the SAM writer fails instead of renaming sequences silently.

### 9.4 PAF

PAF uses its standard 12 columns with zero-based half-open spans and the
original forward query interval. It adds:

```text
cg:Z: extended CIGAR
NM:i: edit distance
AS:i: alignment score
tp:A: P for primary, S for supplementary
```

MAPQ is `255`.

### 9.5 MUMmer-compatible delta

The first line records the normalized absolute reference and query paths, and
the second line is `NUCMER`. Alignment coordinates are one-based and inclusive.
Reference coordinates increase; reverse-query coordinates decrease in original
query space. `D` consumes reference with a query gap, while `I` consumes query
with a reference gap.

The output is designed for MUMmer-compatible consumers, but selection semantics
remain those chosen by RaMA-G.

### 9.6 Pairwise MAF v1

Each `AlignmentRecord` becomes one `a` block containing exactly two `s` rows.
The reference row is on `+`. A reverse query uses the MAF start
`query_size - query_end`, reports `-`, and contains the oriented gapped query
text. `size` counts non-gap bases and `srcSize` is the complete contig length.
The block score is the RaMA-G alignment score.

### 9.7 UCSC chain

The reference is the target and always has strand `+`; the query uses its
actual strand. Consecutive `=`/`X` operations form aligned blocks. `D` becomes
the target gap `dt`, and `I` becomes the query gap `dq`. Terminal gaps that
cannot be represented are removed from the chain span. An alignment with no
aligned block causes the writer to fail rather than disappear silently. Chain
IDs follow final stable record order beginning at 1, and negative alignment
scores are written as chain score 0.

## 10. Progress, signals, and long-running jobs

### 10.1 Progress modes

```text
--progress auto|on|off
--progress-interval SECONDS
```

`auto` is the default. It prints progress to stderr only when stderr is an
interactive terminal. `on` forces newline-delimited status suitable for a log,
and `off` suppresses periodic status. The default interval is 10 seconds and an
explicit value must be at least 1. `--progress off` cannot be combined with an
explicit interval.

Stage changes are reported immediately when progress is enabled. Current stage
names include configuration, input, index build or load, seed enumeration,
seed merge, chaining, extension, conflict resolution, individual writers,
validation, manifest, publication, completion-marker publication, and complete.
When a third-party index routine exposes no precise work total, RaMA-G reports
a heartbeat, elapsed time, thread budget, and current Linux RSS without
inventing a percentage.

Example with a captured log:

```bash
./build-release/ramag align \
  --reference reference.fa.gz \
  --reference-index reference.sufidx \
  --query query.fa.gz \
  --output result.paf \
  --work-dir result-work \
  --threads 32 \
  --progress on \
  --progress-interval 30 \
  >result.stdout.log 2>result.progress.log
```

### 10.2 Signals on POSIX systems

Given a running PID:

```bash
kill -USR1 PID
```

requests one immediate status snapshot without terminating the run.

```bash
kill -INT PID
kill -TERM PID
```

request cooperative interruption. `SIGINT` exits with 130 and `SIGTERM` with
143 after normal-control-flow cleanup. The signal handler itself only sets a
small flag; it does not perform file I/O or cleanup asynchronously.

Input loops, seed callbacks, chaining and extension task boundaries, writers,
validators, and publication check the interruption flag. Sufkit,
libdivsufsort, or CaPS code without a cancellation callback can react only after
it returns control, so the first signal may not stop an index build
immediately.

There is no checkpoint or `--resume`. A rerun receives a new run ID. Reusing a
previously completed explicit reference index avoids rebuilding that index but
does not resume the interrupted alignment.

## 11. Result integrity and provenance

### 11.1 What constitutes success

For prefix `result`, success requires:

- every requested alignment file;
- `result.manifest.json`;
- `result.complete`, published last.

An alignment file, manifest, directory, running or exited process, stdout
message, or temporary file is not sufficient completion evidence.

### 11.2 Transactional publication

RaMA-G writes each artifact to a private temporary sibling, closes it, validates
it, and then publishes it without replacing an existing target. After
publication it verifies sizes and parses each final output again. It publishes
the completion marker only after all requested outputs and the manifest have
survived final validation.

If a writer, validator, publication step, or interruption fails, temporary
files are removed and any files published by this run are rolled back. A
competitor's pre-existing or concurrently published file is not deleted.

### 11.3 Manifest field guide

The JSON manifest uses schema `ramag.run-manifest.v1`, schema version 1. It
records:

- run ID, timestamps, status, and exit code;
- binary path, RaMA-G version/commit, compiler, build type, flags, and OpenMP;
- exact Sufkit/SeqPro identities and dependency mode;
- effective paths, formats, seed/selection mode, thread budget, and numerical
  parameters;
- input catalogs, lengths, and ambiguous-base counts;
- detected compression and input/FAI provenance;
- whether the suffix array was built or loaded and its verified properties;
- actual input, index, seed, chain, and extension routes;
- requested and observed worker counts;
- MAM tile/boundary and chaining resource statistics;
- seed, chain, candidate, rejected-record, alignment, and extension counts;
- per-stage wall time;
- every final output path, byte count, validation state, and format.

Interrupted or failed runs do not publish a success manifest. Once the pipeline
has created its private run directory, its diagnostic is retained below (an
earlier CLI/preflight rejection may occur before a run ID exists):

```text
WORK_DIR/runs/RUN_ID/failure.json
```

It records status, run ID, observed UTC time, failed or interrupted stage,
message, and, for interrupted alignment runs, signal and exit code.

## 12. Large-genome operational recipe

Use separate, explicit paths:

```bash
REF=/data/reference.fa.gz
QUERY=/data/query.fa.gz
INDEX="$PWD/run-assets/reference.sufidx"
PREFIX="$PWD/results/reference-vs-query"
WORK="$PWD/work/reference-vs-query"
```

Build and test a fresh production binary:

```bash
cmake -S . -B build-release \
  -DCMAKE_BUILD_TYPE=Release \
  -DRAMAG_WARNINGS_AS_ERRORS=ON
cmake --build build-release -j
ctest --test-dir build-release --output-on-failure
./build-release/ramag --version
```

Build the reference index once:

```bash
mkdir -p "$PWD/run-assets" "$PWD/results" "$PWD/work/reference-index"
./build-release/ramag index \
  --reference "$REF" \
  --output "$INDEX" \
  --work-dir "$PWD/work/reference-index" \
  --threads 32 \
  --progress on \
  --progress-interval 30
test -s "${INDEX}.complete"
```

Validate the intended alignment configuration:

```bash
./build-release/ramag align \
  --reference "$REF" \
  --reference-index "$INDEX" \
  --query "$QUERY" \
  --output "${PREFIX}.paf" \
  --output "${PREFIX}.maf" \
  --output "${PREFIX}.chain" \
  --work-dir "$WORK" \
  --threads 32 \
  --seed-mode fast \
  --selection-mode all \
  --progress on \
  --progress-interval 30 \
  --print-effective-config
```

Remove only `--print-effective-config` to perform the run. Keep stdout and
stderr logs in paths that do not collide with the output prefix. After launch,
`kill -USR1 PID` requests a snapshot. Accept the result only when:

```bash
test -s "${PREFIX}.complete"
```

To align another query against the same reference, choose a new prefix and work
directory while reusing `--reference-index "$INDEX"`. RaMA-G revalidates the
index/reference binding for every run.

## 13. Exit codes and troubleshooting

| Exit | Category | Typical cause |
|---:|---|---|
| 0 | Success | Requested result bundle or reference index was validated and completed |
| 2 | CLI/configuration | Missing option, invalid number, conflicting output interfaces, invalid suffix or progress combination, or an existing alignment-output target detected during preflight |
| 3 | Input | Malformed FASTA, duplicate ID, empty record, invalid gzip, CRC failure, or unreadable input content |
| 4 | Unsupported | Unsupported seed/input capability or unsupported compression format |
| 5 | Dependency/index | Dependency identity/capability mismatch or invalid persistent index bundle |
| 6 | Compute/resource | Chaining/DP/resource bound, allocation failure, or alignment invariant failure |
| 7 | Output/validation | Writer, serialization, validator, publication failure, or refusal to overwrite an index artifact |
| 8 | Internal | An unexpected uncategorized exception |
| 130 | Interrupted | Cooperative `SIGINT` termination |
| 143 | Interrupted | Cooperative `SIGTERM` termination |

### Dependency checkout is rejected

Verify that the override is a Git checkout, its HEAD equals the complete pin,
and `git status --porcelain` is empty. A matching directory name is not proof.

### Persistent index is rejected

Keep the `.sufidx`, `.sufidx.manifest.json`, and `.sufidx.complete` together.
Verify that the reference has the same normalized record order, names, complete
descriptions, and bases used to build the bundle. A copied payload without its
companions is intentionally rejected.

### Output already exists

RaMA-G never overwrites a final artifact. Use a new prefix. Do not delete a
result merely because one component exists; first determine whether its
`.complete` marker exists and preserve incomplete evidence when needed.

### Delta rejects an input path

Use whitespace-free reference and query paths, or omit delta and select another
format.

### SMEM option is rejected

`--smem-min-occurrences` is accepted only when `--seed-mode smem` is explicitly
selected, and its value must be at least 1.

### A compute/resource limit is reached

The failure is deliberate: RaMA-G does not truncate seeds, edges, or DP work and
return an incomplete success. Inspect stderr and the run's `failure.json`, then
change data/parameters or available resources with full awareness that doing so
defines a different run.

### No completion marker exists

Treat the run as unsuccessful or interrupted, even when one or more output
files exist. Inspect `WORK_DIR/runs/*/failure.json` and retain the diagnostic.

## 14. Current limitations

- There is no automatic filtering/report pipeline beyond
  `--selection-mode all|one-to-one`.
- There is no checkpoint, resume, automatic retry, or automatic rerun.
- Persistent indexes are explicit; there is no automatic cache discovery.
- Gzip/BGZF-compatible input is sequential and has no `.gzi` random-access
  route.
- Alignment outputs cannot currently be gzip/BGZF compressed by RaMA-G.
- BAM, CRAM, VCF, coverage tables, and HTML reports are not produced.
- MAPQ is not calibrated and is reported as 255.
- Strict MUM uniqueness is checked within one complete query record, not across
  the collection of all query records.
- MUM and SMEM Sufkit enumeration are serial in the current implementation.
- Reference and query records are materialized in memory; the current pipeline
  has no general bounded spill/checkpoint path for them.
- The current extension implementation uses exact and ungapped fast paths plus
  bounded scalar affine-gap DP; it does not use KSW2.
- RaMA-G makes no unqualified speed, memory, coverage, or accuracy superiority
  claim without an accepted, version-bound comparison.

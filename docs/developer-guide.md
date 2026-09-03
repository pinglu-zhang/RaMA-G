# RaMA-G Developer Guide

This guide describes the architecture and algorithms implemented in RaMA-G
0.1.0 at commit `d2f46338ca89e54bdbeb5f0f2462b0d81635be54`. It is an implementation guide,
not a roadmap or a function-by-function API reference. For command-line use,
see the [User Guide](user-guide.md).

## 1. Architectural goals and invariants

RaMA-G is a deterministic pairwise whole-genome aligner. The design is built
around a small set of invariants that apply across seed modes, thread counts,
and output formats:

- Internal coordinates are unsigned 64-bit, zero-based, and half-open.
- Query coordinates always refer to the original forward query record, even
  when an alignment is on the reverse strand.
- Reference-contig, query-record, and strand boundaries are hard boundaries.
  `N` is a hard break for exact-match enumeration, uniqueness, containment,
  and diagonal seed merging; a later bounded gap alignment may represent a
  short ambiguous region between two legal anchors as mismatches.
- A seed mode must provide its complete promised occurrence semantics. The
  implementation does not silently sample or truncate occurrences.
- Parallel execution may change scheduling, never result ordering or
  tie-breaking.
- One finalized, immutable alignment-record collection feeds every requested
  writer. A writer cannot make a second alignment-selection decision.
- Explicit resource limits fail the run instead of converting an exact
  operation into a partial-success result.
- Production dependency identity is exact, and a successful output set is
  transactionally bound by its manifest and final completion marker.

These rules are checked at subsystem boundaries. Coordinate arithmetic uses
checked conversions and additions; seed and alignment collections receive
canonical sorting; persistent indexes are validated before use; and output
artifacts are validated before the run is declared complete.

## 2. Layers and module responsibilities

The implementation is deliberately split so that input, exact-match
enumeration, alignment construction, record selection, and serialization do
not redefine one another's semantics.

| Layer | Main responsibility | Useful source entry points |
|---|---|---|
| CLI and configuration | Parse the `align` and `index` commands, reject incompatible options, select outputs, print effective configuration, and configure OpenMP | `include/ramag/cli.hpp`, `src/cli.cpp`, `src/main.cpp` |
| Runtime control | Track stages and counters, report RSS, process snapshot requests, and translate signals into cooperative interruption | `include/ramag/runtime.hpp`, `src/runtime.cpp` |
| Input | Detect compression, parse plain or gzip-compatible FASTA, recover headers, normalize bases, and record provenance | `include/ramag/fasta.hpp`, `src/fasta.cpp` |
| Reference index | Build, save, load, validate, and transactionally publish a full Sufkit suffix array and its RaMA-G companions | `include/ramag/reference_index.hpp`, `src/reference_index.cpp` |
| Sufkit adapter | Enforce the pinned API contract and enumerate complete MAM, MEM, MUM, or SMEM seeds | `include/ramag/sufkit_adapter.hpp`, `src/sufkit_adapter.cpp` |
| Alignment core | Normalize seeds, merge exact diagonals, build sparse exact chains, and extend gaps into canonical alignments | `include/ramag/alignment.hpp`, `src/alignment.cpp` |
| Record resolution | Remove exact duplicate records and apply `all` or reciprocal `one-to-one` selection | `src/alignment.cpp` |
| Writers and validators | Serialize and independently validate SAM, PAF, delta, MAF, and chain | `include/ramag/writers.hpp`, `src/writers.cpp` |
| Manifest and pipeline | Orchestrate stages, collect provenance and statistics, manage private files, publish atomically, and write the completion marker | `include/ramag/pipeline.hpp`, `src/pipeline.cpp`, `src/manifest.cpp` |

The CLI owns policy such as option compatibility and requested paths. The
alignment library owns biological and coordinate semantics. The pipeline owns
lifecycle and publication. This separation is important: a formatter failure
can abort publication, but it cannot alter the computed alignments.

## 3. Canonical data model

The public model is defined in `include/ramag/model.hpp`.

### `SequenceRecord`

A sequence record contains a numeric `SequenceId`, the first
whitespace-delimited header token as `name`, the complete header without `>` as
`header`, and normalized `bases`. Numeric IDs are assigned in input order and
are independent of FASTA names. Names are retained for user-facing formats;
IDs provide compact and unambiguous internal joins.

### `Seed`

A seed is an exact match described by reference ID and begin, query ID and
begin, length, and strand. Both begins are contig-local forward coordinates.
For a reverse-strand seed, `query_begin` is still the begin on the original
forward query. Algorithms that require monotone traversal convert it into
oriented query space temporarily.

### `CigarOp` and `AlignmentRecord`

The canonical CIGAR alphabet is extended CIGAR: `=` for a match, `X` for a
mismatch, `I` for query-only bases, and `D` for reference-only bases. Ambiguous
`M` is intentionally excluded, so matching-base counts, edit distance, and the
SAM `MD` tag can be derived without reinterpreting an operation.

An `AlignmentRecord` holds half-open reference and forward-query spans, strand,
canonical CIGAR, signed score, edit distance, and a primary flag. CIGAR
traversal follows the alignment strand, whereas the stored query interval
always satisfies `query_begin <= query_end` in the original query coordinate
system.

Every record is checked for in-range IDs, non-empty and internally consistent
spans, valid operations, checked length accumulation, and agreement between
CIGAR consumption and coordinates. Stable comparison keys include sequence
IDs, coordinates, strand, CIGAR, score, and edit distance; the precise key for
each stage is centralized rather than inherited from callback order.

### `RunStatistics`

`RunStatistics` carries algorithmic counts and timings: normalized input
sizes, MEM/MAM/MUM/SMEM counts, selected and merged seeds, chaining work,
candidate and finalized alignments, extension paths, thread use, route names,
resource-limit observations, and phase time. CLI, index, writer, and artifact
provenance are added by the pipeline when it creates the run manifest.

## 4. Input and normalization architecture

Input format is detected from suffix and bytes before parsing. A `.gz` name
must be a valid gzip stream. A file with gzip magic is accepted even without a
gzip suffix. Other recognized compressed suffixes are rejected explicitly.

Plain FASTA uses the SeqPro/FAI path. RaMA-G preserves or creates only the
read-only-adjacent index behavior supported by that layer and recovers the full
header text needed for provenance and output. Gzip and BGZF-compatible streams
use the linked zlib reader and are consumed sequentially; this route does not
use `.gzi` random access. Concatenated gzip members are accepted, while zlib
errors, truncation, and checksum failures are fatal input errors.

Both routes produce the same normalized records. ASCII bases are uppercased;
`A`, `C`, `G`, and `T` remain unchanged; every other accepted ambiguity is
represented as `N`. The parser rejects duplicate identifiers, missing
identifiers, empty records, and malformed sequence structure. It also records
source path, source size, detected compression, record lengths, and ambiguous
base counts.

Reference and query loading can run concurrently within the single requested
thread budget. The normalized vectors are then immutable for indexing,
enumeration, extension, and serialization.

## 5. Reference-index architecture

### Ephemeral and persistent paths

Without `--reference-index`, the alignment pipeline builds an in-memory index
for that run. The temporary index is not searched for later and is not saved
automatically. The `ramag index` command builds the same kind of index, saves
it, reloads it for self-validation, and publishes a reusable three-file bundle:

- the Sufkit standalone suffix array (`.sufidx`);
- the RaMA-G companion manifest (`.sufidx.manifest.json`);
- the completion marker (`.sufidx.complete`).

An align run with `--reference-index` requires all three files and loads rather
than rebuilds the index.

### Sufkit construction contract

The production adapter requires Sufkit 0.3.0 at exact commit
`bdb67c6de5daddd8a005640de73d96549d2575f4`. It requests the Fast resource
profile, a standalone full suffix array with sampling rate 1, and the ISA, LCP,
and suffix-link/query capabilities required by every seed mode. Learned-index
construction is disabled.

For a multi-threaded reference of at least 64 MiB, the current policy chooses
the shared-memory CaPS constructor. Smaller inputs, or a one-thread build, use
the bundled divsufsort constructor. The resulting complete query structure and
seed semantics are the same; the choice is recorded as provenance.

### Persistent validation

Sufkit `Save()` and `Load()` provide the binary format, section CRC checks,
suffix-array permutation validation, and internal semantic checks. RaMA-G adds
binding checks for the companion manifest and completion marker, exact Sufkit
version/commit and capabilities, index kind and sampling rate, and the
reference catalog.

The reference catalog includes record count, order, names, descriptions,
lengths, and normalized ambiguous-base statistics. A SHA-256 digest over the
normalized reference establishes algorithm-level identity. Raw source path,
size, and related file facts remain separate provenance. Thus an equivalent
reference representation can reuse an index only when the normalized content
and catalog satisfy the full contract.

The suffix array is not a replacement for reference sequences. RaMA-G keeps
the normalized reference bases available for maximality checks, chaining-gap
extension, CIGAR validation, and sequence-bearing output such as SAM and MAF.

## 6. Seed algorithms

Let a seed identify equal substrings of the reference and an oriented query.
The principal seed concepts are:

- A **MEM** is exact and cannot be extended by one base on either side while
  preserving equality for that occurrence.
- A **MAM** is a maximal exact match whose matched reference string is unique
  in the combined reference; query repetition is allowed.
- A strict **MUM** is maximal and unique in both the combined reference and
  the current complete query record. Overlapping query occurrences count
  against uniqueness.
- An **`(l,c)`-SMEM** is a query interval of length at least `l` with at least
  `c` reference occurrences, and is not properly contained in another query
  interval satisfying those conditions. Every retained interval is expanded
  into every reference occurrence.

RaMA-G enumerates each query in forward and reverse-complement orientation.
Callbacks in reverse-oriented coordinates are checked for range and converted
back to original-forward query coordinates before a `Seed` is stored.

### `maxmatch`

The adapter calls Sufkit MEM enumeration over each complete oriented query
record and preserves the complete occurrence set. It then canonicalizes,
sorts, and removes only exact duplicate seed tuples.

### `mumreference`

Reference-unique matching uses stable 4 MiB non-overlapping query tiles to
provide parallel work without changing global semantics. Tile MAM callbacks
are rechecked for maximality against the complete query record. Exact matches
that cross an internal tile boundary are recovered by targeted whole-record
MEM enumeration, followed by global reference-uniqueness verification.

Tasks occupy deterministic slots for every query, orientation, tile, and
boundary recovery region. Results are merged in stable task order and then
canonically sorted. The adapter has hard limits for boundary MEM callbacks and
its query-proportional workspace; exceeding them fails rather than drops
occurrences.

### `mum`

Strict MUM enumeration calls Sufkit on one complete oriented query record at a
time. No tiles are used, because tile-scoped uniqueness would not be
query-record uniqueness. The current Sufkit enumeration layer for this mode is
serial. The rest of the alignment pipeline may still use its configured
parallel stages.

### `smem`

SMEM enumeration likewise operates serially on each complete oriented query
record. Each callback is checked for a nonzero and internally consistent
reference-occurrence count meeting `c`; the interval is then expanded without
sampling into contig-local `Seed` coordinates. Interval statistics and expanded
coordinate-seed statistics remain distinct.

### `fast`

The default route first obtains the reference-unique MAM skeleton described
above. It also enumerates whole-query MEMs and retains a MEM only when its
forward-query interval is not already covered by the skeleton. This preserves
an occurrence-complete source for uncovered regions while allowing the MAM
skeleton to reduce redundant downstream work.

### Independent oracle route

When the pinned adapters are deliberately disabled, small exhaustive tests can
use an independent seed oracle. It directly checks exact matches, maximality,
reference/query uniqueness, and SMEM containment. This is a development and
correctness oracle only; it is not a production substitute or benchmark
configuration.

## 7. Seed normalization and diagonal merge

Seed callbacks are not trusted to arrive in a stable order. RaMA-G validates
coordinates and lengths, then sorts by reference ID and position, query ID,
strand, query position in the appropriate orientation, length, and stable
tie-break fields. Exact duplicate tuples are removed.

For chaining, a reverse-strand seed beginning at forward coordinate `q` with
length `k` on a query of length `Q` has oriented begin:

```text
q_oriented = Q - (q + k)
```

The signed diagonal is `reference_begin - q_oriented`. Overlapping or adjacent
seeds are merged only when reference contig, query contig, strand, and this
exact diagonal agree. The merged interval is revalidated against both
sequences. Contig, strand, and `N` boundaries therefore cannot be bridged by
the merge.

## 8. Sparse exact chaining

Seeds are partitioned by reference contig, query contig, and strand. Inside a
group, coordinates are monotone in reference and oriented-query space.

For two ordered seeds, define the intervening nonnegative gaps as `r_gap` and
`q_gap`. A legal edge requires:

```text
r_gap <= min(max_gap, break_length)
q_gap <= min(max_gap, break_length)
abs(r_gap - q_gap) <= diag_diff + diag_factor * max(r_gap, q_gap)
(r_gap + 1) * (q_gap + 1) <= max_dp_cells
    when both gaps are nonzero and unequal
```

Zero-length gaps use direct indel paths. Equal-length gaps normally use the
exact or ungapped path and are therefore not rejected by this chaining test;
the limited equal-gap affine-realignment case rechecks the cell bound before
allocating its matrix.

The production chainer discovers spatially plausible predecessor pairs using
endpoint buckets, but every candidate is passed through the exact `CanChain`
predicate. Legal edges form connected components via union-find. Within each
component, deterministic dynamic programming ranks paths, recovers the best
chain, removes its seeds, and repeats until the component is exhausted. The
chain score and all ties use stable coordinate/input-order rules.

The implementation accounts for candidate-pair checks, legal edges, DP edge
relaxations, and a deterministic working-set estimate. Default hard limits are
500,000,000 candidate checks, 100,000,000 legal edges, 1,000,000,000
relaxations, and 6 GiB of chainer-owned working data. Crossing a limit raises a
computation/resource error; it does not invoke an approximate fallback.

A quadratic exact chainer is retained behind test-only access. Differential
tests compare the sparse route against it on small inputs. Production cannot
select it as a fallback.

Independent contig-pair/strand groups run in parallel. Each task writes to a
stable result slot, periodically checks interruption at safe boundaries, and
is merged only after all preceding deterministic ordering rules are applied.
Nested OpenMP teams are avoided.

## 9. Extension and canonical alignment construction

Each recovered chain is converted into one candidate alignment. Exact seed
segments become `=` operations. Intervening regions use three paths:

1. an exact-string fast path when the two gaps are identical;
2. an ungapped path for suitable equal-length gaps, emitting `=` and `X`;
3. bounded scalar affine-gap dynamic programming for the remaining gaps.

The current scoring scheme is match `+2`, mismatch `-4`, gap open `-4`, and
gap extension `-2`. The public tuning interface does not currently expose
these four scoring constants. A two-dimensional unequal gap requiring more
than `--max-dp-cells` is not admitted by chaining. The limited equal-gap
realignment path checks the same bound in the extension layer, so no scalar DP
matrix exceeds the configured cell limit.

The candidate starts at the first seed and ends at the last seed in its chain;
the current implementation does not extend beyond those chain endpoints.
`break-length` is an inter-seed chain split threshold, not an endpoint-extension
radius. Reverse-strand query access uses reverse-complement sequence material
while the resulting interval is converted back to forward coordinates.

Adjacent identical CIGAR operations are coalesced. The constructor then checks
reference/query consumption, span closure, score, edit distance, and sequence
agreement for every `=`/`X` column. KSW2 and other external extension engines
are not active in this implementation.

## 10. Conflict resolution and selection

Candidate records first receive canonical content ordering and exact-record
deduplication.

In `all` mode, all remaining distinct records are retained. For each query
contig, one deterministic best record is marked primary and other records are
secondary/supplementary. Primary ranking uses score, aligned query length,
identity, and stable reference/query coordinate and strand tie-breaks.

In reciprocal `one-to-one` mode, RaMA-G constructs elementary intervals from
all alignment endpoints independently on every reference and query sequence.
For each covered interval, a priority queue selects the best supporting record
by matching bases, identity, covered span on that side, and stable canonical
alignment order. It forms a set of records selected somewhere on the reference
and a separate set selected somewhere on the query, then retains their
intersection. Every retained record is marked primary.

Conceptually:

```text
reference_selected = best records covering each reference elementary interval
query_selected     = best records covering each query elementary interval
retained           = reference_selected intersect query_selected
```

This is reciprocal interval support, not a global collinearity or syntenic-LIS
optimizer. It is also not claimed to be algorithmically equivalent to
MUMmer4's `delta-filter -1`.

## 11. Writer architecture

Writers receive the same finalized, read-only `AlignmentRecord` vector and the
same normalized sequences. They may convert coordinates and CIGAR syntax, but
they may not filter or reprioritize records.

| Writer | Principal projection |
|---|---|
| SAM | Reference position becomes 1-based `POS`; query orientation is represented by flag `0x10`; non-primary records use `0x800`; terminal unaligned query is hard-clipped; `NM`, `MD`, and `AS` are derived from the canonical record; `MAPQ` is 255 |
| PAF | Both spans remain zero-based and half-open; strand is explicit; matching-base count and block length come from CIGAR; `cg`, `NM`, `AS`, and `tp:P/S` tags preserve record details; `MAPQ` is 255 |
| delta | Header paths identify reference and query; coordinates are projected to 1-based inclusive NUCMER conventions, including decreasing query coordinates on the reverse strand; indels are encoded from canonical CIGAR |
| MAF | Each record becomes one pairwise `a` block; gapped reference/query rows are generated directly from `=/X/I/D`; the reference row is `+`; reverse query uses start `query_size - query_end` |
| UCSC chain | Reference is target on `+`; query carries the actual strand; `=`/`X` form aligned blocks, `D` becomes target gap, and `I` becomes query gap; terminal gaps are trimmed from the represented span; IDs follow stable record order |

SAM hard clipping and strand handling preserve the original query length and
coordinates. MAF `size` fields count nongap symbols. Delta serialization
checks that every signed indel displacement fits its format. Chain
serialization normalizes consecutive gap operations between aligned blocks; a
record with no aligned block is an output error rather than a silently omitted
record. Chain scores are `max(0, alignment.score)`.

Each format has an independent structural validator. Because requested files
are one transaction, failure of any writer or validator invalidates the whole
requested output set.

## 12. Pipeline, publication, and diagnostics

The alignment pipeline advances through these observable stages:

1. configuration and private run-directory creation;
2. input;
3. index build or index load;
4. seed enumeration;
5. seed merge;
6. chaining;
7. extension;
8. conflict resolution;
9. one writer stage per requested format;
10. validation;
11. publication;
12. complete.

Each run receives a unique run ID and a private directory beneath the selected
work directory. Writers create private files first. After format validation,
the pipeline publishes output files to their sibling final paths with
no-replace semantics. It revalidates published files, publishes the manifest,
and writes the `.complete` marker last.

If any write, validation, or publication step fails, files published by this
run are rolled back and no completion marker is produced. Existing destination
files are not overwritten. A structured `failure.json` remains in the private
run directory with status, exception category, stage, time, and message.

`SIGINT` and `SIGTERM` set asynchronous flags only. Normal control flow checks
those flags during input, Sufkit callbacks, chaining/extension task boundaries,
writers, validation, and publication, performs cleanup, records an interrupted
diagnostic, and returns 130 or 143. Sufkit, libdivsufsort, and CaPS constructor
regions without cancellation callbacks can react only after library control
returns.

## 13. Concurrency and memory model

`--threads` is one shared budget. RaMA-G configures the OpenMP runtime and uses
that budget sequentially across stages rather than creating independent nested
teams. Parallel regions include eligible input work, large-reference CaPS
construction, tiled reference-MAM work, independent chaining groups, and
independent chain-extension tasks.

MUM and SMEM Sufkit callback enumeration are currently serial and operate on a
whole query record. Small references use divsufsort even when more threads are
requested; the manifest records requested, scheduled, and observed worker
facts rather than implying every stage used all threads.

Normalized reference and query sequences and the reference suffix array are
shared and immutable after construction. Tasks own local callback buffers,
edges, DP workspaces, alignments, and stable output slots. Merges occur in
canonical order, so scheduling does not leak into results.

The current memory model materializes all normalized reference and query
records, the full standalone suffix array and auxiliaries, selected seed
vectors, and bounded chaining/extension workspaces. Hard limits cover MAM
workspace, MAM boundary callbacks, chaining work, and per-gap DP cells. There
is no disk-spill, checkpoint, or partial-result mode.

## 14. Provenance and observability

The run manifest uses conceptual schema `ramag.run-manifest.v1`, version 1. It
binds the result set to:

- RaMA-G version and binary identity, build type, compiler, and OpenMP facts;
- exact Sufkit and SeqPro commits and the selected index backend/profile;
- command configuration and requested outputs;
- input paths, detected encodings, normalized record statistics, and available
  content identities;
- persistent-index identity and whether the action was `built` or `loaded`;
- actual seed, query-parallel, chaining, and extension routes;
- requested and observed thread counts;
- seed, chain, alignment, selection, and resource-limit statistics;
- phase timings and output artifact paths, sizes, validation states, and
  publication states.

The completion marker is the final binding for the requested output set and
manifest. Neither process disappearance nor the presence of an individual
format proves successful completion.

Progress reporting is a view of the same run state, not a second source of
truth. A session stores the run ID, stage, elapsed times, known completed/total
units, counters, active threads, and Linux current RSS from `/proc/self/status`.
Periodic output follows `auto|on|off`; stage changes are immediate when progress
is enabled. `SIGUSR1` sets an async-signal-safe snapshot flag, and a normal
watcher/control boundary prints the requested state without terminating the
run.

## 15. Current correctness boundary

The internal validation baseline used before preparing the public source tree
contains 14 CTest tests covering core units, required OpenMP configuration,
exhaustive seed oracles, Sufkit adapter differential behavior,
sparse-versus-quadratic chaining checks, clean and deliberately invalid
dependency gates, end-to-end integration, progress/signals, and external SAM
validation when the validator is available. Those internal test sources are
not distributed in the public tree. CMake detects their absence and skips test
targets without affecting the `ramag` library or executable build. This
retained validation evidence supports the invariants described here; it is not
a universal accuracy or performance claim.

The implementation boundary is currently:

- one reference multi-FASTA and one query multi-FASTA per align run;
- in-memory normalized genomes and a full standalone reference suffix array;
- the five documented exact-seed modes;
- bounded scalar affine-gap extension;
- `all` or reciprocal interval-based `one-to-one` record selection;
- SAM, PAF, delta, pairwise MAF, and UCSC chain serialization;
- transactional output, provenance, progress, and cooperative interruption;
- no checkpoint/resume, automatic index discovery, calibrated MAPQ, automatic
  coverage/comparison report, global syntenic filtering, BAM/CRAM, or VCF.

Correctness and performance claims beyond these implemented and tested
boundaries require a separately specified and accepted evaluation.

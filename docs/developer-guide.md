# RaMA-G developer guide

## Architecture and invariants

RaMA-G separates normalized sequence input, Sufkit seed enumeration, the standalone
pairwise KSW2 core and output serialization. Internal coordinates are 64-bit,
zero-based and half-open. Query coordinates refer to the original forward record
on either strand; oriented query begin is `query_length - query_end` for minus
records. Contig boundaries and `N` break seeds. Occurrences are not sampled away.

| Layer | Responsibility |
|---|---|
| CLI/runtime | Parameters, launch CPU affinity, thread budget and signal flags |
| FASTA | kseq plus a strict streaming validation wrapper over zlib |
| Index/adapter | Full Sufkit SA, persistence, capabilities and seed conversion |
| Pairwise core | Seed grouping, clustering, KSW2 extension/linking, bilateral selection |
| Writers/pipeline | Final immutable records, validation and no-replace publication |
| Logging | Optional CLI-owned spdlog sinks and retained run statistics |

`SequenceRecord` owns name, full header and normalized bases. `Seed` carries
contig IDs, strand and exact coordinates. `AlignmentRecord` owns coordinates,
score, edit distance, status and canonical `=/X/I/D` CIGAR. Ambiguous `M` is resolved
against bases before serialization; `N/N` is not a canonical match. Overflow,
consumption, bounds and nonzero operations are checked.

## Input and dependencies

The validated read callback checks raw decompressed bytes before kseq can skip
invalid FASTA or accept FASTQ. It tracks line state and a bounded read-ahead queue
of complete headers, preserving descriptions independently of kseq tokenization.
Errors from gzip CRC/truncation and normal EOF are distinguished. File and parser
storage are RAII-owned. Plain and gzip input take the same normalization route;
no SeqPro/FAI path remains. The frozen parser's `int` length boundary is checked
explicitly for an individual record. Record collections and resulting alignment
coordinates remain 64-bit.

Reference and query reading retain stable result/error slots and may use two
OpenMP workers within the requested budget. No exception crosses an OpenMP
structured block. The full normalized inputs remain materialized.

Sufkit stays at `028075e6f2d622fcbdcf76b153bbde0f069ca64f`. Its exact clean checkout
is checked at configuration; its optional SeqPro integration remains disabled.
KSW2, kseq and spdlog are vendored; spdlog uses bundled fmt. kseq/spdlog snapshots
come from the frozen source revision recorded in the root third-party notices.

## Index and seeds

Fast standalone SA retains SA, ISA, LCP and suffix links at sampling rate one.
New indexes use byte-coded LCP, preserving every decoded value. CaPS is selected
for references at or above the existing 64 MiB threshold when parallel building
is applicable; smaller construction uses bundled divsufsort. The executable
captures launch affinity before libgomp initialization and restores only that
authorized set. Library calls do not widen their caller's affinity.

Ephemeral construction does not serialize. Saving calls `Save()` on the existing
index, self-validates the temporary file using `Load()`, then publishes it without
replacement. Loading validates the Sufkit format, CRC, SA structure, sampling and
capabilities. Contig catalog and Sufkit reference fingerprint are compared against
the supplied reference. No SHA-256 or sidecar identity is used. Raw and byte-coded
files retain their stored representation; loading never rebuilds implicitly.

RaMA-G supplies the caller's thread budget to `SuffixArray::Load(path, options)`.
The original Sufkit `Load(path)` overload remains single-threaded by default.
On Linux the loader checks that the requested budget fits the allowed CPU set.
Bounded section readers accumulate CRC while decoding; all sections, including
unconsumed ones, must pass CRC before an index is returned. Linux reads refer to
the same opened inode. Other platforms retain the full-validation fallback.

For a full ISA, range checks and the relation `SA[ISA[p]] == p` establish a
bijection between positions and ranks, replacing the redundant permutation
bitmap. Consistency checks use deterministic parallel chunks and report the
smallest failing position. The no-ISA compatibility path retains permutation
validation. Fast prefix-directory reconstruction receives the loading thread
budget and preserves the directory layout. Byte-coded LCP stays compressed.

Reference matching uses read-only sequence views instead of constructing another
owning `GenomeReference`. Metadata, ambiguity and reference fingerprint checks
remain active. Loading statistics and callbacks report section phases, logical
bytes, CRC CPU time and elapsed load time. RaMA-G logs Sufkit load and reference
validation separately; nested measurements must not be added to their totals.

Both query orientations are searched. MEMs are bilateral maximal exact matches;
MAMs additionally require reference uniqueness. Strict MUMs require uniqueness in
the reference collection and current query orientation/record. Generalized SMEMs
retain eligible query intervals not properly contained in another eligible
interval and expand all reference occurrences.

`maxmatch` enumerates all MEM occurrences. `mumreference` uses stable 4 MiB tiles,
whole-record maximality checks, boundary MEM recovery and reference-uniqueness
verification. `mum` and `smem` use whole-query-record serial enumeration.
`fast` uses the existing MAM skeleton followed by whole-query MEM filtering.
Task results merge in deterministic order before final sorting/deduplication.
Small independent oracles remain internal tests, not production fallback routes.

## Pairwise algorithm and memory

The core groups and filters anchors by reference/query contig and direction,
merges according to the existing diagonal rules, forms clusters, recovers the
best supporting chain and applies minimum support. It materializes alignment
anchors, extends components and attempts connections under the current windows,
gap limits, quality and closure rules.

KSW2 uses the existing scaled HOXD70 matrix and gap convention. Exact and simple
gap paths precede global KSW2. Global alignment uses the existing bandwidth
optimality check, widening or full calculation when necessary. Endpoint extension
and per-worker scratch reuse remain unchanged. Scores used for DP, selection and
serialized records have distinct purposes and are not interchangeable.

Reference-side and query-side DP selection retain their current scores and
tie-breaking. `one-to-one` intersects the two selected candidate sets; `all`
exposes valid preselection records. No former legacy reciprocal interval selector
is applied afterward. Optional residual recovery clips eligible unused portions
of rejected packed-CIGAR candidates. Optional guarded gap fill uses reliable
original flanks, fixed nonrecursive candidates, exact or current KSW2 alignment
and bilateral occupancy checks. Their existing build defaults remain unchanged.

The index is destroyed after seed enumeration. The CLI transfers Seed ownership
into the core; the const-span library entry copies once. Task seed chunks are
released during final collection. Only output candidates expand packed CIGAR
into final records. Capacity estimates and RSS observations describe individual
stages, not additive peaks. Resource failure is explicit; there is no spill or
automatic checkpoint path.

## Output and logs

All writers consume the same immutable finalized record set. SAM performs reverse
flag/hard-clip and MD/NM conversion; PAF carries original-forward intervals and
extended CIGAR. Delta encodes one-based inclusive endpoints and indel displacements.
MAF reconstructs paired gapped rows, using `query_size-query_end` on minus strands.
Chain converts aligned runs to blocks and I/D to query/target gaps, removing
unrepresentable terminal gaps. Writer validation checks coordinate/consumption
closure and format constraints before publication.

Output publication uses temporary siblings and no-replace hard links. A scope
guard rolls back names still referring to this run's inodes on handled failure.
The independent saved-index transaction survives later alignment failure.
No manifest, SHA-256, JSON failure diagnostic or completion marker is generated.
There is no multi-file crash-atomic guarantee: callers must retain process exit
status and validate expected outputs after a forced termination.

`RunLogger` owns unregistered spdlog file and stderr loggers. The CLI creates it
in a fresh private run directory and passes it explicitly to high-level pipelines.
Default library calls create no logger and do not change host global logging or
install signal handlers. The standalone installed pairwise target has no logging
dependency. Effective config may still be printed as JSON to stdout on request.

Stage changes, timing summaries, input/index metadata, actual worker statistics,
memory observations, pairwise counters and artifacts are written to `run.log`.
Periodic terminal progress obeys auto/on/off; explicit snapshots remain available.
The progress watcher has bounded waits and reports failures back to normal flow.
Only flag updates occur in signal handlers. Sufkit's noncancellable internal work
can delay response; normal boundaries check interruption and flush diagnostics.

## Library use and verification

```bash
cmake --install build-release --prefix "$PWD/install"
cmake -S examples/pairwise -B external-build -DCMAKE_PREFIX_PATH="$PWD/install"
cmake --build external-build -j 4
```

Consumers use `find_package(RaMAGPairwise CONFIG REQUIRED)` and link
`RaMAG::pairwise`. Sequence/index inputs are read-only shared data, workspaces
belong to each call and returned records own their storage. Existing high-level
alignment functions remain. Removed manifest, SHA-256 and SeqPro headers are a
source-compatibility change; use `ReadFasta`, returned `RunStatistics` and optional
logging instead.

Internal tests cover independent seed/DP oracles, adapter differential behavior,
ownership paths, public entry-point equivalence, OpenMP and dependency contracts,
all five serializers, input/persistence/logging failure paths and signal rollback.
Release/Werror and ASan/UBSan are validated separately. A non-PIE sanitizer build
is used on the server to avoid the previously observed startup/ASLR issue; this
does not imply sanitizer coverage of an uninstrumented third-party binary.
Public release builds remain supported without internal tests or planning files.
Runtime performance and biological accuracy require separately frozen evidence.

The pending 0.1.1 version is defined by CMake `PROJECT_VERSION`; it supplies
`RAMAG_VERSION` for CLI/SAM identification and the installed package version.
Version preparation changes only that definition and documentation. Historical
Release/Werror and sanitizer results cover the previously accepted loading
implementation, not a fresh execution of the pending version. See the
[changelog](../CHANGELOG.md) for the release boundary.

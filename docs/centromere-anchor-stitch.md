# Anchored path base alignment

Stitching consumes the current primary anchors and accepted copy branches.
It does not rerun direction identification, anchor search or copy acceptance.
Each primary `region_id` and each copy `branch_id` is an independent path.
Different paths may reuse physical intervals. Ambiguous candidates are not
input paths. Only the intervals between two confirmed consecutive anchors are
aligned. An isolated anchor is exported unchanged; regions without anchors,
sequence tails, and unknown directions remain unresolved.

## Build and command line

The normal Linux build runs stitching as part of centromere alignment.

```sh
ramag align -c --reference reference.fa --query query.fa \
  --output-prefix result --threads 1
```

The default standard outputs are SAM, PAF and MAF. Use `--formats LIST` to
select a subset, or provide explicit paths with repeated `--output`. Each
selection keeps the segment, gap, direction, anchor and review TSVs.
`--print-effective-config` reports the scoring model, budgets and selected
formats. The [format contract](centromere-alignment-formats.md) describes
coordinate and copy-path conventions.

## Exact objective and engine selection

The default `centrolign-3piece` objective is

```text
20 * matches - 80 * mismatches
  - sum over same-direction gap runs min(60+30*l, 800+5*l, 2500+l)
```

A gap of length one pays open plus extension. Adjacent I and D runs are charged
separately; adjacent operations of the same type are coalesced. Only canonical
ACGT equality is an exact match: N/N is a mismatch. The pinned implementation
is Centrolign commit `a604cbdcff2e22e643fc9e17d4679296ee5e42fe`.
Its private minimal sources retain MIT attribution and source manifests.

Empty/empty gaps require no alignment. A single nonempty side produces its
unique I or D. Both-nonempty gaps with `(m+1)*(n+1) <= 4,000,000` use full
`po_poa<3>` dynamic programming. Larger gaps use unpruned `wfa_po_poa<3>` on
linear graphs with fixed endpoints. No greedy, pruned, or single-long-deletion
fallback is used. An independent CIGAR scorer checks coordinate consumption,
characters and score. Different optimal CIGARs are allowed when scores tie.

The library also exposes the two-piece scoring control for matched-boundary
comparisons. It uses match 0, mismatch -3 and gap penalty `min(4+2*l,12+l)`.

## Coordinates and output contracts

All positions are original FASTA coordinates, zero-based and half-open.
CIGAR traverses increasing reference coordinates. I consumes query and D
consumes reference. For negative paths the query view is reverse complemented;
its original interval coordinates remain increasing in the tables. An emitted
`=` is verified canonical equality and `X` is a mismatch. No CIGAR is invented
across an unresolved gap or between different paths.

| Output suffix | Content |
| --- | --- |
| `.alignment-segments.tsv` | Contiguous path segments, role, path/region/branch IDs, coordinates, strand, CIGAR, score/model, anchor/gap IDs |
| `.alignment-gaps.tsv` | Every internal gap, flanking anchor IDs, coordinates, complete/unresolved status, reason, engine, CIGAR, score, work counts, measured kernel memory and time |
| `.alignment-summary.tsv` | Internal-gap completion, paths, anchors, engine calls, union coverage and timing |
| `.paf`, `.sam`, `.maf` | Standard per-segment alignments with original coordinates, direction and path provenance; all three by default when stitching is enabled |

`primary` and `copy` roles are explicit. Consumers must not collapse shared
coordinates into a single one-to-one chain or discard copy paths as secondary
SAM/PAF alignments. In the standard SAM export, overlap can require a secondary
flag even for an accepted copy: retain such records and inspect `pr` and `pi`,
as described in [the format contract](centromere-alignment-formats.md).
Coverage uses physical coordinate unions rather than sums
across repeated paths. Paired coverage differs from alignment consumption:
an I/D consumes one axis but introduces no paired bases.

`status=complete` means all eligible internal gaps completed. It does not
mean a whole-input global alignment was performed. `status=partial` preserves
all exact anchors and successful gap segments while listing unresolved gaps.
The CLI may publish such valid partial output and return success; consumers
must inspect the summary and gap status. The `.directed-run.tsv` explicitly distinguishes anchor and stitch status;
its legacy WFA count is labelled as anchor-stage only when stitching is enabled.
New kernel calls are recorded separately in the alignment summary.

## Resource and failure semantics

The default limits are 512 MiB live kernel working allocation and 120 seconds
per gap, with 1,800 seconds for the whole stitch stage. Gaps run serially and
kernels use one thread. Allocations and cooperative interrupt checks cover
the private C++ kernels, not only old WFA2 C allocations. Reported kernel work
memory is separate from the process-tree peak RSS measured by the runner.

Memory/time exhaustion produces an unresolved gap and splits the output path.
An invalid input geometry, invalid CIGAR, or inconsistent score is a program
error. User interruption is propagated; the library observation callback can
retain already computed snapshots. The CLI also retains frozen-anchor diagnostic files and periodic stitch snapshots
in its run directory. An aborted snapshot is explicitly labelled and includes
the failure reason; pending gaps are not completed predictions. These files
are not formal published output. CLI publication remains atomic and rejects
existing files instead of overwriting them.

## Public API and validation

Include `ramag/centromere_stitch.hpp` and link `RaMAG::centromere`:
`AlignCentromereAnchoredPaths(reference, query, primary, optional_copy, options)`.
Inputs are immutable. `CentromereStitchOptions` exposes the scoring model,
cell threshold, memory/time budgets, interruption, progress and observation
callbacks. The output includes path segments, per-gap status and statistics.
The private kernel dependency is embedded in the installed library.

Correctness checks compare independent seven-state dynamic programming with
production DP/WFA and independently replay/rescore each CIGAR. Path tests
cover negative coordinates, I/D direction, gaps, copied intervals, failures
and reversal under input exchange. Regression evaluation separates old anchor
TP/FP, newly paired primary/copy gap bases, and deduplicated union metrics.
Truth is used only after alignment output is frozen. The existing 96 cases
are regression data, not a newly held-out test set. Any decrease is recorded;
this feature remains explicitly enabled and does not become the default.

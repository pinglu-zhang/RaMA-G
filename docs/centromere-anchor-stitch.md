# Anchored-path base alignment

`AlignCentromereAnchoredPaths` is a pure consumer of primary and accepted copy
anchors. It aligns only consecutive-anchor internal gaps within each path. It
does not predict direction, extend unanchored tails or select ambiguous copies.
The complete `align -c` command orchestrates direction, anchors and reviews around
this operation and exports SAM, PAF, MAF and diagnostic TSVs.

## Scoring and kernels

The default three-piece score is

```text
20 * matches - 80 * mismatches
  - sum(min(60 + 30*l, 800 + 5*l, 2500 + l))
```

Each maximal same-direction gap is charged independently. A one-base gap costs
open plus extension; adjacent insertion/deletion runs open separate gaps. N/N
is not an exact match. This objective is distinct from genome HOXD70 alignment
and the low-level source-derived RaMA pair aligner.

An empty gap needs no work. A one-sided gap generates its actual I or D. For two
nonempty sides, at most 4,000,000 grid cells use complete three-piece DP; larger
fixed-end problems use unpruned linear-graph WFA. An independent CIGAR replay
checks score, characters, boundaries and both consumption lengths.

Default kernel working-allocation budget is 512 MiB per gap, with 120 seconds per
gap and 1800 seconds for stitching. These are not whole-process RSS limits.
Resource exhaustion leaves an unresolved gap; invalid CIGAR or inconsistent
score is an error. Caller cancellation propagates rather than pretending success.

## Paths and outputs

Main paths are grouped by directed region, copy paths by branch. Each path must
have one direction and ordered nonoverlapping anchors. Different paths may share
physical intervals. Negative CIGAR traverses reference increasingly and query
decreasingly, while published coordinates remain on the original FASTA.

Only genuinely continuous pieces are joined. A failed gap splits the result;
there is no fabricated I+D replacement. A singleton path outputs its exact anchor.
Completed gaps do not imply whole-input global alignment: tails, zero-anchor
regions, unknown direction and ambiguous copies can remain unaligned.

`.alignment-segments.tsv` is the path-segment representation,
`.alignment-gaps.tsv` records gap status/engine/reason, and
`.alignment-summary.tsv` records counts and completion scope. Standard formats
are equivalent views with path provenance; see [format contracts](centromere-alignment-formats.md).

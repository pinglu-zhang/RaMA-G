# Alignment workflows

## Whole genomes

The program normalizes FASTA input, constructs or loads a pinned Sufkit reference
index, enumerates exact seeds on both query orientations, and forms candidate
chains. KSW2 extends candidates with the genome score model. Reference-side and
query-side DP selections are computed independently; ordinary `align` retains
their intersection and links compatible neighboring records.

Linking retains existing pairs. Only same-contig, same-strand mutual neighbors
are eligible. Supported bounded gaps can add exact sequence, one-sided gaps or
a fixed-end KSW2 alignment. Invalid geometry, insufficient support, N, exhausted
work or failed replay leaves the original split. See [main-path linking](mainpath-linking.md).

The genome scoring model is HOXD70 divided by 10 with affine gap cost
`40 + 3 * length`. It is not the centromere stitching score. Seed modes and
threading are described in the [user guide](user-guide.md).

`batch` reuses a reference/index and processes query assemblies serially. Its
candidate-output policy is deliberately retained. Library clients can explicitly
request candidates, an unlinked bilateral selection, or a linked selection.
Retired coverage and optional recovery experiments are not product capabilities.

## Centromeres

The default Linux `align -c` workflow performs:

1. Direction evidence and V9 planning.
2. RaMA recursive anchors in the planned directions.
3. Accepted copy branches and explicit ambiguity reporting.
4. Root-endpoint review and internal-gap stitching.
5. Bounded block, cross-region and joint-path reviews, followed by supported
   unplanned-core recovery and final path output.

Current direction protections remain enabled; directory organization does not
change their order or acceptance rules. Initial/final anchor updates, trusted
internal support, copy footprints and fixed-gap replay remain coordinated.

Stitching uses `20` for an exact canonical match, `-80` for a mismatch, and the
minimum of `60 + 30*l`, `800 + 5*l`, and `2500 + l` for each same-direction gap.
The DP and unpruned WFA kernels use this same three-piece score. Unknown bases do
not become exact matches. The original RaMA backend's pair-alignment API has a
different role and score; see [upstream integration](centromere-rama-integration.md).

Only supported paths and internal gaps are exported. Unknown directions,
ambiguous copies, unanchored tails and budget-limited gaps may remain unresolved.
Accepted copy paths can share a physical interval. They are not flattened into
one globally monotone correspondence.

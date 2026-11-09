# Copy-aware sparse anchors (development)

The copy layer runs after primary anchor discovery in `ramag align -c`. It
uses those anchors as fixed evidence and records each accepted copy branch
and ambiguity group separately. The stitching stage then fills supported gaps
in accepted branches and primary paths.

## Relations and coordinates

All coordinates refer to the original input sequences and are zero-based,
half-open. Negative anchors relate a query interval to the reverse complement of
the reference interval. Exact anchors contain only A/C/G/T. No CIGAR is emitted
for gaps between anchors. Different accepted branches may reuse an interval;
within a branch the aligned intervals have one orientation and ordered,
nonoverlapping coordinates.

The original `initial-anchors.tsv` and `final-anchors.tsv` remain primary-only.
Additional outputs have explicit roles:

| Suffix | Role |
| --- | --- |
| `.copy-anchors.tsv` | Accepted supplementary exact coordinate relations |
| `.copy-ambiguities.tsv` | Unresolved alternatives; not confirmed alignments |
| `.copy-self-repeats.tsv` | Physical sequence-internal repeat evidence |
| `.copy-tasks.tsv` | Search bounds, candidate counts, timing and task status |
| `.copy-summary.tsv` | Copy-layer totals and resource truncation |

A self-repeat is a sequence relationship, not proof of a historical duplication.
Only original primary anchors can connect self-repeat evidence across the two
sequences. Newly inferred branches cannot seed further inference. The copy layer
can therefore preserve ambiguous ancestral copies instead of choosing a
reference coordinate from identical strings alone.

## Resource and scoring contracts

The public `CentromereCopyAnchorOptions` defaults are rare level 100, 3,000
deduplicated candidates per local problem, and 100,000 candidates per sequence
pair. Resource-exhausted tasks are reported explicitly; a truncated task cannot
publish an accepted branch. The CLI uses these defaults.

Copy-chain match weights are length divided by the square root of the product
of the two occurrence counts. Gap discrepancy is a real penalty; negative
incremental scores are not clamped to the primary algorithm's positive floor.
Competition ties and insufficient independent context are ambiguity reasons,
not permission to choose a unique counterpart. A fixed primary direction does
not identify a unique repeat copy.

## Evaluation

The standalone copy benchmark keeps three results: primary, newly accepted
copy relations excluding duplicate primary pairs, and the deduplicated union.
The strict ancestral key is `(source_id, source_offset, copy_at_M)`. Post-split
duplication may yield valid one-to-many relations; different pre-split copies
remain distinct. Query, reference, and ancestral-site recovery have separate
denominators. A Cartesian product of all allowed copy pairs is not a recall
denominator.

Ambiguity coverage, candidate correctness, truncation, and failures accompany
accepted-anchor scores. Exact sequence identity alone is not correctness.
These outputs are sparse anchors, not complete sequence alignments.

The primary baseline for this development is
`28ee06b505150dec2d74c2fb3f3d91535016e1b8`. Correctness and empirical acceptance
are recorded in the separate experiment report; this interface document does
not itself assert that benchmark gates passed.

## Conditional direction input

The primary route uses a query partition with one known strand per interval.
A query base may have both forward and reverse homologous partners when the
other sequence contains an inverted duplicate. A truth-direction diagnostic
that omits such ambiguous bases can leave no original primary anchor inside
the repeated unit. The copy layer cannot manufacture a primary projection root
from the omitted region. Such a diagnostic must report input-limited recovery
explicitly; it is not equivalent to testing the reciprocal relation with the
same anchor information.

Likewise, exact sequence strings alone do not reveal whether identical bases
were born by duplication or by independent insertion. The benchmark retains
strict source truth and separate ambiguity diagnostics. A long exact match is
not proof of shared ancestry.

## Evidence field scopes

Accepted anchor weights always use occurrences in the two complete input
sequences. Self-repeat occurrence fields include `count_scope` to distinguish
local LCP counts from complete same-sequence counts. They are evidence fields,
not interchangeable estimates of global uniqueness. Ambiguity exports separate
anchor weight from group best/runner chain scores; unscored budget-truncated
records use NA. Branch selection iterations and actual recursive gap-search
nodes are counted separately. Duplicate retained evidence is removed under the
same unique candidate budgets.

Local placement competition is evaluated over the same target interval. Support
outside that interval cannot outweigh a conflicting local correspondence.
Competing self-repeat paths must have positioning support exceeding the existing
logarithmic displacement cost; this avoids treating incidental short collisions
as localized alternatives. Equal or stronger qualifying local support leaves the
whole proposed branch ambiguous. The ambiguity output records the evaluation
interval and distinguishes local scores from whole-branch scores.
## Completion and diagnostic counters

`complete=1` means processing finished; it does not mean every local problem
was resolved. It can coexist with `budget_exhausted`. Consumers must inspect
both the summary and each task status before describing a run as complete
copy recovery.

The candidate budget counts distinct coordinate relations from self-repeat
search, projection, veto-only alternatives, and clipped relations. It is not
an accepted-anchor limit. A late budget failure rolls back all proposed
acceptances for that task, while retaining diagnostics for unfinished work.

`recursion_node_count` includes the depth-zero re-search of the whole branch
window. Inspect `max_recursion_depth` and `recursive_added_count` to distinguish
that verification from searching internal gaps. A positive node count alone
does not imply that recursion recovered new anchors.

Ambiguity `evaluation_target_begin` and `evaluation_target_end` use the axis
identified by the associated task's `side`. Join the ambiguity's task ID to
`.copy-tasks.tsv` before interpreting these coordinates as query or reference.

## Equivalent branches outside the local window

Before accepting a proposed chain, the copy layer also checks whether its
actual supported intervals can be reproduced at another physical location in
the full opposite sequence. Locations already covered by the primary chain
participate in this veto check. Translation and reverse-complement reflection
preserve the proposed chain's anchor lengths and spacing; characters inside
unaligned gaps are not required to match and receive no alignment credit.

These alternative relations are veto-only evidence, not new accepted branches
or projection roots. They use the same candidate budgets. An equal supported
placement leaves the proposed correspondence ambiguous. This check addresses
complete equivalent chains outside local search bounds; it does not claim an
exhaustive search for every distant partial or inexact alternative.

Even exact chains with matching local context do not identify a unique
ancestral history in all cases. A pre-existing repeat followed by deletion can
produce the same final strings as a post-split duplication. Source-truth
accuracy therefore remains an empirical benchmark result, not a consequence
of exact sequence equality or of a completed search status.

## Exported competition witnesses

Ambiguity rows distinguish `evidence_role` and `evidence_route`. A competitor
is veto-only evidence, not an accepted copy or a primary projection root.
`witness_order` and `witness_size` identify the exact intervals of the actual
runner path already used for the comparison. Zero proof IDs mean that the row
has no original primary/self-repeat projection certificate.

Whole-window runner witnesses are appended after existing ambiguity records
with the same group ID. Existing candidate IDs and their original fields are
preserved. Each group exports at most one actual runner path; it does not
export all combinations of repeat positions. Coordinates and weights come from
the existing evaluated path, without another search or score calculation.

## Count memoization

A private fixed-capacity cache can reuse an occurrence count when the target
input identity and actual pattern characters match. Hash hits are checked by
exact character comparison, including the logical reverse complement. Cache
misses call the original count operation; candidate order, budgets and floating
point score expressions do not depend on cache hits or replacement. The default
slot storage is 2.5 MiB on the validated Linux build, plus one reverse-complement
scratch buffer. This allocation figure is not a process RSS limit.

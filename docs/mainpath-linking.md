# Whole-genome main-path linking

This stage is integrated into ramag align, after the existing two-sided
selection and before any alignment writer. It consumes in-memory records and
FASTA sequences. It does not re-read a candidate PAF.

Ordinary `ramag align` now runs bilateral selection and bounded main-path
linking automatically, with additional vacant/chain recovery disabled:

    ramag align --reference reference.fa --query query.fa \
      --output result.paf --work-dir work --threads 16

The former selection, link-gap and link-recovery CLI options have been removed.
Old commands must omit them. The fixed workflow equals the preceding version's
explicit `linked-one-to-one` / `bounded` / `none` configuration. Compared with
the old `all` default it deliberately omits nonselected alternative placements.

Library callers can still select `All`, `ReciprocalOneToOne` or `LinkedOneToOne`
and `LinkGapMode::None` or `Bounded` through `AlignmentOptions`. The library,
batch candidate-output workflow and `align -c` defaults are unchanged.

All requested writers (SAM, PAF, MAF, delta, chain) consume the same final
records. Linking does not force a candidate PAF.

## Record and coordinate contract

- Inputs are the existing reference/query selection intersection.
- Consumed intervals, including insertions and deletions, must be disjoint on
  both axes. Invalid selected inputs fail before publishing output.
- Links require the same contig pair and strand and immediate adjacency in
  both reference and query orders. An intervening record, crossing,
  inversion boundary, or changed partner stops the path.
- CIGAR traverses reference coordinates in increasing order. Reverse query
  gaps use canonical reverse-complement coordinates; published query intervals
  stay on the original FASTA.
- Each selected record's entire CIGAR is retained in order. Newly accepted gaps
  cannot be used to propose another round of links.
- Input records remain unchanged on cancellation. A failed gap leaves the
  original split. A failed merged-CIGAR check retains its original records.
- Coalesced CIGARs are replayed to recompute HOXD70/10 scores, with gap cost
  40 + 3 * length for each maximal same-direction gap, and edit distance.
  The merged path inherits its existing primary designation.

PAF tp:A:P / tp:A:S is a writer designation, not selection evidence. For a
selected or linked run, every emitted record belongs to the fixed accepted
main path. Do not use a generic secondary-filtering PAF adapter to discard
its other records.

## Bounded gap admission

Zero gaps concatenate without a support threshold or DP.

Nonzero gaps need both original flanks to have at least min_cluster exact
bases (default 65), at least min_match inward consecutive exact bases
(default 20), and at least 98% column identity. A gap containing N is rejected.

| Gap | Bound and operation |
|---|---|
| Identical equal-length sequences | At most 10,000 bp on each side; real = |
| One empty side | At most 1,000 bp on the nonempty side; real I or D |
| Both sides nonempty, nonidentical | At most 1,000 bp each; existing fixed-end KSW2 |

KSW2 must consume both complete gaps, contain paired bases, have at least
min_cluster exact matches, and reach max(0.98, min(left_identity, right_identity)).
These are the existing whole-genome gap-fill support thresholds. There is no
clipped or I+D fallback for a rejected alignment.

DP reserves 3 * (reference_gap + 1) * (query_gap + 1) grid positions. The global
budget is the minimum of 500,000,000 and the reference/query input lengths;
library callers may tighten it. Stable support, cost, and coordinate order
determines reservations before parallel work. Work is never refunded.
Proposals run in batches of at most 64. Final assembly follows fixed reference
order.

## Diagnostics and validation

The run log separately reports statistics.candidate_alignment_count,
statistics.conflict_rejected_alignment_count, and final statistics.alignment_count.
The mainpath_link.* fields report input/output counts, zero/exact/single-sided/DP
joins, rejection categories, reserved work, added consumed bases, added paired
columns, KSW2 calls, and phase timings.

All five writers consume the same finalized records. Gap failure, insufficient
support, exhausted work budget or ambiguous geometry leaves the original split.
Fewer PAF records does not necessarily imply better homology recovery or fewer
minipoa calls in a downstream graph builder.

## Library-only recovery

Ordinary CLI alignment does not run vacant-chain or conservative chain recovery.
The library retains `LinkRecoveryMode::None`, `VacantChains`, and `ChainRecover`
for explicit callers and algorithm tests; `ChainRecover` requires bounded
linking. See [Conservative chain recovery](chain-recovery.md) and
[vacant-chain recovery](vacant-chain-recovery.md).

The experimental local-realign and chain-exchange drivers remain retired.
They are not alternate values of the ordinary CLI workflow.

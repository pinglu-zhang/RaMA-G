# Experimental coverage refinement

This experiment extends fragment re-selection. Normal builds keep the existing
baseline policy. The options below are internal CMake settings, not public CLI
arguments. They require RAMAG_INTERNAL_POST_SELECTION=fragment-reselect-v1.

| Option | Default | Operation |
| --- | --- | --- |
| RAMAG_INTERNAL_FRAGMENT_LOCAL_EXCHANGE | OFF | At most two deterministic local-exchange passes |
| RAMAG_INTERNAL_FRAGMENT_SHORT_EXACT | OFF | Add available exact remnants between min_match and min_cluster |
| RAMAG_INTERNAL_SUPPLEMENT_PRUNED_CLUSTERS | OFF | Extend best chains changed or deleted by greedy cleanup |
| RAMAG_INTERNAL_SUPPLEMENT_ALTERNATIVE_CHAIN | OFF | Extend one additional chain per original cluster |

The explicit PairwiseCoreOptions fields provide the same operations for library
callers. High-level align and serial batch apply them only to one-to-one.
Existing RAMAG_INTERNAL_GAP_FILL modes can follow fragment re-selection;
legacy residual recovery remains mutually exclusive.

## Local exchange

The input is the complete post-link candidate collection. The first pass is
unchanged fragment-reselect-v1. Candidate parent identity and priority remain
frozen. A challenger may displace one to four selected fragments. The local
pool contains parents intersecting either axis of the challenger or a displaced
fragment. A pool exceeding 256 distinct parents is skipped. Once a parent is
in that pool, its complete valid unoccupied fragments can compete; restricting
them to the original bounding rectangle would miss the intended one-to-many
exchange.

The challenger is provisionally accepted first, followed by greedy re-selection
from the local pool against all other selected fragments. Commit requires both
span totals not to decrease, strictly more paired columns, and no decrease in
canonical matches. Otherwise occupancy is restored. Two stable passes are
allowed; a pass without an accepted exchange terminates the search.

These are bounded greedy experiments, not a global optimality claim. Quality
must be evaluated independently. No acceptance decision uses MUMmer4 output or
simulation truth.

## Short exact remnants and candidate sources

Short-remnant recovery inspects only parents with at least min_cluster canonical
matches. A remnant must be completely exact, canonical, unoccupied on both axes,
at least min_match long, and shorter than min_cluster. Occupied stretches are
skipped using intervals. Accepted remnants do not seed another search.

The production chain/link path runs unchanged when supplemental candidates are
enabled. Changed or deleted pre-cleanup chains and the one alternative chain
are extended separately and do not enter cross-cluster linking. Complete
coordinate/strand/CIGAR duplicates are removed with earliest-source precedence.
The alternative chain removes only indices selected by the original best-chain
traceback, preserving multiplicity of equal-coordinate seeds.

Internal capture assigns source 0 to production, 1 to pruned-chain supplements,
and 2 to alternative-chain supplements. Its RGA00002 file has a checked record
count and trailer and is renamed from a private partial file. Source IDs are
written alongside it. This is an internal diagnostic format, not an alignment
manifest, product checkpoint, or index format change.

## Evidence and limits

Statistics distinguish exchange attempts, caps, rollbacks, paired-column gains,
short remnants and supplemental sources. Descriptor capacities are estimates,
not process peak RSS. Whole-process sampling remains necessary.

The experiment requires double-sided zero overlap, legal CIGAR consumption and
determinism. Adoption requires median simulated F1 at least the frozen baseline,
median Precision no more than 0.0005 below it, and real end-to-end time and peak
RSS at most 130% of the same-batch baseline. Four real coverage measures must
improve over fragment-reselect-v1. Exceeding MUMmer4 requires all four accepted
coverage counts to be greater; coverage alone does not establish correctness.

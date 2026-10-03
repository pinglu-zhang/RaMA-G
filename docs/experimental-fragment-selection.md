# Experimental fragment-level re-selection

This opt-in policy replaces the final two-sided whole-record DP intersection.
It does not change seeding, clustering, linking, KSW2 scoring or index handling.

Configure with `-DRAMAG_INTERNAL_POST_SELECTION=fragment-reselect-v1`.
The default remains `baseline`. It is mutually exclusive with the other post-
selection policy, legacy recovery, guarded recovery and gap fill. The CLI uses
it only with `--selection-mode one-to-one`; `all` is unchanged.
Library callers can explicitly set `PostSelectionPolicy::FragmentReselectV1`.
Recompile callers after changing the public options/statistics headers.

## Frozen selection contract

All post-link candidates are eligible, independently of their old DP flags.
Validation and the initial sequence scan are fused. Initial slices break at N
and trim terminal gaps; at least `min_cluster` canonical matches are required.
There is no new universal 98-percent identity threshold.

Parent identity is canonical matching columns divided by all alignment columns.
A slice's support is its canonical match count times its immutable parent identity.
The calculation uses double precision in that fixed order. Ties use decreasing
match count, increasing parent ID, increasing starting column, then decreasing
ending column. Clipping cannot increase this priority.

Reference and original-forward query occupancy start empty. A free slice is
accepted; conflicting slices are cut along CIGAR columns on both axes and their
qualified remnants re-enter the queue. Remnants are never joined across a
conflict. The result is deterministic greedy packing, not a global optimum.

Whole accepted parents reuse their packed CIGAR; clipped pieces acquire CIGAR
storage only when accepted. All accepted paired columns come from original
candidates. No new KSW2 or seed search is performed. Original selected records
and primary flags need not survive; the final sorted set receives normal
deterministic primary assignment.

## Interfaces and evidence

The public core still exposes original candidates with selection flags plus
accepted clipped records. Consumers of that low-level result must honor those
flags. The CLI bridge serializes only the selected set. The sequence-free
record-selection API retains its existing behavior and does not run this policy.

Logs report the actual strategy, validated parents, proposed slices, re-clipping,
whole/clipped acceptance, support rejection, N breaks, selected spans and paired
columns, time and memory observations. Capacity accounting covers descriptor,
identity and output-vector storage; it is not a complete allocator or RSS estimate.
Process RSS must be measured separately.

The private replay harness verifies captured online baseline flags before using
the same production selector. Its clipped PAF is not an added-pair audit:
re-selection may both gain and lose records. Pair-coordinate set differences and
independent truth membership are required to judge changes.

No quality, coverage or performance improvement is asserted by this document.
Acceptance requires frozen simulation Precision/F1 nonregression, four real
coverage gains, and time/memory limits. No default is switched automatically.

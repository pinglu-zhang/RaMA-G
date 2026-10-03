# Experimental one-sided guarded post-selection

This opt-in experiment addresses whole-record rejection after the existing
reference/query dynamic-programming selections. It does not change seeding,
clustering, linking, KSW2 scoring, selection weights, or the 5,000-candidate window.

Configure with `-DRAMAG_INTERNAL_POST_SELECTION=one-sided-guarded-v1`.
The default is `baseline`. No public command-line option is added.
The experiment applies only to high-level `one-to-one`; `all` is unchanged.
It is mutually exclusive with legacy residual recovery, guarded recovery, and
gap-fill. Direct C++ callers use `PairwiseCoreOptions::post_selection` or the
existing `ConfiguredPairwiseCoreOptions` conversion. Recompile clients after
updating the public options/statistics structures; binary ABI stability is not claimed.

Only candidates selected on exactly one axis may contribute. Original
dual-selected records are immutable. Both reliable flanks must be original
dual-selected records on the same contig pair and strand, with co-linear,
unoccupied 1–10,000-base gaps on both axes. Each flank requires at least
`min_cluster` canonical matches, 98% real identity, and `min_match` inward exact
bases. Windows containing N are rejected; windows are frozen and never recurse.

Packed CIGAR columns are intersected with both window coordinates and both
occupied-interval indexes. Either-axis conflicts break a fragment. Terminal
gaps are removed without stitching disjoint column ranges. A fragment requires
`min_cluster` matches and identity at least the greater of 98% and the lower
flank identity. Identity is canonical matches divided by all columns, not the
packed-M support proxy used by the unchanged main selector.

Competition uses the existing deterministic residual-fragment ordering.
Accepted fragments reserve both axes; stale competitors are clipped and ranked
again. All new fragments receive validated coordinates, score, NM and CIGAR.
Original primary status is preserved. No new KSW2 or suffix-array operation is
performed by this step.

Logs report `post_selection_*` counters and existing `recovery_*` totals.
Rejection counters count scan events and may include rechecks, not distinct
parents. Parent accounting separately partitions XOR-selected parents into
those intersecting a reliable window and those without one. Paired columns
exclude I/D; span additions include consumed I/D bases. Auxiliary capacity
reports tracked window/index vectors and queue capacities, not total RSS;
RSS is separately sampled and in-process samples are not whole-command peaks.

Production success still means exit code zero and validated published outputs.
There are no runtime manifest, hash-list, or completion-marker additions.
The private replay tool and tests are not part of the installed library.

Correctness or coverage gains alone do not establish biological accuracy.
Adoption requires frozen simulated Precision/F1 non-regression, preservation
of original mappings, deterministic output, and separately measured real-data
coverage/time/memory gates. Until those gates pass this remains an experiment.

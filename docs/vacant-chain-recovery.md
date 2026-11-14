# Supported vacant-chain recovery

This stage is available to explicit library callers through
`AlignmentOptions::link_recovery_mode = LinkRecoveryMode::VacantChains`, with
`selection = AlignmentSelection::LinkedOneToOne`. Ordinary CLI `align` uses
bounded linking without additional recovery; there is no recovery CLI switch.
Centromere and batch defaults are unchanged.

## Protected result

The existing bilateral selector and bounded main-path linker run first. Their
records, including all filled gap pairs, are frozen. Original rejected packed
CIGARs remain alive until recovery completes; only selected records and bounded
local proposals are materialized.

A recovery window must have two existing flanks that are mutual direct neighbors
on both physical axes, with the same contig pair and direction. Both axes must
be vacant. Each flank needs at least min_cluster exact bases (default 65), an
inward exact run of at least gap_fill_min_match (default 20), and identity 0.98.
Reverse query coordinates are normalized before ordering and restored at output.

## Evidence and placement

Donor evidence comes only from existing candidate CIGARs. It is clipped by actual
CIGAR consumption, not bounding-box interpolation. Query-overlapping alternatives
at other reference loci, on other strands, or on occupied partners remain in the
competition. Equal physical alignments are deduplicated. A different placement
covering the same complete query interval with an equal or higher HOXD70 score
makes that fragment ambiguous. This is a local conservative placement check,
not an optimal genome-wide chain assignment.

A donor needs inward exact runs and identity at least the smaller of the two
flank identities (never below 0.98). Candidate chains are ordered and disjoint.
The deterministic DAG prioritizes donor scores minus gap imbalance penalties;
that priority is not reported as an alignment score. Equal best proposals are
rejected. Chosen seams must independently pass the existing fixed-endpoint
linker and its true whole-record score must strictly improve on the two flanks.

## Atomic attachment and resources

Only a complete left-flank to donor-chain to right-flank join is published.
Unresolved seams discard the complete proposal. No isolated fragment is emitted.
New records do not seed another round. Adjacent windows can share an old flank;
it is emitted once. The final record count cannot exceed the frozen baseline,
and the original pair set and both consumed coverage sets must be retained.

The existing global DP work budget is shared with the first linker; attempted
recovery work is charged even if its chain fails. Metadata is bounded to 512 MiB
plus fixed window structures, at most 4096 deduplicated descriptors per window.
Competition and DAG pair checks each have a 2,000,000 bound per window. An
incomplete evidence window is rejected, never called unique after truncation.
Local linking reuses the validated immutable sequence lookups instead of rescanning
the complete FASTAs for each proposal. Cancellation checks precede publication;
output replacement is transactional.

Logs include window coordinates, original flank and donor IDs, accepted/rejected
reasons, counts, budget use and time. This mode does not serialize and re-read a
candidate PAF, nor use simulation truth or downstream graph output for selection.

## Relation to conservative chain recovery

`VacantChains` remains a supported library mode with its original admission and
scheduling rules. chain-recover starts from that accepted result and performs
one additional bounded conservative recovery transaction. The vacant-chain transaction completes before the additional conservative
transaction begins.
See [Conservative chain recovery](chain-recovery.md).

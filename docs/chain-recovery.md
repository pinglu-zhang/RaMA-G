# Conservative chain recovery

Conservative recovery is retained as a library capability, not a CLI option.
Ordinary `ramag align` always uses bilateral main-path selection and bounded
linking with additional recovery disabled. The centromere and batch workflows
keep their existing behavior.

Explicit library callers may set `AlignmentOptions::selection` to
`AlignmentSelection::LinkedOneToOne`, `link_gap_mode` to `LinkGapMode::Bounded`,
and `link_recovery_mode` to `LinkRecoveryMode::ChainRecover`. The command-line
selection/link-gap/link-recovery switches have been removed.

## Protected baseline and candidates

The existing bilateral selection, bounded linker and accepted vacant-chain
recovery run first. Their result is the protected baseline B0. Conservative
chain recovery keeps every original coordinate pair and both consumed spans;
it may add pairs or merge records, and it cannot increase the record count.

Candidates are clipped from original parent CIGARs. Competing placements,
including occupied reference loci and opposite strands, remain in the evidence
universe. Incomplete competition is unresolved. All decisions use the input
FASTA and algorithm evidence; simulation truth is only used for later scoring.
No newly accepted path creates another search round.

A full core is admitted only after its actual CIGAR, placement and protected
consumption checks pass. Equal-score affine-gap translations mark an ambiguous
core envelope. Stable outer portions may be attached separately while the
uncertain interval remains an actual unresolved gap. No synthetic I+D patch is
used to hide a failed gap. Reverse-query work uses oriented coordinates and
returns original-FASTA half-open coordinates.

## Shared tools and resources

The private repair tools contain metadata lookup and bounded gap alignment
shared by the retained recovery implementation. They no longer require the
retired broad local-realignment driver or the exchange driver. Scoring remains
HOXD70/10 with maximal same-direction gap cost 40 + 3 * length.

One invocation-local ledger controls all optional recovery work. The limits
remain 512 MiB of auxiliary capacity reservation, a 10-second task watchdog,
and a stage time cap of min(300 seconds, 25% of seed enumeration plus B0 core
work). The existing cumulative candidate, local-base, seed-hit, edge, replay
and remaining DP limits also apply. The auxiliary reservation is not the full
program RSS. Exhaustion retains B0 and records the reason; caller cancellation
and malformed inputs propagate without publishing a partial transaction.

Candidate precollection consumes the cumulative local-base budget. A window
with incomplete placement competition remains unresolved if this budget is
exhausted before recovery.

## Library compatibility

The retired local-realign and chain-exchange values remain invalid in the library.
The public experimental functions RealignUncoveredPaths and
ExchangeChainAlignmentPaths are removed; the retained public entry is
RestoreChainAlignmentPaths. Retained LinkRecoveryMode numeric values stay
0 (None), 1 (VacantChains) and 3 (ChainRecover); retired values 2 and 4
remain invalid.

The legacy LocalRealignmentStatistics type and local_realign.* log keys
continue to describe the retained chain stage for compatibility with existing
analysis outputs. Keeping those names does not retain either retired driver.

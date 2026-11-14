# Whole-genome main-path linking

Ordinary `ramag align` runs bilateral one-to-one selection and bounded linking
before any output writer. The same final records feed SAM, PAF, MAF, delta and
chain. It does not write or reread a candidate PAF as an intermediate step.

Only mutual direct neighbors on both physical axes are eligible. They must use
the same contig pair and direction. Overlaps, crossings, intervening mappings,
changed partners or inversion boundaries leave separate records. Reverse query
work uses oriented coordinates and restores original-FASTA coordinates at output.

Zero gaps concatenate existing CIGARs. Nonzero gaps require each flank to have
at least `min_cluster` exact bases (default 65), `min_match` inward consecutive
exact bases (default 20), and at least 98% column identity. Gaps containing N are
not filled.

| Gap | Bound |
|---|---|
| Both equal-length sequences identical | At most 10,000 bp per side; actual `=` |
| One side empty | At most 1,000 bp on the nonempty side; actual I or D |
| Both nonempty and different | At most 1,000 bp per side; fixed-end KSW2 |

KSW2 must consume both gaps completely, contain paired bases, have at least
`min_cluster` exact matches and satisfy the existing flank-derived identity
threshold. Rejection leaves the original split; there is no clipped or I+D
fallback. DP reservations keep the 500M upper bound and input-size cap; work is
charged deterministically before parallel execution.

Whole-record CIGAR replay recomputes score, NM and consumption, including seam
gaps. Original selected pairs are preserved. Normal logs retain candidate,
selection, output and linking statistics. Fewer records do not necessarily imply
higher accuracy or fewer downstream minipoa calls.

Explicit library callers may choose `LinkGapMode::None` for zero-gap joining or
`Bounded`; ordinary CLI uses bounded linking. Additional vacant/chain recovery
and older coverage experiments are retired, as described in [migration](migration.md).

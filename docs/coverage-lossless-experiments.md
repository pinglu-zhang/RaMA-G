# Coverage and lossless optimization experiments

These are internal, independently selectable experiments, not new production
defaults. The ordinary build keeps the historical selection and connection
policies. The public CLI, KSW2 scoring, suffix-array format and dependency pin
are unchanged. `batch` remains serial and reuses its reference index.

## Coverage policies

| Build option | Default | Effect |
| --- | --- | --- |
| `RAMAG_INTERNAL_PRESERVE_LINK_CANDIDATES` | `OFF` | Keep valid overlapping candidates not consumed by a successful connection. Do not re-align skipped candidates. |
| `RAMAG_INTERNAL_LINK_PRECHECK` | `OFF` | Reject incompatible IDs, strand or gaps before ranking within the existing candidate window. |
| `RAMAG_INTERNAL_GAP_FILL` | `off` | Independently select `exact-gap` or `ksw2-gap` with reliable original flanks; residual recovery is not required. |
| `RAMAG_INTERNAL_GUARDED_RECOVERY` | `OFF` | Recover only unused CIGAR fragments within reliable original-flank intervals and satisfying the same identity threshold. |
| `RAMAG_INTERNAL_COVERAGE_RECOVERY` | `OFF` | Historical unrestricted residual policy, retained solely as a separate comparison. |

The guarded policy, when explicitly enabled, takes precedence over the
historical recovery policy. Recovery and gap filling apply only to `one-to-one`
in the CLI. Candidate preservation and legality prechecking apply before
selection. No policy increases the 10 kb connection limit or changes KSW2
scoring. Added candidates can change the independent selections and their
intersection; more candidates are not by themselves proof of better coverage.

`--print-effective-config` and the beginning of `run.log` report the actual
compiled policies. Library callers may call `ConfiguredPairwiseCoreOptions()`
to obtain the same conversion, or explicitly configure `PairwiseCoreOptions`.
Direct core defaults remain off. Library calls do not install logging or signal
handlers. Neither operation produces a manifest or completion marker.

Guarded recovery uses only original dual-selected flanks, both at least 98%
canonical identity, with `min_cluster` match support and inward `min_match`
consecutive exact bases. Both free gaps are 1–10,000 bases and contain no N.
A recovered fragment must lie inside both gaps, have `min_cluster` canonical
matches, and have identity at least the larger of 98% and the lower flank
identity. Occupancy is checked again before acceptance. Added fragments never
become new flanks.

## Diagnostics and adoption

Private coverage captures retain stable group-local candidate IDs and link
dispositions. The fixed-candidate replay tool also supports `canonical` weights;
this is a diagnostic only, not a production scoring change. The old proxy counts
paired KSW2 M columns, including mismatches, whereas the diagnostic counts true
canonical matches. Other ordering, window and floating-point rules stay fixed.

Coverage candidates require same-batch, unrounded median Precision and F1 not
below the baseline, better reference and query span and paired-column coverage,
and no more than 10% additional wall time or peak memory. Real assemblies have
no validated base-level truth and must not be assigned an F1 from another
aligner's output. Small correctness tests are not full-genome acceptance.

After a coverage baseline is accepted and frozen, performance changes must
preserve complete seed tuples, finalized records, CIGAR, primary flags and
alignment payload. Raw callback counters may change only for a documented,
proven scan short-circuit; discarded occurrence counts must not be invented.

Experimental outputs and internal tests remain outside the public source tree.
No strategy is automatically enabled or released based on these descriptions.

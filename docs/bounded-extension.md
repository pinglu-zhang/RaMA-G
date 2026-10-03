# Bounded short-chain and endpoint experiments

These are default-off scientific candidates, not adopted production policies.
Both require `RAMAG_INTERNAL_POST_SELECTION=fragment-reselect-v1` and are applied
by the CLI only to `one-to-one`. `all` keeps its existing semantics.

* `RAMAG_INTERNAL_SUPPLEMENT_SHORT_CHAIN=ON` enables S. Main and alternative
  chains retain their existing minimum support. Previously rejected best chains
  with 50–64 exact supporting bases may enter an independent extension path.
  Their existing seed-to-seed alignment remains in the middle of the CIGAR.
  Each end can inspect at most 500 bases per axis, without crossing `N` or a
  contig boundary. The resulting candidate requires 65 canonical matches and
  95% actual identity before fragment selection. This is not an extra seed search.
* `RAMAG_INTERNAL_EXTEND_SELECTED_ENDPOINTS=ON` enables E. Selected flanks with
  65 matches, 98% identity and 20 exact terminal bases may extend into free
  sequence. Flanks and windows are frozen once. Exact new fragments require
  20 bases; imperfect fragments require 65 matches and 95% identity. Occupancy
  is rechecked during stable acceptance; new fragments cannot become flanks.
  E and existing gap filling are deliberately separate experiments.

Each strategy reserves up to `min(5e9,2*min(reference_bp,query_bp))` grid positions
per input query file. Each window reserves `3*(r+1)*(q+1)`. Tasks are ordered by
support quality, work and stable coordinates. Reserved work is never refunded
based on execution speed. This models work, not seconds or total CPU instructions.
New DP uses the existing HOXD70/10, gap 40/3 and connector endpoint Z-drop 400.
Results are batched in groups of at most 64 with at most four computing workers.
Allocation pressure retries the same frozen batch at lower concurrency; failure
at one worker propagates rather than silently discarding candidates.

Public `PairwiseCoreOptions` exposes the same default-off strategies. The 50/65
thresholds are independent of `min_cluster`; window, work cap and worker count
can be reduced. Invalid ranges fail rather than overflow. `PairwiseAlignment`
adds endpoint/source annotations without changing the public Seed layout.
Source IDs are 0 main, 1 pre-cleanup, 2 alternative, 3 short-chain and 4 endpoint.
Private replay accepts old 0–2 sidecars and new 0–4 sidecars with full validation.

Statistics separately report candidates, budget skips, admitted work, calls,
quality/conflict rejections, accepted matches/paired columns, elapsed time,
estimated descriptor capacity and sampled RSS. These values are not summed to
estimate the full process peak. Accepted S candidates can still lose columns
at fragment selection; accepted S support is not final output coverage.

Acceptance requires median simulated F1 no lower than the same-batch B+D,
strictly increasing four-axis real coverage, time at most 130%, and both GNU and
process-tree peak RSS at most 103%. Shared-load interference is reported. No
performance or quality improvement follows merely from enabling these options.
No manifest, checksum inventory, completion marker or product CLI is introduced.

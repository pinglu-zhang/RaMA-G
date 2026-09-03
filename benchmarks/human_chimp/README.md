# Fixed simulated Human–Chimp comparison harness

This directory drives the one-run, 16-thread, reference-unique MAM comparison
against the independently packaged Alignathon benchmark.  It does not alter the
benchmark inputs, truth, or evaluator.  New evidence is written only below the
package's `results/preliminary-mam16-single-<UTC timestamp>` directory.

The fixed protocol is:

- reference `simHuman.fa`, query `simChimp.fa`;
- RaMA-G `--seed-mode mumreference`, MUMmer4's default MUMREFERENCE mode;
- `min-match=20`, `min-cluster=65`, `max-gap=90`, `diag-diff=5`,
  `diag-factor=0.12`, `break-length=200`;
- one MUMmer4 run followed by one RaMA-G run, each with a requested budget of
  16 threads;
- both commands receive the recorded environment `OMP_NUM_THREADS=16`,
  `OMP_DYNAMIC=FALSE`, `OMP_PROC_BIND=SPREAD`, `OMP_PLACES=CORES`,
  `OMP_THREAD_LIMIT=16`, and `OMP_MAX_ACTIVE_LEVELS=1`; `OMP_WAIT_POLICY` is
  left unset.  The binding policy is part of RaMA-G's recorded OpenMP runtime
  configuration; MUMmer4 does not use OpenMP but runs in the same environment;
- immediately before **each** full aligner attempt, read `simHuman.fa` and then
  `simChimp.fa` sequentially through EOF while recomputing their accepted
  SHA-256 values; a mismatch refuses launch, and the per-tool evidence remains
  inside that run's `attempt-NNN` directory;
- delta-only alignment timing under GNU `/usr/bin/time -v`, paired
  `CLOCK_MONOTONIC` and UTC/`CLOCK_REALTIME` endpoints, plus one-second Linux
  process-tree RSS/thread sampling with paired clock samples;
- `all-homology`, one million samples, `near=0`, and seeds
  `20260830,20260831,20260832` for quality scoring.

The harness never invokes `delta-filter -1`.  `delta-filter` and `show-coords`
are used only as raw-delta parser/readability gates.

Preflight rejects a RaMA-G binary unless `--version` reports OpenMP enabled,
the 16-thread runtime maximum, and OpenMP-enabled pinned divsufsort support;
Linux `ldd` must also expose `libgomp` or `libomp`.  A RaMA-G run is accepted
only when its manifest independently records the same runtime, requested thread
budget, disabled dynamic teams, build capabilities, complete stable 4 MiB
`tiled-mam-boundary-mem-v1` tasks, global maximality/reference-uniqueness
revalidation, the seed OpenMP team size, and active workers in chaining and
extension.  Smoke accepts the task-derived seed team: chrD must contain exactly
six tile tasks and four boundary tasks, for 10 completed tasks and a 10-thread
team.  The full four-contig query must contain exactly 92 tile tasks and 84
boundary tasks, for 176 completed tasks and a 16-thread team.  At the fixed
64 MiB reference threshold, the actual smoke backend must be `divsufsort32`
and the full Human reference backend must be `caps32`.  Both smoke and full
runs pass the runtime external-FAI acceptance gate.  Final summarization repeats
the full run's source/copy/hash/size/verification checks against the live source
FAI and retained active copies, and requires the reference/query input sections
to record 2 requested and 2 actual OpenMP workers.  These counts and routes are
independently rechecked before publication.
The one budget is reused by sequential stages and does not create nested
16-thread teams.

## Current implementation evidence

The production `mumreference` path no longer uses one whole-contig task per
orientation.  Each oriented query is split into stable non-overlapping 4 MiB
MAM tiles.  Every internal tile boundary receives a separate
occurrence-complete MEM recovery task; candidates are mapped back to the whole
reference/query, extended, checked for global maximality and reference
uniqueness, then deterministically deduplicated.  This route is named
`tiled-mam-boundary-mem-v1` in the manifest.  It does not claim sufkit's future
query-range or caller-owned-workspace APIs.

Correctness evidence currently includes 28,800 exhaustive single-contig A/C/N
cases, targeted multi-contig/repeat/reverse cases, and identical 1/2/8/16-thread
seed sets and semantic statistics.  A chrD development smoke produced the same
delta SHA-256
`b94dd59c10225b5f0b3c0142725b0b213ae429467df47aeddc9507f352c1abaf` at one
thread and a 16-thread request.  GNU time was 16.61 s and 8.88 s respectively;
the latter had six tile tasks plus four boundary tasks, so the manifest
correctly reported 10 tasks and a 10-thread seed team.  This field records the
OpenMP team size; process sampling and CPU time independently show whether the
team performed parallel work.  These are single development observations, not
the full Human-Chimp comparison.

Focused chrD A/B retained sufkit's automatic algorithms.  MAM `auto` and
forced suffix-link had identical counts/checksums, and this Fast index already
resolves MAM auto to suffix-link.  Boundary MEM `auto` and forced LCP also had
identical counts/checksums; for 4,000 calls their five-run medians were about
0.002560 s and 0.002580 s, while forced suffix-link was about 0.006813 s.
These microbenchmarks justify no hard-coded override: production keeps auto,
while the manifest records the pinned index acceleration and tiled seed route;
this A/B evidence documents how that pinned sufkit revision resolves auto.

The full four-contig run has stricter acceptance.  MUMmer4 runs first and must
publish a non-empty parser-validated raw delta plus accepted GNU/monotonic/
realtime evidence.  Its same-root monotonic elapsed is then frozen in
`speed-budget.json` and becomes the RaMA-G hard deadline.  RaMA-G must create
the complete 16-thread seed team and independently expose parallel threads to
process-tree sampling before that deadline.  If it is still running when
`CLOCK_MONOTONIC elapsed > cutoff`, the runner terminates only its own process
group (SIGTERM, ten-second grace, then SIGKILL if necessary), retains every
partial artifact, writes `RUN_TIMED_OUT`, and stops before MAF/F1.  Even a
passing run remains a preliminary single observation; it is not a stable
speedup or release-level claim.

Clock handling has two layers.  At attempt level, GNU elapsed, monotonic
elapsed, and UTC/realtime elapsed are always recorded.  Non-finite or
non-positive elapsed values, contradictory endpoints, malformed samples, or a
backwards/stopped realtime sample hard-quarantine the attempt.  Pairwise drift
above 1% and positive realtime steps are retained as `CLOCK_WARNING`, because
WSL can apply the same periodic wall-clock correction to both tools.  The same
resource loop writes `timing/resources.clock.tsv`; two or more positive jumps
are additionally labelled `periodic_realtime_jump_detected`.

If RaMA-G completes before the hard deadline, the driver immediately requires
both RaMA-G GNU elapsed <= MUMmer4 GNU elapsed and RaMA-G monotonic elapsed <=
MUMmer4 monotonic elapsed; equality passes.  Only then is
`ALIGNMENT_SPEED_ACCEPTED` created and delta-to-MAF/strict validation/F1
allowed.  A completed clock failure creates `ALIGNMENT_SPEED_REJECTED` and
also stops before MAF/F1.  The final report recomputes `GNU/monotonic` and
`realtime/monotonic`; cross-tool distortion differences above 1% and positive
realtime jumps are retained as warnings and diagnostics, but do not block
`COMPARISON_COMPLETE`.  Any hard-quarantined or timed-out attempt keeps raw GNU
time, stdout/stderr, resource samples, delta/work files, and the clock/timeout
report, but never receives `RUN_ACCEPTED`.

## Safe staged use

First run only the non-alignment preflight with the fresh, already tested
Release binary.  The binary is copied into the immutable experiment bundle:

```bash
python3 benchmarks/human_chimp/run_human_chimp_benchmark.py \
  --stage preflight \
  --ramag-binary /absolute/path/to/fresh-release/ramag
```

The command prints the exclusive result root.  Continue that exact root with
the copied, hash-accepted driver and explicit resume stages:

```bash
RESULT=/absolute/path/to/benchmark-results/preliminary-mam16-single-YYYYMMDDTHHMMSSZ

python3 "$RESULT/scripts/run_human_chimp_benchmark.py" \
  --resume \
  --result-root "$RESULT" \
  --stage smoke

python3 "$RESULT/scripts/run_human_chimp_benchmark.py" \
  --resume \
  --result-root "$RESULT" \
  --stage full

python3 "$RESULT/scripts/run_human_chimp_benchmark.py" \
  --resume \
  --result-root "$RESULT" \
  --stage evaluate

python3 "$RESULT/scripts/run_human_chimp_benchmark.py" \
  --resume \
  --result-root "$RESULT" \
  --stage summarize
```

`--stage all` performs every remaining stage and is intentionally not the
default.  The existing incomplete
`preliminary-ramag-vs-mummer4-20260830T045046Z` directory is explicitly
forbidden as a target.

Every failed run or evaluation is retained as `attempt-NNN`; resuming creates a
new attempt and never overwrites it.  `COMPARISON_COMPLETE` is published last,
only after `MUMMER_BASELINE_ACCEPTED`, `ALIGNMENT_SPEED_ACCEPTED`, both accepted
deltas, canonical MAF files, three XML files per tool, structured metrics, and
the comparison JSON/TSV/Markdown exist and pass their gates.  The MAF conversion,
source canonicalization, strict validation, and F1 stages each retain separate
monotonic timing evidence and do not enter the alignment speed gate.  Immediately
before the final marker, the summarizer re-reads the selected
RaMA-G full and smoke manifests and independently revalidates full=`caps32`,
smoke=`divsufsort32`, the 64 MiB threshold, external-FAI source/copy hashes and
sizes, structure verification, and the 2/2 OpenMP input route.  The backend
wording in `comparison.md` is generated from this accepted manifest evidence,
not from a hard-coded narrative.  F1 is report-only: missing, non-finite,
out-of-range, pair-count-inconsistent, or non-reproducible scores block
completion, but no minimum F1 or RaMA-G-vs-MUMmer4 F1 winner is required.  If
either alignment clock layer fails, `ALIGNMENT_SPEED_REJECTED` is retained,
evaluation and comparison files are not created, and a replacement observation
must use a new exclusive experiment root.

Immediately before each new full aligner attempt, the driver repeats the
resource gate: no other `ramag`/`nucmer` process, at least 20 GiB available
memory, and at least 50 GiB free space on the experiment filesystem.  The
timestamped evidence is retained under `runs/launch-gates/<tool>/attempt-*`.
An accepted run is skipped on resume; a failed alignment is preceded by a new
resource check before its replacement attempt.

After that resource gate creates a new full-run attempt, the intended command
is frozen as `planned-command.json`.  The driver then performs the per-tool
reference-then-query full read and writes `warm-cache.json` plus
`WARM_CACHE_COMPLETE`; `run_with_metrics` is the next operation on the code
path.  `warm-cache-launch.json` records the monotonic interval between the end
of that read and the runner start.  Hash failure creates `WARM_CACHE_FAILED`
and refuses to start the aligner.  Clock assessment writes
`clock-consistency.json` and `CLOCK_ACCEPTED`; warning-quality attempts also
receive `CLOCK_WARNING`, while hard failures receive `CLOCK_QUARANTINED` and
remain unselected.  Resume skips an already selected run, but every replacement
attempt repeats both the resource gate and its own input read, so failed or
quarantined evidence is never repurposed.

## MAF source-name gate

`canonicalize_maf_sources.py` recognizes only the eight canonical contigs and
their eight exact doubled-prefix aliases, for example:

```text
simHuman.simHuman.chrA -> simHuman.chrA
simChimp.simChimp.chrD -> simChimp.chrD
```

Unknown, cross-species, tripled-prefix, and normalization-collision names are
rejected.  The JSON audit records input/output SHA-256, every replacement,
canonical sources seen, and independent hashes proving that non-`s` lines and
all `s`-line fields after the source name were unchanged.

## Local harness tests

These tests use temporary fixture files and a millisecond-scale process sample;
they do not start either aligner or touch the benchmark package:

```bash
python3 -m unittest discover -s benchmarks/human_chimp/tests -v
```

The final report deliberately calls the measurements a preliminary single
observation.  Runtime has no repetition median, and the sampled F1 is not an
exhaustive exact score.

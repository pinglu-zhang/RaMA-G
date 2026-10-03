# Coverage candidate: B+D+C+S+E, min-cluster 40 / max-gap 500

This named profile preserves the configuration used in the parameter screen.
It is opt-in: ordinary CMake builds and CLI parameter defaults are unchanged.
The algorithm implementation was already present in main before this profile
was added. No new scoring, seed search, dependency or index format is introduced.

## Build

From the repository root, using CMake 3.22 or newer:

```bash
cmake --preset coverage-40-500
cmake --build --preset coverage-40-500 -j 8
```

The preset enables fragment-level re-selection, short exact residuals (B),
an alternative chain (D), pre-cleanup cluster supplements (C), bounded short
chain extension (S), and bounded selected-end extension (E). Other coverage
experiments remain off. Dependencies retain the repository's pinned identity.
An existing clean pinned Sufkit checkout can be supplied at configuration time
using the existing `-DRAMAG_SUFKIT_SOURCE_DIR=PATH` option.

## Run

The parameter values are CLI options, not CMake cache variables. Both must be
passed explicitly, along with the seed and selection modes:

```bash
./build/coverage-40-500/ramag align \
  --reference reference.fa --query query.fa \
  --reference-index reference.sufidx \
  --output result.paf --work-dir work \
  --threads 16 --seed-mode mumreference --selection-mode one-to-one \
  --min-match 20 --min-cluster 40 --max-gap 500 \
  --diag-diff 5 --diag-factor 0.12 --progress off
```

Omit `--reference-index` only if an in-memory build is intended. Explicit index
loading continues to validate the index and reference without rebuilding.
Output and work paths should be unique per run. For a non-executing configuration
check, append `--print-effective-config` to the command above. It must show
`min_cluster=40`, `max_gap=500`, `fragment-reselect-v1`, and all five B/D/C/S/E
flags enabled. The same build is usable with serial `batch`; use the same explicit
alignment parameters with batch's own query and output-directory arguments.

## Validation status

The complete Human–Chimp parameter screen measured this configuration using
16 threads, a fixed 16-CPU mask, full inputs, a byte-coded reference index, and PAF-only
output. Coverage used the full reference/query denominators (3,298,430,636 and
3,177,756,316 bases). Paired-column coverage includes matching and mismatching
two-sided columns but excludes single-sided insertion/deletion columns.

The completed September 30, 2026 parameter screen measured:

| Configuration | Reference span | Query span | Reference paired-column | Query paired-column |
|---|---:|---:|---:|---:|
| 40/500 | 82.636777% | 85.762682% | 82.354406% | 85.481789% |
| 50/250 | 82.362983% | 85.485125% | 82.122354% | 85.240926% |
| Historical MUMmer4, filtered -1 | 82.547809% | 85.504488% | 82.288551% | 85.250706% |

On the complete simulated dataset, all-homology, near=0, one million samples,
and seeds 20260830/20260831/20260832 gave median F1
0.9912596925953074 for 40/500 versus 0.99127754492142 for 50/250.
Thus 40/500 did not pass the strict F1-nondecrease gate. Its inclusion records
the user's explicit choice of this higher-coverage configuration, not a claim
that the quality gate passed.

The real 40/500 command took 1260.05 seconds with GNU time peak RSS 54.952 GiB.
Shared server load left the formal performance conclusion unconfirmed.
These are single-run observations, not evidence of universal quality or speed.

To reproduce the previous comparison, use the preserved `coverage-50-250`
configure/build preset and explicitly pass `--min-cluster 50 --max-gap 250`.
Both presets enable the same algorithm flags; the parameter pair must still
be supplied on the command line.

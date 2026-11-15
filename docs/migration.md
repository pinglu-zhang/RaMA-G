# Migration from experimental interfaces

## Fixed ordinary alignment workflow

Ordinary `ramag align` selects the bilateral main path and links it with bounded
gap filling. Extra recovery is disabled. Remove these former command-line options:

- `--selection-mode`
- `--link-gap-mode`
- `--link-recovery-mode`

They report an unknown-option error. The ordinary workflow matches the preceding
explicit linked-one-to-one/bounded/no-recovery configuration, not the older all-
candidate default. `batch` retains candidate output, and `align -c` retains its
complete centromere workflow.

## Removed experimental builds and APIs

The two coverage presets and their residual recovery, fragment re-selection,
one-sided re-selection, short exact/chain supplements, endpoint extensions,
alternative-chain and geometry experiments are retired. Old experimental cache
variables are rejected; configure a new directory with `release`, `genome-only`
or `testing`.

`VacantChains` and `ChainRecover`, their `LinkRecoveryMode` enum and recovery
fields, `RestoreVacantAlignmentChains`, `RestoreChainAlignmentPaths`, recovery-
specific statistics and the repair-budget header have been removed. Library
consumers using them must remove those calls and recompile. The retained main-
path linker, its normal work cap and centromere resource limits are unchanged.

The product keeps basic candidate/unlinked/linked library selection and
centromere direction, anchor, copy, stitch and format APIs. Private workflow and
review headers have moved into source modules; do not include them from an
external consumer. Supported installed public include paths remain unchanged.

## Research scripts and evidence

The top-level benchmark harness is no longer part of the product tree. Historical
source remains in Git history; research results and internal benchmark records
are not rewritten or presented as new product measurements. Real-data benchmarks
are not required for installation or public CI.

## Stage-oriented source layout

The ten supported public header names, library entry points and install targets
are unchanged. Private headers have moved from `src` to
`include/ramag/internal`; unsupported source-tree consumers must update their
includes. Private headers are excluded from installation. Internal tests now
link the compiled genome implementation instead of including its `.cpp` file.

Related runtime, I/O, index and refinement declarations are consolidated.
Direction diagnostic writers are separated from the algorithms, and the genome
extension and main-path linker have independent translation units. The source
reorganization does not introduce CLI options or change alignment defaults,
scoring, tie order, resource limits or output coordinate/CIGAR contracts.

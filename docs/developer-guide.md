# Development and troubleshooting

## Source map

Project headers live under `include`; implementation files live under `src`.
The ten existing `include/ramag/*.hpp` paths remain the supported public API.
`include/ramag/internal/` contains private declarations and required templates;
it is shipped with source packages but is never installed or used by public headers.

```text
include/ramag/              supported public headers
include/ramag/internal/     private app, I/O, index, genome and centromere headers
src/app/                   CLI, orchestration, progress, logging and cancellation
src/io/                    FASTA input and common format writers
src/index/                 reference-index persistence and Sufkit adaptation
src/genome/                alignment orchestration, core chaining/selection,
                           KSW2 extension and bounded main-path linking
src/centromere/            backend, direction, anchors, copies, stitching,
                           refinement and consolidated result writers
third_party/               pinned dependencies, licenses and source attribution
```

The genome entry implementation delegates to the core, extension and linking
translation units. Centromere direction processing separates evidence collection,
planning, orientation assignment and compiled review rules. Recursive anchors,
copy discovery, internal-gap stitching and later path refinement remain distinct
stages. Diagnostic tables and alignment-format rendering are consolidated in the
centromere output implementation rather than mixed into direction search.

Each stage declares cross-file dependencies through private headers. There are
no project headers or implementation fragments in `src`, and no implementation
file is textually included by another source file or test. Ordinary review entry
points are compiled once; callback-dependent search visitors retain statically
typed templates without runtime type erasure or extra data copies.

The small `centromere_overlap_compatibility.hpp` header remains separate because
the pinned RaMA sources include that exact name. Its include directory is private
to the dependency build. This preserves the third-party source and license files.
Genome HOXD70 and centromere three-piece kernels remain separate scoring models.

Public headers must be self-contained and must not include private headers.
I/O and seed-index declarations do not include the CLI parser. Application-only
configuration conversion stays in the application layer; forward declarations
avoid pulling command parsing into lower-level interfaces. A source file's
presence in a checkout does not make it a supported API.

## Build and public tests

```bash
cmake --preset testing
cmake --build --preset testing --parallel 4
ctest --preset testing --output-on-failure --no-tests=error
```

Python 3 drives small artificial-input CLI and format tests. They validate
reference/query coordinates, actual CIGAR pairs, score/consumption consistency,
index persistence, gzip, cancellation, rollback and installed library consumers.
They do not download genomes or measure real-data accuracy. CI tests Release
builds with the centromere backend on and off, and treats zero registered tests
as an error. Internal research tests, when present, are registered separately.

For memory checks, use a separate build with
`-DRAMAG_ENABLE_SANITIZERS=ON` and run the relevant small tests. Do not reuse a
production cache containing retired experimental options.

## Behavioral contracts

Preserve sequence ownership and stable indices during seed/candidate processing.
Returned records own their CIGAR storage. Validate reference and query consumption,
original-forward half-open coordinates and reverse orientation before publishing.
Do not merge a failed gap using an invented I+D patch.

Preserve scoring, tie ordering, selection, exact-match semantics and the existing
resource limits during refactoring. The genome affine score and centromere
three-piece score are distinct. Unknown bases are never exact anchors. Copy
branches, unknown directions and unresolved gaps must remain explicit.

For retained behavior, compare ordered output records and independently expanded
coordinate pairs, not only total coverage or record count. Command-path provenance
may differ; source offsets, directions and CIGAR may not be silently discarded.
Time budgets can affect a computation, so report any budget-sensitive difference
instead of widening the budget to force a test to pass.

## Logs and execution

Each command writes a run log under its work directory. `--print-effective-config`
prints `key=value` text and exits without alignment. It is not a JSON API.
Stage capacities and RSS observations have different scopes and must not be added
together as a process peak. `--threads` is a worker budget, not a memory cap.

SIGUSR1 requests progress on POSIX systems; SIGINT/SIGTERM request cooperative
cancellation. Existing outputs and independently completed saved indexes remain
protected. A forced kill is outside the handled rollback guarantee. There is no
alignment-resume interface.

## Troubleshooting and contributions

Start with [installation](installation.md), then check input validity, supported
formats, fresh output paths, pinned dependency state and the full run log.
A genome-only binary cannot run `align -c`. BED mixed alignment is not supported.
When reporting a problem, include the version, exact command, OS/compiler,
relevant log and a small reproducible input if shareable. Remove private paths
and sequence data before posting. See [CONTRIBUTING](../CONTRIBUTING.md).

## Prepare a product source candidate

`public-source-files.txt` is the explicit export list. New files need to be
classified before inclusion. The export helper copies those files to a new
directory and writes a sibling SHA-256 manifest:

```bash
python3 scripts/prepare_public_source.py --output /path/to/new-candidate
```

It never rewrites versions, updates Git branches, commits or pushes. Internal
tests, experiment reports, benchmark harnesses, plans and local dependency
checkouts are excluded. A candidate must independently configure, build, test
and install using only its own contents before publication. Public headers and
CMake source/include/install paths must all be represented in the export list.

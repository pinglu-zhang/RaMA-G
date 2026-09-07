# RaMA-G

RaMA-G is a deterministic C++20 aligner for pairwise whole-genome and
assembly-to-assembly nucleotide alignment. A run has an explicitly directed
reference multi-FASTA and query multi-FASTA, constructs exact seeds on both
query strands, chains and extends them, resolves the resulting alignment
records, and publishes validated output files as one transactional result set.

RaMA-G currently provides:

- plain FASTA and sequential gzip/BGZF-compatible FASTA input;
- reusable Sufkit full suffix-array reference indexes;
- `fast`, `mumreference`, strict `mum`, generalized `smem`, and `maxmatch`
  seed modes;
- `all` and reciprocal interval-based `one-to-one` alignment selection;
- SAM 1.6, PAF, MUMmer-compatible delta, pairwise MAF v1, and UCSC chain;
- deterministic multi-threaded stages under one OpenMP worker budget;
- terminal or log-friendly progress, `SIGUSR1` status snapshots, and safe
  `SIGINT`/`SIGTERM` interruption;
- provenance manifests, output validation, no-replace publication, rollback,
  and a final completion marker.

## Documentation

- [User guide](docs/user-guide.md): installation, commands, inputs, indexes,
  seed and selection modes, outputs, long-running jobs, result integrity, and
  troubleshooting.
- [Developer guide](docs/developer-guide.md): implemented algorithms, data
  model, subsystem boundaries, concurrency, persistence, publication, and
  correctness invariants.

## Build

RaMA-G's production baseline is Linux x86-64. It requires CMake 3.22 or newer,
a C++20 compiler, Git, zlib, and OpenMP support for both C and C++.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

The public source tree omits the project's internal test suite. CMake detects
that omission and skips test targets; it does not prevent configuring, building,
installing, or running RaMA-G.

The default build resolves exact, clean Sufkit and SeqPro Git snapshots. See
the [user guide](docs/user-guide.md#2-installation-and-build) for fresh release,
offline-development, and dependency-override details.

## Quick start

```bash
./build/ramag align \
  --reference reference.fa.gz \
  --query query.fa.gz \
  --output-prefix result \
  --work-dir work \
  --threads 8
```

The legacy prefix interface defaults to exactly `result.sam`, `result.paf`,
and `result.delta`. A successful run also publishes `result.manifest.json` and
finally `result.complete`. The completion marker, not an individual alignment
file or process exit message, is the authority that the requested result set
was validated and published.

Any supported subset can instead be selected by repeating `--output`:

```bash
./build/ramag align \
  --reference reference.fa.gz \
  --query query.fa.gz \
  --output result.paf \
  --output result.maf \
  --output result.chain \
  --work-dir work \
  --threads 8
```

All explicit output paths must have the same directory and base prefix.

## Current boundaries

RaMA-G does not currently provide an automatic comparison/report pipeline
beyond `--selection-mode all|one-to-one`. It does not call variants, calculate
a coverage report, emit BAM/CRAM, VCF, compressed alignment output, or HTML,
and it does not automatically discover an index. Gzip input is sequential and
does not use `.gzi`. Interrupted runs are safely rolled back but are not
resumable; a previously completed explicit reference index can be reused.
MAPQ is currently reported as unknown (`255`), and strict MUM uniqueness is
scoped to one complete query record at a time.

No unqualified speed or accuracy claim against MUMmer4 is made without a
separately accepted, version-bound benchmark.

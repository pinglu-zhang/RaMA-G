# RaMA-G

RaMA-G is a C++20 library and command-line application for deterministic pairwise
whole-genome nucleotide alignment. It aligns one reference multi-FASTA with one
query multi-FASTA using a complete Sufkit suffix array and a KSW2 pairwise core.

- Plain and gzip/BGZF-compatible sequential FASTA input through kseq and zlib.
- Five seed modes: `fast`, `mumreference`, `mum`, `smem`, `maxmatch`.
- `all` or pairwise bilateral `one-to-one` selection.
- SAM, PAF, delta, pairwise MAF and UCSC chain from the same final records.
- Optional persistent compressed-LCP reference index, progress and signal handling.
- Thread-safe console/file logging using spdlog.

See the [user guide](docs/user-guide.md) for commands and input/output contracts,
and the [developer guide](docs/developer-guide.md) for the algorithms and library.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 8
./build/ramag align --reference reference.fa --query query.fa \
  --output result.paf --work-dir work --threads 8
```

An exit code of zero means all requested outputs were validated and published.
Each CLI run writes `work/runs/<run-id>/run.log`. Alignment and index commands do
not generate JSON manifests, SHA-256 digests or completion marker files.
Ordinary alignment builds its reference index in memory; use `--save PATH` to
persist it, or `--reference-index PATH` to reuse an existing index.

Sufkit is fetched at an exact commit during configuration. KSW2, kseq and spdlog
are vendored. The `cmake/` directory is required source, not a build artifact.
Production builds do not require the internal `tests/` or `plan/` directories.

RaMA-G currently provides no checkpoint/resume, calibrated MAPQ, compressed
alignment output, variant calling or automatic genome-comparison report.
Multi-file publication rolls back on handled failures; a forced process kill
can leave a subset of output files. Check the process exit code before using them.
No unqualified performance or accuracy claim against MUMmer4 is made.
Historical design and experiment material, when present under `plan/`, describes
the state at its recorded date; the two guides describe current behavior.

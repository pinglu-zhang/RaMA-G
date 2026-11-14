# RaMA-G

RaMA-G aligns assembled nucleotide sequences. It provides a command-line program
and C++20 libraries for pairwise whole-genome and centromere alignment.

- **Whole genomes:** search both orientations, select a bilateral one-to-one
  main path, link compatible records, and fill supported bounded gaps.
- **Centromeres:** predict directions, find recursive RaMA anchors and accepted
  copy branches, review paths, and align supported internal gaps.
- **Outputs:** SAM, PAF, pairwise MAF, delta and UCSC chain for genome alignment;
  path-aware SAM, PAF and MAF with diagnostic tables for centromeres.
- **Inputs:** plain or gzip FASTA; an optional persistent reference index can be
  reused across queries.

## Install

The production platform is **Linux x86-64 with SSE4.2**. Build requirements are
C++20 and C11 compilers, CMake 3.22 or newer, Git, zlib development files, C/C++
OpenMP, and Make or Ninja. On Ubuntu 24.04:

```bash
sudo apt-get update
sudo apt-get install build-essential cmake git zlib1g-dev

git clone https://github.com/pinglu-zhang/RaMA-G.git
cd RaMA-G
cmake --preset release
cmake --build --preset release --parallel 4
./build/release/ramag --version
```

The first configuration downloads pinned Sufkit and WFA2 sources. Rust, Cargo,
Python and an independent Centrolign installation are not needed to build or run
the program. For a genome-only build, use the `genome-only` preset; it omits the
centromere backend and WFA2, but still requires Sufkit.

To install into a user-owned prefix:

```bash
cmake --install build/release --prefix "$HOME/.local"
```

See [installation](docs/installation.md) for offline dependencies, build options
and installation checks.

## Run your first alignment

Run from the repository root. The included inputs are small artificial sequences:

```bash
./build/release/ramag align \
  --reference examples/genome/reference.fa \
  --query examples/genome/query.fa \
  --output results/first.paf \
  --work-dir results/work \
  --threads 2
```

Exit status zero means the requested output was validated and published. Open
`results/first.paf` to inspect the mappings and `results/work/runs/*/run.log` for
stage statistics. Use a new output prefix for another run: existing outputs are
never overwritten.

## Common tasks

### Save and reuse a reference index

```bash
./build/release/ramag index \
  --reference examples/genome/reference.fa \
  --output results/reference.sufidx \
  --work-dir results/index-work

./build/release/ramag align \
  --reference examples/genome/reference.fa \
  --reference-index results/reference.sufidx \
  --query examples/genome/query.fa \
  --output results/reused.paf \
  --work-dir results/reused-work
```

### Align a centromere pair

Each FASTA must contain exactly one record:

```bash
./build/release/ramag align \
  -c \
  --reference examples/centromere/reference.fa \
  --query examples/centromere/query.fa \
  --output-prefix results/centromere \
  --threads 1
```

This writes SAM, PAF and MAF plus direction, anchor, copy, segment and gap TSVs.
Only supported paths and gaps are aligned. Unresolved regions and ambiguous
copies remain explicit; a completed run does not promise whole-input coverage.

### Process several queries

```bash
./build/release/ramag batch \
  --reference examples/genome/reference.fa \
  --seqfile examples/batch/queries.txt \
  --output-dir results/batch \
  --threads 2
```

**Batch retains candidate-output semantics.** It shares a reference index and
processes queries serially, but is not a repeated invocation of the linked
one-to-one `align` workflow. Results are listed in `results/batch/batch.tsv`.

## Choose an output

| Entry point | Default formats | Other supported formats |
|---|---|---|
| `align --output-prefix PREFIX` | SAM, PAF, delta | MAF, chain |
| `align --output FILE.paf` | Requested suffix only | SAM, MAF, delta, chain |
| `align -c --output-prefix PREFIX` | SAM, PAF, MAF and TSV diagnostics | Select standard formats with `--formats` |
| `batch` | PAF | SAM, MAF, delta, chain |

PAF includes `cg` CIGAR. Centromere paths can share physical coordinates when an
accepted copy relationship is present. A SAM primary/secondary designation is
not a complete description of that relationship. See the [output contracts](docs/centromere-alignment-formats.md).

## Documentation

- [Installation](docs/installation.md)
- [Quick start](docs/quick-start.md)
- [Commands, inputs and outputs](docs/user-guide.md)
- [Algorithms](docs/algorithms.md) and [C++ libraries](docs/library.md)
- [Development, testing and troubleshooting](docs/developer-guide.md)
- [Migration from experimental interfaces](docs/migration.md)

## Check a build

Public tests use small synthetic inputs, not real genomes or benchmark datasets.
Python 3 is needed for these tests; samtools is optional:

```bash
cmake --preset testing
cmake --build --preset testing --parallel 4
ctest --preset testing --output-on-failure --no-tests=error
```

## Scope and attribution

RaMA-G is an assembly aligner, not a short-read mapper or variant caller. It does
not provide calibrated MAPQ, compressed alignment output or alignment checkpoint
resume. Handled publication failures roll back the current outputs; an unhandled
process kill can leave a subset of files, so always check the exit status.

RaMA-G is distributed under the [MIT license](LICENSE). Source-derived RaMA,
KSW2, Sufkit, WFA2, Centrolign kernels and other dependencies retain their own
licenses and attribution; see [third-party notices](THIRD_PARTY_NOTICES.md).
Please report reproducible problems through the repository's GitHub issues.
See [contributing](CONTRIBUTING.md) before proposing changes.

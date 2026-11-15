# RaMA-G

RaMA-G aligns assembled nucleotide sequences, including whole genomes and
centromere pairs. It reads plain or gzip FASTA and writes base-level alignments.

## Install

### Conda (recommended)

Install the Linux x86-64 package from [malab](https://anaconda.org/malab/rama-g):

```bash
conda create -n ramag --override-channels \
  -c malab -c conda-forge rama-g=0.1.1
conda activate ramag
ramag --version
ramag align --help
```

The package includes the complete centromere backend and requires **glibc 2.17
or newer and SSE4.2**. Conda supplies the C++ runtime, OpenMP and zlib; a compiler
and Python are not needed to run RaMA-G. Example inputs are installed under
`$CONDA_PREFIX/share/ramag/examples`. In the examples below, Conda users can
replace `./build/release/ramag` with `ramag` and use the installed input paths.

### Build from source

Supported platform: **Linux x86-64 with SSE4.2**. Build requirements are a
C++20/C11 compiler, CMake 3.22+, Git, zlib, OpenMP, and Make or Ninja.
On Ubuntu 24.04:

```bash
sudo apt-get update
sudo apt-get install build-essential cmake git zlib1g-dev
git clone --depth 1 https://github.com/pinglu-zhang/RaMA-G.git
cd RaMA-G
cmake --preset release
cmake --build --preset release --parallel 4
./build/release/ramag --version
```

The first configuration downloads pinned Sufkit and WFA2 sources. Python and a
separate Centrolign installation are not required to build or run RaMA-G.
For a genome-only build, replace `release` with `genome-only`; Sufkit remains
required. See [installation](docs/installation.md) for offline builds and details.

## First alignment

Run from the repository root. Included FASTA files are small artificial examples.

```bash
./build/release/ramag align \
  --reference examples/genome/reference.fa \
  --query examples/genome/query.fa \
  --output results/first.paf \
  --work-dir results/first-work --threads 2
```

Use a new output prefix for each run: existing outputs are never overwritten.
Exit status zero indicates success. Logs are under `WORK/runs/*/run.log`.

## Common commands

Save a reference index and reuse it:

```bash
./build/release/ramag index \
  --reference examples/genome/reference.fa \
  --output results/reference.sufidx --work-dir results/index-work
./build/release/ramag align \
  --reference examples/genome/reference.fa \
  --reference-index results/reference.sufidx \
  --query examples/genome/reverse.fa \
  --output results/reverse.paf --work-dir results/reverse-work
```

Align one centromere record per input:

```bash
./build/release/ramag align -c \
  --reference examples/centromere/reference.fa \
  --query examples/centromere/query.fa \
  --output-prefix results/centromere --threads 1
```

This produces SAM, PAF, MAF and diagnostic TSVs. Unresolved regions and ambiguous
copies remain explicit; successful output does not imply whole-input coverage.

Process several queries with one shared reference index:

```bash
./build/release/ramag batch \
  --reference examples/genome/reference.fa \
  --seqfile examples/batch/queries.txt \
  --output-dir results/batch --threads 2
```

**Ordinary `align` outputs a linked bilateral main path; `batch` retains candidate
alignments.** Batch runs queries serially and reports each outcome in `batch.tsv`.

## Documentation

- [Installation](docs/installation.md)
- [Commands and parameters](docs/user-guide.md)
- [Formats, coordinates and diagnostic tables](docs/formats.md)
- [Troubleshooting](docs/troubleshooting.md)

Genome alignment supports SAM, PAF, pairwise MAF, delta and UCSC chain.
Centromere alignment supports SAM, PAF and MAF with copy/path metadata.
Use `ramag align --help` for options and `--print-effective-config` to inspect a
command without aligning. BAM/CRAM, compressed outputs and alignment resume are
not supported. Please report problems through the repository's GitHub issues.

## License

RaMA-G's original code is MIT-licensed. Bundled and linked components retain
their own licenses and source attribution; see [LICENSE](LICENSE).

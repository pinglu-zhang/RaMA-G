# Installation

## Requirements

RaMA-G supports Linux x86-64 with SSE4.2. Use a C++20/C11 compiler, CMake 3.22 or
newer, Git, zlib development files, C/C++ OpenMP, and Make or Ninja.

On Ubuntu 24.04:

```bash
sudo apt-get update
sudo apt-get install build-essential cmake git zlib1g-dev
```

Python 3 is needed only for optional synthetic installation tests. Samtools can
provide an additional SAM format check. Neither is a core runtime dependency.
Rust, Cargo and a separate Centrolign executable are not required.

## Build

```bash
git clone --depth 1 https://github.com/pinglu-zhang/RaMA-G.git
cd RaMA-G
cmake --preset release
cmake --build --preset release --parallel 4
./build/release/ramag --version
./build/release/ramag align --help
```

Initial configuration fetches exact Sufkit and WFA2 Git revisions. Other required
source components are bundled. The default Release preset does not build tests.

For whole-genome alignment without the centromere backend:

```bash
cmake --preset genome-only
cmake --build --preset genome-only --parallel 4
```

Genome-only supports `align`, `index` and `batch`. It still requires Sufkit;
`align -c` reports an unavailable backend. Use a separate build directory for
each configuration. The complete Linux-64 package is also available from the
`malab` Conda channel; see below. No pip distribution is provided.

## Install

```bash
cmake --install build/release --prefix "$HOME/.local"
"$HOME/.local/bin/ramag" --version
```

Add `$HOME/.local/bin` to your PATH to invoke `ramag` directly. License materials
are installed under `share/ramag`; retained C++ packages and public headers are
also installed. When replacing an older installation, use a clean prefix to
avoid retaining obsolete files from that installation.

## Offline configuration

Prepare clean Git checkouts of Sufkit at
`028075e6f2d622fcbdcf76b153bbde0f069ca64f` and WFA2 at
`99a5879a3a2972eff168cbcc6610b6cf02844ed0`. Supply their absolute paths:

```bash
cmake -S . -B build/offline -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_TESTING=OFF -DRAMAG_BUILD_TESTS=OFF \
  -DRAMAG_SUFKIT_SOURCE_DIR=/path/to/clean/pinned/sufkit \
  -DRAMAG_WFA2_SOURCE_DIR=/path/to/clean/pinned/WFA2-lib
cmake --build build/offline --parallel 4
```

A source ZIP or an arbitrary dependency revision is not a compatible override.
The dependency checkouts must have the exact commits and no uncommitted files.
For an offline genome-only build, add `-DRAMAG_ENABLE_CENTROMERE=OFF`; WFA2 is
then unnecessary. Disabling Sufkit is not an offline installation method.

## Verify an installation

Run the [first alignment](../README.md#first-alignment) with the supplied
artificial FASTA files. Check exit status zero and the resulting PAF.
Optional synthetic tests require Python 3:

```bash
cmake --preset testing
cmake --build --preset testing --parallel 4
ctest --preset testing --output-on-failure --no-tests=error
```

These tests use no real genomes or benchmark downloads. See
[troubleshooting](troubleshooting.md) if a build, dependency check or run fails.

## Conda package

The Linux-64 recipe in `recipe/` builds the full program, including the
centromere backend. It targets glibc 2.17 or newer and x86-64 CPUs with SSE4.2.
It does not use `-march=native` or require AVX. Conda supplies the C++ and OpenMP
runtime libraries and zlib; Python is not a program runtime dependency.
The Conda installation itself must also support the host operating system.
For a new Conda installation, consult the official
[Miniforge requirements](https://github.com/conda-forge/miniforge#requirements-and-installers),
which list glibc 2.17 support for Linux x86-64.

Install from the public [malab channel](https://anaconda.org/malab/rama-g):

```bash
conda create -n ramag --override-channels \
  -c malab -c conda-forge rama-g=0.1.1
conda activate ramag
ramag --version
ramag align --help
```

Use a new environment so Conda can resolve the runtime dependencies. Avoid
`conda install --no-deps` and do not replace system glibc libraries. Small example
inputs are installed in `$CONDA_PREFIX/share/ramag/examples`; call `ramag`
instead of `./build/release/ramag` when following the source-build examples.

For example, run an alignment from any writable working directory:

```bash
ramag align \
  --reference "$CONDA_PREFIX/share/ramag/examples/genome/reference.fa" \
  --query "$CONDA_PREFIX/share/ramag/examples/genome/query.fa" \
  --output results/conda-first.paf --work-dir results/conda-first-work --threads 2
```

If you received a local `.conda` file, you can install its absolute path in place
of `rama-g=0.1.1` with `-c conda-forge` for its dependencies.

To build the local package with Conda Build:

```bash
conda build recipe --override-channels -c conda-forge --no-anaconda-upload
```

The recipe pins the RaMA-G source, Sufkit and WFA2 revisions and a glibc 2.17
sysroot. Building requires network access to fetch the pinned sources and
Conda dependencies. Package construction runs the synthetic public tests;
`--no-anaconda-upload` keeps the package local.

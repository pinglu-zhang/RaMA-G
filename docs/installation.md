# Installation

## Supported production configuration

RaMA-G currently targets Linux x86-64. The CPU must support SSE4.2, which is used
by the pinned Sufkit dependency. Use a compiler supporting C++20 and C11, CMake
3.22 or newer, Git, zlib development files, C/C++ OpenMP, and Make or Ninja.
The pthread and math libraries are supplied by the system toolchain.

For Ubuntu 24.04:

```bash
sudo apt-get update
sudo apt-get install build-essential cmake git zlib1g-dev
```

Python 3 is needed only for public tests and the optional diagnostic script.
Samtools is an optional independent SAM validator. Rust, Cargo, Doxygen and
MUMmer4 are not production dependencies. No Conda or pip distribution is provided.

## Build from source

```bash
git clone https://github.com/pinglu-zhang/RaMA-G.git
cd RaMA-G
cmake --preset release
cmake --build --preset release --parallel 4
./build/release/ramag --version
./build/release/ramag --help
```

The `release` preset builds the complete Linux program without test-only tools.
Its explicit equivalent is:

```bash
cmake -S . -B build/release -DCMAKE_BUILD_TYPE=Release   -DBUILD_TESTING=OFF -DRAMAG_BUILD_TESTS=OFF
cmake --build build/release --parallel 4
```

`genome-only` uses a separate build directory and disables the centromere backend:

```bash
cmake --preset genome-only
cmake --build --preset genome-only --parallel 4
```

A genome-only binary supports genome `align`, `index` and `batch`. A centromere
request reports an unavailable dependency; it never falls back to genome mode.

## Pinned dependencies and offline configuration

Initial configuration fetches exact Sufkit and WFA2 Git revisions over HTTPS.
KSW2, RaMA-derived code, the minimal Centrolign stitching kernels, kseq and spdlog
are already included. Genome-only omits WFA2, but still needs Sufkit.

For offline use, prepare clean Git checkouts of the exact revisions recorded in
`CMakeLists.txt` and the dependency CMake modules, including their required source
files. Supply absolute paths:

```bash
cmake \
  -S . \
  -B build/offline \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_TESTING=OFF \
  -DRAMAG_BUILD_TESTS=OFF \
  -DRAMAG_SUFKIT_SOURCE_DIR=/path/to/clean/pinned/sufkit \
  -DRAMAG_WFA2_SOURCE_DIR=/path/to/clean/pinned/WFA2-lib
cmake --build build/offline --parallel 4
```

An arbitrary source ZIP is not equivalent to a clean pinned Git checkout.
Do not change dependency commits to bypass validation. Disabling production
adapters is not an offline installation method.

## Install and verify

```bash
cmake --install build/release --prefix "$HOME/.local"
"$HOME/.local/bin/ramag" --version
```

Add the prefix's `bin` directory to your shell PATH if desired. Installation also
provides the `RaMAGPairwise` and `RaMAGCentromere` CMake packages. See the
[library guide](library.md) for isolated consumer examples.

## Test a build

```bash
sudo apt-get install python3
cmake --preset testing
cmake --build --preset testing --parallel 4
ctest --preset testing --output-on-failure --no-tests=error
```

These tests use artificial sequences and require no benchmark downloads. A zero
registered-test count is an error. Backend-OFF testing uses a separate directory:

```bash
cmake \
  -S . \
  -B build/testing-genome \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_TESTING=ON \
  -DRAMAG_BUILD_TESTS=ON \
  -DRAMAG_WARNINGS_AS_ERRORS=ON \
  -DRAMAG_ENABLE_CENTROMERE=OFF
cmake --build build/testing-genome --parallel 4
ctest --test-dir build/testing-genome --output-on-failure --no-tests=error
```

## Troubleshooting

- **OpenMP or zlib not found:** install compiler and development packages, then
  configure a new build directory.
- **Dependency revision/dirty checkout error:** restore the required dependency
  revision and clean its tracked changes; do not disable validation.
- **Old experimental cache rejected:** create a fresh directory and use a current
  preset. Retired configuration names are listed in [migration](migration.md).
- **Illegal instruction:** check CPU support for SSE4.2.
- **Output already exists:** use a new output prefix; RaMA-G does not overwrite.
- **Large memory use:** reference indexing, query sequences and candidates consume
  memory. `--threads` controls workers, not a whole-program memory cap.

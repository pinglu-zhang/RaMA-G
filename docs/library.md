# C++ libraries

Install the project before compiling a consumer. Public include paths and CMake
package names for retained capabilities remain stable. Headers private to the
CLI, workflow orchestration or review implementation are not installed APIs.

## Pairwise library

```cmake
cmake_minimum_required(VERSION 3.22)
project(MyAlignment LANGUAGES CXX)
find_package(RaMAGPairwise CONFIG REQUIRED)
add_executable(my-alignment main.cpp)
target_link_libraries(my-alignment PRIVATE RaMAG::pairwise)
```

```bash
cmake -S examples/pairwise -B build/pairwise-example   -DCMAKE_PREFIX_PATH="$HOME/.local"
cmake --build build/pairwise-example
./build/pairwise-example/pairwise-example
```

`AlignPairwiseCore` accepts sequence records and exact canonical seeds and returns
owned records with reference/query selection information. Input storage must
remain alive for the duration of the call. Coordinates are zero-based, half-open
on the original forward inputs. `LinkSelectedAlignmentRecords` consumes an
already selected disjoint main path; it validates actual CIGAR consumption and
preserves the input transaction on cancellation or invalid input.

The high-level `AlignmentOptions` retains `All`, `ReciprocalOneToOne` and
`LinkedOneToOne` semantics for explicit callers. The library default is not
implicitly changed to match the ordinary CLI. `LinkGapMode` remains a library
choice. Optional vacant/chain recovery interfaces have been removed; see
[migration](migration.md).

## Centromere library

```cmake
find_package(RaMAGCentromere CONFIG REQUIRED)
target_link_libraries(my-alignment PRIVATE RaMAG::centromere)
```

Installed headers provide source-derived RaMA anchors, direction evidence,
directed recursive anchors, copy relationships, path stitching and format APIs.
`CentromereBackendAvailable()` checks computation support. A backend-OFF build
retains the public contract and reports `DependencyError` for unavailable
computation rather than silently selecting another algorithm.

`AlignCentromereAnchoredPaths` consumes frozen primary/copy anchor results and
aligns eligible internal gaps. It does not predict direction or resolve an
ambiguous copy. Its result owns segments, gap states and diagnostics. The complete
CLI adds direction planning and coordinated review stages around these library
operations; one low-level call is not a replacement for that workflow.

## Ownership, callbacks and installation

Callbacks may run on computation workers and must be thread-safe. Borrowed views
and observation references remain valid only for the documented call. Caller
cancellation propagates; do not convert it into successful partial publication.
Independent invocations do not share mutable sequence ownership.

Public tests install into an isolated prefix and compile consumers using only
installed headers and exported targets. Use the same approach for applications;
do not add this repository's private source directories to a consumer's includes.
Changing removed experimental structures requires recompilation, not just
replacement of a static library archive.

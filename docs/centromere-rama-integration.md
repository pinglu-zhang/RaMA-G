# Single-pair RaMA integration

## Identity and scientific scope

Current CLI default: `RAMAG_CENTROMERE_ANCHORS_ONLY=ON`. `align -c` calls
`FindCentromereAnchors`, exports initial/final TSVs and returns before interval
partitioning or WFA. `AlignCentromeres` remains a separate complete library API;
the complete-path description below is not a claim that anchor TSVs are full
alignments. No WFA execution is required for the anchor-only export tests.

This backend uses the rules of
[pinglu-zhang/RaMA](https://github.com/pinglu-zhang/RaMA/tree/7bc89c4902bc9254f18cc4b05c4bf491a10bd7cb)
at **7bc89c4902bc9254f18cc4b05c4bf491a10bd7cb**.
The untouched source snapshot is in `third_party/rama/upstream`; the
source-derived internal component is in `third_party/rama/core`.
WFA2 is fixed at **99a5879a3a2972eff168cbcc6610b6cf02844ed0**.
Its clean checkout is verified; four safety-patched translation units are
separately stored in `third_party/rama/wfa_safety`.
Original copyright notices and licenses are retained.

This is **original RaMA rules integrated with necessary safety repairs**, not
an assertion that every compiled source byte is unchanged. It is a complete
alignment constrained by recursively selected RaMA anchors, not an inversion
caller or an unconstrained whole-pair optimum.

| Original component | Integrated behavior |
|---|---|
| `Anchor/gsacak.c`, `RMQ.cpp` | Generalized suffix array, LCP, document array, ISA and RMQ |
| `RareMatchFinder::findRareMatch` | Original LCP windows, first qualifying rare level, maximum level 100 |
| Map grouping and left expansion | Same minimum-position/length grouping and occurrence combinations on ACGT |
| `findOptimalPairs` | Integer occurrence weights; original score, floating calculations, predecessor order and ties; no k cap |
| `AnchorFinder::locateAnchor` | Original recursive intervals, traversal and final adjacent-anchor merge |
| `RareMatchPairs2Intervals` | Original leading/inter-anchor/trailing partition; no reinterpretation of nested tail rules |
| `PairAligner::alignIntervals` | Empty-side I/D; original <=5 versus >100 heuristic; WFA2 for other intervals |
| WFA2 | match 0, mismatch 3, gap open/extension 4/2 and 12/1; medium memory; inherited adaptive 10/50/1 |
| Writers | Existing RaMA-G common writers and transactional publication, not RaMA's original format exporters |

The source default WFA heuristic matters. For `A^100` versus `T^100`,
the pinned default gives cost **225**, while unheuristic WFA and independent
dual-affine DP give **224**. The integration preserves 225. Other recorded
100/101, 101/100 and 101/101 boundaries show the same one-unit gap. Tests
retain these facts rather than silently disabling the heuristic.

## Safety repairs, separate from algorithm changes

1. **LCP workspace initialization.** The original complete constructor used
   uninitialized malloc memory. Sparse PLCP construction reads initial LCP
   workspace entries. Value-initialized storage is checked against an
   independent generalized suffix/LCP oracle, including the original 455 bp
   crash fixture; successful execution alone is not the correctness check.
2. **Signed empty sentinel.** Two recursive induced-sorting loops previously
   subtracted one from the signed EMPTY sentinel before checking eligibility.
   The non-positive sentinel check now precedes subtraction.
3. **Ownership and tree traversal.** Per-instance RAII index storage, bounded
   tree insertion/indexing, captured worker exceptions and drain-before-unwind
   replace raw ownership and dropped packaged-task exceptions. Source `exit`
   becomes exceptions. There is no host logger or global signal registration.
4. **WFA unaligned loads/stores.** Scalar 64-bit character accesses now use
   memcpy with the same eight-byte semantics. An empty packed trace no longer
   shifts a 32-bit value by 32. These were detected under UBSan, not inferred
   from the function names. No sanitizer suppression is used.
5. **WFA allocation failure.** Build-local C allocation/exit hooks track
   allocations per worker. Errors unwind through C compiled with exception
   unwind support; residual allocations are released by a per-call guard.
   Interruption probes also occur at allocation checkpoints. A single gSACA-K
   construction or a long allocation-free kernel is not preempted mid-instruction.
6. **Representability.** Public coordinates remain 64-bit. This pinned
   32-bit backend rejects pairs with combined length above
   `(INT_MAX-64)/3` (715,827,861 bp), conservatively protecting offsets,
   score arithmetic and internal sentinel headroom. No coordinate truncation.
   Output uses 64-bit run lengths rather than the original 28-bit packed runs.

N adaptation is deliberate: normalized N maps to distinct internal bytes on
the two sides, LCP is capped at canonical runs, and left expansion stops at
unknowns. Paired N/N and N/base columns cost 3 and become X. Sequence data
written to files always uses the original normalized bases. The original
random-N preprocessing is not reproduced.

## Library and lifecycle

```cmake
find_package(RaMAGCentromere CONFIG REQUIRED)
target_link_libraries(my_program PRIVATE RaMAG::centromere)
```

```cpp
#include <ramag/centromere.hpp>
ramag::CentromereAlignmentOptions options;
options.worker_threads = 4;
auto result = ramag::AlignCentromeres(reference, query, options);
```

Inputs are single, nonempty, normalized A/C/G/T/N `SequenceRecord` objects;
numeric IDs are preserved. Progress callbacks execute on the calling thread.
Interruption callbacks can execute on workers, must be thread-safe and may
throw; their lifetime ends before the call returns. Concurrent independent
calls have independent state.

Anchor recursion and gap alignment use separate phases of the supplied thread
budget. No nested OpenMP teams are introduced. Index/RMQ/recursion ownership
ends before gap workspaces are allocated. Gaps carry coordinates; at most 64
gap results reside in a batch and only active workers copy interval text.
Results are consumed in interval order, independently of task completion.
There is no arbitrary occurrence truncation, sampling or predecessor cap.

The exact kernel and installed library are Linux production targets. Other
platforms default to a disabled-backend stub. `RAMAG_ENABLE_CENTROMERE=OFF`
also permits genome-only builds without downloading WFA2.

## Direction-region API

The installed headers expose `FindCentromereDirectionRegions` through
`<ramag/centromere_blocks.hpp>`. It reports query direction evidence and unknown
regions; it does not establish reference-copy correspondence or generate CIGAR.
Direction planning and local direction repair are separate from the default
anchor-only CLI route. Experimental direction features are disabled by default.
See [direction configuration](centromere-directions.md) for opt-in settings.

## Validation scope

The integration is checked with independent suffix/LCP oracles, source-derived
anchor differential tests, coordinate and lifecycle tests, and installed-library
consumers. The public source distribution does not include the internal test
harness or experimental datasets. Anchor correctness and direction labels are
separate from complete alignment accuracy.

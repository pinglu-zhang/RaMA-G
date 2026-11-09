# Private Centrolign stitching kernels

Pinned upstream: https://github.com/jeizenga/centrolign at
`a604cbdcff2e22e643fc9e17d4679296ee5e42fe`. MIT copyright and complete permission
notice are retained in LICENSE. source-manifest.json hashes the original inputs;
local-patches.json inventories local edits and new instrumentation.

Only the base graph, unpruned graph-to-graph wavefront, and full-matrix PO-POA
kernels are called by RaMA-G. The caller supplies already oriented linear gaps
between two accepted anchors. No upstream anchoring, stitching heuristic,
pruned WFA, long-deletion shortcut, or cyclic graph flow is invoked.

Local changes isolate names/includes, replace the instantiated dynamic
containers with per-call counting allocators, and add cooperative interruption
checks. The numerical recurrences, ordering and scoring parameters are unchanged.
The recursive Euclidean gcd function was rewritten as its equivalent iterative
form to avoid a hidden std::function allocation.

The cap counts actual requested allocation payload bytes for simultaneously
live graph storage, matrix/frontier/queues, hash buckets and nodes, traceback,
and public CIGAR capacity during conversion. Allocator bookkeeping, process RSS,
caller-owned sequences, callbacks, and returned output after leaving the kernel
are outside this per-gap scope. No global operator new or process allocator is
replaced. Each allocator remembers its originating context for exception-safe
unwinding; caller exceptions retain their original type.

The adapter independently validates every emitted coordinate and extended CIGAR,
fully consumes both strings, and recomputes the score from same-direction gap
runs. It intentionally does not use upstream score_alignment/rescore or the
pure_deletion score helper. N in opposite input graphs receives distinct labels,
so N/N is a mismatch rather than an exact match. Empty sides are handled by the
public assembler, not by the graph kernel.

DP uses upstream int32 cells with an explicit conservative numeric-range check.
Unpruned GGWFA uses its int64 cost/conversion range and is not restricted by the
DP scalar. Time and memory exhaustion are reported as unresolved resources;
there is no heuristic successful fallback.

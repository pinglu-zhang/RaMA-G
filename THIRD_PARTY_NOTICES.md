# Third-party notices

RaMA-G is distributed under the MIT License in [LICENSE](LICENSE).
The root license retains both RaMA-G contributor and Pinglu Zhang copyright.

## Sufkit

The separately maintained library https://github.com/malabz/sufkit is fixed to
`f8c4c386ee883e45ad0f973efc4c8e1148b0068a` and obtained by CMake. Its own license
applies when fetched/linked. RaMA-G does not directly depend on SeqPro.

## KSW2

Vendored `third_party/ksw2/ksw2.h` and `ksw2_extz2_sse.c` retain their upstream
MIT notices; see [KSW2 license](third_party/ksw2/LICENSE).

## kseq and spdlog

`third_party/kseq/kseq.h` retains Attractive Chaos's MIT license in its header.
spdlog headers, bundled fmt and license are in `third_party/spdlog/`; embedded fmt
notices are retained. See [spdlog license](third_party/spdlog/LICENSE).
These exact snapshots were taken from RaMAx revision
`7d08359e0df7f7e6ffcfe67217c3399761cb2129`: `include/kseq.h` and
`third_party/spdlog/include/` respectively. No global/system fmt is substituted.

## Pairwise source history

The pairwise core was extracted from RaMAx revision
`7d08359e0df7f7e6ffcfe67217c3399761cb2129`, copyright 2026 Pinglu Zhang, MIT.
Sources included the anchor, cluster, connection, alignment and CIGAR modules,
plus the pairwise rare-aligner path before graph insertion. The extraction uses
checked 64-bit public coordinates, immutable configuration, in-memory sequence
views and RaMA-G naming. Signed coordinate subtraction and the final floating
DP comparison were corrected. Graph construction and multi-species orchestration
were excluded. Subsequent RaMA-G ownership, recovery and guarded gap-fill changes
are maintained in this repository's history.

The previous separate attribution directory and file-hash inventory were removed.
Required notices are consolidated here and in the root license; third-party
library licenses remain with each library.

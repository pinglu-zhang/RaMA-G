# RaMA and stitching integration

RaMA-G integrates source-derived RaMA anchor search with checked sequence access,
explicit ownership and resource handling. The unchanged upstream source and
license are retained under `third_party/rama/upstream`; project adaptation lives
separately. Sufkit and WFA2 are pinned to the revisions required by CMake.
See [third-party notices](../THIRD_PARTY_NOTICES.md).

## Complete CLI versus low-level APIs

The Linux default `align -c` predicts directions, obtains recursive/copy anchors,
reviews paths and stitches their supported internal gaps. It does not stop after
anchor search. The current stitching score is the Centrolign three-piece model,
with minimal source-attributed kernels under `third_party/centrolign_stitch`.
It does not invoke the full Centrolign executable, graph merger or anchor search.

Low-level RaMA-derived `FindCentromereAnchors` returns sparse exact anchors without
WFA. `AlignCentromeres` is the supplied-orientation source-derived pair alignment
API. Its original RaMA/WFA score and source behavior are separate from complete
CLI direction planning and three-piece path stitching. Do not infer CLI scoring
or inversion support from one low-level API result.

## Safety and builds

The sequence callback verifies access bounds before reading. Anchor ownership,
recursive half-open coordinates and CIGAR consumption are checked. Gap solvers
keep their documented budget and interruption contracts. Reverse coordinates are
returned on the original FASTA; ambiguous bases are never exact anchors.

`RAMAG_ENABLE_CENTROMERE=OFF` builds an unavailable-backend stub. Public headers
remain usable, while unavailable computation raises `DependencyError` instead of
silently selecting another aligner. The complete Linux build compiles the needed
WFA2 C sources without importing its global build settings. No Rust or SeqPro
runtime is required.

Source attribution, upstream manifests and local patch records remain with the
vendored dependencies. Directory cleanup does not rewrite upstream checksums.
Public synthetic tests independently check coordinates, CIGARs, formats and
installed consumers. Large internal experiment records are not installation tests.

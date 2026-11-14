# Directed recursive anchors

`FindCentromereAnchorsInDirections` consumes directed regions and searches for exact
RaMA anchors in their permitted orientation. The complete CLI first predicts the
regions; explicit library callers may supply their own valid region result.

Each anchor is an actual coordinate correspondence, not the bounding rectangle
between two hits. Recursive boundaries are checked as half-open intervals. N
breaks exact matching; substitutions are not exact anchors. Negative query
coordinates are normalized for search and restored to original-forward FASTA
coordinates in results.

Initial roots and final recursive anchors remain distinct. Endpoint, block,
cross-region and joint-path reviews may replace a supported local proposal; their
updates to initial and final tables remain coordinated. Accepted copy footprints,
trusted internal support and fixed-gap CIGAR evidence are protected where required.
A tie, incomplete context or resource limit does not become a confirmed match.

The CLI's corrected-boundary behavior and reference-window policy are part of its
effective configuration. Low-level source-derived APIs keep their documented
boundary semantics; do not assume every library call is the whole CLI workflow.

Anchor outputs are sparse exact correspondences, not a complete base alignment.
Only the stitcher emits internal-gap CIGAR. Unanchored tails and unknown regions
remain in diagnostics. See [stitching](centromere-anchor-stitch.md) and
[copy anchors](centromere-copy-anchors.md).

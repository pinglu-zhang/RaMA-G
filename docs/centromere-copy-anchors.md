# Copy-aware anchors

The centromere workflow keeps its primary anchor paths and discovers additional
copy-supported branches separately. Different branches may share reference or
query coordinates. An accepted branch is an actual correspondence; an ambiguity
group is a set of alternatives, not a confirmed prediction.

Evidence combines within-sequence repetition, projection through fixed primary
anchors, local exact-search verification, direction/order consistency and
independent context. An exact repeated string alone does not establish a unique
reference placement or evolutionary copy history. Competing equal support remains
ambiguous. Newly inferred branches do not become roots for another inference round.

The retained public options bound rare levels to 100, deduplicated local candidates
to 3,000 and total candidates to 100,000 by default. Budget exhaustion is recorded,
not reported as absence of duplication. Identifier assignment and output ordering
are deterministic for a fixed input and configuration.

Copy anchors preserve original-forward half-open coordinates and branch identity.
They are stitched independently from the main paths. Only accepted branches enter
formal alignments; ambiguity tables remain diagnostics. The SAM flag alone is
insufficient to distinguish these roles; use the path/branch provenance tags and
TSV records described in the [format contract](centromere-alignment-formats.md).

This layer predicts sequence correspondence. It does not label a historical
copy event, force all possible copy combinations into output, or claim a unique
ancestral-copy assignment where the sequences contain no discriminating evidence.

# Experimental direction-region API

The default `align -c` command exports RaMA anchors in the supplied orientation.
Automatic direction-region planning is an opt-in library capability. It is not
an inversion caller or complete base-pair alignment, and it does not resolve
which reference repeat copy corresponds to every inferred query interval.

Use `find_package(RaMAGCentromere CONFIG REQUIRED)` and link `RaMAG::centromere`.
The header `<ramag/centromere_blocks.hpp>` provides `CentromereBlockOptions`,
`FindCentromereDirectionRegions` and the direction evidence writers.

## Selecting the local direction repair policy

Set `planning_strategy` to `CentromereBlockStrategy::OrientationEvidenceV9`.
For the full local direction repair pipeline, enable the following option fields:

```text
orientation_connected_gaps
orientation_bridge_unknown_gaps
orientation_audit_known_blocks
recover_long_gap_ends
recover_sparse_pivots
recover_short_gap_ends
recover_mixed_gap_ends
recover_sparse_gap_ends
prefer_stronger_sparse_gap_end
unify_real_anchor_ends
review_local_competitive_gaps
review_isolated_opposite_anchors
review_collinear_fragment_connections
review_collinear_window_gaps
review_alternate_original_endpoints
review_neutral_window_gaps
review_region_chains
supplement_region_anchors
review_unified_gaps
review_gap_interior_opposition
review_known_block_counterchains
repair_inherited_directions
```

`repair_drift_consistent_gaps=false` selects the v1 final repair;
`repair_drift_consistent_gaps=true` additionally enables the v2 long-edge rule.
These flags default to false. Invalid prerequisite combinations raise an error.
Configure `worker_threads`, memory limits and search bounds for the workload;
the public defaults are not a claim of benchmark-equivalent performance.

The CMake equivalent of the v2 flag is
`RAMAG_CENTROMERE_REPAIR_DRIFT_CONSISTENT_GAPS`, also OFF by default. CMake
configuration validates the required preceding options. These settings do not
change the default anchor-only CLI into the direction-region route by themselves.

## Results and interpretation

`local_direction_repair_intervals` is the complete final known-direction set.
`orientation_unknown` records the remaining unknown query intervals. Do not
union final intervals with older intermediate direction tables a second time.
Reviews, witnesses and changes preserve the frozen evidence and decision origin;
changes describe the pre-repair-to-final transition, not arbitrary version deltas.
Coordinates are zero-based, half-open, in the original query orientation.

Both repair policies make one final pass over frozen evidence; repaired output
never becomes new evidence. The final repair stage adds no searches and invokes
no WFA. Earlier enabled stages can perform their configured bounded searches.
For long-edge admission, source reference/query counts of 1/1 describe LCP
source-interval counts, not a proof of global uniqueness in either sequence.
Direction inference does not claim exact matching across inferred gaps.

See [input contract](centromere-interface.md) and
[source integration](centromere-rama-integration.md) for input normalization,
thread ownership and the separately retained complete alignment API.

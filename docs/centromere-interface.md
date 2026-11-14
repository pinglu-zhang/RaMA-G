# Centromere command-line interface

On Linux with the computation backend enabled, `ramag align -c` treats both
FASTA inputs as complete centromere sequences. Each must contain exactly one
nonempty record. The program predicts directions, recursively searches for
anchors, finds accepted copy branches, reviews paths and stitches internal gaps.

```bash
ramag align -c --reference reference.fa --query query.fa   --output-prefix comparison --threads 1
```

Default standard outputs are `comparison.sam`, `comparison.paf` and
`comparison.maf`. Direction, anchor, copy, review, segment, gap and summary TSVs
are also written. `--formats paf` selects only PAF among standard formats;
`--output comparison.paf` is the suffix-driven alternative. Do not mix `--output`
with `--output-prefix`/`--formats`. Delta and chain cannot preserve these path
contracts and are rejected for `-c`.

If omitted, the work directory is `.ramag-work` beside the output prefix. Existing
alignment or diagnostic outputs are not overwritten. Formal outputs are published
only after validation; interrupted snapshots in the work directory are diagnostics.
Check process exit status before consuming results.

The current workflow uses V9 direction planning, corrected recursive boundaries,
copy anchors, endpoint/block/cross-region/joint reviews, supported unplanned-core
recovery and three-piece stitching. These are not independent CLI switches.
`--print-effective-config` prints the resolved `key=value` contract without
performing alignment.

`--reference-index`, `--save`, seed and whole-genome tuning options are not
accepted with `-c`. Backend-OFF builds report dependency error 5 for a computation
request and do not fall back to whole-genome alignment.

The paired BED options `--ref-cen`/`--query-cen` (long aliases
`--reference-centromeres`/`--query-centromeres`) validate zero-based, half-open
annotations. Mixed genome/centromere alignment is not implemented. Both BEDs
are required together, cannot be combined with `-c`, and are not an automatic
region detector. Do not interpret interface validation as a completed alignment.

See [algorithms](algorithms.md), [format contracts](centromere-alignment-formats.md)
and [migration](migration.md).

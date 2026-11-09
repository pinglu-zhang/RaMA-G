# Centromere input interface

On Linux the normal build runs the complete anchored-path workflow with
`ramag align -c`. The command aligns centromere sequences supplied as FASTA.
Each input contains exactly one nonempty record. Plain and gzip FASTA use
the same normalization rules.

```sh
ramag align -c --reference reference_centromere.fa --query query_centromere.fa \
  --output-prefix comparison --threads 4
```

RaMA-G predicts forward and reverse query intervals, searches for recursive
primary anchors and copy branches, reviews endpoints and paths, and fills
supported internal gaps with the `centrolign-3piece` model. In large direction
gaps it first checks the original initial anchors and then searches locally
when the existing evidence is insufficient and the window meets the resource
and context criteria. All coordinates refer to the original input sequences.

The default standard outputs are `comparison.sam`, `comparison.paf` and
`comparison.maf`. Select a subset with `--formats LIST` or provide repeatable
`--output` paths with the same directory and base prefix. The command also
publishes diagnostic TSVs alongside those outputs.

| Suffix | Contents |
| --- | --- |
| `.orientation-local-direction-repair-intervals.tsv` | Direction plan passed to recursive anchor search |
| `.direction-initial-anchors.tsv` | Initial evidence from direction prediction |
| `.initial-anchors.tsv`, `.final-anchors.tsv` | Primary anchors before and after recursion, with the applicable reviews |
| `.copy-anchors.tsv`, `.copy-ambiguities.tsv` | Accepted copy anchors and separately recorded ambiguity groups |
| `.alignment-segments.tsv`, `.alignment-gaps.tsv`, `.alignment-summary.tsv` | Path coordinates, CIGARs, gap status and run totals |
| `.unplanned-core-reviews.tsv` | Local review of direction gaps and their candidate reference locations |

Sequence tails, regions without supported paths, and resource-limited gaps
can remain unresolved. Read the segment coordinates and gap status to determine
the extent of each alignment. See the [format contract](centromere-alignment-formats.md)
and [stitching model](centromere-anchor-stitch.md).

The [direction-region API](centromere-directions.md) also exposes the planner for library use.

Add `--print-effective-config` to inspect the mode, scoring and output formats.
Omitting `--work-dir` places the run logs in `.ramag-work` beside the output prefix.
A build with the centromere backend disabled reports dependency error **5**.

Global genome options (`--seed-mode`, `--selection-mode`, `--min-match`,
`--smem-min-occurrences`, `--min-cluster`, `--max-gap`, `--diag-diff`,
`--diag-factor`, `--break-length`, `--max-dp-cells`, `--save`,
`--reference-index`) are rejected when explicitly combined with `-c`.
The source's rare-level limit is 100, not a predecessor search window.

N is never an exact anchor or an `=` column, even opposite another N. Valid
ambiguous input letters normalize deterministically to N; no random bases are
generated. Anchor rules are source-compatible on ACGT; this N behavior is an
explicit extension, not equivalence to RaMA's random-N preprocessing.

See [source mapping and validation](centromere-rama-integration.md) for the
original RaMA backend dependencies and safety patches.

## Full genomes with known centromere intervals

Use two short, purpose-specific BED options:

```bash
ramag align --reference chm13.fa --query chm1.fa \
  --ref-cen chm13.centromeres.bed --query-cen chm1.centromeres.bed \
  --output comparison.paf --work-dir work --threads 16 \
  --print-effective-config
```

Long aliases are `--reference-centromeres` and `--query-centromeres`. The mode is
`genome-with-centromeres`: the requested scope remains the whole genomes, with
the marked regions routed to a specialised algorithm once that backend is
available. It is **not** a request to extract regions and ignore the rest of the
genomes. Reconciliation of ordinary and specialised alignments remains backend
work; this interface does not simply concatenate overlapping output sets.

Both BED paths are required. Their region counts may differ. BED row order,
fourth-column names and chromosome names never specify cross-genome pairing.
This permits different chromosome numbers without inventing correspondences.
The later backend must establish and report supported correspondences.

## BED contract

Use an uncompressed BED with at least three whitespace-separated columns:

```text
# contig    start    end      optional_name
chr1        100000   120000   cen_chr1
chr2        200000   230000   cen_chr2
```

- Coordinates are zero-based, half-open: `0 <= start < end`, represented with
  64-bit integers. The example interval `[100000,120000)` contains 20,000 bases.
- Contig IDs must match the first token of the corresponding FASTA header.
  Bounds and contig existence are checked against sequence records by the
  reusable validation helper. `--print-effective-config` checks BED syntax
  only: it does not read genomes or verify sequence bounds.
- BED3 and a fourth-column name are supported. Later columns are ignored; a
  BED strand column does not reverse-complement input or force an alignment
  orientation. The algorithm must search orientations rather than infer an
  inversion from a BED label.
- Blank lines, comment lines, LF/CRLF, a first-line UTF-8 BOM, and standard UCSC
  `track`/`browser` declarations are accepted. Unsorted intervals retain their
  input order.
- Empty annotations, negative/invalid/overflowing coordinates, empty intervals,
  duplicates and overlapping intervals on the same contig are errors. Adjacent
  intervals are allowed. Errors identify the file and source line when applicable.
- The software does not infer centromere features from GFF/GTF or repeat-class
  names. Convert your selected annotation intervals to BED before using them.

`-c` cannot be combined with either BED option. Repeating a flag or its alias is
an error. These options currently apply only to `align`, not `batch` or `index`.
All ordinary input/output, workspace, format and thread arguments keep their
existing contracts. With no centromere options, alignment behavior is unchanged.

The mode and backend readiness are visible in `--print-effective-config`.
Accepting a configuration does not imply that a backend is available, that a
biological correspondence is known, or that an inversion has been called.

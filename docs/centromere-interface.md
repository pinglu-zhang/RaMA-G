# Centromere input interface

Centromere input stays within `ramag align`. RaMA-G does not locate centromeres
automatically and does not require a species tree, a region-pairing table, or
matching chromosome names.

**Current stage:** Linux builds default to anchor-only export. `align -c` stops
after original RaMA anchor recursion and before any WFA or gap alignment.
The complete backend is retained but is not invoked on this route.
Inversions, automatic orientation,
region pairing and BED mixed alignment are not implemented. A disabled backend
or BED-mode run fails with dependency error **5**, without genome fallback.
The experimental [direction-region API](centromere-directions.md) is a separate
library capability and does not change these default CLI semantics.

## One centromere sequence on each side

Use `-c` or its identical long alias `--centromere`. Each input must contain
exactly one nonempty record. Plain/gzip FASTA and the existing normalization
rules apply. Multiple records are rejected, never concatenated or implicitly paired.

```bash
ramag align -c --reference reference_centromere.fa --query query_centromere.fa \
  --output-prefix comparison --work-dir work --threads 4
```

The reported mode is `centromere-only`. This flag is a declaration by the user,
not an automatic centromere detector. Both sequences are used in their supplied
orientation; reverse-complement the query externally when required.
Add `--print-effective-config` to inspect the source identities, scoring,
thread budget and N policy without running the alignment.
For `-c` only, omitting `--work-dir` creates per-run logs under
`.ramag-work` beside the output prefix. Ordinary genome commands retain their
existing explicit workspace requirement.

The current default exports `comparison.initial-anchors.tsv` and
`comparison.final-anchors.tsv`. Both contain stable anchor IDs, sequence IDs
and lengths, reference/query begin/end, strand, match length and source weight.
Coordinates are zero-based half-open. Strand is `+` in the supplied orientation;
there is no automatic reverse search. Empty anchor sets retain a header.
Use final anchors for plotting and initial anchors to compare recursion effects.
Do not connect anchors across unaligned gaps as if those gaps were aligned.

`--output` and `--formats` are rejected in this stage. No PAF, SAM, delta,
MAF or chain is produced. Effective configuration reports
`centromere_stage=anchors-only` and `wfa_execution=disabled`.

### Retained complete backend (not the current CLI default)

Only a build explicitly configured with internal
`RAMAG_CENTROMERE_ANCHORS_ONLY=OFF` enables the complete CLI route. Then existing
output options select PAF, SAM, delta, MAF or chain:

```bash
ramag align -c --reference ref.fa.gz --query query.fa.gz --threads 4 \
  --output-prefix result --formats paf,sam,delta,maf,chain
```

The complete internal CIGAR consumes both sequences, including terminal gaps.
PAF/SAM/delta/MAF preserve that consumption. UCSC chain has no terminal-gap
fields: the existing writer projects the paired blocks and reports their
actual coordinates, rather than inventing aligned bases. Paired-coordinate
sets remain identical, but chain spans can exclude terminal unpaired bases.
Chain cannot encode exact/mismatch distinctions or NM. Use PAF/SAM for those.

Global genome options (`--seed-mode`, `--selection-mode`, `--min-match`,
`--smem-min-occurrences`, `--min-cluster`, `--max-gap`, `--diag-diff`,
`--diag-factor`, `--break-length`, `--max-dp-cells`, `--save`,
`--reference-index`) are rejected when explicitly combined with `-c`.
The source's rare-level limit is 100, not a predecessor search window.

N is never an exact anchor or an `=` column, even opposite another N. Valid
ambiguous input letters normalize deterministically to N; no random bases are
generated. Anchor rules are source-compatible on ACGT; this N behavior is an
explicit extension, not equivalence to RaMA's random-N preprocessing.

See [source mapping and validation](centromere-rama-integration.md) for exact
dependencies, safety patches and limitations. The retained complete backend is
RaMA-anchor-constrained, not an unconstrained globally optimal alignment.

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

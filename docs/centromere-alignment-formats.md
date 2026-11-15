# Centromere alignment exports: PAF, SAM and MAF

The fixed-direction anchored-path stitcher exports the same completed segments
in three standard formats. It consumes existing `CentromereStitchResult` data:
export never searches for anchors, aligns another gap, changes copy acceptance,
or uses simulation truth. The authoritative path and gap TSVs remain available.
The format contract is `path-segments-v1`.

## CLI

A normal Linux build exports SAM, PAF and MAF from the same command.

```sh
ramag align -c --reference reference.fa --query query.fa \
  --output-prefix result --threads 1
```

This produces `result.paf`, `result.sam`, `result.maf` and the direction,
anchor, copy and alignment TSVs. `--formats paf` selects PAF among the standard
formats and keeps the diagnostic TSVs. Repeated explicit outputs also work.

```sh
ramag align -c --reference reference.fa --query query.fa \
  --output result.paf --output result.sam --output result.maf --threads 1
```

Explicit outputs share a directory and base prefix. Use either explicit
outputs or `--output-prefix` with optional `--formats`. The centromere path
format supports SAM, PAF and MAF, including accepted copy paths.

All selected outputs participate in the existing no-overwrite publication
transaction. A format failure rolls back this run's published files and removes
its temporary files. Existing files and input files are preserved. Diagnostic
checkpoints in the run directory remain available for failure diagnosis.

## Shared alignment meaning

One segment becomes one PAF record, one mapped SAM record, and one MAF block.
No record bridges a failed gap, a direction change, or a path boundary. Only
accepted copy branches are exported; ambiguity groups remain in their own TSV.
Different segments can share reference or query coordinates. Consumers must
preserve physical coordinates and path identity instead of deduplicating by a
single coordinate or forcing all segments into a one-to-one chain.

Input geometry, original sequence IDs and lengths, normalized ACGTN characters,
CIGAR consumption, edit counts and exact-match characters are validated before
writing a format. CIGAR uses `=`, `X`, `I`, `D`, traversing reference increasingly;
`I` consumes query and `D` consumes reference. N/N is `X`, not an exact match.

PAF and SAM use MAPQ 255 (unavailable); the exporter does not estimate mapping
confidence. Exporting a segment does not increase its scientific confidence.
`complete` refers to eligible internal gaps, not whole-input global coverage.
Partial results retain successful segments and advertise `partial`; aborted
diagnostic snapshots are rejected as formal alignment exports. The gap/summary
TSVs record unresolved regions and the scope of completion.

## PAF

- The 12 standard columns contain full input lengths and original 0-based,
  half-open intervals. Negative query coordinates are still on the original
  forward FASTA, with strand `-`.
- `cg:Z` contains the extended segment CIGAR, without clipping.
- Matching bases are the sum of `=`; alignment block length is the sum of all
  `=XID` runs. NM is the sum of `XID`.
- `tp` is omitted: RaMA-G primary/copy roles are expressed with `pr`, rather than
  being interpreted as minimap2 primary/secondary selection.
- An empty alignment result gives an empty PAF. See the summary TSV for run state.

## SAM

SAM v1.6 uses the original query FASTA name as QNAME, the original reference name
as RNAME, and POS=`reference_begin+1`. The header records the reference length.
SEQ contains only the aligned query fragment (reverse complemented for flag16),
QUAL is `*`, and unaligned flanks use hard clipping. For a negative interval
`[qb,qe)`, left clip=`query_length-qe` and right clip=`qb`. Both the CIGAR and SEQ
are in reference traversal order. Long individual operations, including clips,
are split to fit the BAM 28-bit operation length without changing consumption.

There is exactly one representative mapped record per query. Selection prefers
algorithm role `primary`, then greatest paired-base count, then smallest segment
ID. The remaining records are visited with primary roles before copy roles and
ascending segment ID. A query interval disjoint from the already selected
representative/supplementary group gets flag2048; an overlapping interval gets
flag256. Negative records additionally receive flag16. This stable format-level
classification neither removes records nor changes the algorithm's path roles.

**A secondary record can therefore contain an accepted additional copy relation.**
Do not filter out flags256/2048 when evaluating all RaMA-G relationships. Use
`pr`, `pi`, `bi` and `si` to distinguish roles and paths. SAM flags alone cannot
express the full accepted-copy semantics. No `SA` tags are emitted; metadata
links the records without generating a quadratic cross-reference list. If no
segment exists, SAM contains one unmapped query record (flag4) and the input SEQ.

The exporter rejects names outside the SAM v1.6 QNAME/RNAME grammar and reference
lengths outside the SAM header range; it does not silently rename or truncate.

## MAF

MAF v1 uses independent two-row blocks with the reference row first and always
on `+`. Each block contains its score and actual gapped strings. Source names
are `ref.<original_reference_name>` and `query.<original_query_name>`, so equal
FASTA names do not collide. Remove only the first role prefix when recovering
an original name.

For the negative query row, start=`query_length-query_end`, size is the number
of non-gap bases (`query_end-query_begin`), and the text is the aligned reverse
complement. Source size is the full input length. Blocks are not flattened into
a whole-input multiple alignment; duplicated relationships remain separate
blocks even when they reuse coordinates. Each block's `# ramag` comment stores
path metadata. Readers that discard comments lose that metadata; keep the
segment TSV alongside such transformed MAFs. Empty results contain headers only.

## Provenance fields

PAF/SAM custom tags use lower-case names. MAF writes the corresponding `key=value`
fields in a comment immediately before each block.

| Field | Meaning |
| --- | --- |
| `si` | Segment ID; joins `alignment-segments.tsv` |
| `pi` | Path ID |
| `pr` | Algorithm role: `primary` or `copy` |
| `ri`, `bi` | Original region and accepted-copy branch IDs |
| `sm` | Stitch scoring model |
| `st` | Overall stitch state: `complete` or `partial` |
| `sc` | Exact signed 64-bit segment score |
| `ai`, `gi` | Associated anchor and gap ID lists, when nonempty |

IDs and `sc` use SAM-compatible `Z` strings to avoid integer truncation. Standard
`AS:i` and `NM:i` are included within signed 32-bit range. Larger AS remains in
`sc:Z`; larger NM is represented by `ne:Z`. No clipping or saturation is applied.

## Library interface and checks

Include `ramag/centromere/stitch.hpp`, link `RaMAG::centromere`, and call:

```cpp
WriteCentromereStitchPaf(output, reference, query, stitched);
WriteCentromereStitchSam(output, reference, query, stitched);
WriteCentromereStitchMaf(output, reference, query, stitched);
```

These const-input writers are usable without rerunning the alignment when the
in-memory stitch result is retained. They validate all segments before writing;
callers remain responsible for their own stream/file publication transaction.
The public writers also build without the centromere computation backend.

Tests include independent hand-computed mixed CIGARs, negative coordinates,
copy overlaps, hard clipping, N/N, partial gaps, empty results, malformed inputs,
64-bit metadata, and BAM run encoding. The synthetic CLI test independently
reconstructs every exported coordinate pair and checks it against the TSV. When
SAMtools is installed, it additionally verifies SAM-to-BAM-to-SAM preservation.
These tests establish export correctness, not alignment accuracy improvements.

## Format references

- [PAF specification](https://github.com/lh3/miniasm/blob/master/PAF.md) and
  [minimap2 CIGAR tag](https://github.com/lh3/minimap2/blob/master/minimap2.1).
- [SAM v1.6](https://samtools.github.io/hts-specs/SAMv1.pdf) and
  [optional tags](https://samtools.github.io/hts-specs/SAMtags.pdf).
- [UCSC MAF specification](https://genome.ucsc.edu/FAQ/FAQformat.html#format5).

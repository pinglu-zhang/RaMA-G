# Formats and coordinates

## Select output files

| Command | Default with an output prefix | Supported standard formats |
|---|---|---|
| `align` | SAM, PAF, delta | SAM, PAF, MAF, delta, chain |
| `align -c` | SAM, PAF, MAF plus TSVs | SAM, PAF, MAF |
| `batch` | PAF per query | SAM, PAF, MAF, delta, chain |

Repeated `--output result.paf --output result.maf` selects formats by suffix.
Explicit outputs share one directory and base prefix. Alternatively, use
`--output-prefix result --formats paf,maf`. Do not mix these two styles.
Compressed outputs and BAM/CRAM are not supported. Each format is validated
before publication; existing outputs are never overwritten.

## Coordinate conventions

| Format | Coordinates and orientation |
|---|---|
| PAF | Original forward-input, zero-based half-open intervals; query strand is explicit |
| SAM | One-based reference POS; CIGAR and sequence follow reference traversal; flag16 marks reverse query |
| MAF | Zero-based start on each row's stated strand; reverse-query start is `query_length-query_end` |
| delta | One-based inclusive endpoints; reverse query has reversed endpoints |
| chain | Positive reference target; query coordinates interpreted on its stated strand |
| Centromere TSV | Original forward-input, zero-based half-open intervals with an explicit strand |

CIGAR traverses increasing reference coordinates. `=` is an exact canonical
match, `X` a substitution, `I` consumes query only and `D` consumes reference
only. `N/N` is not an exact match. Insertions/deletions consume bases without
creating a paired column. Identity, edit count and alignment score describe the
emitted alignment; they do not establish a unique historical correspondence.

## PAF

PAF has the standard 12 columns plus tags. `cg:Z` carries the actual extended
CIGAR without clipping. Negative-strand intervals remain on the original forward
query FASTA. Read the strand and CIGAR rather than pairing an entire interval
rectangle. MAPQ is currently 255 (unavailable), not a confidence estimate.
A valid empty result produces an empty PAF; use exit status and logs to distinguish
it from failure.

Centromere PAF does not use `tp` to classify primary/copy paths; use `pr` and the
other path tags below. The matching-base count is the sum of `=` runs and the
block length is the sum of all `=XID` runs.

## SAM

SAM records use original query/reference identifiers and full reference lengths.
NM/MD and score tags describe the alignment, not calibrated mapping confidence.

For centromere output, SEQ contains the aligned query fragment, reverse
complemented on flag16; QUAL is `*`. Unaligned flanks use hard clipping. For a
negative query interval `[qb,qe)`, left clip is `query_length-qe` and right clip
is `qb`. Individual long operations are split to fit the BAM operation-length
representation without changing consumption.

Centromere output chooses one representative mapped record per query, preferring
the primary path role, then paired-base count and segment ID. Further disjoint
query segments receive flag2048 and overlapping segments receive flag256;
negative segments additionally receive flag16. **An accepted copy relationship
can therefore be a secondary record.** Do not discard every flag256/2048 record
when interpreting all accepted relationships. Preserve path metadata as well as
SAM flags. No quadratic list of `SA` cross-references is emitted. An empty
centromere result contains an unmapped query record.

Names and lengths must fit the SAM specification; invalid names are not silently
renamed. If converting to BAM with another tool, retain a copy of the original
PAF/TSV metadata where that conversion cannot preserve every field.

## Pairwise MAF

Each MAF block has two rows: positive reference first, then the oriented query.
Minus-query start is `query_length-query_end`; size counts non-gap bases.
For centromeres, source names are `ref.<original_name>` and
`query.<original_name>` to distinguish equal FASTA names. Remove only that first
role prefix to recover the original name.

Centromere blocks retain duplicated relationships even when coordinates overlap.
A preceding `# ramag` comment carries metadata. Tools that discard comments lose
that metadata, so retain the segment TSV alongside transformed MAF files.
An empty result has headers only. Blocks are pairwise segments, not a whole-input
multiple alignment.

## Delta and UCSC chain

These formats are available for genome alignment. Delta uses indel displacement
encoding; input FASTA paths must not contain whitespace. Chain normalizes
representable blocks and omits unrepresentable terminal gaps. Both formats must
be interpreted with their own endpoint/strand conventions, not PAF coordinates.
MAF and chain preserve the emitted genome record set rather than applying a new
SAM-specific filtering step. Batch candidate output can include alternative
placements; ordinary align outputs the linked bilateral main path.

## Centromere paths and unresolved regions

One completed segment becomes one PAF record, mapped SAM record or MAF block.
Records do not bridge failed gaps, direction changes or path boundaries. Accepted
copy paths can share physical coordinates. Ambiguous copy candidates are kept in
their diagnostic tables rather than treated as accepted alignments.

Keep these files with standard-format output:

| File | Use |
|---|---|
| `PREFIX.alignment-segments.tsv` | Segment/path roles and IDs, physical coordinates, strand, actual CIGAR, score and associated anchor/gap IDs |
| `PREFIX.alignment-gaps.tsv` | Internal gap status/reason, coordinates, engine, CIGAR, score, resource observations and time |
| `PREFIX.alignment-summary.tsv` | Overall state, completion scope and aggregate counts |
| Other direction/anchor/copy TSVs | Inspect supporting evidence, accepted copies, ambiguity and unanchored regions |

Use the table header for exact column names. `complete` refers to eligible
internal-gap work, not whole-input alignment coverage. Partial results retain
completed segments; aborted diagnostic snapshots are not successful final
alignments. Unresolved regions and unaligned tails must not be counted as paired
bases merely because they lie between two reported coordinates.

### Path metadata

PAF/SAM tags use lower-case names; MAF comments use corresponding `key=value`
fields. Tags join records to the diagnostic tables:

| Tag | Meaning |
|---|---|
| `si` | Segment ID |
| `pi` | Path ID |
| `pr` | `primary` or `copy` algorithm role |
| `ri`, `bi` | Original region ID and accepted-copy branch ID |
| `sm` | Stitch scoring model |
| `st` | Overall stitch state: `complete` or `partial` |
| `sc` | Exact signed 64-bit segment score |
| `ai`, `gi` | Anchor and gap ID lists, when present |

IDs and `sc` are string tags to avoid truncation. AS/NM integer tags are included
when representable in signed 32-bit form; larger edit counts use `ne:Z` and the
exact score remains in `sc:Z`. No clipping or saturation of those values occurs.

See [commands](user-guide.md) for inputs and options, or
[troubleshooting](troubleshooting.md) for empty and partial results.

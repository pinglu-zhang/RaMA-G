# Quick start

Run these commands from the repository root after building the `release` preset.
The FASTA files are artificial examples, not a biological dataset. Choose a new
output directory when repeating the commands.

## Genome alignment

```bash
./build/release/ramag align \
  --reference examples/genome/reference.fa \
  --query examples/genome/query.fa \
  --output results/genome.paf \
  --output results/genome.maf \
  --work-dir results/genome-work \
  --threads 2
```

The output is the linked bilateral main path. PAF `cg` describes the actual
base alignment. Inspect the run log below the work directory and check exit
status zero. The complete [user guide](user-guide.md) explains formats and logs.

## Save and load an index

```bash
./build/release/ramag index \
  --reference examples/genome/reference.fa \
  --output results/reference.sufidx \
  --work-dir results/index-work
./build/release/ramag align \
  --reference examples/genome/reference.fa \
  --reference-index results/reference.sufidx \
  --query examples/genome/reverse.fa \
  --output results/reverse.paf \
  --work-dir results/reverse-work
```

The reference FASTA remains required and is validated against the stored index.
The reverse example contains the reverse complement of the artificial query.

## Centromere alignment

```bash
./build/release/ramag align \
  -c \
  --reference examples/centromere/reference.fa \
  --query examples/centromere/query.fa \
  --output-prefix results/centromere
```

Each input has one record. Default outputs include `.sam`, `.paf`, `.maf`,
`.alignment-segments.tsv` and `.alignment-gaps.tsv`, together with direction,
anchor, copy and review diagnostics. `--formats paf` restricts standard alignment
formats while retaining the diagnostic tables. A supported internal gap may be
complete even when input tails or an ambiguous region remain unaligned.

## Batch queries

```bash
./build/release/ramag batch \
  --reference examples/genome/reference.fa \
  --seqfile examples/batch/queries.txt \
  --output-dir results/batch \
  --threads 2
```

The two query paths in the seqfile are relative to its own directory. Batch runs
queries serially using one resident index and writes `batch.tsv`. It retains
candidate-output semantics, whereas ordinary `align` emits a linked one-to-one
main path. Do not treat their records as the same output policy.

## Inspect effective configuration

Append `--print-effective-config` to a valid command to print `key=value` text
and exit before reading FASTA contents or building an index. Use `--help` for the
current command-line options. The old selection/linking/recovery switches are
not accepted; see [migration](migration.md).

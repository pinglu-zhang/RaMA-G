# Commands and parameters

Use the [README](../README.md) for runnable artificial examples and
[installation](installation.md) for build requirements.

## Choose a command

| Command | Use |
|---|---|
| `align` | Pairwise whole-genome alignment with linked bilateral main-path output |
| `align -c` | One centromere FASTA record per side, with direction and copy paths |
| `index` | Save a reusable whole-genome reference index |
| `batch` | Run queries serially using a shared reference index; retain candidate output |

RaMA-G aligns assemblies, not short reads. Genome and centromere commands have
different workflows and output defaults. In particular, `batch` is not a repeated
invocation of ordinary `align` with the same selection policy.

## FASTA inputs

Both sides accept plain and gzip FASTA, including BGZF and concatenated gzip
members. Compression is detected from bytes. Damaged gzip or a gzip-looking
suffix without gzip data is rejected. No external decompressor or FAI/GZI is
required; existing index companions are ignored.

Record IDs are the first whitespace-delimited header token and must be unique
within each input. Empty records, FASTQ, sequence before the first header,
nonletter sequence characters and unsupported compressed formats are rejected.
LF and CRLF are supported. Letters are uppercased; non-ACGT letters become `N`.
`N` and record boundaries break exact seed matching. A single record longer than
`INT_MAX` bases is rejected; collections may total more than that. Inputs are
loaded into memory.

Delta output requires input paths without whitespace. Other validation and
output-name restrictions are described in [formats](formats.md).

## Ordinary genome alignment

```bash
ramag align --reference reference.fa --query query.fa \
  --output result.paf --output result.maf --work-dir work --threads 8
```

Ordinary alignment selects a bilateral one-to-one main path and links compatible
records, including supported bounded gaps. Alternative repeat placements can be
omitted. The former `--selection-mode`, `--link-gap-mode` and
`--link-recovery-mode` switches are no longer accepted; remove them from older
commands. This output policy differs from older all-candidate output.

Use either repeated `--output` or an output prefix:

```bash
ramag align --reference reference.fa --query query.fa \
  --output-prefix result --formats sam,paf,maf --work-dir work
```

Explicit outputs must have the same directory and base prefix, with each
supported suffix appearing once. Do not mix `--output` with `--output-prefix`
or `--formats`. Defaults and file meanings are in [formats](formats.md).

## Centromere pairs

```bash
ramag align -c --reference reference.fa --query query.fa \
  --output-prefix centromere --threads 1
```

Each FASTA must contain exactly one record. The command predicts directions,
finds recursive and accepted-copy anchors, reviews paths and aligns supported
internal gaps. It writes SAM, PAF, MAF and diagnostic TSVs by default. Restrict
standard formats with `--formats paf`, or request explicit `.paf`, `.sam` and
`.maf` files with repeated `--output`.

Unknown directions, ambiguous copy choices, tails and unresolved gaps can remain
unaligned. Read the summary/gap TSVs alongside the standard output. A completed
run does not mean that every input base has a homologous partner or is aligned.
See [formats](formats.md) for copy/path metadata.

When omitted, the centromere work directory is `.ramag-work` beside the output.
Genome seed/tuning and persistent-index options cannot be combined with `-c`.
A genome-only build cannot run this command. BED-based mixed genome/centromere
alignment is not implemented: `--ref-cen` and `--query-cen` can be inspected with
`--print-effective-config`, but do not provide a runnable mixed-alignment mode.

## Persistent reference indexes

Build an index and reuse it:

```bash
ramag index --reference reference.fa.gz --output reference.sufidx \
  --work-dir index-work --threads 8
ramag align --reference reference.fa.gz --reference-index reference.sufidx \
  --query query.fa --output result.paf --work-dir align-work --threads 8
```

Alternatively, add `--save reference.sufidx` to a genome `align` or `batch` command
to save its newly built index. `--save` and `--reference-index` are mutually
exclusive. Index files and outputs are never overwritten.

The reference FASTA is still required and must match the indexed normalized
reference. A mismatch or corrupt index fails rather than silently rebuilding.
The index is one `.sufidx` file. There is no automatic index discovery and no
required companion manifest. Compatible older indexes remain loadable; unrelated
old companion files are left untouched. A successfully saved index survives a
later alignment failure.

## Batch queries

```bash
ramag batch --reference reference.fa --reference-index reference.sufidx \
  --seqfile queries.txt --output-dir results --threads 8
```

A seqfile has one name and FASTA path per line:

```text
query1  data/query1.fa.gz
query2  data/query2.fa.gz
```

Paths are relative to the seqfile directory. After the first whitespace
separator, the remaining trimmed text is the path; internal spaces are allowed.
Blank lines and whole-line `#` comments are supported. Shell variables, globs
and quoting are not expanded. Names must be unique and use letters, digits,
`.`, `_` or `-`; empty names and reserved path/table names are rejected.
The [supplied seqfile](../examples/batch/queries.txt) is a runnable example.

Repeated `--query PATH` is an alternative to `--seqfile`; do not combine them.
Names then come from query filenames. If filenames collide, use a named seqfile.
`--output-dir` is required. The default work directory is `OUTPUT_DIR/.ramag-work`.
Each query writes `OUTPUT_DIR/NAME/alignment.paf` by default; `--formats` selects
other genome formats. `--output` and `--output-prefix` are not batch options.

One reference/index stays resident while queries run serially. Missing or invalid
query inputs are recorded and later queries continue. Initialization errors or
interruption stop the batch. Successfully completed outputs remain available.
There is no query concurrency, prefetch or batch resume. A batch can retain more
memory than standalone alignment because its reference index stays resident.

`batch.tsv` reports order, name, query path, status, exit code, elapsed time,
alignment record count, output directory, log path and message. Status is
`success`, `failed`, `interrupted` or `not-run`. A currently running query has
no terminal row yet. The table is flushed per finished query. Per-query elapsed
time excludes shared initialization. Tabs/newlines/backslashes are escaped.
Exit zero means all queries succeeded; otherwise the first failed query's code
is returned, or 130/143 for interruption. Existing outputs/tables are protected.

## Parameters

| Option | Default | Use |
|---|---:|---|
| `--threads N` | 1 | Worker budget, not a memory limit |
| `--progress auto\|on\|off` | `auto` | Terminal progress; auto requires a TTY |
| `--progress-interval N` | 10 | Seconds between progress reports; minimum 1 |
| `--print-effective-config` | off | Validate and print configuration without aligning |
| `--seed-mode MODE` | `fast` | Genome seed mode |
| `--min-match N` | 20 | Minimum exact genome seed length |
| `--smem-min-occurrences N` | 1 | Reference occurrence threshold, SMEM mode only |
| `--max-gap N` | 90 | Genome seed-chain gap constraint |
| `--diag-diff N` | 5 | Absolute genome diagonal slack |
| `--diag-factor X` | 0.12 | Relative genome diagonal slack |
| `--min-cluster N` | 65 | Minimum genome matching support |

Seed modes are `fast`, `mumreference`, `mum`, `smem`, and `maxmatch`.
`fast` is the default; `mumreference` requires reference uniqueness while
allowing query repetition; `mum` requires strict combined reference/query-record
uniqueness. `smem` and `maxmatch` can enumerate many repeat occurrences and use
more resources. An explicit `--smem-min-occurrences` outside SMEM is rejected.
Use `--help` for the current options.

`--break-length` and `--max-dp-cells` are compatibility-only arguments: only their
default values, 200 and 4000000 respectively, are accepted. They do not select a
solver or change centromere gap budgets. An explicit progress interval cannot
accompany `--progress off`. Effective configuration is `key=value` text, not JSON.

## Logs and exit status

Each run logs under `WORK/runs/RUN_ID/run.log`. Logs include the command,
configuration, input counts, stages, resource observations and final status.
Batch query logs are under `WORK/queries/NAME/runs/RUN_ID/run.log`.
Progress-off does not suppress final errors or the run log. Linux RSS is a
process observation, not the sum of stage memory capacities.

| Exit code | Meaning |
|---|---|
| 0 | Requested outputs validated and published |
| 2 | CLI/configuration error |
| 3 | FASTA/input error |
| 4 | Unsupported input or mode |
| 5 | Dependency/index error |
| 6 | Compute/resource error |
| 7 | Output, validation or logging error |
| 8 | Unexpected internal exception |
| 130 / 143 | SIGINT / SIGTERM interruption |

Check process exit status and logs rather than waiting for a `.complete` file
or `.manifest.json`; these completion artifacts are not produced. For help with
failed or partial runs, see [troubleshooting](troubleshooting.md).

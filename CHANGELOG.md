# Changelog

## 0.1.1 — Unreleased

### Index loading

- Pin Sufkit to `028075e6f2d622fcbdcf76b153bbde0f069ca64f`.
- Accumulate section CRC during index reads instead of scanning the payload
  separately before decoding. All sections still require successful validation.
- Parallelize full SA/ISA consistency validation and avoid its redundant
  permutation bitmap. Indexes without ISA retain permutation validation.
- Pass the requested loading thread budget to Fast prefix-directory construction.
- Validate the normalized reference through read-only sequence views, avoiding
  another owning copy of the reference collection.
- Expose loading stages, logical read bytes and CRC CPU time in logs/statistics.

### Compatibility

The suffix-index format, complete validation requirements, alignment algorithms,
scoring, CLI and defaults are unchanged by this loading update. Compatible raw-
and byte-coded-LCP indexes remain reusable. Sufkit retains `Load(path)` alongside
its options overload. CRC, structure, capability and reference matching checks
remain mandatory in the RaMA-G loading path.

The product version becomes 0.1.1, including CLI identification, installed CMake
package metadata and SAM `@PG VN`; that header change is not an alignment change.

### Existing input and automation contract

The source preceding this patch already uses kseq/zlib for plain and gzip FASTA,
spdlog logging, and single-file persistent indexes. Commands do not generate
SHA-256 manifests or `.complete` markers. Automation migrating from earlier
snapshots should check exit status, parse requested outputs and inspect
`WORK/runs/<run-id>/run.log`. A saved index can survive a later alignment failure.
See the [user guide](docs/user-guide.md#exit-codes-and-migration).

### Previously measured loading performance

Before this version/documentation preparation, three serial paired server runs
loaded the same 32,986,750,296-byte GRCh38 byte-coded index with 16 threads, using
the same warm-cache procedure and alternating baseline/candidate order.

| Load-only metric | Baseline median | Optimized median |
|---|---:|---:|
| RaMA-G adapter load, including reference matching | 522.593 s | 163.754 s |
| Loading process-tree peak RSS | 39.919 GiB | 33.940 GiB |

These are loading measurements on that server, not complete Human–Chimp
alignment results. The server had shared load; optimized adapter times ranged
from 162.872 to 163.771 s. No new coverage or F1 claim follows from this result.

### Verification boundary

The accepted implementation preceding this version change passed Sufkit
Release/Werror and ASan/UBSan (24/24 each), RaMA-G equivalents (19/19 each),
107 runtime-contract scenarios and a 40-run alignment regression matrix.
These are retained implementation results. This pending release preparation
performs static checks only; it does not rerun builds, tests or experiments.

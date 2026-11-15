# Public synthetic tests

These tests exercise the installed product without real genomes or benchmark
artifacts. `cli` checks input, gzip, indexes, publication failures and removed
options. `formats` independently expands genome SAM/PAF/MAF/delta/chain into
original-coordinate pairs and checks a known artificial ancestry. `centromere`
independently checks path-aware formats, scores, unresolved scope, cancellation
and rollback, or the backend-OFF error. `install` compiles retained public APIs
using an isolated install prefix.

The scripts use only Python's standard library. Samtools adds an optional SAM
round trip. CTest registers this small suite directly; CI treats an absent suite as an error.

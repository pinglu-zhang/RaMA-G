# Third-party notices

RaMA-G itself is licensed under the MIT License. The baseline implementation in
this repository does not copy source code from RaMAx, MUMmer4, minibwa, sufkit,
or SeqPro.

The production dependency contract pins the following separately maintained
libraries. Their own licenses and notices continue to apply when fetched or
linked:

- **sufkit**, commit `bdb67c6de5daddd8a005640de73d96549d2575f4`.
- **SeqPro**, commit `6781cadcf81a0da53d7573444594c1484947017c`.
- **zlib**, discovered through CMake's `ZLIB::ZLIB` target and linked for
  in-process streaming gzip/BGZF-compatible FASTA decompression. The zlib
  license and the system/package distributor's notices apply.

MUMmer4 and RaMAx are read-only algorithm/compatibility references and are not
linked into RaMA-G. minibwa commit
`f0e117436c28addc359b67123d2353f0d4a1f9e8` is an evidence and potential
MIT-compatible-code source only. Its optional GPL-2.0 `bwtgen` component is
explicitly excluded from RaMA-G.

KSW2 is planned as the optimized gap-extension implementation. Until that
adapter and its source notice are added, RaMA-G uses its own bounded scalar
dynamic-programming baseline. A future KSW2 import must retain the upstream
copyright and MIT notice in this file and in the vendored source directory.

# Third-party notices

RaMA-G itself is licensed under the MIT License. Vendored source attribution
and separately maintained dependency licenses are described below.

The production dependency contract pins the following separately maintained
libraries. Their own licenses and notices continue to apply when fetched or
linked:

- **sufkit**, commit `f8c4c386ee883e45ad0f973efc4c8e1148b0068a`.
- **SeqPro**, commit `6781cadcf81a0da53d7573444594c1484947017c`.
- **zlib**, discovered through CMake's `ZLIB::ZLIB` target and linked for
  in-process streaming gzip/BGZF-compatible FASTA decompression. The zlib
  license and the system/package distributor's notices apply.

MUMmer4 is a read-only compatibility reference, not a linked dependency. minibwa commit
`f0e117436c28addc359b67123d2353f0d4a1f9e8` is an evidence and potential
MIT-compatible-code source only. Its optional GPL-2.0 `bwtgen` component is
explicitly excluded from RaMA-G.

The opt-in graph-free pairwise core incorporates algorithm code and a
KSW2 snapshot from RaMAx `7d08359e0df7f7e6ffcfe67217c3399761cb2129`.
See `third_party/attribution/pairwise-source.json` for source identities and adaptations and
`third_party/attribution/LICENSE.pairwise` for the original MIT copyright and license. The default legacy scalar
build does not compile this optional snapshot. KSW2's upstream MIT notice is
retained separately in `third_party/ksw2/LICENSE`.

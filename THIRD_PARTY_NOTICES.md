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

Optional compile-time-only extension experiments use the following separately
fetched MIT-licensed implementations. Neither is linked by the default scalar
build:

- **KSW2**, commit `289609bd9e5381a13b16239d0a7703f1ff03f9ca`,
  copyright Dana-Farber Cancer Institute and Broad Institute, Inc. The upstream
  MIT license applies to the selected `ksw2_extz2_sse.c` and `ksw2_gg.c`
  implementations.
- **Block Aligner**, commit
  `4fcf630cf775de5b578fe63971f210e1dc958791`, copyright Daniel Liu. The
  upstream MIT license applies to its Rust implementation and C ABI.

RaMA-G keeps its bounded scalar dynamic-programming implementation as the
default. Selecting either experimental dependency requires an explicit
`RAMAG_INTERNAL_EXTENSION_BACKEND` CMake value and does not add a public CLI
backend switch.

Their common MIT license text is reproduced below with both upstream copyright
notices:

> Copyright (c) 2018- Dana-Farber Cancer Institute<br>
> Copyright (c) 2017-2018 Broad Institute, Inc.<br>
> Copyright (c) 2021 Daniel Liu
>
> Permission is hereby granted, free of charge, to any person obtaining a copy
> of this software and associated documentation files (the "Software"), to deal
> in the Software without restriction, including without limitation the rights
> to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
> copies of the Software, and to permit persons to whom the Software is
> furnished to do so, subject to the following conditions:
>
> The above copyright notice and this permission notice shall be included in all
> copies or substantial portions of the Software.
>
> THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
> IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
> FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
> AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
> LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
> OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
> SOFTWARE.

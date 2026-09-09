# Current whole-genome benchmark runner

This portable runner follows RaMA-G's exit-code, output and run-log contract.
It does not require or read program manifests, hashes or completion markers.
It replaces the old fixed-server harness; old `--package`, `--stage` and `--resume`
arguments are rejected. Historical experiments remain interpretable with the
script revision that created them, not with this new driver.

```bash
python3 benchmarks/human_chimp/run_human_chimp_benchmark.py \
  --binary build-release/ramag --reference reference.fa --query query.fa \
  --reference-index reference.sufidx --result-root benchmark-new \
  --threads 16 --repetitions 1
```

The result root must be new. Runs are serial with `mumreference`, `one-to-one`,
PAF-only and periodic progress off. A saved index is optional; if supplied, each
log must report it as loaded. No implicit index construction is substituted.
The wrapper records `/usr/bin/time -v`, sampled process-tree resources and exit
status in benchmark-owned reports. These are distinct from program outputs.

The summarizer validates nonempty PAF, coordinates, input catalog and CIGAR
consumption; coverage is the interval union divided by complete input lengths.
It compares repetitions directly in bounded chunks. It reports time, sampled
peak memory and coverage, with no F1 claim on real data. It includes input/index
loading in wall time and makes no assumption about cache warmth or server load.
To reread completed runs:

```bash
python3 benchmarks/human_chimp/summarize_comparison.py --result-root benchmark-new
```

The runner does not launch MUMmer4 or alter accepted historical MUMmer4 results.
Use the repository's explicit TSV-plan harness for separately authorized tool
comparisons, and retain the exact endpoint and timing definitions.
`canonicalize_maf_sources.py` remains the frozen simulation alias transformer;
it now validates unchanged payload by paired streaming comparison, without hashes.

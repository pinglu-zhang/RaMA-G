# Reproducible benchmark harness

This directory is infrastructure, not a performance claim. Do not populate the
owner-decision fields or use a convenient local dataset to declare that RaMA-G
has surpassed MUMmer4. Copy `benchmark-plan.example.tsv`, fill every `OWNER_TBD`
field with an owner-approved value, and preserve that accepted file with every
result bundle.

`run_benchmark.sh` records raw `/usr/bin/time -v` output, command lines, binary
versions, host information, completion-marker status, and artifact sizes. It
supports both RaMA-G delta-only headline runs and a separate all-format product
throughput run. MUMmer4 commands are supplied explicitly in the plan so the
harness never guesses a semantically unfair comparison.

Only runs with exit code zero, a non-empty requested artifact, a valid RaMA-G
manifest and `.complete` marker (where applicable), and successful external
format validation are eligible for aggregation. Cold-cache manipulation is not
performed automatically because it usually requires privileged, system-wide
state changes; cold and warm runs must instead be scheduled and recorded by the
experiment owner.

The fixed simulated Human–Chimp MAM protocol has a separate staged driver,
strict MAF source canonicalizer, process-tree resource recorder, and comparison
summarizer under [`human_chimp/`](human_chimp/README.md). Its default action is
preflight only; full alignments and million-sample scoring require explicit
stages.

# Centromere direction evidence

The complete `align -c` workflow predicts direction from input sequence evidence.
Its V9 planner combines initial forward/reverse anchors, local direction support,
reference-placement checks, neighboring context and the retained protection rules.
A direction label is not sufficient evidence for a unique reference copy.

Unknown or competing evidence can remain unresolved. Current gap, placement,
opposition, island and protection rules run in their defined order. Grouping their
implementation by function does not merge their acceptance criteria into one vote.
The planner uses input FASTA and its own search evidence; it does not consume
simulation truth or downstream benchmark scores.

Reverse searches use canonical reverse-complement coordinates internally. Output
query coordinates always refer to the original forward FASTA. All diagnostic
intervals are zero-based and half-open. Known directions constrain the recursive
anchor stage; later reviews update corresponding initial/final tables through the
coordinated replacement operations.

The program preserves direction evidence and planned-region TSVs as diagnostics.
Use the effective configuration and run log to identify the current workflow.
Coverage is not a homology accuracy score, and historical inversion events are not
an input to direction prediction. See [directed anchors](centromere-directed-anchors.md)
and [format contracts](centromere-alignment-formats.md).

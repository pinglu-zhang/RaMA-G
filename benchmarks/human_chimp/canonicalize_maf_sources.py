#!/usr/bin/env python3
"""Canonicalize the eight explicitly permitted Alignathon MAF source aliases.

This utility intentionally is *not* a general species-name rewriter.  It only
accepts the eight canonical Human/Chimp contig identifiers and the eight exact
``species.species.chr`` aliases listed below.  Every non-``s`` line and every
field after the source name is copied byte-for-byte.
"""

from __future__ import annotations

import argparse
import gzip
import json
import itertools
import os
import re
import sys
from collections import Counter
from contextlib import contextmanager
from pathlib import Path
from typing import Any, BinaryIO, Iterator


CANONICAL_SOURCES = tuple(
    f"{species}.chr{chromosome}"
    for species in ("simHuman", "simChimp")
    for chromosome in "ABCD"
)
CANONICAL_SET = frozenset(CANONICAL_SOURCES)
ALIASES = {
    f"{source.split('.', 1)[0]}.{source}": source for source in CANONICAL_SOURCES
}
S_LINE_RE = re.compile(rb"^([ \t]*s[ \t]+)([^ \t\r\n]+)(.*)$", re.DOTALL)


class CanonicalizationError(RuntimeError):
    """A validation failure that must stop benchmark evaluation."""


@contextmanager
def _open_input(path: Path) -> Iterator[BinaryIO]:
    if path.suffix == ".gz":
        with gzip.open(path, "rb") as handle:
            yield handle
    else:
        with path.open("rb") as handle:
            yield handle


def canonicalize(input_path: Path, output_path: Path) -> dict[str, object]:
    """Rewrite allowed aliases and atomically publish ``output_path``.

    A block that would contain the same canonical source twice is rejected as a
    normalization collision.  This catches, for example, one canonical row and
    one doubled-prefix row for the same contig in a two-row prediction block.
    """

    if not input_path.is_file():
        raise CanonicalizationError(f"input MAF does not exist: {input_path}")
    if output_path.exists():
        raise CanonicalizationError(f"refusing to overwrite output MAF: {output_path}")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = output_path.with_name(f".{output_path.name}.tmp.{os.getpid()}")
    if temporary.exists():
        raise CanonicalizationError(f"temporary output already exists: {temporary}")

    replacements: Counter[str] = Counter()
    rows_by_source: Counter[str] = Counter()
    seen_sources: set[str] = set()
    block_sources: set[str] = set()
    s_rows = 0
    line_count = 0
    block_count = 0
    in_block = False

    try:
        with _open_input(input_path) as source, temporary.open("xb") as destination:
            for line_number, line in enumerate(source, 1):
                line_count += 1
                match = S_LINE_RE.match(line)
                if match is None:
                    destination.write(line)
                    stripped = line.strip()
                    if stripped == b"" or stripped == b"a" or stripped.startswith(b"a "):
                        block_sources.clear()
                        in_block = False
                    continue

                if not in_block:
                    block_count += 1
                    in_block = True
                prefix, encoded_source, suffix = match.groups()
                try:
                    original = encoded_source.decode("ascii")
                except UnicodeDecodeError as exc:
                    raise CanonicalizationError(
                        f"non-ASCII MAF source on line {line_number}"
                    ) from exc

                if original in CANONICAL_SET:
                    canonical = original
                elif original in ALIASES:
                    canonical = ALIASES[original]
                    replacements[original] += 1
                else:
                    raise CanonicalizationError(
                        f"unknown or forbidden MAF source on line {line_number}: {original}"
                    )

                if canonical in block_sources:
                    raise CanonicalizationError(
                        "source collision after canonicalization in block "
                        f"{block_count}: {canonical}"
                    )
                block_sources.add(canonical)
                seen_sources.add(canonical)
                rows_by_source[canonical] += 1
                s_rows += 1

                # The suffix contains start/size/strand/srcSize/text and its
                # original whitespace/newline; direct paired re-reading verifies it.
                destination.write(prefix)
                destination.write(canonical.encode("ascii"))
                destination.write(suffix)

        with _open_input(input_path) as before, temporary.open("rb") as after:
            for a,b in itertools.zip_longest(before,after):
                if a is None or b is None:raise CanonicalizationError("line count changed")
                ma=S_LINE_RE.match(a);mb=S_LINE_RE.match(b)
                if ma is None:
                    if a!=b:raise CanonicalizationError("non-s line changed")
                elif mb is None or ma.group(1)!=mb.group(1) or ma.group(3)!=mb.group(3):
                    raise CanonicalizationError("s-line payload changed")
                else:
                    original=ma.group(2).decode("ascii")
                    expected=ALIASES.get(original,original).encode("ascii")
                    if mb.group(2)!=expected:raise CanonicalizationError("canonical source mismatch")
        os.link(temporary, output_path)
        temporary.unlink()
    except BaseException:
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass
        raise

    return {
        "schema": "ramag.maf-source-canonicalization.v1",
        "input": str(input_path.resolve()),
        "output": str(output_path.resolve()),
        "input_bytes": input_path.stat().st_size,
        "output_bytes": output_path.stat().st_size,
        "canonical_whitelist": list(CANONICAL_SOURCES),
        "alias_map": dict(sorted(ALIASES.items())),
        "replacement_counts": {
            alias: replacements.get(alias, 0) for alias in sorted(ALIASES)
        },
        "total_replacements": sum(replacements.values()),
        "line_count": line_count,
        "block_count": block_count,
        "s_rows": s_rows,
        "canonical_sources_seen": sorted(seen_sources),
        "rows_by_source": dict(sorted(rows_by_source.items())),
        "non_s_lines_unchanged": True,
        "s_line_payload_unchanged": True,
        "collision_count": 0,
        "status": "success",
    }


def _write_json_exclusive(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8", newline="\n") as handle:
        json.dump(value, handle, indent=2, sort_keys=True, ensure_ascii=False)
        handle.write("\n")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--stats", required=True, type=Path)
    args = parser.parse_args(argv)
    if args.stats.exists():
        raise CanonicalizationError(f"refusing to overwrite stats JSON: {args.stats}")
    result = canonicalize(args.input, args.output)
    _write_json_exclusive(args.stats, result)
    print(json.dumps(result, indent=2, sort_keys=True, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except CanonicalizationError as error:
        print(f"canonicalize_maf_sources: {error}", file=sys.stderr)
        raise SystemExit(2)

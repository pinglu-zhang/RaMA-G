#!/usr/bin/env python3
"""Stage the fixed preliminary Human--Chimp MAM comparison.

The safe default is ``--stage preflight``.  Full alignments and the expensive
three-seed quality evaluation only start when their stage is named explicitly.
An existing result root is accepted only together with ``--resume``; failed
attempt directories are preserved and the next resume uses a new attempt.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

from canonicalize_maf_sources import CANONICAL_SET, canonicalize
from run_with_metrics import (
    CLOCK_CONSISTENCY_TOLERANCE_FRACTION,
    assess_clock_consistency,
    run as run_with_metrics,
)
from summarize_comparison import summarize


DEFAULT_REPO = Path(__file__).resolve().parents[2]
DEFAULT_PACKAGE = Path(
    os.environ.get(
        "RAMAG_BENCHMARK_ROOT",
        DEFAULT_REPO / "benchmark-results/alignathon-sim-human-chimp-v1",
    )
)
DEFAULT_MUMMER_ROOT = Path(
    os.environ.get("MUMMER4_ROOT", DEFAULT_REPO.parent / "mummer-4.0.1")
)
EXPECTED_SUFKIT_COMMIT = "f8c4c386ee883e45ad0f973efc4c8e1148b0068a"
EXPECTED_SEQPRO_COMMIT = "6781cadcf81a0da53d7573444594c1484947017c"
CONFIG_SCHEMA = "ramag.human-chimp-preliminary-config.v3"
METRICS_SCHEMA = "ramag.command-metrics.v3"
EXISTING_RAMAG_DIAGNOSTIC_F1 = 0.9742776368
EXPECTED_FASTA_SHA256 = {
    "simHuman.fa": "66af39d61a3b9f4ddd0b202c2420dfcfe9521bf68ecbbe580b1debc5e2d7b8f9",
    "simChimp.fa": "fe1e3d5147fd50e832cd921e06963fedb36221c0b79c06b8a3200b884e8a3ed1",
}
EXPECTED_FAI_SHA256 = {
    "simHuman.fa.fai": "e359b4878aa56f189a96eebc9660b68bf0bdb3ae6138ff5607fdb9708cee97ec",
    "simChimp.fa.fai": "edc78a8f78e0a97db7174dd40672b333a190367d332b188b17899c6b69d3de1b",
}
EXPECTED_PACKAGE_MANIFEST_SHA256 = (
    "42cc2ea128f2f78d8d527686b29d6ba077a0f4a6b227e32051aeb47cbe822757"
)
EXPECTED_PACKAGE_SHA256SUMS_SHA256 = (
    "a43e5f34c9e40b3de4282c9650f05ee68b6abb26fff53ffa79789589cb40bd81"
)
EXPECTED_TRUTH_SHA256 = {
    "simHuman-simChimp.all-homology.maf.gz": (
        "5f649b375f56e4dd8419072751a854db8211e86f8c05c7b9c847119352ca1c6a"
    ),
    "simHuman-simChimp.no-paralogy.maf.gz": (
        "05c1297c624a948b7d4dba687fa1f01c68eb9f52b3bf843135396e94654e69b5"
    ),
    "simHuman-simChimp.single-copy-compat.maf.gz": (
        "62463797ff0cf1449a1509945d1d40887dbab7533f360f7097dc44b314a855f9"
    ),
}
EXPECTED_EVALUATOR_SHA256 = {
    "mafComparator": (
        "ded583c6b2c4eba5e408eebfcf8c6a4d307d724ec58b41757338b7f2efb39b16"
    ),
    "mafPairCounter": (
        "48a6930fef3d5e7322cea9f384e4b221b47e7cd27b3a68a1b4ccb31d5a08f207"
    ),
}
EXPECTED_EVALUATOR_VERSION_LINES = {
    "mafComparator": (
        "mafComparator, version 0.9 May 2013",
        "build: 2024-11-14T11:02CST, master, "
        "4e5b5de3f275f61b36b9762824cc1edbead31820",
    ),
    "mafPairCounter": (
        "mafPairCounter, version 0.1 July 2012",
        "build: 2024-11-14T11:02CST, master, "
        "4e5b5de3f275f61b36b9762824cc1edbead31820",
    ),
}
CAPS_MIN_REFERENCE_BASES = 64 * 1024 * 1024
THREADS = 16
REPETITIONS = 1
OPENMP_ENVIRONMENT = {
    "OMP_DYNAMIC": "FALSE",
    "OMP_MAX_ACTIVE_LEVELS": "1",
    "OMP_NUM_THREADS": str(THREADS),
    "OMP_PLACES": "CORES",
    "OMP_PROC_BIND": "SPREAD",
    "OMP_THREAD_LIMIT": str(THREADS),
}
SAMPLES = 1_000_000
NEAR = 0
SEEDS = [20260830, 20260831, 20260832]
MIN_AVAILABLE_MEMORY_BYTES = 20 * 1024**3
MIN_FREE_DISK_BYTES = 50 * 1024**3
FORBIDDEN_OLD_ROOTS = {
    "preliminary-ramag-vs-mummer4-20260830T045046Z",
    "preliminary-mam16-single-20260830T053945Z",
    "preliminary-mam16-single-20260830T083137Z",
    "preliminary-mam16-single-20260830T134638Z",
}
MAM_TILE_BASES = 4 * 1024 * 1024
SMOKE_MAM_TASK_SHAPE = {
    "oriented_queries": 2,
    "tile_tasks": 6,
    "boundary_tasks": 4,
    "total_tasks": 10,
}
FULL_MAM_TASK_SHAPE = {
    "oriented_queries": 8,
    "tile_tasks": 92,
    "boundary_tasks": 84,
    "total_tasks": 176,
}
EXPECTED_SEED_ROUTE = (
    "sufkit-full-sa:mumreference+tiled-mam-boundary-mem-v1+"
    "openmp-dynamic-stable-tile-boundary-tasks"
)
SMOKE_FASTA_CONTRACT = {
    "simHuman.chrD.fa": {
        "id": "simHuman.chrD",
        "length": 10_572_275,
        "sha256": "bd1dd33db481d2b110399bb1702c1b459dd3a3134576bae9edc33e488d283fdb",
    },
    "simChimp.chrD.fa": {
        "id": "simChimp.chrD",
        "length": 10_574_168,
        "sha256": "93b398ad9ed5411c196a0c24be876f05e7fdb48d2193a485e3d0bb9a0ff732ea",
    },
}


class ExperimentError(RuntimeError):
    """A user-facing stage or acceptance failure."""


class ClockQuarantineError(ExperimentError):
    """A successful command whose elapsed-clock evidence is not admissible."""


class AlignmentTimeoutError(ExperimentError):
    """A full alignment stopped at the same-run monotonic speed budget."""

    def __init__(
        self,
        tool: str,
        attempt: Path,
        metrics: dict[str, object],
    ) -> None:
        super().__init__(f"{tool} alignment timed out; preserved {attempt}")
        self.tool = tool
        self.attempt = attempt
        self.metrics = metrics


def fasta_record_lengths(path: Path) -> list[int]:
    """Read record lengths without materializing sequences.

    The immutable benchmark ships authoritative ``.fai`` files.  The chrD
    smoke copies do not, so the same contract has a strict streaming FASTA
    fallback.
    """

    fai = Path(str(path) + ".fai")
    if fai.is_file():
        lengths: list[int] = []
        seen: set[str] = set()
        for line_number, line in enumerate(
            fai.read_text(encoding="utf-8").splitlines(), start=1
        ):
            fields = line.split("\t")
            if len(fields) < 2 or not fields[0] or fields[0] in seen:
                raise ExperimentError(f"invalid FASTA index record at {fai}:{line_number}")
            try:
                length = int(fields[1])
            except ValueError as exc:
                raise ExperimentError(
                    f"invalid FASTA index length at {fai}:{line_number}"
                ) from exc
            if length <= 0:
                raise ExperimentError(
                    f"non-positive FASTA index length at {fai}:{line_number}"
                )
            seen.add(fields[0])
            lengths.append(length)
        if not lengths:
            raise ExperimentError(f"empty FASTA index: {fai}")
        return lengths

    lengths = []
    current_name: str | None = None
    current_length = 0
    seen: set[str] = set()
    with path.open("r", encoding="ascii", newline=None) as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith(">"):
                if current_name is not None:
                    if current_length <= 0:
                        raise ExperimentError(f"empty FASTA record {current_name!r}: {path}")
                    lengths.append(current_length)
                header = line[1:].strip()
                current_name = header.split(maxsplit=1)[0] if header else ""
                if not current_name or current_name in seen:
                    raise ExperimentError(
                        f"empty or duplicate FASTA ID at {path}:{line_number}"
                    )
                seen.add(current_name)
                current_length = 0
                continue
            sequence = line.strip()
            if current_name is None or not sequence or any(ch.isspace() for ch in sequence):
                raise ExperimentError(f"invalid FASTA sequence line at {path}:{line_number}")
            current_length += len(sequence)
    if current_name is None:
        raise ExperimentError(f"FASTA has no records: {path}")
    if current_length <= 0:
        raise ExperimentError(f"empty FASTA record {current_name!r}: {path}")
    lengths.append(current_length)
    return lengths


def expected_mam_task_shape(query: Path) -> dict[str, int]:
    lengths = fasta_record_lengths(query)
    tiles_per_forward_record = [
        (length + MAM_TILE_BASES - 1) // MAM_TILE_BASES for length in lengths
    ]
    tile_tasks = 2 * sum(tiles_per_forward_record)
    boundary_tasks = 2 * sum(tiles - 1 for tiles in tiles_per_forward_record)
    return {
        "oriented_queries": 2 * len(lengths),
        "tile_tasks": tile_tasks,
        "boundary_tasks": boundary_tasks,
        "total_tasks": tile_tasks + boundary_tasks,
    }


def utc_stamp() -> str:
    return dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def resolve_ramag_commit(repo: Path) -> str:
    """Return the clean RaMA-G Git HEAD, or ``unknown`` for a non-Git tree.

    A Git-backed formal run must name the repository root exactly and must not
    mix committed source with dirty or untracked files.  The ``unknown``
    fallback is retained only for preserved historical source snapshots that
    predate the server-side repository.
    """

    repo = repo.resolve()
    if not repo.is_dir():
        raise ExperimentError(f"RaMA-G source repository is not a directory: {repo}")

    try:
        top_level = subprocess.run(
            ["git", "-C", str(repo), "rev-parse", "--show-toplevel"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )
    except FileNotFoundError as exc:
        if (repo / ".git").exists():
            raise ExperimentError("git is unavailable for the RaMA-G repository") from exc
        return "unknown"
    if top_level.returncode != 0:
        return "unknown"

    reported_root = Path(top_level.stdout.strip()).resolve()
    if reported_root != repo:
        raise ExperimentError(
            f"--repo must name the RaMA-G Git worktree root exactly: {reported_root}"
        )

    head = subprocess.run(
        ["git", "-C", str(repo), "rev-parse", "--verify", "HEAD^{commit}"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )
    commit = head.stdout.strip().lower()
    if head.returncode != 0 or re.fullmatch(r"[0-9a-f]{40}", commit) is None:
        raise ExperimentError("RaMA-G Git repository has no valid 40-hex HEAD commit")

    status = subprocess.run(
        ["git", "-C", str(repo), "status", "--porcelain", "--untracked-files=all"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )
    if status.returncode != 0:
        raise ExperimentError("cannot inspect the RaMA-G Git worktree status")
    if status.stdout.strip():
        raise ExperimentError(
            "RaMA-G Git worktree must be clean before creating or resuming a result root"
        )
    return commit


def _read_key_value_marker(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for line_number, line in enumerate(
        path.read_text(encoding="utf-8").splitlines(), start=1
    ):
        if not line or "=" not in line:
            raise ExperimentError(f"invalid key=value marker at {path}:{line_number}")
        key, value = line.split("=", 1)
        if not key or key in values:
            raise ExperimentError(f"duplicate/empty marker key at {path}:{line_number}")
        values[key] = value
    return values


def _read_sha256sums(path: Path) -> dict[str, str]:
    entries: dict[str, str] = {}
    for line_number, line in enumerate(
        path.read_text(encoding="ascii").splitlines(), start=1
    ):
        if (
            len(line) < 69
            or line[64:68] != "  ./"
            or any(character not in "0123456789abcdef" for character in line[:64])
        ):
            raise ExperimentError(f"invalid SHA256SUMS record at {path}:{line_number}")
        relative_text = line[68:]
        relative = Path(relative_text)
        if (
            not relative_text
            or relative.is_absolute()
            or relative_text != relative.as_posix()
            or any(part in {"", ".", ".."} for part in relative.parts)
            or relative_text in entries
        ):
            raise ExperimentError(
                f"unsafe or duplicate SHA256SUMS path at {path}:{line_number}"
            )
        entries[relative_text] = line[:64]
    if not entries:
        raise ExperimentError(f"empty SHA256SUMS: {path}")
    return entries


def audit_benchmark_package_integrity(package: Path) -> dict[str, object]:
    """Revalidate the immutable benchmark package against fixed v1 anchors."""

    package = package.resolve()
    failures: list[str] = []
    marker_path = package / "PACKAGE_COMPLETE"
    manifest_path = package / "MANIFEST.json"
    sums_path = package / "SHA256SUMS"
    marker: dict[str, str] = {}
    sums: dict[str, str] = {}
    manifest: dict[str, object] = {}

    try:
        marker = _read_key_value_marker(marker_path)
    except (OSError, UnicodeError, ExperimentError) as error:
        failures.append(f"PACKAGE_COMPLETE is invalid: {error}")
    try:
        sums = _read_sha256sums(sums_path)
    except (OSError, UnicodeError, ExperimentError) as error:
        failures.append(f"SHA256SUMS is invalid: {error}")
    try:
        manifest = load_json(manifest_path)
    except ExperimentError as error:
        failures.append(f"MANIFEST.json is invalid: {error}")

    try:
        manifest_sha256 = sha256_file(manifest_path) if manifest_path.is_file() else None
    except OSError as error:
        manifest_sha256 = None
        failures.append(f"cannot hash MANIFEST.json: {error}")
    try:
        sums_sha256 = sha256_file(sums_path) if sums_path.is_file() else None
    except OSError as error:
        sums_sha256 = None
        failures.append(f"cannot hash SHA256SUMS: {error}")
    if manifest_sha256 != EXPECTED_PACKAGE_MANIFEST_SHA256:
        failures.append("MANIFEST.json does not match the fixed benchmark-v1 hash")
    if sums_sha256 != EXPECTED_PACKAGE_SHA256SUMS_SHA256:
        failures.append("SHA256SUMS does not match the fixed benchmark-v1 hash")
    if marker:
        if marker.get("status") != "validated":
            failures.append("PACKAGE_COMPLETE status is not validated")
        if marker.get("manifest_sha256") != EXPECTED_PACKAGE_MANIFEST_SHA256:
            failures.append("PACKAGE_COMPLETE records the wrong MANIFEST.json hash")
        if marker.get("sha256sums_sha256") != EXPECTED_PACKAGE_SHA256SUMS_SHA256:
            failures.append("PACKAGE_COMPLETE records the wrong SHA256SUMS hash")
    if sums:
        if sums.get("MANIFEST.json") != EXPECTED_PACKAGE_MANIFEST_SHA256:
            failures.append("SHA256SUMS records the wrong MANIFEST.json hash")
        for relative, expected in sums.items():
            artifact = package / relative
            if artifact.is_symlink() or not artifact.is_file():
                failures.append(f"missing/non-regular package artifact: {relative}")
                continue
            try:
                observed = sha256_file(artifact)
            except OSError as error:
                failures.append(f"cannot hash package artifact {relative}: {error}")
                continue
            if observed != expected:
                failures.append(f"package artifact SHA-256 mismatch: {relative}")

    artifacts = manifest.get("artifacts") if manifest else None
    if manifest:
        if (
            manifest.get("schema_version") != 1
            or manifest.get("status") != "validated"
            or manifest.get("package_name") != "alignathon-sim-human-chimp-v1"
            or manifest.get("headline_truth") != "all-homology"
        ):
            failures.append("MANIFEST.json does not describe the fixed validated package")
    if not isinstance(artifacts, dict):
        failures.append("MANIFEST.json lacks an artifacts object")
        artifacts = {}
    if sums and set(artifacts) != set(sums) - {"MANIFEST.json"}:
        failures.append("MANIFEST.json artifact set disagrees with SHA256SUMS")
    for relative, metadata in artifacts.items():
        expected = sums.get(relative)
        artifact = package / relative
        try:
            observed_bytes = artifact.stat().st_size if artifact.is_file() else -1
        except OSError:
            observed_bytes = -1
        if (
            not isinstance(metadata, dict)
            or metadata.get("sha256") != expected
            or not isinstance(metadata.get("bytes"), int)
            or observed_bytes < 0
            or int(metadata.get("bytes", -1)) != observed_bytes
        ):
            failures.append(f"MANIFEST.json artifact metadata mismatch: {relative}")

    truth_profiles = manifest.get("truth_profiles") if manifest else None
    if not isinstance(truth_profiles, dict):
        failures.append("MANIFEST.json lacks truth_profiles")
        truth_profiles = {}
    profile_to_filename = {
        "all-homology": "simHuman-simChimp.all-homology.maf.gz",
        "no-paralogy": "simHuman-simChimp.no-paralogy.maf.gz",
        "single-copy-compat": "simHuman-simChimp.single-copy-compat.maf.gz",
    }
    for profile, filename in profile_to_filename.items():
        relative = f"truth/{filename}"
        expected = EXPECTED_TRUTH_SHA256[filename]
        metadata = truth_profiles.get(profile)
        if sums.get(relative) != expected:
            failures.append(f"fixed truth hash is absent/wrong in SHA256SUMS: {profile}")
        if (
            not isinstance(metadata, dict)
            or metadata.get("artifact") != relative
            or metadata.get("sha256") != expected
        ):
            failures.append(f"MANIFEST.json truth profile mismatch: {profile}")

    evaluator_versions: dict[str, object] = {}
    for filename, expected_hash in EXPECTED_EVALUATOR_SHA256.items():
        relative = f"tools/bin/{filename}"
        executable = package / relative
        if sums.get(relative) != expected_hash:
            failures.append(f"fixed evaluator hash is absent/wrong: {filename}")
        if (
            not executable.is_file()
            or executable.is_symlink()
            or not os.access(executable, os.X_OK)
        ):
            failures.append(f"fixed evaluator is missing/not executable: {filename}")
            continue
        if sha256_file(executable) != expected_hash:
            # Never execute a binary that failed the fixed-hash gate.
            continue
        try:
            completed = subprocess.run(
                [str(executable), "--version"],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                encoding="utf-8",
                errors="replace",
                env={**os.environ, "LC_ALL": "C"},
                timeout=15,
                check=False,
            )
        except (OSError, subprocess.SubprocessError) as error:
            failures.append(f"cannot execute {filename} --version: {error}")
            continue
        version_lines = tuple(
            line.strip()
            for line in (completed.stdout + completed.stderr).splitlines()
            if line.strip()
        )
        evaluator_versions[filename] = {
            "exit_code": completed.returncode,
            "lines": list(version_lines),
        }
        if (
            completed.returncode != 0
            or version_lines != EXPECTED_EVALUATOR_VERSION_LINES[filename]
        ):
            failures.append(f"fixed evaluator version mismatch: {filename}")

    return {
        "schema": "ramag.human-chimp-package-integrity.v1",
        "package": str(package),
        "status": "failed" if failures else "success",
        "failures": failures,
        "manifest_sha256": manifest_sha256,
        "sha256sums_sha256": sums_sha256,
        "artifact_count": len(sums),
        "truth_sha256": dict(EXPECTED_TRUTH_SHA256),
        "evaluator_sha256": dict(EXPECTED_EVALUATOR_SHA256),
        "evaluator_versions": evaluator_versions,
    }


def audit_fixed_input_integrity(
    config: dict[str, object],
    *,
    reference: Path | None = None,
    query: Path | None = None,
    include_fasta: bool = True,
    include_fai: bool = True,
) -> dict[str, object]:
    """Hash fixed FASTA/FAI inputs against constants and accepted config."""

    failures: list[str] = []
    accepted_fasta = config.get("input_sha256")
    accepted_fai = config.get("input_fai_sha256")
    if accepted_fasta != EXPECTED_FASTA_SHA256:
        failures.append("accepted config does not contain the fixed FASTA hashes")
    if accepted_fai != EXPECTED_FAI_SHA256:
        failures.append("accepted config does not contain the fixed FAI hashes")

    try:
        configured_reference = Path(str(config["reference"])).resolve()
        configured_query = Path(str(config["query"])).resolve()
        package = Path(str(config["package"])).resolve()
    except KeyError as error:
        failures.append(f"accepted config lacks fixed input path: {error.args[0]}")
        configured_reference = reference.resolve() if reference is not None else Path()
        configured_query = query.resolve() if query is not None else Path()
        package = Path()
    expected_reference = package / "inputs/simHuman.fa"
    expected_query = package / "inputs/simChimp.fa"
    if configured_reference != expected_reference:
        failures.append("accepted reference path is not package inputs/simHuman.fa")
    if configured_query != expected_query:
        failures.append("accepted query path is not package inputs/simChimp.fa")
    if reference is not None and reference.resolve() != configured_reference:
        failures.append("runtime reference path differs from accepted config")
    if query is not None and query.resolve() != configured_query:
        failures.append("runtime query path differs from accepted config")

    observations: dict[str, object] = {}
    targets: list[tuple[str, Path, str]] = []
    if include_fasta:
        targets.extend(
            (
                ("simHuman.fa", configured_reference, EXPECTED_FASTA_SHA256["simHuman.fa"]),
                ("simChimp.fa", configured_query, EXPECTED_FASTA_SHA256["simChimp.fa"]),
            )
        )
    if include_fai:
        targets.extend(
            (
                (
                    "simHuman.fa.fai",
                    Path(str(configured_reference) + ".fai"),
                    EXPECTED_FAI_SHA256["simHuman.fa.fai"],
                ),
                (
                    "simChimp.fa.fai",
                    Path(str(configured_query) + ".fai"),
                    EXPECTED_FAI_SHA256["simChimp.fa.fai"],
                ),
            )
        )
    for filename, path, expected_hash in targets:
        if path.is_symlink() or not path.is_file():
            failures.append(f"fixed input is missing/non-regular: {path}")
            continue
        try:
            before = _stable_file_identity(path.stat())
            observed_hash = sha256_file(path)
            after = _stable_file_identity(path.stat())
        except OSError as error:
            failures.append(f"cannot hash fixed input {path}: {error}")
            continue
        observations[filename] = {
            "path": str(path.resolve()),
            "bytes": after["bytes"],
            "expected_sha256": expected_hash,
            "observed_sha256": observed_hash,
            "stable_identity": before == after,
        }
        if before != after:
            failures.append(f"fixed input changed while hashing: {path}")
        if observed_hash != expected_hash:
            failures.append(f"fixed input SHA-256 mismatch: {filename}")
    return {
        "schema": "ramag.human-chimp-fixed-input-integrity.v1",
        "status": "failed" if failures else "success",
        "failures": failures,
        "include_fasta": include_fasta,
        "include_fai": include_fai,
        "observations": observations,
    }


def require_fixed_input_integrity(
    config: dict[str, object],
    *,
    reference: Path | None = None,
    query: Path | None = None,
    include_fasta: bool = True,
    include_fai: bool = True,
) -> dict[str, object]:
    report = audit_fixed_input_integrity(
        config,
        reference=reference,
        query=query,
        include_fasta=include_fasta,
        include_fai=include_fai,
    )
    if report["status"] != "success":
        raise ExperimentError(
            "fixed FASTA/FAI integrity gate failed: " + "; ".join(report["failures"])
        )
    return report


def _utc_from_realtime_ns(realtime_ns: int) -> str:
    seconds, nanoseconds = divmod(realtime_ns, 1_000_000_000)
    instant = dt.datetime.fromtimestamp(seconds, tz=dt.timezone.utc).replace(
        microsecond=nanoseconds // 1_000
    )
    return instant.isoformat(timespec="microseconds").replace("+00:00", "Z")


def _stable_file_identity(stat_result: os.stat_result) -> dict[str, int]:
    return {
        "device": int(stat_result.st_dev),
        "inode": int(stat_result.st_ino),
        "bytes": int(stat_result.st_size),
        "mtime_ns": int(stat_result.st_mtime_ns),
    }


def _sequential_hash_evidence(
    role: str,
    path: Path,
    expected_sha256: str | None,
) -> dict[str, object]:
    started_monotonic_ns = time.monotonic_ns()
    started_realtime_ns = time.time_ns()
    digest = hashlib.sha256()
    bytes_read = 0
    failure: str | None = None
    before: dict[str, int] | None = None
    after: dict[str, int] | None = None
    try:
        before = _stable_file_identity(path.stat())
        with path.open("rb") as handle:
            while chunk := handle.read(8 * 1024 * 1024):
                digest.update(chunk)
                bytes_read += len(chunk)
        after = _stable_file_identity(path.stat())
    except OSError as exc:
        failure = f"{type(exc).__name__}: {exc}"
    finished_monotonic_ns = time.monotonic_ns()
    finished_realtime_ns = time.time_ns()
    observed_sha256 = digest.hexdigest() if failure is None else None
    if expected_sha256 is None:
        failure = failure or "accepted config has no expected SHA-256"
    elif observed_sha256 != expected_sha256:
        failure = failure or "SHA-256 mismatch"
    if before is not None and after is not None:
        if before != after:
            failure = failure or "file identity changed during sequential read"
        elif bytes_read != after["bytes"]:
            failure = failure or "complete byte count does not match file size"
    return {
        "role": role,
        "path": str(path.resolve()),
        "filename": path.name,
        "read_order": 1 if role == "reference" else 2,
        "read_method": "single-open sequential 8-MiB chunks through EOF",
        "expected_sha256": expected_sha256,
        "observed_sha256": observed_sha256,
        "bytes_read": bytes_read,
        "file_identity_before": before,
        "file_identity_after": after,
        "started_monotonic_ns": started_monotonic_ns,
        "finished_monotonic_ns": finished_monotonic_ns,
        "monotonic_elapsed_seconds": (
            finished_monotonic_ns - started_monotonic_ns
        )
        / 1_000_000_000.0,
        "started_utc": _utc_from_realtime_ns(started_realtime_ns),
        "finished_utc": _utc_from_realtime_ns(finished_realtime_ns),
        "failure": failure,
        "status": "success" if failure is None else "failed",
    }


def warm_cache_inputs_for_tool(
    tool: str,
    attempt: Path,
    reference: Path,
    query: Path,
    config: dict[str, object],
) -> dict[str, object]:
    """Sequentially hash both fixed inputs immediately before one full run."""

    if tool not in {"mummer4", "ramag"}:
        raise ExperimentError(f"unknown warm-cache tool: {tool}")
    expected = config.get("input_sha256")
    if not isinstance(expected, dict):
        expected = {}
    reference_expected = expected.get(reference.name)
    query_expected = expected.get(query.name)
    started_monotonic_ns = time.monotonic_ns()
    started_realtime_ns = time.time_ns()
    inputs = [
        _sequential_hash_evidence(
            "reference",
            reference,
            reference_expected if isinstance(reference_expected, str) else None,
        ),
        _sequential_hash_evidence(
            "query",
            query,
            query_expected if isinstance(query_expected, str) else None,
        ),
    ]
    index_integrity = audit_fixed_input_integrity(
        config,
        reference=reference,
        query=query,
        include_fasta=False,
        include_fai=True,
    )
    finished_monotonic_ns = time.monotonic_ns()
    finished_realtime_ns = time.time_ns()
    failures = [
        f"{entry['role']}: {entry['failure']}"
        for entry in inputs
        if entry["status"] != "success"
    ]
    failures.extend(str(failure) for failure in index_integrity["failures"])
    report = {
        "schema": "ramag.human-chimp-per-tool-warm-cache.v1",
        "tool": tool,
        "read_order": ["reference", "query"],
        "purpose": (
            "per-tool verified warm-cache starting condition; no privileged cache drop"
        ),
        "inputs": inputs,
        "input_index_integrity": index_integrity,
        "started_monotonic_ns": started_monotonic_ns,
        "finished_monotonic_ns": finished_monotonic_ns,
        "monotonic_elapsed_seconds": (
            finished_monotonic_ns - started_monotonic_ns
        )
        / 1_000_000_000.0,
        "started_utc": _utc_from_realtime_ns(started_realtime_ns),
        "finished_utc": _utc_from_realtime_ns(finished_realtime_ns),
        "failures": failures,
        "status": "success" if not failures else "failed",
    }
    write_json_exclusive(attempt / "warm-cache.json", report)
    marker = "WARM_CACHE_COMPLETE" if not failures else "WARM_CACHE_FAILED"
    (attempt / marker).write_text(
        "success\n" if not failures else "\n".join(failures) + "\n",
        encoding="utf-8",
    )
    if failures:
        raise ExperimentError(
            f"{tool} per-tool warm-cache hash gate failed; aligner was not launched; "
            f"inspect {attempt / 'warm-cache.json'}"
        )
    return report


def record_warm_cache_launch_relation(
    attempt: Path,
    warm_cache: dict[str, object],
    metrics: dict[str, object],
) -> dict[str, object]:
    measurements = metrics.get("clock_measurements")
    warm_finished = warm_cache.get("finished_monotonic_ns")
    runner_started = (
        measurements.get("monotonic_started_ns")
        if isinstance(measurements, dict)
        else None
    )
    valid = (
        isinstance(warm_finished, int)
        and isinstance(runner_started, int)
        and runner_started >= warm_finished
    )
    gap = (
        (runner_started - warm_finished) / 1_000_000_000.0
        if valid
        else None
    )
    report = {
        "schema": "ramag.human-chimp-warm-cache-launch-relation.v1",
        "warm_cache_finished_monotonic_ns": warm_finished,
        "metrics_runner_started_monotonic_ns": runner_started,
        "warm_cache_to_metrics_runner_seconds": gap,
        "code_path": (
            "warm-cache evidence publication followed directly by run_with_metrics"
        ),
        "status": "passed" if valid else "quarantined",
    }
    write_json_exclusive(attempt / "warm-cache-launch.json", report)
    if not valid:
        (attempt / "WARM_CACHE_LAUNCH_QUARANTINED").write_text(
            "warm-cache and runner monotonic timestamps are inconsistent\n",
            encoding="ascii",
        )
        raise ClockQuarantineError(
            "warm-cache evidence cannot be ordered immediately before the aligner"
        )
    return report


def write_json_exclusive(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8", newline="\n") as handle:
        json.dump(value, handle, indent=2, sort_keys=True, ensure_ascii=False)
        handle.write("\n")


def load_json(path: Path) -> dict[str, object]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ExperimentError(f"cannot load JSON: {path}") from exc
    if not isinstance(value, dict):
        raise ExperimentError(f"expected a JSON object: {path}")
    return value


def run_capture(
    command: list[str],
    stdout: Path,
    stderr: Path,
    environment_overrides: dict[str, str] | None = None,
) -> int:
    stdout.parent.mkdir(parents=True, exist_ok=True)
    environment = os.environ.copy()
    environment["LC_ALL"] = "C"
    environment.update(environment_overrides or {})
    with stdout.open("x", encoding="utf-8", newline="\n") as out, stderr.open(
        "x", encoding="utf-8", newline="\n"
    ) as err:
        completed = subprocess.run(
            command, stdout=out, stderr=err, env=environment, check=False
        )
    return completed.returncode


def command_evidence(
    command: list[str], output_dir: Path, name: str, *, discard_stdout: bool = False
) -> None:
    command_path = output_dir / f"{name}.command.json"
    stdout_path = output_dir / f"{name}.stdout.txt"
    stderr_path = output_dir / f"{name}.stderr.txt"
    write_json_exclusive(command_path, command)
    started_monotonic_ns = time.monotonic_ns()
    started_realtime_ns = time.time_ns()
    if discard_stdout:
        with stderr_path.open("x", encoding="utf-8", newline="\n") as stderr:
            completed = subprocess.run(
                command, stdout=subprocess.DEVNULL, stderr=stderr, check=False
            )
        stdout_path.write_text(
            "stdout intentionally discarded after parser exit-status validation\n",
            encoding="ascii",
        )
    else:
        completed_code = run_capture(command, stdout_path, stderr_path)
        completed = type("Completed", (), {"returncode": completed_code})()
    finished_monotonic_ns = time.monotonic_ns()
    finished_realtime_ns = time.time_ns()
    write_json_exclusive(
        output_dir / f"{name}.timing.json",
        {
            "schema": "ramag.benchmark-postprocess-timing.v1",
            "name": name,
            "clock": "CLOCK_MONOTONIC",
            "started_monotonic_ns": started_monotonic_ns,
            "finished_monotonic_ns": finished_monotonic_ns,
            "monotonic_elapsed_seconds": (
                finished_monotonic_ns - started_monotonic_ns
            )
            / 1_000_000_000.0,
            "started_utc": _utc_from_realtime_ns(started_realtime_ns),
            "finished_utc": _utc_from_realtime_ns(finished_realtime_ns),
            "exit_code": completed.returncode,
            "status": "success" if completed.returncode == 0 else "failed",
        },
    )
    (output_dir / f"{name}.exit-code.txt").write_text(
        f"{completed.returncode}\n", encoding="ascii"
    )
    if completed.returncode != 0:
        raise ExperimentError(
            f"command failed ({name}, exit {completed.returncode}); see {stderr_path}"
        )


def next_attempt(base: Path) -> tuple[str, Path]:
    base.mkdir(parents=True, exist_ok=True)
    indices: list[int] = []
    for child in base.iterdir():
        if child.is_dir() and child.name.startswith("attempt-"):
            try:
                indices.append(int(child.name.removeprefix("attempt-")))
            except ValueError:
                continue
    index = max(indices, default=0) + 1
    name = f"attempt-{index:03d}"
    attempt = base / name
    attempt.mkdir()
    return name, attempt


def selected_attempt(base: Path, required_marker: str) -> Path | None:
    selector = base / "SELECTED_ATTEMPT"
    if not selector.exists():
        return None
    name = selector.read_text(encoding="ascii").strip()
    if not name.startswith("attempt-") or "/" in name or "\\" in name:
        raise ExperimentError(f"invalid selected attempt marker: {selector}")
    attempt = base / name
    if not (attempt / required_marker).is_file():
        raise ExperimentError(f"selected attempt is incomplete: {attempt}")
    if required_marker == "RUN_ACCEPTED" and (attempt / "timing/RUN_TIMED_OUT").is_file():
        raise ExperimentError(f"timed-out attempt cannot be selected: {attempt}")
    return attempt


def accept_attempt(base: Path, name: str, attempt: Path, marker: str) -> None:
    if marker == "RUN_ACCEPTED" and (attempt / "timing/RUN_TIMED_OUT").is_file():
        raise ExperimentError(f"timed-out attempt cannot be accepted: {attempt}")
    (attempt / marker).write_text("success\n", encoding="ascii")
    selector = base / "SELECTED_ATTEMPT"
    with selector.open("x", encoding="ascii", newline="\n") as handle:
        handle.write(name + "\n")


def copy_source_snapshot(repo: Path, destination: Path) -> dict[str, object]:
    ignored_directories = {
        ".git",
        ".agents",
        ".codex",
        ".cache",
        "__pycache__",
        "_deps",
        "deps-src",
        "result",
        "results",
        "temp",
        "tmp",
        "work",
    }
    destination.mkdir(parents=True)
    entries: list[tuple[str, str, int]] = []
    for source in sorted(repo.rglob("*")):
        relative = source.relative_to(repo)
        if any(
            part in ignored_directories or part.startswith("build")
            for part in relative.parts
        ):
            continue
        target = destination / relative
        if source.is_dir():
            target.mkdir(parents=True, exist_ok=True)
        elif source.is_file() and not source.is_symlink():
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, target)
            digest = sha256_file(target)
            entries.append((relative.as_posix(), digest, target.stat().st_size))
    manifest = destination.parent / "ramag-source-files.sha256.tsv"
    with manifest.open("x", encoding="utf-8", newline="\n") as handle:
        handle.write("relative_path\tsha256\tbytes\n")
        for relative, digest, size in entries:
            handle.write(f"{relative}\t{digest}\t{size}\n")
    normalized = hashlib.sha256()
    for relative, digest, size in entries:
        normalized.update(f"{relative}\0{digest}\0{size}\n".encode("utf-8"))
    return {
        "file_count": len(entries),
        "tree_sha256": normalized.hexdigest(),
        "manifest": str(manifest.resolve()),
    }


def fingerprint_tree(source_root: Path, manifest: Path) -> dict[str, object]:
    """Hash a current local tree without copying or modifying it."""

    entries: list[tuple[str, str, int]] = []
    for path in sorted(source_root.rglob("*")):
        if not path.is_file() or path.is_symlink() or ".git" in path.parts:
            continue
        relative = path.relative_to(source_root).as_posix()
        entries.append((relative, sha256_file(path), path.stat().st_size))
    normalized = hashlib.sha256()
    with manifest.open("x", encoding="utf-8", newline="\n") as handle:
        handle.write("relative_path\tsha256\tbytes\n")
        for relative, digest, size in entries:
            handle.write(f"{relative}\t{digest}\t{size}\n")
            normalized.update(f"{relative}\0{digest}\0{size}\n".encode("utf-8"))
    return {
        "root": str(source_root.resolve()),
        "file_count": len(entries),
        "tree_sha256": normalized.hexdigest(),
        "manifest": str(manifest.resolve()),
    }


def read_mem_available() -> int:
    for line in Path("/proc/meminfo").read_text(encoding="ascii").splitlines():
        if line.startswith("MemAvailable:"):
            return int(line.split()[1]) * 1024
    raise ExperimentError("cannot read MemAvailable from /proc/meminfo")


def active_aligners() -> list[dict[str, object]]:
    found: list[dict[str, object]] = []
    for entry in Path("/proc").iterdir():
        if not entry.name.isdigit() or int(entry.name) == os.getpid():
            continue
        try:
            executable = Path(os.readlink(entry / "exe")).name
            command = (entry / "cmdline").read_bytes().replace(b"\0", b" ").decode(
                "utf-8", errors="replace"
            )
        except (FileNotFoundError, PermissionError, ProcessLookupError, OSError):
            continue
        if executable in {"ramag", "nucmer", "lt-nucmer"}:
            found.append({"pid": int(entry.name), "executable": executable, "command": command})
    return sorted(found, key=lambda value: int(value["pid"]))


def record_full_launch_gate(
    root: Path, tool: str, config: dict[str, object]
) -> Path:
    """Recheck and preserve host resources immediately before a full aligner.

    A successful check is deliberately not reused: if an alignment attempt
    fails, the next resume records a new check immediately before launching a
    replacement attempt. Already accepted alignment runs never call this
    function again.
    """

    if tool not in {"mummer4", "ramag"}:
        raise ExperimentError(f"unknown full-run launch-gate tool: {tool}")
    name, attempt = next_attempt(root / "runs/launch-gates" / tool)
    available_memory = read_mem_available()
    free_disk = shutil.disk_usage(root).free
    conflicts = active_aligners()
    failures: list[str] = []
    input_integrity = audit_fixed_input_integrity(config)
    failures.extend(str(failure) for failure in input_integrity["failures"])
    if conflicts:
        failures.append(f"other ramag/nucmer processes are active: {conflicts}")
    if available_memory < MIN_AVAILABLE_MEMORY_BYTES:
        failures.append(
            f"available memory below 20 GiB: {available_memory} bytes"
        )
    if free_disk < MIN_FREE_DISK_BYTES:
        failures.append(f"free disk below 50 GiB: {free_disk} bytes")
    report = {
        "schema": "ramag.human-chimp-full-launch-gate.v1",
        "attempt": name,
        "tool": tool,
        "checked_at_utc": dt.datetime.now(dt.timezone.utc)
        .isoformat(timespec="seconds")
        .replace("+00:00", "Z"),
        "status": "failed" if failures else "success",
        "failures": failures,
        "active_aligners": conflicts,
        "fixed_input_integrity": input_integrity,
        "available_memory_bytes": available_memory,
        "required_available_memory_bytes": MIN_AVAILABLE_MEMORY_BYTES,
        "free_disk_bytes": free_disk,
        "required_free_disk_bytes": MIN_FREE_DISK_BYTES,
    }
    write_json_exclusive(attempt / "launch-gate.json", report)
    marker = "LAUNCH_GATE_FAILED" if failures else "LAUNCH_GATE_ACCEPTED"
    (attempt / marker).write_text(
        ("\n".join(failures) + "\n") if failures else "success\n",
        encoding="utf-8",
    )
    if failures:
        raise ExperimentError(
            f"{tool} full-run launch gate failed; inspect {attempt / 'launch-gate.json'}"
        )
    return attempt


def snapshot_binary(source: Path, destination: Path) -> dict[str, object]:
    if not source.is_file() or not os.access(source, os.X_OK):
        raise ExperimentError(f"missing executable binary: {source}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        raise ExperimentError(f"refusing to overwrite binary snapshot: {destination}")
    shutil.copy2(source, destination)
    return {
        "source": str(source.resolve()),
        "snapshot": str(destination.resolve()),
        "source_sha256": sha256_file(source),
        "snapshot_sha256": sha256_file(destination),
    }


def create_result_root(args: argparse.Namespace) -> tuple[Path, dict[str, object]]:
    package = args.package.resolve()
    output_parent = package / "results"
    source_binary = args.ramag_binary.resolve()
    repo = args.repo.resolve()
    if not source_binary.is_file() or not os.access(source_binary, os.X_OK):
        raise ExperimentError(
            "--ramag-binary must name the fresh, executable Release binary: "
            f"{source_binary}"
        )
    ramag_commit = resolve_ramag_commit(repo)
    if args.result_root is None:
        root = output_parent / f"preliminary-mam16-single-{utc_stamp()}"
    else:
        root = args.result_root.resolve()
    if root.name in FORBIDDEN_OLD_ROOTS:
        raise ExperimentError(f"the preserved incomplete root must not be reused: {root}")
    if root.parent != output_parent:
        raise ExperimentError(f"result root must be a direct child of {output_parent}")
    if root.exists():
        raise ExperimentError(f"refusing to overwrite existing result root: {root}")
    root.mkdir(parents=True)
    for directory in ("provenance", "build", "smoke", "runs", "evaluation", "scripts"):
        (root / directory).mkdir()

    scripts_source = Path(__file__).resolve().parent
    for name in (
        "run_human_chimp_benchmark.py",
        "run_with_metrics.py",
        "canonicalize_maf_sources.py",
        "summarize_comparison.py",
    ):
        shutil.copy2(scripts_source / name, root / "scripts" / name)
    harness_scripts = {
        path.name: sha256_file(path)
        for path in sorted((root / "scripts").iterdir())
        if path.is_file()
    }

    binary_info = snapshot_binary(source_binary, root / "build/ramag")
    config: dict[str, object] = {
        "schema": CONFIG_SCHEMA,
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat().replace("+00:00", "Z"),
        "package": str(package),
        "reference": str(package / "inputs/simHuman.fa"),
        "query": str(package / "inputs/simChimp.fa"),
        "repo": str(repo),
        "ramag_source_binary": str(source_binary),
        "ramag_run_binary": str((root / "build/ramag").resolve()),
        "mummer_root": str(args.mummer_root.resolve()),
        "mummer_nucmer": str((args.mummer_root / "nucmer").resolve()),
        "mummer_show_coords": str((args.mummer_root / "show-coords").resolve()),
        "mummer_delta_filter": str((args.mummer_root / "delta-filter").resolve()),
        "threads": THREADS,
        "repetitions": REPETITIONS,
        "openmp_environment": OPENMP_ENVIRONMENT,
        "seed_semantics": "reference-unique MAM",
        "ramag_seed_mode": "mumreference",
        "ramag_selection_mode": "one-to-one",
        "min_match": 20,
        "max_gap": 90,
        "diag_diff": 5,
        "diag_factor": 0.12,
        "min_cluster": 65,
        "break_length": 200,
        "max_dp_cells": 4_000_000,
        "truth_profile": "all-homology",
        "samples": SAMPLES,
        "near": NEAR,
        "seeds": SEEDS,
        "sample_interval_seconds": args.sample_interval,
        "clock_consistency_tolerance_fraction": (
            CLOCK_CONSISTENCY_TOLERANCE_FRACTION
        ),
        "full_run_warm_cache_policy": (
            "reference-then-query full sequential SHA-256 immediately before each tool"
        ),
        "expected_sufkit_commit": EXPECTED_SUFKIT_COMMIT,
        "expected_seqpro_commit": EXPECTED_SEQPRO_COMMIT,
        "input_sha256": EXPECTED_FASTA_SHA256,
        "input_fai_sha256": EXPECTED_FAI_SHA256,
        "package_manifest_sha256": EXPECTED_PACKAGE_MANIFEST_SHA256,
        "package_sha256sums_sha256": EXPECTED_PACKAGE_SHA256SUMS_SHA256,
        "truth_sha256": EXPECTED_TRUTH_SHA256,
        "evaluator_sha256": EXPECTED_EVALUATOR_SHA256,
        "runtime_endpoint": "FASTA read + index build + alignment + raw delta write",
        "runtime_order": ["mummer4", "ramag"],
        "runtime_is_single_observation": True,
        "ramag_commit": ramag_commit,
        "metrics_schema": METRICS_SCHEMA,
        "early_stop_policy": {
            "enabled": True,
            "cutoff_source": (
                "same result root's accepted MUMmer4 CLOCK_MONOTONIC elapsed"
            ),
            "timeout_condition": "ramag_monotonic_elapsed_seconds > cutoff_seconds",
            "equal_is_accepted": True,
            "termination_scope": "runner-created exact process group only",
            "termination_sequence": "SIGTERM; wait at most 10 seconds; SIGKILL if live",
            "partial_evidence_policy": "preserve all files; never select or resume timed-out attempt",
        },
        "alignment_speed_policy": {
            "gnu_requirement": "ramag_gnu_elapsed_seconds <= mummer4_gnu_elapsed_seconds",
            "monotonic_requirement": (
                "ramag_monotonic_elapsed_seconds <= mummer4_monotonic_elapsed_seconds"
            ),
            "equal_is_accepted": True,
            "observations_per_tool": 1,
            "claim_scope": "preliminary single observation",
        },
        "quality_policy": {
            "mode": "report-only",
            "minimum_f1": None,
            "comparative_f1_gate": False,
            "format_coordinate_pair_count_reproducibility_are_hard_gates": True,
            "existing_ramag_diagnostic_f1": EXISTING_RAMAG_DIAGNOSTIC_F1,
        },
        "binary_snapshot": binary_info,
        "harness_scripts_sha256": harness_scripts,
    }
    write_json_exclusive(root / "accepted-config.json", config)
    return root, config


def load_resume_root(args: argparse.Namespace) -> tuple[Path, dict[str, object]]:
    if args.result_root is None:
        raise ExperimentError("--resume requires --result-root")
    root = args.result_root.resolve()
    if root.name in FORBIDDEN_OLD_ROOTS:
        raise ExperimentError(f"the preserved incomplete root must not be reused: {root}")
    if not root.is_dir():
        raise ExperimentError(f"resume root does not exist: {root}")
    config = load_json(root / "accepted-config.json")
    if config.get("schema") != CONFIG_SCHEMA:
        raise ExperimentError("resume config predates the early alignment-speed gate")
    package = Path(str(config["package"])).resolve()
    if root.parent != package / "results":
        raise ExperimentError("resume root is outside its accepted package results directory")
    if int(config.get("threads", 0)) != THREADS or int(config.get("repetitions", 0)) != 1:
        raise ExperimentError("resume config is not the fixed 16-thread, one-run protocol")
    if config.get("openmp_environment") != OPENMP_ENVIRONMENT:
        raise ExperimentError("resume config does not preserve the fixed OpenMP environment")
    if config.get("input_sha256") != EXPECTED_FASTA_SHA256 or config.get(
        "input_fai_sha256"
    ) != EXPECTED_FAI_SHA256:
        raise ExperimentError("resume config does not preserve fixed FASTA/FAI hashes")
    if (
        config.get("package_manifest_sha256") != EXPECTED_PACKAGE_MANIFEST_SHA256
        or config.get("package_sha256sums_sha256")
        != EXPECTED_PACKAGE_SHA256SUMS_SHA256
        or config.get("truth_sha256") != EXPECTED_TRUTH_SHA256
        or config.get("evaluator_sha256") != EXPECTED_EVALUATOR_SHA256
    ):
        raise ExperimentError("resume config does not preserve fixed package anchors")
    if float(config.get("sample_interval_seconds", 0.0)) != 1.0:
        raise ExperimentError("resume config does not preserve one-second process sampling")
    if float(config.get("clock_consistency_tolerance_fraction", -1.0)) != (
        CLOCK_CONSISTENCY_TOLERANCE_FRACTION
    ):
        raise ExperimentError("resume config does not preserve the 1% clock gate")
    if (
        config.get("expected_sufkit_commit") != EXPECTED_SUFKIT_COMMIT
        or config.get("expected_seqpro_commit") != EXPECTED_SEQPRO_COMMIT
        or config.get("metrics_schema") != METRICS_SCHEMA
    ):
        raise ExperimentError("resume config does not preserve source/dependency identity")
    accepted_ramag_commit = config.get("ramag_commit")
    if not isinstance(accepted_ramag_commit, str) or (
        accepted_ramag_commit != "unknown"
        and re.fullmatch(r"[0-9a-f]{40}", accepted_ramag_commit) is None
    ):
        raise ExperimentError("resume config has an invalid RaMA-G commit identity")
    accepted_repo = config.get("repo")
    if accepted_repo is None:
        if accepted_ramag_commit != "unknown":
            raise ExperimentError("resume config lacks the RaMA-G repository path")
    else:
        current_ramag_commit = resolve_ramag_commit(Path(str(accepted_repo)))
        if current_ramag_commit != accepted_ramag_commit:
            raise ExperimentError(
                "resume RaMA-G repository HEAD differs from the accepted commit"
            )
    early_stop = config.get("early_stop_policy")
    speed_policy = config.get("alignment_speed_policy")
    quality_policy = config.get("quality_policy")
    if (
        not isinstance(early_stop, dict)
        or early_stop.get("enabled") is not True
        or early_stop.get("equal_is_accepted") is not True
        or not isinstance(speed_policy, dict)
        or speed_policy.get("equal_is_accepted") is not True
        or speed_policy.get("observations_per_tool") != 1
        or not isinstance(quality_policy, dict)
        or quality_policy.get("mode") != "report-only"
        or quality_policy.get("minimum_f1") is not None
        or quality_policy.get("comparative_f1_gate") is not False
    ):
        raise ExperimentError("resume config does not preserve speed/F1 publication policy")
    require_fixed_input_integrity(config)
    accepted_scripts = config.get("harness_scripts_sha256")
    if not isinstance(accepted_scripts, dict):
        raise ExperimentError("resume config lacks accepted harness script hashes")
    current_scripts = Path(__file__).resolve().parent
    for filename, expected in accepted_scripts.items():
        current = current_scripts / str(filename)
        if not current.is_file() or sha256_file(current) != expected:
            raise ExperimentError(
                "resume harness differs from the accepted snapshot; execute "
                f"{root / 'scripts/run_human_chimp_benchmark.py'}"
            )
    return root, config


def collect_provenance(root: Path, config: dict[str, object]) -> dict[str, object]:
    capture_base = root / "provenance/captures"
    selected = selected_attempt(capture_base, "PROVENANCE_ACCEPTED")
    if selected is not None:
        return load_json(selected / "summary.json")
    name, provenance = next_attempt(capture_base)
    commands = {
        "host.uname.txt": ["uname", "-a"],
        "host.lscpu.txt": ["lscpu"],
        "host.free.txt": ["free", "-b"],
        "host.df.txt": ["df", "-B1", str(root)],
        "ramag.version.txt": [str(config["ramag_run_binary"]), "--version"],
        "ramag.ldd.txt": ["ldd", str(config["ramag_run_binary"])],
        "ramag.build-id.txt": ["readelf", "-n", str(config["ramag_run_binary"])],
        "mummer4.version.txt": [str(config["mummer_nucmer"]), "--version"],
        "mummer4.ldd.txt": ["ldd", str(Path(str(config["mummer_root"])) / ".libs/nucmer")],
        "mummer4.build-id.txt": [
            "readelf",
            "-n",
            str(Path(str(config["mummer_root"])) / ".libs/nucmer"),
        ],
    }
    command_status: dict[str, int] = {}
    for filename, command in commands.items():
        stdout = provenance / filename
        stderr = provenance / (filename + ".stderr")
        command_status[filename] = run_capture(
            command,
            stdout,
            stderr,
            OPENMP_ENVIRONMENT if filename == "ramag.version.txt" else None,
        )

    ramag_version = (provenance / "ramag.version.txt").read_text(
        encoding="utf-8", errors="replace"
    )
    expected_ramag_commit = str(config.get("ramag_commit", ""))
    version_commit = re.search(
        r"\(commit ([0-9a-f]{40}|unknown(?:-status-unknown)?)(-dirty)?;",
        ramag_version,
    )
    if (
        version_commit is None
        or version_commit.group(1) != expected_ramag_commit
        or version_commit.group(2) is not None
    ):
        command_status["ramag-version-commit-mismatch"] = 2
    for required in (
        "openmp enabled",
        f"max_threads {THREADS}",
        "sufkit_divsufsort_openmp enabled",
    ):
        if required not in ramag_version:
            command_status[f"ramag-version-missing:{required}"] = 2
    if "runtime unavailable" in ramag_version:
        command_status["ramag-version-openmp-runtime-unavailable"] = 2
    ramag_ldd = (provenance / "ramag.ldd.txt").read_text(
        encoding="utf-8", errors="replace"
    )
    if "libgomp" not in ramag_ldd and "libomp" not in ramag_ldd:
        command_status["ramag-ldd-openmp-runtime-missing"] = 2

    mummer_root = Path(str(config["mummer_root"]))
    files = {
        "wrapper": Path(str(config["mummer_nucmer"])),
        "elf": mummer_root / ".libs/nucmer",
        "library": mummer_root / "src/tigr/libumdmummer/.libs/libumdmummer.so.0.0.0",
    }
    # Accommodate the layout of the local 4.0.1 snapshot without guessing a
    # replacement library if the canonical location is absent.
    if not files["library"].is_file():
        candidates = list(mummer_root.rglob("libumdmummer.so.0.0.0"))
        if len(candidates) == 1:
            files["library"] = candidates[0]
    mummer_hashes: dict[str, object] = {}
    for role, path in files.items():
        mummer_hashes[role] = {
            "path": str(path.resolve()),
            "sha256": sha256_file(path) if path.is_file() else None,
            "bytes": path.stat().st_size if path.is_file() else None,
        }
        if not path.is_file():
            command_status[f"mummer4-required-{role}"] = 2
    write_json_exclusive(provenance / "mummer4-binaries.json", mummer_hashes)
    mummer_tree = fingerprint_tree(
        mummer_root, provenance / "mummer4-current-tree.sha256.tsv"
    )
    write_json_exclusive(provenance / "mummer4-current-tree.json", mummer_tree)
    source_snapshot = copy_source_snapshot(
        Path(str(config["repo"])), provenance / "ramag-source"
    )
    write_json_exclusive(provenance / "ramag-source-tree.json", source_snapshot)

    build_source = Path(str(config["ramag_source_binary"])).parent
    build_files: dict[str, object] = {}
    for filename in ("CMakeCache.txt", "compile_commands.json"):
        source = build_source / filename
        if source.is_file():
            destination = provenance / ("ramag-build-" + filename)
            shutil.copy2(source, destination)
            build_files[filename] = {
                "source": str(source.resolve()),
                "snapshot": str(destination.resolve()),
                "sha256": sha256_file(destination),
                "bytes": destination.stat().st_size,
            }
        else:
            build_files[filename] = None
    if build_files["CMakeCache.txt"] is None:
        command_status["ramag-build-CMakeCache.txt"] = 2
    write_json_exclusive(provenance / "ramag-build-files.json", build_files)
    failed_commands = {
        command: exit_code for command, exit_code in command_status.items() if exit_code != 0
    }
    summary = {
        "command_exit_codes": command_status,
        "failed_commands": failed_commands,
        "mummer4": mummer_hashes,
        "mummer4_current_tree": mummer_tree,
        "source": source_snapshot,
        "ramag_build_files": build_files,
    }
    write_json_exclusive(provenance / "summary.json", summary)
    if failed_commands:
        raise ExperimentError(f"provenance commands failed: {failed_commands}")
    accept_attempt(capture_base, name, provenance, "PROVENANCE_ACCEPTED")
    return summary


def run_preflight(root: Path, config: dict[str, object]) -> None:
    name, _attempt = next_attempt(root / "provenance/preflight")
    report_path = root / "provenance/preflight" / f"{name}.json"
    failures: list[str] = []
    package = Path(str(config["package"]))
    required = [
        package / "PACKAGE_COMPLETE",
        package / "tools/normalize_prediction.py",
        package / "tools/evaluate.py",
        Path(str(config["ramag_run_binary"])),
        Path(str(config["mummer_nucmer"])),
        Path(str(config["mummer_show_coords"])),
        Path(str(config["mummer_delta_filter"])),
        Path(str(config["reference"])),
        Path(str(config["query"])),
        Path(str(config["reference"]) + ".fai"),
        Path(str(config["query"]) + ".fai"),
    ]
    for path in required:
        if not path.is_file():
            failures.append(f"missing required file: {path}")

    package_integrity = audit_benchmark_package_integrity(package)
    failures.extend(str(failure) for failure in package_integrity["failures"])
    input_integrity = audit_fixed_input_integrity(config)
    failures.extend(str(failure) for failure in input_integrity["failures"])

    available_memory = read_mem_available()
    free_disk = shutil.disk_usage(root).free
    if available_memory < MIN_AVAILABLE_MEMORY_BYTES:
        failures.append(
            f"available memory below 20 GiB: {available_memory} bytes"
        )
    if free_disk < MIN_FREE_DISK_BYTES:
        failures.append(f"free disk below 50 GiB: {free_disk} bytes")
    conflicts = active_aligners()
    if conflicts:
        failures.append(f"other ramag/nucmer processes are active: {conflicts}")

    provenance: dict[str, object] | None = None
    try:
        provenance = collect_provenance(root, config)
    except (OSError, subprocess.SubprocessError, ExperimentError) as error:
        failures.append(f"provenance collection failed: {error}")

    report = {
        "schema": "ramag.human-chimp-preflight.v2",
        "attempt": name,
        "status": "failed" if failures else "success",
        "failures": failures,
        "package_integrity": package_integrity,
        "fixed_input_integrity": input_integrity,
        "available_memory_bytes": available_memory,
        "required_available_memory_bytes": MIN_AVAILABLE_MEMORY_BYTES,
        "free_disk_bytes": free_disk,
        "required_free_disk_bytes": MIN_FREE_DISK_BYTES,
        "active_aligners": conflicts,
        "provenance": provenance,
    }
    write_json_exclusive(report_path, report)
    if failures:
        (root / "PREFLIGHT_COMPLETE").unlink(missing_ok=True)
        (root / "PREFLIGHT_FAILED").write_text(
            "\n".join(failures) + "\n", encoding="utf-8"
        )
        raise ExperimentError("preflight failed; inspect " + str(report_path))
    (root / "PREFLIGHT_COMPLETE").write_text("success\n", encoding="ascii")


def extract_fasta_contig(source: Path, name: str, output: Path) -> None:
    if output.exists() or Path(str(output) + ".fai").exists():
        raise ExperimentError(f"refusing to overwrite smoke FASTA: {output}")
    sequence: list[str] = []
    active = False
    found = False
    with source.open("r", encoding="ascii") as handle:
        for line in handle:
            if line.startswith(">"):
                identifier = line[1:].split(None, 1)[0]
                if active:
                    break
                active = identifier == name
                found = found or active
            elif active:
                sequence.append(line.strip())
    if not found:
        raise ExperimentError(f"contig {name} not found in {source}")
    bases = "".join(sequence)
    if not bases:
        raise ExperimentError(f"contig {name} is empty in {source}")
    line_bases = 60
    temporary_output = output.with_name(f".{output.name}.{os.getpid()}.tmp")
    final_fai = Path(str(output) + ".fai")
    temporary_fai = output.with_name(f".{output.name}.{os.getpid()}.fai.tmp")
    with temporary_output.open("x", encoding="ascii", newline="\n") as handle:
        handle.write(f">{name}\n")
        for start in range(0, len(bases), line_bases):
            handle.write(bases[start : start + line_bases] + "\n")
    offset = len(name) + 2
    with temporary_fai.open("x", encoding="ascii", newline="\n") as handle:
        handle.write(f"{name}\t{len(bases)}\t{offset}\t{line_bases}\t{line_bases + 1}\n")
    os.replace(temporary_fai, final_fai)
    os.replace(temporary_output, output)


def validate_smoke_fasta(path: Path) -> None:
    contract = SMOKE_FASTA_CONTRACT.get(path.name)
    if contract is None:
        raise ExperimentError(f"no fixed smoke FASTA contract for {path.name}")
    fai = Path(str(path) + ".fai")
    if not path.is_file() or not fai.is_file():
        raise ExperimentError(f"smoke FASTA or index is incomplete: {path}")
    if sha256_file(path) != contract["sha256"]:
        raise ExperimentError(f"smoke FASTA content hash is not the fixed chrD extract: {path}")
    if fasta_record_lengths(path) != [int(contract["length"])]:
        raise ExperimentError(f"smoke FASTA length is not the fixed chrD length: {path}")
    expected_fai = (
        f"{contract['id']}\t{contract['length']}\t15\t60\t61\n"
    )
    if fai.read_text(encoding="ascii") != expected_fai:
        raise ExperimentError(f"smoke FASTA index is not canonical: {fai}")


def format_validate(
    attempt: Path,
    delta: Path,
    config: dict[str, object],
) -> None:
    if not delta.is_file() or delta.stat().st_size == 0:
        raise ExperimentError(f"raw delta is missing or empty: {delta}")
    command_evidence(
        [str(config["mummer_show_coords"]), "-rcl", str(delta)],
        attempt,
        "show-coords",
        discard_stdout=True,
    )
    command_evidence(
        [str(config["mummer_delta_filter"]), str(delta)],
        attempt,
        "delta-filter",
        discard_stdout=True,
    )
    write_json_exclusive(
        attempt / "format-validation.json",
        {
            "delta": str(delta.resolve()),
            "bytes": delta.stat().st_size,
            "sha256": sha256_file(delta),
            "show_coords_exit_code": 0,
            "delta_filter_exit_code": 0,
            "delta_filter_mode": "unfiltered parser/readability check; no -1",
            "status": "success",
        },
    )


def validate_accepted_metrics(
    metrics_path: Path,
    expected_environment: dict[str, str] | None = None,
    *,
    require_clock_consistency: bool = False,
    clock_report_path: Path | None = None,
) -> dict[str, object]:
    metrics = load_json(metrics_path)
    if metrics.get("schema") != METRICS_SCHEMA:
        raise ExperimentError(
            f"timing metrics schema is not {METRICS_SCHEMA}: {metrics.get('schema')!r}"
        )
    if metrics.get("status") != "success" or int(metrics.get("exit_code", 1)) != 0:
        raise ExperimentError("timing metrics do not report a successful command")
    if metrics.get("gnu_time_parse_error") is not None:
        raise ExperimentError(
            f"GNU time output was not parseable: {metrics.get('gnu_time_parse_error')}"
        )
    gnu_time = metrics.get("gnu_time")
    if not isinstance(gnu_time, dict):
        raise ExperimentError("timing metrics lack parsed GNU time fields")
    for field in ("elapsed_seconds", "user_seconds", "system_seconds"):
        value = gnu_time.get(field)
        if not isinstance(value, (int, float)) or value < 0:
            raise ExperimentError(f"GNU time field is invalid: {field}={value!r}")
    peak_rss = gnu_time.get("maximum_resident_set_kbytes")
    if not isinstance(peak_rss, int) or peak_rss <= 0:
        raise ExperimentError(f"GNU time maximum RSS is invalid: {peak_rss!r}")
    if gnu_time.get("exit_status") != 0:
        raise ExperimentError(
            f"GNU time exit status is not zero: {gnu_time.get('exit_status')!r}"
        )
    for field in (
        "sample_count",
        "observed_max_total_threads",
        "observed_max_single_process_threads",
        "observed_max_aligner_tree_threads",
        "observed_max_aligner_single_process_threads",
        "observed_process_tree_peak_rss_bytes",
    ):
        value = metrics.get(field)
        if not isinstance(value, int) or value <= 0:
            raise ExperimentError(f"process-tree metric is invalid: {field}={value!r}")
    total_threads = int(metrics["observed_max_total_threads"])
    total_single = int(metrics["observed_max_single_process_threads"])
    aligner_tree = int(metrics["observed_max_aligner_tree_threads"])
    aligner_single = int(metrics["observed_max_aligner_single_process_threads"])
    if not (
        aligner_single <= aligner_tree <= total_threads
        and aligner_single <= total_single <= total_threads
    ):
        raise ExperimentError(
            "process-tree thread metrics violate aligner/single-process/total hierarchy"
        )
    if expected_environment is not None and metrics.get(
        "environment_overrides"
    ) != expected_environment:
        raise ExperimentError("timing metrics do not record the accepted environment")
    timeout = metrics.get("timeout")
    if (
        not isinstance(timeout, dict)
        or timeout.get("clock") != "CLOCK_MONOTONIC"
        or timeout.get("exceeded") is not False
        or timeout.get("elapsed_at_signal_seconds") is not None
        or timeout.get("sigterm_sent") is not False
        or timeout.get("sigkill_sent") is not False
        or timeout.get("deadline_overshoot_seconds") is not None
    ):
        raise ExperimentError("successful timing metrics contain invalid timeout evidence")
    timeout_enabled = timeout.get("enabled")
    timeout_limit = timeout.get("limit_seconds")
    if (
        not isinstance(timeout_enabled, bool)
        or not isinstance(timeout_limit, (int, float))
        or not math.isfinite(float(timeout_limit))
        or (timeout_enabled and float(timeout_limit) <= 0.0)
        or (not timeout_enabled and float(timeout_limit) != 0.0)
        or (metrics_path.parent / "RUN_TIMED_OUT").exists()
    ):
        raise ExperimentError("successful run has contradictory timeout state")
    if require_clock_consistency:
        report = assess_clock_consistency(
            metrics,
            metrics_path.parent / "resources.clock.tsv",
            tolerance_fraction=CLOCK_CONSISTENCY_TOLERANCE_FRACTION,
        )
        if clock_report_path is not None:
            write_json_exclusive(clock_report_path, report)
        if report.get("status") == "quarantined":
            raise ClockQuarantineError(
                "elapsed-clock evidence is non-finite/non-positive, internally "
                "inconsistent, or contains a backwards CLOCK_REALTIME sample"
            )
    return metrics


def canonicalize_with_timing(
    source: Path,
    destination: Path,
    timing_path: Path,
) -> dict[str, object]:
    started_monotonic_ns = time.monotonic_ns()
    started_realtime_ns = time.time_ns()
    status = "failed"
    try:
        result = canonicalize(source, destination)
        status = "success"
        return result
    finally:
        finished_monotonic_ns = time.monotonic_ns()
        finished_realtime_ns = time.time_ns()
        write_json_exclusive(
            timing_path,
            {
                "schema": "ramag.benchmark-postprocess-timing.v1",
                "name": "canonicalize-sources",
                "clock": "CLOCK_MONOTONIC",
                "started_monotonic_ns": started_monotonic_ns,
                "finished_monotonic_ns": finished_monotonic_ns,
                "monotonic_elapsed_seconds": (
                    finished_monotonic_ns - started_monotonic_ns
                )
                / 1_000_000_000.0,
                "started_utc": _utc_from_realtime_ns(started_realtime_ns),
                "finished_utc": _utc_from_realtime_ns(finished_realtime_ns),
                "status": status,
            },
        )


def validate_ramag_manifest(
    manifest: dict[str, object],
    config: dict[str, object],
    *,
    expected_task_shape: dict[str, int],
    expected_reference: Path,
    expected_query: Path,
    verify_fixed_input_routes: bool = False,
) -> None:
    dependencies = manifest.get("dependencies")
    if manifest.get("status") != "success" or int(manifest.get("exit_code", 1)) != 0:
        raise ExperimentError("RaMA-G manifest does not report success")
    if (
        not isinstance(dependencies, dict)
        or dependencies.get("sufkit_commit") != config["expected_sufkit_commit"]
        or dependencies.get("seqpro_commit") != config["expected_seqpro_commit"]
    ):
        raise ExperimentError("RaMA-G manifest records the wrong dependency commit")
    effective = manifest.get("effective_config")
    if (
        not isinstance(effective, dict)
        or effective.get("reference") != str(expected_reference.resolve())
        or effective.get("query") != str(expected_query.resolve())
        or effective.get("formats") != "delta"
        or effective.get("seed_mode") != "mumreference"
        or effective.get("selection_mode") != config["ramag_selection_mode"]
        or int(effective.get("threads", 0)) != int(config["threads"])
        or int(effective.get("min_match", -1)) != int(config["min_match"])
        or int(effective.get("max_gap", -1)) != int(config["max_gap"])
        or int(effective.get("diag_diff", -1)) != int(config["diag_diff"])
        or float(effective.get("diag_factor", -1.0)) != float(config["diag_factor"])
        or int(effective.get("min_cluster", -1)) != int(config["min_cluster"])
        or int(effective.get("break_length", -1)) != int(config["break_length"])
        or int(effective.get("max_dp_cells", -1)) != int(config["max_dp_cells"])
    ):
        raise ExperimentError("RaMA-G manifest records the wrong effective alignment config")
    actual_routes = manifest.get("actual_routes")
    chaining = manifest.get("chaining")
    expected_chain_route = "sparse-exact-edge-components-v1"
    if (
        not isinstance(actual_routes, dict)
        or actual_routes.get("seed") != EXPECTED_SEED_ROUTE
        or actual_routes.get("chain") != expected_chain_route
        or not isinstance(chaining, dict)
        or chaining.get("route") != expected_chain_route
    ):
        raise ExperimentError(
            "RaMA-G manifest does not record sparse-exact-edge-components-v1"
        )
    adapter = manifest.get("adapter_provenance")
    expected_backend = (
        "caps32"
        if expected_task_shape == FULL_MAM_TASK_SHAPE
        else "divsufsort32"
    )
    if (
        not isinstance(adapter, dict)
        or adapter.get("sufkit.commit") != config["expected_sufkit_commit"]
        or adapter.get("sufkit.source_state") != "exact-commit-clean-at-configure"
        or adapter.get("sufkit.backend") != expected_backend
        or int(adapter.get("sufkit.build.caps_min_reference_bases", 0))
        != CAPS_MIN_REFERENCE_BASES
    ):
        raise ExperimentError(
            "RaMA-G adapter provenance does not record the accepted clean sufkit "
            "source and fixed SA-backend policy"
        )
    openmp = manifest.get("openmp")
    build = manifest.get("build")
    threading = manifest.get("threading")
    mam_tiling = manifest.get("mam_tiling")
    threads = int(config.get("threads", 0))
    seed_scheduled = (
        int(threading.get("seed_scheduled_threads", 0))
        if isinstance(threading, dict)
        else 0
    )
    seed_workers = (
        int(threading.get("seed_worker_threads", 0))
        if isinstance(threading, dict)
        else 0
    )
    seed_tasks = (
        int(threading.get("seed_task_count", 0))
        if isinstance(threading, dict)
        else 0
    )
    seed_completed = (
        int(threading.get("seed_tasks_completed", 0))
        if isinstance(threading, dict)
        else 0
    )
    tile_tasks = (
        int(mam_tiling.get("tile_tasks", 0))
        if isinstance(mam_tiling, dict)
        else 0
    )
    tile_completed = (
        int(mam_tiling.get("tile_tasks_completed", 0))
        if isinstance(mam_tiling, dict)
        else 0
    )
    boundary_tasks = (
        int(mam_tiling.get("boundary_tasks", 0))
        if isinstance(mam_tiling, dict)
        else 0
    )
    boundary_completed = (
        int(mam_tiling.get("boundary_tasks_completed", 0))
        if isinstance(mam_tiling, dict)
        else 0
    )
    mam_limits = (
        mam_tiling.get("resource_limits", {})
        if isinstance(mam_tiling, dict)
        else {}
    )
    expected_oriented = int(expected_task_shape["oriented_queries"])
    expected_tile_tasks = int(expected_task_shape["tile_tasks"])
    expected_boundary_tasks = int(expected_task_shape["boundary_tasks"])
    expected_total_tasks = int(expected_task_shape["total_tasks"])
    if expected_total_tasks != expected_tile_tasks + expected_boundary_tasks:
        raise ExperimentError("internal expected tiled-MAM task shape is inconsistent")
    expected_seed_workers = min(expected_total_tasks, threads, 16)
    if (
        not isinstance(openmp, dict)
        or openmp.get("enabled") is not True
        or openmp.get("runtime") in {None, "", "unavailable"}
        or int(openmp.get("runtime_max_threads", 0)) != threads
        or int(openmp.get("run_requested_threads", 0)) != threads
        or int(openmp.get("configured_requested_threads", 0)) != threads
        or openmp.get("dynamic") is not False
        or openmp.get("sufkit_divsufsort_openmp") is not True
        or not isinstance(build, dict)
        or build.get("openmp_enabled") is not True
        or build.get("sufkit_divsufsort_openmp") is not True
        or not isinstance(threading, dict)
        or int(threading.get("openmp_actual_requested_threads", 0)) != threads
        or int(threading.get("openmp_runtime_max_threads", 0)) != threads
        or int(threading.get("input_requested_workers", 0)) != 2
        or int(threading.get("input_actual_workers", 0)) != 2
        or threading.get("input_parallel_route")
        != "openmp-sections-reference-query-v1"
        or actual_routes.get("input_parallel")
        != "openmp-sections-reference-query-v1"
        or int(threading.get("seed_requested_threads", 0)) != threads
        or not 1 < seed_scheduled <= threads
        or seed_workers != seed_scheduled
        or seed_scheduled != expected_seed_workers
        or seed_tasks < seed_scheduled
        or seed_completed != seed_tasks
        or threading.get("seed_parallel_route")
        != "openmp-dynamic-stable-tile-boundary-tasks"
        or not isinstance(mam_tiling, dict)
        or int(mam_tiling.get("worker_cap", 0)) != 16
        or int(mam_tiling.get("tile_bases", 0)) != MAM_TILE_BASES
        or int(mam_tiling.get("oriented_queries", 0)) != expected_oriented
        or tile_tasks != expected_tile_tasks
        or boundary_tasks != expected_boundary_tasks
        or tile_completed != tile_tasks
        or boundary_completed != boundary_tasks
        or seed_tasks != expected_total_tasks
        or seed_tasks != tile_tasks + boundary_tasks
        or not isinstance(mam_limits, dict)
        or int(mam_tiling.get("workspace_peak_bytes", 0))
        > int(mam_limits.get("workspace_bytes", -1))
        or int(threading.get("chaining_requested_threads", 0)) != threads
        or int(threading.get("chaining_worker_threads", 0)) <= 1
        or int(threading.get("chaining_worker_threads", 0)) > threads
        or int(threading.get("extension_worker_threads", 0)) <= 1
        or int(threading.get("extension_worker_threads", 0)) > threads
        or threading.get("seed_chain_extension_parallel") is not True
    ):
        raise ExperimentError(
            "RaMA-G manifest does not record the accepted OpenMP runtime/build contract"
        )

    if verify_fixed_input_routes:
        fixed_input_integrity: dict[str, object] | None = None
        if expected_task_shape == FULL_MAM_TASK_SHAPE:
            fixed_input_integrity = require_fixed_input_integrity(
                config,
                reference=expected_reference,
                query=expected_query,
            )
        if actual_routes.get("input") != (
            "seqpro-external-fai-mmap+caller-buffer-copy"
        ):
            raise ExperimentError(
                "RaMA-G full/smoke run did not use the accepted external-FAI input route"
            )
        for role, fasta in (
            ("reference", expected_reference),
            ("query", expected_query),
        ):
            source_fai = Path(str(fasta) + ".fai")
            if not source_fai.is_file():
                raise ExperimentError(f"accepted {role} FAI is missing: {source_fai}")
            fai_filename = source_fai.name
            if fixed_input_integrity is not None:
                expected_fai_sha256 = EXPECTED_FAI_SHA256.get(fai_filename)
                if not isinstance(expected_fai_sha256, str):
                    raise ExperimentError(
                        f"fixed post-run input has an unexpected FAI name: {fai_filename}"
                    )
                observation = fixed_input_integrity["observations"].get(fai_filename)
                if not isinstance(observation, dict):
                    raise ExperimentError(
                        f"fixed post-run input evidence lacks {fai_filename}"
                    )
                expected_fai_bytes = int(observation["bytes"])
            else:
                # The chrD smoke inputs are deterministic derived artifacts,
                # validated separately by validate_smoke_fasta().
                expected_fai_sha256 = sha256_file(source_fai)
                expected_fai_bytes = source_fai.stat().st_size
            prefix = f"seqpro.{role}."
            active_fai = Path(str(adapter.get(prefix + "fai", "")))
            if (
                adapter.get(prefix + "source_fai") != str(source_fai.resolve())
                or int(adapter.get(prefix + "source_fai_bytes", -1))
                != expected_fai_bytes
                or adapter.get(prefix + "source_fai_sha256")
                != expected_fai_sha256
                or adapter.get(prefix + "source_fai.copy_status") != "copied"
                or adapter.get(prefix + "external_fai.adoption_status")
                != "adopted"
                or adapter.get(prefix + "build_action") != "reused"
                or adapter.get(prefix + "index_origin")
                != "external-standard-fai"
                or adapter.get(prefix + "verification") != "structure-validated"
                or adapter.get(prefix + "metadata") != ""
                or not active_fai.is_file()
                or active_fai.stat().st_size != expected_fai_bytes
                or sha256_file(active_fai) != expected_fai_sha256
            ):
                raise ExperimentError(
                    f"RaMA-G {role} external-FAI provenance/active copy is invalid"
                )

    counts = manifest.get("counts")
    if not isinstance(counts, dict):
        raise ExperimentError("RaMA-G manifest lacks tiled-MAM count cross-checks")
    statistic_names = (
        "short_query_tasks",
        "tile_raw_mams",
        "tile_globally_maximal_mams",
        "boundary_raw_mems",
        "boundary_patterns",
        "boundary_reference_unique_patterns",
        "boundary_recovered_mams",
        "workspace_baseline_bytes",
        "workspace_peak_bytes",
    )
    statistics = {name: int(mam_tiling.get(name, -1)) for name in statistic_names}
    if any(value < 0 for value in statistics.values()):
        raise ExperimentError("RaMA-G manifest has missing or negative tiled-MAM statistics")
    boundary_limit = int(mam_limits.get("boundary_mem_occurrences", -1))
    workspace_limit = int(mam_limits.get("workspace_bytes", -1))
    if boundary_limit <= 0 or workspace_limit <= 0:
        raise ExperimentError("RaMA-G manifest has invalid tiled-MAM resource limits")
    if not (
        statistics["short_query_tasks"] <= expected_oriented
        and statistics["short_query_tasks"] <= tile_tasks
        and statistics["tile_globally_maximal_mams"]
        <= statistics["tile_raw_mams"]
        and statistics["boundary_patterns"] <= statistics["boundary_raw_mems"]
        and statistics["boundary_reference_unique_patterns"]
        <= statistics["boundary_patterns"]
        and statistics["boundary_recovered_mams"]
        <= statistics["boundary_reference_unique_patterns"]
        and statistics["boundary_raw_mems"] <= boundary_limit
        and statistics["workspace_baseline_bytes"]
        <= statistics["workspace_peak_bytes"]
        <= workspace_limit
    ):
        raise ExperimentError("RaMA-G manifest violates tiled-MAM statistical invariants")
    mem_seeds = int(counts.get("mem_seeds", -1))
    mam_seeds = int(counts.get("mam_seeds", -1))
    selected_seeds = int(counts.get("selected_seeds", -1))
    if (
        mem_seeds != statistics["boundary_raw_mems"]
        or mam_seeds
        != statistics["tile_globally_maximal_mams"]
        + statistics["boundary_recovered_mams"]
        or not 0 <= selected_seeds <= mam_seeds
    ):
        raise ExperimentError("RaMA-G manifest tiled-MAM counts disagree with top-level counts")
    if expected_task_shape == FULL_MAM_TASK_SHAPE and (
        seed_scheduled != threads or seed_workers != threads
    ):
        raise ExperimentError(
            "full Human-Chimp RaMA-G run did not use the complete 16-worker "
            "tiled-MAM budget"
        )


def alignment_command(
    tool: str,
    attempt: Path,
    reference: Path,
    query: Path,
    config: dict[str, object],
) -> tuple[list[str], Path]:
    if tool == "mummer4":
        delta = attempt / "mummer4.mumreference.raw.delta"
        command = [
            str(config["mummer_nucmer"]),
            f"--threads={config['threads']}",
            f"--minmatch={config['min_match']}",
            f"--mincluster={config['min_cluster']}",
            f"--maxgap={config['max_gap']}",
            f"--diagdiff={config['diag_diff']}",
            f"--diagfactor={config['diag_factor']}",
            f"--breaklen={config['break_length']}",
            f"--delta={delta}",
            str(reference),
            str(query),
        ]
        return command, delta
    prefix = attempt / "ramag.mumreference.raw"
    delta = Path(str(prefix) + ".delta")
    command = [
        str(config["ramag_run_binary"]),
        "align",
        "--reference",
        str(reference),
        "--query",
        str(query),
        "--output-prefix",
        str(prefix),
        "--work-dir",
        str(attempt / "work"),
        "--threads",
        str(config["threads"]),
        "--formats",
        "delta",
        "--seed-mode",
        "mumreference",
        "--selection-mode",
        str(config["ramag_selection_mode"]),
        "--min-match",
        str(config["min_match"]),
        "--max-gap",
        str(config["max_gap"]),
        "--diag-diff",
        str(config["diag_diff"]),
        "--diag-factor",
        str(config["diag_factor"]),
        "--min-cluster",
        str(config["min_cluster"]),
        "--break-length",
        str(config["break_length"]),
        "--max-dp-cells",
        str(config["max_dp_cells"]),
    ]
    return command, delta


def run_alignment_attempt(
    base: Path,
    tool: str,
    reference: Path,
    query: Path,
    config: dict[str, object],
    *,
    expected_task_shape: dict[str, int] | None = None,
    full_acceptance: bool = False,
    timeout_seconds: float | None = None,
) -> Path:
    accepted = selected_attempt(base, "RUN_ACCEPTED")
    if accepted is not None:
        return accepted
    name, attempt = next_attempt(base)
    command, delta = alignment_command(tool, attempt, reference, query, config)
    write_json_exclusive(attempt / "planned-command.json", command)
    warm_cache: dict[str, object] | None = None
    if full_acceptance:
        warm_cache = warm_cache_inputs_for_tool(
            tool, attempt, reference, query, config
        )
    exit_code = run_with_metrics(
        attempt / "timing",
        command,
        float(config["sample_interval_seconds"]),
        dict(config["openmp_environment"]),
        timeout_seconds,
    )
    raw_metrics = load_json(attempt / "timing/metrics.json")
    if full_acceptance:
        assert warm_cache is not None
        record_warm_cache_launch_relation(
            attempt,
            warm_cache,
            raw_metrics,
        )
    if raw_metrics.get("status") == "timed_out":
        if (
            exit_code == 0
            or raw_metrics.get("schema") != METRICS_SCHEMA
            or not (attempt / "timing/RUN_TIMED_OUT").is_file()
            or (attempt / "timing/RUN_COMPLETE").exists()
            or (attempt / "RUN_ACCEPTED").exists()
        ):
            raise ExperimentError(f"{tool} timeout evidence is internally inconsistent")
        raise AlignmentTimeoutError(tool, attempt, raw_metrics)
    if exit_code != 0:
        raise ExperimentError(f"{tool} alignment failed; preserved {attempt}")
    try:
        metrics = validate_accepted_metrics(
            attempt / "timing/metrics.json",
            dict(config["openmp_environment"]),
            require_clock_consistency=full_acceptance,
            clock_report_path=(
                attempt / "clock-consistency.json" if full_acceptance else None
            ),
        )
    except ClockQuarantineError as error:
        (attempt / "CLOCK_QUARANTINED").write_text(
            str(error) + "\n", encoding="utf-8"
        )
        raise ExperimentError(
            f"{tool} full alignment was quarantined and cannot be selected; "
            f"raw evidence is preserved at {attempt}"
        ) from error
    if full_acceptance:
        clock_report = load_json(attempt / "clock-consistency.json")
        if clock_report.get("status") == "warning":
            (attempt / "CLOCK_WARNING").write_text(
                "1% drift or positive realtime jump retained for cross-tool review\n",
                encoding="ascii",
            )
        (attempt / "CLOCK_ACCEPTED").write_text(
            str(clock_report.get("status")) + "\n", encoding="ascii"
        )
    if tool == "ramag" and int(
        metrics.get("observed_max_aligner_single_process_threads", 0)
    ) <= 1:
        raise ExperimentError(
            "RaMA-G process sampling did not independently observe parallel threads"
        )
    format_validate(attempt, delta, config)
    if tool == "ramag":
        if expected_task_shape is None:
            raise ExperimentError("RaMA-G alignment lacks an expected tiled-MAM task shape")
        prefix = attempt / "ramag.mumreference.raw"
        manifest_path = Path(str(prefix) + ".manifest.json")
        complete_path = Path(str(prefix) + ".complete")
        if not manifest_path.is_file() or not complete_path.is_file():
            raise ExperimentError("RaMA-G did not publish manifest and complete marker")
        manifest = load_json(manifest_path)
        validate_ramag_manifest(
            manifest,
            config,
            expected_task_shape=expected_task_shape,
            expected_reference=reference,
            expected_query=query,
            verify_fixed_input_routes=True,
        )
    accept_attempt(base, name, attempt, "RUN_ACCEPTED")
    return attempt


def normalize_and_validate_smoke(
    tool: str,
    alignment_attempt: Path,
    validation_attempt: Path,
    reference: Path,
    query: Path,
    config: dict[str, object],
) -> None:
    package = Path(str(config["package"]))
    delta = alignment_attempt / (
        "mummer4.mumreference.raw.delta" if tool == "mummer4" else "ramag.mumreference.raw.delta"
    )
    converted = validation_attempt / "prediction.converter.maf"
    converter_stats = validation_attempt / "prediction.converter.stats.json"
    command_evidence(
        [
            sys.executable,
            str(package / "tools/normalize_prediction.py"),
            "--prediction",
            str(delta),
            "--format",
            "delta",
            "--reference",
            str(reference),
            "--query",
            str(query),
            "--output",
            str(converted),
            "--stats",
            str(converter_stats),
        ],
        validation_attempt,
        "normalize-prediction",
    )
    canonical = validation_attempt / "prediction.canonical.maf"
    canonical_stats = canonicalize_with_timing(
        converted,
        canonical,
        validation_attempt / "canonicalize-sources.timing.json",
    )
    write_json_exclusive(validation_attempt / "canonicalization.json", canonical_stats)
    command_evidence(
        [
            sys.executable,
            str(package / "tools/validate_maf.py"),
            "--maf",
            str(canonical),
            "--reference",
            str(reference),
            "--query",
            str(query),
            "--output-json",
            str(validation_attempt / "prediction.validation.json"),
            "--allow-singletons",
        ],
        validation_attempt,
        "validate-maf",
    )


def run_smoke(root: Path, config: dict[str, object]) -> None:
    if (root / "SMOKE_COMPLETE").is_file():
        return
    if not (root / "PREFLIGHT_COMPLETE").is_file():
        raise ExperimentError("smoke requires a successful preflight")
    inputs = root / "smoke/inputs"
    inputs.mkdir(parents=True, exist_ok=True)
    reference = inputs / "simHuman.chrD.fa"
    query = inputs / "simChimp.chrD.fa"
    if not reference.exists():
        extract_fasta_contig(Path(str(config["reference"])), "simHuman.chrD", reference)
    if not query.exists():
        extract_fasta_contig(Path(str(config["query"])), "simChimp.chrD", query)
    validate_smoke_fasta(reference)
    validate_smoke_fasta(query)
    smoke_task_shape = expected_mam_task_shape(query)
    if smoke_task_shape != SMOKE_MAM_TASK_SHAPE:
        raise ExperimentError(
            "fixed chrD smoke query does not yield the accepted 6 tile + 4 "
            "boundary task shape"
        )
    for tool in ("mummer4", "ramag"):
        base = root / "smoke" / tool
        attempt = run_alignment_attempt(
            base,
            tool,
            reference,
            query,
            config,
            expected_task_shape=(
                smoke_task_shape if tool == "ramag" else None
            ),
        )
        validation_base = base / "maf-validation"
        if selected_attempt(validation_base, "SMOKE_MAF_VALIDATED") is None:
            validation_name, validation_attempt = next_attempt(validation_base)
            normalize_and_validate_smoke(
                tool, attempt, validation_attempt, reference, query, config
            )
            accept_attempt(
                validation_base,
                validation_name,
                validation_attempt,
                "SMOKE_MAF_VALIDATED",
            )
    (root / "SMOKE_COMPLETE").write_text("success\n", encoding="ascii")


def _accepted_alignment_clocks(
    attempt: Path,
    config: dict[str, object],
) -> tuple[dict[str, object], float, float]:
    metrics = validate_accepted_metrics(
        attempt / "timing/metrics.json",
        dict(config["openmp_environment"]),
        require_clock_consistency=True,
    )
    gnu_time = metrics.get("gnu_time")
    measurements = metrics.get("clock_measurements")
    if not isinstance(gnu_time, dict) or not isinstance(measurements, dict):
        raise ExperimentError(f"accepted attempt lacks elapsed clocks: {attempt}")
    gnu_elapsed = gnu_time.get("elapsed_seconds")
    monotonic_elapsed = measurements.get("monotonic_elapsed_seconds")
    if (
        not isinstance(gnu_elapsed, (int, float))
        or not math.isfinite(float(gnu_elapsed))
        or float(gnu_elapsed) <= 0.0
        or not isinstance(monotonic_elapsed, (int, float))
        or not math.isfinite(float(monotonic_elapsed))
        or float(monotonic_elapsed) <= 0.0
    ):
        raise ExperimentError(f"accepted attempt has invalid elapsed clocks: {attempt}")
    return metrics, float(gnu_elapsed), float(monotonic_elapsed)


def _write_marker_exclusive(path: Path, value: str = "success\n") -> None:
    with path.open("x", encoding="ascii", newline="\n") as handle:
        handle.write(value)


def freeze_mummer_speed_budget(
    root: Path,
    mummer_attempt: Path,
    config: dict[str, object],
) -> dict[str, object]:
    metrics, gnu_elapsed, monotonic_elapsed = _accepted_alignment_clocks(
        mummer_attempt, config
    )
    metrics_path = mummer_attempt / "timing/metrics.json"
    report = {
        "schema": "ramag.human-chimp-speed-budget.v1",
        "status": "accepted",
        "result_root": str(root.resolve()),
        "source_tool": "mummer4",
        "source_attempt": str(mummer_attempt.relative_to(root)),
        "source_metrics": str(metrics_path.relative_to(root)),
        "source_metrics_sha256": sha256_file(metrics_path),
        "source_metrics_schema": metrics.get("schema"),
        "source_gnu_elapsed_seconds": gnu_elapsed,
        "cutoff_clock": "CLOCK_MONOTONIC",
        "cutoff_seconds": monotonic_elapsed,
        "timeout_condition": "ramag_monotonic_elapsed_seconds > cutoff_seconds",
        "equal_is_accepted": True,
        "observation_policy": "same-root single accepted MUMmer4 observation",
    }
    path = root / "speed-budget.json"
    marker = root / "MUMMER_BASELINE_ACCEPTED"
    if marker.exists() and not path.is_file():
        raise ExperimentError("MUMMER_BASELINE_ACCEPTED exists without speed-budget.json")
    if path.exists():
        if load_json(path) != report:
            raise ExperimentError("frozen MUMmer4 speed budget no longer matches its metrics")
    else:
        write_json_exclusive(path, report)
    if not marker.exists():
        _write_marker_exclusive(marker)
    return report


def _publish_alignment_speed_gate(root: Path, report: dict[str, object]) -> None:
    status = report.get("status")
    if status not in {"accepted", "rejected"}:
        raise ExperimentError(f"invalid alignment speed-gate status: {status!r}")
    accepted_marker = root / "ALIGNMENT_SPEED_ACCEPTED"
    rejected_marker = root / "ALIGNMENT_SPEED_REJECTED"
    marker = accepted_marker if status == "accepted" else rejected_marker
    opposite = rejected_marker if status == "accepted" else accepted_marker
    if opposite.exists():
        raise ExperimentError("alignment speed gate has contradictory completion markers")
    path = root / "alignment-speed-gate.json"
    if path.exists():
        if load_json(path) != report:
            raise ExperimentError("alignment speed-gate evidence changed within one root")
    else:
        write_json_exclusive(path, report)
    if not marker.exists():
        _write_marker_exclusive(marker, status + "\n")


def build_completed_alignment_speed_gate(
    root: Path,
    mummer_attempt: Path,
    ramag_attempt: Path,
    budget: dict[str, object],
    config: dict[str, object],
) -> dict[str, object]:
    mummer_metrics, mummer_gnu, mummer_monotonic = _accepted_alignment_clocks(
        mummer_attempt, config
    )
    ramag_metrics, ramag_gnu, ramag_monotonic = _accepted_alignment_clocks(
        ramag_attempt, config
    )
    cutoff = float(budget["cutoff_seconds"])
    if (
        budget.get("result_root") != str(root.resolve())
        or budget.get("source_attempt") != str(mummer_attempt.relative_to(root))
        or cutoff != mummer_monotonic
    ):
        raise ExperimentError("speed budget is not tied to this root's accepted MUMmer4 run")
    conditions = {
        "ramag_completed_before_hard_deadline": ramag_monotonic <= cutoff,
        "ramag_gnu_not_slower": ramag_gnu <= mummer_gnu,
        "ramag_monotonic_not_slower": ramag_monotonic <= mummer_monotonic,
    }
    accepted = all(conditions.values())
    return {
        "schema": "ramag.human-chimp-alignment-speed-gate.v1",
        "status": "accepted" if accepted else "rejected",
        "reason": "dual_clock_passed" if accepted else "dual_clock_failed",
        "scope": "preliminary single observation",
        "endpoint": config["runtime_endpoint"],
        "equal_is_accepted": True,
        "hard_deadline_clock": "CLOCK_MONOTONIC",
        "hard_deadline_seconds": cutoff,
        "budget_sha256": sha256_file(root / "speed-budget.json"),
        "conditions": conditions,
        "mummer4": {
            "attempt": str(mummer_attempt.relative_to(root)),
            "metrics_sha256": sha256_file(mummer_attempt / "timing/metrics.json"),
            "metrics_schema": mummer_metrics.get("schema"),
            "gnu_elapsed_seconds": mummer_gnu,
            "monotonic_elapsed_seconds": mummer_monotonic,
        },
        "ramag": {
            "attempt": str(ramag_attempt.relative_to(root)),
            "metrics_sha256": sha256_file(ramag_attempt / "timing/metrics.json"),
            "metrics_schema": ramag_metrics.get("schema"),
            "status": "success",
            "gnu_elapsed_seconds": ramag_gnu,
            "monotonic_elapsed_seconds": ramag_monotonic,
        },
    }


def build_timeout_alignment_speed_gate(
    root: Path,
    mummer_attempt: Path,
    timeout_error: AlignmentTimeoutError,
    budget: dict[str, object],
    config: dict[str, object],
) -> dict[str, object]:
    mummer_metrics, mummer_gnu, mummer_monotonic = _accepted_alignment_clocks(
        mummer_attempt, config
    )
    cutoff = float(budget["cutoff_seconds"])
    metrics = timeout_error.metrics
    timeout = metrics.get("timeout")
    measurements = metrics.get("clock_measurements")
    if not isinstance(timeout, dict) or not isinstance(measurements, dict):
        raise ExperimentError("RaMA-G timeout lacks structured clock evidence")
    elapsed_at_signal = timeout.get("elapsed_at_signal_seconds")
    if (
        budget.get("result_root") != str(root.resolve())
        or budget.get("source_attempt") != str(mummer_attempt.relative_to(root))
        or cutoff != mummer_monotonic
        or metrics.get("schema") != METRICS_SCHEMA
        or metrics.get("status") != "timed_out"
        or timeout.get("enabled") is not True
        or timeout.get("exceeded") is not True
        or timeout.get("clock") != "CLOCK_MONOTONIC"
        or float(timeout.get("limit_seconds", -1.0)) != cutoff
        or not isinstance(elapsed_at_signal, (int, float))
        or float(elapsed_at_signal) <= cutoff
        or timeout.get("sigterm_sent") is not True
        or timeout.get("process_group_alive_after_termination") is not False
        or not (timeout_error.attempt / "timing/RUN_TIMED_OUT").is_file()
        or (timeout_error.attempt / "timing/RUN_COMPLETE").exists()
        or (timeout_error.attempt / "RUN_ACCEPTED").exists()
    ):
        raise ExperimentError("RaMA-G timeout evidence violates the hard-deadline contract")
    return {
        "schema": "ramag.human-chimp-alignment-speed-gate.v1",
        "status": "rejected",
        "reason": "ramag_timed_out",
        "scope": "preliminary single observation",
        "endpoint": config["runtime_endpoint"],
        "equal_is_accepted": True,
        "hard_deadline_clock": "CLOCK_MONOTONIC",
        "hard_deadline_seconds": cutoff,
        "budget_sha256": sha256_file(root / "speed-budget.json"),
        "conditions": {
            "ramag_completed_before_hard_deadline": False,
            "ramag_gnu_not_slower": None,
            "ramag_monotonic_not_slower": False,
        },
        "mummer4": {
            "attempt": str(mummer_attempt.relative_to(root)),
            "metrics_sha256": sha256_file(mummer_attempt / "timing/metrics.json"),
            "metrics_schema": mummer_metrics.get("schema"),
            "gnu_elapsed_seconds": mummer_gnu,
            "monotonic_elapsed_seconds": mummer_monotonic,
        },
        "ramag": {
            "attempt": str(timeout_error.attempt.relative_to(root)),
            "metrics_sha256": sha256_file(
                timeout_error.attempt / "timing/metrics.json"
            ),
            "metrics_schema": metrics.get("schema"),
            "status": "timed_out",
            "gnu_elapsed_seconds": (
                metrics.get("gnu_time", {}).get("elapsed_seconds")
                if isinstance(metrics.get("gnu_time"), dict)
                else None
            ),
            "monotonic_elapsed_seconds": measurements.get(
                "monotonic_elapsed_seconds"
            ),
            "timeout": timeout,
        },
    }


def run_full(root: Path, config: dict[str, object]) -> None:
    if (root / "ALIGNMENT_SPEED_REJECTED").is_file():
        raise ExperimentError(
            "this experiment root already failed the alignment speed gate; "
            "use a new exclusive root"
        )
    if (root / "FULL_RUNS_COMPLETE").is_file():
        if not (root / "ALIGNMENT_SPEED_ACCEPTED").is_file():
            raise ExperimentError("FULL_RUNS_COMPLETE lacks ALIGNMENT_SPEED_ACCEPTED")
        return
    if not (root / "SMOKE_COMPLETE").is_file():
        raise ExperimentError("full runs require the chrD smoke gate")
    reference = Path(str(config["reference"]))
    query = Path(str(config["query"]))
    require_fixed_input_integrity(config, reference=reference, query=query)
    full_task_shape = expected_mam_task_shape(query)
    if full_task_shape != FULL_MAM_TASK_SHAPE:
        raise ExperimentError(
            "fixed Human-Chimp query does not yield the accepted 92 tile + 84 "
            "boundary task shape"
        )
    # Fixed fairness order: MUMmer4 first, then RaMA-G.
    mummer_base = root / "runs/mummer4"
    mummer_attempt = selected_attempt(mummer_base, "RUN_ACCEPTED")
    if mummer_attempt is None:
        record_full_launch_gate(root, "mummer4", config)
        mummer_attempt = run_alignment_attempt(
            mummer_base,
            "mummer4",
            reference,
            query,
            config,
            full_acceptance=True,
        )
    budget = freeze_mummer_speed_budget(root, mummer_attempt, config)
    ramag_base = root / "runs/ramag"
    ramag_attempt = selected_attempt(ramag_base, "RUN_ACCEPTED")
    if ramag_attempt is None:
        record_full_launch_gate(root, "ramag", config)
        try:
            ramag_attempt = run_alignment_attempt(
                ramag_base,
                "ramag",
                reference,
                query,
                config,
                expected_task_shape=full_task_shape,
                full_acceptance=True,
                timeout_seconds=float(budget["cutoff_seconds"]),
            )
        except AlignmentTimeoutError as error:
            report = build_timeout_alignment_speed_gate(
                root, mummer_attempt, error, budget, config
            )
            _publish_alignment_speed_gate(root, report)
            raise ExperimentError(
                "RaMA-G exceeded this root's accepted MUMmer4 monotonic time; "
                "MAF conversion and F1 were not started"
            ) from error
    report = build_completed_alignment_speed_gate(
        root, mummer_attempt, ramag_attempt, budget, config
    )
    _publish_alignment_speed_gate(root, report)
    if report["status"] != "accepted":
        raise ExperimentError(
            "RaMA-G completed but failed the GNU/monotonic alignment speed gate; "
            "MAF conversion and F1 were not started"
        )
    if not (root / "FULL_RUNS_COMPLETE").exists():
        _write_marker_exclusive(root / "FULL_RUNS_COMPLETE")


def run_evaluation_attempt(
    root: Path, tool: str, config: dict[str, object]
) -> Path:
    base = root / "evaluation" / tool
    accepted = selected_attempt(base, "EVALUATION_ACCEPTED")
    if accepted is not None:
        return accepted
    run_attempt = selected_attempt(root / "runs" / tool, "RUN_ACCEPTED")
    if run_attempt is None:
        raise ExperimentError(f"cannot evaluate {tool}: no accepted full run")
    name, attempt = next_attempt(base)
    package = Path(str(config["package"]))
    delta = run_attempt / (
        "mummer4.mumreference.raw.delta" if tool == "mummer4" else "ramag.mumreference.raw.delta"
    )
    converted = attempt / "prediction.converter.maf"
    command_evidence(
        [
            sys.executable,
            str(package / "tools/normalize_prediction.py"),
            "--prediction",
            str(delta),
            "--format",
            "delta",
            "--reference",
            str(config["reference"]),
            "--query",
            str(config["query"]),
            "--output",
            str(converted),
            "--stats",
            str(attempt / "prediction.converter.stats.json"),
        ],
        attempt,
        "normalize-prediction",
    )
    canonical_maf = attempt / "prediction.canonical.maf"
    canonical_stats = canonicalize_with_timing(
        converted,
        canonical_maf,
        attempt / "canonicalize-sources.timing.json",
    )
    if not set(canonical_stats["canonical_sources_seen"]) <= CANONICAL_SET:
        raise ExperimentError("canonical prediction contains sources outside the whitelist")
    write_json_exclusive(attempt / "canonicalization.json", canonical_stats)
    command_evidence(
        [
            sys.executable,
            str(package / "tools/validate_maf.py"),
            "--maf",
            str(canonical_maf),
            "--reference",
            str(config["reference"]),
            "--query",
            str(config["query"]),
            "--output-json",
            str(attempt / "prediction.validation.json"),
            "--allow-singletons",
        ],
        attempt,
        "validate-canonical-maf",
    )
    score = attempt / "score"
    command_evidence(
        [
            sys.executable,
            str(package / "tools/evaluate.py"),
            "--prediction",
            str(canonical_maf),
            "--format",
            "maf",
            "--name",
            tool,
            "--truth-profile",
            "all-homology",
            "--output",
            str(score),
            "--samples",
            str(config["samples"]),
            "--near",
            str(config["near"]),
            "--seeds",
            ",".join(str(seed) for seed in config["seeds"]),
        ],
        attempt,
        "evaluate",
    )
    summary = load_json(score / "summary.json")
    if summary.get("seeds") != config["seeds"] or int(summary.get("samples", 0)) != int(
        config["samples"]
    ) or int(summary.get("near", -1)) != int(config["near"]):
        raise ExperimentError(f"{tool} evaluation did not preserve fixed sample parameters")
    if not (score / "EVALUATION_COMPLETE").is_file():
        raise ExperimentError(f"{tool} package evaluator did not publish completion")
    xml_files = sorted((score / "raw").glob("mafComparator.seed-*.xml"))
    if len(xml_files) != len(config["seeds"]) or any(path.stat().st_size == 0 for path in xml_files):
        raise ExperimentError(f"{tool} evaluation lacks three non-empty comparator XML files")
    accept_attempt(base, name, attempt, "EVALUATION_ACCEPTED")
    return attempt


def run_evaluate(root: Path, config: dict[str, object]) -> None:
    if (root / "EVALUATIONS_COMPLETE").is_file():
        return
    if not (root / "FULL_RUNS_COMPLETE").is_file():
        raise ExperimentError("evaluation requires two accepted full runs")
    if not (root / "ALIGNMENT_SPEED_ACCEPTED").is_file() or (
        root / "ALIGNMENT_SPEED_REJECTED"
    ).exists():
        raise ExperimentError("evaluation requires the accepted early alignment speed gate")
    speed_gate = load_json(root / "alignment-speed-gate.json")
    conditions = speed_gate.get("conditions")
    if (
        speed_gate.get("schema")
        != "ramag.human-chimp-alignment-speed-gate.v1"
        or speed_gate.get("status") != "accepted"
        or speed_gate.get("reason") != "dual_clock_passed"
        or not isinstance(conditions, dict)
        or any(value is not True for value in conditions.values())
    ):
        raise ExperimentError("evaluation speed-gate evidence is invalid")
    for tool in ("mummer4", "ramag"):
        run_evaluation_attempt(root, tool, config)
    (root / "EVALUATIONS_COMPLETE").write_text("success\n", encoding="ascii")


def execute(root: Path, config: dict[str, object], stage: str) -> None:
    stages = ["preflight", "smoke", "full", "evaluate", "summarize"]
    target = "summarize" if stage == "all" else stage
    for current in stages[: stages.index(target) + 1]:
        if current == "preflight":
            run_preflight(root, config)
        elif current == "smoke":
            run_smoke(root, config)
        elif current == "full":
            run_full(root, config)
        elif current == "evaluate":
            run_evaluate(root, config)
        elif current == "summarize":
            if (root / "COMPARISON_COMPLETE").is_file():
                return
            result = summarize(root)
            performance_gate = result.get("performance_gate")
            clock_gate = result.get("clock_robustness_gate")
            if (
                not isinstance(performance_gate, dict)
                or performance_gate.get("status") != "passed"
                or not isinstance(clock_gate, dict)
                or clock_gate.get("status") != "passed"
            ):
                raise ExperimentError(
                    "GNU-time performance or cross-tool clock-robustness gate failed; "
                    "comparison artifacts were preserved without COMPARISON_COMPLETE. "
                    "Use a new exclusive experiment root for a replacement observation."
                )


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--stage",
        choices=("preflight", "smoke", "full", "evaluate", "summarize", "all"),
        default="preflight",
        help="furthest stage to execute; default is the non-alignment preflight",
    )
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--result-root", type=Path)
    parser.add_argument("--package", type=Path, default=DEFAULT_PACKAGE)
    parser.add_argument("--repo", type=Path, default=DEFAULT_REPO)
    parser.add_argument("--mummer-root", type=Path, default=DEFAULT_MUMMER_ROOT)
    parser.add_argument(
        "--ramag-binary",
        type=Path,
        default=DEFAULT_REPO / "build/ramag",
        help="fresh validated Release binary to snapshot into the result root",
    )
    parser.add_argument("--sample-interval", type=float, default=1.0)
    args = parser.parse_args(argv)
    if args.sample_interval != 1.0:
        parser.error("--sample-interval is fixed at 1.0 seconds for this protocol")
    return args


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if args.resume:
        root, config = load_resume_root(args)
    else:
        root, config = create_result_root(args)
    print(f"result_root={root}", flush=True)
    execute(root, config, args.stage)
    print(f"completed_stage={args.stage}", flush=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ExperimentError as error:
        print(f"run_human_chimp_benchmark: {error}", file=sys.stderr)
        raise SystemExit(2)

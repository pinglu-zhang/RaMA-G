#!/usr/bin/env python3

from __future__ import annotations

import copy
import datetime as dt
import gzip
import hashlib
import json
import subprocess
import sys
import tempfile
import time
import unittest
from argparse import Namespace
from pathlib import Path
from unittest import mock


MODULE_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(MODULE_DIR))

from canonicalize_maf_sources import (  # noqa: E402
    ALIASES,
    CANONICAL_SOURCES,
    CanonicalizationError,
    canonicalize,
)
from run_with_metrics import (  # noqa: E402
    assess_clock_consistency,
    parse_elapsed,
    parse_gnu_time,
    run,
)
from run_human_chimp_benchmark import (  # noqa: E402
    AlignmentTimeoutError,
    ClockQuarantineError,
    ExperimentError,
    FULL_MAM_TASK_SHAPE,
    OPENMP_ENVIRONMENT,
    audit_benchmark_package_integrity,
    build_completed_alignment_speed_gate,
    copy_source_snapshot,
    create_result_root,
    expected_mam_task_shape,
    load_resume_root,
    parse_args,
    record_full_launch_gate,
    resolve_ramag_commit,
    run_alignment_attempt,
    run_evaluate,
    run_full,
    validate_accepted_metrics,
    validate_ramag_manifest,
    warm_cache_inputs_for_tool,
)
from summarize_comparison import SummaryError, summarize  # noqa: E402

FIXTURE_REFERENCE = Path("/fixture/reference.fa")
FIXTURE_QUERY = Path("/fixture/query.fa")


def make_fixed_input_fixture(
    root: Path,
) -> tuple[Path, Path, dict[str, object], dict[str, str], dict[str, str]]:
    package = root / "package"
    inputs = package / "inputs"
    (package / "results").mkdir(parents=True)
    inputs.mkdir()
    reference = inputs / "simHuman.fa"
    query = inputs / "simChimp.fa"
    reference.write_bytes(b">simHuman.chrA\nACGT\n")
    query.write_bytes(b">simChimp.chrA\nTGCA\n")
    Path(str(reference) + ".fai").write_text(
        "simHuman.chrA\t4\t16\t4\t5\n", encoding="ascii"
    )
    Path(str(query) + ".fai").write_text(
        "simChimp.chrA\t4\t15\t4\t5\n", encoding="ascii"
    )
    fasta_hashes = {
        path.name: hashlib.sha256(path.read_bytes()).hexdigest()
        for path in (reference, query)
    }
    fai_paths = (Path(str(reference) + ".fai"), Path(str(query) + ".fai"))
    fai_hashes = {
        path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in fai_paths
    }
    config: dict[str, object] = {
        "package": str(package.resolve()),
        "reference": str(reference.resolve()),
        "query": str(query.resolve()),
        "input_sha256": fasta_hashes,
        "input_fai_sha256": fai_hashes,
    }
    return reference, query, config, fasta_hashes, fai_hashes


def make_package_integrity_fixture(root: Path) -> dict[str, object]:
    package = root / "package"
    truth_dir = package / "truth"
    bin_dir = package / "tools/bin"
    truth_dir.mkdir(parents=True)
    bin_dir.mkdir(parents=True)
    truth_names = (
        "simHuman-simChimp.all-homology.maf.gz",
        "simHuman-simChimp.no-paralogy.maf.gz",
        "simHuman-simChimp.single-copy-compat.maf.gz",
    )
    for index, filename in enumerate(truth_names, start=1):
        (truth_dir / filename).write_bytes(f"truth-{index}\n".encode("ascii"))
    version_lines = {
        "mafComparator": ("fixture mafComparator 0.9", "fixture build commit"),
        "mafPairCounter": ("fixture mafPairCounter 0.1", "fixture build commit"),
    }
    for filename, lines in version_lines.items():
        executable = bin_dir / filename
        executable.write_text(
            "#!/bin/sh\nprintf '%s\\n' "
            + " ".join(repr(line) for line in lines)
            + "\n",
            encoding="ascii",
        )
        executable.chmod(0o755)

    artifact_paths = [
        *(truth_dir / filename for filename in truth_names),
        *(bin_dir / filename for filename in version_lines),
    ]
    artifacts: dict[str, object] = {}
    for path in artifact_paths:
        relative = path.relative_to(package).as_posix()
        artifacts[relative] = {
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "bytes": path.stat().st_size,
        }
    profile_names = ("all-homology", "no-paralogy", "single-copy-compat")
    manifest = {
        "schema_version": 1,
        "status": "validated",
        "package_name": "alignathon-sim-human-chimp-v1",
        "headline_truth": "all-homology",
        "artifacts": artifacts,
        "truth_profiles": {
            profile: {
                "artifact": f"truth/{filename}",
                "sha256": artifacts[f"truth/{filename}"]["sha256"],
            }
            for profile, filename in zip(profile_names, truth_names, strict=True)
        },
    }
    manifest_path = package / "MANIFEST.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    manifest_hash = hashlib.sha256(manifest_path.read_bytes()).hexdigest()
    sums_entries = {
        "MANIFEST.json": manifest_hash,
        **{
            relative: str(metadata["sha256"])
            for relative, metadata in artifacts.items()
        },
    }
    sums_path = package / "SHA256SUMS"
    sums_path.write_text(
        "".join(
            f"{digest}  ./{relative}\n"
            for relative, digest in sorted(sums_entries.items())
        ),
        encoding="ascii",
    )
    sums_hash = hashlib.sha256(sums_path.read_bytes()).hexdigest()
    (package / "PACKAGE_COMPLETE").write_text(
        "status=validated\n"
        f"manifest_sha256={manifest_hash}\n"
        f"sha256sums_sha256={sums_hash}\n",
        encoding="ascii",
    )
    return {
        "package": package,
        "manifest_hash": manifest_hash,
        "sums_hash": sums_hash,
        "truth_hashes": {
            filename: artifacts[f"truth/{filename}"]["sha256"]
            for filename in truth_names
        },
        "evaluator_hashes": {
            filename: artifacts[f"tools/bin/{filename}"]["sha256"]
            for filename in version_lines
        },
        "version_lines": version_lines,
    }


def utc_from_ns(realtime_ns: int) -> str:
    seconds, nanoseconds = divmod(realtime_ns, 1_000_000_000)
    value = dt.datetime.fromtimestamp(seconds, tz=dt.timezone.utc).replace(
        microsecond=nanoseconds // 1_000
    )
    return value.isoformat(timespec="microseconds").replace("+00:00", "Z")


def write_resource_clocks(
    path: Path,
    monotonic_values: list[int],
    realtime_values: list[int],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    started = monotonic_values[0]
    lines = ["sample_index\telapsed_seconds\tmonotonic_ns\trealtime_ns\tutc"]
    for index, (monotonic_ns, realtime_ns) in enumerate(
        zip(monotonic_values, realtime_values, strict=True), start=1
    ):
        lines.append(
            f"{index}\t{(monotonic_ns - started) / 1e9:.9f}\t{monotonic_ns}\t"
            f"{realtime_ns}\t{utc_from_ns(realtime_ns)}"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def consistent_clock_fields(
    elapsed_seconds: float = 100.0,
) -> dict[str, object]:
    monotonic_started = 1_000_000_000_000
    monotonic_finished = monotonic_started + round(elapsed_seconds * 1e9)
    realtime_started = 1_700_000_000_000_000_000
    realtime_finished = realtime_started + round(elapsed_seconds * 1e9)
    started_utc = utc_from_ns(realtime_started)
    finished_utc = utc_from_ns(realtime_finished)
    return {
        "started_utc": started_utc,
        "finished_utc": finished_utc,
        "runner_wall_seconds": elapsed_seconds,
        "timeout": {
            "enabled": False,
            "clock": "CLOCK_MONOTONIC",
            "limit_seconds": 0.0,
            "exceeded": False,
            "elapsed_at_signal_seconds": None,
            "sigterm_sent": False,
            "sigkill_sent": False,
            "deadline_overshoot_seconds": None,
        },
        "clock_measurements": {
            "monotonic_clock": "CLOCK_MONOTONIC",
            "monotonic_started_ns": monotonic_started,
            "monotonic_finished_ns": monotonic_finished,
            "monotonic_elapsed_seconds": elapsed_seconds,
            "realtime_clock": "CLOCK_REALTIME",
            "realtime_started_ns": realtime_started,
            "realtime_finished_ns": realtime_finished,
            "realtime_elapsed_seconds": elapsed_seconds,
            "started_utc": started_utc,
            "finished_utc": finished_utc,
        },
    }

class CanonicalizerTests(unittest.TestCase):
    def test_exact_aliases_and_canonical_names(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "input.maf"
            destination = root / "output.maf"
            lines = ["##maf version=1 scoring=test\n", "# untouched comment\n"]
            for index, canonical in enumerate(CANONICAL_SOURCES):
                name = next(alias for alias, target in ALIASES.items() if target == canonical)
                # One row per block deliberately avoids a source collision while
                # exercising every exact alias.
                lines.extend(
                    [
                        f"a score={index}\n",
                        f"s  {name}\t0 4 + 4 ACGT\n",
                        "\n",
                    ]
                )
            lines.extend(
                [
                    "a\n",
                    "s simHuman.chrA 0 4 + 4 ACGT\n",
                    "s simChimp.chrA 0 4 + 4 ACGT\n",
                    "\n",
                ]
            )
            source.write_text("".join(lines), encoding="ascii", newline="")

            stats = canonicalize(source, destination)
            output = destination.read_text(encoding="ascii")
            self.assertNotIn("simHuman.simHuman", output)
            self.assertNotIn("simChimp.simChimp", output)
            self.assertEqual(stats["total_replacements"], 8)
            self.assertEqual(stats["s_rows"], 10)
            self.assertTrue(stats["non_s_lines_unchanged"])
            self.assertTrue(stats["s_line_payload_unchanged"])
            self.assertEqual(set(stats["canonical_sources_seen"]), set(CANONICAL_SOURCES))
            self.assertTrue(output.endswith("\n\n"), "line endings must be preserved")

    def test_gzip_input(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "input.maf.gz"
            with gzip.open(source, "wb") as handle:
                handle.write(
                    b"##maf version=1\n\na\ns simHuman.simHuman.chrA 0 1 + 1 A\n\n"
                )
            destination = root / "output.maf"
            stats = canonicalize(source, destination)
            self.assertEqual(stats["total_replacements"], 1)
            self.assertIn("s simHuman.chrA", destination.read_text(encoding="ascii"))

    def test_unknown_and_tripled_names_are_rejected_atomically(self) -> None:
        for source_name in ("unknown.chrA", "simHuman.simHuman.simHuman.chrA"):
            with self.subTest(source_name=source_name), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                source = root / "input.maf"
                output = root / "output.maf"
                source.write_text(
                    f"##maf version=1\n\na\ns {source_name} 0 1 + 1 A\n\n",
                    encoding="ascii",
                )
                with self.assertRaises(CanonicalizationError):
                    canonicalize(source, output)
                self.assertFalse(output.exists())
                self.assertEqual(list(root.glob(".*.tmp.*")), [])

    def test_normalization_collision_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "input.maf"
            output = root / "output.maf"
            source.write_text(
                "##maf version=1\n\n"
                "a\n"
                "s simHuman.chrA 0 1 + 1 A\n"
                "s simHuman.simHuman.chrA 1 1 + 2 C\n\n",
                encoding="ascii",
            )
            with self.assertRaisesRegex(CanonicalizationError, "collision"):
                canonicalize(source, output)
            self.assertFalse(output.exists())

    def test_refuses_overwrite(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "input.maf"
            output = root / "output.maf"
            source.write_text("a\ns simHuman.chrA 0 1 + 1 A\n", encoding="ascii")
            output.write_text("keep\n", encoding="ascii")
            with self.assertRaisesRegex(CanonicalizationError, "overwrite"):
                canonicalize(source, output)
            self.assertEqual(output.read_text(encoding="ascii"), "keep\n")


class MetricsTests(unittest.TestCase):
    def test_elapsed_parser(self) -> None:
        self.assertAlmostEqual(parse_elapsed("1:02.50"), 62.5)
        self.assertAlmostEqual(parse_elapsed("2:03:04"), 7384.0)

    def test_time_fixture_parser(self) -> None:
        fixture = MODULE_DIR / "tests/fixtures/gnu-time-v.txt"
        parsed = parse_gnu_time(fixture)
        self.assertEqual(parsed["user_seconds"], 12.5)
        self.assertEqual(parsed["system_seconds"], 1.25)
        self.assertEqual(parsed["elapsed_seconds"], 14.0)
        self.assertEqual(parsed["maximum_resident_set_kbytes"], 123456)
        self.assertEqual(parsed["exit_status"], 0)

    def test_runner_records_structured_metrics(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "run"
            environment = {"OMP_DYNAMIC": "FALSE", "OMP_NUM_THREADS": "3"}
            code = run(
                output,
                [
                    sys.executable,
                    "-c",
                    (
                        "import time\n"
                        "payload = bytearray(16 * 1024 * 1024)\n"
                        "for offset in range(0, len(payload), 4096):\n"
                        "    payload[offset] = 1\n"
                        "time.sleep(0.5)\n"
                    ),
                ],
                0.01,
                environment,
            )
            self.assertEqual(code, 0)
            metrics = json.loads((output / "metrics.json").read_text(encoding="utf-8"))
            self.assertEqual(metrics["status"], "success")
            self.assertGreaterEqual(metrics["sample_count"], 1)
            self.assertIsNotNone(metrics["gnu_time"])
            self.assertEqual(metrics["environment_overrides"], environment)
            self.assertEqual(metrics["schema"], "ramag.command-metrics.v3")
            self.assertEqual(
                metrics["timeout"],
                {
                    "clock": "CLOCK_MONOTONIC",
                    "deadline_overshoot_seconds": None,
                    "elapsed_at_signal_seconds": None,
                    "enabled": False,
                    "exceeded": False,
                    "limit_seconds": 0.0,
                    "process_group_alive_after_termination": None,
                    "process_group_id": metrics["timeout"]["process_group_id"],
                    "sigkill_sent": False,
                    "signal_monotonic_ns": None,
                    "sigterm_sent": False,
                    "termination": None,
                    "termination_grace_seconds": 10.0,
                },
            )
            clocks = metrics["clock_measurements"]
            self.assertGreater(clocks["monotonic_elapsed_seconds"], 0.0)
            self.assertGreater(clocks["realtime_elapsed_seconds"], 0.0)
            self.assertTrue((output / "resources.clock.tsv").is_file())
            self.assertTrue((output / "RUN_COMPLETE").is_file())
            accepted = validate_accepted_metrics(output / "metrics.json")
            self.assertEqual(accepted["exit_code"], 0)
            header = (output / "resources.tree.tsv").read_text(encoding="utf-8").splitlines()[0]
            self.assertEqual(
                header,
                "elapsed_seconds\tutc\tpid\tppid\tprocess\tthread_count\t"
                "cpu_percent\trss_bytes\tvms_bytes",
            )

    def test_runner_timeout_uses_monotonic_deadline_and_preserves_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "run"
            code = run(
                output,
                ["/bin/sh", "-c", "sleep 30"],
                0.01,
                timeout_seconds=0.05,
            )
            self.assertEqual(code, 124)
            metrics = json.loads((output / "metrics.json").read_text(encoding="utf-8"))
            timeout = metrics["timeout"]
            self.assertEqual(metrics["schema"], "ramag.command-metrics.v3")
            self.assertEqual(metrics["status"], "timed_out")
            self.assertEqual(timeout["clock"], "CLOCK_MONOTONIC")
            self.assertTrue(timeout["enabled"])
            self.assertTrue(timeout["exceeded"])
            self.assertGreater(timeout["elapsed_at_signal_seconds"], 0.05)
            self.assertGreaterEqual(timeout["deadline_overshoot_seconds"], 0.0)
            self.assertTrue(timeout["sigterm_sent"])
            self.assertFalse(timeout["sigkill_sent"])
            self.assertFalse(timeout["process_group_alive_after_termination"])
            self.assertTrue((output / "RUN_TIMED_OUT").is_file())
            self.assertTrue((output / "timeout.json").is_file())
            self.assertTrue((output / "exit-code.txt").is_file())
            self.assertTrue((output / "stdout.txt").is_file())
            self.assertTrue((output / "stderr.txt").is_file())
            self.assertTrue((output / "resources.tree.tsv").is_file())
            self.assertTrue((output / "resources.clock.tsv").is_file())
            self.assertFalse((output / "RUN_COMPLETE").exists())
            self.assertFalse((output / "RUN_ACCEPTED").exists())

    def test_runner_short_command_completes_before_enabled_deadline(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "run"
            code = run(
                output,
                ["/bin/sh", "-c", "sleep 0.02"],
                0.01,
                timeout_seconds=1.0,
            )
            self.assertEqual(code, 0)
            metrics = json.loads((output / "metrics.json").read_text(encoding="utf-8"))
            self.assertEqual(metrics["status"], "success")
            self.assertTrue(metrics["timeout"]["enabled"])
            self.assertFalse(metrics["timeout"]["exceeded"])
            self.assertTrue((output / "RUN_COMPLETE").is_file())
            self.assertFalse((output / "RUN_TIMED_OUT").exists())

    def test_runner_timeout_kills_the_runner_owned_descendant_group(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "run"
            code = run(
                output,
                [
                    "/bin/sh",
                    "-c",
                    "sleep 30 & child=$!; printf '%s\\n' \"$child\"; wait",
                ],
                0.01,
                timeout_seconds=0.1,
            )
            self.assertEqual(code, 124)
            child_pid = int((output / "stdout.txt").read_text(encoding="ascii").strip())
            self.assertFalse(Path(f"/proc/{child_pid}").exists())
            metrics = json.loads((output / "metrics.json").read_text(encoding="utf-8"))
            self.assertFalse(
                metrics["timeout"]["process_group_alive_after_termination"]
            )

    def test_runner_timeout_escalates_to_sigkill_after_ten_seconds(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "run"
            code = run(
                output,
                [
                    "/bin/sh",
                    "-c",
                    "trap '' TERM; while :; do sleep 30; done",
                ],
                0.01,
                timeout_seconds=0.05,
            )
            self.assertEqual(code, 124)
            metrics = json.loads((output / "metrics.json").read_text(encoding="utf-8"))
            timeout = metrics["timeout"]
            self.assertTrue(timeout["sigterm_sent"])
            self.assertTrue(timeout["sigkill_sent"])
            self.assertGreaterEqual(
                timeout["termination"]["grace_elapsed_seconds"], 10.0
            )
            self.assertFalse(timeout["process_group_alive_after_termination"])

    def test_clock_gate_accepts_three_consistent_clocks(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            path = root / "metrics.json"
            metrics = {
                "schema": "ramag.command-metrics.v3",
                "status": "success",
                "exit_code": 0,
                "gnu_time_parse_error": None,
                "gnu_time": {
                    "elapsed_seconds": 100.0,
                    "user_seconds": 10.0,
                    "system_seconds": 1.0,
                    "maximum_resident_set_kbytes": 1024,
                    "exit_status": 0,
                },
                "sample_count": 3,
                "observed_max_total_threads": 2,
                "observed_max_single_process_threads": 2,
                "observed_max_aligner_tree_threads": 1,
                "observed_max_aligner_single_process_threads": 1,
                "observed_process_tree_peak_rss_bytes": 4096,
                **consistent_clock_fields(100.0),
            }
            write_json(path, metrics)
            write_resource_clocks(
                root / "resources.clock.tsv",
                [1_000_000_000_000, 1_001_000_000_000, 1_002_000_000_000],
                [
                    1_700_000_000_000_000_000,
                    1_700_000_001_000_000_000,
                    1_700_000_002_000_000_000,
                ],
            )
            accepted = validate_accepted_metrics(
                path, require_clock_consistency=True
            )
            self.assertEqual(accepted["exit_code"], 0)
            report = assess_clock_consistency(metrics, root / "resources.clock.tsv")
            self.assertEqual(report["status"], "passed")
            self.assertEqual(report["resource_samples"]["realtime_jump_count"], 0)

    def test_clock_gate_retains_endpoint_drift_as_warning(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            path = root / "metrics.json"
            metrics = {
                "schema": "ramag.command-metrics.v3",
                "status": "success",
                "exit_code": 0,
                "gnu_time_parse_error": None,
                "gnu_time": {
                    "elapsed_seconds": 102.0,
                    "user_seconds": 10.0,
                    "system_seconds": 1.0,
                    "maximum_resident_set_kbytes": 1024,
                    "exit_status": 0,
                },
                "sample_count": 2,
                "observed_max_total_threads": 2,
                "observed_max_single_process_threads": 2,
                "observed_max_aligner_tree_threads": 1,
                "observed_max_aligner_single_process_threads": 1,
                "observed_process_tree_peak_rss_bytes": 4096,
                **consistent_clock_fields(100.0),
            }
            write_json(path, metrics)
            write_resource_clocks(
                root / "resources.clock.tsv",
                [1_000_000_000_000, 1_001_000_000_000],
                [1_700_000_000_000_000_000, 1_700_000_001_000_000_000],
            )
            accepted = validate_accepted_metrics(
                path, require_clock_consistency=True
            )
            self.assertEqual(accepted["exit_code"], 0)
            report = assess_clock_consistency(metrics, root / "resources.clock.tsv")
            self.assertEqual(report["status"], "warning")
            self.assertTrue(report["warnings"])

    def test_clock_gate_warns_on_periodic_resource_realtime_jumps(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            metrics = {
                "gnu_time": {"elapsed_seconds": 100.0},
                "sample_count": 3,
                **consistent_clock_fields(100.0),
            }
            write_resource_clocks(
                root / "resources.clock.tsv",
                [1_000_000_000_000, 1_001_000_000_000, 1_002_000_000_000],
                [
                    1_700_000_000_000_000_000,
                    1_700_000_001_200_000_000,
                    1_700_000_002_000_000_000,
                ],
            )
            report = assess_clock_consistency(metrics, root / "resources.clock.tsv")
            self.assertEqual(report["status"], "warning")
            self.assertEqual(report["resource_samples"]["realtime_jump_count"], 2)
            self.assertTrue(
                report["resource_samples"]["periodic_realtime_jump_detected"]
            )

    def test_clock_gate_hard_quarantines_realtime_backwards(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            path = root / "metrics.json"
            metrics = {
                "schema": "ramag.command-metrics.v3",
                "status": "success",
                "exit_code": 0,
                "gnu_time_parse_error": None,
                "gnu_time": {
                    "elapsed_seconds": 100.0,
                    "user_seconds": 10.0,
                    "system_seconds": 1.0,
                    "maximum_resident_set_kbytes": 1024,
                    "exit_status": 0,
                },
                "sample_count": 2,
                "observed_max_total_threads": 2,
                "observed_max_single_process_threads": 2,
                "observed_max_aligner_tree_threads": 1,
                "observed_max_aligner_single_process_threads": 1,
                "observed_process_tree_peak_rss_bytes": 4096,
                **consistent_clock_fields(100.0),
            }
            write_json(path, metrics)
            write_resource_clocks(
                root / "resources.clock.tsv",
                [1_000_000_000_000, 1_001_000_000_000],
                [1_700_000_001_000_000_000, 1_700_000_000_000_000_000],
            )
            with self.assertRaises(ClockQuarantineError):
                validate_accepted_metrics(path, require_clock_consistency=True)

    def test_acceptance_rejects_unparsed_or_empty_metrics(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "metrics.json"
            write_json(
                path,
                {
                    "schema": "ramag.command-metrics.v3",
                    "status": "success",
                    "exit_code": 0,
                    "gnu_time_parse_error": "broken fixture",
                    "gnu_time": None,
                    "sample_count": 0,
                    "observed_max_total_threads": 0,
                    "observed_max_single_process_threads": 0,
                    "observed_process_tree_peak_rss_bytes": 0,
                },
            )
            with self.assertRaisesRegex(ExperimentError, "not parseable"):
                validate_accepted_metrics(path)

    def test_acceptance_rejects_impossible_thread_hierarchy(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "metrics.json"
            write_json(
                path,
                {
                    "schema": "ramag.command-metrics.v3",
                    "status": "success",
                    "exit_code": 0,
                    "gnu_time_parse_error": None,
                    "gnu_time": {
                        "elapsed_seconds": 1.0,
                        "user_seconds": 1.0,
                        "system_seconds": 0.0,
                        "maximum_resident_set_kbytes": 1024,
                        "exit_status": 0,
                    },
                    "sample_count": 1,
                    "observed_max_total_threads": 2,
                    "observed_max_single_process_threads": 2,
                    "observed_max_aligner_tree_threads": 3,
                    "observed_max_aligner_single_process_threads": 3,
                    "observed_process_tree_peak_rss_bytes": 4096,
                },
            )
            with self.assertRaisesRegex(ExperimentError, "hierarchy"):
                validate_accepted_metrics(path)


def make_ramag_manifest(
    shape: dict[str, int] = FULL_MAM_TASK_SHAPE,
) -> tuple[dict[str, object], dict[str, object]]:
    commit = "bdb67c6de5daddd8a005640de73d96549d2575f4"
    seqpro_commit = "6781cadcf81a0da53d7573444594c1484947017c"
    workers = min(shape["total_tasks"], 16)
    backend = "caps32" if shape == FULL_MAM_TASK_SHAPE else "divsufsort32"
    config: dict[str, object] = {
        "expected_sufkit_commit": commit,
        "expected_seqpro_commit": seqpro_commit,
        "ramag_selection_mode": "one-to-one",
        "threads": 16,
        "min_match": 20,
        "max_gap": 90,
        "diag_diff": 5,
        "diag_factor": 0.12,
        "min_cluster": 65,
        "break_length": 200,
        "max_dp_cells": 4_000_000,
    }
    manifest: dict[str, object] = {
        "status": "success",
        "exit_code": 0,
        "dependencies": {
            "sufkit_commit": commit,
            "seqpro_commit": seqpro_commit,
        },
        "effective_config": {
            "reference": str(FIXTURE_REFERENCE),
            "query": str(FIXTURE_QUERY),
            "threads": 16,
            "formats": "delta",
            "seed_mode": "mumreference",
            "selection_mode": "one-to-one",
            "min_match": 20,
            "max_gap": 90,
            "diag_diff": 5,
            "diag_factor": 0.12,
            "min_cluster": 65,
            "break_length": 200,
            "max_dp_cells": 4_000_000,
        },
        "actual_routes": {
            "seed": "sufkit-full-sa:mumreference+tiled-mam-boundary-mem-v1+"
            "openmp-dynamic-stable-tile-boundary-tasks",
            "chain": "sparse-exact-edge-components-v1",
            "input_parallel": "openmp-sections-reference-query-v1",
        },
        "chaining": {"route": "sparse-exact-edge-components-v1"},
        "adapter_provenance": {
            "sufkit.commit": commit,
            "sufkit.source_state": "exact-commit-clean-at-configure",
            "sufkit.backend": backend,
            "sufkit.build.caps_min_reference_bases": str(64 * 1024 * 1024),
        },
        "build": {
            "openmp_enabled": True,
            "sufkit_divsufsort_openmp": True,
        },
        "threading": {
            "openmp_actual_requested_threads": 16,
            "openmp_runtime_max_threads": 16,
            "input_requested_workers": 2,
            "input_actual_workers": 2,
            "input_parallel_route": "openmp-sections-reference-query-v1",
            "seed_requested_threads": 16,
            "seed_scheduled_threads": workers,
            "seed_worker_threads": workers,
            "seed_task_count": shape["total_tasks"],
            "seed_tasks_completed": shape["total_tasks"],
            "seed_parallel_route": "openmp-dynamic-stable-tile-boundary-tasks",
            "chaining_requested_threads": 16,
            "chaining_worker_threads": 16,
            "extension_worker_threads": 16,
            "seed_chain_extension_parallel": True,
        },
        "mam_tiling": {
            "worker_cap": 16,
            "tile_bases": 4 * 1024 * 1024,
            "oriented_queries": shape["oriented_queries"],
            "tile_tasks": shape["tile_tasks"],
            "tile_tasks_completed": shape["tile_tasks"],
            "short_query_tasks": 0,
            "boundary_tasks": shape["boundary_tasks"],
            "boundary_tasks_completed": shape["boundary_tasks"],
            "tile_raw_mams": 120,
            "tile_globally_maximal_mams": 100,
            "boundary_raw_mems": 30,
            "boundary_patterns": 20,
            "boundary_reference_unique_patterns": 5,
            "boundary_recovered_mams": 5,
            "workspace_baseline_bytes": 1024,
            "workspace_peak_bytes": 1536,
            "resource_limits": {
                "boundary_mem_occurrences": 100,
                "workspace_bytes": 2048,
            },
        },
        "counts": {"mem_seeds": 30, "mam_seeds": 105, "selected_seeds": 100},
        "openmp": {
            "enabled": True,
            "runtime": "gomp",
            "runtime_max_threads": 16,
            "run_requested_threads": 16,
            "configured_requested_threads": 16,
            "dynamic": False,
            "sufkit_divsufsort_openmp": True,
        },
    }
    return manifest, config


class DriverTests(unittest.TestCase):
    def test_package_integrity_rehashes_manifest_truth_and_evaluators(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            fixture = make_package_integrity_fixture(Path(temporary))
            module = sys.modules["run_human_chimp_benchmark"]
            with mock.patch.multiple(
                module,
                EXPECTED_PACKAGE_MANIFEST_SHA256=fixture["manifest_hash"],
                EXPECTED_PACKAGE_SHA256SUMS_SHA256=fixture["sums_hash"],
                EXPECTED_TRUTH_SHA256=fixture["truth_hashes"],
                EXPECTED_EVALUATOR_SHA256=fixture["evaluator_hashes"],
                EXPECTED_EVALUATOR_VERSION_LINES=fixture["version_lines"],
            ):
                report = audit_benchmark_package_integrity(fixture["package"])
                self.assertEqual(report["status"], "success")
                self.assertEqual(report["artifact_count"], 6)
                self.assertEqual(
                    set(report["evaluator_versions"]),
                    {"mafComparator", "mafPairCounter"},
                )

                truth = (
                    fixture["package"]
                    / "truth/simHuman-simChimp.all-homology.maf.gz"
                )
                truth.write_bytes(b"tampered\n")
                tampered = audit_benchmark_package_integrity(fixture["package"])
                self.assertEqual(tampered["status"], "failed")
                self.assertTrue(
                    any("artifact SHA-256 mismatch" in failure for failure in tampered["failures"])
                )

    def test_package_integrity_rejects_wrong_executable_version(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            fixture = make_package_integrity_fixture(Path(temporary))
            module = sys.modules["run_human_chimp_benchmark"]
            wrong_versions = dict(fixture["version_lines"])
            wrong_versions["mafComparator"] = ("wrong version",)
            with mock.patch.multiple(
                module,
                EXPECTED_PACKAGE_MANIFEST_SHA256=fixture["manifest_hash"],
                EXPECTED_PACKAGE_SHA256SUMS_SHA256=fixture["sums_hash"],
                EXPECTED_TRUTH_SHA256=fixture["truth_hashes"],
                EXPECTED_EVALUATOR_SHA256=fixture["evaluator_hashes"],
                EXPECTED_EVALUATOR_VERSION_LINES=wrong_versions,
            ):
                report = audit_benchmark_package_integrity(fixture["package"])
            self.assertEqual(report["status"], "failed")
            self.assertIn(
                "fixed evaluator version mismatch: mafComparator", report["failures"]
            )

    def test_resume_rehashes_fixed_fasta_and_fai(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root_path = Path(temporary)
            reference, query, config, fasta_hashes, fai_hashes = (
                make_fixed_input_fixture(root_path)
            )
            package = Path(str(config["package"]))
            result_root = package / "results/preliminary-mam16-single-resume-fixture"
            result_root.mkdir()
            module = sys.modules["run_human_chimp_benchmark"]
            config.update(
                {
                    "schema": "ramag.human-chimp-preliminary-config.v3",
                    "threads": 16,
                    "repetitions": 1,
                    "openmp_environment": OPENMP_ENVIRONMENT,
                    "sample_interval_seconds": 1.0,
                    "clock_consistency_tolerance_fraction": 0.01,
                    "package_manifest_sha256": module.EXPECTED_PACKAGE_MANIFEST_SHA256,
                    "package_sha256sums_sha256": module.EXPECTED_PACKAGE_SHA256SUMS_SHA256,
                    "truth_sha256": module.EXPECTED_TRUTH_SHA256,
                    "evaluator_sha256": module.EXPECTED_EVALUATOR_SHA256,
                    "harness_scripts_sha256": {},
                    "expected_sufkit_commit": module.EXPECTED_SUFKIT_COMMIT,
                    "expected_seqpro_commit": module.EXPECTED_SEQPRO_COMMIT,
                    "metrics_schema": "ramag.command-metrics.v3",
                    "ramag_commit": "unknown",
                    "early_stop_policy": {
                        "enabled": True,
                        "equal_is_accepted": True,
                    },
                    "alignment_speed_policy": {
                        "equal_is_accepted": True,
                        "observations_per_tool": 1,
                    },
                    "quality_policy": {
                        "mode": "report-only",
                        "minimum_f1": None,
                        "comparative_f1_gate": False,
                    },
                }
            )
            write_json(result_root / "accepted-config.json", config)
            args = Namespace(result_root=result_root)
            with mock.patch.multiple(
                module,
                EXPECTED_FASTA_SHA256=fasta_hashes,
                EXPECTED_FAI_SHA256=fai_hashes,
            ):
                resumed, _ = load_resume_root(args)
                self.assertEqual(resumed, result_root.resolve())
                Path(str(query) + ".fai").write_text("tampered\n", encoding="ascii")
                with self.assertRaisesRegex(ExperimentError, "FASTA/FAI integrity"):
                    load_resume_root(args)

    def test_per_tool_warm_cache_reads_reference_then_query(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            reference, query, config, fasta_hashes, fai_hashes = (
                make_fixed_input_fixture(root)
            )
            attempt = root / "attempt-001"
            attempt.mkdir()
            module = sys.modules["run_human_chimp_benchmark"]
            with mock.patch.multiple(
                module,
                EXPECTED_FASTA_SHA256=fasta_hashes,
                EXPECTED_FAI_SHA256=fai_hashes,
            ):
                report = warm_cache_inputs_for_tool(
                    "mummer4", attempt, reference, query, config
                )
            self.assertEqual(report["status"], "success")
            self.assertEqual(report["read_order"], ["reference", "query"])
            self.assertEqual(
                [entry["role"] for entry in report["inputs"]],
                ["reference", "query"],
            )
            self.assertEqual(report["input_index_integrity"]["status"], "success")
            self.assertTrue((attempt / "WARM_CACHE_COMPLETE").is_file())
            self.assertFalse((attempt / "WARM_CACHE_FAILED").exists())

    def test_warm_cache_hash_mismatch_refuses_aligner_launch(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            reference, query, config, fasta_hashes, fai_hashes = (
                make_fixed_input_fixture(root)
            )
            config["input_sha256"] = {
                reference.name: "0" * 64,
                query.name: fasta_hashes[query.name],
            }
            config.update({
                "mummer_nucmer": "/bin/false",
                "threads": 16,
                "min_match": 20,
                "min_cluster": 65,
                "max_gap": 90,
                "diag_diff": 5,
                "diag_factor": 0.12,
                "break_length": 200,
                "sample_interval_seconds": 1.0,
                "openmp_environment": OPENMP_ENVIRONMENT,
            })
            module = sys.modules["run_human_chimp_benchmark"]
            with (
                mock.patch.multiple(
                    module,
                    EXPECTED_FASTA_SHA256=fasta_hashes,
                    EXPECTED_FAI_SHA256=fai_hashes,
                ),
                mock.patch.object(module, "run_with_metrics") as runner,
                self.assertRaisesRegex(ExperimentError, "aligner was not launched"),
            ):
                run_alignment_attempt(
                    root / "runs/mummer4",
                    "mummer4",
                    reference,
                    query,
                    config,
                    full_acceptance=True,
                )
            runner.assert_not_called()
            attempt = root / "runs/mummer4/attempt-001"
            report = json.loads((attempt / "warm-cache.json").read_text())
            self.assertEqual(report["status"], "failed")
            self.assertTrue((attempt / "WARM_CACHE_FAILED").is_file())
            self.assertFalse((attempt / "timing").exists())
            self.assertFalse((attempt / "RUN_ACCEPTED").exists())
            self.assertFalse((attempt.parent / "SELECTED_ATTEMPT").exists())

    def test_clock_quarantine_preserves_attempt_without_selection(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            reference, query, config, fasta_hashes, fai_hashes = (
                make_fixed_input_fixture(root)
            )
            config.update({
                "mummer_nucmer": "/bin/true",
                "threads": 16,
                "min_match": 20,
                "min_cluster": 65,
                "max_gap": 90,
                "diag_diff": 5,
                "diag_factor": 0.12,
                "break_length": 200,
                "sample_interval_seconds": 1.0,
                "openmp_environment": OPENMP_ENVIRONMENT,
            })

            def fake_runner(
                output: Path,
                _command: list[str],
                _interval: float,
                environment: dict[str, str],
                _timeout_seconds: float | None = None,
            ) -> int:
                output.mkdir(parents=True)
                monotonic_started = time.monotonic_ns()
                realtime_started = time.time_ns()
                elapsed = 2.0
                monotonic_finished = monotonic_started + 2_000_000_000
                realtime_finished = realtime_started + 2_000_000_000
                started_utc = utc_from_ns(realtime_started)
                finished_utc = utc_from_ns(realtime_finished)
                metrics = {
                    "schema": "ramag.command-metrics.v3",
                    "status": "success",
                    "exit_code": 0,
                    "gnu_time_parse_error": None,
                    "gnu_time": {
                        "elapsed_seconds": 2.1,
                        "user_seconds": 1.0,
                        "system_seconds": 0.1,
                        "maximum_resident_set_kbytes": 1024,
                        "exit_status": 0,
                    },
                    "sample_count": 2,
                    "observed_max_total_threads": 2,
                    "observed_max_single_process_threads": 2,
                    "observed_max_aligner_tree_threads": 1,
                    "observed_max_aligner_single_process_threads": 1,
                    "observed_process_tree_peak_rss_bytes": 4096,
                    "environment_overrides": environment,
                    "started_utc": started_utc,
                    "finished_utc": finished_utc,
                    "runner_wall_seconds": elapsed,
                    "timeout": {
                        "enabled": False,
                        "clock": "CLOCK_MONOTONIC",
                        "limit_seconds": 0.0,
                        "exceeded": False,
                        "elapsed_at_signal_seconds": None,
                        "sigterm_sent": False,
                        "sigkill_sent": False,
                        "deadline_overshoot_seconds": None,
                    },
                    "clock_measurements": {
                        "monotonic_started_ns": monotonic_started,
                        "monotonic_finished_ns": monotonic_finished,
                        "monotonic_elapsed_seconds": elapsed,
                        "realtime_started_ns": realtime_started,
                        "realtime_finished_ns": realtime_finished,
                        "realtime_elapsed_seconds": elapsed,
                        "started_utc": started_utc,
                        "finished_utc": finished_utc,
                    },
                }
                write_json(output / "metrics.json", metrics)
                write_resource_clocks(
                    output / "resources.clock.tsv",
                    [monotonic_started, monotonic_started + 1_000_000_000],
                    [realtime_started, realtime_started - 1_000_000_000],
                )
                return 0

            module = sys.modules["run_human_chimp_benchmark"]
            with (
                mock.patch.multiple(
                    module,
                    EXPECTED_FASTA_SHA256=fasta_hashes,
                    EXPECTED_FAI_SHA256=fai_hashes,
                ),
                mock.patch.object(module, "run_with_metrics", side_effect=fake_runner),
                self.assertRaisesRegex(ExperimentError, "quarantined"),
            ):
                run_alignment_attempt(
                    root / "runs/mummer4",
                    "mummer4",
                    reference,
                    query,
                    config,
                    full_acceptance=True,
                )
            attempt = root / "runs/mummer4/attempt-001"
            self.assertTrue((attempt / "CLOCK_QUARANTINED").is_file())
            self.assertTrue((attempt / "clock-consistency.json").is_file())
            self.assertFalse((attempt / "CLOCK_ACCEPTED").exists())
            self.assertFalse((attempt / "RUN_ACCEPTED").exists())
            self.assertFalse((attempt.parent / "SELECTED_ATTEMPT").exists())

    def test_clock_warning_attempt_remains_selectable(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            reference, query, config, fasta_hashes, fai_hashes = (
                make_fixed_input_fixture(root)
            )
            config.update({
                "mummer_nucmer": "/bin/true",
                "threads": 16,
                "min_match": 20,
                "min_cluster": 65,
                "max_gap": 90,
                "diag_diff": 5,
                "diag_factor": 0.12,
                "break_length": 200,
                "sample_interval_seconds": 1.0,
                "openmp_environment": OPENMP_ENVIRONMENT,
            })

            def warning_runner(
                output: Path,
                _command: list[str],
                _interval: float,
                environment: dict[str, str],
                _timeout_seconds: float | None = None,
            ) -> int:
                output.mkdir(parents=True)
                monotonic_started = time.monotonic_ns()
                realtime_started = time.time_ns()
                monotonic_finished = monotonic_started + 2_000_000_000
                realtime_finished = realtime_started + 2_200_000_000
                started_utc = utc_from_ns(realtime_started)
                finished_utc = utc_from_ns(realtime_finished)
                metrics = {
                    "schema": "ramag.command-metrics.v3",
                    "status": "success",
                    "exit_code": 0,
                    "gnu_time_parse_error": None,
                    "gnu_time": {
                        "elapsed_seconds": 2.2,
                        "user_seconds": 1.0,
                        "system_seconds": 0.1,
                        "maximum_resident_set_kbytes": 1024,
                        "exit_status": 0,
                    },
                    "sample_count": 2,
                    "observed_max_total_threads": 2,
                    "observed_max_single_process_threads": 2,
                    "observed_max_aligner_tree_threads": 1,
                    "observed_max_aligner_single_process_threads": 1,
                    "observed_process_tree_peak_rss_bytes": 4096,
                    "environment_overrides": environment,
                    "started_utc": started_utc,
                    "finished_utc": finished_utc,
                    "runner_wall_seconds": 2.0,
                    "timeout": {
                        "enabled": False,
                        "clock": "CLOCK_MONOTONIC",
                        "limit_seconds": 0.0,
                        "exceeded": False,
                        "elapsed_at_signal_seconds": None,
                        "sigterm_sent": False,
                        "sigkill_sent": False,
                        "deadline_overshoot_seconds": None,
                    },
                    "clock_measurements": {
                        "monotonic_started_ns": monotonic_started,
                        "monotonic_finished_ns": monotonic_finished,
                        "monotonic_elapsed_seconds": 2.0,
                        "realtime_started_ns": realtime_started,
                        "realtime_finished_ns": realtime_finished,
                        "realtime_elapsed_seconds": 2.2,
                        "started_utc": started_utc,
                        "finished_utc": finished_utc,
                    },
                }
                write_json(output / "metrics.json", metrics)
                write_resource_clocks(
                    output / "resources.clock.tsv",
                    [monotonic_started, monotonic_started + 1_000_000_000],
                    [realtime_started, realtime_started + 1_100_000_000],
                )
                return 0

            module = sys.modules["run_human_chimp_benchmark"]
            with (
                mock.patch.multiple(
                    module,
                    EXPECTED_FASTA_SHA256=fasta_hashes,
                    EXPECTED_FAI_SHA256=fai_hashes,
                ),
                mock.patch.object(
                    module, "run_with_metrics", side_effect=warning_runner
                ),
                mock.patch.object(module, "format_validate"),
            ):
                selected = run_alignment_attempt(
                    root / "runs/mummer4",
                    "mummer4",
                    reference,
                    query,
                    config,
                    full_acceptance=True,
                )
            self.assertEqual(selected.name, "attempt-001")
            self.assertTrue((selected / "CLOCK_WARNING").is_file())
            self.assertTrue((selected / "CLOCK_ACCEPTED").is_file())
            self.assertTrue((selected / "RUN_ACCEPTED").is_file())
            self.assertFalse((selected / "CLOCK_QUARANTINED").exists())
            self.assertEqual(
                (selected.parent / "SELECTED_ATTEMPT").read_text().strip(),
                "attempt-001",
            )

    def test_full_launch_gate_records_each_check_and_all_limits(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            _reference, _query, config, fasta_hashes, fai_hashes = (
                make_fixed_input_fixture(root)
            )
            disk = type("DiskUsage", (), {"free": 60 * 1024**3})()
            module = sys.modules["run_human_chimp_benchmark"]
            with (
                mock.patch.multiple(
                    module,
                    EXPECTED_FASTA_SHA256=fasta_hashes,
                    EXPECTED_FAI_SHA256=fai_hashes,
                ),
                mock.patch.object(module, "active_aligners", return_value=[]),
                mock.patch.object(
                    module, "read_mem_available", return_value=24 * 1024**3
                ),
                mock.patch.object(module.shutil, "disk_usage", return_value=disk),
            ):
                first = record_full_launch_gate(root, "mummer4", config)
                second = record_full_launch_gate(root, "mummer4", config)
            self.assertEqual(first.name, "attempt-001")
            self.assertEqual(second.name, "attempt-002")
            report = json.loads((second / "launch-gate.json").read_text())
            self.assertEqual(report["tool"], "mummer4")
            self.assertEqual(report["status"], "success")
            self.assertEqual(report["available_memory_bytes"], 24 * 1024**3)
            self.assertEqual(report["free_disk_bytes"], 60 * 1024**3)
            self.assertEqual(report["fixed_input_integrity"]["status"], "success")
            self.assertTrue(report["checked_at_utc"].endswith("Z"))
            self.assertTrue((second / "LAUNCH_GATE_ACCEPTED").is_file())
            self.assertFalse((second.parent / "SELECTED_ATTEMPT").exists())

    def test_full_launch_gate_preserves_failure_without_acceptance(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            _reference, _query, config, fasta_hashes, fai_hashes = (
                make_fixed_input_fixture(root)
            )
            disk = type("DiskUsage", (), {"free": 49 * 1024**3})()
            conflicts = [{"pid": 123, "executable": "nucmer", "command": "nucmer"}]
            module = sys.modules["run_human_chimp_benchmark"]
            with (
                mock.patch.multiple(
                    module,
                    EXPECTED_FASTA_SHA256=fasta_hashes,
                    EXPECTED_FAI_SHA256=fai_hashes,
                ),
                mock.patch.object(module, "active_aligners", return_value=conflicts),
                mock.patch.object(
                    module, "read_mem_available", return_value=19 * 1024**3
                ),
                mock.patch.object(module.shutil, "disk_usage", return_value=disk),
                self.assertRaisesRegex(ExperimentError, "launch gate failed"),
            ):
                record_full_launch_gate(root, "ramag", config)
            attempt = root / "runs/launch-gates/ramag/attempt-001"
            report = json.loads((attempt / "launch-gate.json").read_text())
            self.assertEqual(report["status"], "failed")
            self.assertEqual(len(report["failures"]), 3)
            self.assertTrue((attempt / "LAUNCH_GATE_FAILED").is_file())
            self.assertFalse((attempt / "LAUNCH_GATE_ACCEPTED").exists())
            self.assertFalse((attempt.parent / "SELECTED_ATTEMPT").exists())

    def test_full_launch_gate_rejects_changed_fixed_fai(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            _reference, query, config, fasta_hashes, fai_hashes = (
                make_fixed_input_fixture(root)
            )
            Path(str(query) + ".fai").write_text("tampered\n", encoding="ascii")
            disk = type("DiskUsage", (), {"free": 60 * 1024**3})()
            module = sys.modules["run_human_chimp_benchmark"]
            with (
                mock.patch.multiple(
                    module,
                    EXPECTED_FASTA_SHA256=fasta_hashes,
                    EXPECTED_FAI_SHA256=fai_hashes,
                ),
                mock.patch.object(module, "active_aligners", return_value=[]),
                mock.patch.object(
                    module, "read_mem_available", return_value=24 * 1024**3
                ),
                mock.patch.object(module.shutil, "disk_usage", return_value=disk),
                self.assertRaisesRegex(ExperimentError, "launch gate failed"),
            ):
                record_full_launch_gate(root, "ramag", config)
            attempt = root / "runs/launch-gates/ramag/attempt-001"
            report = json.loads((attempt / "launch-gate.json").read_text())
            self.assertEqual(report["fixed_input_integrity"]["status"], "failed")
            self.assertTrue(
                any("FAI" in failure or ".fai" in failure for failure in report["failures"])
            )

    def test_ramag_manifest_gate_requires_sparse_route_and_clean_sufkit(self) -> None:
        manifest, config = make_ramag_manifest()
        validate_ramag_manifest(
            manifest,
            config,
            expected_task_shape=FULL_MAM_TASK_SHAPE,
            expected_reference=FIXTURE_REFERENCE,
            expected_query=FIXTURE_QUERY,
        )
        manifest["actual_routes"] = {"chain": "baseline-diagonal"}
        with self.assertRaisesRegex(ExperimentError, "sparse-exact"):
            validate_ramag_manifest(
                manifest,
                config,
                expected_task_shape=FULL_MAM_TASK_SHAPE,
                expected_reference=FIXTURE_REFERENCE,
                expected_query=FIXTURE_QUERY,
            )

        manifest, config = make_ramag_manifest()
        manifest["adapter_provenance"]["sufkit.backend"] = "divsufsort32"
        with self.assertRaisesRegex(ExperimentError, "SA-backend policy"):
            validate_ramag_manifest(
                manifest,
                config,
                expected_task_shape=FULL_MAM_TASK_SHAPE,
                expected_reference=FIXTURE_REFERENCE,
                expected_query=FIXTURE_QUERY,
            )

    def test_ramag_manifest_gate_verifies_external_fai_route(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            package = root / "package"
            inputs = package / "inputs"
            inputs.mkdir(parents=True)
            reference = inputs / "simHuman.fa"
            query = inputs / "simChimp.fa"
            reference.write_text(">simHuman.chrA\nACGT\n", encoding="ascii")
            query.write_text(">simChimp.chrA\nACGT\n", encoding="ascii")
            for fasta, name in (
                (reference, "simHuman.chrA"),
                (query, "simChimp.chrA"),
            ):
                Path(str(fasta) + ".fai").write_text(
                    f"{name}\t4\t{len(name) + 2}\t4\t5\n", encoding="ascii"
                )

            manifest, config = make_ramag_manifest()
            fasta_hashes = {
                path.name: hashlib.sha256(path.read_bytes()).hexdigest()
                for path in (reference, query)
            }
            fai_hashes = {
                Path(str(path) + ".fai").name: hashlib.sha256(
                    Path(str(path) + ".fai").read_bytes()
                ).hexdigest()
                for path in (reference, query)
            }
            config.update(
                {
                    "package": str(package.resolve()),
                    "reference": str(reference.resolve()),
                    "query": str(query.resolve()),
                    "input_sha256": fasta_hashes,
                    "input_fai_sha256": fai_hashes,
                }
            )
            manifest["effective_config"]["reference"] = str(reference.resolve())
            manifest["effective_config"]["query"] = str(query.resolve())
            manifest["actual_routes"]["input"] = (
                "seqpro-external-fai-mmap+caller-buffer-copy"
            )
            adapter = manifest["adapter_provenance"]
            work = root / "work"
            work.mkdir()
            for role, fasta in (("reference", reference), ("query", query)):
                source_fai = Path(str(fasta) + ".fai")
                active_fai = work / f"{role}.fai"
                active_fai.write_bytes(source_fai.read_bytes())
                prefix = f"seqpro.{role}."
                adapter[prefix + "fai"] = str(active_fai.resolve())
                adapter[prefix + "metadata"] = ""
                adapter[prefix + "source_fai"] = str(source_fai.resolve())
                adapter[prefix + "source_fai_bytes"] = str(source_fai.stat().st_size)
                adapter[prefix + "source_fai_sha256"] = hashlib.sha256(
                    source_fai.read_bytes()
                ).hexdigest()
                adapter[prefix + "source_fai.copy_status"] = "copied"
                adapter[prefix + "external_fai.adoption_status"] = "adopted"
                adapter[prefix + "build_action"] = "reused"
                adapter[prefix + "index_origin"] = "external-standard-fai"
                adapter[prefix + "verification"] = "structure-validated"

            module = sys.modules["run_human_chimp_benchmark"]
            with mock.patch.multiple(
                module,
                EXPECTED_FASTA_SHA256=fasta_hashes,
                EXPECTED_FAI_SHA256=fai_hashes,
            ):
                validate_ramag_manifest(
                    manifest,
                    config,
                    expected_task_shape=FULL_MAM_TASK_SHAPE,
                    expected_reference=reference,
                    expected_query=query,
                    verify_fixed_input_routes=True,
                )
                adapter["seqpro.query.verification"] = "metadata-validated"
                with self.assertRaisesRegex(ExperimentError, "external-FAI"):
                    validate_ramag_manifest(
                        manifest,
                        config,
                        expected_task_shape=FULL_MAM_TASK_SHAPE,
                        expected_reference=reference,
                        expected_query=query,
                        verify_fixed_input_routes=True,
                    )
                adapter["seqpro.query.verification"] = "structure-validated"
                Path(str(query) + ".fai").write_text("tampered\n", encoding="ascii")
                with self.assertRaisesRegex(ExperimentError, "FASTA/FAI integrity"):
                    validate_ramag_manifest(
                        manifest,
                        config,
                        expected_task_shape=FULL_MAM_TASK_SHAPE,
                        expected_reference=reference,
                        expected_query=query,
                        verify_fixed_input_routes=True,
                    )

        manifest, config = make_ramag_manifest()
        manifest["dependencies"]["seqpro_commit"] = "wrong"
        with self.assertRaisesRegex(ExperimentError, "dependency commit"):
            validate_ramag_manifest(
                manifest,
                config,
                expected_task_shape=FULL_MAM_TASK_SHAPE,
                expected_reference=FIXTURE_REFERENCE,
                expected_query=FIXTURE_QUERY,
            )

        manifest, config = make_ramag_manifest()
        manifest["actual_routes"]["seed"] = "whole-contig-mam"
        with self.assertRaisesRegex(ExperimentError, "sparse-exact"):
            validate_ramag_manifest(
                manifest,
                config,
                expected_task_shape=FULL_MAM_TASK_SHAPE,
                expected_reference=FIXTURE_REFERENCE,
                expected_query=FIXTURE_QUERY,
            )

    def test_tiled_mam_shape_is_computed_and_exactly_enforced(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            query = Path(temporary) / "query.fa"
            query.write_text(
                ">q\n" + "A" * (2 * 4 * 1024 * 1024 + 1) + "\n",
                encoding="ascii",
            )
            smoke_shape = {
                "oriented_queries": 2,
                "tile_tasks": 6,
                "boundary_tasks": 4,
                "total_tasks": 10,
            }
            self.assertEqual(expected_mam_task_shape(query), smoke_shape)
            smoke_manifest, config = make_ramag_manifest(smoke_shape)
            validate_ramag_manifest(
                smoke_manifest,
                config,
                expected_task_shape=smoke_shape,
                expected_reference=FIXTURE_REFERENCE,
                expected_query=FIXTURE_QUERY,
            )
            smoke_manifest["adapter_provenance"]["sufkit.backend"] = "caps32"
            with self.assertRaisesRegex(ExperimentError, "OpenMP runtime/build contract"):
                validate_ramag_manifest(
                    smoke_manifest,
                    config,
                    expected_task_shape=FULL_MAM_TASK_SHAPE,
                    expected_reference=FIXTURE_REFERENCE,
                    expected_query=FIXTURE_QUERY,
                )

            full_manifest, config = make_ramag_manifest()
            full_manifest["mam_tiling"]["tile_tasks"] = 91
            with self.assertRaisesRegex(ExperimentError, "OpenMP runtime/build contract"):
                validate_ramag_manifest(
                    full_manifest,
                    config,
                    expected_task_shape=FULL_MAM_TASK_SHAPE,
                    expected_reference=FIXTURE_REFERENCE,
                    expected_query=FIXTURE_QUERY,
                )

    def test_tiled_mam_statistical_invariants_are_enforced(self) -> None:
        mutations = (
            ("tile_globally_maximal_mams", 121),
            ("boundary_patterns", 31),
            ("boundary_reference_unique_patterns", 21),
            ("boundary_recovered_mams", 6),
            ("workspace_baseline_bytes", 1537),
        )
        for field, value in mutations:
            with self.subTest(field=field):
                manifest, config = make_ramag_manifest()
                manifest = copy.deepcopy(manifest)
                manifest["mam_tiling"][field] = value
                with self.assertRaisesRegex(ExperimentError, "statistical invariants"):
                    validate_ramag_manifest(
                        manifest,
                        config,
                        expected_task_shape=FULL_MAM_TASK_SHAPE,
                        expected_reference=FIXTURE_REFERENCE,
                        expected_query=FIXTURE_QUERY,
                    )

        manifest, config = make_ramag_manifest()
        manifest["counts"]["mam_seeds"] = 104
        with self.assertRaisesRegex(ExperimentError, "top-level counts"):
            validate_ramag_manifest(
                manifest,
                config,
                expected_task_shape=FULL_MAM_TASK_SHAPE,
                expected_reference=FIXTURE_REFERENCE,
                expected_query=FIXTURE_QUERY,
            )

    def test_sample_interval_is_fixed(self) -> None:
        with self.assertRaises(SystemExit):
            parse_args(["--sample-interval", "0.5"])

    def test_new_root_is_exclusive_and_fixed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            package = base / "package"
            (package / "results").mkdir(parents=True)
            repo = base / "repo"
            repo.mkdir()
            mummer = base / "mummer"
            mummer.mkdir()
            root = package / "results/preliminary-mam16-single-fixture"
            args = Namespace(
                package=package,
                result_root=root,
                ramag_binary=Path("/bin/true"),
                repo=repo,
                mummer_root=mummer,
                sample_interval=1.0,
            )
            created, config = create_result_root(args)
            self.assertEqual(created, root)
            self.assertEqual(config["threads"], 16)
            self.assertEqual(config["repetitions"], 1)
            self.assertEqual(config["ramag_seed_mode"], "mumreference")
            self.assertEqual(config["ramag_selection_mode"], "one-to-one")
            self.assertEqual(config["openmp_environment"], OPENMP_ENVIRONMENT)
            self.assertEqual(config["schema"], "ramag.human-chimp-preliminary-config.v3")
            self.assertEqual(config["ramag_commit"], "unknown")
            self.assertEqual(config["clock_consistency_tolerance_fraction"], 0.01)
            self.assertIn("reference-then-query", config["full_run_warm_cache_policy"])
            self.assertTrue((root / "accepted-config.json").is_file())
            self.assertTrue((root / "build/ramag").is_file())
            with self.assertRaisesRegex(ExperimentError, "overwrite"):
                create_result_root(args)

    def test_git_commit_is_recorded_and_dirty_tree_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            repo = base / "repo"
            repo.mkdir()
            subprocess.run(["git", "init", "--quiet", str(repo)], check=True)
            (repo / "tracked.txt").write_text("committed\n", encoding="ascii")
            subprocess.run(["git", "-C", str(repo), "add", "tracked.txt"], check=True)
            subprocess.run(
                [
                    "git",
                    "-C",
                    str(repo),
                    "-c",
                    "user.name=RaMA-G test",
                    "-c",
                    "user.email=ramag-test@example.invalid",
                    "commit",
                    "--quiet",
                    "-m",
                    "fixture",
                ],
                check=True,
            )
            expected = subprocess.check_output(
                ["git", "-C", str(repo), "rev-parse", "HEAD"], text=True
            ).strip()
            self.assertEqual(resolve_ramag_commit(repo), expected)
            (repo / "untracked.txt").write_text("dirty\n", encoding="ascii")
            with self.assertRaisesRegex(ExperimentError, "must be clean"):
                resolve_ramag_commit(repo)

    def test_source_snapshot_excludes_build_results_and_temp(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            repo = root / "repo"
            (repo / "src").mkdir(parents=True)
            (repo / "src/main.cpp").write_text("int main() {}\n", encoding="ascii")
            for ignored in (
                "build-release",
                "deps-src",
                "results",
                "temp",
                "tmp",
                "work",
            ):
                (repo / ignored).mkdir()
                (repo / ignored / "artifact.txt").write_text("ignored\n", encoding="ascii")
            destination = root / "snapshot"
            result = copy_source_snapshot(repo, destination)
            self.assertEqual(result["file_count"], 1)
            self.assertTrue((destination / "src/main.cpp").is_file())
            self.assertFalse((destination / "results").exists())

    def test_preserved_incomplete_root_is_forbidden(self) -> None:
        for name in (
            "preliminary-ramag-vs-mummer4-20260830T045046Z",
            "preliminary-mam16-single-20260830T083137Z",
        ):
            with self.subTest(name=name), tempfile.TemporaryDirectory() as temporary:
                base = Path(temporary)
                package = base / "package"
                (package / "results").mkdir(parents=True)
                args = Namespace(
                    package=package,
                    result_root=package / "results" / name,
                    ramag_binary=Path("/bin/true"),
                    repo=base,
                    mummer_root=base,
                    sample_interval=1.0,
                )
                with self.assertRaisesRegex(ExperimentError, "must not be reused"):
                    create_result_root(args)

    def test_full_does_not_start_ramag_without_an_accepted_mummer_run(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "SMOKE_COMPLETE").write_text("success\n", encoding="ascii")
            config = {
                "reference": str(root / "reference.fa"),
                "query": str(root / "query.fa"),
            }
            module = sys.modules["run_human_chimp_benchmark"]
            with (
                mock.patch.object(module, "require_fixed_input_integrity"),
                mock.patch.object(
                    module, "expected_mam_task_shape", return_value=FULL_MAM_TASK_SHAPE
                ),
                mock.patch.object(module, "record_full_launch_gate"),
                mock.patch.object(
                    module,
                    "run_alignment_attempt",
                    side_effect=ExperimentError("MUMmer4 failed"),
                ) as alignment,
                self.assertRaisesRegex(ExperimentError, "MUMmer4 failed"),
            ):
                run_full(root, config)
            self.assertEqual(alignment.call_count, 1)
            self.assertEqual(alignment.call_args.args[1], "mummer4")
            self.assertFalse((root / "runs/ramag").exists())
            self.assertFalse((root / "MUMMER_BASELINE_ACCEPTED").exists())

    def test_mummer_accepted_metrics_freeze_the_same_root_speed_budget(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            attempt = root / "runs/mummer4/attempt-001"
            timing = attempt / "timing"
            timing.mkdir(parents=True)
            metrics = {
                "schema": "ramag.command-metrics.v3",
                "status": "success",
                "exit_code": 0,
                "gnu_time_parse_error": None,
                "gnu_time": {
                    "elapsed_seconds": 100.0,
                    "user_seconds": 10.0,
                    "system_seconds": 1.0,
                    "maximum_resident_set_kbytes": 1024,
                    "exit_status": 0,
                },
                "sample_count": 2,
                "observed_max_total_threads": 2,
                "observed_max_single_process_threads": 2,
                "observed_max_aligner_tree_threads": 1,
                "observed_max_aligner_single_process_threads": 1,
                "observed_process_tree_peak_rss_bytes": 4096,
                "environment_overrides": OPENMP_ENVIRONMENT,
                **consistent_clock_fields(100.0),
            }
            write_json(timing / "metrics.json", metrics)
            measurements = metrics["clock_measurements"]
            write_resource_clocks(
                timing / "resources.clock.tsv",
                [
                    measurements["monotonic_started_ns"],
                    measurements["monotonic_started_ns"] + 1_000_000_000,
                ],
                [
                    measurements["realtime_started_ns"],
                    measurements["realtime_started_ns"] + 1_000_000_000,
                ],
            )
            config = {"openmp_environment": OPENMP_ENVIRONMENT}
            module = sys.modules["run_human_chimp_benchmark"]
            budget = module.freeze_mummer_speed_budget(root, attempt, config)
            self.assertEqual(budget["cutoff_seconds"], 100.0)
            self.assertEqual(budget["source_attempt"], "runs/mummer4/attempt-001")
            self.assertTrue((root / "speed-budget.json").is_file())
            self.assertTrue((root / "MUMMER_BASELINE_ACCEPTED").is_file())
            self.assertEqual(
                module.freeze_mummer_speed_budget(root, attempt, config), budget
            )

    def test_completed_speed_gate_requires_both_clocks_and_accepts_equality(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            mummer_attempt = root / "runs/mummer4/attempt-001"
            ramag_attempt = root / "runs/ramag/attempt-001"
            for attempt in (mummer_attempt, ramag_attempt):
                (attempt / "timing").mkdir(parents=True)
                write_json(attempt / "timing/metrics.json", {"fixture": True})
            write_json(root / "speed-budget.json", {"cutoff_seconds": 100.0})
            budget = {
                "result_root": str(root.resolve()),
                "source_attempt": "runs/mummer4/attempt-001",
                "cutoff_seconds": 100.0,
            }
            config = {"runtime_endpoint": "alignment raw delta"}
            module = sys.modules["run_human_chimp_benchmark"]

            def report_for(ramag_gnu: float, ramag_monotonic: float) -> dict[str, object]:
                with mock.patch.object(
                    module,
                    "_accepted_alignment_clocks",
                    side_effect=[
                        ({"schema": "ramag.command-metrics.v3"}, 100.0, 100.0),
                        (
                            {"schema": "ramag.command-metrics.v3"},
                            ramag_gnu,
                            ramag_monotonic,
                        ),
                    ],
                ):
                    return build_completed_alignment_speed_gate(
                        root,
                        mummer_attempt,
                        ramag_attempt,
                        budget,
                        config,
                    )

            equal = report_for(100.0, 100.0)
            self.assertEqual(equal["status"], "accepted")
            self.assertTrue(all(equal["conditions"].values()))
            self.assertEqual(report_for(100.001, 100.0)["status"], "rejected")
            self.assertEqual(report_for(100.0, 100.001)["status"], "rejected")

    def test_ramag_timeout_rejects_root_before_any_evaluation(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "SMOKE_COMPLETE").write_text("success\n", encoding="ascii")
            mummer_attempt = root / "runs/mummer4/attempt-001"
            ramag_attempt = root / "runs/ramag/attempt-001"
            mummer_attempt.mkdir(parents=True)
            ramag_attempt.mkdir(parents=True)
            timeout = AlignmentTimeoutError(
                "ramag",
                ramag_attempt,
                {"schema": "ramag.command-metrics.v3", "status": "timed_out"},
            )
            config = {
                "reference": str(root / "reference.fa"),
                "query": str(root / "query.fa"),
            }
            rejected = {
                "schema": "ramag.human-chimp-alignment-speed-gate.v1",
                "status": "rejected",
                "reason": "ramag_timed_out",
            }
            module = sys.modules["run_human_chimp_benchmark"]
            with (
                mock.patch.object(module, "require_fixed_input_integrity"),
                mock.patch.object(
                    module, "expected_mam_task_shape", return_value=FULL_MAM_TASK_SHAPE
                ),
                mock.patch.object(module, "record_full_launch_gate"),
                mock.patch.object(
                    module,
                    "run_alignment_attempt",
                    side_effect=[mummer_attempt, timeout],
                ) as alignment,
                mock.patch.object(
                    module,
                    "freeze_mummer_speed_budget",
                    return_value={"cutoff_seconds": 100.0},
                ),
                mock.patch.object(
                    module,
                    "build_timeout_alignment_speed_gate",
                    return_value=rejected,
                ),
                self.assertRaisesRegex(ExperimentError, "MAF conversion and F1 were not started"),
            ):
                run_full(root, config)
            self.assertEqual(alignment.call_count, 2)
            self.assertEqual(alignment.call_args_list[1].kwargs["timeout_seconds"], 100.0)
            self.assertTrue((root / "ALIGNMENT_SPEED_REJECTED").is_file())
            self.assertFalse((root / "ALIGNMENT_SPEED_ACCEPTED").exists())
            self.assertFalse((root / "FULL_RUNS_COMPLETE").exists())
            self.assertFalse((root / "evaluation").exists())

    def test_evaluation_refuses_missing_or_rejected_speed_gate(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "FULL_RUNS_COMPLETE").write_text("success\n", encoding="ascii")
            module = sys.modules["run_human_chimp_benchmark"]
            with (
                mock.patch.object(module, "run_evaluation_attempt") as evaluation,
                self.assertRaisesRegex(ExperimentError, "alignment speed gate"),
            ):
                run_evaluate(root, {})
            evaluation.assert_not_called()
            (root / "ALIGNMENT_SPEED_REJECTED").write_text(
                "rejected\n", encoding="ascii"
            )
            with (
                mock.patch.object(module, "run_evaluation_attempt") as evaluation,
                self.assertRaisesRegex(ExperimentError, "alignment speed gate"),
            ):
                run_evaluate(root, {})
            evaluation.assert_not_called()

    def test_resume_selector_rejects_timed_out_attempt(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary) / "runs/ramag"
            attempt = base / "attempt-001"
            (attempt / "timing").mkdir(parents=True)
            (base / "SELECTED_ATTEMPT").write_text(
                "attempt-001\n", encoding="ascii"
            )
            (attempt / "RUN_ACCEPTED").write_text("success\n", encoding="ascii")
            (attempt / "timing/RUN_TIMED_OUT").write_text(
                "timed_out\n", encoding="ascii"
            )
            module = sys.modules["run_human_chimp_benchmark"]
            with self.assertRaisesRegex(ExperimentError, "timed-out"):
                module.selected_attempt(base, "RUN_ACCEPTED")


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value) + "\n", encoding="utf-8")


class SummaryTests(unittest.TestCase):
    def create_bundle(self, root: Path) -> None:
        ramag_binary = root / "build/ramag"
        ramag_binary.parent.mkdir(parents=True)
        ramag_binary.write_bytes(b"fixture-ramag-binary\n")
        ramag_binary_sha = hashlib.sha256(ramag_binary.read_bytes()).hexdigest()
        reference = root / "inputs/simHuman.fa"
        query = root / "inputs/simChimp.fa"
        mummer_wrapper = Path("/fixture/mummer4/nucmer")
        write_json(
            root / "accepted-config.json",
            {
                "schema": "ramag.human-chimp-preliminary-config.v3",
                "threads": 16,
                "repetitions": 1,
                "samples": 1_000_000,
                "near": 0,
                "seeds": [20260830, 20260831, 20260832],
                "expected_sufkit_commit": "bdb67c6de5daddd8a005640de73d96549d2575f4",
                "expected_seqpro_commit": "6781cadcf81a0da53d7573444594c1484947017c",
                "ramag_selection_mode": "one-to-one",
                "openmp_environment": OPENMP_ENVIRONMENT,
                "reference": str(reference.resolve()),
                "query": str(query.resolve()),
                "input_sha256": {
                    "simHuman.fa": "1" * 64,
                    "simChimp.fa": "2" * 64,
                },
                "clock_consistency_tolerance_fraction": 0.01,
                "metrics_schema": "ramag.command-metrics.v3",
                "ramag_commit": "unknown",
                "runtime_endpoint": (
                    "FASTA read + index build + alignment + raw delta write"
                ),
                "early_stop_policy": {
                    "enabled": True,
                    "equal_is_accepted": True,
                },
                "alignment_speed_policy": {
                    "equal_is_accepted": True,
                    "observations_per_tool": 1,
                },
                "quality_policy": {
                    "mode": "report-only",
                    "minimum_f1": None,
                    "comparative_f1_gate": False,
                    "existing_ramag_diagnostic_f1": 0.9742776368,
                },
                "mummer_nucmer": str(mummer_wrapper),
                "min_match": 20,
                "max_gap": 90,
                "diag_diff": 5,
                "diag_factor": 0.12,
                "min_cluster": 65,
                "break_length": 200,
                "max_dp_cells": 4_000_000,
                "binary_snapshot": {
                    "source_sha256": "b" * 64,
                    "snapshot": str(ramag_binary.resolve()),
                    "snapshot_sha256": ramag_binary_sha,
                },
            },
        )
        provenance = root / "provenance/captures/attempt-001"
        provenance.mkdir(parents=True)
        (provenance.parent / "SELECTED_ATTEMPT").write_text(
            "attempt-001\n", encoding="ascii"
        )
        (provenance / "PROVENANCE_ACCEPTED").write_text("success\n", encoding="ascii")
        write_json(
            provenance / "summary.json",
            {
                "failed_commands": {},
                "mummer4": {
                    "wrapper": {"sha256": "a" * 64},
                    "elf": {"sha256": "c" * 64},
                    "library": {"sha256": "d" * 64},
                },
                "mummer4_current_tree": {"tree_sha256": "e" * 64},
                "source": {"tree_sha256": "f" * 64},
            },
        )
        for tool in ("mummer4", "ramag"):
            run_base = root / "runs" / tool
            run_attempt = run_base / "attempt-001"
            run_attempt.mkdir(parents=True)
            (run_base / "SELECTED_ATTEMPT").write_text("attempt-001\n", encoding="ascii")
            (run_attempt / "RUN_ACCEPTED").write_text("success\n", encoding="ascii")
            metrics = {
                "schema": "ramag.command-metrics.v3",
                "status": "success",
                "exit_code": 0,
                "gnu_time_parse_error": None,
                "sample_count": 2,
                "observed_max_total_threads": 4,
                "observed_max_single_process_threads": 4,
                "observed_max_aligner_tree_threads": 3,
                "observed_max_aligner_single_process_threads": 3,
                "observed_process_tree_peak_rss_bytes": 2048,
                "binary_path": str(
                    ramag_binary.resolve() if tool == "ramag" else mummer_wrapper
                ),
                "binary_sha256": ramag_binary_sha if tool == "ramag" else "a" * 64,
                "environment_overrides": OPENMP_ENVIRONMENT,
                "gnu_time": {
                    "elapsed_seconds": 2.0,
                    "user_seconds": 3.0,
                    "system_seconds": 0.5,
                    "cpu_percent": 175.0,
                    "maximum_resident_set_kbytes": 1000,
                    "exit_status": 0,
                },
                **consistent_clock_fields(2.0),
            }
            write_json(
                run_attempt / "timing/metrics.json",
                metrics,
            )
            monotonic_started = metrics["clock_measurements"]["monotonic_started_ns"]
            realtime_started = metrics["clock_measurements"]["realtime_started_ns"]
            write_resource_clocks(
                run_attempt / "timing/resources.clock.tsv",
                [monotonic_started, monotonic_started + 1_000_000_000],
                [realtime_started, realtime_started + 1_000_000_000],
            )
            clock_report = assess_clock_consistency(
                metrics, run_attempt / "timing/resources.clock.tsv"
            )
            write_json(run_attempt / "clock-consistency.json", clock_report)
            (run_attempt / "CLOCK_ACCEPTED").write_text(
                "success\n", encoding="ascii"
            )
            warm_finished = monotonic_started - 100_000_000
            file_identities = {
                "reference": {"device": 1, "inode": 10, "bytes": 100, "mtime_ns": 1},
                "query": {"device": 1, "inode": 11, "bytes": 100, "mtime_ns": 1},
            }
            warm_inputs = []
            for order, (role, path, digest) in enumerate(
                (
                    ("reference", reference, "1" * 64),
                    ("query", query, "2" * 64),
                ),
                start=1,
            ):
                identity = file_identities[role]
                warm_inputs.append(
                    {
                        "role": role,
                        "path": str(path.resolve()),
                        "filename": path.name,
                        "read_order": order,
                        "read_method": "single-open sequential 8-MiB chunks through EOF",
                        "expected_sha256": digest,
                        "observed_sha256": digest,
                        "bytes_read": identity["bytes"],
                        "file_identity_before": identity,
                        "file_identity_after": identity,
                        "failure": None,
                        "status": "success",
                    }
                )
            warm_cache = {
                "schema": "ramag.human-chimp-per-tool-warm-cache.v1",
                "tool": tool,
                "read_order": ["reference", "query"],
                "inputs": warm_inputs,
                "finished_monotonic_ns": warm_finished,
                "monotonic_elapsed_seconds": 0.05,
                "failures": [],
                "status": "success",
            }
            write_json(run_attempt / "warm-cache.json", warm_cache)
            (run_attempt / "WARM_CACHE_COMPLETE").write_text(
                "success\n", encoding="ascii"
            )
            write_json(
                run_attempt / "warm-cache-launch.json",
                {
                    "schema": "ramag.human-chimp-warm-cache-launch-relation.v1",
                    "warm_cache_finished_monotonic_ns": warm_finished,
                    "metrics_runner_started_monotonic_ns": monotonic_started,
                    "warm_cache_to_metrics_runner_seconds": 0.1,
                    "status": "passed",
                },
            )
            (run_attempt / "timing/RUN_COMPLETE").write_text(
                "success\n", encoding="ascii"
            )
            delta_name = (
                "mummer4.mumreference.raw.delta"
                if tool == "mummer4"
                else "ramag.mumreference.raw.delta"
            )
            (run_attempt / delta_name).write_text("delta\n", encoding="ascii")
            delta = run_attempt / delta_name
            write_json(
                run_attempt / "format-validation.json",
                {
                    "status": "success",
                    "bytes": delta.stat().st_size,
                    "sha256": hashlib.sha256(delta.read_bytes()).hexdigest(),
                    "show_coords_exit_code": 0,
                    "delta_filter_exit_code": 0,
                },
            )
            if tool == "ramag":
                write_json(
                    run_attempt / "ramag.mumreference.raw.manifest.json",
                    {
                        "status": "success",
                        "exit_code": 0,
                        "dependencies": {
                            "sufkit_commit": "bdb67c6de5daddd8a005640de73d96549d2575f4",
                            "seqpro_commit": "6781cadcf81a0da53d7573444594c1484947017c",
                        },
                        "effective_config": {
                            "reference": str(reference.resolve()),
                            "query": str(query.resolve()),
                            "threads": 16,
                            "formats": "delta",
                            "seed_mode": "mumreference",
                            "selection_mode": "one-to-one",
                            "min_match": 20,
                            "max_gap": 90,
                            "diag_diff": 5,
                            "diag_factor": 0.12,
                            "min_cluster": 65,
                            "break_length": 200,
                            "max_dp_cells": 4_000_000,
                        },
                        "actual_routes": {
                            "index": "divsufsort32",
                            "seed": "sufkit-full-sa:mumreference+"
                            "tiled-mam-boundary-mem-v1+"
                            "openmp-dynamic-stable-tile-boundary-tasks",
                            "chain": "sparse-exact-edge-components-v1",
                        },
                        "chaining": {"route": "sparse-exact-edge-components-v1"},
                        "adapter_provenance": {
                            "sufkit.backend": "divsufsort32",
                            "sufkit.commit": "bdb67c6de5daddd8a005640de73d96549d2575f4",
                            "sufkit.source_state": "exact-commit-clean-at-configure",
                        },
                        "build": {
                            "openmp_enabled": True,
                            "sufkit_divsufsort_openmp": True,
                        },
                        "openmp": {
                            "enabled": True,
                            "runtime": "gomp",
                            "runtime_max_threads": 16,
                            "configured_requested_threads": 16,
                            "dynamic": False,
                            "sufkit_divsufsort_openmp": True,
                        },
                        "threading": {
                            "seed_requested_threads": 16,
                            "seed_scheduled_threads": 16,
                            "seed_worker_threads": 16,
                            "seed_task_count": 176,
                            "seed_tasks_completed": 176,
                            "seed_parallel_route": "openmp-dynamic-stable-tile-boundary-tasks",
                            "chaining_requested_threads": 16,
                            "chaining_worker_threads": 16,
                            "extension_worker_threads": 16,
                            "seed_chain_extension_parallel": True,
                        },
                        "mam_tiling": {
                            "worker_cap": 16,
                            "tile_bases": 4 * 1024 * 1024,
                            "oriented_queries": 8,
                            "tile_tasks": 92,
                            "tile_tasks_completed": 92,
                            "short_query_tasks": 0,
                            "boundary_tasks": 84,
                            "boundary_tasks_completed": 84,
                            "tile_raw_mams": 120,
                            "tile_globally_maximal_mams": 100,
                            "boundary_raw_mems": 30,
                            "boundary_patterns": 20,
                            "boundary_reference_unique_patterns": 5,
                            "boundary_recovered_mams": 5,
                            "workspace_baseline_bytes": 1024,
                            "workspace_peak_bytes": 1536,
                            "resource_limits": {
                                "boundary_mem_occurrences": 100,
                                "workspace_bytes": 2048,
                            },
                        },
                        "counts": {
                            "mem_seeds": 30,
                            "mam_seeds": 105,
                            "selected_seeds": 100,
                        },
                    },
                )
                (run_attempt / "ramag.mumreference.raw.complete").write_text(
                    "success\n", encoding="ascii"
                )

            evaluation_base = root / "evaluation" / tool
            evaluation = evaluation_base / "attempt-001"
            (evaluation / "score").mkdir(parents=True)
            (evaluation_base / "SELECTED_ATTEMPT").write_text(
                "attempt-001\n", encoding="ascii"
            )
            (evaluation / "EVALUATION_ACCEPTED").write_text("success\n", encoding="ascii")
            (evaluation / "score/EVALUATION_COMPLETE").write_text(
                "success\n", encoding="ascii"
            )
            (evaluation / "prediction.canonical.maf").write_text(
                "##maf version=1\n", encoding="ascii"
            )
            (evaluation / "score/summary.tsv").write_text(
                "seed\tprecision\trecall\tf1\n", encoding="ascii"
            )
            aggregate = {"median": 0.9, "min": 0.89, "max": 0.91}
            raw = evaluation / "score/raw"
            raw.mkdir()
            runs = []
            for seed, score_value in zip(
                (20260830, 20260831, 20260832),
                (0.89, 0.9, 0.91),
                strict=True,
            ):
                xml = raw / f"mafComparator.seed-{seed}.xml"
                xml.write_text(
                    '<alignmentComparisons seed="{seed}" numberOfSamples="1000000" '
                    'near="0" numberOfPairsInMaf1="100" '
                    'numberOfPairsInMaf2="100"/>\n'.format(seed=seed),
                    encoding="ascii",
                )
                runs.append(
                    {
                        "xml": str(xml.resolve()),
                        "seed": seed,
                        "requested_samples": 1_000_000,
                        "near": 0,
                        "tp_truth_to_prediction": 89,
                        "false_negative": 11,
                        "tp_prediction_to_truth": 89,
                        "false_positive": 11,
                        "precision": score_value,
                        "recall": score_value,
                        "f1": score_value,
                    }
                )
            write_json(
                evaluation / "score/summary.json",
                {
                    "precision": aggregate,
                    "recall": aggregate,
                    "f1": aggregate,
                    "prediction_pair_count": 100,
                    "truth_profile": "all-homology",
                    "run_count": 3,
                    "runs": runs,
                    "seeds": [20260830, 20260831, 20260832],
                    "samples": 1_000_000,
                    "near": 0,
                },
            )
            write_json(
                evaluation / "score/prediction-validation.json",
                {
                    "total_pair_count": 100,
                    "blocks": 1,
                    "forward_rows": 1,
                    "reverse_rows": 1,
                    "coverage": {
                        "simHuman.chrA": {
                            "covered_bases": 5,
                            "sequence_length": 10,
                        },
                        "simChimp.chrA": {
                            "covered_bases": 4,
                            "sequence_length": 10,
                        },
                    },
                },
            )
            write_json(
                evaluation / "prediction.validation.json",
                {
                    "total_pair_count": 100,
                    "status": "success",
                },
            )
            for timing_name in (
                "normalize-prediction",
                "canonicalize-sources",
                "validate-canonical-maf",
                "evaluate",
            ):
                write_json(
                    evaluation / f"{timing_name}.timing.json",
                    {
                        "schema": "ramag.benchmark-postprocess-timing.v1",
                        "name": timing_name,
                        "clock": "CLOCK_MONOTONIC",
                        "monotonic_elapsed_seconds": 0.25,
                        "status": "success",
                    },
                )
            write_json(
                evaluation / "canonicalization.json",
                {
                    "status": "success",
                    "non_s_lines_unchanged": True,
                    "s_line_payload_unchanged": True,
                    "canonical_sources_seen": ["simHuman.chrA", "simChimp.chrA"],
                    "canonical_whitelist": list(CANONICAL_SOURCES),
                },
            )

        # Keep the general comparison fixture aligned with the independently
        # revalidated final-manifest policy: full Human--Chimp uses CaPS,
        # chrD smoke uses divsufsort, and the full run carries fixed external
        # FAI plus 2/2 input-OpenMP evidence.
        config_path = root / "accepted-config.json"
        config = json.loads(config_path.read_text(encoding="utf-8"))
        full_attempt = root / "runs/ramag/attempt-001"
        manifest_path = full_attempt / "ramag.mumreference.raw.manifest.json"
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        routes = manifest["actual_routes"]
        routes.update(
            {
                "index": "sufkit-full-sa:caps32:sampling=1:acceleration=suffix-link",
                "input": "seqpro-external-fai-mmap+caller-buffer-copy",
                "input_parallel": "openmp-sections-reference-query-v1",
            }
        )
        adapter = manifest["adapter_provenance"]
        adapter.update(
            {
                "sufkit.backend": "caps32",
                "sufkit.build.caps_min_reference_bases": str(64 * 1024 * 1024),
                "sufkit.reference_bases": "185145446",
            }
        )
        manifest["threading"].update(
            {
                "input_requested_workers": 2,
                "input_actual_workers": 2,
                "input_parallel_route": "openmp-sections-reference-query-v1",
            }
        )
        fai_hashes: dict[str, str] = {}
        for role, fasta_key, content in (
            ("reference", "reference", b"simHuman.chrA\t10\t15\t10\t11\n"),
            ("query", "query", b"simChimp.chrA\t10\t15\t10\t11\n"),
        ):
            fasta = Path(str(config[fasta_key]))
            source_fai = Path(str(fasta) + ".fai")
            source_fai.parent.mkdir(parents=True, exist_ok=True)
            source_fai.write_bytes(content)
            digest = hashlib.sha256(content).hexdigest()
            fai_hashes[source_fai.name] = digest
            active_fai = full_attempt / "work/runs/final/seqpro" / f"{role}.fai"
            active_fai.parent.mkdir(parents=True, exist_ok=True)
            active_fai.write_bytes(content)
            prefix = f"seqpro.{role}."
            adapter.update(
                {
                    prefix + "fai": str(active_fai.resolve()),
                    prefix + "source_fai": str(source_fai.resolve()),
                    prefix + "source_fai_bytes": str(source_fai.stat().st_size),
                    prefix + "source_fai_sha256": digest,
                    prefix + "source_fai.copy_status": "copied",
                    prefix + "external_fai.adoption_status": "adopted",
                    prefix + "build_action": "reused",
                    prefix + "index_origin": "external-standard-fai",
                    prefix + "verification": "structure-validated",
                    prefix + "metadata": "",
                }
            )
        config["input_fai_sha256"] = fai_hashes
        write_json(config_path, config)
        write_json(manifest_path, manifest)

        smoke_base = root / "smoke/ramag"
        smoke_attempt = smoke_base / "attempt-001"
        smoke_attempt.mkdir(parents=True)
        (smoke_base / "SELECTED_ATTEMPT").write_text(
            "attempt-001\n", encoding="ascii"
        )
        (smoke_attempt / "RUN_ACCEPTED").write_text("success\n", encoding="ascii")
        write_json(
            smoke_attempt / "ramag.mumreference.raw.manifest.json",
            {
                "status": "success",
                "exit_code": 0,
                "actual_routes": {
                    "index": (
                        "sufkit-full-sa:divsufsort32:sampling=1:"
                        "acceleration=suffix-link"
                    )
                },
                "adapter_provenance": {
                    "sufkit.backend": "divsufsort32",
                    "sufkit.build.caps_min_reference_bases": str(64 * 1024 * 1024),
                    "sufkit.reference_bases": "10572275",
                },
            },
        )

        mummer_attempt = root / "runs/mummer4/attempt-001"
        ramag_attempt = root / "runs/ramag/attempt-001"
        mummer_metrics_path = mummer_attempt / "timing/metrics.json"
        ramag_metrics_path = ramag_attempt / "timing/metrics.json"
        speed_budget = {
            "schema": "ramag.human-chimp-speed-budget.v1",
            "status": "accepted",
            "result_root": str(root.resolve()),
            "source_tool": "mummer4",
            "source_attempt": "runs/mummer4/attempt-001",
            "source_metrics": "runs/mummer4/attempt-001/timing/metrics.json",
            "source_metrics_sha256": hashlib.sha256(
                mummer_metrics_path.read_bytes()
            ).hexdigest(),
            "source_metrics_schema": "ramag.command-metrics.v3",
            "source_gnu_elapsed_seconds": 2.0,
            "cutoff_clock": "CLOCK_MONOTONIC",
            "cutoff_seconds": 2.0,
            "timeout_condition": (
                "ramag_monotonic_elapsed_seconds > cutoff_seconds"
            ),
            "equal_is_accepted": True,
            "observation_policy": (
                "same-root single accepted MUMmer4 observation"
            ),
        }
        write_json(root / "speed-budget.json", speed_budget)
        speed_gate = {
            "schema": "ramag.human-chimp-alignment-speed-gate.v1",
            "status": "accepted",
            "reason": "dual_clock_passed",
            "scope": "preliminary single observation",
            "endpoint": config["runtime_endpoint"],
            "equal_is_accepted": True,
            "hard_deadline_clock": "CLOCK_MONOTONIC",
            "hard_deadline_seconds": 2.0,
            "budget_sha256": hashlib.sha256(
                (root / "speed-budget.json").read_bytes()
            ).hexdigest(),
            "conditions": {
                "ramag_completed_before_hard_deadline": True,
                "ramag_gnu_not_slower": True,
                "ramag_monotonic_not_slower": True,
            },
            "mummer4": {
                "attempt": "runs/mummer4/attempt-001",
                "metrics_sha256": hashlib.sha256(
                    mummer_metrics_path.read_bytes()
                ).hexdigest(),
                "metrics_schema": "ramag.command-metrics.v3",
                "gnu_elapsed_seconds": 2.0,
                "monotonic_elapsed_seconds": 2.0,
            },
            "ramag": {
                "attempt": "runs/ramag/attempt-001",
                "metrics_sha256": hashlib.sha256(
                    ramag_metrics_path.read_bytes()
                ).hexdigest(),
                "metrics_schema": "ramag.command-metrics.v3",
                "status": "success",
                "gnu_elapsed_seconds": 2.0,
                "monotonic_elapsed_seconds": 2.0,
            },
        }
        write_json(root / "alignment-speed-gate.json", speed_gate)
        for marker in (
            "MUMMER_BASELINE_ACCEPTED",
            "ALIGNMENT_SPEED_ACCEPTED",
            "FULL_RUNS_COMPLETE",
            "EVALUATIONS_COMPLETE",
        ):
            (root / marker).write_text("success\n", encoding="ascii")

    def test_complete_bundle_publishes_marker_last(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.create_bundle(root)
            result = summarize(root)
            self.assertEqual(result["status"], "complete")
            self.assertEqual(result["performance_gate"]["status"], "passed")
            self.assertEqual(result["clock_robustness_gate"]["status"], "passed")
            self.assertTrue((root / "comparison.tsv").is_file())
            self.assertTrue((root / "comparison.json").is_file())
            self.assertTrue((root / "comparison.md").is_file())
            self.assertTrue((root / "COMPARISON_COMPLETE").is_file())
            header = (root / "comparison.tsv").read_text(encoding="utf-8").splitlines()[0]
            self.assertIn("source_tree_sha256", header)
            self.assertIn("ramag_elf_sha256", header)
            self.assertIn("mummer_elf_sha256", header)
            self.assertIn("openmp_runtime", header)
            rows = result["tools"]
            ramag = next(row for row in rows if row["tool_id"] == "ramag")
            self.assertEqual(ramag["source_tree_sha256"], "f" * 64)
            self.assertTrue(ramag["sufkit_divsufsort_openmp"])

    def test_tampered_slower_ramag_invalidates_the_frozen_speed_gate(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.create_bundle(root)
            ramag_metrics = (
                root / "runs/ramag/attempt-001/timing/metrics.json"
            )
            metrics = json.loads(ramag_metrics.read_text(encoding="utf-8"))
            metrics["gnu_time"]["elapsed_seconds"] = 2.001
            ramag_metrics.write_text(json.dumps(metrics) + "\n", encoding="utf-8")
            write_json(
                root / "runs/ramag/attempt-001/clock-consistency.json",
                assess_clock_consistency(
                    metrics,
                    root / "runs/ramag/attempt-001/timing/resources.clock.tsv",
                ),
            )
            with self.assertRaisesRegex(SummaryError, "alignment speed gate"):
                summarize(root)
            self.assertFalse((root / "COMPARISON_COMPLETE").exists())
            self.assertFalse((root / "comparison.json").exists())

    def test_summary_rejects_clock_quarantine_before_headline(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.create_bundle(root)
            metrics_path = root / "runs/ramag/attempt-001/timing/metrics.json"
            metrics = json.loads(metrics_path.read_text(encoding="utf-8"))
            clocks_path = root / "runs/ramag/attempt-001/timing/resources.clock.tsv"
            monotonic_started = metrics["clock_measurements"]["monotonic_started_ns"]
            realtime_started = metrics["clock_measurements"]["realtime_started_ns"]
            write_resource_clocks(
                clocks_path,
                [monotonic_started, monotonic_started + 1_000_000_000],
                [realtime_started, realtime_started - 1_000_000_000],
            )
            write_json(
                root / "runs/ramag/attempt-001/clock-consistency.json",
                assess_clock_consistency(metrics, clocks_path),
            )
            with self.assertRaisesRegex(SummaryError, "clock-quarantined"):
                summarize(root)
            self.assertFalse((root / "comparison.json").exists())
            self.assertFalse((root / "COMPARISON_COMPLETE").exists())

    def test_cross_tool_clock_distortion_is_a_nonblocking_warning(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.create_bundle(root)
            metrics_path = root / "runs/ramag/attempt-001/timing/metrics.json"
            metrics = json.loads(metrics_path.read_text(encoding="utf-8"))
            measurements = metrics["clock_measurements"]
            realtime_started = measurements["realtime_started_ns"]
            realtime_finished = realtime_started + 2_100_000_000
            measurements["realtime_finished_ns"] = realtime_finished
            measurements["realtime_elapsed_seconds"] = 2.1
            measurements["finished_utc"] = utc_from_ns(realtime_finished)
            metrics["finished_utc"] = measurements["finished_utc"]
            metrics_path.write_text(json.dumps(metrics) + "\n", encoding="utf-8")
            speed_gate_path = root / "alignment-speed-gate.json"
            speed_gate = json.loads(speed_gate_path.read_text(encoding="utf-8"))
            speed_gate["ramag"]["metrics_sha256"] = hashlib.sha256(
                metrics_path.read_bytes()
            ).hexdigest()
            write_json(speed_gate_path, speed_gate)
            report = assess_clock_consistency(
                metrics,
                root / "runs/ramag/attempt-001/timing/resources.clock.tsv",
            )
            self.assertEqual(report["status"], "warning")
            write_json(
                root / "runs/ramag/attempt-001/clock-consistency.json", report
            )
            (root / "runs/ramag/attempt-001/CLOCK_WARNING").write_text(
                "warning\n", encoding="ascii"
            )
            result = summarize(root)
            self.assertEqual(result["status"], "complete")
            self.assertEqual(result["performance_gate"]["status"], "passed")
            self.assertEqual(result["clock_robustness_gate"]["status"], "passed")
            self.assertEqual(
                result["clock_robustness_gate"]["diagnostic_status"], "warning"
            )
            self.assertTrue((root / "COMPARISON_COMPLETE").is_file())

    def test_summary_rejects_incomplete_full_tiled_mam_shape(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.create_bundle(root)
            path = root / "runs/ramag/attempt-001/ramag.mumreference.raw.manifest.json"
            manifest = json.loads(path.read_text(encoding="utf-8"))
            manifest["mam_tiling"]["tile_tasks"] = 91
            manifest["mam_tiling"]["tile_tasks_completed"] = 91
            manifest["threading"]["seed_task_count"] = 175
            manifest["threading"]["seed_tasks_completed"] = 175
            path.write_text(json.dumps(manifest) + "\n", encoding="utf-8")
            with self.assertRaisesRegex(SummaryError, "OpenMP contract"):
                summarize(root)

    def test_summary_rejects_impossible_tiled_mam_statistics(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.create_bundle(root)
            path = root / "runs/ramag/attempt-001/ramag.mumreference.raw.manifest.json"
            manifest = json.loads(path.read_text(encoding="utf-8"))
            manifest["mam_tiling"]["boundary_recovered_mams"] = 999
            path.write_text(json.dumps(manifest) + "\n", encoding="utf-8")
            with self.assertRaisesRegex(SummaryError, "statistical invariants"):
                summarize(root)

    def test_summary_rejects_impossible_thread_metric_hierarchy(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.create_bundle(root)
            path = root / "runs/ramag/attempt-001/timing/metrics.json"
            metrics = json.loads(path.read_text(encoding="utf-8"))
            metrics["observed_max_total_threads"] = 2
            metrics["observed_max_aligner_tree_threads"] = 3
            path.write_text(json.dumps(metrics) + "\n", encoding="utf-8")
            with self.assertRaisesRegex(SummaryError, "hierarchy"):
                summarize(root)

    def test_summary_rejects_missing_fixed_seed_xml(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.create_bundle(root)
            xml = (
                root
                / "evaluation/ramag/attempt-001/score/raw/mafComparator.seed-20260832.xml"
            )
            xml.unlink()
            with self.assertRaisesRegex(SummaryError, "fixed-seed XML"):
                summarize(root)

    def test_summary_rejects_changed_scoring_parameters(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.create_bundle(root)
            path = root / "evaluation/ramag/attempt-001/score/summary.json"
            summary = json.loads(path.read_text(encoding="utf-8"))
            summary["seeds"] = [9]
            path.write_text(json.dumps(summary) + "\n", encoding="utf-8")
            with self.assertRaisesRegex(SummaryError, "fixed scoring contract"):
                summarize(root)

    def test_summary_rejects_missing_ramag_complete_marker(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.create_bundle(root)
            (root / "runs/ramag/attempt-001/ramag.mumreference.raw.complete").unlink()
            with self.assertRaisesRegex(SummaryError, "atomic complete"):
                summarize(root)

    def test_summary_rejects_runtime_binary_hash_mismatch(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.create_bundle(root)
            path = root / "runs/ramag/attempt-001/timing/metrics.json"
            metrics = json.loads(path.read_text(encoding="utf-8"))
            metrics["binary_sha256"] = "0" * 64
            path.write_text(json.dumps(metrics) + "\n", encoding="utf-8")
            with self.assertRaisesRegex(SummaryError, "accepted snapshot"):
                summarize(root)

    def test_incomplete_bundle_does_not_publish_marker(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            root.mkdir(exist_ok=True)
            write_json(
                root / "accepted-config.json",
                {
                    "threads": 16,
                    "repetitions": 1,
                    "samples": 1,
                    "near": 0,
                    "seeds": [1],
                    "expected_sufkit_commit": "x",
                },
            )
            with self.assertRaises(SummaryError):
                summarize(root)
            self.assertFalse((root / "COMPARISON_COMPLETE").exists())


if __name__ == "__main__":
    unittest.main()

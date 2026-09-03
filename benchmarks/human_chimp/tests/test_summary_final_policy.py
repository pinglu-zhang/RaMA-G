#!/usr/bin/env python3

from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path


TEST_DIR = Path(__file__).resolve().parent
MODULE_DIR = TEST_DIR.parent
sys.path.insert(0, str(TEST_DIR))
sys.path.insert(0, str(MODULE_DIR))

import test_harness as fixture_support  # noqa: E402
from run_with_metrics import assess_clock_consistency  # noqa: E402
from summarize_comparison import (  # noqa: E402
    CAPS_MIN_REFERENCE_BASES,
    SummaryError,
    build_clock_robustness_gate,
    sha256_file,
    summarize,
    validate_final_ramag_manifest_gate,
)


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value) + "\n", encoding="utf-8")


def install_formal_manifest_evidence(root: Path) -> tuple[dict[str, object], Path]:
    config_path = root / "accepted-config.json"
    config = json.loads(config_path.read_text(encoding="utf-8"))
    config["schema"] = "ramag.human-chimp-preliminary-config.v3"

    full_attempt = root / "runs/ramag/attempt-001"
    full_manifest_path = full_attempt / "ramag.mumreference.raw.manifest.json"
    manifest = json.loads(full_manifest_path.read_text(encoding="utf-8"))
    adapter = manifest["adapter_provenance"]
    routes = manifest["actual_routes"]
    threading = manifest["threading"]
    routes.update(
        {
            "index": "sufkit-full-sa:caps32:sampling=1:acceleration=suffix-link",
            "input": "seqpro-external-fai-mmap+caller-buffer-copy",
            "input_parallel": "openmp-sections-reference-query-v1",
        }
    )
    adapter.update(
        {
            "sufkit.backend": "caps32",
            "sufkit.build.caps_min_reference_bases": str(
                CAPS_MIN_REFERENCE_BASES
            ),
            "sufkit.reference_bases": "185145446",
        }
    )
    threading.update(
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
        digest = sha256_file(source_fai)
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
    write_json(full_manifest_path, manifest)

    smoke_base = root / "smoke/ramag"
    smoke_attempt = smoke_base / "attempt-001"
    smoke_attempt.mkdir(parents=True, exist_ok=True)
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
                "sufkit.build.caps_min_reference_bases": str(
                    CAPS_MIN_REFERENCE_BASES
                ),
                "sufkit.reference_bases": "10572275",
            },
        },
    )
    return config, full_attempt


class FinalManifestGateTests(unittest.TestCase):
    def create_formal_bundle(self, root: Path) -> tuple[dict[str, object], Path]:
        fixture_support.SummaryTests().create_bundle(root)
        return install_formal_manifest_evidence(root)

    def test_revalidates_full_smoke_external_fai_and_input_openmp(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config, full_attempt = self.create_formal_bundle(root)
            evidence = validate_final_ramag_manifest_gate(
                root, config, full_attempt
            )
            self.assertEqual(evidence["status"], "passed")
            self.assertEqual(evidence["full"]["backend"], "caps32")
            self.assertEqual(evidence["smoke"]["backend"], "divsufsort32")
            self.assertEqual(
                evidence["full"]["caps_min_reference_bases"],
                CAPS_MIN_REFERENCE_BASES,
            )
            self.assertEqual(evidence["full"]["input_actual_workers"], 2)
            self.assertEqual(
                set(evidence["external_fai"]), {"reference", "query"}
            )

    def test_rejects_tampered_external_fai_copy(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config, full_attempt = self.create_formal_bundle(root)
            manifest = json.loads(
                (
                    full_attempt / "ramag.mumreference.raw.manifest.json"
                ).read_text(encoding="utf-8")
            )
            active = Path(
                manifest["adapter_provenance"]["seqpro.query.fai"]
            )
            active.write_text("tampered\n", encoding="ascii")
            with self.assertRaisesRegex(SummaryError, "query external-FAI"):
                validate_final_ramag_manifest_gate(root, config, full_attempt)

    def test_rejects_wrong_full_or_smoke_backend(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config, full_attempt = self.create_formal_bundle(root)
            manifest_path = full_attempt / "ramag.mumreference.raw.manifest.json"
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            manifest["adapter_provenance"]["sufkit.backend"] = "divsufsort32"
            write_json(manifest_path, manifest)
            with self.assertRaisesRegex(SummaryError, "full caps32"):
                validate_final_ramag_manifest_gate(root, config, full_attempt)

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config, full_attempt = self.create_formal_bundle(root)
            smoke_path = (
                root
                / "smoke/ramag/attempt-001/ramag.mumreference.raw.manifest.json"
            )
            smoke = json.loads(smoke_path.read_text(encoding="utf-8"))
            smoke["adapter_provenance"]["sufkit.backend"] = "caps32"
            write_json(smoke_path, smoke)
            with self.assertRaisesRegex(SummaryError, "smoke divsufsort32"):
                validate_final_ramag_manifest_gate(root, config, full_attempt)


class FrozenClockPolicyTests(unittest.TestCase):
    @staticmethod
    def row(
        gnu: float,
        monotonic: float,
        realtime: float,
        jumps: int = 0,
    ) -> dict[str, object]:
        return {
            "wall_seconds": gnu,
            "monotonic_wall_seconds": monotonic,
            "realtime_utc_wall_seconds": realtime,
            "clock_consistency_status": "warning" if jumps else "passed",
            "resource_realtime_jump_count": jumps,
        }

    def test_distortion_and_positive_jump_are_nonblocking_warnings(self) -> None:
        gate = build_clock_robustness_gate(
            {
                "mummer4": self.row(100.0, 100.0, 100.0),
                "ramag": self.row(90.0, 90.0, 120.0, jumps=2),
            }
        )
        self.assertEqual(gate["schema"], "ramag.human-chimp-clock-robustness-gate.v2")
        self.assertEqual(gate["status"], "passed")
        self.assertEqual(gate["diagnostic_status"], "warning")
        self.assertFalse(
            gate["diagnostic_conditions"][
                "realtime_over_monotonic_distortion_relative_difference_at_most_1pct"
            ]
        )
        self.assertTrue(gate["warnings"])

    def test_monotonic_winner_remains_a_hard_gate(self) -> None:
        gate = build_clock_robustness_gate(
            {
                "mummer4": self.row(100.0, 100.0, 100.0),
                "ramag": self.row(90.0, 101.0, 101.0),
            }
        )
        self.assertEqual(gate["status"], "failed")
        self.assertFalse(
            gate["hard_conditions"]["monotonic_ramag_not_slower"]
        )

    def test_summary_publishes_marker_with_diagnostic_warning(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fixture_support.SummaryTests().create_bundle(root)
            install_formal_manifest_evidence(root)
            metrics_path = root / "runs/ramag/attempt-001/timing/metrics.json"
            metrics = json.loads(metrics_path.read_text(encoding="utf-8"))
            measurements = metrics["clock_measurements"]
            realtime_started = measurements["realtime_started_ns"]
            realtime_finished = realtime_started + 2_100_000_000
            measurements["realtime_finished_ns"] = realtime_finished
            measurements["realtime_elapsed_seconds"] = 2.1
            measurements["finished_utc"] = fixture_support.utc_from_ns(
                realtime_finished
            )
            metrics["finished_utc"] = measurements["finished_utc"]
            write_json(metrics_path, metrics)
            speed_gate_path = root / "alignment-speed-gate.json"
            speed_gate = json.loads(speed_gate_path.read_text(encoding="utf-8"))
            speed_gate["ramag"]["metrics_sha256"] = sha256_file(metrics_path)
            write_json(speed_gate_path, speed_gate)
            clock_path = root / "runs/ramag/attempt-001/timing/resources.clock.tsv"
            report = assess_clock_consistency(metrics, clock_path)
            self.assertEqual(report["status"], "warning")
            write_json(
                root / "runs/ramag/attempt-001/clock-consistency.json", report
            )
            (root / "runs/ramag/attempt-001/CLOCK_WARNING").write_text(
                "warning\n", encoding="ascii"
            )

            result = summarize(root)
            self.assertEqual(result["schema"], "ramag.human-chimp-preliminary-comparison.v4")
            self.assertEqual(result["status"], "complete")
            self.assertEqual(result["clock_robustness_gate"]["status"], "passed")
            self.assertEqual(
                result["clock_robustness_gate"]["diagnostic_status"], "warning"
            )
            self.assertTrue((root / "COMPARISON_COMPLETE").is_file())
            markdown = (root / "comparison.md").read_text(encoding="utf-8")
            self.assertIn("full=`caps32`", markdown)
            self.assertIn("smoke=`divsufsort32`", markdown)
            self.assertIn("仅作为 warning/诊断", markdown)


if __name__ == "__main__":
    unittest.main()

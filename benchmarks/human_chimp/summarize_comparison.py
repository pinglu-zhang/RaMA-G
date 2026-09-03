#!/usr/bin/env python3
"""Create the auditable RaMA-G/MUMmer4 Human--Chimp comparison tables."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

from run_with_metrics import (
    CLOCK_CONSISTENCY_TOLERANCE_FRACTION,
    assess_clock_consistency,
)


class SummaryError(RuntimeError):
    """An incomplete or inconsistent result bundle."""


TOOLS = ("mummer4", "ramag")
MAM_TILE_BASES = 4 * 1024 * 1024
CAPS_MIN_REFERENCE_BASES = 64 * 1024 * 1024
FULL_SA_BACKEND = "caps32"
SMOKE_SA_BACKEND = "divsufsort32"
EXTERNAL_FAI_INPUT_ROUTE = "seqpro-external-fai-mmap+caller-buffer-copy"
INPUT_PARALLEL_ROUTE = "openmp-sections-reference-query-v1"
CONFIG_SCHEMA = "ramag.human-chimp-preliminary-config.v3"
METRICS_SCHEMA = "ramag.command-metrics.v3"
EXISTING_RAMAG_DIAGNOSTIC_F1 = 0.9742776368
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


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, object]:
    if not path.is_file():
        raise SummaryError(f"missing JSON artifact: {path}")
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise SummaryError(f"cannot parse JSON artifact: {path}") from exc
    if not isinstance(value, dict):
        raise SummaryError(f"expected JSON object: {path}")
    return value


def selected_attempt(root: Path, category: str, tool: str) -> Path:
    base = root / category / tool
    selector = base / "SELECTED_ATTEMPT"
    if not selector.is_file():
        raise SummaryError(f"missing selected attempt marker: {selector}")
    name = selector.read_text(encoding="ascii").strip()
    if not name.startswith("attempt-") or "/" in name or "\\" in name:
        raise SummaryError(f"invalid selected attempt name in {selector}: {name!r}")
    attempt = base / name
    required_marker = (
        "RUN_ACCEPTED"
        if category in {"runs", "smoke"}
        else "EVALUATION_ACCEPTED"
    )
    if not (attempt / required_marker).is_file():
        raise SummaryError(f"selected attempt lacks {required_marker}: {attempt}")
    if category == "runs" and (attempt / "timing/RUN_TIMED_OUT").exists():
        raise SummaryError(f"timed-out attempt cannot be selected: {attempt}")
    return attempt


def load_successful_timing(path: Path, expected_name: str) -> dict[str, object]:
    report = load_json(path)
    elapsed = report.get("monotonic_elapsed_seconds")
    if (
        report.get("schema") != "ramag.benchmark-postprocess-timing.v1"
        or report.get("name") != expected_name
        or report.get("clock") != "CLOCK_MONOTONIC"
        or report.get("status") != "success"
        or not isinstance(elapsed, (int, float))
        or not math.isfinite(float(elapsed))
        or float(elapsed) < 0.0
    ):
        raise SummaryError(f"invalid postprocessing timing evidence: {path}")
    return report


def metric_value(metrics: dict[str, object], key: str) -> object:
    gnu_time = metrics.get("gnu_time")
    if not isinstance(gnu_time, dict):
        return None
    return gnu_time.get(key)


def total_coverage(validation: dict[str, object], species: str) -> dict[str, object]:
    coverage = validation.get("coverage")
    if not isinstance(coverage, dict):
        raise SummaryError("prediction validation lacks coverage")
    covered = 0
    length = 0
    per_contig: dict[str, object] = {}
    for source, value in coverage.items():
        if not str(source).startswith(species + "."):
            continue
        if not isinstance(value, dict):
            raise SummaryError(f"invalid coverage entry for {source}")
        observed = int(value.get("covered_bases", 0))
        source_length = int(value.get("sequence_length", 0))
        covered += observed
        length += source_length
        per_contig[str(source)] = value
    return {
        "covered_bases": covered,
        "sequence_length": length,
        "fraction": covered / length if length else None,
        "per_contig": per_contig,
    }


def load_provenance(root: Path) -> tuple[dict[str, object], Path]:
    base = root / "provenance/captures"
    selector = base / "SELECTED_ATTEMPT"
    if not selector.is_file():
        raise SummaryError("missing selected provenance capture")
    name = selector.read_text(encoding="ascii").strip()
    attempt = base / name
    if not name.startswith("attempt-") or not (attempt / "PROVENANCE_ACCEPTED").is_file():
        raise SummaryError("invalid selected provenance capture")
    summary = load_json(attempt / "summary.json")
    if summary.get("failed_commands"):
        raise SummaryError("selected provenance capture contains failed commands")
    return summary, attempt


def required_int(mapping: dict[str, object], key: str, context: str) -> int:
    try:
        return int(mapping[key])
    except (KeyError, TypeError, ValueError) as exc:
        raise SummaryError(f"{context} lacks an integer {key}") from exc


def validate_final_ramag_manifest_gate(
    root: Path,
    config: dict[str, object],
    full_run_attempt: Path,
) -> dict[str, object]:
    """Re-read final manifests and independently enforce fixed I/O/backend routes.

    The alignment-stage gate already checks this contract before selecting an
    attempt.  Final summarization intentionally repeats the checks against the
    still-present source and copied FAI files, so a stale selector, edited
    manifest, or damaged work directory cannot receive COMPARISON_COMPLETE.
    """

    manifest_path = full_run_attempt / "ramag.mumreference.raw.manifest.json"
    manifest = load_json(manifest_path)
    adapter = manifest.get("adapter_provenance")
    routes = manifest.get("actual_routes")
    threading = manifest.get("threading")
    dependencies = manifest.get("dependencies")
    effective = manifest.get("effective_config")
    if (
        manifest.get("status") != "success"
        or int(manifest.get("exit_code", 1)) != 0
        or not isinstance(adapter, dict)
        or not isinstance(routes, dict)
        or not isinstance(threading, dict)
        or not isinstance(dependencies, dict)
        or not isinstance(effective, dict)
    ):
        raise SummaryError("final RaMA-G manifest gate lacks a successful full manifest")

    reference = Path(str(config.get("reference", ""))).resolve()
    query = Path(str(config.get("query", ""))).resolve()
    threads = required_int(config, "threads", "accepted config")
    if (
        dependencies.get("sufkit_commit") != config.get("expected_sufkit_commit")
        or dependencies.get("seqpro_commit") != config.get("expected_seqpro_commit")
        or effective.get("reference") != str(reference)
        or effective.get("query") != str(query)
        or effective.get("formats") != "delta"
        or effective.get("seed_mode") != "mumreference"
        or effective.get("selection_mode") != config.get("ramag_selection_mode")
        or required_int(effective, "threads", "RaMA-G effective config") != threads
    ):
        raise SummaryError("final RaMA-G manifest gate found the wrong full-run contract")

    full_backend = adapter.get("sufkit.backend")
    full_index_route = routes.get("index")
    threshold = required_int(
        adapter,
        "sufkit.build.caps_min_reference_bases",
        "RaMA-G adapter provenance",
    )
    full_reference_bases = required_int(
        adapter, "sufkit.reference_bases", "RaMA-G adapter provenance"
    )
    if (
        full_backend != FULL_SA_BACKEND
        or not isinstance(full_index_route, str)
        or f":{FULL_SA_BACKEND}:" not in full_index_route
        or threshold != CAPS_MIN_REFERENCE_BASES
        or full_reference_bases < threshold
    ):
        raise SummaryError(
            "final RaMA-G manifest gate requires full caps32 at the fixed 64 MiB threshold"
        )

    if (
        routes.get("input") != EXTERNAL_FAI_INPUT_ROUTE
        or routes.get("input_parallel") != INPUT_PARALLEL_ROUTE
        or required_int(threading, "input_requested_workers", "RaMA-G threading") != 2
        or required_int(threading, "input_actual_workers", "RaMA-G threading") != 2
        or threading.get("input_parallel_route") != INPUT_PARALLEL_ROUTE
    ):
        raise SummaryError(
            "final RaMA-G manifest gate requires the external-FAI and 2/2 OpenMP input routes"
        )

    accepted_fai_hashes = config.get("input_fai_sha256")
    if not isinstance(accepted_fai_hashes, dict):
        raise SummaryError("accepted config lacks fixed input FAI hashes")
    fai_evidence: dict[str, object] = {}
    for role, fasta in (("reference", reference), ("query", query)):
        source_fai = Path(str(fasta) + ".fai")
        if not source_fai.is_file():
            raise SummaryError(f"final RaMA-G {role} source FAI is missing: {source_fai}")
        expected_hash = accepted_fai_hashes.get(source_fai.name)
        source_hash = sha256_file(source_fai)
        source_bytes = source_fai.stat().st_size
        prefix = f"seqpro.{role}."
        active_fai = Path(str(adapter.get(prefix + "fai", "")))
        if (
            not isinstance(expected_hash, str)
            or source_hash != expected_hash
            or adapter.get(prefix + "source_fai") != str(source_fai.resolve())
            or required_int(
                adapter, prefix + "source_fai_bytes", f"RaMA-G {role} FAI provenance"
            )
            != source_bytes
            or adapter.get(prefix + "source_fai_sha256") != expected_hash
            or adapter.get(prefix + "source_fai.copy_status") != "copied"
            or adapter.get(prefix + "external_fai.adoption_status") != "adopted"
            or adapter.get(prefix + "build_action") != "reused"
            or adapter.get(prefix + "index_origin") != "external-standard-fai"
            or adapter.get(prefix + "verification") != "structure-validated"
            or adapter.get(prefix + "metadata") != ""
            or not active_fai.is_file()
            or active_fai.resolve() == source_fai.resolve()
            or full_run_attempt.resolve() not in active_fai.resolve().parents
            or active_fai.stat().st_size != source_bytes
            or sha256_file(active_fai) != expected_hash
        ):
            raise SummaryError(
                f"final RaMA-G {role} external-FAI source/copy evidence is invalid"
            )
        fai_evidence[role] = {
            "source": str(source_fai.resolve()),
            "active_copy": str(active_fai.resolve()),
            "bytes": source_bytes,
            "sha256": expected_hash,
            "copy_status": "copied",
            "adoption_status": "adopted",
            "verification": "structure-validated",
        }

    smoke_attempt = selected_attempt(root, "smoke", "ramag")
    smoke_manifest = load_json(smoke_attempt / "ramag.mumreference.raw.manifest.json")
    smoke_adapter = smoke_manifest.get("adapter_provenance")
    smoke_routes = smoke_manifest.get("actual_routes")
    if (
        smoke_manifest.get("status") != "success"
        or int(smoke_manifest.get("exit_code", 1)) != 0
        or not isinstance(smoke_adapter, dict)
        or not isinstance(smoke_routes, dict)
    ):
        raise SummaryError("final RaMA-G manifest gate lacks a successful smoke manifest")
    smoke_backend = smoke_adapter.get("sufkit.backend")
    smoke_index_route = smoke_routes.get("index")
    smoke_threshold = required_int(
        smoke_adapter,
        "sufkit.build.caps_min_reference_bases",
        "RaMA-G smoke adapter provenance",
    )
    smoke_reference_bases = required_int(
        smoke_adapter, "sufkit.reference_bases", "RaMA-G smoke adapter provenance"
    )
    if (
        smoke_backend != SMOKE_SA_BACKEND
        or not isinstance(smoke_index_route, str)
        or f":{SMOKE_SA_BACKEND}:" not in smoke_index_route
        or smoke_threshold != CAPS_MIN_REFERENCE_BASES
        or smoke_reference_bases >= smoke_threshold
    ):
        raise SummaryError(
            "final RaMA-G manifest gate requires smoke divsufsort32 below the fixed 64 MiB threshold"
        )

    return {
        "schema": "ramag.human-chimp-final-manifest-gate.v1",
        "status": "passed",
        "manifest": str(manifest_path.relative_to(root)),
        "full": {
            "backend": full_backend,
            "index_route": full_index_route,
            "reference_bases": full_reference_bases,
            "caps_min_reference_bases": threshold,
            "input_route": routes.get("input"),
            "input_parallel_route": routes.get("input_parallel"),
            "input_requested_workers": 2,
            "input_actual_workers": 2,
        },
        "smoke": {
            "manifest": str(
                (smoke_attempt / "ramag.mumreference.raw.manifest.json").relative_to(root)
            ),
            "backend": smoke_backend,
            "index_route": smoke_index_route,
            "reference_bases": smoke_reference_bases,
            "caps_min_reference_bases": smoke_threshold,
        },
        "external_fai": fai_evidence,
    }


def validate_full_run_fairness_evidence(
    run_attempt: Path,
    config: dict[str, object],
    metrics: dict[str, object],
    tool: str,
) -> tuple[dict[str, object], dict[str, object], dict[str, object]]:
    if float(config.get("clock_consistency_tolerance_fraction", -1.0)) != (
        CLOCK_CONSISTENCY_TOLERANCE_FRACTION
    ):
        raise SummaryError("accepted config lacks the fixed 1% clock gate")
    if not (run_attempt / "WARM_CACHE_COMPLETE").is_file():
        raise SummaryError(f"selected {tool} run lacks per-tool warm-cache evidence")
    if not (run_attempt / "CLOCK_ACCEPTED").is_file() or (
        run_attempt / "CLOCK_QUARANTINED"
    ).exists():
        raise SummaryError(f"selected {tool} run is not clock-accepted")

    warm_cache = load_json(run_attempt / "warm-cache.json")
    expected_hashes = config.get("input_sha256")
    if (
        warm_cache.get("status") != "success"
        or warm_cache.get("tool") != tool
        or warm_cache.get("read_order") != ["reference", "query"]
        or not isinstance(expected_hashes, dict)
    ):
        raise SummaryError(f"selected {tool} warm-cache contract is invalid")
    entries = warm_cache.get("inputs")
    if not isinstance(entries, list) or len(entries) != 2:
        raise SummaryError(f"selected {tool} warm-cache input evidence is incomplete")
    expected_paths = {
        "reference": Path(str(config["reference"])).resolve(),
        "query": Path(str(config["query"])).resolve(),
    }
    for expected_order, (role, entry) in enumerate(
        zip(("reference", "query"), entries, strict=True), start=1
    ):
        if not isinstance(entry, dict):
            raise SummaryError(f"selected {tool} warm-cache entry is invalid")
        expected_path = expected_paths[role]
        expected_hash = expected_hashes.get(expected_path.name)
        before = entry.get("file_identity_before")
        after = entry.get("file_identity_after")
        if (
            entry.get("role") != role
            or int(entry.get("read_order", 0)) != expected_order
            or Path(str(entry.get("path", ""))).resolve() != expected_path
            or not isinstance(expected_hash, str)
            or entry.get("expected_sha256") != expected_hash
            or entry.get("observed_sha256") != expected_hash
            or entry.get("status") != "success"
            or entry.get("failure") is not None
            or not isinstance(before, dict)
            or before != after
            or int(entry.get("bytes_read", -1)) != int(before.get("bytes", -2))
            or int(entry.get("bytes_read", 0)) <= 0
        ):
            raise SummaryError(
                f"selected {tool} warm-cache {role} hash/read evidence is invalid"
            )

    launch = load_json(run_attempt / "warm-cache-launch.json")
    measurements = metrics.get("clock_measurements")
    warm_finished = warm_cache.get("finished_monotonic_ns")
    runner_started = (
        measurements.get("monotonic_started_ns")
        if isinstance(measurements, dict)
        else None
    )
    expected_gap = (
        (runner_started - warm_finished) / 1_000_000_000.0
        if isinstance(warm_finished, int)
        and isinstance(runner_started, int)
        and runner_started >= warm_finished
        else None
    )
    if (
        launch.get("status") != "passed"
        or launch.get("warm_cache_finished_monotonic_ns") != warm_finished
        or launch.get("metrics_runner_started_monotonic_ns") != runner_started
        or launch.get("warm_cache_to_metrics_runner_seconds") != expected_gap
        or expected_gap is None
    ):
        raise SummaryError(
            f"selected {tool} warm-cache evidence is not ordered before launch"
        )

    recorded_clock = load_json(run_attempt / "clock-consistency.json")
    recomputed_clock = assess_clock_consistency(
        metrics,
        run_attempt / "timing/resources.clock.tsv",
        tolerance_fraction=CLOCK_CONSISTENCY_TOLERANCE_FRACTION,
    )
    if (
        recomputed_clock.get("status") == "quarantined"
        or recorded_clock != recomputed_clock
    ):
        raise SummaryError(
            f"selected {tool} run is clock-quarantined and cannot enter headline results"
        )
    if recomputed_clock.get("status") == "warning" and not (
        run_attempt / "CLOCK_WARNING"
    ).is_file():
        raise SummaryError(f"selected {tool} warning-quality clock run lacks its marker")
    return warm_cache, launch, recomputed_clock


def build_row(
    root: Path,
    config: dict[str, object],
    provenance: dict[str, object],
    provenance_attempt: Path,
    tool: str,
) -> dict[str, object]:
    run_attempt = selected_attempt(root, "runs", tool)
    evaluation_attempt = selected_attempt(root, "evaluation", tool)
    metrics = load_json(run_attempt / "timing/metrics.json")
    summary = load_json(evaluation_attempt / "score/summary.json")
    validation = load_json(evaluation_attempt / "score/prediction-validation.json")
    canonical = load_json(evaluation_attempt / "canonicalization.json")
    score_dir = evaluation_attempt / "score"
    canonical_maf = evaluation_attempt / "prediction.canonical.maf"
    summary_tsv = score_dir / "summary.tsv"
    if not canonical_maf.is_file() or canonical_maf.stat().st_size == 0:
        raise SummaryError(f"{tool} evaluation lacks a non-empty canonical MAF")
    if not summary_tsv.is_file() or summary_tsv.stat().st_size == 0:
        raise SummaryError(f"{tool} evaluation lacks a non-empty summary TSV")
    if not (run_attempt / "timing/RUN_COMPLETE").is_file():
        raise SummaryError(f"selected {tool} run lacks its timing completion marker")
    if (
        metrics.get("schema") != METRICS_SCHEMA
        or metrics.get("status") != "success"
        or int(metrics.get("exit_code", 1)) != 0
    ):
        raise SummaryError(f"selected {tool} run was not successful")
    timeout = metrics.get("timeout")
    if (
        not isinstance(timeout, dict)
        or not isinstance(timeout.get("enabled"), bool)
        or timeout.get("clock") != "CLOCK_MONOTONIC"
        or timeout.get("exceeded") is not False
        or timeout.get("elapsed_at_signal_seconds") is not None
        or timeout.get("sigterm_sent") is not False
        or timeout.get("sigkill_sent") is not False
        or timeout.get("deadline_overshoot_seconds") is not None
        or not isinstance(timeout.get("limit_seconds"), (int, float))
        or not math.isfinite(float(timeout.get("limit_seconds", float("nan"))))
        or (
            timeout.get("enabled") is True
            and float(timeout["limit_seconds"]) <= 0.0
        )
        or (
            timeout.get("enabled") is False
            and float(timeout["limit_seconds"]) != 0.0
        )
    ):
        raise SummaryError(f"selected {tool} run has contradictory timeout evidence")
    gnu_time = metrics.get("gnu_time")
    if metrics.get("gnu_time_parse_error") is not None or not isinstance(gnu_time, dict):
        raise SummaryError(f"selected {tool} run lacks valid GNU time metrics")
    if gnu_time.get("exit_status") != 0:
        raise SummaryError(f"selected {tool} GNU time exit status is not zero")
    for field in ("elapsed_seconds", "user_seconds", "system_seconds"):
        value = gnu_time.get(field)
        if not isinstance(value, (int, float)) or value < 0:
            raise SummaryError(f"selected {tool} GNU time field is invalid: {field}")
    if not isinstance(gnu_time.get("maximum_resident_set_kbytes"), int) or int(
        gnu_time["maximum_resident_set_kbytes"]
    ) <= 0:
        raise SummaryError(f"selected {tool} GNU time peak RSS is invalid")
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
            raise SummaryError(f"selected {tool} process-tree metric is invalid: {field}")
    total_threads = int(metrics["observed_max_total_threads"])
    total_single = int(metrics["observed_max_single_process_threads"])
    aligner_tree = int(metrics["observed_max_aligner_tree_threads"])
    aligner_single = int(metrics["observed_max_aligner_single_process_threads"])
    if not (
        aligner_single <= aligner_tree <= total_threads
        and aligner_single <= total_single <= total_threads
    ):
        raise SummaryError(
            f"selected {tool} thread metrics violate aligner/single-process/total hierarchy"
        )
    accepted_environment = config.get("openmp_environment")
    if not isinstance(accepted_environment, dict) or metrics.get(
        "environment_overrides"
    ) != accepted_environment:
        raise SummaryError(f"selected {tool} run lacks the accepted OpenMP environment")
    warm_cache, warm_launch, clock_consistency = validate_full_run_fairness_evidence(
        run_attempt, config, metrics, tool
    )
    if canonical.get("status") != "success":
        raise SummaryError(f"selected {tool} canonicalization was not successful")
    if not canonical.get("non_s_lines_unchanged") or not canonical.get(
        "s_line_payload_unchanged"
    ):
        raise SummaryError(f"{tool} MAF canonicalization altered non-name content")
    seen = set(canonical.get("canonical_sources_seen", []))
    whitelist = set(canonical.get("canonical_whitelist", []))
    if not seen <= whitelist:
        raise SummaryError(f"{tool} canonical MAF contains forbidden sources")
    if not (evaluation_attempt / "score/EVALUATION_COMPLETE").is_file():
        raise SummaryError(f"{tool} evaluation lacks EVALUATION_COMPLETE")

    expected_seeds = list(config["seeds"])
    expected_samples = int(config["samples"])
    expected_near = int(config["near"])
    runs = summary.get("runs")
    if (
        summary.get("truth_profile") != "all-homology"
        or summary.get("seeds") != expected_seeds
        or int(summary.get("samples", -1)) != expected_samples
        or int(summary.get("near", -1)) != expected_near
        or int(summary.get("run_count", -1)) != len(expected_seeds)
        or not isinstance(runs, list)
        or len(runs) != len(expected_seeds)
    ):
        raise SummaryError(f"{tool} evaluation summary violates the fixed scoring contract")
    expected_xml_paths = [
        score_dir / "raw" / f"mafComparator.seed-{seed}.xml" for seed in expected_seeds
    ]
    observed_xml_paths = sorted((score_dir / "raw").glob("mafComparator.seed-*.xml"))
    if observed_xml_paths != expected_xml_paths or any(
        path.stat().st_size == 0 for path in observed_xml_paths
    ):
        raise SummaryError(f"{tool} evaluation lacks exactly three fixed-seed XML files")
    per_run_metrics: dict[str, list[float]] = {
        "precision": [],
        "recall": [],
        "f1": [],
    }
    quality_runs: list[dict[str, object]] = []
    for expected_seed, run, xml_path in zip(
        expected_seeds, runs, expected_xml_paths, strict=True
    ):
        if not isinstance(run, dict):
            raise SummaryError(f"{tool} evaluation contains a non-object seed run")
        if (
            int(run.get("seed", -1)) != expected_seed
            or int(run.get("requested_samples", -1)) != expected_samples
            or int(run.get("near", -1)) != expected_near
            or Path(str(run.get("xml", ""))).resolve() != xml_path.resolve()
        ):
            raise SummaryError(f"{tool} evaluation seed run violates the fixed contract")
        for count_name in (
            "tp_truth_to_prediction",
            "false_negative",
            "tp_prediction_to_truth",
            "false_positive",
        ):
            value = run.get(count_name)
            if not isinstance(value, int) or value < 0:
                raise SummaryError(f"{tool} evaluation has invalid {count_name}")
        try:
            root_element = ET.parse(xml_path).getroot()
        except (OSError, ET.ParseError) as exc:
            raise SummaryError(f"{tool} evaluation XML is not parseable: {xml_path}") from exc
        if (
            int(root_element.attrib.get("seed", -1)) != expected_seed
            or int(root_element.attrib.get("numberOfSamples", -1)) != expected_samples
            or int(root_element.attrib.get("near", -1)) != expected_near
        ):
            raise SummaryError(f"{tool} evaluation XML violates the fixed scoring contract")
        for metric_name in per_run_metrics:
            value = run.get(metric_name)
            if not isinstance(value, (int, float)) or not math.isfinite(float(value)):
                raise SummaryError(f"{tool} evaluation has non-finite {metric_name}")
            numeric = float(value)
            if not 0.0 <= numeric <= 1.0:
                raise SummaryError(f"{tool} evaluation has out-of-range {metric_name}")
            per_run_metrics[metric_name].append(numeric)
        quality_runs.append(
            {
                "seed": expected_seed,
                "precision": float(run["precision"]),
                "recall": float(run["recall"]),
                "f1": float(run["f1"]),
                "requested_samples": int(run["requested_samples"]),
                "near": int(run["near"]),
                "xml": str(xml_path.relative_to(root)),
                "xml_sha256": sha256_file(xml_path),
            }
        )

    final_manifest_gate: dict[str, object] | None = None
    if tool == "ramag":
        delta = run_attempt / "ramag.mumreference.raw.delta"
        complete = run_attempt / "ramag.mumreference.raw.complete"
        if not complete.is_file():
            raise SummaryError("RaMA-G run lacks its atomic complete marker")
        manifest = load_json(run_attempt / "ramag.mumreference.raw.manifest.json")
        if manifest.get("status") != "success" or int(manifest.get("exit_code", 1)) != 0:
            raise SummaryError("RaMA-G manifest is not successful")
        dependencies = manifest.get("dependencies", {})
        if (
            not isinstance(dependencies, dict)
            or dependencies.get("sufkit_commit")
            != str(config["expected_sufkit_commit"])
            or dependencies.get("seqpro_commit")
            != str(config["expected_seqpro_commit"])
        ):
            raise SummaryError("RaMA-G manifest does not record accepted dependency commits")
        adapter_provenance = manifest.get("adapter_provenance", {})
        routes = manifest.get("actual_routes", {})
        backend = (
            adapter_provenance.get("sufkit.backend")
            if isinstance(adapter_provenance, dict)
            else None
        ) or (routes.get("index") if isinstance(routes, dict) else None)
        openmp = manifest.get("openmp")
        build = manifest.get("build")
        effective = manifest.get("effective_config")
        chaining = manifest.get("chaining")
        threading = manifest.get("threading")
        mam_tiling = manifest.get("mam_tiling")
        counts = manifest.get("counts")
        if (
            not isinstance(effective, dict)
            or effective.get("reference") != str(Path(str(config["reference"])).resolve())
            or effective.get("query") != str(Path(str(config["query"])).resolve())
            or effective.get("formats") != "delta"
            or effective.get("seed_mode") != "mumreference"
            or effective.get("selection_mode")
            != config.get("ramag_selection_mode")
            or int(effective.get("threads", 0)) != int(config["threads"])
            or int(effective.get("min_match", -1)) != int(config["min_match"])
            or int(effective.get("max_gap", -1)) != int(config["max_gap"])
            or int(effective.get("diag_diff", -1)) != int(config["diag_diff"])
            or float(effective.get("diag_factor", -1.0))
            != float(config["diag_factor"])
            or int(effective.get("min_cluster", -1)) != int(config["min_cluster"])
            or int(effective.get("break_length", -1))
            != int(config["break_length"])
            or int(effective.get("max_dp_cells", -1))
            != int(config["max_dp_cells"])
            or not isinstance(routes, dict)
            or routes.get("seed") != EXPECTED_SEED_ROUTE
            or routes.get("chain") != "sparse-exact-edge-components-v1"
            or not isinstance(chaining, dict)
            or chaining.get("route") != "sparse-exact-edge-components-v1"
            or not isinstance(adapter_provenance, dict)
            or adapter_provenance.get("sufkit.commit")
            != str(config["expected_sufkit_commit"])
            or adapter_provenance.get("sufkit.source_state")
            != "exact-commit-clean-at-configure"
            or not isinstance(openmp, dict)
            or openmp.get("enabled") is not True
            or openmp.get("runtime") in {None, "", "unavailable"}
            or int(openmp.get("runtime_max_threads", 0)) != int(config["threads"])
            or int(openmp.get("configured_requested_threads", 0))
            != int(config["threads"])
            or openmp.get("dynamic") is not False
            or openmp.get("sufkit_divsufsort_openmp") is not True
            or not isinstance(build, dict)
            or build.get("openmp_enabled") is not True
            or build.get("sufkit_divsufsort_openmp") is not True
            or not isinstance(threading, dict)
            or int(threading.get("seed_requested_threads", 0))
            != int(config["threads"])
            or int(threading.get("seed_scheduled_threads", 0))
            != int(config["threads"])
            or int(threading.get("seed_worker_threads", 0))
            != int(config["threads"])
            or int(threading.get("chaining_requested_threads", 0))
            != int(config["threads"])
            or int(threading.get("seed_task_count", 0))
            != int(threading.get("seed_tasks_completed", -1))
            or threading.get("seed_parallel_route")
            != "openmp-dynamic-stable-tile-boundary-tasks"
            or not isinstance(mam_tiling, dict)
            or int(mam_tiling.get("worker_cap", 0)) != int(config["threads"])
            or int(mam_tiling.get("tile_bases", 0)) != MAM_TILE_BASES
            or int(mam_tiling.get("oriented_queries", 0))
            != FULL_MAM_TASK_SHAPE["oriented_queries"]
            or int(mam_tiling.get("tile_tasks", 0))
            != FULL_MAM_TASK_SHAPE["tile_tasks"]
            or int(mam_tiling.get("tile_tasks_completed", -1))
            != FULL_MAM_TASK_SHAPE["tile_tasks"]
            or int(mam_tiling.get("boundary_tasks", 0))
            != FULL_MAM_TASK_SHAPE["boundary_tasks"]
            or int(mam_tiling.get("boundary_tasks_completed", -1))
            != FULL_MAM_TASK_SHAPE["boundary_tasks"]
            or int(threading.get("seed_task_count", 0))
            != FULL_MAM_TASK_SHAPE["total_tasks"]
            or int(threading.get("seed_tasks_completed", -1))
            != FULL_MAM_TASK_SHAPE["total_tasks"]
            or int(threading.get("chaining_worker_threads", 0)) <= 1
            or int(threading.get("extension_worker_threads", 0)) <= 1
            or threading.get("seed_chain_extension_parallel") is not True
            or int(metrics.get("observed_max_aligner_single_process_threads", 0)) <= 1
        ):
            raise SummaryError("RaMA-G manifest lacks the accepted OpenMP contract")
        mam_limits = mam_tiling.get("resource_limits")
        if not isinstance(mam_limits, dict) or not isinstance(counts, dict):
            raise SummaryError("RaMA-G manifest lacks tiled-MAM count/resource evidence")
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
        boundary_limit = int(mam_limits.get("boundary_mem_occurrences", -1))
        workspace_limit = int(mam_limits.get("workspace_bytes", -1))
        if (
            any(value < 0 for value in statistics.values())
            or boundary_limit <= 0
            or workspace_limit <= 0
            or statistics["short_query_tasks"]
            > FULL_MAM_TASK_SHAPE["oriented_queries"]
            or statistics["short_query_tasks"]
            > FULL_MAM_TASK_SHAPE["tile_tasks"]
            or statistics["tile_globally_maximal_mams"]
            > statistics["tile_raw_mams"]
            or statistics["boundary_patterns"] > statistics["boundary_raw_mems"]
            or statistics["boundary_reference_unique_patterns"]
            > statistics["boundary_patterns"]
            or statistics["boundary_recovered_mams"]
            > statistics["boundary_reference_unique_patterns"]
            or statistics["boundary_raw_mems"] > boundary_limit
            or statistics["workspace_baseline_bytes"]
            > statistics["workspace_peak_bytes"]
            or statistics["workspace_peak_bytes"] > workspace_limit
            or int(counts.get("mem_seeds", -1))
            != statistics["boundary_raw_mems"]
            or int(counts.get("mam_seeds", -1))
            != statistics["tile_globally_maximal_mams"]
            + statistics["boundary_recovered_mams"]
            or not 0
            <= int(counts.get("selected_seeds", -1))
            <= int(counts.get("mam_seeds", -1))
        ):
            raise SummaryError("RaMA-G manifest violates tiled-MAM statistical invariants")
        # Re-read the selected full manifest and the selected smoke manifest.
        # This is deliberately separate from the alignment-stage acceptance
        # check and runs before any comparison artifact or completion marker is
        # published.
        final_manifest_gate = validate_final_ramag_manifest_gate(
            root, config, run_attempt
        )
        full_backend = final_manifest_gate["full"]
        if not isinstance(full_backend, dict):
            raise SummaryError("final RaMA-G manifest gate produced invalid evidence")
        backend = full_backend.get("backend")
        snapshot = config.get("binary_snapshot")
        if not isinstance(snapshot, dict):
            raise SummaryError("accepted config lacks the RaMA-G binary snapshot")
        snapshot_path = Path(str(snapshot.get("snapshot", "")))
        expected_binary_sha = snapshot.get("snapshot_sha256")
        if (
            not snapshot_path.is_file()
            or sha256_file(snapshot_path) != expected_binary_sha
            or metrics.get("binary_sha256") != expected_binary_sha
            or Path(str(metrics.get("binary_path", ""))).resolve()
            != snapshot_path.resolve()
        ):
            raise SummaryError("RaMA-G runtime binary does not match the accepted snapshot")
        binary_sha = expected_binary_sha
    else:
        delta = run_attempt / "mummer4.mumreference.raw.delta"
        manifest = None
        backend = "local-mummer4-current-build"
        mummer_binaries = provenance.get("mummer4")
        wrapper = (
            mummer_binaries.get("wrapper") if isinstance(mummer_binaries, dict) else None
        )
        if (
            not isinstance(wrapper, dict)
            or metrics.get("binary_sha256") != wrapper.get("sha256")
            or Path(str(metrics.get("binary_path", ""))).resolve()
            != Path(str(config["mummer_nucmer"])).resolve()
        ):
            raise SummaryError("MUMmer4 runtime wrapper does not match accepted provenance")
        binary_sha = metrics.get("binary_sha256")
    if not delta.is_file() or delta.stat().st_size == 0:
        raise SummaryError(f"missing non-empty raw delta: {delta}")
    format_validation = load_json(run_attempt / "format-validation.json")
    if (
        format_validation.get("status") != "success"
        or int(format_validation.get("bytes", -1)) != delta.stat().st_size
        or format_validation.get("sha256") != sha256_file(delta)
        or int(format_validation.get("show_coords_exit_code", -1)) != 0
        or int(format_validation.get("delta_filter_exit_code", -1)) != 0
    ):
        raise SummaryError(f"{tool} raw delta no longer matches parser validation evidence")

    precision = summary.get("precision")
    recall = summary.get("recall")
    f1 = summary.get("f1")
    if not all(isinstance(value, dict) for value in (precision, recall, f1)):
        raise SummaryError(f"{tool} evaluation summary lacks aggregate scores")
    for metric_name, aggregate_values in (
        ("precision", precision),
        ("recall", recall),
        ("f1", f1),
    ):
        observed = [float(value) for value in per_run_metrics[metric_name]]
        expected_min = min(observed)
        expected_max = max(observed)
        ordered = sorted(observed)
        expected_median = ordered[len(ordered) // 2]
        aggregate_min = aggregate_values.get("min")
        aggregate_median = aggregate_values.get("median")
        aggregate_max = aggregate_values.get("max")
        if (
            not all(
                isinstance(value, (int, float)) and math.isfinite(float(value))
                for value in (aggregate_min, aggregate_median, aggregate_max)
            )
            or float(aggregate_min) != expected_min
            or float(aggregate_median) != expected_median
            or float(aggregate_max) != expected_max
        ):
            raise SummaryError(f"{tool} {metric_name} aggregate disagrees with seed runs")
    strict_validation = load_json(evaluation_attempt / "prediction.validation.json")
    parser_pair_count = validation.get("total_pair_count")
    strict_pair_count = strict_validation.get("total_pair_count")
    prediction_pair_count = summary.get("prediction_pair_count")
    if (
        not isinstance(parser_pair_count, int)
        or parser_pair_count < 0
        or strict_pair_count != parser_pair_count
        or prediction_pair_count != parser_pair_count
    ):
        raise SummaryError(
            f"{tool} parser and mafPairCounter prediction pair counts disagree"
        )
    postprocessing = {
        "delta_to_maf": load_successful_timing(
            evaluation_attempt / "normalize-prediction.timing.json",
            "normalize-prediction",
        ),
        "source_canonicalization": load_successful_timing(
            evaluation_attempt / "canonicalize-sources.timing.json",
            "canonicalize-sources",
        ),
        "strict_maf_validation": load_successful_timing(
            evaluation_attempt / "validate-canonical-maf.timing.json",
            "validate-canonical-maf",
        ),
        "f1_evaluation": load_successful_timing(
            evaluation_attempt / "evaluate.timing.json",
            "evaluate",
        ),
    }
    row: dict[str, object] = {
        "tool": "RaMA-G" if tool == "ramag" else "local MUMmer4 4.0.1 current build",
        "tool_id": tool,
        "seed_semantics": "reference-unique MAM",
        "requested_threads": int(config["threads"]),
        "run_environment_overrides": metrics.get("environment_overrides"),
        "observed_max_total_threads": metrics.get("observed_max_total_threads"),
        "observed_max_single_process_threads": metrics.get(
            "observed_max_single_process_threads"
        ),
        "observed_max_aligner_tree_threads": metrics.get(
            "observed_max_aligner_tree_threads"
        ),
        "observed_max_aligner_single_process_threads": metrics.get(
            "observed_max_aligner_single_process_threads"
        ),
        "actual_backend": backend,
        "wall_seconds": metric_value(metrics, "elapsed_seconds"),
        "monotonic_wall_seconds": clock_consistency["clock_values_seconds"].get(
            "monotonic_elapsed_seconds"
        ),
        "realtime_utc_wall_seconds": clock_consistency[
            "clock_values_seconds"
        ].get("realtime_utc_elapsed_seconds"),
        "clock_consistency_status": clock_consistency.get("status"),
        "resource_realtime_jump_count": clock_consistency[
            "resource_samples"
        ].get("realtime_jump_count"),
        "warm_cache_seconds": warm_cache.get("monotonic_elapsed_seconds"),
        "warm_cache_to_launch_seconds": warm_launch.get(
            "warm_cache_to_metrics_runner_seconds"
        ),
        "user_seconds": metric_value(metrics, "user_seconds"),
        "system_seconds": metric_value(metrics, "system_seconds"),
        "gnu_time_cpu_percent": metric_value(metrics, "cpu_percent"),
        "gnu_time_peak_rss_kbytes": metric_value(
            metrics, "maximum_resident_set_kbytes"
        ),
        "process_tree_peak_rss_bytes": metrics.get(
            "observed_process_tree_peak_rss_bytes"
        ),
        "delta_bytes": delta.stat().st_size,
        "prediction_pair_count": prediction_pair_count,
        "alignment_blocks": validation.get("blocks"),
        "forward_rows": validation.get("forward_rows"),
        "reverse_rows": validation.get("reverse_rows"),
        "human_coverage": total_coverage(validation, "simHuman"),
        "chimp_coverage": total_coverage(validation, "simChimp"),
        "precision_median": precision.get("median"),
        "precision_min": precision.get("min"),
        "precision_max": precision.get("max"),
        "recall_median": recall.get("median"),
        "recall_min": recall.get("min"),
        "recall_max": recall.get("max"),
        "f1_median": f1.get("median"),
        "f1_min": f1.get("min"),
        "f1_max": f1.get("max"),
        "quality_runs": quality_runs,
        "score_seeds": summary.get("seeds"),
        "samples": summary.get("samples"),
        "near": summary.get("near"),
        "binary_sha256": binary_sha,
        "raw_delta_sha256": sha256_file(delta),
        "canonical_maf_bytes": canonical_maf.stat().st_size,
        "canonical_maf_sha256": sha256_file(canonical_maf),
        "evaluation_summary_sha256": sha256_file(score_dir / "summary.json"),
        "postprocessing": postprocessing,
        "delta_to_maf_seconds": postprocessing["delta_to_maf"][
            "monotonic_elapsed_seconds"
        ],
        "source_canonicalization_seconds": postprocessing[
            "source_canonicalization"
        ]["monotonic_elapsed_seconds"],
        "strict_maf_validation_seconds": postprocessing["strict_maf_validation"][
            "monotonic_elapsed_seconds"
        ],
        "f1_evaluation_seconds": postprocessing["f1_evaluation"][
            "monotonic_elapsed_seconds"
        ],
        "provenance_capture": str(provenance_attempt.relative_to(root)),
        "run_attempt": str(run_attempt.relative_to(root)),
        "evaluation_attempt": str(evaluation_attempt.relative_to(root)),
        "completion_state": "complete",
    }
    if manifest is not None:
        row["ramag_manifest"] = str(
            (run_attempt / "ramag.mumreference.raw.manifest.json").relative_to(root)
        )
        snapshot = config.get("binary_snapshot")
        source_tree = provenance.get("source")
        row["binary_provenance"] = snapshot
        row["ramag_elf_sha256"] = (
            snapshot.get("snapshot_sha256") if isinstance(snapshot, dict) else None
        )
        row["mummer_wrapper_sha256"] = None
        row["mummer_elf_sha256"] = None
        row["mummer_library_sha256"] = None
        row["source_tree_sha256"] = (
            source_tree.get("tree_sha256") if isinstance(source_tree, dict) else None
        )
        openmp = manifest["openmp"]
        row["openmp_enabled"] = openmp.get("enabled")
        row["openmp_runtime"] = openmp.get("runtime")
        row["openmp_runtime_max_threads"] = openmp.get("runtime_max_threads")
        row["openmp_configured_requested_threads"] = openmp.get(
            "configured_requested_threads"
        )
        row["sufkit_divsufsort_openmp"] = openmp.get(
            "sufkit_divsufsort_openmp"
        )
        if final_manifest_gate is None:
            raise SummaryError("final RaMA-G manifest evidence was not retained")
        row["ramag_final_manifest_gate"] = final_manifest_gate
        full_backend = final_manifest_gate["full"]
        smoke_backend = final_manifest_gate["smoke"]
        if not isinstance(full_backend, dict) or not isinstance(smoke_backend, dict):
            raise SummaryError("final RaMA-G backend evidence is invalid")
        row["full_backend"] = full_backend.get("backend")
        row["smoke_backend"] = smoke_backend.get("backend")
        row["caps_min_reference_bases"] = full_backend.get(
            "caps_min_reference_bases"
        )
        row["input_route"] = full_backend.get("input_route")
        row["input_parallel_route"] = full_backend.get("input_parallel_route")
    else:
        binaries = provenance.get("mummer4")
        current_tree = provenance.get("mummer4_current_tree")
        row["binary_provenance"] = binaries
        row["ramag_elf_sha256"] = None
        row["mummer_wrapper_sha256"] = (
            binaries.get("wrapper", {}).get("sha256")
            if isinstance(binaries, dict) and isinstance(binaries.get("wrapper"), dict)
            else None
        )
        row["mummer_elf_sha256"] = (
            binaries.get("elf", {}).get("sha256")
            if isinstance(binaries, dict) and isinstance(binaries.get("elf"), dict)
            else None
        )
        row["mummer_library_sha256"] = (
            binaries.get("library", {}).get("sha256")
            if isinstance(binaries, dict) and isinstance(binaries.get("library"), dict)
            else None
        )
        row["source_tree_sha256"] = (
            current_tree.get("tree_sha256") if isinstance(current_tree, dict) else None
        )
        row["openmp_enabled"] = None
        row["openmp_runtime"] = None
        row["openmp_runtime_max_threads"] = None
        row["openmp_configured_requested_threads"] = None
        row["sufkit_divsufsort_openmp"] = None
        row["ramag_final_manifest_gate"] = None
        row["full_backend"] = row["actual_backend"]
        row["smoke_backend"] = None
        row["caps_min_reference_bases"] = None
        row["input_route"] = None
        row["input_parallel_route"] = None
    return row


def format_cell(value: object) -> str:
    if value is None:
        return "NA"
    if isinstance(value, float):
        return "NA" if math.isnan(value) else f"{value:.9f}"
    return str(value)


def write_exclusive(path: Path, text: str) -> None:
    with path.open("x", encoding="utf-8", newline="\n") as handle:
        handle.write(text)


def relative_difference(left: float, right: float) -> float:
    denominator = max(abs(left), abs(right))
    return 0.0 if denominator == 0.0 else abs(left - right) / denominator


def build_clock_robustness_gate(
    rows_by_id: dict[str, dict[str, object]],
) -> dict[str, object]:
    mummer = rows_by_id["mummer4"]
    ramag = rows_by_id["ramag"]
    values = {
        "mummer4_gnu_seconds": float(mummer["wall_seconds"]),
        "ramag_gnu_seconds": float(ramag["wall_seconds"]),
        "mummer4_monotonic_seconds": float(mummer["monotonic_wall_seconds"]),
        "ramag_monotonic_seconds": float(ramag["monotonic_wall_seconds"]),
        "mummer4_realtime_seconds": float(mummer["realtime_utc_wall_seconds"]),
        "ramag_realtime_seconds": float(ramag["realtime_utc_wall_seconds"]),
    }
    if any(not math.isfinite(value) or value <= 0.0 for value in values.values()):
        raise SummaryError("clock robustness inputs must be finite and positive")
    distortions = {
        "mummer4_gnu_over_monotonic": (
            values["mummer4_gnu_seconds"]
            / values["mummer4_monotonic_seconds"]
        ),
        "ramag_gnu_over_monotonic": (
            values["ramag_gnu_seconds"] / values["ramag_monotonic_seconds"]
        ),
        "mummer4_realtime_over_monotonic": (
            values["mummer4_realtime_seconds"]
            / values["mummer4_monotonic_seconds"]
        ),
        "ramag_realtime_over_monotonic": (
            values["ramag_realtime_seconds"]
            / values["ramag_monotonic_seconds"]
        ),
    }
    gnu_distortion_difference = relative_difference(
        distortions["ramag_gnu_over_monotonic"],
        distortions["mummer4_gnu_over_monotonic"],
    )
    realtime_distortion_difference = relative_difference(
        distortions["ramag_realtime_over_monotonic"],
        distortions["mummer4_realtime_over_monotonic"],
    )
    hard_conditions = {
        "gnu_headline_ramag_not_slower": (
            values["ramag_gnu_seconds"] <= values["mummer4_gnu_seconds"]
        ),
        "monotonic_ramag_not_slower": (
            values["ramag_monotonic_seconds"]
            <= values["mummer4_monotonic_seconds"]
        ),
    }
    diagnostic_conditions = {
        "gnu_over_monotonic_distortion_relative_difference_at_most_1pct": (
            gnu_distortion_difference <= CLOCK_CONSISTENCY_TOLERANCE_FRACTION
        ),
        "realtime_over_monotonic_distortion_relative_difference_at_most_1pct": (
            realtime_distortion_difference
            <= CLOCK_CONSISTENCY_TOLERANCE_FRACTION
        ),
    }
    warnings: list[str] = []
    if not diagnostic_conditions[
        "gnu_over_monotonic_distortion_relative_difference_at_most_1pct"
    ]:
        warnings.append(
            "cross-tool GNU/monotonic distortion differs by more than 1%"
        )
    if not diagnostic_conditions[
        "realtime_over_monotonic_distortion_relative_difference_at_most_1pct"
    ]:
        warnings.append(
            "cross-tool realtime/monotonic distortion differs by more than 1%"
        )
    for tool in TOOLS:
        jump_count = int(rows_by_id[tool]["resource_realtime_jump_count"])
        if jump_count > 0:
            warnings.append(
                f"{tool} resource samples contain {jump_count} positive realtime jump(s)"
            )
        if rows_by_id[tool]["clock_consistency_status"] == "warning":
            warnings.append(f"{tool} selected attempt has warning-quality clock evidence")
    passed = all(hard_conditions.values())
    return {
        "schema": "ramag.human-chimp-clock-robustness-gate.v2",
        "status": "passed" if passed else "failed",
        "diagnostic_status": "warning" if warnings else "clear",
        "tolerance_fraction": CLOCK_CONSISTENCY_TOLERANCE_FRACTION,
        "formulas": {
            "distortion": "tool_clock_elapsed_seconds / tool_monotonic_elapsed_seconds",
            "relative_difference": "abs(left-right) / max(abs(left),abs(right))",
            "hard_runtime_requirements": (
                "ramag_gnu<=mummer4_gnu and "
                "ramag_monotonic<=mummer4_monotonic"
            ),
            "diagnostic_policy": (
                "distortion parity and positive realtime jumps are warnings only"
            ),
        },
        "elapsed_values_seconds": values,
        "distortion_values": distortions,
        "gnu_distortion_relative_difference": gnu_distortion_difference,
        "realtime_distortion_relative_difference": realtime_distortion_difference,
        "conditions": hard_conditions,
        "hard_conditions": hard_conditions,
        "diagnostic_conditions": diagnostic_conditions,
        "warnings": warnings,
        "tool_attempt_clock_quality": {
            tool: {
                "status": rows_by_id[tool]["clock_consistency_status"],
                "resource_realtime_jump_count": rows_by_id[tool][
                    "resource_realtime_jump_count"
                ],
            }
            for tool in TOOLS
        },
    }


def validate_early_alignment_speed_gate(
    root: Path,
    config: dict[str, object],
    rows_by_id: dict[str, dict[str, object]],
) -> dict[str, object]:
    required_markers = (
        "MUMMER_BASELINE_ACCEPTED",
        "ALIGNMENT_SPEED_ACCEPTED",
        "FULL_RUNS_COMPLETE",
        "EVALUATIONS_COMPLETE",
    )
    missing = [name for name in required_markers if not (root / name).is_file()]
    if missing:
        raise SummaryError(f"missing root completion prerequisite(s): {missing}")
    if (root / "ALIGNMENT_SPEED_REJECTED").exists():
        raise SummaryError("alignment speed gate has contradictory accepted/rejected markers")
    if config.get("schema") != CONFIG_SCHEMA or config.get("metrics_schema") != METRICS_SCHEMA:
        raise SummaryError("accepted config predates the early alignment speed gate")
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
        raise SummaryError("accepted speed/F1 publication policy is invalid")

    budget_path = root / "speed-budget.json"
    gate_path = root / "alignment-speed-gate.json"
    budget = load_json(budget_path)
    gate = load_json(gate_path)
    mummer = rows_by_id["mummer4"]
    ramag = rows_by_id["ramag"]
    mummer_attempt = root / str(mummer["run_attempt"])
    ramag_attempt = root / str(ramag["run_attempt"])
    mummer_metrics = mummer_attempt / "timing/metrics.json"
    ramag_metrics = ramag_attempt / "timing/metrics.json"
    mummer_monotonic = float(mummer["monotonic_wall_seconds"])
    mummer_gnu = float(mummer["wall_seconds"])
    ramag_monotonic = float(ramag["monotonic_wall_seconds"])
    ramag_gnu = float(ramag["wall_seconds"])
    budget_cutoff = budget.get("cutoff_seconds")
    gate_mummer = gate.get("mummer4")
    gate_ramag = gate.get("ramag")
    conditions = gate.get("conditions")
    if (
        budget.get("schema") != "ramag.human-chimp-speed-budget.v1"
        or budget.get("status") != "accepted"
        or budget.get("result_root") != str(root.resolve())
        or budget.get("source_tool") != "mummer4"
        or budget.get("source_attempt") != str(mummer_attempt.relative_to(root))
        or budget.get("source_metrics") != str(mummer_metrics.relative_to(root))
        or budget.get("source_metrics_sha256") != sha256_file(mummer_metrics)
        or budget.get("source_metrics_schema") != METRICS_SCHEMA
        or not isinstance(budget_cutoff, (int, float))
        or float(budget_cutoff) != mummer_monotonic
        or budget.get("cutoff_clock") != "CLOCK_MONOTONIC"
        or budget.get("equal_is_accepted") is not True
    ):
        raise SummaryError("frozen MUMmer4 speed budget is invalid")
    if (
        gate.get("schema") != "ramag.human-chimp-alignment-speed-gate.v1"
        or gate.get("status") != "accepted"
        or gate.get("reason") != "dual_clock_passed"
        or gate.get("scope") != "preliminary single observation"
        or gate.get("endpoint") != config.get("runtime_endpoint")
        or gate.get("equal_is_accepted") is not True
        or gate.get("hard_deadline_clock") != "CLOCK_MONOTONIC"
        or gate.get("hard_deadline_seconds") != budget_cutoff
        or gate.get("budget_sha256") != sha256_file(budget_path)
        or not isinstance(gate_mummer, dict)
        or not isinstance(gate_ramag, dict)
        or not isinstance(conditions, dict)
        or conditions
        != {
            "ramag_completed_before_hard_deadline": True,
            "ramag_gnu_not_slower": True,
            "ramag_monotonic_not_slower": True,
        }
        or gate_mummer.get("attempt") != str(mummer_attempt.relative_to(root))
        or gate_mummer.get("metrics_sha256") != sha256_file(mummer_metrics)
        or gate_mummer.get("metrics_schema") != METRICS_SCHEMA
        or gate_mummer.get("gnu_elapsed_seconds") != mummer_gnu
        or gate_mummer.get("monotonic_elapsed_seconds") != mummer_monotonic
        or gate_ramag.get("attempt") != str(ramag_attempt.relative_to(root))
        or gate_ramag.get("metrics_sha256") != sha256_file(ramag_metrics)
        or gate_ramag.get("metrics_schema") != METRICS_SCHEMA
        or gate_ramag.get("status") != "success"
        or gate_ramag.get("gnu_elapsed_seconds") != ramag_gnu
        or gate_ramag.get("monotonic_elapsed_seconds") != ramag_monotonic
        or ramag_gnu > mummer_gnu
        or ramag_monotonic > mummer_monotonic
    ):
        raise SummaryError("published early alignment speed gate is invalid")
    return gate


def summarize(root: Path) -> dict[str, object]:
    if (root / "COMPARISON_COMPLETE").exists():
        raise SummaryError(f"comparison is already complete: {root}")
    for name in ("comparison.json", "comparison.tsv", "comparison.md"):
        if (root / name).exists():
            raise SummaryError(f"refusing to overwrite partial summary: {root / name}")
    config = load_json(root / "accepted-config.json")
    provenance, provenance_attempt = load_provenance(root)
    rows = [
        build_row(root, config, provenance, provenance_attempt, tool) for tool in TOOLS
    ]
    rows_by_id = {str(row["tool_id"]): row for row in rows}
    early_alignment_speed_gate = validate_early_alignment_speed_gate(
        root, config, rows_by_id
    )
    mummer_elapsed = float(rows_by_id["mummer4"]["wall_seconds"])
    ramag_elapsed = float(rows_by_id["ramag"]["wall_seconds"])
    performance_passed = ramag_elapsed <= mummer_elapsed
    performance_gate = {
        "status": "passed" if performance_passed else "failed",
        "metric": "gnu_time.elapsed_seconds",
        "requirement": "ramag_elapsed_seconds <= mummer4_elapsed_seconds",
        "ramag_elapsed_seconds": ramag_elapsed,
        "mummer4_elapsed_seconds": mummer_elapsed,
        "equal_is_accepted": True,
    }
    clock_robustness_gate = build_clock_robustness_gate(rows_by_id)
    clock_robustness_passed = clock_robustness_gate["status"] == "passed"
    completion_passed = (
        early_alignment_speed_gate["status"] == "accepted"
        and performance_passed
        and clock_robustness_passed
    )
    mummer_f1 = float(rows_by_id["mummer4"]["f1_median"])
    ramag_f1 = float(rows_by_id["ramag"]["f1_median"])
    quality_comparison = {
        "policy": "report-only; no minimum or comparative F1 completion threshold",
        "mummer4_f1_median": mummer_f1,
        "ramag_f1_median": ramag_f1,
        "ramag_minus_mummer4_f1_median": ramag_f1 - mummer_f1,
        "existing_ramag_diagnostic_f1": EXISTING_RAMAG_DIAGNOSTIC_F1,
        "ramag_minus_existing_diagnostic_f1": (
            ramag_f1 - EXISTING_RAMAG_DIAGNOSTIC_F1
        ),
    }
    for row in rows:
        row["performance_gate_status"] = performance_gate["status"]
        row["clock_robustness_gate_status"] = clock_robustness_gate["status"]
        row["clock_diagnostic_status"] = clock_robustness_gate[
            "diagnostic_status"
        ]
        row["f1_median_minus_mummer4"] = float(row["f1_median"]) - mummer_f1
        row["f1_median_minus_existing_ramag_diagnostic"] = (
            float(row["f1_median"]) - EXISTING_RAMAG_DIAGNOSTIC_F1
            if row["tool_id"] == "ramag"
            else None
        )
    result_status = (
        "complete"
        if completion_passed
        else "performance_gate_failed"
        if not performance_passed
        else "monotonic_performance_gate_failed"
    )
    ramag_manifest_gate = rows_by_id["ramag"].get("ramag_final_manifest_gate")
    if not isinstance(ramag_manifest_gate, dict) or ramag_manifest_gate.get(
        "status"
    ) != "passed":
        raise SummaryError("final RaMA-G manifest gate did not pass")
    result = {
        "schema": "ramag.human-chimp-preliminary-comparison.v4",
        "status": result_status,
        "scope": "preliminary single observation",
        "runtime_repetitions": int(config["repetitions"]),
        "runtime_aggregation": "none; one observation per tool",
        "quality_truth_profile": "all-homology",
        "quality_is_sampled": True,
        "quality_samples": int(config["samples"]),
        "quality_near": int(config["near"]),
        "quality_seeds": config["seeds"],
        "quality_policy": config["quality_policy"],
        "quality_comparison": quality_comparison,
        "openmp_environment": config["openmp_environment"],
        "attempt_clock_acceptance": {
            "status": "hard-gates-passed",
            "tolerance_fraction": CLOCK_CONSISTENCY_TOLERANCE_FRACTION,
            "required_clocks": [
                "gnu_time_elapsed_seconds",
                "monotonic_elapsed_seconds",
                "realtime_utc_elapsed_seconds",
            ],
            "positive_realtime_jumps": "warning; diagnostic only",
            "realtime_backward_or_nonpositive": "hard quarantine",
        },
        "clock_robustness_gate": clock_robustness_gate,
        "early_alignment_speed_gate": early_alignment_speed_gate,
        "performance_gate": performance_gate,
        "ramag_final_manifest_gate": ramag_manifest_gate,
        "warning": (
            "This single run cannot establish a stable speedup or that RaMA-G "
            "outperforms MUMmer4. Requested and observed parallelism differ."
        ),
        "provenance_capture": str(provenance_attempt.relative_to(root)),
        "tools": rows,
    }
    write_exclusive(
        root / "comparison.json",
        json.dumps(result, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
    )

    columns = [
        "tool",
        "performance_gate_status",
        "clock_robustness_gate_status",
        "clock_diagnostic_status",
        "seed_semantics",
        "requested_threads",
        "observed_max_aligner_tree_threads",
        "observed_max_aligner_single_process_threads",
        "observed_max_total_threads",
        "openmp_enabled",
        "openmp_runtime",
        "openmp_runtime_max_threads",
        "openmp_configured_requested_threads",
        "sufkit_divsufsort_openmp",
        "actual_backend",
        "full_backend",
        "smoke_backend",
        "caps_min_reference_bases",
        "input_route",
        "input_parallel_route",
        "wall_seconds",
        "monotonic_wall_seconds",
        "realtime_utc_wall_seconds",
        "clock_consistency_status",
        "resource_realtime_jump_count",
        "warm_cache_seconds",
        "warm_cache_to_launch_seconds",
        "user_seconds",
        "system_seconds",
        "gnu_time_peak_rss_kbytes",
        "process_tree_peak_rss_bytes",
        "delta_bytes",
        "prediction_pair_count",
        "raw_delta_sha256",
        "canonical_maf_bytes",
        "canonical_maf_sha256",
        "evaluation_summary_sha256",
        "delta_to_maf_seconds",
        "source_canonicalization_seconds",
        "strict_maf_validation_seconds",
        "f1_evaluation_seconds",
        "precision_median",
        "recall_median",
        "f1_median",
        "f1_min",
        "f1_max",
        "f1_median_minus_mummer4",
        "f1_median_minus_existing_ramag_diagnostic",
        "binary_sha256",
        "source_tree_sha256",
        "ramag_elf_sha256",
        "mummer_wrapper_sha256",
        "mummer_elf_sha256",
        "mummer_library_sha256",
        "completion_state",
    ]
    tsv = ["\t".join(columns)]
    for row in rows:
        tsv.append("\t".join(format_cell(row.get(column)) for column in columns))
    write_exclusive(root / "comparison.tsv", "\n".join(tsv) + "\n")

    ramag_row = rows_by_id["ramag"]
    hard_clock_conditions = clock_robustness_gate["hard_conditions"]
    if not isinstance(hard_clock_conditions, dict):
        raise SummaryError("clock gate lacks hard conditions")
    monotonic_passed = bool(hard_clock_conditions["monotonic_ramag_not_slower"])
    clock_warnings = clock_robustness_gate["warnings"]
    if not isinstance(clock_warnings, list):
        raise SummaryError("clock gate warnings are invalid")
    md = [
        "# RaMA-G 与 MUMmer4 Human–Chimp 初步对比",
        "",
        (
            "> 性能门禁：**通过**（RaMA-G GNU time elapsed <= MUMmer4）。"
            if performance_passed
            else "> 性能门禁：**失败**（RaMA-G GNU time elapsed > MUMmer4）；不发布 `COMPARISON_COMPLETE`。"
        ),
        "",
        f"Headline 仍比较同一字段 `gnu_time.elapsed_seconds`：RaMA-G={ramag_elapsed:.9f}s，"
        f"MUMmer4={mummer_elapsed:.9f}s；相等视为通过。",
        "",
        (
            "> monotonic 性能门禁：**通过**（RaMA-G monotonic elapsed <= MUMmer4）。"
            if monotonic_passed
            else "> monotonic 性能门禁：**失败**；不发布 `COMPARISON_COMPLETE`。"
        ),
        "发布硬门禁只比较两个同名时钟：RaMA-G 的 GNU elapsed 和 monotonic elapsed "
        "都不得慢于 MUMmer4。GNU/monotonic 与 realtime/monotonic 的跨工具膨胀差异，"
        "以及正向 realtime 跳变，仅作为 warning/诊断，不参与完成标记判定。"
        f"本次 GNU distortion 相对差={clock_robustness_gate['gnu_distortion_relative_difference']:.9%}，"
        f"realtime distortion 相对差={clock_robustness_gate['realtime_distortion_relative_difference']:.9%}。",
        (
            "时钟诊断：无 warning。"
            if not clock_warnings
            else "时钟诊断 warning：" + "；".join(str(item) for item in clock_warnings) + "。"
        ),
        "单次 attempt 的正向 realtime 跳变或三时钟超过 1% 仅标 warning；"
        "非有限/非正、端点矛盾或 realtime 倒退才硬隔离。"
        "因此绝对 GNU wall 可能受 WSL 校时膨胀影响，诊断仍随结果完整保留。",
        "",
        "> 这是同一主机、16-thread 请求预算下各一次的初步观测。它不能证明稳定加速，"
        "也不能据此宣称 RaMA-G 已超过 MUMmer4。",
        "",
        "| 工具 | GNU wall (s) | monotonic wall (s) | clock quality | GNU time peak RSS (KiB) | observed threads | Precision | Recall | F1 |",
        "|---|---:|---:|---|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        md.append(
            "| {tool} | {wall} | {monotonic} | {clock} | {rss} | {threads} | {precision} | {recall} | {f1} |".format(
                tool=row["tool"],
                wall=format_cell(row["wall_seconds"]),
                monotonic=format_cell(row["monotonic_wall_seconds"]),
                clock=format_cell(row["clock_consistency_status"]),
                rss=format_cell(row["gnu_time_peak_rss_kbytes"]),
                threads=format_cell(row["observed_max_aligner_tree_threads"]),
                precision=format_cell(row["precision_median"]),
                recall=format_cell(row["recall_median"]),
                f1=format_cell(row["f1_median"]),
            )
        )
    md.extend(
        [
            "",
            "质量分数使用 `all-homology` 真值、1,000,000 samples、`near=0` 和三个固定 seed；"
            "表中质量指标为三次抽样的中位数。F1 只报告，不设最低值或胜负门禁；"
            "运行时间没有取中位数。",
            f"RaMA-G median F1 - 本轮 MUMmer4 median F1 = "
            f"{quality_comparison['ramag_minus_mummer4_f1_median']:.10f}；"
            f"RaMA-G median F1 - 既有诊断值 {EXISTING_RAMAG_DIAGNOSTIC_F1:.10f} = "
            f"{quality_comparison['ramag_minus_existing_diagnostic_f1']:.10f}。",
            "",
            "## 固定 seed 质量分数",
            "",
            "| 工具 | seed | Precision | Recall | F1 | XML SHA-256 |",
            "|---|---:|---:|---:|---:|---|",
        ]
    )
    for row in rows:
        quality_runs = row["quality_runs"]
        if not isinstance(quality_runs, list):
            raise SummaryError("quality run details are invalid")
        for quality_run in quality_runs:
            if not isinstance(quality_run, dict):
                raise SummaryError("quality run detail is invalid")
            md.append(
                "| {tool} | {seed} | {precision} | {recall} | {f1} | `{xml_sha}` |".format(
                    tool=row["tool"],
                    seed=quality_run["seed"],
                    precision=format_cell(quality_run["precision"]),
                    recall=format_cell(quality_run["recall"]),
                    f1=format_cell(quality_run["f1"]),
                    xml_sha=quality_run["xml_sha256"],
                )
            )
    md.extend(
        [
            "",
            "## 派生产物与后处理时间",
            "",
            "| 工具 | raw delta SHA-256 | canonical MAF SHA-256 | evaluation summary SHA-256 | delta→MAF (s) | canonicalize (s) | strict validate (s) | F1 evaluation (s) |",
            "|---|---|---|---|---:|---:|---:|---:|",
        ]
    )
    for row in rows:
        md.append(
            "| {tool} | `{delta}` | `{maf}` | `{summary}` | {convert} | {canonicalize} | {validate} | {evaluate} |".format(
                tool=row["tool"],
                delta=row["raw_delta_sha256"],
                maf=row["canonical_maf_sha256"],
                summary=row["evaluation_summary_sha256"],
                convert=format_cell(row["delta_to_maf_seconds"]),
                canonicalize=format_cell(row["source_canonicalization_seconds"]),
                validate=format_cell(row["strict_maf_validation_seconds"]),
                evaluate=format_cell(row["f1_evaluation_seconds"]),
            )
        )
    md.extend(
        [
            "",
            "两方运行均记录 `OMP_NUM_THREADS=16`、`OMP_DYNAMIC=FALSE`、"
            "`OMP_PROC_BIND=SPREAD`、`OMP_PLACES=CORES`、`OMP_THREAD_LIMIT=16` 和 "
            "`OMP_MAX_ACTIVE_LEVELS=1`；未设置 `OMP_WAIT_POLICY`。"
            "MUMmer4 不使用 OpenMP，但仍在同一记录环境中运行。",
            "",
            "RaMA-G 后端文字直接来自最终重验的 manifest："
            f"full=`{ramag_row['full_backend']}`，smoke=`{ramag_row['smoke_backend']}`，"
            f"切换阈值={ramag_row['caps_min_reference_bases']} bases（64 MiB）。"
            f"full 输入路由=`{ramag_row['input_route']}`，并行路由="
            f"`{ramag_row['input_parallel_route']}`（2/2 workers）。",
            "RaMA-G 顺序复用一个 16-thread OpenMP 预算：full 索引构建使用上述实际后端，"
            "reference-MAM seed 以 4 MiB oriented-query tiles 加内部边界 MEM 恢复形成"
            "稳定任务并动态调度，随后 chaining groups 与 chain extension 最多使用 16 个 worker。"
            "各阶段的实际 worker 数由 manifest 和进程树观测共同记录；"
            "MUMmer4 的 query-level 有效并行度仍受 4 条 query contig 限制。",
        ]
    )
    write_exclusive(root / "comparison.md", "\n".join(md) + "\n")
    if completion_passed:
        # The marker is published last and with O_EXCL semantics.
        fd = os.open(
            root / "COMPARISON_COMPLETE", os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644
        )
        with os.fdopen(fd, "w", encoding="ascii") as handle:
            handle.write("success\n")
    return result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result-root", required=True, type=Path)
    args = parser.parse_args(argv)
    result = summarize(args.result_root.resolve())
    print(json.dumps(result, indent=2, sort_keys=True, ensure_ascii=False))
    return 0 if result["status"] == "complete" else 3


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except SummaryError as error:
        print(f"summarize_comparison: {error}", file=sys.stderr)
        raise SystemExit(2)

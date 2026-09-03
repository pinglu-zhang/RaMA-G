#!/usr/bin/env python3
"""Run one command under GNU time and sample its Linux process tree."""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import math
import os
import shlex
import signal
import subprocess
import sys
import time
from pathlib import Path


class MetricsError(RuntimeError):
    """A runner setup or result parsing error."""


CLOCK_CONSISTENCY_TOLERANCE_FRACTION = 0.01
TIMEOUT_TERMINATION_GRACE_SECONDS = 10.0


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def parse_elapsed(value: str) -> float:
    parts = value.strip().split(":")
    try:
        if len(parts) == 2:
            minutes, seconds = parts
            return float(minutes) * 60.0 + float(seconds)
        if len(parts) == 3:
            hours, minutes, seconds = parts
            return float(hours) * 3600.0 + float(minutes) * 60.0 + float(seconds)
    except ValueError as exc:
        raise MetricsError(f"invalid GNU time elapsed value: {value!r}") from exc
    raise MetricsError(f"invalid GNU time elapsed value: {value!r}")


def parse_gnu_time(path: Path) -> dict[str, object]:
    fields: dict[str, str] = {}
    for raw_line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw_line.lstrip()
        if ": " not in line:
            continue
        key, value = line.split(": ", 1)
        fields[key] = value.strip()

    def float_field(name: str) -> float | None:
        value = fields.get(name)
        return None if value is None else float(value)

    def int_field(name: str) -> int | None:
        value = fields.get(name)
        return None if value is None else int(value)

    elapsed_raw = fields.get("Elapsed (wall clock) time (h:mm:ss or m:ss)")
    cpu_raw = fields.get("Percent of CPU this job got")
    return {
        "raw_fields": fields,
        "user_seconds": float_field("User time (seconds)"),
        "system_seconds": float_field("System time (seconds)"),
        "elapsed_seconds": parse_elapsed(elapsed_raw) if elapsed_raw else None,
        "cpu_percent": float(cpu_raw.rstrip("%")) if cpu_raw else None,
        "maximum_resident_set_kbytes": int_field("Maximum resident set size (kbytes)"),
        "major_page_faults": int_field("Major (requiring I/O) page faults"),
        "minor_page_faults": int_field("Minor (reclaiming a frame) page faults"),
        "voluntary_context_switches": int_field("Voluntary context switches"),
        "involuntary_context_switches": int_field("Involuntary context switches"),
        "filesystem_inputs": int_field("File system inputs"),
        "filesystem_outputs": int_field("File system outputs"),
        "exit_status": int_field("Exit status"),
    }


def _proc_stat(pid: int) -> tuple[int, int] | None:
    try:
        text = Path(f"/proc/{pid}/stat").read_text(encoding="ascii")
    except (FileNotFoundError, PermissionError, ProcessLookupError):
        return None
    close = text.rfind(")")
    if close < 0:
        return None
    fields = text[close + 2 :].split()
    try:
        return int(fields[1]), int(fields[11]) + int(fields[12])
    except (IndexError, ValueError):
        return None


def _status(pid: int) -> dict[str, int | str] | None:
    values: dict[str, int | str] = {}
    try:
        lines = Path(f"/proc/{pid}/status").read_text(encoding="ascii").splitlines()
    except (FileNotFoundError, PermissionError, ProcessLookupError):
        return None
    for line in lines:
        if ":" not in line:
            continue
        key, value = line.split(":", 1)
        value = value.strip()
        if key == "Name":
            values["name"] = value
        elif key in {"VmRSS", "VmSize"}:
            try:
                values[key] = int(value.split()[0]) * 1024
            except (IndexError, ValueError):
                values[key] = 0
        elif key == "Threads":
            try:
                values["threads"] = int(value)
            except ValueError:
                values["threads"] = 0
    values.setdefault("name", "unknown")
    values.setdefault("VmRSS", 0)
    values.setdefault("VmSize", 0)
    values.setdefault("threads", 0)
    return values


def _descendants(root_pid: int) -> list[tuple[int, int, int]]:
    records: dict[int, tuple[int, int]] = {}
    try:
        proc_entries = list(Path("/proc").iterdir())
    except OSError:
        return []
    for entry in proc_entries:
        if not entry.name.isdigit():
            continue
        pid = int(entry.name)
        stat = _proc_stat(pid)
        if stat is not None:
            records[pid] = stat
    children: dict[int, list[int]] = {}
    for pid, (ppid, _ticks) in records.items():
        children.setdefault(ppid, []).append(pid)
    selected: list[int] = []
    pending = [root_pid]
    visited: set[int] = set()
    while pending:
        pid = pending.pop()
        if pid in visited:
            continue
        visited.add(pid)
        if pid in records:
            selected.append(pid)
        pending.extend(children.get(pid, ()))
    return [(pid, records[pid][0], records[pid][1]) for pid in sorted(selected)]


def _utc_from_realtime_ns(realtime_ns: int) -> str:
    seconds, nanoseconds = divmod(realtime_ns, 1_000_000_000)
    instant = dt.datetime.fromtimestamp(seconds, tz=dt.timezone.utc).replace(
        microsecond=nanoseconds // 1_000
    )
    return instant.isoformat(timespec="microseconds").replace("+00:00", "Z")


def _relative_difference(left: float, right: float) -> float:
    denominator = max(abs(left), abs(right))
    if denominator == 0.0:
        return 0.0
    return abs(left - right) / denominator


def _parse_utc_ns(value: str) -> int:
    if not isinstance(value, str):
        raise MetricsError(f"invalid UTC resource-sample timestamp: {value!r}")
    try:
        parsed = dt.datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError as exc:
        raise MetricsError(f"invalid UTC resource-sample timestamp: {value!r}") from exc
    if parsed.tzinfo is None:
        raise MetricsError(f"UTC resource-sample timestamp lacks timezone: {value!r}")
    return int(parsed.timestamp() * 1_000_000_000)


def analyze_resource_clock_samples(
    path: Path,
    *,
    expected_sample_count: int,
    tolerance_fraction: float = CLOCK_CONSISTENCY_TOLERANCE_FRACTION,
) -> dict[str, object]:
    """Audit paired CLOCK_MONOTONIC/CLOCK_REALTIME resource-loop samples."""

    if tolerance_fraction <= 0.0:
        raise MetricsError("clock consistency tolerance must be positive")
    if not path.is_file():
        raise MetricsError(f"missing resource clock samples: {path}")
    required = {
        "sample_index",
        "elapsed_seconds",
        "monotonic_ns",
        "realtime_ns",
        "utc",
    }
    samples: list[dict[str, object]] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or not required <= set(reader.fieldnames):
            raise MetricsError(
                "resource clock sample header is incomplete: "
                f"{reader.fieldnames!r}"
            )
        for line_number, row in enumerate(reader, start=2):
            try:
                sample_index = int(row["sample_index"])
                elapsed_seconds = float(row["elapsed_seconds"])
                monotonic_ns = int(row["monotonic_ns"])
                realtime_ns = int(row["realtime_ns"])
                utc = row["utc"]
                utc_ns = _parse_utc_ns(utc)
            except (KeyError, TypeError, ValueError, MetricsError) as exc:
                raise MetricsError(
                    f"invalid resource clock sample at {path}:{line_number}"
                ) from exc
            if (
                sample_index <= 0
                or not math.isfinite(elapsed_seconds)
                or elapsed_seconds < 0.0
                or monotonic_ns <= 0
                or realtime_ns <= 0
            ):
                raise MetricsError(
                    f"non-positive resource clock sample at {path}:{line_number}"
                )
            # UTC is rendered at microsecond precision from the exact realtime
            # nanoseconds; allow only that deterministic sub-microsecond loss.
            if abs(utc_ns - realtime_ns) > 2_000:
                raise MetricsError(
                    f"UTC/realtime mismatch in resource clock sample at "
                    f"{path}:{line_number}"
                )
            samples.append(
                {
                    "sample_index": sample_index,
                    "elapsed_seconds": elapsed_seconds,
                    "monotonic_ns": monotonic_ns,
                    "realtime_ns": realtime_ns,
                    "utc": utc,
                }
            )

    if len(samples) != expected_sample_count:
        raise MetricsError(
            "resource clock sample count disagrees with metrics: "
            f"{len(samples)} != {expected_sample_count}"
        )
    if not samples:
        raise MetricsError("resource clock sample file is empty")
    for expected_index, sample in enumerate(samples, start=1):
        if int(sample["sample_index"]) != expected_index:
            raise MetricsError("resource clock sample indices are not contiguous")

    jumps: list[dict[str, object]] = []
    hard_failures: list[str] = []
    warnings: list[str] = []
    for previous, current in zip(samples, samples[1:]):
        monotonic_delta = (
            int(current["monotonic_ns"]) - int(previous["monotonic_ns"])
        ) / 1_000_000_000.0
        realtime_delta = (
            int(current["realtime_ns"]) - int(previous["realtime_ns"])
        ) / 1_000_000_000.0
        relative_difference = _relative_difference(monotonic_delta, realtime_delta)
        classification: str | None = None
        if monotonic_delta <= 0.0:
            classification = "hard_failure"
            hard_failures.append(
                "CLOCK_MONOTONIC did not advance between resource samples "
                f"{previous['sample_index']} and {current['sample_index']}"
            )
        elif realtime_delta <= 0.0:
            classification = "hard_failure"
            hard_failures.append(
                "CLOCK_REALTIME moved backwards or stopped between resource samples "
                f"{previous['sample_index']} and {current['sample_index']}"
            )
        elif relative_difference > tolerance_fraction:
            classification = "warning"
            warnings.append(
                "positive CLOCK_REALTIME step differs from CLOCK_MONOTONIC by "
                f"{relative_difference:.9%} between samples "
                f"{previous['sample_index']} and {current['sample_index']}"
            )
        if classification is not None:
            jumps.append(
                {
                    "previous_sample_index": int(previous["sample_index"]),
                    "sample_index": int(current["sample_index"]),
                    "monotonic_delta_seconds": monotonic_delta,
                    "realtime_delta_seconds": realtime_delta,
                    "absolute_difference_seconds": abs(
                        monotonic_delta - realtime_delta
                    ),
                    "relative_difference": relative_difference,
                    "classification": classification,
                }
            )

    forward_jumps = [
        jump for jump in jumps if jump["classification"] == "warning"
    ]
    backwards = [
        jump for jump in jumps if jump["classification"] == "hard_failure"
    ]
    status = "quarantined" if hard_failures else "warning" if warnings else "passed"
    return {
        "path": str(path.resolve()),
        "sample_count": len(samples),
        "interval_count": max(0, len(samples) - 1),
        "tolerance_fraction": tolerance_fraction,
        "realtime_jump_count": len(jumps),
        "realtime_jump_detected": bool(jumps),
        "forward_realtime_jump_count": len(forward_jumps),
        "realtime_backward_or_stopped_count": len(backwards),
        "periodic_realtime_jump_detected": len(forward_jumps) >= 2,
        "jumps": jumps,
        "hard_failures": hard_failures,
        "warnings": warnings,
        "status": status,
    }


def assess_clock_consistency(
    metrics: dict[str, object],
    resource_clock_path: Path,
    *,
    tolerance_fraction: float = CLOCK_CONSISTENCY_TOLERANCE_FRACTION,
) -> dict[str, object]:
    """Compare GNU elapsed, monotonic elapsed, realtime elapsed, and samples."""

    hard_failures: list[str] = []
    warnings: list[str] = []
    clock_values: dict[str, float] = {}
    gnu_time = metrics.get("gnu_time")
    measurements = metrics.get("clock_measurements")
    if isinstance(gnu_time, dict):
        gnu_elapsed = gnu_time.get("elapsed_seconds")
        if (
            isinstance(gnu_elapsed, (int, float))
            and math.isfinite(float(gnu_elapsed))
            and float(gnu_elapsed) > 0.0
        ):
            clock_values["gnu_time_elapsed_seconds"] = float(gnu_elapsed)
        else:
            hard_failures.append(
                "GNU elapsed is missing, non-finite, or non-positive"
            )
    else:
        hard_failures.append("GNU time measurements are missing")
    if isinstance(measurements, dict):
        endpoint_specs = (
            (
                "monotonic",
                "monotonic_started_ns",
                "monotonic_finished_ns",
                "monotonic_elapsed_seconds",
            ),
            (
                "realtime",
                "realtime_started_ns",
                "realtime_finished_ns",
                "realtime_elapsed_seconds",
            ),
        )
        for label, started_field, finished_field, elapsed_field in endpoint_specs:
            started = measurements.get(started_field)
            finished = measurements.get(finished_field)
            declared = measurements.get(elapsed_field)
            if (
                not isinstance(started, int)
                or not isinstance(finished, int)
                or finished <= started
                or not isinstance(declared, (int, float))
            ):
                hard_failures.append(f"{label} clock endpoints are invalid")
                continue
            derived = (finished - started) / 1_000_000_000.0
            if _relative_difference(derived, float(declared)) > 1e-12:
                hard_failures.append(
                    f"{label} elapsed disagrees with its clock endpoints"
                )
        realtime_started = measurements.get("realtime_started_ns")
        realtime_finished = measurements.get("realtime_finished_ns")
        started_utc = measurements.get("started_utc")
        finished_utc = measurements.get("finished_utc")
        if (
            not isinstance(realtime_started, int)
            or not isinstance(realtime_finished, int)
            or not isinstance(started_utc, str)
            or not isinstance(finished_utc, str)
        ):
            hard_failures.append("UTC endpoints are missing or invalid")
        else:
            try:
                rendered_started_ns = _parse_utc_ns(started_utc)
                rendered_finished_ns = _parse_utc_ns(finished_utc)
            except MetricsError as exc:
                hard_failures.append(str(exc))
            else:
                if (
                    abs(rendered_started_ns - realtime_started) > 2_000
                    or abs(rendered_finished_ns - realtime_finished) > 2_000
                ):
                    hard_failures.append(
                        "UTC endpoints disagree with CLOCK_REALTIME"
                    )
        if (
            metrics.get("started_utc") != measurements.get("started_utc")
            or metrics.get("finished_utc") != measurements.get("finished_utc")
        ):
            hard_failures.append(
                "top-level UTC endpoints disagree with clock measurements"
            )
        runner_wall = metrics.get("runner_wall_seconds")
        monotonic_declared = measurements.get("monotonic_elapsed_seconds")
        if (
            not isinstance(runner_wall, (int, float))
            or not isinstance(monotonic_declared, (int, float))
            or _relative_difference(
                float(runner_wall), float(monotonic_declared)
            )
            > 1e-12
        ):
            hard_failures.append(
                "legacy runner wall value disagrees with monotonic elapsed"
            )
        for source, field in (
            ("monotonic_elapsed_seconds", "monotonic_elapsed_seconds"),
            ("realtime_utc_elapsed_seconds", "realtime_elapsed_seconds"),
        ):
            value = measurements.get(field)
            if (
                isinstance(value, (int, float))
                and math.isfinite(float(value))
                and float(value) > 0.0
            ):
                clock_values[source] = float(value)
            else:
                hard_failures.append(
                    f"{source} is missing, non-finite, or non-positive"
                )
    else:
        hard_failures.append("monotonic/realtime clock measurements are missing")

    comparisons: list[dict[str, object]] = []
    names = tuple(clock_values)
    for left_index, left_name in enumerate(names):
        for right_name in names[left_index + 1 :]:
            left = clock_values[left_name]
            right = clock_values[right_name]
            relative_difference = _relative_difference(left, right)
            passed = relative_difference <= tolerance_fraction
            comparisons.append(
                {
                    "left": left_name,
                    "right": right_name,
                    "absolute_difference_seconds": abs(left - right),
                    "relative_difference": relative_difference,
                    "status": "passed" if passed else "quarantined",
                }
            )
            if not passed:
                warnings.append(
                    f"{left_name} and {right_name} differ by "
                    f"{relative_difference:.9%} (> {tolerance_fraction:.2%})"
                )

    sample_count = metrics.get("sample_count")
    resource_report: dict[str, object]
    try:
        if not isinstance(sample_count, int) or sample_count <= 0:
            raise MetricsError(f"invalid metrics sample count: {sample_count!r}")
        resource_report = analyze_resource_clock_samples(
            resource_clock_path,
            expected_sample_count=sample_count,
            tolerance_fraction=tolerance_fraction,
        )
    except MetricsError as exc:
        resource_report = {
            "path": str(resource_clock_path.resolve()),
            "status": "quarantined",
            "error": str(exc),
            "realtime_jump_detected": None,
            "periodic_realtime_jump_detected": None,
        }
        hard_failures.append(str(exc))
    else:
        resource_hard = resource_report.get("hard_failures")
        resource_warnings = resource_report.get("warnings")
        if isinstance(resource_hard, list):
            hard_failures.extend(str(value) for value in resource_hard)
        if isinstance(resource_warnings, list):
            warnings.extend(str(value) for value in resource_warnings)

    if len(clock_values) != 3:
        hard_failures.append("all three elapsed clocks are required")
    status = (
        "quarantined"
        if hard_failures
        else "warning"
        if warnings
        else "passed"
    )
    return {
        "schema": "ramag.clock-consistency.v1",
        "tolerance_fraction": tolerance_fraction,
        "clock_values_seconds": clock_values,
        "pairwise_comparisons": comparisons,
        "resource_samples": resource_report,
        "hard_failures": hard_failures,
        "warnings": warnings,
        "reasons": [*hard_failures, *warnings],
        "status": status,
    }


def write_json(path: Path, value: object) -> None:
    with path.open("x", encoding="utf-8", newline="\n") as handle:
        json.dump(value, handle, indent=2, sort_keys=True, ensure_ascii=False)
        handle.write("\n")


def _normalize_environment_overrides(
    values: dict[str, str] | None,
) -> dict[str, str]:
    normalized: dict[str, str] = {}
    for raw_key, raw_value in (values or {}).items():
        key = str(raw_key)
        value = str(raw_value)
        if not key or "=" in key or "\x00" in key or "\x00" in value:
            raise MetricsError(f"invalid environment override: {key!r}")
        normalized[key] = value
    return dict(sorted(normalized.items()))


def _process_group_exists(process_group_id: int) -> bool:
    """Return whether the runner-owned process group still has a member."""

    try:
        os.killpg(process_group_id, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def _terminate_process_group(
    process: subprocess.Popen[bytes],
    *,
    grace_seconds: float = TIMEOUT_TERMINATION_GRACE_SECONDS,
) -> dict[str, object]:
    """Terminate only ``process``'s new-session process group, then reap it."""

    if grace_seconds < 0.0 or not math.isfinite(grace_seconds):
        raise MetricsError("termination grace must be finite and non-negative")
    process_group_id = process.pid
    sigterm_sent = False
    sigkill_sent = False
    try:
        os.killpg(process_group_id, signal.SIGTERM)
        sigterm_sent = True
    except ProcessLookupError:
        pass

    grace_started_ns = time.monotonic_ns()
    grace_deadline_ns = grace_started_ns + int(grace_seconds * 1_000_000_000)
    while True:
        # poll() reaps the /usr/bin/time wrapper if it has already exited.  That
        # prevents a wrapper zombie from making killpg(..., 0) look live while
        # still allowing us to detect an independently surviving descendant.
        process.poll()
        if not _process_group_exists(process_group_id):
            break
        now_ns = time.monotonic_ns()
        if now_ns >= grace_deadline_ns:
            try:
                os.killpg(process_group_id, signal.SIGKILL)
                sigkill_sent = True
            except ProcessLookupError:
                pass
            break
        remaining_seconds = (grace_deadline_ns - now_ns) / 1_000_000_000.0
        time.sleep(min(0.05, remaining_seconds))

    if process.poll() is None:
        process.wait()
    if sigkill_sent:
        # Give init/subreapers a short bounded opportunity to reap any killed
        # descendants before recording the final process-group state.
        reap_deadline_ns = time.monotonic_ns() + 1_000_000_000
        while _process_group_exists(process_group_id) and time.monotonic_ns() < reap_deadline_ns:
            time.sleep(0.01)
    grace_finished_ns = time.monotonic_ns()
    return {
        "process_group_id": process_group_id,
        "sigterm_sent": sigterm_sent,
        "sigkill_sent": sigkill_sent,
        "grace_started_monotonic_ns": grace_started_ns,
        "grace_finished_monotonic_ns": grace_finished_ns,
        "grace_elapsed_seconds": (
            grace_finished_ns - grace_started_ns
        )
        / 1_000_000_000.0,
        "process_group_alive_after_termination": _process_group_exists(
            process_group_id
        ),
        "process_exit_code": process.returncode,
    }


def run(
    output_dir: Path,
    command: list[str],
    sample_interval: float,
    environment_overrides: dict[str, str] | None = None,
    timeout_seconds: float | None = None,
) -> int:
    if output_dir.exists():
        raise MetricsError(f"refusing to overwrite run directory: {output_dir}")
    if not command:
        raise MetricsError("no command supplied")
    time_binary = Path("/usr/bin/time")
    if not time_binary.is_file() or not os.access(time_binary, os.X_OK):
        raise MetricsError("GNU /usr/bin/time is required")
    if sample_interval <= 0:
        raise MetricsError("sample interval must be positive")
    if timeout_seconds is not None and (
        not math.isfinite(timeout_seconds) or timeout_seconds <= 0.0
    ):
        raise MetricsError("timeout must be finite and positive")

    output_dir.mkdir(parents=True)
    (output_dir / "command.txt").write_text(
        shlex.join(command) + "\n", encoding="utf-8", newline="\n"
    )
    (output_dir / "command.json").write_text(
        json.dumps(command, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    stdout_path = output_dir / "stdout.txt"
    stderr_path = output_dir / "stderr.txt"
    time_path = output_dir / "time.raw.txt"
    resources_path = output_dir / "resources.tree.tsv"
    resource_clocks_path = output_dir / "resources.clock.tsv"
    started_monotonic_ns = time.monotonic_ns()
    started_realtime_ns = time.time_ns()
    started_utc = _utc_from_realtime_ns(started_realtime_ns)
    clock_ticks = os.sysconf("SC_CLK_TCK")
    previous: dict[int, tuple[int, float]] = {}
    peak_tree_rss = 0
    peak_tree_vms = 0
    max_total_threads = 0
    max_single_process_threads = 0
    max_aligner_tree_threads = 0
    max_aligner_single_process_threads = 0
    sample_count = 0
    timed_out = False
    timeout_signal_monotonic_ns: int | None = None
    timeout_elapsed_at_signal: float | None = None
    termination: dict[str, object] | None = None

    normalized_environment = _normalize_environment_overrides(environment_overrides)
    environment = os.environ.copy()
    environment["LC_ALL"] = "C"
    environment.update(normalized_environment)
    wrapped = [str(time_binary), "-v", "-o", str(time_path), "--", *command]
    with stdout_path.open("xb") as stdout, stderr_path.open("xb") as stderr, resources_path.open(
        "x", encoding="utf-8", newline="\n"
    ) as resources, resource_clocks_path.open(
        "x", encoding="utf-8", newline="\n"
    ) as resource_clocks:
        resources.write(
            "elapsed_seconds\tutc\tpid\tppid\tprocess\tthread_count\t"
            "cpu_percent\trss_bytes\tvms_bytes\n"
        )
        resource_clocks.write(
            "sample_index\telapsed_seconds\tmonotonic_ns\trealtime_ns\tutc\n"
        )
        process = subprocess.Popen(
            wrapped,
            stdout=stdout,
            stderr=stderr,
            env=environment,
            start_new_session=True,
        )
        interrupted = False
        try:
            while True:
                sample_monotonic_ns = time.monotonic_ns()
                sample_realtime_ns = time.time_ns()
                now = sample_monotonic_ns / 1_000_000_000.0
                elapsed = (
                    sample_monotonic_ns - started_monotonic_ns
                ) / 1_000_000_000.0
                return_code = process.poll()
                if return_code is not None:
                    break
                if timeout_seconds is not None and elapsed > timeout_seconds:
                    timed_out = True
                    timeout_signal_monotonic_ns = time.monotonic_ns()
                    timeout_elapsed_at_signal = (
                        timeout_signal_monotonic_ns - started_monotonic_ns
                    ) / 1_000_000_000.0
                    termination = _terminate_process_group(process)
                    break
                timestamp = _utc_from_realtime_ns(sample_realtime_ns)
                sample_index = sample_count + 1
                resource_clocks.write(
                    f"{sample_index}\t{elapsed:.9f}\t{sample_monotonic_ns}\t"
                    f"{sample_realtime_ns}\t{timestamp}\n"
                )
                resource_clocks.flush()
                snapshot = _descendants(process.pid)
                tree_rss = 0
                tree_vms = 0
                total_threads = 0
                aligner_tree_threads = 0
                for pid, ppid, ticks in snapshot:
                    status = _status(pid)
                    if status is None:
                        continue
                    old = previous.get(pid)
                    cpu_percent = 0.0
                    if old is not None and now > old[1] and ticks >= old[0]:
                        cpu_percent = (ticks - old[0]) / clock_ticks / (now - old[1]) * 100.0
                    previous[pid] = (ticks, now)
                    rss = int(status["VmRSS"])
                    vms = int(status["VmSize"])
                    threads = int(status["threads"])
                    tree_rss += rss
                    tree_vms += vms
                    total_threads += threads
                    max_single_process_threads = max(max_single_process_threads, threads)
                    if pid != process.pid:
                        aligner_tree_threads += threads
                        max_aligner_single_process_threads = max(
                            max_aligner_single_process_threads, threads
                        )
                    name = str(status["name"]).replace("\t", " ").replace("\n", " ")
                    resources.write(
                        f"{elapsed:.6f}\t{timestamp}\t{pid}\t{ppid}\t{name}\t"
                        f"{threads}\t{cpu_percent:.3f}\t{rss}\t{vms}\n"
                    )
                resources.flush()
                sample_count += 1
                peak_tree_rss = max(peak_tree_rss, tree_rss)
                peak_tree_vms = max(peak_tree_vms, tree_vms)
                max_total_threads = max(max_total_threads, total_threads)
                max_aligner_tree_threads = max(
                    max_aligner_tree_threads, aligner_tree_threads
                )
                wait_seconds = sample_interval
                if timeout_seconds is not None:
                    wait_seconds = min(
                        wait_seconds,
                        max(0.0, timeout_seconds - elapsed),
                    )
                try:
                    process.wait(timeout=wait_seconds)
                except subprocess.TimeoutExpired:
                    pass
                else:
                    break
        except KeyboardInterrupt:
            interrupted = True
            termination = _terminate_process_group(process)
        process_exit_code = process.returncode
        exit_code = 124 if timed_out else process_exit_code
        if interrupted and exit_code == 0:
            exit_code = 130

    finished_monotonic_ns = time.monotonic_ns()
    finished_realtime_ns = time.time_ns()
    finished_utc = _utc_from_realtime_ns(finished_realtime_ns)
    monotonic_elapsed = (
        finished_monotonic_ns - started_monotonic_ns
    ) / 1_000_000_000.0
    realtime_elapsed = (
        finished_realtime_ns - started_realtime_ns
    ) / 1_000_000_000.0
    parsed_time: dict[str, object] | None = None
    parse_error: str | None = None
    try:
        parsed_time = parse_gnu_time(time_path)
    except (OSError, MetricsError, ValueError) as error:
        parse_error = str(error)

    binary = Path(command[0])
    binary_sha = sha256_file(binary) if binary.is_file() else None
    timeout_evidence = {
        "enabled": timeout_seconds is not None,
        "clock": "CLOCK_MONOTONIC",
        "limit_seconds": float(timeout_seconds) if timeout_seconds is not None else 0.0,
        "exceeded": timed_out,
        "elapsed_at_signal_seconds": timeout_elapsed_at_signal,
        "signal_monotonic_ns": timeout_signal_monotonic_ns,
        "sigterm_sent": bool(termination and termination.get("sigterm_sent"))
        if timed_out
        else False,
        "sigkill_sent": bool(termination and termination.get("sigkill_sent"))
        if timed_out
        else False,
        "termination_grace_seconds": TIMEOUT_TERMINATION_GRACE_SECONDS,
        "deadline_overshoot_seconds": (
            max(0.0, timeout_elapsed_at_signal - float(timeout_seconds))
            if timed_out
            and timeout_elapsed_at_signal is not None
            and timeout_seconds is not None
            else None
        ),
        "process_group_id": process.pid,
        "process_group_alive_after_termination": (
            termination.get("process_group_alive_after_termination")
            if timed_out and termination is not None
            else None
        ),
        "termination": termination if timed_out else None,
    }
    status = "timed_out" if timed_out else "success" if exit_code == 0 else "failed"
    metrics = {
        "schema": "ramag.command-metrics.v3",
        "command": command,
        "environment_overrides": normalized_environment,
        "started_utc": started_utc,
        "finished_utc": finished_utc,
        "runner_wall_seconds": monotonic_elapsed,
        "clock_measurements": {
            "monotonic_clock": "CLOCK_MONOTONIC",
            "monotonic_started_ns": started_monotonic_ns,
            "monotonic_finished_ns": finished_monotonic_ns,
            "monotonic_elapsed_seconds": monotonic_elapsed,
            "realtime_clock": "CLOCK_REALTIME",
            "realtime_started_ns": started_realtime_ns,
            "realtime_finished_ns": finished_realtime_ns,
            "realtime_elapsed_seconds": realtime_elapsed,
            "started_utc": started_utc,
            "finished_utc": finished_utc,
        },
        "exit_code": exit_code,
        "process_exit_code": process_exit_code,
        "sample_interval_seconds": sample_interval,
        "sample_count": sample_count,
        "observed_process_tree_peak_rss_bytes": peak_tree_rss,
        "observed_process_tree_peak_vms_bytes": peak_tree_vms,
        "observed_max_total_threads": max_total_threads,
        "observed_max_single_process_threads": max_single_process_threads,
        "observed_max_aligner_tree_threads": max_aligner_tree_threads,
        "observed_max_aligner_single_process_threads": max_aligner_single_process_threads,
        "binary_path": str(binary.resolve()) if binary.exists() else command[0],
        "binary_sha256": binary_sha,
        "gnu_time": parsed_time,
        "gnu_time_parse_error": parse_error,
        "status": status,
        "timeout": timeout_evidence,
    }
    write_json(output_dir / "metrics.json", metrics)
    (output_dir / "exit-code.txt").write_text(f"{exit_code}\n", encoding="ascii")
    if timed_out:
        write_json(output_dir / "timeout.json", timeout_evidence)
    marker = (
        "RUN_TIMED_OUT"
        if timed_out
        else "RUN_COMPLETE"
        if exit_code == 0
        else "RUN_FAILED"
    )
    (output_dir / marker).write_text(metrics["status"] + "\n", encoding="ascii")
    return int(exit_code)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--sample-interval", type=float, default=1.0)
    parser.add_argument(
        "--timeout-seconds",
        type=float,
        help="hard CLOCK_MONOTONIC deadline; equality is accepted",
    )
    parser.add_argument(
        "--env",
        action="append",
        default=[],
        metavar="KEY=VALUE",
        help="recorded environment override; may be repeated",
    )
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args(argv)
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    environment_overrides: dict[str, str] = {}
    for item in args.env:
        if "=" not in item:
            raise MetricsError(f"--env requires KEY=VALUE, got {item!r}")
        key, value = item.split("=", 1)
        if key in environment_overrides:
            raise MetricsError(f"duplicate --env key: {key}")
        environment_overrides[key] = value
    return run(
        args.output_dir,
        command,
        args.sample_interval,
        environment_overrides,
        args.timeout_seconds,
    )


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except MetricsError as error:
        print(f"run_with_metrics: {error}", file=sys.stderr)
        raise SystemExit(2)

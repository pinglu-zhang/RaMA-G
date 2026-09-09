#!/usr/bin/env python3
"""Summarize current exit/output/log benchmark runs; never infer biological F1."""
import argparse, gzip, json, pathlib, re, statistics

from pathlib import Path
import math
from run_with_metrics import assess_clock_consistency
ExperimentError=ValueError
ClockQuarantineError=ValueError
METRICS_SCHEMA="ramag.command-metrics.v3"
CLOCK_CONSISTENCY_TOLERANCE_FRACTION=0.05
def load_json(path):return json.loads(path.read_text())
def write_json_exclusive(path,value):
    with path.open('x') as h:json.dump(value,h,indent=2)

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


def catalog(path):
    path=pathlib.Path(path)
    with path.open('rb') as h: compressed=h.read(2)==b'\x1f\x8b'
    result={}; name=None
    with (gzip.open(path,'rt') if compressed else path.open()) as h:
        for line in h:
            if line.startswith('>'):
                name=line[1:].split()[0]
                if name in result:raise ValueError('duplicate FASTA identifier')
                result[name]=0
            elif line.strip():
                if name is None:raise ValueError('sequence before header')
                result[name]+=len(line.strip())
    if not result or any(n==0 for n in result.values()):raise ValueError('empty FASTA')
    return result

def union(intervals):
    total=0; end=0
    for a,b in sorted(intervals):
        total+=max(0,b-max(a,end));end=max(end,b)
    return total

def paf_coverage(path,refs,queries):
    ri={k:[] for k in refs};qi={k:[] for k in queries};count=0
    with pathlib.Path(path).open() as h:
        for line in h:
            f=line.rstrip('\n').split('\t')
            if len(f)<12:raise ValueError('invalid PAF columns')
            q,ql,qs,qe,strand,r,rl,rs,re_,matches,columns,mapq=f[:12]
            ql,qs,qe,rl,rs,re_,matches,columns,mapq=map(int,[ql,qs,qe,rl,rs,re_,matches,columns,mapq])
            if q not in queries or r not in refs or ql!=queries[q] or rl!=refs[r]:raise ValueError('PAF/input catalog mismatch')
            if not(0<=qs<qe<=ql and 0<=rs<re_<=rl and strand in ('+','-') and 0<=matches<=columns and 0<=mapq<=255):raise ValueError('invalid PAF coordinates')
            cg=next((v[5:] for v in f[12:] if v.startswith('cg:Z:')),None)
            if cg is None:raise ValueError('PAF lacks CIGAR')
            ops=re.findall(r'(\d+)([=XID])',cg)
            if ''.join(n+op for n,op in ops)!=cg or any(int(n)==0 for n,op in ops):raise ValueError('invalid CIGAR')
            if sum(int(n) for n,op in ops if op!='I')!=re_-rs or sum(int(n) for n,op in ops if op!='D')!=qe-qs:raise ValueError('CIGAR span mismatch')
            ri[r].append((rs,re_));qi[q].append((qs,qe));count+=1
    if count==0:raise ValueError('empty PAF')
    rb=sum(map(union,ri.values()));qb=sum(map(union,qi.values()))
    return dict(records=count,reference_covered_bases=rb,query_covered_bases=qb,reference_coverage=rb/sum(refs.values()),query_coverage=qb/sum(queries.values()))

def summarize(root):
    root=pathlib.Path(root);config=json.loads((root/'config.json').read_text())
    refs=catalog(config['reference']);queries=catalog(config['query']);rows=[]
    for i in range(1,config['repetitions']+1):
        d=root/f'run-{i}';metrics=validate_accepted_metrics(d/'metrics.json')
        if metrics.get('exit_code')!=0 or metrics.get('status')!='success':raise ValueError(f'run {i} did not succeed')
        logs=list((d/'work').glob('runs/*/run.log'))
        if len(logs)!=1:raise ValueError('expected one run log')
        text=logs[0].read_text()
        if 'status=success exit_code=0' not in text or 'status=failed' in text or 'status=interrupted' in text:raise ValueError('run log does not describe a successful run')
        action='loaded' if config.get('reference_index') else 'built'
        if 'sufkit.index.action='+action not in text:raise ValueError('index action mismatch')
        row=paf_coverage(d/'result.paf',refs,queries)
        row.update(trial=i,wall_seconds=metrics['runner_wall_seconds'],peak_tree_rss_bytes=metrics['observed_process_tree_peak_rss_bytes'],gnu_time=metrics['gnu_time'],stage_seconds={k:float(v) for k,v in re.findall(r'timing stage=(\S+) seconds=([\d.eE+-]+)',text)})
        rows.append(row)
    # Compare complete output bytes in bounded chunks, without digest computation.
    for i in range(2,len(rows)+1):
        with (root/'run-1/result.paf').open('rb') as a,(root/f'run-{i}/result.paf').open('rb') as b:
            while True:
                x=a.read(1048576);y=b.read(1048576)
                if x!=y:raise ValueError('PAF differs between repetitions')
                if not x:break
    result={'status':'success','definition':'alignment span union / complete input length','f1':None,'runs':rows,'median_wall_seconds':statistics.median(x['wall_seconds'] for x in rows)}
    (root/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    lines=['# RaMA-G benchmark','', 'Real data has no validated base-pair truth; F1 is not reported.','', '| Trial | Wall seconds | Peak tree GiB | Reference coverage | Query coverage | Records |','|---|---:|---:|---:|---:|---:|']
    for x in rows:lines.append(f"| {x['trial']} | {x['wall_seconds']:.3f} | {x['peak_tree_rss_bytes']/2**30:.3f} | {x['reference_coverage']:.6%} | {x['query_coverage']:.6%} | {x['records']} |")
    lines+=['','Index loading/validation is included in wall time. Cache state and competing server load are uncontrolled unless recorded externally.','Use the separate explicit-plan harness for MUMmer4; different seed/filter/index policies require qualified comparison.']
    (root/'comparison.md').write_text('\n'.join(lines)+'\n');return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--result-root',required=True);a=p.parse_args()
    summarize(a.result_root)

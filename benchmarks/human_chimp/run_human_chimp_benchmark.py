#!/usr/bin/env python3
"""Portable serial runner for the current RaMA-G exit/output/log contract.

Historical fixed-server CLI and resume bundles are intentionally not accepted.
This script never starts work merely by being imported.
"""
import argparse,json,pathlib,sys
from run_with_metrics import run
from summarize_comparison import summarize

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary',type=pathlib.Path,required=True)
    p.add_argument('--reference',type=pathlib.Path,required=True)
    p.add_argument('--query',type=pathlib.Path,required=True)
    p.add_argument('--reference-index',type=pathlib.Path)
    p.add_argument('--result-root',type=pathlib.Path,required=True)
    p.add_argument('--threads',type=int,default=16)
    p.add_argument('--repetitions',type=int,default=1)
    p.add_argument('--sample-interval',type=float,default=1.0)
    a=p.parse_args()
    if a.threads<1 or a.repetitions<1 or a.sample_interval<=0:p.error('counts and interval must be positive')
    for path in [a.binary,a.reference,a.query]+([a.reference_index] if a.reference_index else []):
        if not path.is_file():p.error(f'missing file: {path}')
    a.result_root.mkdir(parents=True,exist_ok=False)
    config={k:str(v.resolve()) if isinstance(v,pathlib.Path) else v for k,v in vars(a).items()}
    (a.result_root/'config.json').write_text(json.dumps(config,indent=2)+'\n')
    for i in range(1,a.repetitions+1):
        d=(a.result_root/f'run-{i}').resolve()
        cmd=[str(a.binary.resolve()),'align','--reference',str(a.reference.resolve()),'--query',str(a.query.resolve()),'--output',str(d/'result.paf'),'--work-dir',str(d/'work'),'--threads',str(a.threads),'--seed-mode','mumreference','--selection-mode','one-to-one','--progress','off']
        if a.reference_index:cmd+=['--reference-index',str(a.reference_index.resolve())]
        code=run(d,cmd,a.sample_interval)
        if code:return code
    summarize(a.result_root)
    return 0
if __name__=='__main__':sys.exit(main())

#!/usr/bin/env python3
"""Exit/output/log contract, gzip, persistence and rollback integration gates."""
import argparse, gzip, os, pathlib, random, signal, subprocess, tempfile, time, zlib
p=argparse.ArgumentParser();p.add_argument('--binary',required=True);p.add_argument('--root',required=True);a=p.parse_args()
root=pathlib.Path(a.root);root.mkdir(parents=True,exist_ok=True)
root=pathlib.Path(tempfile.mkdtemp(prefix='contract-',dir=root))
rng=random.Random(17031)
bases=''.join(rng.choice('ACGT') for _ in range(8000))
ref=root/'reference.fa';query=root/'query.fa'
ref.write_text('>ref description  preserved\n'+bases+'\n>other\nACGTNNNNACGT\n')
q=bases[:2500]+bases[2501:5200]+'ACT'+bases[5200:]
query.write_text('>query\n'+q+'\n')
checks=[]
def run(name, extra=(), expected=0, reference=ref, qry=query, env=None, binary=None, formats='sam,paf,delta,maf,chain'):
 d=root/name;d.mkdir();prefix=d/'result'
 cmd=[binary or a.binary,'align','--reference',str(reference),'--query',str(qry),'--output-prefix',str(prefix),'--formats',formats,'--work-dir',str(d/'work'),'--min-match','8','--min-cluster','20','--progress','off',*map(str,extra)]
 r=subprocess.run(cmd,env={**os.environ,**(env or {})},capture_output=True,text=True,timeout=30)
 (d/'stdout').write_text(r.stdout);(d/'stderr').write_text(r.stderr)
 assert r.returncode==expected,(name,r.returncode,expected,r.stderr)
 if not binary:
  assert not list(d.rglob('*.manifest.json')) and not list(d.rglob('*.complete')) and not list(d.rglob('failure.json')),name
  assert not r.stdout,name
 if expected==0:
  for f in formats.split(','):assert (d/('result.'+f)).stat().st_size>0,(name,f)
  if not binary:assert 'status=success exit_code=0' in logs(d),name
 else:
  assert not any((d/('result.'+f)).exists() for f in formats.split(',')),name
 checks.append(name);return d
def logs(d):return '\n'.join(x.read_text() for x in d.rglob('run.log'))
def payload(d,f):
 b=(d/('result.'+f)).read_bytes()
 if f=='sam':b=b'\n'.join(x for x in b.splitlines() if not x.startswith(b'@PG'))
 if f=='delta':b=b'\n'.join(b.splitlines()[1:])
 return b
def same(x,y):
 for f in ['sam','paf','delta','maf','chain']:assert payload(x,f)==payload(y,f),(x,y,f)
idx=root/'saved.sufidx'
e=run('ephemeral');s=run('save',['--save',idx]);l=run('load',['--reference-index',idx]);same(e,s);same(s,l)
assert idx.is_file() and not pathlib.Path(str(idx)+'.manifest.json').exists() and not pathlib.Path(str(idx)+'.complete').exists()
assert 'sufkit.index.action=loaded' in logs(l) and 'stage=index-build' not in logs(l)
assert logs(s).count('event=stage stage=index-build ')==1
for mode in ['fast','mumreference','mum','smem','maxmatch']:
 previous=None
 for threads in [1,4]:
  tag=f'{mode}-linked-default-{threads}';args=['--seed-mode',mode,'--threads',threads]
  d=run(tag,args)
  if previous:same(previous,d)
  previous=d
  same(d,run('loaded-'+tag,[*args,'--reference-index',idx]))
for path in [ref,query]:path.with_suffix('.fa.gz').write_bytes(gzip.compress(path.read_bytes()))
for name,r,qf in [('gz-ref',ref.with_suffix('.fa.gz'),query),('gz-query',ref,query.with_suffix('.fa.gz')),('gz-both',ref.with_suffix('.fa.gz'),query.with_suffix('.fa.gz'))]:same(e,run(name,reference=r,qry=qf))
magic=root/'magic-input';magic.write_bytes(gzip.compress(ref.read_bytes()));same(e,run('magic',reference=magic))
concat=root/'concat.gz';b=ref.read_bytes();split=b.index(b'>other');concat.write_bytes(gzip.compress(b[:split])+gzip.compress(b[split:]));same(e,run('concat',reference=concat))
crlf=root/'crlf.fa';crlf.write_bytes(ref.read_bytes().replace(b'\n',b'\r\n'));same(e,run('crlf',reference=crlf))
malformed={'empty.fa':b'','blank.fa':b'\n>r\nACGT\n','fastq.fa':b'@r\nACGT\n+\nIIII\n','duplicate.fa':b'>r\nACGT\n>r\nACGT\n','no-bases.fa':b'>r\n>q\nACGT\n','space.fa':b'>r\nAC GT\n','mixed.fa':b'>r\nACGT\n@q\nACGT\n+\nIIII\n','fake.gz':b'>r\nACGT\n','truncated.gz':gzip.compress(b'>r\nACGT\n')[:-5]}
crc=bytearray(gzip.compress(b'>r\nACGT\n'));crc[-8]^=1;malformed['crc.gz']=crc
for name,data in malformed.items():
 f=root/name;f.write_bytes(data);run('bad-'+name,reference=f,expected=3)
bad=root/'bad.xz';bad.write_bytes(b'not supported');run('unsupported',reference=bad,expected=4)
changed=root/'changed.fa';changed.write_text(ref.read_text().replace(bases[100:120],'N'+bases[101:120],1));run('wrong-ref',['--reference-index',idx],reference=changed,expected=5)
damaged=root/'damaged.sufidx';bb=bytearray(idx.read_bytes());bb[len(bb)//2]^=1;damaged.write_bytes(bb);run('corrupt-load',['--reference-index',damaged],expected=5)
truncated=root/'truncated.sufidx';truncated.write_bytes(idx.read_bytes()[:100]);run('truncated-index',['--reference-index',truncated],expected=5)
for name,args in [('missing',['--save']),('duplicate-save',['--save',root/'d1','--save',root/'d2']),('mixed-save',['--save',root/'d3','--reference-index',idx]),('existing-save',['--save',idx]),('smem-zero',['--seed-mode','smem','--smem-min-occurrences','0']),('smem-misuse',['--smem-min-occurrences','2'])]:run(name,args,expected=2)
for phase in ['save','validation','publication']:
 target=root/('failure-'+phase+'.sufidx');run('index-failure-'+phase,['--save',target],env={'RAMAG_TEST_FAIL_INDEX_STAGE':phase},expected=7);assert not target.exists()
for f in ['sam','paf','delta','maf','chain']:
 for hook in ['WRITER','VALIDATOR']:run(f'failure-{hook}-{f}',env={'RAMAG_TEST_FAIL_'+hook:f},expected=7)
retained=root/'retained.sufidx';d=run('after-save',['--save',retained],env={'RAMAG_TEST_FAIL_WRITER':'paf'},expected=7)
assert retained.is_file() and 'retained_reference_index='+str(retained) in logs(d);same(e,run('retained-load',['--reference-index',retained]))
# Standalone command publishes only the index; loading ignores historical companions.
standalone=root/'standalone.sufidx'
r=subprocess.run([a.binary,'index','--reference',str(ref),'--output',str(standalone),'--work-dir',str(root/'index-work'),'--threads','1'],capture_output=True,text=True)
assert r.returncode==0,(r.returncode,r.stderr)
assert standalone.exists() and not pathlib.Path(str(standalone)+'.complete').exists()
pathlib.Path(str(standalone)+'.manifest.json').write_text('obsolete unrelated content')
pathlib.Path(str(standalone)+'.complete').write_text('obsolete unrelated content')
same(e,run('old-companions-ignored',['--reference-index',standalone]))
assert pathlib.Path(str(standalone)+'.manifest.json').read_text()=='obsolete unrelated content'
(pathlib.Path(a.root)/'result.sam').write_bytes((e/'result.sam').read_bytes())
(root/'report.txt').write_text('PASS '+str(len(checks))+' scenarios\n'+'\n'.join(checks))
print('PASS',len(checks),'runtime contract scenarios',root)

for flag in ['--selection-mode','--link-gap-mode','--link-recovery-mode']:
 d=run('removed-'+flag[2:],[flag,'none'],expected=2)
 assert not any((d/('result.'+f)).exists() for f in ['sam','paf','delta','maf','chain'])
help_text=subprocess.check_output([a.binary,'align','--help'],text=True)
assert all(flag not in help_text for flag in ['--selection-mode','--link-gap-mode','--link-recovery-mode'])
batch=root/'batch'
cp=subprocess.run([a.binary,'batch','--reference',str(ref),'--query',str(query),
 '--output-dir',str(batch),'--reference-index',str(idx),'--threads','1','--progress','off'],capture_output=True,text=True,timeout=60)
assert cp.returncode==0,cp.stderr
assert (batch/'query/alignment.paf').is_file() and '\tsuccess\t0\t' in (batch/'batch.tsv').read_text()
print('Removed-option rejection and serial batch passed')

# Cooperative cancellation is exercised in both backend configurations.
if os.name == 'posix':
 for stage, sig, expected in [('input',signal.SIGINT,130),('chaining',signal.SIGTERM,143),
                              ('writer-paf',signal.SIGINT,130),('publication',signal.SIGTERM,143)]:
  d=root/('signal-'+stage);d.mkdir();err=d/'stderr';output=d/'result.paf'
  cmd=[a.binary,'align','--reference',str(ref),'--query',str(query),'--output',str(output),
       '--work-dir',str(d/'work'),'--threads','1','--progress','on']
  with err.open('w') as stream:
   child=subprocess.Popen(cmd,stdout=subprocess.DEVNULL,stderr=stream,
                          env={**os.environ,'RAMAG_TEST_PAUSE_STAGE':stage})
   try:
    deadline=time.monotonic()+10
    while 'stage='+stage not in err.read_text():
     assert child.poll() is None and time.monotonic()<deadline,(stage,err.read_text())
     time.sleep(.02)
    if stage=='writer-paf':
     child.send_signal(signal.SIGUSR1);deadline=time.monotonic()+5
     while 'event=snapshot' not in err.read_text():
      assert child.poll() is None and time.monotonic()<deadline,err.read_text()
      time.sleep(.02)
    child.send_signal(sig);assert child.wait(timeout=10)==expected,(stage,err.read_text())
   finally:
    if child.poll() is None:child.kill();child.wait()
  assert not output.exists(),stage
  assert 'status=interrupted' in logs(d),stage
 print('Four genome cancellation/publication stages and SIGUSR1 snapshot passed')

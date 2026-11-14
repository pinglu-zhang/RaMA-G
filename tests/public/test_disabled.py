"""A backend-OFF binary rejects centromere work before publishing output."""
from pathlib import Path
import subprocess, sys, tempfile
with tempfile.TemporaryDirectory(prefix='ramag-public-disabled-') as d:
 p=Path(d);r=p/'r.fa';q=p/'q.fa'
 r.write_text('>r\nACGTACGT\n');q.write_text('>q\nACGTACGT\n')
 cp=subprocess.run([sys.argv[1],'align','-c','--reference',str(r),'--query',str(q),
  '--output-prefix',str(p/'out')],capture_output=True,text=True,timeout=20)
 assert cp.returncode==5,(cp.returncode,cp.stderr)
 assert not list(p.glob('out.*'))
print('Disabled centromere request rejected without output')

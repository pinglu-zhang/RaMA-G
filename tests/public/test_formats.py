"""Independent genome coordinate/CIGAR checks for five output formats."""
from pathlib import Path
import random
import re
import subprocess
import sys
import tempfile

BINARY = sys.argv[1]
FORMATS = ('paf', 'sam', 'maf', 'delta', 'chain')
COMPLEMENT = str.maketrans('ACGTN', 'TGCAN')

def rc(sequence):
    return sequence.translate(COMPLEMENT)[::-1]

def operations(text):
    result = [(int(n), op) for n, op in re.findall(r'([1-9][0-9]*)([MIDNSHP=X])', text)]
    assert ''.join(str(n) + op for n, op in result) == text, text
    return result

def replay(reference, query, rb, qb, strand, ops):
    ri, qi = rb, qb
    direction = 1 if strand == '+' else -1
    pairs = set()
    for n, op in ops:
        if op in '=XM':
            for _ in range(n):
                assert 0 <= ri < len(reference) and 0 <= qi < len(query)
                a, b = reference[ri], query[qi]
                if direction < 0:
                    b = b.translate(COMPLEMENT)
                exact = a in 'ACGT' and a == b
                if op != 'M':
                    assert (op == '=') == exact, (ri, qi, op, a, b)
                pairs.add((ri, qi, strand))
                ri += 1
                qi += direction
        elif op == 'D' or op == 'N':
            ri += n
        elif op in 'ISH':
            qi += direction * n
        else:
            assert op == 'P'
    return ri, qi, pairs

def paf_pairs(path, reference, query):
    all_pairs = set()
    intervals = []
    for line in path.read_text().splitlines():
        f = line.split('\t')
        assert len(f) >= 12 and f[0] == 'q' and f[5] == 'r'
        assert int(f[1]) == len(query) and int(f[6]) == len(reference)
        tags = {x.split(':', 2)[0]: x.split(':', 2)[2] for x in f[12:]}
        ops = operations(tags['cg'])
        rb, re_, qb, qe = map(int, (f[7], f[8], f[2], f[3]))
        strand = f[4]
        rend, qend, pairs = replay(reference, query, rb, qb if strand == '+' else qe - 1, strand, ops)
        assert rend == re_ and qend == (qe if strand == '+' else qb - 1)
        assert int(f[9]) == sum(n for n, op in ops if op == '=')
        assert int(f[10]) == sum(n for n, _ in ops)
        assert int(tags['NM']) == sum(n for n, op in ops if op in 'XID')
        for ar, br, aq, bq in intervals:
            assert re_ <= ar or br <= rb, ('reference overlap', intervals, f)
            assert qe <= aq or bq <= qb, ('query overlap', intervals, f)
        intervals.append((rb, re_, qb, qe))
        all_pairs.update(pairs)
    return all_pairs

def sam_pairs(path, reference, query):
    result = set()
    for line in path.read_text().splitlines():
        if line.startswith('@'):
            continue
        f = line.split('\t')
        flag = int(f[1])
        if flag & 4:
            continue
        assert f[0] == 'q' and f[2] == 'r'
        strand = '-' if flag & 16 else '+'
        _, _, pairs = replay(reference, query, int(f[3]) - 1,
                             0 if strand == '+' else len(query) - 1,
                             strand, operations(f[5]))
        result.update(pairs)
    return result

def maf_pairs(path, reference, query):
    result, rows = set(), []
    def flush():
        if not rows:
            return
        assert len(rows) == 2
        a, b = rows
        assert a[1] == 'r' and b[1] == 'q' and a[4] == '+'
        assert int(a[5]) == len(reference) and int(b[5]) == len(query)
        ri = int(a[2])
        qi = int(b[2]) if b[4] == '+' else len(query) - int(b[2]) - 1
        step = 1 if b[4] == '+' else -1
        assert len(a[6]) == len(b[6])
        assert sum(x != '-' for x in a[6]) == int(a[3])
        assert sum(x != '-' for x in b[6]) == int(b[3])
        for x, y in zip(a[6], b[6]):
            assert x != '-' or y != '-'
            if x != '-':
                assert x == reference[ri]
            if y != '-':
                expected = query[qi] if step > 0 else query[qi].translate(COMPLEMENT)
                assert y == expected
            if x != '-' and y != '-':
                result.add((ri, qi, b[4]))
            ri += x != '-'
            qi += step * (y != '-')
        rows.clear()
    for line in path.read_text().splitlines() + ['']:
        if not line.strip():
            flush()
        elif line.startswith('s '):
            rows.append(line.split())
    return result

def delta_pairs(path, reference, query):
    lines = path.read_text().splitlines()
    assert len(lines) >= 2 and lines[1] == 'NUCMER'
    result = set()
    i = 2
    while i < len(lines):
        if lines[i].startswith('>'):
            names = lines[i][1:].split()
            assert names == ['r', 'q', str(len(reference)), str(len(query))]
            i += 1
            continue
        rstart, rend, qstart, qend, *_ = map(int, lines[i].split())
        ri, qi = rstart - 1, qstart - 1
        step = 1 if qstart <= qend else -1
        strand = '+' if step > 0 else '-'
        i += 1
        while True:
            delta = int(lines[i]);i += 1
            if delta == 0:
                break
            for _ in range(abs(delta) - 1):
                result.add((ri, qi, strand));ri += 1;qi += step
            if delta > 0:
                ri += 1
            else:
                qi += step
        while ri < rend:
            result.add((ri, qi, strand));ri += 1;qi += step
        assert qi == qend if step > 0 else qi == qend - 2
    return result

def chain_pairs(path, reference, query):
    result = set()
    ri = qi = None
    for line in path.read_text().splitlines():
        if not line:
            continue
        f = line.split()
        if f[0] == 'chain':
            assert f[2] == 'r' and f[7] == 'q' and f[4] == '+'
            assert int(f[3]) == len(reference) and int(f[8]) == len(query)
            strand = f[9]
            step = 1 if strand == '+' else -1
            ri = int(f[5])
            qi = int(f[10]) if step > 0 else len(query) - int(f[10]) - 1
        else:
            assert ri is not None and len(f) in (1, 3)
            for _ in range(int(f[0])):
                result.add((ri, qi, strand));ri += 1;qi += step
            if len(f) == 3:
                ri += int(f[1]);qi += step * int(f[2])
    return result

with tempfile.TemporaryDirectory(prefix='ramag-public-formats-') as tmp:
    root = Path(tmp)
    rng = random.Random(419)
    ref = ''.join(rng.choice('ACGT') for _ in range(2400))
    insertion = next(x for x in 'ACGT' if x not in (ref[699], ref[700])) * 7
    deletion = next(x for x in range(1300, 1400)
                    if ref[x - 1] != ref[x + 4] and ref[x] != ref[x + 5])
    query = ref[:700] + insertion + ref[700:deletion] + ref[deletion + 5:]
    query = query[:300] + next(x for x in 'ACGT' if x != query[300]) + query[301:]
    source = list(range(700)) + [None] * 7 + list(range(700, deletion)) + list(range(deletion + 5, len(ref)))
    for name, reference, sequence, ancestry in [
            ('identical', ref, ref, list(range(len(ref)))),
            ('forward', ref, query, source),
            ('reverse', ref, rc(query), list(reversed(source))),
            ('unknown', 'N' * 97, 'N' * 103, [None] * 103)]:
        d = root / name;d.mkdir()
        r, q = d / 'r.fa', d / 'q.fa'
        r.write_text('>r\n' + reference + '\n');q.write_text('>q\n' + sequence + '\n')
        cp = subprocess.run([BINARY, 'align', '--reference', str(r), '--query', str(q),
                             '--output-prefix', str(d / 'a'), '--formats', ','.join(FORMATS),
                             '--work-dir', str(d / 'work'), '--threads', '2', '--progress', 'off'],
                            capture_output=True, text=True, timeout=60)
        assert cp.returncode == 0, cp.stderr
        got = paf_pairs(d / 'a.paf', reference, sequence)
        for ext, parse in [('sam', sam_pairs), ('maf', maf_pairs), ('delta', delta_pairs), ('chain', chain_pairs)]:
            assert parse(d / ('a.' + ext), reference, sequence) == got, (name, ext)
        if name == 'unknown':
            assert not got
        else:
            strand = '-' if name == 'reverse' else '+'
            expected = {(rpos, qpos, strand) for qpos, rpos in enumerate(ancestry) if rpos is not None}
            assert got <= expected, (name, list(got - expected)[:3])
            assert got, name
            if name == 'identical':
                assert got == expected, (name, len(got), len(expected))
            # The product promises supported segments, not whole-input global
            # coverage. Edited examples test correctness of every emitted pair.
print('Five formats preserve independent pairs, consumption, strand and known artificial ancestry')

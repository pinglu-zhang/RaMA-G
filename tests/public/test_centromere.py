#!/usr/bin/env python3
"""Small synthetic CLI export regression; no real samples or benchmark runs.

Each format is parsed independently into original-coordinate paired bases and
path identities, then compared with the authoritative segments TSV. SAMtools
round-trip is additional validation when installed, not the coordinate oracle.
"""
import collections
import csv
import os
from pathlib import Path
import random
import re
import shutil
import signal
import subprocess
import sys
import tempfile
import time

BINARY = str(Path(sys.argv[1]).resolve())
SAMTOOLS = shutil.which('samtools') or ('/usr/bin/samtools' if Path('/usr/bin/samtools').is_file() else None)


def rc(sequence):
    return sequence.translate(str.maketrans('ACGTN', 'TGCAN'))[::-1]


def table(path):
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def cigar_ops(cigar, allowed='=XID'):
    ops = [(int(n), op) for n, op in re.findall(r'([1-9][0-9]*)([' + allowed + r'])', cigar)]
    assert ''.join(str(n) + op for n, op in ops) == cigar, cigar
    return ops


def canonical_cigar(ops):
    merged = []
    for n, op in ops:
        if merged and merged[-1][1] == op:
            merged[-1] = (merged[-1][0] + n, op)
        else:
            merged.append((n, op))
    return ''.join(str(n) + op for n, op in merged)


def metadata(tags):
    return tuple(tags[k] for k in ('si', 'pi', 'ri', 'bi', 'pr', 'sm', 'st', 'sc'))


def aux_tags(fields):
    result = {}
    for field in fields:
        key, kind, value = field.split(':', 2)
        assert key not in result, field
        assert kind in 'AifZHB', field
        result[key] = value
    return result


def paired_bases(rb, qb, qe, strand, ops, ref, query):
    ri = rb
    qi = qb if strand == '+' else qe - 1
    direction = 1 if strand == '+' else -1
    pairs = []
    for n, op in ops:
        if op in '=X':
            for _ in range(n):
                assert 0 <= ri < len(ref) and 0 <= qi < len(query)
                a, b = ref[ri], query[qi]
                if direction == -1:
                    b = rc(b)
                assert (a in 'ACGT' and a == b) == (op == '='), (ri, qi, op, a, b)
                pairs.append((ri, qi, strand))
                ri += 1
                qi += direction
        elif op == 'I':
            qi += direction * n
        elif op == 'D':
            ri += n
        else:
            raise AssertionError(op)
    assert qi == (qe if direction == 1 else qb - 1)
    return ri, tuple(pairs)


def signature(rb, re_, qb, qe, strand, ops, tags, ref, query):
    end, pairs = paired_bases(rb, qb, qe, strand, ops, ref, query)
    assert end == re_
    return (rb, re_, qb, qe, strand, canonical_cigar(ops), metadata(tags), pairs)


def expected_records(prefix, ref, query):
    summary = {x['field']: x['value'] for x in table(Path(str(prefix) + '.alignment-summary.tsv'))}
    rows = table(Path(str(prefix) + '.alignment-segments.tsv'))
    expected = []
    for row in rows:
        tags = {k: row[v] for k, v in (
            ('si', 'segment_id'), ('pi', 'path_id'), ('ri', 'region_id'),
            ('bi', 'branch_id'), ('pr', 'path_role'), ('sm', 'score_model'), ('sc', 'score'))}
        tags['st'] = summary['status']
        assert row['score_model'] == 'centrolign-3piece'
        score = 0
        for length, op in cigar_ops(row['cigar']):
            if op == '=':
                score += 20 * length
            elif op == 'X':
                score -= 80 * length
            else:
                score -= min(60 + 30 * length, 800 + 5 * length, 2500 + length)
        assert score == int(row['score']), row
        expected.append(signature(int(row['reference_begin']), int(row['reference_end']),
                                  int(row['query_begin']), int(row['query_end']), row['strand'],
                                  cigar_ops(row['cigar']), tags, ref, query))
    return collections.Counter(expected)


def parse_paf(path, ref, query):
    records = []
    for line in path.read_text().splitlines():
        fields = line.split('\t')
        assert len(fields) >= 12
        assert fields[0] == 'q' and fields[5] == 'r'
        assert int(fields[1]) == len(query) and int(fields[6]) == len(ref)
        assert fields[4] in ('+', '-') and fields[11] == '255'
        tags = aux_tags(fields[12:])
        assert 'tp' not in tags, 'Do not imply alternatives from path roles'
        ops = cigar_ops(tags['cg'])
        assert int(fields[9]) == sum(n for n, op in ops if op == '=')
        assert int(fields[10]) == sum(n for n, _ in ops)
        assert int(tags['NM']) == sum(n for n, op in ops if op != '=')
        assert tags['AS'] == tags['sc']
        records.append(signature(int(fields[7]), int(fields[8]), int(fields[2]), int(fields[3]),
                                 fields[4], ops, tags, ref, query))
    return collections.Counter(records)


def parse_sam_text(text, ref, query):
    records, flags, unmapped = [], [], 0
    assert '@HD\tVN:1.6' in text
    assert '@SQ\tSN:r\tLN:' + str(len(ref)) in text
    for line in text.splitlines():
        if line.startswith('@'):
            continue
        fields = line.split('\t')
        assert len(fields) >= 11
        flag = int(fields[1])
        if flag & 4:
            unmapped += 1
            assert fields[0] == 'q' and fields[2:4] == ['*', '0']
            assert fields[5] == '*' and fields[9] == query
            continue
        flags.append(flag)
        assert fields[0] == 'q' and fields[2] == 'r'
        assert fields[4] == '255' and fields[6:9] == ['*', '0', '0']
        assert fields[10] == '*'
        tags = aux_tags(fields[11:])
        assert 'SA' not in tags
        ops = cigar_ops(fields[5], '=XIDH')
        leading = ops.pop(0)[0] if ops and ops[0][1] == 'H' else 0
        trailing = ops.pop()[0] if ops and ops[-1][1] == 'H' else 0
        assert all(op != 'H' for _, op in ops)
        consumed_q = sum(n for n, op in ops if op != 'D')
        consumed_r = sum(n for n, op in ops if op != 'I')
        assert leading + consumed_q + trailing == len(query)
        strand = '-' if flag & 16 else '+'
        qb, qe = ((trailing, len(query) - leading) if strand == '-'
                  else (leading, len(query) - trailing))
        sequence = query[qb:qe]
        if strand == '-':
            sequence = rc(sequence)
        assert fields[9] == sequence
        assert int(tags['NM']) == sum(n for n, op in ops if op != '=')
        assert tags['AS'] == tags['sc']
        rb = int(fields[3]) - 1
        records.append(signature(rb, rb + consumed_r, qb, qe, strand, ops, tags, ref, query))
    if records:
        assert sum((flag & 0x900) == 0 for flag in flags) == 1
        assert not unmapped
    else:
        assert unmapped == 1
    return collections.Counter(records)


def parse_maf(path, ref, query):
    text = path.read_text()
    assert text.startswith('##maf version=1')
    records, tags, rows, score = [], None, [], None

    def flush():
        nonlocal tags, rows, score
        if not rows:
            return
        assert tags is not None and len(rows) == 2 and score is not None
        a, b = rows
        assert a[0] == 'ref.r' and b[0] == 'query.q'
        assert a[3] == '+' and b[3] in ('+', '-')
        assert int(a[4]) == len(ref) and int(b[4]) == len(query)
        rb, re_ = int(a[1]), int(a[1]) + int(a[2])
        qb, qe = int(b[1]), int(b[1]) + int(b[2])
        if b[3] == '-':
            qb, qe = len(query) - qe, len(query) - qb
        assert len(a[5]) == len(b[5])
        assert len(a[5].replace('-', '')) == int(a[2])
        assert len(b[5].replace('-', '')) == int(b[2])
        assert a[5].replace('-', '') == ref[rb:re_]
        assert b[5].replace('-', '') == (query[qb:qe] if b[3] == '+' else rc(query[qb:qe]))
        ops = []
        for x, y in zip(a[5], b[5]):
            assert x != '-' or y != '-'
            op = 'I' if x == '-' else 'D' if y == '-' else '=' if x in 'ACGT' and x == y else 'X'
            ops.append((1, op))
        assert int(score) == int(tags['sc'])
        records.append(signature(rb, re_, qb, qe, b[3], ops, tags, ref, query))
        tags, rows, score = None, [], None

    for line in text.splitlines() + ['']:
        if not line:
            flush()
        elif line.startswith('# ramag '):
            assert not rows
            tags = dict(item.split('=', 1) for item in line.split()[2:])
        elif line.startswith('a '):
            score = dict(item.split('=', 1) for item in line.split()[1:])['score']
        elif line.startswith('s '):
            fields = line.split()
            assert len(fields) == 7
            rows.append(fields[1:])
        else:
            assert line.startswith('#'), line
    return collections.Counter(records)


def validate(prefix, ref, query, selected=('paf', 'sam', 'maf')):
    expected = expected_records(prefix, ref, query)
    for fmt in ('paf', 'sam', 'maf'):
        path = Path(str(prefix) + '.' + fmt)
        if fmt not in selected:
            assert not path.exists(), path
            continue
        if fmt == 'paf':
            got = parse_paf(path, ref, query)
        elif fmt == 'sam':
            got = parse_sam_text(path.read_text(), ref, query)
            if SAMTOOLS:
                bam = Path(str(prefix) + '.roundtrip.bam')
                result = subprocess.run([SAMTOOLS, 'view', '-b', '-o', str(bam), str(path)],
                                        capture_output=True, text=True, errors='replace')
                assert result.returncode == 0, result.stderr
                result = subprocess.run([SAMTOOLS, 'view', '-h', str(bam)],
                                        capture_output=True, text=True, errors='replace')
                assert result.returncode == 0, result.stderr
                assert parse_sam_text(result.stdout, ref, query) == expected
        else:
            got = parse_maf(path, ref, query)
        assert got == expected, (fmt, list((got - expected).elements())[:1],
                                 list((expected - got).elements())[:1])
    return expected


with tempfile.TemporaryDirectory(prefix='ramag-export-cli-') as temporary:
    root = Path(temporary)
    rng = random.Random(918)
    ref = ''.join(rng.choice('ACGT') for _ in range(18000))
    ref = ref[:4500] + 'NNN' + ref[4503:]
    alternate = 'A' if ref[3301] != 'A' else 'T'
    query = ref[:3301] + alternate + ref[3302:6801] + 'AA' + ref[6801:12003] + ref[12006:]
    query = 'N' * 11 + query + 'N' * 13
    r, q = root / 'r.fa', root / 'q.fa'

    def inputs(reference, sequence):
        r.write_text('>r\n' + reference + '\n')
        q.write_text('>q\n' + sequence + '\n')

    inputs(ref, query)
    original = [p.read_bytes() for p in (r, q)]

    def command(name, extra=(), explicit=None, reference=r, query_path=q, progress='off'):
        args = [BINARY, 'align', '-c', '--reference', str(reference), '--query', str(query_path),
                '--progress', progress, '--threads', '1']
        if explicit is None:
            args += ['--output-prefix', str(root / name)]
        else:
            for fmt in explicit:
                args += ['--output', str(root / (name + '.' + fmt))]
        return args + list(extra)

    def run(name, env=None, **kwargs):
        return subprocess.run(command(name, **kwargs), capture_output=True, text=True, env=env,
                              timeout=30)

    result = run('default')
    assert result.returncode == 0, result.stderr
    assert (root / '.ramag-work' / 'runs').is_dir()
    expected = validate(root / 'default', ref, query)
    assert expected and all(item[4] == '+' for item in expected)
    assert any('X' in item[5] for item in expected)
    assert any('I' in item[5] for item in expected)
    assert any('D' in item[5] for item in expected)
    assert any(item[2] > 0 and item[3] < len(query) for item in expected)

    for name in ('orientation-local-direction-repair-intervals', 'direction-initial-anchors',
                 'directed-regions', 'initial-anchors', 'final-anchors', 'alignment-gaps'):
        assert table(root / ('default.' + name + '.tsv')), name
    for name in ('copy-anchors', 'copy-ambiguities', 'endpoint-reviews', 'block-reviews',
                 'cross-region-reviews', 'joint-path-reviews', 'unplanned-core-reviews'):
        assert (root / ('default.' + name + '.tsv')).is_file(), name
    status = {row['field']: row['value'] for row in table(root / 'default.directed-run.tsv')}
    assert status['orientation_model'] == 'orientation-evidence-v9'
    assert status['wfa_calls_scope'] == 'anchor_stage' and status['wfa_calls'] == '0'
    assert status['stitch_status'] == 'complete'

    result = run('selected', extra=('--formats', 'paf'))
    assert result.returncode == 0, result.stderr
    assert validate(root / 'selected', ref, query, ('paf',)) == expected
    result = run('explicit', explicit=('paf', 'sam', 'maf'))
    assert result.returncode == 0, result.stderr
    assert validate(root / 'explicit', ref, query) == expected
    result = run('config', extra=('--print-effective-config',))
    assert result.returncode == 0, result.stderr
    config = dict(line.split('=', 1) for line in result.stdout.splitlines() if '=' in line)
    assert config['centromere_stage'] == 'complete-alignment'
    assert config['orientation'] == 'predicted'
    assert config['stitch_score_model'] == 'centrolign-3piece'
    assert config['unplanned_core_recovery'] in ('initial-anchors-and-local-search', 'initial-anchors-only')
    assert set(config['formats'].split(',')) == {'sam', 'paf', 'maf'}
    assert not list(root.glob('config.*'))

    removed = (
        '--centromere-fixed-directions', '--centromere-copy-anchors', '--centromere-stitch',
        '--centromere-stitch-scoring', '--centromere-refine-endpoints', '--centromere-refine-blocks',
        '--centromere-refine-cross-regions', '--centromere-refine-joint-paths',
        '--centromere-initial-anchor-evidence', '--centromere-search-unplanned-cores')
    help_text = subprocess.check_output([BINARY, 'align', '--help'], text=True)
    for flag in removed:
        result = run('removed', extra=(flag,))
        assert result.returncode == 2 and 'unknown option: ' + flag in result.stderr
        assert flag not in help_text
    assert not list(root.glob('removed.*'))
    for fmt in ('delta', 'chain'):
        assert run('bad-' + fmt, extra=('--formats', fmt)).returncode == 2
        assert run('bad-explicit-' + fmt, explicit=(fmt,)).returncode == 2
    assert run('conflict', extra=('--output', str(root / 'conflict.sam'))).returncode == 2
    assert run('conflict-explicit', explicit=('paf',), extra=('--formats', 'sam')).returncode == 2
    assert run('invalid-seeds', extra=('--seed-mode', 'fast')).returncode == 2

    for suffix in ('paf', 'sam', 'maf', 'orientation-local-direction-repair-intervals.tsv',
                   'direction-initial-anchors.tsv'):
        name = 'exists-' + suffix.split('.')[0]
        protected = root / (name + '.' + suffix)
        protected.write_text('preexisting artifact\n')
        result = run(name)
        assert result.returncode != 0, (suffix, result.stderr)
        assert list(root.glob(name + '.*')) == [protected]
        assert protected.read_text() == 'preexisting artifact\n'
    collision = root / 'input-collision.paf'
    collision.write_bytes(r.read_bytes())
    result = run('input-collision', reference=collision)
    assert result.returncode != 0, result.stderr
    assert list(root.glob('input-collision.*')) == [collision]
    assert collision.read_bytes() == r.read_bytes()
    for stage in ('stitch-paf', 'stitch-sam', 'stitch-maf', 'stitch-segments',
                  'directed-final-anchors', 'endpoint-review', 'block-review',
                  'cross-region-review', 'direction-evidence'):
        result = run('failure', env=dict(os.environ, RAMAG_TEST_FAIL_WRITER=stage))
        assert result.returncode == 7, (stage, result.stderr)
        assert not list(root.glob('failure.*')), (stage, list(root.glob('failure.*')))
    assert [p.read_bytes() for p in (r, q)] == original

    multi = root / 'multi.fa'
    multi.write_bytes(r.read_bytes() + b'>extra\nACGT\n')
    result = run('multiple-records', reference=multi)
    assert result.returncode == 3 and not list(root.glob('multiple-records.*')), result.stderr

    cancel_log = root / 'cancel.log'
    with cancel_log.open('w') as err:
        process = subprocess.Popen(command('cancelled', progress='on'), stderr=err,
                                   stdout=subprocess.DEVNULL, env=dict(os.environ,
                                       RAMAG_TEST_PAUSE_STAGE='centromere-direction-plan'))
        try:
            deadline = time.monotonic() + 10
            while 'centromere-direction-plan' not in cancel_log.read_text():
                assert process.poll() is None and time.monotonic() < deadline
                time.sleep(0.02)
            process.send_signal(signal.SIGINT)
            assert process.wait(timeout=10) == 130
        finally:
            if process.poll() is None:
                process.kill()
                process.wait()
    assert not list(root.glob('cancelled.*'))

    inputs(ref, rc(query))
    result = run('reverse')
    assert result.returncode == 0, result.stderr
    reverse = validate(root / 'reverse', ref, rc(query))
    assert reverse and all(item[4] == '-' for item in reverse)
    inputs('N' * 80, 'N' * 93)
    result = run('empty')
    assert result.returncode == 0, result.stderr
    assert not validate(root / 'empty', 'N' * 80, 'N' * 93)

if SAMTOOLS:
    result = subprocess.run([SAMTOOLS, '--version'], capture_output=True, text=True, errors='replace')
    print(result.stdout.splitlines()[0] if result.stdout else SAMTOOLS)
else:
    print('SAMtools unavailable: optional external round-trip skipped; independent parsers passed')
print('stitch standard exports: coordinate/CIGAR/path round-trip, routing, preservation, and atomic failure passed')

"""Compile retained public APIs using only an isolated installed prefix."""
import argparse
import json
from pathlib import Path
import subprocess
import shlex
import tempfile
p = argparse.ArgumentParser()
p.add_argument('--build', type=Path, required=True)
p.add_argument('--source', type=Path, required=True)
p.add_argument('--root', type=Path, required=True)
p.add_argument('--cmake', required=True)
p.add_argument('--backend', required=True)
a = p.parse_args()
a.root.mkdir(parents=True, exist_ok=True)
with tempfile.TemporaryDirectory(prefix='consumer-', dir=a.root) as d:
    root = Path(d)
    prefix, build = root / 'prefix', root / 'build'
    commands = [
        [a.cmake, '--install', str(a.build), '--prefix', str(prefix)],
        [a.cmake, '-S', str(a.source / 'install'), '-B', str(build),
         '-DCMAKE_PREFIX_PATH=' + str(prefix), '-DEXPECT_BACKEND=' + a.backend],
        [a.cmake, '--build', str(build), '--parallel', '2'],
        [str(build / 'public-consumer')]]
    for command in commands:
        cp = subprocess.run(command, text=True, capture_output=True, timeout=90)
        assert cp.returncode == 0, (command, cp.stdout, cp.stderr)
    commands_json = json.loads((build / 'compile_commands.json').read_text())
    assert commands_json
    assert not (prefix / 'include' / 'ramag' / 'internal').exists()
    first = commands_json[0]
    compiler = shlex.split(first['command'])[0] if 'command' in first else first['arguments'][0]
    headers = sorted((prefix / 'include' / 'ramag').rglob('*.hpp'))
    expected_headers = {'alignment.hpp', 'model.hpp', 'pairwise_core.hpp',
                        'centromere.hpp', 'centromere/direction.hpp',
                        'centromere/anchors.hpp', 'centromere/stitch.hpp'}
    header_root = prefix / 'include' / 'ramag'
    assert {str(header.relative_to(header_root)) for header in headers} == expected_headers
    for header in headers:
        unit = root / ('standalone-' + header.stem + '.cpp')
        unit.write_text('#include <ramag/' + str(header.relative_to(header_root)) + '>\n')
        cp = subprocess.run([compiler, '-std=c++20', '-fsyntax-only',
                             '-I' + str(prefix / 'include'), str(unit)],
                            text=True, capture_output=True, timeout=30)
        assert cp.returncode == 0, (header.name, cp.stdout, cp.stderr)
    repository = a.source.resolve().parents[1]
    for entry in commands_json:
        command = entry.get('command', ' '.join(entry.get('arguments', [])))
        assert str(prefix / 'include') in command, command
        assert str(repository / 'include') not in command, command
print('Installed-only pairwise/centromere consumer passed')

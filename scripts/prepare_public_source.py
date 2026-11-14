#!/usr/bin/env python3
"""Export an explicit file list into a new reviewable product source directory.

This command never commits, merges or pushes. It never rewrites versions. The
manifest lists individual relative files, not globs; new files need explicit
classification. Internal tests and research artifacts cannot be exported.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import shutil
import tempfile

PRODUCT_DOCS = {
    'README.md', 'installation.md', 'quick-start.md', 'user-guide.md',
    'algorithms.md', 'library.md', 'developer-guide.md', 'migration.md',
    'mainpath-linking.md', 'centromere-interface.md', 'centromere-directions.md',
    'centromere-directed-anchors.md', 'centromere-copy-anchors.md',
    'centromere-anchor-stitch.md', 'centromere-rama-integration.md',
    'centromere-alignment-formats.md'}
ROOT_FILES = {'.gitignore', 'README.md', 'CHANGELOG.md', 'LICENSE',
              'THIRD_PARTY_NOTICES.md', 'CONTRIBUTING.md', 'CMakeLists.txt',
              'CMakePresets.json', 'public-source-files.txt'}
PRODUCT_CMAKE = {'RaMACentromere.cmake', 'RaMADependencies.cmake',
                 'RaMAGitSnapshot.cmake', 'RaMAGPairwiseCore.cmake',
                 'RaMAGCentromereConfig.cmake.in', 'RaMAGPairwiseConfig.cmake.in'}

def allowed(name):
    parts = PurePosixPath(name).parts
    if len(parts) == 1:
        return name in ROOT_FILES
    if parts[0] in ('src', 'include', 'third_party', 'examples'):
        return True
    if name.startswith('tests/public/') or name.startswith('.github/workflows/'):
        return True
    if parts[0] == 'docs':
        return len(parts) == 2 and parts[1] in PRODUCT_DOCS
    if parts[0] == 'cmake':
        return len(parts) == 2 and parts[1] in PRODUCT_CMAKE
    return name in ('scripts/plot_centromere_directions.py', 'scripts/prepare_public_source.py')

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--manifest', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    source = args.source.resolve(strict=True)
    manifest = args.manifest or source / 'public-source-files.txt'
    output = args.output.resolve()
    evidence = output.with_name(output.name + '.manifest.json')
    if output.exists() or evidence.exists() or source == output or source in output.parents:
        parser.error('output and evidence must be new paths outside the source tree')
    names = [line.strip() for line in manifest.read_text().splitlines()
             if line.strip() and not line.lstrip().startswith('#')]
    if not names or len(names) != len(set(names)):
        parser.error('manifest must list unique individual files')
    checked = []
    for name in names:
        rel = PurePosixPath(name)
        if rel.is_absolute() or '..' in rel.parts or str(rel) != name or not allowed(name):
            parser.error('unclassified or forbidden manifest path: ' + name)
        original = source.joinpath(*rel.parts)
        if original.is_symlink() or not original.is_file() or source not in original.resolve().parents:
            parser.error('missing, symlink or outside-source manifest file: ' + name)
        checked.append((name, original))
    output.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix='.' + output.name + '-staging-', dir=output.parent))
    records = []
    for name, original in checked:
        target = staging / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(original, target)
        digest = hashlib.sha256(target.read_bytes()).hexdigest()
        if digest != hashlib.sha256(original.read_bytes()).hexdigest():
            raise RuntimeError('source changed during export: ' + name)
        records.append({'path': name, 'size': target.stat().st_size, 'sha256': digest})
    os.rename(staging, output)
    evidence.write_text(json.dumps({'format': 1, 'files': records,
                                    'version_rewritten': False}, indent=2) + '\n')
    print('Exported', len(records), 'explicit files to', output)

if __name__ == '__main__':
    main()

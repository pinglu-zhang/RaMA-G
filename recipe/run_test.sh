#!/usr/bin/env bash
set -euo pipefail
ramag --version | grep -F "ramag ${PKG_VERSION} "
ramag align --help >/dev/null
python rama-g/tests/public/test_cli.py --binary "${PREFIX}/bin/ramag" --root "${PWD}/cli-check"
python rama-g/tests/public/test_formats.py "${PREFIX}/bin/ramag"
python rama-g/tests/public/test_centromere.py "${PREFIX}/bin/ramag"
report="$(ldd "${PREFIX}/bin/ramag")"
printf '%s\n' "${report}"
if printf '%s\n' "${report}" | grep -qE 'not found|_build_env|_h_env'; then
  echo 'Unresolved or build-private runtime dependency' >&2
  exit 1
fi
test -s "${PREFIX}/share/ramag/LICENSE"
test -s "${PREFIX}/share/ramag/third_party/sufkit/LICENSE"
test -s "${PREFIX}/share/ramag/third_party/rama/LICENSE"
test -s "${PREFIX}/share/ramag/docs/user-guide.md"

python - <<'PY'
import os, re, subprocess
from pathlib import Path
binary = Path(os.environ['PREFIX']) / 'bin/ramag'
info = subprocess.check_output(['readelf', '--version-info', binary], text=True)
versions = {tuple(map(int, v.split('.'))) for v in re.findall(r'\bGLIBC_(\d+\.\d+(?:\.\d+)?)\b', info)}
assert versions and max(versions) <= (2, 17), sorted(versions)
print('Executable maximum required glibc:', '.'.join(map(str, max(versions))))
PY

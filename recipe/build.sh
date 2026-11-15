#!/usr/bin/env bash
set -euo pipefail

export CFLAGS="${CFLAGS:-} -march=x86-64 -mtune=generic"
export CXXFLAGS="${CXXFLAGS:-} -march=x86-64 -mtune=generic"
source_dir="${SRC_DIR}/rama-g"
cmake -S "${source_dir}" -B build-conda -G Ninja \
  ${CMAKE_ARGS:-} \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
  -DCMAKE_INSTALL_PREFIX="${PREFIX}" \
  -DCMAKE_INSTALL_LIBDIR=lib \
  -DCMAKE_PREFIX_PATH="${PREFIX}" \
  -DCMAKE_INSTALL_RPATH="${PREFIX}/lib" \
  -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=OFF \
  -DBUILD_SHARED_LIBS=OFF \
  -DBUILD_TESTING=ON -DRAMAG_BUILD_TESTS=ON \
  -DRAMAG_ENABLE_CENTROMERE=ON \
  -DRAMAG_SUFKIT_SOURCE_DIR="${SRC_DIR}/sufkit" \
  -DRAMAG_WFA2_SOURCE_DIR="${SRC_DIR}/wfa2"
cmake --build build-conda --parallel "${CPU_COUNT:-4}"
python - <<'PY'
import json, shlex
from pathlib import Path
units = json.loads(Path('build-conda/compile_commands.json').read_text())
assert units, 'Missing compile commands'
for unit in units:
    flags = shlex.split(unit['command']) if 'command' in unit else unit['arguments']
    arches = [flag for flag in flags if flag.startswith('-march=')]
    assert arches and arches[-1] == '-march=x86-64', unit['file']
    assert not any(flag.startswith(('-mavx', '-mfma')) for flag in flags), unit['file']
print('Validated portable CPU flags for', len(units), 'translation units')
PY
ctest --test-dir build-conda --output-on-failure --no-tests=error -j 1
cmake --install build-conda

install -d "${PREFIX}/share/ramag/docs" "${PREFIX}/share/ramag/examples"
install -m 0644 "${source_dir}/README.md" "${PREFIX}/share/ramag/README.md"
cp -R "${source_dir}/docs/." "${PREFIX}/share/ramag/docs/"
cp -R "${source_dir}/examples/." "${PREFIX}/share/ramag/examples/"
"${PREFIX}/bin/ramag" --version | grep -F "ramag ${PKG_VERSION} "

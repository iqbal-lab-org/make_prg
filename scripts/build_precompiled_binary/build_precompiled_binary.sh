#!/usr/bin/env bash
set -eu

version="0.4.0"

CURRENT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
SCRIPTS_DIR="$(dirname "${CURRENT_DIR}")"
MAKE_PRG_DIR="$(dirname "${SCRIPTS_DIR}")"
PORTABLE_EXECUTABLE_BUILD_DIR="${MAKE_PRG_DIR}/precompiled_binary"

cd "$MAKE_PRG_DIR"

if [ -d "${PORTABLE_EXECUTABLE_BUILD_DIR}" ]; then
  echo "Please remove ${PORTABLE_EXECUTABLE_BUILD_DIR} before proceeding."
  exit 1
fi

docker run --rm -v "$(pwd)":/make_prg make_prg:${version} /bin/bash -c \
  "pyinstaller \
  --noconfirm \
  --clean \
  --onefile \
  --distpath precompiled_binary/dist \
  --workpath precompiled_binary/build \
  --specpath precompiled_binary/spec \
  --name make_prg_${version} \
  --hidden-import=\"sklearn.utils._weight_vector\" \
  make_prg/__main__.py"

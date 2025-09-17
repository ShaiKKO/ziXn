#!/usr/bin/env bash
##!
# \file wsl-determinism.sh
# \brief Run two deterministic scene executions and diff telemetry for parity.
# \author Colin Macritchie (Ripple Group, LLC)
# \version 1.0.0
# \date 2025-09-16
##

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
BUILD_DIR="${ROOT_DIR}/wsl-build"
# Default multi-config (Windows/MSVC) uses configuration subfolder; single-config (Makefiles/Ninja) does not.
CONFIG="${CONFIG:-Release}"
# Detect candidate binary locations cross-platform
BIN_CANDIDATES=(
  "${BUILD_DIR}/zx_cli"
  "${BUILD_DIR}/${CONFIG}/zx_cli.exe"
  "${BUILD_DIR}/zx_cli.exe"
)
BIN=""
for c in "${BIN_CANDIDATES[@]}"; do
  if [ -f "$c" ]; then BIN="$c"; break; fi
done
OUT_A="${BUILD_DIR}/telem_a.json"
OUT_B="${BUILD_DIR}/telem_b.json"
SCENE="dambreak"
SEED="42"
FALLBACK="off"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --scene) SCENE="$2"; shift 2 ;;
    --seed) SEED="$2"; shift 2 ;;
    --fallback) FALLBACK="$2"; shift 2 ;;
    *) echo "Unknown arg: $1" >&2; exit 2 ;;
  esac
done

# Ensure binary exists; if not, configure and build, passing --config for multi-config generators
if [[ -z "${BIN}" ]]; then
  echo "zx_cli not found, configuring & building fresh..." >&2
  rm -rf "${BUILD_DIR}"
  cmake -S "${ROOT_DIR}" -B "${BUILD_DIR}" -DZX_ENABLE_TESTS=ON -DZX_BUILD_SHARED=ON
  cmake --build "${BUILD_DIR}" --config "${CONFIG}" -j
  # Re-resolve after build
  for c in "${BIN_CANDIDATES[@]}"; do
    if [ -f "$c" ]; then BIN="$c"; break; fi
  done
fi

if [[ -z "${BIN}" ]]; then
  echo "ERROR: zx_cli not found after build in ${BUILD_DIR} (checked ${BIN_CANDIDATES[*]})" >&2
  exit 1
fi

"${BIN}" --mode scene --scene "${SCENE}" --deterministic 1 --seed "${SEED}" --fallback "${FALLBACK}" --telemetry-out "${OUT_A}"
"${BIN}" --mode scene --scene "${SCENE}" --deterministic 1 --seed "${SEED}" --fallback "${FALLBACK}" --telemetry-out "${OUT_B}"

# Sanitize time-variant fields (e.g., *_ms) before diff to avoid false negatives
SAN_A="${BUILD_DIR}/telem_a.sanitized.json"
SAN_B="${BUILD_DIR}/telem_b.sanitized.json"
sed -E 's/("[A-Za-z_]*_ms"\s*:\s*)[-0-9eE\.]+/\10/g' "${OUT_A}" > "${SAN_A}"
sed -E 's/("[A-Za-z_]*_ms"\s*:\s*)[-0-9eE\.]+/\10/g' "${OUT_B}" > "${SAN_B}"

hash_cmd=$(command -v sha256sum || command -v sha1sum || true)
size_a=$(stat -c%s "${SAN_A}" 2>/dev/null || wc -c < "${SAN_A}")
size_b=$(stat -c%s "${SAN_B}" 2>/dev/null || wc -c < "${SAN_B}")
hash_a=$([ -n "$hash_cmd" ] && $hash_cmd "${SAN_A}" | awk '{print $1}' || echo "nohash")
hash_b=$([ -n "$hash_cmd" ] && $hash_cmd "${SAN_B}" | awk '{print $1}' || echo "nohash")

if diff -q "${SAN_A}" "${SAN_B}" >/dev/null 2>&1; then
  echo "Determinism OK: telemetry identical | sizeA=${size_a} sizeB=${size_b} hashA=${hash_a} hashB=${hash_b}"
else
  echo "Determinism MISMATCH: printing unified diff" >&2
  diff -u "${OUT_A}" "${OUT_B}" || true
  exit 1
fi



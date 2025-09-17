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
# Prefer Windows .exe on Windows runners; else ELF
if [[ "${OS:-}" == "Windows_NT" ]]; then
  BIN="${BUILD_DIR}/Release/zx_cli.exe"
else
  BIN="${BUILD_DIR}/zx_cli"
fi
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

if [[ ! -x "${BIN}" ]]; then
  echo "zx_cli not found, configuring & building fresh..." >&2
  rm -rf "${BUILD_DIR}"
  cmake -S "${ROOT_DIR}" -B "${BUILD_DIR}" -DZX_ENABLE_TESTS=ON -DZX_BUILD_SHARED=ON
  if [[ "${OS:-}" == "Windows_NT" ]]; then
    cmake --build "${BUILD_DIR}" --config Release -j
    BIN="${BUILD_DIR}/Release/zx_cli.exe"
  else
    cmake --build "${BUILD_DIR}" -j
    BIN="${BUILD_DIR}/zx_cli"
  fi
fi

if [[ ! -x "${BIN}" ]]; then
  echo "ERROR: zx_cli not found at ${BIN}" >&2
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



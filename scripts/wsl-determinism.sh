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
BIN="${BUILD_DIR}/zx_cli"
OUT_A="${BUILD_DIR}/telem_a.json"
OUT_B="${BUILD_DIR}/telem_b.json"

if [[ ! -x "${BIN}" ]]; then
  echo "zx_cli not found, building..." >&2
  cmake --build "${BUILD_DIR}" -j
fi

"${BIN}" --mode scene --scene dambreak --deterministic 1 --seed 42 --fallback off --telemetry-out "${OUT_A}"
"${BIN}" --mode scene --scene dambreak --deterministic 1 --seed 42 --fallback off --telemetry-out "${OUT_B}"

if diff -q "${OUT_A}" "${OUT_B}" >/dev/null 2>&1; then
  echo "Determinism OK: telemetry files identical"
else
  echo "Determinism MISMATCH: printing unified diff" >&2
  diff -u "${OUT_A}" "${OUT_B}" || true
  exit 1
fi



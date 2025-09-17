#!/usr/bin/env bash
##!
# \file wsl-run-tests.sh
# \brief CI-like test runner for ziXn under WSL (configures, builds, runs ctest, summarizes).
# \author Colin Macritchie (Ripple Group, LLC)
# \version 1.0.0
# \date 2025-09-16
# \license Proprietary — Copyright (c) 2025 Colin Macritchie / Ripple Group, LLC.
##

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
BUILD_DIR="${ROOT_DIR}/wsl-build"
BUILD_TYPE="${CMAKE_BUILD_TYPE:-Release}"

mkdir -p "${BUILD_DIR}"

if [[ ! -f "${BUILD_DIR}/CMakeCache.txt" ]] || [[ "${1:-}" == "--reconfigure" ]]; then
  cmake -S "${ROOT_DIR}" -B "${BUILD_DIR}" \
    -DCMAKE_BUILD_TYPE="${BUILD_TYPE}" \
    -DCMAKE_C_COMPILER=/usr/bin/gcc \
    -DCMAKE_CXX_COMPILER=/usr/bin/g++
fi

cmake --build "${BUILD_DIR}" -j

ctest --test-dir "${BUILD_DIR}" -j --output-on-failure | tee "${BUILD_DIR}/ctest.log"

echo
passed=$(grep -Eo "[0-9]+ tests passed" "${BUILD_DIR}/ctest.log" | tail -n1 || true)
failed=$(grep -Eo "[0-9]+ tests failed" "${BUILD_DIR}/ctest.log" | tail -n1 || true)
echo "Summary: ${passed:-unknown}, ${failed:-0 tests failed}"



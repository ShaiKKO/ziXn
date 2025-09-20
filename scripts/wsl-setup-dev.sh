#!/usr/bin/env bash

# WSL development bootstrap for ziXn
# - Installs toolchain (cmake, ninja, clang-format, clang-tidy)
# - Configures Ninja Release build with compile_commands.json
# - Builds and runs tests
# - Runs clang-format dry-run check and clang-tidy over core/src

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
BUILD_DIR="${ROOT_DIR}/wsl-build"

export DEBIAN_FRONTEND=noninteractive

echo "[1/5] Installing toolchain (requires sudo)"
sudo apt-get update -qq
sudo apt-get install -y -qq build-essential pkg-config curl git jq unzip wget cmake ninja-build python3 python3-pip clang-format || true

# Prefer distro clang-tidy; fall back to versioned package and symlink
if ! command -v clang-tidy >/dev/null 2>&1; then
  sudo apt-get install -y -qq clang-tidy || true
fi
if ! command -v clang-tidy >/dev/null 2>&1; then
  sudo apt-get install -y -qq clang-tidy-15 || true
  if [ -x /usr/bin/clang-tidy-15 ] && [ ! -x /usr/bin/clang-tidy ]; then
    sudo ln -sf /usr/bin/clang-tidy-15 /usr/bin/clang-tidy
  fi
fi

echo "Tool versions:"
(cmake --version | head -n1) || true
(ninja --version) || true
(g++ --version | head -n1) || true
(clang-format --version) || true
(clang-tidy --version) || true

echo "[2/5] Configuring Ninja Release build at ${BUILD_DIR}"
rm -rf "${BUILD_DIR}"
cmake -S "${ROOT_DIR}" -B "${BUILD_DIR}" \
  -DZX_ENABLE_TESTS=ON -DZX_BUILD_SHARED=ON \
  -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_EXPORT_COMPILE_COMMANDS=ON

echo "[3/5] Building"
cmake --build "${BUILD_DIR}" --config Release --parallel

echo "[4/5] Running tests"
ctest --test-dir "${BUILD_DIR}" --output-on-failure -C Release -E "determinism_script"

echo "[5/5] Running format and tidy checks"
# clang-format dry-run (matches CI)
if command -v clang-format >/dev/null 2>&1; then
  git -C "${ROOT_DIR}" ls-files "core/src/*.{c,cpp}" "core/include/**/*.{h,hpp}" "kernels/*.hlsl" \
    | tr '\n' '\0' | xargs -0 -I{} clang-format --dry-run -Werror -style=file {}
fi

# clang-tidy over core/src (same scope as CI)
if command -v clang-tidy >/dev/null 2>&1; then
  files=$(git -C "${ROOT_DIR}" ls-files 'core/src/*.cpp')
  clang-tidy -p "${BUILD_DIR}" ${files} | tee "${BUILD_DIR}/clang-tidy.log" || true
  echo "clang-tidy errors: $(grep -n '\berror:\b' "${BUILD_DIR}/clang-tidy.log" | wc -l || true)"
  echo "clang-tidy warnings: $(grep -n '\bwarning:\b' "${BUILD_DIR}/clang-tidy.log" | wc -l || true)"
fi

echo "WSL setup complete."



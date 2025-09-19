#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
MODE="nd"
COUNT=60
SEED=123
SCENE="dambreak"
if [[ $# -ge 1 ]]; then MODE="$1"; fi
if [[ "$MODE" == "nd" ]]; then DET=0; else DET=1; fi
for i in $(seq 1 ${COUNT}); do
  ./wsl-build/zx_cli --mode scene --scene "$SCENE" --deterministic ${DET} --seed ${SEED} >/dev/null
done

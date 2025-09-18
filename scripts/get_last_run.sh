#!/usr/bin/env bash
set -euo pipefail

repo_root="/mnt/c/Users/cisco/Desktop/dev/research/ziXn"
cd "$repo_root"
mkdir -p docs/logs

# Fetch latest run metadata for our workflow
gh run list --workflow "ci.yml" -L 1 --json databaseId,conclusion,status,displayTitle > docs/logs/last_run.json

# Extract run id (handle possible UTF-8 BOM)
python3 - <<'PY' > docs/logs/last_run_id.txt
import json, io
with io.open('docs/logs/last_run.json','r',encoding='utf-8-sig') as f:
    data = json.load(f)
rid = data[0]['databaseId']
print(rid)
PY

rid=$(tr -d '\r\n' < docs/logs/last_run_id.txt)

# Fetch full logs for that run
gh run view "$rid" --log > docs/logs/gha_ci_latest.log
echo "Saved docs/logs/gha_ci_latest.log (run $rid)" >&2



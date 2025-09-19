$ErrorActionPreference = 'Stop'
$repoRoot = Resolve-Path (Join-Path $PSScriptRoot '..')
Set-Location -LiteralPath $repoRoot

# Inputs
$BaseBranch   = 'main'
$SourceBranch = 'feature/p3-streaming-fallback'
$NewBranch    = 'feature/p3-streaming-fallback-v2'

Write-Host ("Repo: {0}" -f $repoRoot)

# 1) Close old PR (ignore failure if already closed)
try {
  gh pr close 11 --comment 'Superseded by v2 PR' --delete-branch=false | Out-Null
  Write-Host 'Closed PR #11'
} catch {
  Write-Host 'PR #11 close skipped'
}

# 2) Clean and fetch
git reset --hard HEAD | Out-Null
git clean -fdx -e .git -e .gitignore | Out-Null
git fetch origin --prune | Out-Null

# Delete mistaken local branch named "origin/main" if present
if (git branch --list "origin/main") {
  git branch -D "origin/main" | Out-Null
}

# 3) Create fresh branch from base
git switch -C $NewBranch origin/$BaseBranch | Out-Null

# 4) Squash-merge source branch content
git fetch origin $SourceBranch | Out-Null
git merge --squash origin/$SourceBranch | Out-Null
if ($LASTEXITCODE -ne 0) {
  Write-Host 'Merge conflicts during squash. Aborting.'
  exit 1
}

# 5) Commit and push
if (-not (git diff --cached --quiet)) {
  git commit -m "feat: Windows tidy clean + CI readiness (v2, squashed)" | Out-Null
} else {
  git commit --allow-empty -m "chore: v2 PR scaffolding (no-op)" | Out-Null
}
git push -u origin $NewBranch | Out-Null

# 6) Create PR
gh pr create -B $BaseBranch -H $NewBranch -t "Windows tidy clean + CI readiness (v2, squashed)" -b "Follow-up to PR #11; squashed for CI stability." | Out-Null

# 7) Trigger CI and show first diagnostics
gh workflow run ci --ref $NewBranch | Out-Null
Start-Sleep -Seconds 10
$fetch = Join-Path $PSScriptRoot 'fetch-workflow.ps1'
Write-Host ("Fetch script: {0}" -f $fetch)
& $fetch -Workflow ci -PullRequestOnly | Out-Host



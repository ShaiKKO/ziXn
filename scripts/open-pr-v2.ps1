$ErrorActionPreference = 'Stop'
$repoRoot = Resolve-Path (Join-Path $PSScriptRoot '..')
Set-Location -LiteralPath $repoRoot

$newBranch = 'feature/p3-streaming-fallback-v2'

Write-Host ("Repo: {0}" -f $repoRoot)
$cur = (git rev-parse --abbrev-ref HEAD).Trim()
Write-Host ("Current branch: {0}" -f $cur)

# Create new branch from current HEAD (no file changes)
if (-not (git branch --list $newBranch)) {
  git switch -c $newBranch | Out-Null
} else {
  git switch $newBranch | Out-Null
}

# Push branch
git push -u origin $newBranch | Out-Host

# Create PR targeting main
gh pr create -B main -H $newBranch -t "Windows tidy clean + CI readiness (v2)" -b "Follow-up to #11; continuing from current HEAD." | Out-Host

# Trigger CI for this branch and print initial log excerpt
gh workflow run ci --ref $newBranch | Out-Host
Start-Sleep -Seconds 10

$run = gh run list --workflow ci --branch $newBranch --limit 1 --json databaseId | ConvertFrom-Json
if ($run -and $run[0] -and $run[0].databaseId) {
  $rid = $run[0].databaseId
  $log = "ci-$rid.log"
  gh run view $rid --log | Out-File -FilePath $log -Encoding utf8
  Write-Host '=== First 120 lines ==='
  Get-Content $log -TotalCount 120 | ForEach-Object { $_ }
  Write-Host '=== clang-tidy mentions ==='
  (Select-String -Path $log -Pattern 'clang-tidy' | Select-Object -First 20) | ForEach-Object { "{0}:{1}" -f $_.LineNumber, $_.Line }
  Write-Host '=== : error: lines (first 60) ==='
  (Select-String -Path $log -Pattern ': error:' | Select-Object -First 60) | ForEach-Object { "{0}:{1}" -f $_.LineNumber, $_.Line }
} else {
  Write-Host 'No CI run yet.'
}



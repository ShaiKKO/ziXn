param(
  [string]$Branch = 'feature/p3-streaming-fallback-v2',
  [string]$Workflow = 'ci'
)
$ErrorActionPreference = 'Stop'

$repoRoot = Resolve-Path (Join-Path $PSScriptRoot '..')
Set-Location -LiteralPath $repoRoot
$env:GH_PAGER = 'cat'

$runs = gh run list --workflow $Workflow --branch $Branch --limit 1 --json databaseId,status,conclusion | ConvertFrom-Json
if (-not $runs -or -not $runs[0]) { Write-Host 'No runs found'; exit 0 }
$rid = $runs[0].databaseId
Write-Host ("RUN_ID={0} STATUS={1} CONC={2}" -f $rid, $runs[0].status, $runs[0].conclusion)

$log = "${Workflow}-${rid}.log"
gh run view $rid --log | Out-File -FilePath $log -Encoding utf8

Write-Host '=== First 120 lines ==='
Get-Content -Path $log -TotalCount 120 | ForEach-Object { $_ }

Write-Host '=== failure markers ==='
Select-String -Path $log -Pattern '##\[error\]|invalid path|fatal:|The process .* failed with exit code|CMake Error|MSB\d{4}' | Select-Object -First 80 | ForEach-Object { "{0}:{1}" -f $_.LineNumber, $_.Line }



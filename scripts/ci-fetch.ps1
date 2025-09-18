param(
    [string]$Workflow = 'ci',
    [string]$Branch = ''
)

$ErrorActionPreference = 'Stop'
$repoRoot = Resolve-Path (Join-Path $PSScriptRoot '..')
Set-Location -LiteralPath $repoRoot
$env:GH_PAGER = 'cat'

function Get-HeadBranch {
    param([int]$PrNumber)
    try {
        $pr = gh pr view $PrNumber --json headRefName | ConvertFrom-Json
        if ($pr -and $pr.headRefName) { return $pr.headRefName }
    } catch {}
    try {
        $b = (git rev-parse --abbrev-ref HEAD).Trim()
        if ($b) { return $b }
    } catch {}
    return ''
}

if (-not $Branch -or $Branch -eq '') {
    $Branch = Get-HeadBranch -PrNumber 11
}

$runs = @()
try {
    if ($Branch) {
        $runs = gh run list --workflow $Workflow --branch $Branch --limit 5 --json databaseId,headBranch,status,conclusion,createdAt | ConvertFrom-Json
    }
} catch {}
if (-not $runs -or $runs.Count -eq 0) {
    try {
        $runs = gh run list --workflow $Workflow --limit 5 --json databaseId,headBranch,status,conclusion,createdAt | ConvertFrom-Json
    } catch {}
}
if (-not $runs -or $runs.Count -eq 0) {
    Write-Host "No $Workflow runs found."
    exit 0
}

$run = ($runs | Sort-Object -Property createdAt -Descending | Select-Object -First 1)
$rid = $run.databaseId
Write-Host ("RUN_ID={0} BRANCH={1} STATUS={2} CONC={3}" -f $rid, $run.headBranch, $run.status, $run.conclusion)

$log = "ci-$rid.log"
gh run view $rid --log | Out-File -FilePath $log -Encoding utf8

Write-Host '=== First 120 lines ==='
Get-Content -Path $log -TotalCount 120 | ForEach-Object { $_ }

Write-Host '=== clang-tidy mentions ==='
$ct = Select-String -Path $log -Pattern 'clang-tidy' | Select-Object -First 20
$ct | ForEach-Object { "{0}:{1}" -f $_.LineNumber, $_.Line }

Write-Host '=== : error: lines (first 100) ==='
$errs = Select-String -Path $log -Pattern ': error:' | Select-Object -First 100
$errs | ForEach-Object { "{0}:{1}" -f $_.LineNumber, $_.Line }



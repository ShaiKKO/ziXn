param(
  [string]$Branch = 'feature/p3-streaming-fallback-v2'
)

$ErrorActionPreference = 'Stop'

$repoRoot = Resolve-Path (Join-Path $PSScriptRoot '..')
Set-Location -LiteralPath $repoRoot

git switch $Branch | Out-Null

# Ensure we are aligned with remote tip in the index without touching the working tree
git fetch origin $Branch | Out-Null
git reset --mixed ("origin/" + $Branch) | Out-Null

# Detect the stray invalid filename by searching the remote tree and index using NUL-safe lists
$rawRemote = git ls-tree -r -z --name-only --full-tree ("origin/" + $Branch)
$remoteList = @()
if ($rawRemote) { $remoteList = $rawRemote -split "`0" }

$rawIndex = git ls-files -z
$indexList = @()
if ($rawIndex) { $indexList = $rawIndex -split "`0" }

$bad = @()
if ($remoteList) { $bad += ($remoteList | Where-Object { $_ -like 'hell -NoProfile*' -or $_ -match 'ExecutionPolicy Bypass|ErrorActionPreference' }) }
if ($indexList) { $bad += ($indexList | Where-Object { $_ -like 'hell -NoProfile*' -or $_ -match 'ExecutionPolicy Bypass|ErrorActionPreference' }) }
$bad = $bad | Where-Object { $_ -and $_.Length -gt 0 } | Sort-Object -Unique

if ($bad -and $bad.Count -gt 0) {
  Write-Host "Found suspicious tracked paths:" -ForegroundColor Yellow
  $bad | ForEach-Object { Write-Host "  $_" }

  # Use a pathspec file to robustly remove from index even if not present in the working tree
  $tmp = Join-Path $env:TEMP "paths-to-remove.txt"
  # Write NUL-separated for exact matching of special characters
  $bytes = New-Object System.Collections.Generic.List[byte]
  foreach ($p in $bad) {
    $pBytes = [System.Text.Encoding]::UTF8.GetBytes($p)
    [void]$bytes.AddRange($pBytes)
    $bytes.Add(0)
  }
  [System.IO.File]::WriteAllBytes($tmp, $bytes.ToArray())

  # Remove from index only; working tree may not contain the file on Windows
  git rm -f -r --cached --pathspec-from-file=$tmp --pathspec-file-nul | Out-Null
  Remove-Item $tmp -Force -ErrorAction SilentlyContinue

  git commit -m "ci: remove invalid path that breaks Windows checkout" | Out-Null
  git push -u origin $Branch | Out-Host
} else {
  Write-Host 'No invalid path found in remote or index.'
}



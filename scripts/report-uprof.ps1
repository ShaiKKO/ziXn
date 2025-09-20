<#
/**
 * @file report-uprof.ps1
 * @brief Generate AMD uProf CSV report for an existing TBP session and print a summary.
 * @details
 *   Processes a previously collected AMD uProf session directory and emits
 *   a CSV report (`report.csv`) using AMDuProfCLI. Prints the first lines for
 *   quick inspection and attempts to extract the HOT FUNCTIONS section.
 *   Deterministic, non-interactive, and suitable for CI.
 *
 *   Owner: Colin Macritchie / Ripple Group, LLC. Engineering Standards (2025)
 *   apply: clear errors, deterministic outputs, documented behavior.
 *
 * @copyright
 *   (c) 2025 Colin Macritchie / Ripple Group, LLC. All rights reserved.
 *   Licensed for use within the ziXn project under project terms.
 */
#>

[CmdletBinding(PositionalBinding=$false)]
param(
  [Parameter(Mandatory=$true)]
  [string] $Root,
  [string] $OutputCsv = 'report.csv',
  [int]    $Head = 80,
  [string] $ExtraSymbolPath = ''
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

function Write-Info([string] $msg) { Write-Host "[report-uprof] $msg" -ForegroundColor Cyan }
function Fail([string] $msg) { Write-Error "[report-uprof] $msg"; exit 1 }

function Resolve-AMDuProfCLI {
  $path = (Get-Command AMDuProfCLI -ErrorAction SilentlyContinue | Select-Object -First 1).Path
  if (-not $path) {
    $path = (Get-ChildItem 'C:\Program Files','C:\Program Files (x86)' -Recurse -Filter 'AMDuProfCLI.exe' -ErrorAction SilentlyContinue |
      Select-Object -First 1 -ExpandProperty FullName)
  }
  if (-not $path) { Fail 'AMDuProfCLI.exe not found. Install AMD uProf and re-run.' }
  return $path
}

function Find-SessionRoot([string] $root) {
  # Use provided root if it contains AMDuProf session markers; else, pick latest subfolder.
  $markers = @(Get-ChildItem -LiteralPath $root -Recurse -Directory -ErrorAction SilentlyContinue |
    Where-Object { $_.Name -like 'AMDuProf-*' })
  if ($markers -and $markers.Count -gt 0) {
    return ($markers | Sort-Object LastWriteTime -Descending | Select-Object -First 1 -ExpandProperty FullName)
  }
  return $root
}

try {
  if (-not (Test-Path -LiteralPath $Root)) { Fail ("Root not found: " + $Root) }
  $baseSym = 'srv;C:\Symbols;srv*https://msdl.microsoft.com/download/symbols'
  $repoRoot = Split-Path -Parent $PSScriptRoot
  $relBuild = Join-Path $repoRoot 'build-msvc-rel\Release'
  if (-not $ExtraSymbolPath -and (Test-Path -LiteralPath $relBuild)) { $ExtraSymbolPath = $relBuild }
  if ($ExtraSymbolPath) { $env:_NT_SYMBOL_PATH = $baseSym + ';' + $ExtraSymbolPath } else { $env:_NT_SYMBOL_PATH = $baseSym }
  $uProf = Resolve-AMDuProfCLI
  $session = Find-SessionRoot -root $Root
  $csv = Join-Path $Root $OutputCsv

  Write-Info ("uProf: " + $uProf)
  Write-Info ("session: " + $session)
  Write-Info ("csv: " + $csv)

  # AMDuProfCLI: emit CSV in-place (no -o). Then discover most recent CSV.
  & $uProf report -i $session -f csv | Out-Null
  $csvCand = @(Get-ChildItem -Path $session -Recurse -Filter *.csv -File -ErrorAction SilentlyContinue |
    Sort-Object LastWriteTime -Descending | Select-Object -First 1 -ExpandProperty FullName)
  if (-not $csvCand -or $csvCand.Count -eq 0) {
    $csvCand = @(Get-ChildItem -Path $Root -Recurse -Filter *.csv -File -ErrorAction SilentlyContinue |
      Sort-Object LastWriteTime -Descending | Select-Object -First 1 -ExpandProperty FullName)
  }
  if ($csvCand -and $csvCand.Count -gt 0) { $csv = $csvCand[0] }
  if (-not (Test-Path -LiteralPath $csv)) {
    Fail 'Report not generated: no CSV found.'
  }

  Write-Info 'HEAD:'
  Get-Content -LiteralPath $csv -Head $Head | Write-Output

  # Attempt to extract top functions section
  Write-Info 'HOT FUNCTIONS (CPU_TIME):'
  $lines = @(Get-Content -LiteralPath $csv)
  $hotA = ($lines | Select-String -SimpleMatch 'HOTTEST FUNCTIONS (Sort Event - CPU_TIME)')
  if (-not $hotA) { $hotA = ($lines | Select-String -SimpleMatch 'HOT FUNCTIONS (Sort Event - CPU_TIME)') }
  if ($hotA) {
    $start = $hotA[0].LineNumber
    $idx = [int]$start - 1
    $end = [Math]::Min($idx + 40, $lines.Count - 1)
    $lines[$idx..$end] | Write-Output
  } else {
    Write-Output 'No explicit HOT FUNCTIONS section found in CSV.'
    # For memory sessions, emit a short cacheline summary if present
    $memHdr = ($lines | Select-String -SimpleMatch 'SHARED DATA CACHELINE SUMMARY')
    if ($memHdr) {
      $mstart = $memHdr[0].LineNumber
      $mid = [int]$mstart
      $mend = [Math]::Min($mid + 20, $lines.Count - 1)
      Write-Info 'CACHELINE SUMMARY:'
      $lines[$mid..$mend] | Write-Output
    }
  }

} catch {
  Fail ("Unhandled error: " + $_.Exception.Message)
}



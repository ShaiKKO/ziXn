<#
/**
 * @file profile-uprof.ps1
 * @brief Collect AMD uProf TBP with call stacks for zx_cli and emit CSV report.
 * @details
 *   Windows-first profiling helper for ziXn. Collects TBP (Time-Based Profiling)
 *   with FPO call stacks for release builds of `zx_cli`, writes outputs under
 *   `traces/`, and produces a `report.csv` via `AMDuProfCLI report`. Designed for
 *   CI and local usage: deterministic filenames, explicit error handling, no
 *   interactive prompts. Honors _NT_SYMBOL_PATH for symbolization.
 *
 *   Owner: Colin Macritchie / Ripple Group, LLC. Engineering Standards (2025)
 *   apply: no placeholders, clear errors, deterministic outputs, and docs.
 *
 * @copyright
 *   (c) 2025 Colin Macritchie / Ripple Group, LLC. All rights reserved.
 *   Licensed for use within the ziXn project under project terms.
 */
#>

[CmdletBinding(PositionalBinding=$false)]
param(
  [string] $AppPath,
  [string] $AppArgs = "--mode bench --bench-up2 --size 2048 --iters 100000",
  [int]    $DurationSec = 90,
  [int]    $CallGraphDepth = 128,
  
  [int]    $CallGraphIntervalMs = 1,
  [string] $OutputDir = "",
  [switch] $VerboseLog
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

function Write-Info([string] $msg) { Write-Host "[profile-uprof] $msg" -ForegroundColor Cyan }
function Fail([string] $msg) { Write-Error "[profile-uprof] $msg"; exit 1 }

function Resolve-AMDuProfCLI {
  $path = (Get-Command AMDuProfCLI -ErrorAction SilentlyContinue | Select-Object -First 1).Path
  if (-not $path) {
    $path = (Get-ChildItem 'C:\Program Files','C:\Program Files (x86)' -Recurse -Filter 'AMDuProfCLI.exe' -ErrorAction SilentlyContinue |
      Select-Object -First 1 -ExpandProperty FullName)
  }
  if (-not $path) { Fail 'AMDuProfCLI.exe not found. Install AMD uProf and re-run.' }
  return $path
}

function Resolve-ZxCli([string] $hint) {
  if ($hint) {
    $p = Resolve-Path -LiteralPath $hint -ErrorAction SilentlyContinue
    if ($p) { return $p.Path }
  }
  $exe = Get-ChildItem -Path . -Recurse -File -Include 'zx_cli.exe' -ErrorAction SilentlyContinue |
    Sort-Object LastWriteTime -Descending | Select-Object -First 1 -ExpandProperty FullName
  if (-not $exe) { Fail 'zx_cli.exe not found. Build the project first.' }
  return $exe
}

function Ensure-OutputDir([string] $dir) {
  if (-not $dir) {
    $ts = Get-Date -Format 'yyyyMMdd_HHmmss'
    $dir = Join-Path (Resolve-Path '.').Path (Join-Path 'traces' ("uprof_tbp_" + $ts))
  }
  New-Item -ItemType Directory -Path $dir -Force | Out-Null
  return $dir
}

function Generate-ReportCsv([string] $uProf, [string] $root) {
  # Some AMDuProfCLI versions do not accept -o/--output-dir. Emit in-place and discover CSV.
  & $uProf report -i $root -f csv | Out-Null
  $csv = (Get-ChildItem -Path $root -Recurse -Filter *.csv -File -ErrorAction SilentlyContinue |
    Sort-Object LastWriteTime -Descending | Select-Object -First 1 -ExpandProperty FullName)
  if (-not $csv) { Fail ("CSV report not generated under: " + $root) }
  return $csv
}

try {
  if (-not $env:_NT_SYMBOL_PATH) {
    $env:_NT_SYMBOL_PATH = 'srv;C:\Symbols;srv*https://msdl.microsoft.com/download/symbols'
  }
  $uProf = Resolve-AMDuProfCLI
  $exe   = Resolve-ZxCli -hint $AppPath
  $out   = Ensure-OutputDir -dir $OutputDir

  Write-Info ("uProf: " + $uProf)
  Write-Info ("zx_cli: " + $exe)
  Write-Info ("output: " + $out)

  $args = @(
    'collect',
    '--config','tbp',
    '--call-graph-mode','fpo',
    '--call-graph-type','user',
    '--call-graph-interval', [string]$CallGraphIntervalMs,
    '--call-graph-depth',    [string]$CallGraphDepth,
    '-d', [string]$DurationSec,
    '-o', $out,
    $exe
  )

  if ($AppArgs) {
    $args += ($AppArgs -split '\s+')
  }

  if ($VerboseLog) {
    Write-Info ("AMDuProfCLI " + ($args -join ' '))
  }

  & $uProf @args | Write-Output
  Write-Info 'Collection completed.'

  # Find root output directory created by uProf (it nests a session folder)
  $root = (Get-ChildItem $out -Directory | Sort-Object LastWriteTime -Descending | Select-Object -First 1 -ExpandProperty FullName)
  if (-not $root) { $root = $out }

  $csv = Generate-ReportCsv -uProf $uProf -root $root
  Write-Info ("CSV: " + $csv)

  # Print first lines as a quick sanity check
  Get-Content -Path $csv -Head 60 | Write-Output

} catch {
  Fail ("Unhandled error: " + $_.Exception.Message)
}



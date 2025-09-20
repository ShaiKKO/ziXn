<#
/**
 * @file profile-uprof-memory.ps1
 * @brief Collect AMD uProf Memory/Cache profiling for zx_cli and emit CSV report.
 * @details
 *   Runs AMD uProf with the memory/config to capture cache behavior (incl. L3)
 *   for a sustained bench run of `zx_cli`. Outputs are written under `traces/`
 *   in a timestamped folder. Generates a CSV via `AMDuProfCLI report` and
 *   prints its head for quick inspection.
 *
 *   Owner: Colin Macritchie / Ripple Group, LLC. Engineering Standards (2025)
 *   apply: deterministic outputs, clear errors, Doxygen header, no prompts.
 */
#>

[CmdletBinding(PositionalBinding=$false)]
param(
  [string] $AppPath,
  [string] $AppArgs = "--mode bench --bench-up2 --size 2048 --iters 100000",
  [int]    $DurationSec = 60,
  [string] $Cpu = "",
  [string] $OutputDir = "",
  [switch] $VerboseLog
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

function Write-Info([string] $msg) { Write-Host "[profile-uprof-memory] $msg" -ForegroundColor Cyan }
function Fail([string] $msg) { Write-Error "[profile-uprof-memory] $msg"; exit 1 }

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
    $dir = Join-Path (Resolve-Path '.').Path (Join-Path 'traces' ("uprof_mem_" + $ts))
  }
  New-Item -ItemType Directory -Path $dir -Force | Out-Null
  return $dir
}

function Generate-ReportCsv([string] $uProf, [string] $root) {
  # Emit CSV in-place and locate it
  & $uProf report -i $root -f csv | Out-Null
  $csv = (Get-ChildItem -Path $root -Recurse -Filter *.csv -File -ErrorAction SilentlyContinue |
    Sort-Object LastWriteTime -Descending | Select-Object -First 1 -ExpandProperty FullName)
  if (-not $csv) { Fail ("CSV report not generated under: " + $root) }
  return $csv
}

try {
  $uProf = Resolve-AMDuProfCLI
  $exe   = Resolve-ZxCli -hint $AppPath
  $out   = Ensure-OutputDir -dir $OutputDir

  Write-Info ("uProf: " + $uProf)
  Write-Info ("zx_cli: " + $exe)
  Write-Info ("output: " + $out)

  $args = @(
    'collect',
    '--config','memory',
    $(if($Cpu){'--cpu'})
    $(if($Cpu){$Cpu})
    $(if($Cpu){'--affinity'})
    $(if($Cpu){$Cpu})
    '-d', [string]$DurationSec,
    '-o', $out,
    $exe
  )
  if ($AppArgs) { $args += ($AppArgs -split '\s+') }
  if ($VerboseLog) { Write-Info ("AMDuProfCLI " + ($args -join ' ')) }

  & $uProf @args | Write-Output
  if ($LASTEXITCODE -ne 0) { Fail ("uProf collect failed with exit code: " + $LASTEXITCODE) }
  Write-Info 'Memory collection completed.'

  $root = (Get-ChildItem $out -Directory | Sort-Object LastWriteTime -Descending | Select-Object -First 1 -ExpandProperty FullName)
  if (-not $root) { $root = $out }

  $csv = Generate-ReportCsv -uProf $uProf -root $root
  Write-Info ("CSV: " + $csv)
  Get-Content -Path $csv -Head 80 | Write-Output

} catch {
  Fail ("Unhandled error: " + $_.Exception.Message)
}



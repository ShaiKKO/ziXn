<#
/**
 * @file analyze-perfview.ps1
 * @brief Export fully symbolized top stacks from a PerfView ETL(.zip) headlessly.
 * @details
 *   Finds PerfView.exe, locates an ETL/ETL.ZIP trace under `traces/` (or a
 *   provided path), sets _NT_SYMBOL_PATH to include local PDBs and MS server,
 *   and invokes PerfView without GUI to save top stacks to a text file.
 *   Non-interactive; suitable for CI and local dev.
 *
 *   Owner: Colin Macritchie / Ripple Group, LLC. Engineering Standards (2025)
 *   apply: deterministic outputs, clear errors, Doxygen header.
 */
#>

[CmdletBinding(PositionalBinding=$false)]
param(
  [string] $InputTrace = "",
  [string] $OutputPath = "",
  [string] $ExtraSymbolPath = "",
  [string] $FocusProcess = "zx_cli"
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

function Write-Info([string] $msg) { Write-Host "[analyze-perfview] $msg" -ForegroundColor Cyan }
function Fail([string] $msg) { Write-Error "[analyze-perfview] $msg"; exit 1 }

function Resolve-PerfView {
  $pv = (Get-Command PerfView.exe -ErrorAction SilentlyContinue | Select-Object -First 1).Path
  if (-not $pv) {
    $pv = (Get-ChildItem 'C:\Users','C:\Program Files','C:\Program Files (x86)' -Recurse -Filter 'PerfView.exe' -ErrorAction SilentlyContinue |
      Select-Object -First 1 -ExpandProperty FullName)
  }
  if (-not $pv) { Fail 'PerfView.exe not found. Install PerfView and retry.' }
  return $pv
}

function Resolve-Trace([string] $hint) {
  if ($hint -and (Test-Path -LiteralPath $hint)) { return (Resolve-Path -LiteralPath $hint).Path }
  $traces = Join-Path (Resolve-Path '.').Path 'traces'
  $cand = @(Get-ChildItem -Path $traces -Recurse -File -Include '*.etl','*.etl.zip' -ErrorAction SilentlyContinue |
    Sort-Object LastWriteTime -Descending | Select-Object -First 1 -ExpandProperty FullName)
  if ($cand -and $cand.Count -gt 0) { return $cand[0] }
  Fail 'No ETL/ETL.ZIP trace found under traces/. Provide -InputTrace.'
}

try {
  $pv = Resolve-PerfView
  $trace = Resolve-Trace -hint $InputTrace
  if (-not $OutputPath) {
    $ts = Get-Date -Format 'yyyyMMdd_HHmmss'
    $outDir = Join-Path (Split-Path -Parent $trace) ('perfview_' + $ts)
    New-Item -ItemType Directory -Path $outDir -Force | Out-Null
    $OutputPath = Join-Path $outDir 'topstacks.txt'
  }
  $baseSym = 'srv;C:\Symbols;srv*https://msdl.microsoft.com/download/symbols'
  $repoRoot = Split-Path -Parent $PSScriptRoot
  $relBuild = Join-Path $repoRoot 'build-msvc-rel\Release'
  if (-not $ExtraSymbolPath -and (Test-Path -LiteralPath $relBuild)) { $ExtraSymbolPath = $relBuild }
  if ($ExtraSymbolPath) { $env:_NT_SYMBOL_PATH = $baseSym + ';' + $ExtraSymbolPath } else { $env:_NT_SYMBOL_PATH = $baseSym }

  Write-Info ("PerfView: " + $pv)
  Write-Info ("Trace: " + $trace)
  Write-Info ("Output: " + $OutputPath)

  # Try known CLI variants to save stacks headlessly
  $ok = $false
  $symArg = @()
  if ($env:_NT_SYMBOL_PATH) { $symArg = @('/SymbolPath:' + $env:_NT_SYMBOL_PATH) }
  $logPath = (Join-Path (Split-Path -Parent $OutputPath) 'perfview.log.txt')
  $cmds = @(
    @('/NoGui','/AcceptEula') + $symArg + @('/LogFile:' + $logPath, '/ThreadTime',                '/ProcessFilter:' + $FocusProcess, '/SaveThreadTimeStacks', $OutputPath, $trace),
    @('/NoGui','/AcceptEula') + $symArg + @('/LogFile:' + $logPath, '/ThreadTime',                '/OnlyProcess:'  + $FocusProcess, '/SaveThreadTimeStacks', $OutputPath, $trace),
    @('/NoGui','/AcceptEula') + $symArg + @('/LogFile:' + $logPath, '/StackViewer:ThreadTime',    '/ProcessFilter:' + $FocusProcess, '/SaveThreadTimeStacks', $OutputPath, $trace),
    @('/NoGui','/AcceptEula') + $symArg + @('/LogFile:' + $logPath, '/StackViewer:ThreadTime',    '/OnlyProcess:'  + $FocusProcess, '/SaveThreadTimeStacks', $OutputPath, $trace),
    @('/NoGui','/AcceptEula') + $symArg + @('/LogFile:' + $logPath,                                '/ProcessFilter:' + $FocusProcess, '/SaveThreadTimeStacks', $OutputPath, $trace),
    @('/NoGui','/AcceptEula') + $symArg + @('/LogFile:' + $logPath, '/StackViewer:CPUStacks',                                      '/SaveStacks',           $OutputPath, $trace)
  )
  foreach($args in $cmds) {
    try {
      & $pv @args | Out-Null
      if (Test-Path -LiteralPath $OutputPath) { $ok = $true; break }
    } catch {
      continue
    }
  }
  if (-not $ok) {
    if (Test-Path -LiteralPath $logPath) {
      Write-Info 'PerfView log (first 120 lines):'
      Get-Content -LiteralPath $logPath -Head 120 | Write-Output
    }
    Fail 'Failed to export top stacks via PerfView CLI.'
  }

  Write-Info 'TOP STACKS (head):'
  Get-Content -LiteralPath $OutputPath -Head 100 | Write-Output

} catch {
  Fail ("Unhandled error: " + $_.Exception.Message)
}



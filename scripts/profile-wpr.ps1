Param(
  [ValidateSet("cpu","heap","both")] [string]$Mode = 'cpu',
  [string]$OutDir = 'traces',
  [string]$Inst = 'zx_profile',
  [string]$TestCmd = 'build-rel\\zx_cli.exe --mode scene --scene dambreak --telemetry-out traces\\dambreak.json',
  [string]$WorkingDir = ''
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

function Resolve-WprPath {
  $wpr = Get-Command wpr.exe -ErrorAction SilentlyContinue
  if (-not $wpr) {
    $candidate = 'C:\Program Files (x86)\Windows Kits\10\Windows Performance Toolkit\wpr.exe'
    if (Test-Path $candidate) { return $candidate }
    throw 'WPR (wpr.exe) not found. Install Windows Performance Toolkit (ADK).'
  }
  return $wpr.Source
}

function Assert-Admin {
  $id  = [Security.Principal.WindowsIdentity]::GetCurrent()
  $pri = New-Object Security.Principal.WindowsPrincipal $id
  if (-not $pri.IsInRole([Security.Principal.WindowsBuiltinRole]::Administrator)) {
    throw 'Administrator privileges required. Re-run this script from an elevated PowerShell.'
  }
}

function Start-WprCpu([string]$wpr, [string]$instance) {
  $p = Start-Process -FilePath $wpr -ArgumentList ('-start','CPU','-filemode','-instancename', $instance) -PassThru -Wait -NoNewWindow
  if ($p.ExitCode -ne 0) {
    # Fallback to GeneralProfile if CPU profile fails
    $p2 = Start-Process -FilePath $wpr -ArgumentList ('-start','GeneralProfile','-filemode','-instancename', $instance) -PassThru -Wait -NoNewWindow
    if ($p2.ExitCode -ne 0) { throw "Failed to start WPR CPU/GeneralProfile (exit=$($p.ExitCode)/$($p2.ExitCode))" }
  }
}

function Start-WprHeap([string]$wpr, [string]$instance) {
  & $wpr -start Heap -heapuser -filemode -instancename $instance | Out-Null
}

function Stop-Wpr([string]$wpr, [string]$outFile) {
  $p = Start-Process -FilePath $wpr -ArgumentList ('-stop', $outFile) -PassThru -Wait -NoNewWindow
  if ($p.ExitCode -ne 0) { Write-Warning "WPR stop returned exit code $($p.ExitCode). Output may be missing." }
}

function Stop-WprRecording([string]$wpr) {
  try { & $wpr -cancel | Out-Null } catch { }
}

if ([string]::IsNullOrWhiteSpace($WorkingDir)) {
  # Default working dir to repo root (script's parent directory)
  $WorkingDir = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
}
Push-Location $WorkingDir
try { New-Item -ItemType Directory -Force -Path $OutDir | Out-Null } catch {}
$wprPath = Resolve-WprPath
Assert-Admin
Stop-WprRecording $wprPath

$cpuEtl = Join-Path $OutDir 'cpu.etl'
$heapEtl = Join-Path $OutDir 'heap.etl'

$startCpu = ($Mode -eq 'cpu' -or $Mode -eq 'both')
$startHeap = ($Mode -eq 'heap' -or $Mode -eq 'both')

try {
  if ($startCpu) { Start-WprCpu -wpr $wprPath -instance $Inst }
  if ($startHeap) { Start-WprHeap -wpr $wprPath -instance "$Inst-heap" }

  $psArgs = @('-NoProfile','-ExecutionPolicy','Bypass','-Command', $TestCmd)
  $p = Start-Process -FilePath 'powershell.exe' -ArgumentList $psArgs -WorkingDirectory $WorkingDir -PassThru
  Wait-Process -Id $p.Id

  if ($startCpu) { Stop-Wpr -wpr $wprPath -outFile $cpuEtl }
  if ($startHeap) { Stop-Wpr -wpr $wprPath -outFile $heapEtl }
}
catch {
  Stop-WprRecording $wprPath
  throw
}

Write-Host "Trace(s) written to:"
if ($startCpu) {
  if (Test-Path $cpuEtl) { Write-Host "  $cpuEtl" } else { Write-Warning "CPU ETL not found at $cpuEtl" }
}
if ($startHeap) {
  if (Test-Path $heapEtl) { Write-Host "  $heapEtl" } else { Write-Warning "Heap ETL not found at $heapEtl" }
}
Write-Host "Open ETL files with Windows Performance Analyzer (wpa.exe) for analysis."

Pop-Location


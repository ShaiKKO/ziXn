Param(
  [ValidateSet("cpu","heap","both")] [string]$Mode = 'cpu',
  [string]$OutDir = 'traces',
  [string]$Inst = 'zx_profile',
  [string]$TestCmd = 'ctest --test-dir build-rel -C Release -R dambreak_integration_test -V'
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
  & $wpr -start CPU -filemode -instancename $instance | Out-Null
}

function Start-WprHeap([string]$wpr, [string]$instance) {
  & $wpr -start Heap -heapuser -filemode -instancename $instance | Out-Null
}

function Stop-Wpr([string]$wpr, [string]$outFile) {
  & $wpr -stop $outFile | Out-Null
}

function Stop-WprRecording([string]$wpr) {
  try { & $wpr -cancel | Out-Null } catch { }
}

New-Item -ItemType Directory -Force -Path $OutDir | Out-Null
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
  $p = Start-Process -FilePath 'powershell.exe' -ArgumentList $psArgs -PassThru
  Wait-Process -Id $p.Id

  if ($startCpu) { Stop-Wpr -wpr $wprPath -outFile $cpuEtl }
  if ($startHeap) { Stop-Wpr -wpr $wprPath -outFile $heapEtl }
}
catch {
  Stop-WprRecording $wprPath
  throw
}

Write-Host "Trace(s) written to:"
if ($startCpu) { Write-Host "  $cpuEtl" }
if ($startHeap) { Write-Host "  $heapEtl" }
Write-Host "Open ETL files with Windows Performance Analyzer (wpa.exe) for analysis."



Param(
  [string]$Command = 'build-rel\\zx_cli.exe --mode bench --bench-up2 --size 2048 --iters 400',
  [string]$OutDir = 'traces',
  [int]$DurationSeconds = 30,
  [int]$SampleInterval = 1
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

function Get-InstanceNameByPid([int]$procId) {
  $cs = Get-Counter '\\Process(*)\\ID Process'
  foreach ($s in $cs.CounterSamples) {
    if ([int]$s.CookedValue -eq $procId) { return $s.InstanceName }
  }
  throw "Process instance for PID $procId not found in performance counters."
}

New-Item -ItemType Directory -Force -Path $OutDir | Out-Null
$csv = Join-Path $OutDir 'zx_cli_typeperf.csv'

$p = Start-Process -FilePath 'powershell.exe' -ArgumentList @('-NoProfile','-ExecutionPolicy','Bypass','-Command', $Command) -PassThru
Start-Sleep -Seconds 1

$inst = Get-InstanceNameByPid -procId $p.Id
$samples = [Math]::Max(1, [int]([double]$DurationSeconds / [double]$SampleInterval))

# CPU is reported per logical processor; Process % >100 on multi-core machines is expected
$counters = @(
  "\\Process($inst)\\% Processor Time",
  "\\Process($inst)\\Working Set",
  "\\Process($inst)\\Private Bytes",
  "\\Processor(_Total)\\% Processor Time",
  "\\Memory\\Available MBytes"
)

& typeperf $counters -si $SampleInterval -sc $samples -f CSV -o $csv | Out-Null

try {
  if (-not $p.HasExited) { $p.WaitForExit() }
} catch {}

Write-Host "Typeperf CSV written:" (Resolve-Path $csv).Path


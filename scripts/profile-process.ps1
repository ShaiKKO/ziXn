Param(
  [string]$Exe = 'build-rel\\zx_cli.exe',
  [string]$ProcArgs = '--mode bench --bench-up2 --size 2048 --iters 400',
  [string]$Out = 'traces\\zx_cli_process_profile.csv',
  [int]$IntervalMs = 500,
  [int]$MaxMs = 60000
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

New-Item -ItemType Directory -Force -Path (Split-Path $Out) | Out-Null

$psi = New-Object System.Diagnostics.ProcessStartInfo
$psi.FileName = $Exe
$psi.Arguments = $ProcArgs
$psi.UseShellExecute = $false
$psi.RedirectStandardOutput = $true
$psi.RedirectStandardError = $true
$proc = [System.Diagnostics.Process]::Start($psi)

$sw = [System.Diagnostics.Stopwatch]::StartNew()
$rows = New-Object System.Collections.Generic.List[string]
$rows.Add('ms,cpu_ms,working_set,private_bytes,threads') | Out-Null

$last_total = $proc.TotalProcessorTime
Start-Sleep -Milliseconds $IntervalMs
while (-not $proc.HasExited -and $sw.ElapsedMilliseconds -lt $MaxMs) {
  try {
    $proc.Refresh()
    $now = $sw.ElapsedMilliseconds
    $total = $proc.TotalProcessorTime
    $delta_ms = ($total - $last_total).TotalMilliseconds
    $ws = $proc.WorkingSet64
    $pb = $proc.PrivateMemorySize64
    $thr = $proc.Threads.Count
    $rows.Add("$now,$delta_ms,$ws,$pb,$thr") | Out-Null
    $last_total = $total
  } catch {}
  Start-Sleep -Milliseconds $IntervalMs
}

$rows -join "`n" | Set-Content -NoNewline -Path $Out
try { if (-not $proc.HasExited) { $proc.WaitForExit() } } catch {}

Write-Host "Process profile CSV written: $Out"


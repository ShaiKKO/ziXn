param(
  [ValidateSet('nd','det','both')] [string]$Mode = 'both',
  [int]$DurationSec = 45,
  [int]$Iters = 60
)
$ErrorActionPreference = 'Stop'

# Find AMD uProf CLI
$uprof = Get-Command 'AMDuProfCLI.exe' -ErrorAction SilentlyContinue
if (-not $uprof) {
  $uprof = 'C:\\Program Files\\AMD\\uProf\\bin\\AMDuProfCLI.exe'
}
if (-not (Test-Path $uprof)) {
  throw "AMDuProfCLI.exe not found. Install AMD uProf or add it to PATH."
}

# Build Release zx_cli if not present
$root = Split-Path -Parent $MyInvocation.MyCommand.Path
$ws = Resolve-Path (Join-Path $root '..')
Set-Location $ws
if (-not (Test-Path 'build-win')) { cmake -S . -B build-win -G "Visual Studio 17 2022" -A x64 | Out-Null }
cmake --build build-win --config Release --target zx_cli | Out-Null
$exe = Resolve-Path 'build-win/Release/zx_cli.exe'

function Invoke-Profiling([string]$DetLabel){
  $detFlag = if ($DetLabel -eq 'det') { 1 } else { 0 }
  $outDir = Join-Path $ws "uprof_$DetLabel"
  if (Test-Path $outDir) { Remove-Item -Recurse -Force $outDir }
  New-Item -ItemType Directory -Path $outDir | Out-Null
  & "$uprof" collect --config tbp --duration $DurationSec --output $outDir --working-directory (Split-Path $exe) -- "$exe" --mode scene --scene dambreak --deterministic $detFlag --seed 123 | Write-Host
  Write-Host "uProf TBP trace saved to: $outDir"
}

switch ($Mode) {
  'nd'  { Invoke-Profiling 'nd' }
  'det' { Invoke-Profiling 'det' }
  'both' { Invoke-Profiling 'nd'; Invoke-Profiling 'det' }
}

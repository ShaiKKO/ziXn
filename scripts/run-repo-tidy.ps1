Param(
  [string]$BuildDir = 'build-rel',
  [string]$Out = 'ci-errors.txt',
  [string[]]$Focus = @('core/src/zx_api.cpp','core/src/zx_apic_ref.cpp','core/src/zx_residency.cpp'),
  [string]$Checks = '*',
  [switch]$NoQuiet,
  [string[]]$ExtraArgs = @(),
  [string]$ExportFixes = ''
)
$ErrorActionPreference = 'Stop'
$ct = Get-Command clang-tidy -ErrorAction SilentlyContinue
if(-not $ct){ $ct = Get-Item 'C:\Program Files\LLVM\bin\clang-tidy.exe' }
if(-not $ct){ throw 'clang-tidy not found' }
$compdb = Join-Path $BuildDir 'compile_commands.json'
if(-not (Test-Path $compdb)){ throw "compile_commands.json not found in $BuildDir" }
$files = if($Focus -and $Focus.Count -gt 0){ $Focus } else { git ls-files | Where-Object { $_ -match '^(core\\src|core\\include|kernels)\\.*\.(cpp|c|h|hpp)$' } }
$files = $files | Where-Object { $_ -notmatch '^tools\\' }
$errs = New-Object System.Collections.Generic.List[string]
foreach($f in $files){
  try {
    $args = @($f, '-p', $BuildDir, "-checks=$Checks", "-warnings-as-errors=*")
    foreach($ea in $ExtraArgs){ $args += "--extra-arg=$ea" }
    if(-not $NoQuiet){ $args += '-quiet' }
    if($ExportFixes -and $ExportFixes.Length -gt 0){ $args += "-export-fixes=$ExportFixes" }
    $r = & $ct @args 2>&1 | Out-String
  } catch {
    $r = $_ | Out-String
  }
  $errs.Add("==== $f ====") | Out-Null
  $errs.Add($r) | Out-Null
}
$errs -join "`n" | Set-Content -NoNewline -Path $Out
Write-Host "Tidy run complete. See $Out"

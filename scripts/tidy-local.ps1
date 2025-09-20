param(
  [string]$BuildDir = '.build/tidy',
  [int]$Jobs = 8
)

$ErrorActionPreference = 'Stop'

function Test-Exe($name) {
  $cmd = Get-Command $name -ErrorAction SilentlyContinue
  if ($null -eq $cmd) { return $false }
  try { & $name --version 1>$null 2>$null; return $true } catch { return $true }
}

$repoRoot = Resolve-Path (Join-Path $PSScriptRoot '..')
Set-Location -LiteralPath $repoRoot

if (-not (Test-Exe 'clang-format')) {
  Write-Host 'clang-format not found in PATH. Install LLVM (choco install llvm -y) or add to PATH.' -ForegroundColor Yellow
}
if (-not (Test-Exe 'clang-tidy')) {
  Write-Host 'clang-tidy not found in PATH. Install LLVM (choco install llvm -y) or add to PATH.' -ForegroundColor Yellow
}

# Collect tracked files
$tracked = git ls-files
$formatFiles = $tracked | Where-Object { $_ -match '^(core/|kernels/)' } |
  Where-Object { $_ -match '\.(cpp|c|h|hpp|hlsl)$' }

# 1) clang-format (in-place)
if (Test-Exe 'clang-format') {
  Write-Host 'Running clang-format on tracked sources/headers/shaders...'
  foreach ($f in $formatFiles) {
    & clang-format -i -style=file -- "$f"
  }
}

# 2) clang-tidy passes (apply auto-fixes), no compile DB required for many checks
$cppFiles = $tracked | Where-Object { $_ -match '^core/src/.*\.cpp$' }
if ($cppFiles.Count -gt 0 -and (Test-Exe 'clang-tidy')) {
  New-Item -ItemType Directory -Force -Path (Join-Path $repoRoot '.local-tidy-logs') | Out-Null

  $common = @('-format-style=file', '-header-filter', '^(core/(include|src)|kernels)/.*', '-fix', '-fix-errors')
  $pp     = @('--', '-Icore/include')

  # Pass A: braces, bool conversion, math parentheses, uppercase literal suffix, isolate decl, qualified auto
  $checksA = 'readability-braces-around-statements,readability-implicit-bool-conversion,readability-math-missing-parentheses,readability-uppercase-literal-suffix,readability-isolate-declaration,readability-qualified-auto'
  Write-Host "clang-tidy pass A: $checksA"
  & clang-tidy -checks $checksA @common @cppFiles @pp 2>&1 | Tee-Object '.local-tidy-logs/passA.log' | Out-Null

  # Pass B: loop convert, use-std-array, use-nullptr, performance-trivially-destructible
  $checksB = 'modernize-loop-convert,modernize-use-std-array,modernize-use-nullptr,performance-trivially-destructible'
  Write-Host "clang-tidy pass B: $checksB"
  & clang-tidy -checks $checksB @common @cppFiles @pp 2>&1 | Tee-Object '.local-tidy-logs/passB.log' | Out-Null

  # Pass C: member-init, redundant-void-arg, bugprone-easily-swappable-parameters (annotate via NOLINT where needed)
  $checksC = 'cppcoreguidelines-pro-type-member-init,modernize-redundant-void-arg,bugprone-easily-swappable-parameters'
  Write-Host "clang-tidy pass C: $checksC"
  & clang-tidy -checks $checksC @common @cppFiles @pp 2>&1 | Tee-Object '.local-tidy-logs/passC.log' | Out-Null
}

# Summary
git status --porcelain | Out-Host
Write-Host 'Local tidy/format passes complete. Review changes, build locally, and commit when ready.' -ForegroundColor Green



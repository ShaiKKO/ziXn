Param()
$ErrorActionPreference = 'Stop'
$paths = git ls-files | Where-Object { $_ -match '\.(c|cpp|h|hpp)$' }
$found = @()
foreach($p in $paths){
  $m = Select-String -Path $p -Pattern 'NOLINT' -SimpleMatch -Quiet
  if($m){ $found += $p }
}
if($found.Count -gt 0){
  Write-Host "'NOLINT' markers found in:" -ForegroundColor Red
  $found | ForEach-Object { Write-Host "  $_" }
  exit 1
}
Write-Host "No NOLINT markers found." -ForegroundColor Green

<#
/**
 * @file parse-perfview-xml.ps1
 * @brief Parse PerfView .perfView.xml and emit top-N frames as CSV.
 * @details
 *   Reads a PerfView saved view XML, resolves Frames, Stacks, and Samples,
 *   aggregates per-frame metric (CPU_TIME or sample count), and writes a CSV
 *   with columns: Rank, Frame, Samples, Metric, Percent.
 *   Owner: Colin Macritchie / Ripple Group, LLC. Engineering Standards (2025).
 */
#>

[CmdletBinding(PositionalBinding=$false)]
param(
  [Parameter(Mandatory=$true)] [string] $InputXml,
  [string] $OutputCsv = "",
  [int] $TopN = 15
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

function Fail([string] $msg) { Write-Error $msg; exit 1 }

if (-not (Test-Path -LiteralPath $InputXml)) {
  Fail "Input XML not found: $InputXml"
}

[xml]$doc = Get-Content -LiteralPath $InputXml -Raw

$frames = @{}
$frameIdxToName = @{}
$metricFrameId = -1

# Build frame map and detect metric frame id (CPU_TIME, BLOCKED_TIME, etc.)
$frameNodes = $doc.SelectNodes('//StackWindow/StackSource/Frames/Frame')
foreach ($f in $frameNodes) {
  $id = [int]$f.ID
  $name = [string]$f.'#text'
  $frameIdxToName[$id] = $name
  if ($name -eq 'CPU_TIME') {
    $metricFrameId = $id
  }
}

# Build parent links for stacks for quick upward walk
$stackToParent = @{}
$stackToFrame = @{}
$stackNodes = $doc.SelectNodes('//StackWindow/StackSource/Stacks/Stack')
foreach ($s in $stackNodes) {
  $sid = [int]$s.ID
  $caller = [int]$s.CallerID
  $fid = [int]$s.FrameID
  $stackToParent[$sid] = $caller
  $stackToFrame[$sid] = $fid
}

# Aggregate samples by walking up to CPU_TIME and counting the next frame
# Strategy: For each Sample(StackID, Metric?), walk up until we see frame==CPU_TIME.
# Then take its child frame (the next up) as the top-of-interest for accumulation.
function Is-GenericFrameName([string] $name) {
  if (-not $name) { return $true }
  if ($name -eq 'ROOT') { return $true }
  if ($name -eq 'BROKEN') { return $true }
  if ($name -eq '?!?') { return $true }
  if ($name -eq 'OVERHEAD') { return $true }
  if ($name -eq 'CPU_TIME') { return $true }
  if ($name -eq 'BLOCKED_TIME') { return $true }
  if ($name.StartsWith('Thread (')) { return $true }
  if ($name.StartsWith('Process ')) { return $true }
  if ($name.StartsWith('Process64 ')) { return $true }
  if ($name.StartsWith('Processor (')) { return $true }
  if ($name.StartsWith('READIED BY')) { return $true }
  return $false
}

$accum = @{}
$samples = $doc.SelectNodes('//StackWindow/StackSource/Samples/Sample')
foreach ($s in $samples) {
  $sid = [int]$s.StackID
  $metricVal = 1.0
  if ($s.Attributes['Metric']) { $metricVal = [double]$s.Metric }

  $current = $sid
  $candidateFrameId = -1
  $lastNonGenericFrameId = -1

  while ($current -ne -1 -and $stackToFrame.ContainsKey($current)) {
    $fid = $stackToFrame[$current]
    $fname = $frameIdxToName[$fid]

    if (-not (Is-GenericFrameName $fname)) {
      $lastNonGenericFrameId = $fid
    }

    if ($fid -eq $metricFrameId) {
      if ($lastNonGenericFrameId -ne -1) { $candidateFrameId = $lastNonGenericFrameId }
      break
    }

    $current = $stackToParent[$current]
  }

  if ($candidateFrameId -eq -1) {
    if ($lastNonGenericFrameId -ne -1) { $candidateFrameId = $lastNonGenericFrameId }
    else {
      # Fallback to direct leaf frame
      if ($stackToFrame.ContainsKey($sid)) { $candidateFrameId = $stackToFrame[$sid] }
    }
  }

  if ($candidateFrameId -ne -1) {
    if (-not $accum.ContainsKey($candidateFrameId)) { $accum[$candidateFrameId] = @{ Samples = 0; Metric = 0.0 } }
    $accum[$candidateFrameId].Samples += 1
    $accum[$candidateFrameId].Metric += $metricVal
  }
}

# If CPU_TIME not present or no accumulation happened, fall back: attribute top by sample's direct frame
if ($accum.Count -eq 0) {
  foreach ($s in $samples) {
    $sid = [int]$s.StackID
    $metricVal = 1.0
    if ($s.Attributes['Metric']) { $metricVal = [double]$s.Metric }
    $fid = $stackToFrame[$sid]
    if (-not $accum.ContainsKey($fid)) { $accum[$fid] = @{ Samples = 0; Metric = 0.0 } }
    $accum[$fid].Samples += 1
    $accum[$fid].Metric += $metricVal
  }
}

# Compute totals and ranking
$totalMetric = 0.0
foreach ($kv in $accum.GetEnumerator()) { $totalMetric += $kv.Value.Metric }
if ($totalMetric -le 0.0) { $totalMetric = [double]($samples.Count) }

$rows = @()
foreach ($fid in $accum.Keys) {
  $name = $frameIdxToName[$fid]
  if (-not $name) { $name = "Frame_$fid" }
  $metric = [double]$accum[$fid].Metric
  $samplesCnt = [int]$accum[$fid].Samples
  $pct = if ($totalMetric -gt 0.0) { 100.0 * $metric / $totalMetric } else { 0.0 }
  $rows += [pscustomobject]@{
    Frame = $name
    Samples = $samplesCnt
    Metric = '{0:F6}' -f $metric
    Percent = '{0:F2}' -f $pct
  }
}

$rows = $rows | Sort-Object {[double]$_.Metric} -Descending | Select-Object -First $TopN

# Prepare output path
if (-not $OutputCsv -or $OutputCsv -eq '') {
  $dir = Split-Path -Parent (Resolve-Path -LiteralPath $InputXml).Path
  $base = [System.IO.Path]::GetFileNameWithoutExtension($InputXml)
  $OutputCsv = Join-Path $dir ($base + '.top.csv')
}

# Write CSV
"Rank,Frame,Samples,Metric,Percent" | Out-File -LiteralPath $OutputCsv -Encoding UTF8
$rank = 1
foreach ($r in $rows) {
  $line = ($rank.ToString() + ',' + '"' + ($r.Frame -replace '"','''') + '"' + ',' + $r.Samples + ',' + $r.Metric + ',' + $r.Percent)
  $line | Out-File -LiteralPath $OutputCsv -Append -Encoding UTF8
  $rank++
}

Write-Host "[parse-perfview-xml] Wrote: $OutputCsv" -ForegroundColor Cyan

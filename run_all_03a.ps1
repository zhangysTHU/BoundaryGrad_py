param(
  [int]$FromStep = 1,
  [int]$ToStep = 11,
  [string]$SampleName = "CRC1",
  [string]$Python = "python"
)

$ErrorActionPreference = "Stop"
$root = $PSScriptRoot

$previousSampleName = [Environment]::GetEnvironmentVariable("COTTRAZM_SAMPLE_NAME", "Process")
$previousRouteSuffix = [Environment]::GetEnvironmentVariable("COTTRAZM_ROUTE_SUFFIX", "Process")

try {
  [Environment]::SetEnvironmentVariable("COTTRAZM_SAMPLE_NAME", $SampleName, "Process")
  [Environment]::SetEnvironmentVariable("COTTRAZM_ROUTE_SUFFIX", "03a", "Process")

  & (Join-Path $root "run_all.ps1") `
    -FromStep $FromStep `
    -ToStep $ToStep `
    -InfercnvMode "external_r" `
    -SampleName $SampleName `
    -RouteSuffix "03a" `
    -Python $Python
}
finally {
  [Environment]::SetEnvironmentVariable("COTTRAZM_SAMPLE_NAME", $previousSampleName, "Process")
  [Environment]::SetEnvironmentVariable("COTTRAZM_ROUTE_SUFFIX", $previousRouteSuffix, "Process")
}

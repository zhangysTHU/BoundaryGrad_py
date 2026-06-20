param(
  [int]$FromStep = 1,
  [int]$ToStep = 11,
  [ValidateSet("infercnvpy", "external_r")]
  [string]$InfercnvMode = "infercnvpy",
  [string]$SampleName = "CRC1",
  [string]$RouteSuffix = "",
  [string]$Python = "python"
)

$ErrorActionPreference = "Stop"
$root = $PSScriptRoot
$previousSampleName = [Environment]::GetEnvironmentVariable("COTTRAZM_SAMPLE_NAME", "Process")
$previousRouteSuffix = [Environment]::GetEnvironmentVariable("COTTRAZM_ROUTE_SUFFIX", "Process")

if (-not (Test-Path -LiteralPath $Python)) {
  throw "Python not found: $Python"
}

if ([string]::IsNullOrWhiteSpace($RouteSuffix) -and $InfercnvMode -eq "external_r") {
  $RouteSuffix = "03a"
}

[Environment]::SetEnvironmentVariable("COTTRAZM_SAMPLE_NAME", $SampleName, "Process")
[Environment]::SetEnvironmentVariable("COTTRAZM_ROUTE_SUFFIX", $RouteSuffix, "Process")

$defaults = @{
  COTTRAZM_FAST_CNV = ""
  COTTRAZM_INFERCNVPY_N_JOBS = "4"
  COTTRAZM_INFERCNVPY_CHUNKSIZE = "1000"
  COTTRAZM_PROGRESS_SECONDS = "120"
}

foreach ($key in $defaults.Keys) {
  if ([string]::IsNullOrEmpty([Environment]::GetEnvironmentVariable($key, "Process"))) {
    [Environment]::SetEnvironmentVariable($key, $defaults[$key], "Process")
  }
}

$infercnvScript = if ($InfercnvMode -eq "external_r") {
  "03a_run_infercnv_external_r.py"
} else {
  "03b_run_infercnvpy.py"
}

$steps = @(
  @{ Step = 1; Script = "01_preprocess_st.py" },
  @{ Step = 2; Script = "02_morphology_adjusted_cluster.py" },
  @{ Step = 3; Script = $infercnvScript },
  @{ Step = 4; Script = "04_score_cnv.py" },
  @{ Step = 5; Script = "05_define_boundary.py" },
  @{ Step = 6; Script = "06_prepare_single_cell_reference.py" },
  @{ Step = 7; Script = "07_spatial_deconvolution.py" },
  @{ Step = 8; Script = "08_spatial_reconstruction.py" },
  @{ Step = 9; Script = "09_diff_and_enrichment.py" },
  @{ Step = 10; Script = "10_plot_results.py" },
  @{ Step = 11; Script = "11_lsgi_gradient.py" }
)

Push-Location $root
try {
  foreach ($item in $steps) {
    if ($item.Step -lt $FromStep -or $item.Step -gt $ToStep) {
      continue
    }
    Write-Host ("[{0}] Running {1}" -f (Get-Date -Format "yyyy-MM-dd HH:mm:ss"), $item.Script)
    & $Python $item.Script
    if ($LASTEXITCODE -ne 0) {
      throw "$($item.Script) failed with exit code $LASTEXITCODE"
    }
  }
  Write-Host ("[{0}] Python workflow completed." -f (Get-Date -Format "yyyy-MM-dd HH:mm:ss"))
}
finally {
  [Environment]::SetEnvironmentVariable("COTTRAZM_SAMPLE_NAME", $previousSampleName, "Process")
  [Environment]::SetEnvironmentVariable("COTTRAZM_ROUTE_SUFFIX", $previousRouteSuffix, "Process")
  Pop-Location
}

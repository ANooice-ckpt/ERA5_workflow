$ErrorActionPreference = 'Stop'
$maskRoot = $PSScriptRoot
$maskTemp = Join-Path $maskRoot '_tmp'
New-Item -ItemType Directory -Path $maskTemp -Force | Out-Null
$env:TMPDIR = $maskTemp
$env:TMP = $maskTemp
$env:TEMP = $maskTemp
$env:OMP_NUM_THREADS = '1'
$env:OPENBLAS_NUM_THREADS = '1'
$env:MKL_NUM_THREADS = '1'
$env:ERA5_MASK_SCRIPT = Join-Path $maskRoot 'test_population_mask.R'
$rscript = 'C:\Program Files\R\R-4.5.0\bin\Rscript.exe'
& $rscript --vanilla -e 'Sys.setlocale("LC_CTYPE","English_United States.utf8"); source(Sys.getenv("ERA5_MASK_SCRIPT"),encoding="UTF-8")'
if ($LASTEXITCODE -ne 0) { throw "Mask test failed with exit code $LASTEXITCODE" }

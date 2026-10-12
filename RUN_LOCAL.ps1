param([ValidateSet('pack','verify','test')][string]$Mode='pack')
$ErrorActionPreference='Stop'
$bundle=Join-Path $PSScriptRoot '03_server_bundle'
$local=Join-Path $PSScriptRoot '02_local'
$runtime=Join-Path $local 'runtime'
$env:ERA5_BUNDLE_ROOT=$bundle.Replace('\','/')
$env:TMPDIR=Join-Path $runtime 'tmp'
$env:TMP=$env:TMPDIR
$env:TEMP=$env:TMPDIR
$env:R_LIBS_USER=(Join-Path $runtime 'Rlibrary')+';'+(Join-Path $env:LOCALAPPDATA 'R/win-library/4.5')
$env:OMP_NUM_THREADS='1'
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:LC_ALL='C'
$rscript='C:\Program Files\R\R-4.5.0\bin\Rscript.exe'
New-Item -ItemType Directory -Path $env:TMPDIR -Force | Out-Null
if (!(Test-Path -LiteralPath $rscript)) { throw 'R 4.5.0 is not installed at the configured path.' }
$config=Join-Path $bundle 'config.R'
if ($Mode -eq 'test') {
    & $rscript --vanilla (Join-Path $bundle 'code/06_test_regressions.R')
    if ($LASTEXITCODE -ne 0) { throw 'Regression tests failed.' }
    & $rscript --vanilla (Join-Path $bundle 'code/10_test_parallel_runtime.R')
    if ($LASTEXITCODE -ne 0) { throw 'Parallel runtime tests failed.' }
    exit 0
}
if ($Mode -eq 'pack') {
    & $rscript --vanilla (Join-Path $bundle 'code/01_pack_local.R') $config
    if ($LASTEXITCODE -ne 0) { throw 'Local pack failed; inspect the log. Committed years are resumable.' }
}
& $rscript --vanilla (Join-Path $bundle 'code/04_verify_outputs.R') $config pack
if ($LASTEXITCODE -ne 0) { throw 'Weather pack verification failed.' }
& $rscript --vanilla (Join-Path $bundle 'code/08_summarize_pack.R') $config (Join-Path $local 'audit')
if ($LASTEXITCODE -ne 0) { throw 'Pack summary failed.' }

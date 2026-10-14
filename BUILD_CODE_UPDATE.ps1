param([string]$OutputPath='02_local/audit/era5_memory_v2_code_update.tar.gz')
$ErrorActionPreference='Stop'
$projectRoot=[IO.Path]::GetFullPath($PSScriptRoot).TrimEnd('\')
$bundleRoot=Join-Path $projectRoot '03_server_bundle'
$archivePath=[IO.Path]::GetFullPath((Join-Path $projectRoot $OutputPath))
if (!$archivePath.StartsWith($projectRoot+'\',[StringComparison]::OrdinalIgnoreCase)) {
    throw 'The update archive must be inside this project.'
}
$payload=@('run_server.sh','config_memory.R','config_pilot.R','README_SERVER.txt',
    'DEPLOY_MEMORY.txt','SHA256SUMS.txt','WEATHER_MD5SUMS.txt','code/model_environment_versions.csv',
    'code/13_test_server_launcher.sh')
$payload+=@(Get-ChildItem -LiteralPath (Join-Path $bundleRoot 'code') -Filter '*.R' -File |
    ForEach-Object { 'code/'+$_.Name })
if ($payload -contains 'config.R' -or @($payload | Where-Object { $_ -match '^(models|weather|results)/' }).Count) {
    throw 'Unexpected scientific configuration or input payload.'
}
foreach ($entry in Get-Content -LiteralPath (Join-Path $bundleRoot 'SHA256SUMS.txt')) {
    if ($entry -notmatch '^([0-9a-f]{64})\s+(.+)$') { throw 'Invalid static checksum entry.' }
    $expected=$Matches[1]; $relativePath=$Matches[2]
    if ($relativePath -in $payload) {
        $actual=(Get-FileHash -LiteralPath (Join-Path $bundleRoot $relativePath) -Algorithm SHA256).Hash.ToLowerInvariant()
        if ($expected -ne $actual) { throw "Update SHA256SUMS.txt for $relativePath before building." }
    }
}
$stageRoot=Join-Path $projectRoot ('02_local/audit/code_update_stage_'+[guid]::NewGuid().ToString('N'))
New-Item -ItemType Directory -Path $stageRoot,(Split-Path -Parent $archivePath) -Force | Out-Null
try {
    foreach ($relativePath in $payload) {
        $destination=Join-Path $stageRoot $relativePath
        New-Item -ItemType Directory -Path (Split-Path -Parent $destination) -Force | Out-Null
        Copy-Item -LiteralPath (Join-Path $bundleRoot $relativePath) -Destination $destination
    }
    $hashLines=@($payload | Sort-Object | ForEach-Object {
        (Get-FileHash -LiteralPath (Join-Path $stageRoot $_) -Algorithm SHA256).Hash.ToLowerInvariant()+'  '+$_
    })
    [IO.File]::WriteAllText((Join-Path $stageRoot 'UPDATE_SHA256SUMS.txt'),
        ($hashLines -join "`n")+"`n",[Text.UTF8Encoding]::new($false))
    & tar -czf $archivePath -C $stageRoot -- @payload UPDATE_SHA256SUMS.txt
    if ($LASTEXITCODE -ne 0) { throw 'Code update archive failed.' }
    Get-Item -LiteralPath $archivePath | Select-Object FullName,Length
} finally {
    $resolvedStage=[IO.Path]::GetFullPath($stageRoot).TrimEnd('\')
    if ($resolvedStage.StartsWith($projectRoot+'\02_local\audit\code_update_stage_',[StringComparison]::OrdinalIgnoreCase)) {
        Remove-Item -LiteralPath $resolvedStage -Recurse -Force
    } else { throw 'Unexpected staging cleanup path.' }
}

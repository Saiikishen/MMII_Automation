[CmdletBinding()]
param(
    [Parameter(Mandatory = $true, Position = 0)]
    [ValidateNotNullOrEmpty()]
    [string]$Images,

    [ValidateNotNullOrEmpty()]
    [string]$PythonExe = 'python',

    [ValidateNotNullOrEmpty()]
    [string]$Config = (Join-Path $PSScriptRoot 'mmii_config.json'),

    [string]$Output,
    [switch]$WithQC,
    [switch]$CheckInputs
)

$ErrorActionPreference = 'Stop'
$imageItem = Get-Item -LiteralPath $Images
if (-not $imageItem.PSIsContainer -or $imageItem.PSProvider.Name -ne 'FileSystem') {
    throw '-Images must name a filesystem folder containing t1.nii and the available modality files.'
}
$imagePath = $imageItem.FullName
$configPath = (Get-Item -LiteralPath $Config).FullName
$builderPath = Join-Path $PSScriptRoot 'build_mmii.py'
$builderArguments = @($builderPath, '--config', $configPath, '--images', $imagePath)

if ($CheckInputs) {
    & $PythonExe @builderArguments --check-inputs
    if ($LASTEXITCODE -ne 0) {
        throw "MMII input check failed (exit code $LASTEXITCODE). See the message above."
    }
    return
}

if ([string]::IsNullOrWhiteSpace($Output)) {
    $patientName = $imageItem.Name -replace '(?i)_Images$', ''
    if ([string]::IsNullOrWhiteSpace($patientName)) { $patientName = 'patient' }
    $outputName = $patientName + '_MMII_' + (Get-Date -Format 'yyyyMMdd_HHmmss_fff')
    $Output = Join-Path $imageItem.Parent.FullName $outputName
} else {
    $Output = $ExecutionContext.SessionState.Path.GetUnresolvedProviderPathFromPSPath($Output)
}
$builderArguments += @('--output', $Output)
& $PythonExe @builderArguments
if ($LASTEXITCODE -ne 0) {
    throw "MMII build failed (exit code $LASTEXITCODE). See the message above and any output logs."
}
if ($WithQC) {
    & $PythonExe (Join-Path $PSScriptRoot 'qc_registration.py') --output $Output
    if ($LASTEXITCODE -ne 0) {
        throw "Workspace built, but QC generation failed (exit code $LASTEXITCODE)."
    }
}
Write-Host ('Workspace: ' + (Join-Path $Output 'MMII.itksnap'))

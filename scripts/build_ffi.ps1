# Explicit distribution build for rssn C/C++ FFI libraries on Windows; keeps
# ordinary Cargo dependency builds lean.
[CmdletBinding()]
param(
    [string]$Target = "",
    [string]$Profile = $(if ($env:RSSN_FFI_PROFILE) { $env:RSSN_FFI_PROFILE } else { "release" }),
    [switch]$Package,
    [string]$DistName = "",
    [string]$DistDir = "",
    [string]$Version = "",
    [Parameter(ValueFromRemainingArguments = $true)]
    [string[]]$CargoArgs = @()
)

$ErrorActionPreference = "Stop"

$RssnRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
if ([string]::IsNullOrWhiteSpace($DistDir)) {
    $DistDir = Join-Path $RssnRoot "dist"
}

$ManifestPath = Join-Path $RssnRoot "Cargo.toml"
$CmdArgs = @(
    "rustc",
    "--manifest-path", $ManifestPath,
    "--lib",
    "--profile", $Profile,
    "--crate-type", "staticlib,cdylib"
)
if (-not [string]::IsNullOrWhiteSpace($Target)) {
    $CmdArgs += @("--target", $Target)
}
if ($CargoArgs.Count -gt 0) {
    $CmdArgs += $CargoArgs
}

& cargo @CmdArgs
if ($LASTEXITCODE -ne 0) {
    exit $LASTEXITCODE
}

if (-not $Package -and [string]::IsNullOrWhiteSpace($DistName)) {
    exit 0
}

if (-not [string]::IsNullOrWhiteSpace($Target)) {
    $ArtifactDir = Join-Path (Join-Path (Join-Path $RssnRoot "target") $Target) $Profile
} else {
    $ArtifactDir = Join-Path (Join-Path $RssnRoot "target") $Profile
}

if ([string]::IsNullOrWhiteSpace($DistName)) {
    if ([string]::IsNullOrWhiteSpace($Version)) {
        $CargoLines = Get-Content -Path $ManifestPath -Encoding UTF8
        $InPackage = $false
        foreach ($Line in $CargoLines) {
            if ($Line -match '^\[package\]') {
                $InPackage = $true
                continue
            }
            if ($InPackage -and $Line -match '^\[') {
                break
            }
            if ($InPackage -and $Line -match '^version\s*=\s*"([^"]+)"') {
                $Version = $Matches[1]
                break
            }
        }
    }
    if ([string]::IsNullOrWhiteSpace($Version)) {
        $Version = "0.0.0"
    }
    $TargetLabel = if ([string]::IsNullOrWhiteSpace($Target)) { "host" } else { $Target }
    $DistName = "rssn-$Version-$TargetLabel"
}

$StageDir = Join-Path $DistDir $DistName
if (Test-Path $StageDir) {
    Remove-Item -Path $StageDir -Recurse -Force
}
New-Item -Path $StageDir -ItemType Directory -Force | Out-Null

$Copied = 0
foreach ($Lib in @("rssn.dll", "rssn.dll.lib", "rssn.lib", "librssn.a", "librssn.so", "librssn.dylib")) {
    $Src = Join-Path $ArtifactDir $Lib
    if (Test-Path $Src) {
        Copy-Item -Path $Src -Destination $StageDir -Force
        $Copied++
    }
}

if ($Copied -eq 0) {
    throw "No FFI library artifacts found in $ArtifactDir"
}

foreach ($Doc in @("rssn.h", "rssn.hpp", "LICENSE", "README.md")) {
    $Src = Join-Path $RssnRoot $Doc
    if (Test-Path $Src) {
        Copy-Item -Path $Src -Destination $StageDir -Force
    }
}

$ArchiveName = "$DistName.tar.gz"
$ArchivePath = Join-Path $DistDir $ArchiveName
$ShaPath = Join-Path $DistDir "$ArchiveName.sha256"

Push-Location $DistDir
try {
    & tar -czf $ArchiveName $DistName
    if ($LASTEXITCODE -ne 0) {
        throw "tar failed with exit code $LASTEXITCODE"
    }
    $Hash = (Get-FileHash -Path $ArchivePath -Algorithm SHA256).Hash.ToLowerInvariant()
    Set-Content -Path $ShaPath -Value "$Hash  $ArchiveName" -Encoding Ascii -NoNewline
    Add-Content -Path $ShaPath -Value "" -Encoding Ascii
} finally {
    Pop-Location
}

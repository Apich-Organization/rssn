# Explicit distribution build for rssn C/C++ FFI libraries on Windows; keeps
# ordinary Cargo dependency builds lean.
[CmdletBinding()]
param(
    [string]$Target = "",
    [string]$Profile = $(if ($env:RSSN_FFI_PROFILE) { $env:RSSN_FFI_PROFILE } else { "ffi" }),
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

$ArtifactProfile = switch ($Profile) {
    "dev" { "debug" }
    "test" { "debug" }
    "bench" { "release" }
    default { $Profile }
}
$ArtifactRoot = if ($env:CARGO_TARGET_DIR) {
    $env:CARGO_TARGET_DIR
} else {
    Join-Path $RssnRoot "target"
}
if (-not [string]::IsNullOrWhiteSpace($Target)) {
    $ArtifactRoot = Join-Path $ArtifactRoot $Target
}
$ArtifactDir = Join-Path $ArtifactRoot $ArtifactProfile

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

$Documents = @(
    "rssn.h", "rssn.hpp", "LICENSE", "ARCHITECTURE.md", "ATTRIBUTIONS.md",
    "CHANGELOG.md", "CODE_OF_CONDUCT.md", "CONTRIBUTING.md", "SECURITY.md"
)
foreach ($Document in $Documents) {
    Copy-Item -LiteralPath (Join-Path $RssnRoot $Document) -Destination $StageDir
}
Get-ChildItem -LiteralPath (Join-Path $PSScriptRoot "ffi-package") -Filter "README*.md" -File |
    ForEach-Object { Copy-Item -LiteralPath $_.FullName -Destination $StageDir }
$DocsDir = Join-Path $StageDir "docs"
$ExamplesDir = Join-Path $StageDir "examples"
$CExamplesDir = Join-Path $ExamplesDir "c"
New-Item -Path $DocsDir, $ExamplesDir, $CExamplesDir -ItemType Directory | Out-Null
Get-ChildItem -LiteralPath (Join-Path $RssnRoot "docs") -Filter "*.md" -File |
    ForEach-Object { Copy-Item -LiteralPath $_.FullName -Destination $DocsDir }
Get-ChildItem -LiteralPath (Join-Path $RssnRoot "examples") -Filter "*.rs" -File |
    ForEach-Object { Copy-Item -LiteralPath $_.FullName -Destination $ExamplesDir }
Copy-Item -LiteralPath (Join-Path $RssnRoot "examples/README.md") -Destination $ExamplesDir
Get-ChildItem -LiteralPath (Join-Path $RssnRoot "examples/c") -Filter "*.c" -File |
    ForEach-Object { Copy-Item -LiteralPath $_.FullName -Destination $CExamplesDir }

# Match build_ffi.sh: use the checkout commit for files outside the archive.
$SourceRef = "v0.3"
if ((Test-Path -LiteralPath (Join-Path $RssnRoot ".git")) -and (Get-Command git -ErrorAction SilentlyContinue)) {
    $Revision = & git -C $RssnRoot rev-parse --verify HEAD 2>$null
    if ($LASTEXITCODE -eq 0) {
        $SourceRef = $Revision
    }
}
$SourceRepository = "https://github.com/Apich-Organization/rssn"
$LinkPrefix = "$SourceRepository/blob/$SourceRef/"
$Utf8 = [System.Text.UTF8Encoding]::new($false)
Get-ChildItem -LiteralPath $StageDir -Filter "*.md" -Recurse -File | ForEach-Object {
    $Text = [System.IO.File]::ReadAllText($_.FullName)
    $Text = $Text.Replace("@RSSN_SOURCE@", "$SourceRepository/tree/$SourceRef")
    $Text = $Text -creplace '(\]\()(\.\./)?((src|tools)/[^\s)]+)', ('$1' + $LinkPrefix + '$3')
    $Text = $Text.Replace("](README.md#testing)", "](${LinkPrefix}README.md#testing)")
    [System.IO.File]::WriteAllText($_.FullName, $Text, $Utf8)
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

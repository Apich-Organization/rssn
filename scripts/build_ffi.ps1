# Compatibility entry point for the RSSN Bazel FFI targets.
[CmdletBinding()]
param(
    [string]$Target = "",
    [string]$Profile = $(if ($env:RSSN_FFI_PROFILE) { $env:RSSN_FFI_PROFILE } else { "ffi" }),
    [switch]$Package,
    [switch]$LlvmOverlay,
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

$BazelCommand = Get-Command bazel, bazelisk -CommandType Application -ErrorAction SilentlyContinue |
    Select-Object -First 1
if (-not $BazelCommand) { throw "Install Bazel or Bazelisk to build FFI libraries." }
$BazelExecutable = $BazelCommand.Source
$BazelItem = Get-Item -LiteralPath $BazelExecutable
if ($BazelItem.LinkType -eq "SymbolicLink") {
    $LinkTarget = @($BazelItem.Target)[0]
    if (-not [System.IO.Path]::IsPathRooted($LinkTarget)) {
        $LinkTarget = Join-Path $BazelItem.DirectoryName $LinkTarget
    }
    $BazelExecutable = (Resolve-Path -LiteralPath $LinkTarget).Path
}

$BazelFlags = @()
foreach ($Variable in @("RSSN_TARGET_ENV", "RUSTC", "CARGO", "RSSN_LLVM_CONFIG", "RSSN_LLVM_IMPORT_LIBRARY")) {
    $Value = [Environment]::GetEnvironmentVariable($Variable)
    if ($Value) { $BazelFlags += "--repo_env=$Variable=$Value" }
}
if ($Target) { $BazelFlags += "--define=rust_target=$Target" }
if (-not [System.IO.Path]::IsPathRooted($DistDir)) {
    $DistDir = Join-Path (Get-Location) $DistDir
}
$DistFlags = @("--dist-dir", $DistDir)
if ($DistName) { $DistFlags += @("--dist-name", $DistName) }
if ($Version) { $DistFlags += @("--version", $Version) }

Push-Location $RssnRoot
try {
    if ($Profile -in @("ffi", "release")) {
        foreach ($Arg in $CargoArgs) {
            if ($Arg -notin @("--locked", "--all-features")) {
                throw "Unsupported release Cargo argument $Arg."
            }
        }
        if ($Package -or $DistName) {
            $ArchiveTarget = if ($LlvmOverlay) { "//:package_ffi_llvm" } else { "//:package_ffi" }
            & $BazelExecutable run @BazelFlags $ArchiveTarget -- @DistFlags
        } else {
            & $BazelExecutable build @BazelFlags //:ffi_distribution
        }
    } elseif ($Profile -in @("dev", "test")) {
        if ($LlvmOverlay) { throw "The LLVM option requires profile ffi or release." }
        if ($Package -or $DistName) { $DistFlags += "--package" }
        & $BazelExecutable run @BazelFlags //:ffi_legacy_dev -- --profile $Profile @DistFlags @CargoArgs
    } else {
        throw "Use profile ffi, release, dev, or test."
    }
    exit $LASTEXITCODE
} finally {
    Pop-Location
}

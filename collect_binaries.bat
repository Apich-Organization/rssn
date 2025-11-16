@echo off
setlocal enabledelayedexpansion

:: --- Configuration ---

:: 1. Output directory name
set OUTPUT_DIR=dist

:: 2. Project Name / Executable Name (PLEASE CHANGE THIS MANUALLY)
set EXECUTABLE_NAME=my_rust_project
:: !!! ⚠️ Change the line above to match the [package].name in your Cargo.toml !!!
echo ----------------------------------------------------
echo Please verify EXECUTABLE_NAME (Current: %EXECUTABLE_NAME%)
echo ----------------------------------------------------
pause

:: 3. Target List
set ALL_TARGETS=^
x86_64-apple-darwin ^
aarch64-apple-darwin ^
x86_64-pc-windows-msvc ^
aarch64-pc-windows-msvc ^
i686-pc-windows-msvc ^
x86_64-unknown-linux-gnu ^
aarch64-unknown-linux-gnu ^
x86_64-unknown-linux-musl ^
armv7-unknown-linux-gnueabihf ^
i686-unknown-linux-gnu ^
x86_64-unknown-freebsd ^
x86_64-unknown-netbsd ^
aarch64-linux-android ^
x86_64-linux-android ^
wasm32-unknown-unknown

:: --- Initialization ---
echo --- 📦 Collecting Release Artifacts into %OUTPUT_DIR%\ directory ---
if not exist "%OUTPUT_DIR%" mkdir "%OUTPUT_DIR%"

:: --- Loop through targets ---
for %%T in (%ALL_TARGETS%) do (
    set "TARGET=%%T"
    set "EXT="

    :: Determine file extension
    if "!TARGET:windows=!" NEQ "!TARGET!" (
        set "EXT=.exe"
    ) else if "!TARGET:wasm=!" NEQ "!TARGET!" (
        set "EXT=.wasm"
    )

    set "SOURCE_FILE=target\%%T\release\%EXECUTABLE_NAME%!EXT!"
    set "DEST_FILE=%OUTPUT_DIR%\%EXECUTABLE_NAME%-%%T!EXT!"

    if exist "!SOURCE_FILE!" (
        copy "!SOURCE_FILE!" "!DEST_FILE!" > NUL
        echo ✅ Copied successfully: %%T -> %OUTPUT_DIR%\%EXECUTABLE_NAME%-%%T!EXT!
    ) else (
        echo ⚠️ Warning: File not found at "!SOURCE_FILE!" (Not compiled or failed compilation).
    )
)

echo.
echo === 🎉 All release artifacts collected (%OUTPUT_DIR%) ===

endlocal
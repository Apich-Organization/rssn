@echo off
setlocal enabledelayedexpansion

echo --- Rust Cross-Platform Operations Script ---
echo Please select the operation type:
echo 1) Debug Build (cargo build)
echo 2) Release Build (cargo build --release)
echo 3) Test (cargo test)
echo 4) Bench (cargo bench)
set /p CHOICE="Enter option (1-4): "

set CARGO_COMMAND=
set MODE_DESCRIPTION=

if "%CHOICE%"=="1" (
    set CARGO_COMMAND=build --features full
    set MODE_DESCRIPTION=Debug Build
) else if "%CHOICE%"=="2" (
    set CARGO_COMMAND=build --release --features full
    set MODE_DESCRIPTION=Release Build
) else if "%CHOICE%"=="3" (
    set CARGO_COMMAND=test --features full
    set MODE_DESCRIPTION=Test
) else if "%CHOICE%"=="4" (
    set CARGO_COMMAND=bench --features full
    set MODE_DESCRIPTION=Bench
) else (
    echo Invalid option. Exiting.
    goto :EOF
)

echo.
echo You selected: %MODE_DESCRIPTION%
echo.

:: --- Configure Targets ---
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

:: List of targets that cannot execute 'cargo test' or 'cargo bench'
set NO_TEST_OR_BENCH=^
wasm32-unknown-unknown ^
aarch64-linux-android ^
x86_64-linux-android

:: --- Execution Loop ---
for %%T in (%ALL_TARGETS%) do (
    set "TARGET=%%T"
    
    set "SKIP=0"
    if "%CARGO_COMMAND%"=="test" (
        for %%N in (%NO_TEST_OR_BENCH%) do (
            if "!TARGET!"=="%%N" (
                echo --- ⚠️ Skipping !TARGET!: test is not supported for this target.---
                set "SKIP=1"
            )
        )
    )
    if "%CARGO_COMMAND%"=="bench" (
        for %%N in (%NO_TEST_OR_BENCH%) do (
            if "!TARGET!"=="%%N" (
                echo --- ⚠️ Skipping !TARGET!: bench is not supported for this target.---
                set "SKIP=1"
            )
        )
    )

    if "!SKIP!"=="0" (
        echo ================================================================
        echo 🚀 Running 'cargo %CARGO_COMMAND%' for Target: !TARGET!
        echo ================================================================
        
        :: Execute Cargo Command
        cargo %CARGO_COMMAND% --target "!TARGET!"
        
        if errorlevel 1 (
            echo !!! ❌ !TARGET! Build/Test failed. Continuing to the next target...
        ) else (
            echo ✅ !TARGET! execution successful.
        )
    )
)

echo.
echo === All targets operations completed (%MODE_DESCRIPTION%) ===

endlocal
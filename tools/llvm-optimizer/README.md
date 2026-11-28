# RSSN LLVM Optimizer

[English](README.md) | [简体中文](README.zh-CN.md)

Explicitly loaded LLVM New Pass Manager plugin for RSSN release builds. Validated with `rustc` 1.100.0 (`0ef1bc386059596dada8fb65b50425aab2ea637a`) and LLVM 23.1.5, running at the end of the `O3` optimization pipeline for fat LTO (and verified with ThinLTO and no-LTO on Linux x86-64).

## Build and Test

Run from the RSSN repository root. Bazel builds the host plugin using the LLVM SDK matching `rustc`, a C++17 compiler, and CMake. Tests require Python 3.11 and SDK `opt` and `llc`.

```bash
rustc -vV
llvm-config-23 --version
bazel build --repo_env=RSSN_LLVM_CONFIG="$(command -v llvm-config-23)" //:llvm_optimizer
bazel test --repo_env=RSSN_LLVM_CONFIG="$(command -v llvm-config-23)" //:llvm_optimizer_test
bazel test --config=release-max --repo_env=RSSN_LLVM_CONFIG="$(command -v llvm-config-23)" //:release_pass_test
```

- `//:llvm_optimizer_tools_test` runs standalone `opt`/`llc` IR, MIR, and object-file tests (usable with `-DRSSN_TEST_RUSTC=OFF` when the host `rustc` lacks exported LLVM C++ symbols, such as official Windows x64 or ARMv7 `rustc`).
- `//:llvm_optimizer_rustc_test` verifies `rustc` plugin loading across no-LTO, ThinLTO, and fat LTO.
- For Windows builds, set `RSSN_LLVM_IMPORT_LIBRARY` to the host LLVM C++ import library.
- For ELF linker-plugin LTO (e.g. ARMv7 with `lld`):

```bash
plugin="$PWD/bazel-bin/llvm_optimizer/RssnLLVMOptimizer.so"
RSSN_SPILL_MODE=off RSSN_LIBM_UNIT_LOG=transform \
  rustc -O -Clinker-plugin-lto -Clink-arg=-fuse-ld=lld \
  "-Clink-arg=-Wl,--load-pass-plugin=$plugin" -Clink-arg=-Wl,--lto-O3 main.rs
```

To run the full release pass pipeline on already-optimized IR or fat LTO bitcode:

```bash
opt -load-pass-plugin="$plugin" -passes='rssn-release,verify' input.bc -o out.bc
```

## Passes and Environment Controls

All optional passes default to `off` when loading the plugin directly unless configured by the Bazel release runner (`bazel/runner.py`). Set `RSSN_SPILL_FILTER=<substring>` to restrict passes to matching LLVM symbol names.

| Environment Variable | Modes | `opt` / `llc` Pass Name | Purpose |
|---|---|---|---|
| `RSSN_SPILL_MODE` | `off` \| `report` \| `transform` | `function(rssn-spill-pressure)` | Joint live-range motion, cross-block rematerialization, and load fusion/merging. |
| `RSSN_SPILL_PLAN_WIDTH` | `1`..`4` (default `4`) | — | Max members moved together in a joint motion plan. |
| `RSSN_SPILL_ABI` | `off` \| `report` \| `transform` | `rssn-call-abi` | SysV local call `x86_regcallcc` rewriting and indirect scalar field load detection. |
| `RSSN_SPILL_LIBM` | `off` \| `report` \| `transform` | `function(rssn-libm-pressure)` | Vector `exp`/`log` branch recovery and SysV scalar libm live-range motion. |
| `RSSN_SPILL_LIBM_REUSE` | `off` \| `report` \| `transform` | `function(rssn-libm-reuse)` | Bitwise-checked reuse of mirrored `exp` calls within a loop (x86-64 SysV & Win64). |
| `RSSN_LIBM_UNIT_LOG` | `off` \| `report` \| `transform` | `function(rssn-libm-unit-log)` | Reuse `exp(C*log(x)-x)` at unit nodes in exp-sinh quadrature. |
| `RSSN_LIBM_TAIL_CUTOFF` | `off` \| `report` \| `transform` | `function(rssn-libm-tail-cutoff)` | Early cutoff / zero-numerator proof for far exp-sinh quadrature nodes. |
| `RSSN_POW_EXPONENT_PACK` | `off` \| `report` \| `transform` | `function(rssn-pow-exponent-pack)` | Pack loop-invariant exponents across adjacent `pow` calls on Win64 MSVC AVX. |
| `RSSN_LOOP_EXIT_COMPRESSION` | `off` \| `report` \| `transform` | `function(rssn-loop-exit-compression)` | Precompute pure FP expressions at loop exit across calls (x86-64 SysV & Win64). |
| `RSSN_LOOP_PREHEADER_PRESSURE` | `off` \| `report` \| `transform` | `function(rssn-loop-preheader-pressure)` | Sink branch-only integer/address/vector builds from loop preheaders (x86-64 & AArch64). |
| `RSSN_FP_ENDPOINT_GUARDS` | `off` \| `report` \| `transform` | `function(rssn-fp-endpoint-guards)` | Merge endpoint/weight checks in double-exponential quadrature to skip dead `llvm.exp`. |
| `RSSN_FP_CLASS_INTERVAL` | `off` \| `report` \| `transform` | `function(rssn-fp-class-interval)` | Lower binary64 finite-nonzero branch checks in loops to a single unsigned interval check. |
| `RSSN_SCALAR_LOAD_PACK` | `off` \| `report` \| `transform` | `function(rssn-scalar-load-pack)` | Merge adjacent narrow integer loads in high-GP-pressure loops (4B Linux, 8B Win64). |
| `RSSN_SWITCH_LOAD_REUSE` | `off` \| `report` \| `transform` | `function(rssn-switch-load-reuse)` | Merge wide loads across switch/branch arms and reuse low bytes for narrow loads. |
| `RSSN_FP_REGISTER_BANK` | `off` \| `report` \| `transform` | `function(rssn-fp-register-bank)` | Carry scalar `double` values across libm calls in non-volatile GPRs on Win64 AVX. |
| `RSSN_POW_UNIT_BASE` | `off` \| `report` \| `transform` | `function(rssn-pow-unit-base)` | Skip `llvm.pow` groups when the base is `+1.0` within rounding bounds. |
| — | MIR (`llc -run-pass`) | `rssn-edge-restore` | Post-RA x86-64 MIR pass merging or sinking identical callee-saved restore sequences. |
| `RSSN_MACHINE_STACK_DETAIL` | `0` \| `1` | `rssn-machine-stack-traffic` | Post-RA MIR diagnostic reporting spill, fixed-slot, and constant-pool stack traffic. |

Append `-report` to any IR or `rssn-edge-restore` pass name for analysis-only remarks (`-pass-remarks-analysis=<pass>`).

## Architecture Overview

- **SSA Liveness & Cost Model** (`Pressure.cpp`, `Planning.h`): Computes backward `LiveBefore = Use ∪ (LiveAfter − Def)` across GP and SIMD register banks using TTI legalization counts, ABI preserved-register budgets, and bounded quadratic over-limit penalties.
- **Joint Motion & Rematerialization** (`Motion.cpp`): Bucketed beam search (width 8, up to 12 candidates per bucket, up to 4 moved members) evaluating net live-range deltas across calls and local regions, plus bounded cross-block rematerialization (up to 3 same-depth blocks).
- **Load Merging & Fusion** (`LoadFusion.cpp`, `MemoryModel.cpp`): Eliminates duplicate-address loads, fuses adjacent loads into 128/256-bit vectors, and builds 128-bit vector carriers for scalars that must cross modifying calls—strictly gated by alias analysis, `nounwind`/`nofree`/`nosync`, and non-volatile/non-atomic checks.

## Benchmark Reproduction

```bash
bazel run --repo_env=RSSN_LLVM_CONFIG="$(command -v llvm-config-23)" \
  //:llvm_benchmark -- --cpu 3
```

Supports `--opt`, `--llc`, `--jobs`, `--pairs`, `--work-dir`, and `--cpu`. Runs interleaved A/B and B/A builds in isolated directories with whole-array bitwise hash verification and LLVM `Spill`/`Reload` remark accounting.

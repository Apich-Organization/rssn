# RSSN LLVM 优化器

[English](README.md) | [简体中文](README.zh-CN.md)

面向 RSSN release 构建的显式加载 LLVM New Pass Manager 插件。基于 `rustc` 1.100.0（`0ef1bc386059596dada8fb65b50425aab2ea637a`）与 LLVM 23.1.5 验证，在 `O3` 末尾阶段运行以适配 fat LTO（已在 Linux x86-64 上通过 ThinLTO 与无 LTO 验证）。

## 构建与测试

在 RSSN 仓库根目录运行。Bazel 使用与 `rustc` 匹配的 LLVM SDK、C++17 编译器和 CMake 构建宿主插件。测试需要 Python 3.11 以及 SDK 的 `opt` 和 `llc`。

```bash
rustc -vV
llvm-config-23 --version
bazel build --repo_env=RSSN_LLVM_CONFIG="$(command -v llvm-config-23)" //:llvm_optimizer
bazel test --repo_env=RSSN_LLVM_CONFIG="$(command -v llvm-config-23)" //:llvm_optimizer_test
bazel test --config=release-max --repo_env=RSSN_LLVM_CONFIG="$(command -v llvm-config-23)" //:release_pass_test
```

- `//:llvm_optimizer_tools_test` 运行独立的 `opt`/`llc` IR、MIR 与目标对象测试（当宿主 `rustc` 未导出 LLVM C++ 符号时，如官方 Windows x64 或 ARMv7 `rustc`，可配合 `-DRSSN_TEST_RUSTC=OFF` 使用）。
- `//:llvm_optimizer_rustc_test` 验证 `rustc` 在无 LTO、ThinLTO 与 fat LTO 下的插件加载。
- Windows 构建需将 `RSSN_LLVM_IMPORT_LIBRARY` 指向宿主 LLVM C++ 导入库。
- ELF 链接器插件 LTO 用法（例如 ARMv7 配合 `lld`）：

```bash
plugin="$PWD/bazel-bin/llvm_optimizer/RssnLLVMOptimizer.so"
RSSN_SPILL_MODE=off RSSN_LIBM_UNIT_LOG=transform \
  rustc -O -Clinker-plugin-lto -Clink-arg=-fuse-ld=lld \
  "-Clink-arg=-Wl,--load-pass-plugin=$plugin" -Clink-arg=-Wl,--lto-O3 main.rs
```

对已完成优化的 IR 或 fat LTO bitcode 显式运行完整发布 pass 管线：

```bash
opt -load-pass-plugin="$plugin" -passes='rssn-release,verify' input.bc -o out.bc
```

## Pass 与环境变量控制

直接加载插件时，所有可选 pass 默认为 `off`（Bazel 发布配置通过 `bazel/runner.py` 启用已验证的变换）。设置 `RSSN_SPILL_FILTER=<子串>` 可按 LLVM 符号名过滤目标函数。

| 环境变量 | 模式 | `opt` / `llc` Pass 名称 | 用途 |
|---|---|---|---|
| `RSSN_SPILL_MODE` | `off` \| `report` \| `transform` | `function(rssn-spill-pressure)` | 联合活跃区间移动、跨块重物化与读取合并/融合。 |
| `RSSN_SPILL_PLAN_WIDTH` | `1`..`4`（默认 `4`） | — | 单次联合移动计划允许的最大成员数。 |
| `RSSN_SPILL_ABI` | `off` \| `report` \| `transform` | `rssn-call-abi` | SysV 本地调用 `x86_regcallcc` 改写与间接标量字段读取检测。 |
| `RSSN_SPILL_LIBM` | `off` \| `report` \| `transform` | `function(rssn-libm-pressure)` | 向量 `exp`/`log` 分支恢复与 SysV 标量 libm 活跃区间移动。 |
| `RSSN_SPILL_LIBM_REUSE` | `off` \| `report` \| `transform` | `function(rssn-libm-reuse)` | 循环内镜像 `exp` 调用的按位校验复用（x86-64 SysV 与 Win64）。 |
| `RSSN_LIBM_UNIT_LOG` | `off` \| `report` \| `transform` | `function(rssn-libm-unit-log)` | 在 exp-sinh 求积单位节点复用 `exp(C*log(x)-x)`。 |
| `RSSN_LIBM_TAIL_CUTOFF` | `off` \| `report` \| `transform` | `function(rssn-libm-tail-cutoff)` | 为 exp-sinh 远端求积节点提前证明截断或零分子。 |
| `RSSN_POW_EXPONENT_PACK` | `off` \| `report` \| `transform` | `function(rssn-pow-exponent-pack)` | 为 Win64 MSVC AVX 相邻 `pow` 打包循环不变指数。 |
| `RSSN_LOOP_EXIT_COMPRESSION` | `off` \| `report` \| `transform` | `function(rssn-loop-exit-compression)` | 提前计算循环退出处的纯浮点表达式以减少跨调用活跃值。 |
| `RSSN_LOOP_PREHEADER_PRESSURE` | `off` \| `report` \| `transform` | `function(rssn-loop-preheader-pressure)` | 将循环前置块中分支专用的整数/地址/向量构造下沉到目标分支（x86-64 与 AArch64）。 |
| `RSSN_FP_ENDPOINT_GUARDS` | `off` \| `report` \| `transform` | `function(rssn-fp-endpoint-guards)` | 合并双指数求积端点与权重检查，跳过无效节点的 `llvm.exp`。 |
| `RSSN_FP_CLASS_INTERVAL` | `off` \| `report` \| `transform` | `function(rssn-fp-class-interval)` | 将循环内 binary64 有限非零检查降低为单次无符号区间判断。 |
| `RSSN_SCALAR_LOAD_PACK` | `off` \| `report` \| `transform` | `function(rssn-scalar-load-pack)` | 在高 GP 压力循环中合并相邻窄整数读取（Linux 4 字节，Win64 8 字节）。 |
| `RSSN_SWITCH_LOAD_REUSE` | `off` \| `report` \| `transform` | `function(rssn-switch-load-reuse)` | 合并多分支共用的宽读取并复用低位字节替代窄读取。 |
| `RSSN_FP_REGISTER_BANK` | `off` \| `report` \| `transform` | `function(rssn-fp-register-bank)` | 在 Win64 AVX 上使用非易失通用寄存器跨 libm 调用保存标量 `double`。 |
| `RSSN_POW_UNIT_BASE` | `off` \| `report` \| `transform` | `function(rssn-pow-unit-base)` | 当底数在舍入边界内等于 `+1.0` 时跳过 `llvm.pow` 组。 |
| — | MIR（`llc -run-pass`） | `rssn-edge-restore` | 寄存器分配后 x86-64 MIR pass，合并或下沉相同的被调用者保存恢复序列。 |
| `RSSN_MACHINE_STACK_DETAIL` | `0` \| `1` | `rssn-machine-stack-traffic` | 寄存器分配后 MIR 诊断，统计 spill、固定栈槽与常量池栈流量。 |

在任意 IR pass 或 `rssn-edge-restore` 名称后附加 `-report` 可仅输出分析诊断（`-pass-remarks-analysis=<pass>`）。

## 架构概览

- **SSA 活跃度与代价模型**（`Pressure.cpp`、`Planning.h`）：基于 TTI legalization 寄存器占用、ABI 保存寄存器预算与有界二次超限惩罚，沿 GP 与 SIMD 寄存器组反向计算 `LiveBefore = Use ∪ (LiveAfter − Def)`。
- **联合移动与重物化**（`Motion.cpp`）：分桶束搜索（宽度 8，每桶最多 12 个候选，最多联合移动 4 个成员），评估跨调用与局部区域的净活跃区间变化，并支持有界跨块重物化（最多 3 个同层循环块）。
- **读取合并与融合**（`LoadFusion.cpp`、`MemoryModel.cpp`）：消除重复地址读取、将相邻读取融合为 128/256 位向量，并为必须跨修改性调用的标量建立 128 位向量载体；所有变换均受别名分析、`nounwind`/`nofree`/`nosync` 与非 volatile/非原子约束保护。

## 复现实验

```bash
bazel run --repo_env=RSSN_LLVM_CONFIG="$(command -v llvm-config-23)" \
  //:llvm_benchmark -- --cpu 3
```

支持 `--opt`、`--llc`、`--jobs`、`--pairs`、`--work-dir` 与 `--cpu`。脚本在独立目录中交替执行 A/B 与 B/A 构建，进行全数组位级哈希校验并统计 LLVM 的 `Spill`/`Reload` 诊断。

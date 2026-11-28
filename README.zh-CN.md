# RSSN：Rust 科学计算库

[English](README.md) | [简体中文](README.zh-CN.md)

[![Crates.io](https://img.shields.io/crates/v/rssn.svg)](https://crates.io/crates/rssn)
[![Docs.rs](https://docs.rs/rssn/badge.svg)](https://docs.rs/rssn)
[![License](https://img.shields.io/crates/l/rssn)](LICENSE)
[![Discord Server](https://img.shields.io/discord/1459399539403522074.svg?label=Discord&logo=discord&color=blue)](https://discord.gg/D5e2czMTT9)
[![DOI](https://zenodo.org/badge/DOI/10.6084/m9.figshare.31044715.svg)](https://doi.org/10.6084/m9.figshare.31044715)

**rssn** 是面向 Rust 的开源科学计算库：在一个简洁的 API 背后集成计算机代数系统、数值内核与物理仿真。

所有符号运算均视为*恒等变换*。`diff(f, x)` 等请求在哈希驻留的**等价图（e-graph DAG）**上由规则集重写，直至提取出等价闭式或数值证据。详见 [ARCHITECTURE.zh-CN.md](ARCHITECTURE.zh-CN.md)。

---

## 核心特性

- **统一入口**：`Session::compute(term, &config)` 根据目标阶段返回符号闭式或数值结果。
- **量纲与单位**：标量项自动推导物理量纲，支持精确有理数/符号指数、国际单位制与仿射单位，以及建议/严格模式检查（[docs/DIMENSIONS.zh.md](docs/DIMENSIONS.zh.md)）。
- **开放等价图引擎**：哈希驻留 DAG、可扩展算子注册表、分层规则调度、基于代价的提取、假设推理（`positive`、`integer`、`real`）及 `cse` 公共子表达式提取。
- **微积分与微分方程**：求导、不定/定/反常积分（含 Hermite 归约）、极限、级数、Meijer G 函数、Gosper 求和、积分变换、ODE（精确求解、李对称、方程组）、1D/2D/3D PDE 与变分法。
- **代数与分析**：多元因式分解、Gröbner 基（Buchberger/F4/F5）、CAD 量词消去、符号/数值线性代数、数论、组合数学、特殊函数、复分析、泛函分析与统计学。
- **数值内核与信号处理**：自适应/双指数积分、刚性与非刚性 ODE 求解器、求根、稠密/稀疏/Krylov 线性代数、锥规划/MILP 优化、FFT/小波、图算法及任意精度运算（`kernels::mp`）（[docs/MATH_COVERAGE.zh.md](docs/MATH_COVERAGE.zh.md)）。
- **物理与仿真**：经典力学、电磁学、相对论、量子力学、量子场论、凝聚态与热力学规则集，以及 FDM/FEM/FVM/BEM/CFD、分子动力学、SPH、谱方法与多重网格求解器。
- **编译与设备执行**：通过 `Backend` 将闭式编译为 SSA 磁带（支持反向自动微分与区间包络）、Cranelift JIT（`jit`）或 CPU/wgpu 设备（`gpu-wgpu`，`f64` 采用双单精度模拟）（[docs/BACKENDS.zh.md](docs/BACKENDS.zh.md)）。
- **数据 I/O 与可视化**：可往返文本、LaTeX、Typst、Unicode 排版输出、绘图、`.npy`（`npy`）、DataFrame、Parquet、HDF5 与 NetCDF。
- **嵌入与 C ABI**：确定性*内核清单（Kernel Manifest）*、动态 ABI（`Session::call_value`）、宿主回调（[docs/CALLBACKS.zh.md](docs/CALLBACKS.zh.md)）、动态库插件（`plugin-cdylib`）与 C/C++ 头文件（`rssn.h`/`rssn.hpp`）（[docs/KERNEL.zh.md](docs/KERNEL.zh.md)）。

---

## 快速开始（3.0 Beta API）

在现有 Cargo 项目中从 `v0.3` 分支的本地检出添加 RSSN：

```bash
git clone --branch v0.3 https://github.com/Apich-Organization/rssn.git ../rssn
cargo add --path ../rssn
```

保存为 `src/main.rs` 并运行 `cargo run`：

```rust
use rssn::api::*;

fn main() -> Result<(), ComputeError> {
    let s = Session::new();
    let x = s.sym("x");

    // 符号阶段：求取闭式
    let d = s.compute(diff(sin(x) * cos(x), x), &Config::new())?;
    assert_eq!(d.term.to_string(), "cos(2*x)");

    // 数值阶段：在 x = 0.0 处求值
    let n = s.compute(
        diff(sin(x) * cos(x), x),
        &Config::new().numeric(1e-12).bind("x", 0.0),
    )?;
    assert_eq!(n.value, Some(1.0));

    // 也可直接从文本解析请求
    let sum = s.compute(s.parse("sum(1/(k*(k+1)), k, 1, n)")?, &Config::new())?;
    assert_eq!(sum.term.to_string(), "1 - 1/(n + 1)");

    Ok(())
}
```

更多完整示例见[示例索引](examples/README.zh-CN.md)：
- `bazel run //:example -- --profile release quickstart`
- `bazel run //:example -- --profile release simulation`
- `examples/c/smoke.c`（C 客户端）
- `bazel run //:benchmark -- engine`（性能基准）

### Cargo 特性

默认启用 `jit`、`npy` 与 `gpu-wgpu`。使用 `--no-default-features` 构建最小 Rust 库，或运行 `bazel build //:ffi_distribution //:headers` 构建 C/C++ 分发库与头文件（[docs/BUILDING.zh.md](docs/BUILDING.zh.md)）。

| 特性 | 启用内容 |
|---|---|
| `jit` | Cranelift JIT 后端（`rssn::backend::jit`） |
| `npy` | `rssn::io` 中的 NumPy `.npy` 读写 |
| `plugin-cdylib` | 从共享库加载插件（在宿主进程内运行本地代码） |
| `dev-headers` | 使用 cbindgen 重新生成 `rssn.h`/`rssn.hpp`（需 `DEV=1`） |
| `gpu-wgpu` | wgpu GPU 设备（`rssn::backend::device`，支持 `f32` 与双单精度 `f64`；不包含在 `full` 中） |
| `full` | `npy`、`jit` 与 `plugin-cdylib` |

---

## 测试

```bash
bazel test //:rust_tests
bazel test //:rust_tests --test_arg=--test --test_arg=dimensions
bazel test //:rust_tests --test_arg=--test --test_arg=sim --test_arg=acoustics
bazel test //:gpu_test
bazel run //:check -- --lib --tests   # 快速编译检查（约 2 分钟）
```

完整单元测试与 `sim` 套件计算量较大（约 19–25 分钟），建议日常使用 `--test_arg` 过滤运行。CI 会执行 `//:rust_all_features_test`、`//:api_generator_test`、`//:contracts_test` 与 `//:headers_check_test`。详见 [CONTRIBUTING.zh-CN.md](CONTRIBUTING.zh-CN.md)。

---

## 文档导航

| 文档 | 主题 |
|---|---|
| [ARCHITECTURE.zh-CN.md](ARCHITECTURE.zh-CN.md) | 引擎、规则集、内核、仿真、C ABI 与 JIT 架构 |
| [docs/README.zh.md](docs/README.zh.md) | 完整文档索引 |
| [docs/KERNEL.zh.md](docs/KERNEL.zh.md) | 嵌入内核：清单、动态 ABI、插件与单位 |
| [docs/DIMENSIONS.zh.md](docs/DIMENSIONS.zh.md) | 量纲与单位层 |
| [docs/CALLBACKS.zh.md](docs/CALLBACKS.zh.md) · [docs/NUMERICAL_CALLBACKS.zh.md](docs/NUMERICAL_CALLBACKS.zh.md) | 宿主回调与数值回调算子 |
| [docs/BUILDING.zh.md](docs/BUILDING.zh.md) | 构建 Rust 与 C 产物 |
| [docs/BACKENDS.zh.md](docs/BACKENDS.zh.md) | 执行后端：CPU、wgpu、双单精度 `f64` 与启动融合 |
| [docs/MATH_COVERAGE.zh.md](docs/MATH_COVERAGE.zh.md) | 数值约定、覆盖范围与已知边界 |
| [CHANGELOG.md](CHANGELOG.md) | 版本变更日志 |
| [API 参考](https://docs.rs/rssn) | rustdoc 文档 |

---

## 社区与规范

- **参与贡献**：[CONTRIBUTING.zh-CN.md](CONTRIBUTING.zh-CN.md)
- **开源许可**：Apache 2.0（[LICENSE](LICENSE)）· [ATTRIBUTIONS.md](ATTRIBUTIONS.md)
- **安全策略**：[SECURITY.zh-CN.md](SECURITY.zh-CN.md)
- **行为准则**：[CODE_OF_CONDUCT.md](CODE_OF_CONDUCT.md)

# RSSN 示例

[English](README.md) | [简体中文](README.zh-CN.md)

在仓库根目录运行 `bazel run //:example -- NAME` 执行 Rust 示例。

| 示例 | 内容 |
|---|---|
| [quickstart](quickstart.rs) | 会话、符号与数值目标 |
| [symbolic_calculus](symbolic_calculus.rs)、[solve_equations](solve_equations.rs) | 重写、微积分与方程求解 |
| [units_and_dimensions](units_and_dimensions.rs) | 物理单位与量纲检查 |
| [linear_algebra](linear_algebra.rs)、[numeric_kernels](numeric_kernels.rs) | 矩阵分解、数值求解与实数 FFT 重建 |
| [ode_pde](ode_pde.rs)、[simulation_field](simulation_field.rs) | 演化模型与场 |
| [simulation](simulation.rs)、[scenario_catalog](scenario_catalog.rs) | 按名称调度仿真与查看场景类型元数据 |
| [compile_expression](compile_expression.rs) | 可复用的编译求值 |
| [weighted_calibration](weighted_calibration.rs) | 加权分位数与阈值置信度 |
| [engine_timing](engine_timing.rs)、[backend_timing](backend_timing.rs) | 符号负载与 CPU 后端计时 |
| [io_formats](io_formats.rs) | HDF5（分块/过滤/变长/软链接）、Parquet（Brotli/Bloom/PageIndex）、VTK、WAV、PNG、GIF 与 CSV |
| [sim_postprocess_vis](sim_postprocess_vis.rs) | 2D 仿真结果后处理导出为 HDF5、VTK、PNG、SVG、GIF 与 WAV |
| [C 冒烟示例](c/smoke.c) | C ABI 会话句柄与所有权管理 |

计时示例需使用 `tools` 配置，后端对比还需启用 `jit`：

```sh
bazel run //:example -- --profile tools engine_timing
bazel run //:example -- --profile tools backend_timing
```

特性与 C 工具链详见[构建指南](../docs/BUILDING.zh.md)。

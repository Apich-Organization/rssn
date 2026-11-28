# RSSN examples · 示例

Run a Rust example with `bazel run //:example -- NAME` from the repository root.
在仓库根目录运行 `bazel run //:example -- NAME`。

| Example | What it exercises · 内容 |
|---|---|
| [quickstart](quickstart.rs) | Sessions, symbolic and numeric targets · 会话、符号与数值目标 |
| [symbolic_calculus](symbolic_calculus.rs), [solve_equations](solve_equations.rs) | Rewrites, calculus and equations · 重写、微积分和方程 |
| [units_and_dimensions](units_and_dimensions.rs) | Physical units and dimension policies · 单位和量纲检查 |
| [linear_algebra](linear_algebra.rs), [numeric_kernels](numeric_kernels.rs) | Matrix factorizations, solvers and real FFT reconstruction · 矩阵分解、数值求解与实数 FFT 重建 |
| [ode_pde](ode_pde.rs), [simulation_field](simulation_field.rs) | Evolution models and fields · 演化模型与场 |
| [simulation](simulation.rs), [scenario_catalog](scenario_catalog.rs) | Named simulation dispatch and typed scenario metadata · 按名称运行仿真与查看场景类型信息 |
| [compile_expression](compile_expression.rs) | Reusable compiled evaluation · 可复用的编译求值 |
| [weighted_calibration](weighted_calibration.rs) | Weighted quantiles and threshold confidence · 加权分位数和阈值判断 |
| [engine_timing](engine_timing.rs), [backend_timing](backend_timing.rs) | Symbolic workload and CPU backend timing · 符号任务和 CPU 后端计时 |
| [io_formats](io_formats.rs) | HDF5 (chunked/filters/vlen/soft-links), Parquet (Brotli/Bloom/PageIndex), VTK, WAV, PNG, GIF & CSV · 全协议科学数据与媒体 IO |
| [sim_postprocess_vis](sim_postprocess_vis.rs) | 2D simulation post-processing to HDF5, VTK, PNG, SVG, GIF & WAV · 仿真结果后处理与可视化 |
| [C smoke example](c/smoke.c) | Session handles and ownership across the C ABI · C ABI 会话句柄和所有权 |

Use the `tools` profile for timing. The backend comparison also needs `jit`.
计时示例使用 `tools` 配置，后端对比还需要启用 `jit`。

```sh
bazel run //:example -- --profile tools engine_timing
bazel run //:example -- --profile tools backend_timing
```

The [build guide](../docs/BUILDING.md) specifies features and the C toolchain.
特性与 C 工具链见[构建指南](../docs/BUILDING.md)。

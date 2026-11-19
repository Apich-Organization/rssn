# RSSN examples · 示例

Run a Rust example with `cargo run --locked --example NAME` from the repository root.
在仓库根目录运行 `cargo run --locked --example NAME`。

| Example | What it exercises · 内容 |
|---|---|
| [quickstart](quickstart.rs) | Sessions, symbolic and numeric targets · 会话、符号与数值目标 |
| [symbolic_calculus](symbolic_calculus.rs), [solve_equations](solve_equations.rs) | Rewrites, calculus and equations · 重写、微积分和方程 |
| [units_and_dimensions](units_and_dimensions.rs) | Physical units and dimension policies · 单位和量纲检查 |
| [linear_algebra](linear_algebra.rs), [numeric_kernels](numeric_kernels.rs) | Matrix factorizations and numerical algorithms · 矩阵分解与数值算法 |
| [ode_pde](ode_pde.rs), [simulation_field](simulation_field.rs) | Evolution models and fields · 演化模型与场 |
| [compile_expression](compile_expression.rs) | Reusable compiled evaluation · 可复用的编译求值 |
| [weighted_calibration](weighted_calibration.rs) | Weighted quantiles and threshold confidence · 加权分位数和阈值判断 |
| [C smoke example](c/smoke.c) | Session handles and ownership across the C ABI · C ABI 会话句柄和所有权 |

The [build guide](../docs/BUILDING.md) specifies features and the C toolchain.
特性与 C 工具链见[构建指南](../docs/BUILDING.md)。

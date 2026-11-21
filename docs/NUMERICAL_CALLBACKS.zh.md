# 数值宿主回调算子

[English](NUMERICAL_CALLBACKS.md) · [宿主回调契约](CALLBACKS.zh.md)

`numerical_callback` 规则集通过内核清单公开有宿主副作用的串行数值运算。末尾的函数实参是注册在该次 `CallCtx` 中的 `CallbackRef`。回调库在注册时保存其具体签名。按名称调用和按 `OpKey` 调用使用同一求值器。

| 算子 | 普通实参 | 回调 | 结果字段 |
|---|---|---|---|
| `num_integrate` | `a`、`b`、绝对容差、相对容差、子区间预算 | `Float[D] -> Float[E]` | `value: Float[D*E]`、`error: Float[D*E]`、`evaluations`、`converged` |
| `num_brent` | `a`、`b`、正的自变量容差、正的迭代预算 | `Float[D] -> Float[E]` | `root: Float[D]`、`evaluations` |
| `num_brent_with_bracket` | 同一夹逼区间实参 | 同一标量残差 | `root: Float[D]`、`evaluations`、`bracket: Tensor<[2], Float[D]>` |
| `num_least_squares` | 初始参数向量、残差容差、梯度容差、步长容差、初始阻尼、迭代预算、求值预算 | `Tensor<1, Float[D]> -> Tensor<1, Float[E]>` | `x: Tensor<1, Float[D]>`、`residual_norm: Float[E]`、`iterations`、`evaluations`、`jacobian_evaluations`、`accepted_steps`、`rejected_steps`、`converged`、`status` |
| `num_ode_dopri5` | 初始时间、终止时间、状态向量、状态绝对容差、相对容差、最大时间步长、尝试预算 | `(Float[D], Tensor<1, Float[E]>) -> Tensor<1, Float[E/D]>` | `times: Tensor<1, Float[D]>`、`states: Tensor<2, Float[E]>`、`evaluations`、`rejected` |
| `num_ode_rosenbrock23`、`num_ode_radau5`、`num_ode_bdf` | 同一初值问题实参 | 同一向量右端项 | 同一初值问题结果 |

表中数值系数要求有限浮点数。线性物理量在进入求解器前转为相干数值。仿射点会被拒绝。定义域、状态和导数的量纲按实际回调签名检查。绝对容差恰为零时可匹配任意量纲。两种容差均非负且不能同时为零。ODE 最大步长必须为正且有限，起始步长在指定区间内自动选择。

`num_least_squares` 通过一个作用域内残差回调实现 Levenberg–Marquardt 型阻尼非线性最小二乘。参数向量与残差向量的量纲 `D` 和 `E` 独立。残差容差具有 `E`，梯度容差、步长容差和初始阻尼无量纲。参数向量非空。该算子通过有限差分计算雅可比。使用解析雅可比的 Rust 接口见 `kernels::rootfind::least_squares`。`status` 为 `residual_converged`、`gradient_converged`、`iteration_limit`、`evaluation_limit` 或 `stagnated`。前两种状态的 `converged` 为真。

ODE 状态向量必须非空。每次导数的宽度与状态一致。调取宿主代码前，求解器会检查签名给出的固定形状。每次返回的张量还要检查秩、形状与存储的一致性、宽度和物理量纲。结果状态矩阵每行对应一个接受的时间点。四个方法均支持反向时间。需要事件、质量矩阵、微分代数方程或稠密输出时，可使用专用的原生算子与 Rust 求解器接口。

首个宿主错误或回调结果转换错误原样传递，包括 `OpError::HostDiagnostic` 的错误码。之后求解器停止调用宿主。ODE 适配器借助内部取消检查停止自适应重试。`evaluations` 记录实际调用次数。Brent 不会为报告额外求一次根的残差，其预算统计两端采样后的迭代。ODE 的预算统计所有接受及拒绝的步长尝试。预算耗尽产生显式 `NotReduced`。数值积分保留自己的 `converged` 字段。

数值积分与 ODE 的估计均为启发式浮点误差估计。Brent 需要两端异号或某端点恰为根，自变量容差使用定义域单位。
`num_brent_with_bracket` 按升序返回最终区间的端点。找到精确零点时，两端点均为该坐标。

`num_native_hessian_diagonal_scaled` 接收带量纲的坐标点、正的坐标与输出参考量、求值预算和标量回调。`num_native_hessian_diagonal_with_steps_scaled` 在回调之前接收每个坐标方向的正步长。回调输入量纲为 `D`、输出量纲为 `E` 时，对角项的量纲为 `E/D²`。拉普拉斯使用同一组采样对并对相关对角项求和。

## 共用契约

`numerical_callbacks::common::CallbackRun` 统一保留首个错误并计数。注册表为这些算子声明有状态、副作用、会话约束与解释器执行能力。宿主句柄只在注册其的调用上下文内有效。

`tools/rssn-api-gen` 的 structured 模式根据公开 Rust 结构体生成记录编解码器。生成代码使用穷尽字段解构，字段变更会在审查策略前触发编译错误。积分的 `Integral` 与 SVD 记录使用同一字段策略生成清单类型和运行时量纲。编解码器覆盖与可调用算子覆盖分别报告。求解器的量纲契约由其签名定义。

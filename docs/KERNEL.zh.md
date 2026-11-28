# RSSN 内核嵌入接口

[English](KERNEL.md) · [量纲与单位](DIMENSIONS.zh.md)

RSSN 为宿主程序提供两个嵌入模块。它们共用表达式图、算子表和量纲域。

| 模块 | 接口 |
|---|---|
| `rssn::kernel` | 清单、动态值、签名、插件、编码和结构化结果 |
| `rssn::dim` | 量纲表达式、量纲域、单位注册表和物理常数 |

## 内核契约

### 1. 编码

- `encode_binary` 和 `decode_binary` 使用 bincode-next 3.2 原生确定性编码。它用于内核二进制存储与清单哈希。解码需要完整且仅含一个值的缓冲区。
- 使用 Serde 类型的宿主可以调用 `encode_deterministic` 和 `decode_deterministic`。对支持的内核记录，原生与 Serde 二进制编码遵循相同线格式。
- `encode_cbor` 和 `decode_cbor` 使用确定性的自描述 CBOR，用于跨系统交换。`digest` 和 `short_digest` 对确定性编码计算 BLAKE3。
- `KernelManifest::to_json` 与 C ABI 的 JSON 接口提供可读的元数据和结果。

构建会嵌入版本控制中的 `signatures.overlay.json` 和 `signatures.witnessed.json`。需要声明更精确的契约时改写 overlay。见证表由实际运算重放生成。缺少其中任一文件都会导致构建失败。[构建指南](BUILDING.zh.md)

### 2. 算子标识

`OpDescriptor` 是算子的注册描述。它产生 `Signature`、`OpKey` 和清单中的 `OpInfo`。`OpKey` 包含 `(rule_set, name, signature_hash)`，对于相同签名不随注册顺序改变。一个算子的每个重载有自己的键。

```rust
let manifest = session.manifest();
let sin = manifest.op("sin").expect("sin is a core operator");
assert!(!sin.overloads.is_empty());
let key = &sin.keys[0];
```

`OpInfo` 的每个字段对每个算子都有值，宿主生成接口时不必为缺口单独处理。注册未声明时有两类字段会取默认值：

* `category` —— 宿主分组用的族名。注册可用 `.category(..)` 指定更细的类别，否则取规则集名，因此不会为空。
* `doc`、`doc_source` 和 `doc_path` —— 说明文本、来源，以及提供 rustdoc 的函数（见 3b）。
* `since` —— 算子出现的版本。注册未记录时为空；宿主应把空值读作"未记录"而非"未发布"，算子是否稳定由 `stability` 表示。

`tests/manifest_completeness.rs` 固化这些不变量，任一算子违反即失败。

### 3. 签名推导

注册信息足以推导可调用算子的基础签名。缺少专门类型的非叶算子有符号项签名。声明数值求值器的算子还可有数值签名。例如 `sin(0.5)` 可以折叠为数值，`sin(x)` 保持符号项。

| 注册信息 | 推导结果 |
|---|---|
| 标量求值器 `eval` | 实数参数和实数结果 |
| `OpFlags::PREDICATE` | `Bool` 结果 |
| 第 `k` 个参数带 `Binder` | 第 `k` 个参数是绑定 `Symbol` |
| `Arity::Variadic` | 尾部参数可重复 |
| `OpFlags::LEAF` | 直接记录 `dynamic_reason` |

`SignatureSource` 区分直接声明 `Declared`、人工覆盖 `Overlay`、推导 `Derived`、见证 `Witnessed` 和缺少静态签名 `None`。

overlay 只作为未声明签名的算子的类型来源：声明优先，因此给已声明签名的算子再写 overlay 条目没有作用，还会在声明变化时与代码悄悄分叉。文件因此恰好包含清单报告为 `overlay` 来源的算子；两侧任一方向漂移（登记了不存在的算子，或清单从不读取的条目）都会让 `tests/signature_overlay.rs` 失败。

### 3a. 见证签名

`factorint` 或 `groebner` 这类算子如果只声明参数个数，就会得到宽泛的 `(Term, ..) -> Term`。当 `factorint(360)` 与 `list(list(2, 3), list(3, 2), list(5, 1))` 进入同一等价类时，计算结果可以为更具体的签名提供证据。

1. `kernel::sort` 将等价类归入由 `Int`、`Rat`、`Float`、`Sym`、`Expr`、`Bool`、`Prop`、具名对象、列表和元组等组成的格。推导依据是字面量、列表构造器、谓词标记和量纲规则。
2. 设置 `RSSN_SORT_TRACE` 后，每次顶层运算结束时记录算子、实参、结果类别和可重放调用。
3. `signature-mine` 在完整标准内核上重放去重后的调用。它保留有回答的调用，按实参类别分组，再合并同组结果。仍留在原等价类中的未求值请求不构成答案证据。
4. 候选签名会再用扰动后的数值实参检查，保留原形状。比如 `sin(0) = 0` 单独无法证明整数签名。证据稀少且含符号实参时，结果中的精确数值也会拓宽为表达式。
5. 只有实参签名与其他重载可比较或不相交的候选才保留。同一实参签名下更精确的结果单独记录为 `refinements`，不参与重载选择。

```text
RSSN_SORT_TRACE=sorts.tsv cargo test --lib
cargo run --profile tools --bin signature-mine -- sorts.tsv
cargo run --profile tools --bin signature-mine -- --check
```

生成的 `signatures.witnessed.json` 随源码保存。清单用 `signature_source` 和 `overload_sources` 记录来源，用 `witnesses` 提供可重放调用。见证不证明对所有实参都能求出同类结果，宿主还要验证解码结果。两个可比较的重载同时匹配时，逐参数比较其具体程度。

### 4. 动态算子的原因

| `DynamicReason` | 含义 |
|---|---|
| `request` | 由规则求值，使用实际触发规则产生的结果类型 |
| `variadic` | 参数个数可变 |
| `fully_polymorphic` | 接受任意类型 |
| `shape_dependent` | 容器形状由实参决定 |
| `opaque_payload` | 携带内核不解释的宿主数据 |
| `not_yet_specified` | 缺失声明，覆盖检查报错 |

`cargo run --locked --bin signature-audit` 按签名来源汇总当前算子。若算子既无签名也无获认可的原因，命令返回非零退出码。

### 5. 确定性

相同的算子表、规则集合和构建能力生成逐字节一致的清单。清单对各列表排序去重，不包含进程内编号、指针或时间戳。计算 `manifest_hash` 时，将该字段置零后以确定性二进制编码计算 BLAKE3。宿主可用 `(manifest_hash, OpKey, target)` 作为缓存键。

## 预计算模算术

`num_modular_montgomery_new(modulus)` 接受大于 1 的奇数无符号 64 位模数，
`num_modular_barrett_new(modulus)` 也支持偶数模数。`_new_prime` 构造函数会检查素性。
运算使用精确整数，支持输入范围内的 64 位与 128 位数据。求逆或除法在分母与模数
不互素时返回 `None`，参数超出范围、句柄类型不符或已释放时返回错误。

Montgomery 编码运算要求输入小于对应的模数。应使用同一个计划转换、运算和解码，
原始 REDC 接受小于 `modulus * 2^64` 的无符号整数。Barrett 支持整个无符号 128 位范围。
计划是会话所属的不可变对象，可重复使用。完成后调用 `num_native_release` 释放。

## 动态调用与重载选择

```rust
use rssn::kernel::{CallCtx, Value};

let session = rssn::api::Session::new();
let args = [Value::Float(0.5)];
let key = session.resolve_overload("sin", &args)?;
let value = session.call_key(&key, &args, &CallCtx::new())?;
```

`call_value` 按名称、参数个数和类型选择重载，也可返回仍需符号求值的项句柄。`call_key` 直接选定重载，适合已生成静态绑定的宿主。`Session::overloads(name)` 列出签名及键。无签名算子通过名称调用，签名哈希为全零。

匹配优先级是 `Exact`、`Widened`、`Shaped`、`Generic`，调用按最弱的实参匹配评分。等价候选产生 `OpError::Ambiguous`。参数个数或找不到算子分别有 `OpError::Arity` 与 `OpError::UnknownOp`。`CallCtx` 携带符号或数值目标、搜索预算、量纲模式、假设、变量绑定、代入、宿主回调和可选跟踪。算子可向跟踪追加事件。

### 拟合和复用多项式

便携多项式模型用降幂系数表示归一化坐标 `x / coordinate_scale` 的多项式。
拟合一次后，可以把模型传回内核，一次求出多个坐标的值。结果保持查询顺序。

```python
import rssn

with rssn.Session() as session:
    model = session.call_value(
        "num_polynomial_model_interpolate_lagrange", [0.0, 2.0, 4.0], [1.0, 5.0, 17.0], 2.0
    )
    values = session.call_value("num_polynomial_model_eval_batch", model, [3.0, 1.0])
    print([v.float() for v in values])  # [10.0, 2.0]
```

对于带单位的数据，正坐标标度和查询坐标与拟合横轴同量纲，预测值保留纵坐标的量纲。
`num_polynomial_model_from_roots` 从根构造的多项式在归一化坐标下首一。
调用 `num_polynomial_model_compose` 时先传外层模型，其坐标量纲须与内层模型的输出量纲相同。

### 批量计算特殊函数

向 `num_special_gamma_batch` 传入一维 Float 张量，结果保持输入顺序。
批量数据达到内核的并行阈值后会按执行策略分配工作。对数伽马、双伽马、
误差函数、余误差函数及贝塞尔函数的批量接口沿用同一输入和结果形式。

```python
import rssn

with rssn.Session() as session:
    result = session.call_value("num_special_gamma_batch", [1.0, 2.0, 3.0, 4.0])
    print([item.float() for item in result])  # [1.0, 1.0, 2.0, 6.0]
```

输入和输出均无量纲。如果标量内核产生非有限值，例如实数伽马函数在极点的结果，
值接口会报告定义域错误。

### 双曲坐标与魏尔斯特拉斯周期

用 `num_uhp_to_disk` 和 `num_disk_to_uhp` 在开上半平面与庞加莱开圆盘之间
转换复坐标。`num_uhp_distance` 与 `num_poincare_disk_distance` 分别在两种
模型中计算曲率为 −1 的测地距离。例如，上半平面中的 `i` 和 `4i` 相距
`ln(4)`。参数与结果均无量纲。边界点以及舍入后落在目标开区域之外的转换
结果会报定义域错误。

`num_weierstrass_p_prime_complex(z, omega1, omega2)` 接收两个复数半周期，
对应的基本周期为 `2 * omega1` 和 `2 * omega2`。
`num_weierstrass_invariants_complex(omega1, omega2)` 返回带有 `g2` 与
`g3` 字段的记录。半周期须有限且不共线。导数遇到格点极点时报定义域错误。
计算采用 binary64，底层 theta 函数求值若产生非有限结果也会报定义域错误。

### 向量更新和批量扫描

`num_vector_axpy(alpha, x, y)` 计算 `alpha*x + y`，增量必须与原向量同量纲。
`num_vector_axpby` 对两个向量分别缩放，要求缩放结果的量纲相同。
两个接口都返回新向量。

```python
import rssn

with rssn.Session() as session:
    updated = session.call_value("num_vector_axpy", 2.0, [1.0, 3.0], [5.0, 7.0])
    print([x.float() for x in updated])  # [7.0, 13.0]
```

逐元素积和外积将输入量纲相乘，`num_vector_hadamard_div` 将输入量纲相除。
固定分块的二范数和包含当前项的前缀和保留输入量纲。
`num_par_dot_fixed` 对输入向量的共同前缀求点积，结果为两者的量纲乘积。
分块大小决定浮点数的累加分组，`grain` 是启动并行计算的阈值。
比较不同线程数的结果时应使用相同的分块大小。

### 精确整数和模运算

`num_finite_field_` 系列的向量算子使用 `Tensor<1, Int>`，接受完整的
`u64` 范围和大于等于 2 的模数，也支持合数模。加减、AXPY 和点积要求输入
向量等长。批量求逆的严格版本在存在不可逆元素时报告定义域错误，
`num_finite_field_batch_inverse_skip_zeros` 则将不可逆位置置零。
空向量可参与向量运算，空向量的点积为 0。

`num_number_theory_gcd_u128`、`num_number_theory_mul_mod_u128` 和
`num_number_theory_mod_pow_u128` 接受完整的无符号 128 位范围。
模运算要求正模数。所有参数和结果均为无量纲精确整数。
绑定的宽整数转换方法见[语言绑定契约](../bindings/api/CONTRACT.md)。

### 通过 C ABI 传入带量纲浮点值

用 `rssn_value_with_coherent_dimension` 将具体量纲附到相干单位下的 Float
或稠密 binary64 张量。对于米，传入 UTF-8 JSON `{"base":"Length"}`
及其字节数。结果在同一会话中使用独立的值单元，保留张量的阶、空维长度及 Float 表示。
要求带量纲参数的算子也须接收带量纲的零值。

JSON 沿用清单的量纲格式。为带单位名称的值添加量纲之前，先换算为相干单位下的值。
可用 `rssn_value_json` 查看返回值的量纲。

### 带量纲三对角方程组

对称正定三对角方程组使用 `num_dense_solve_spd_tridiagonal(diag, off, rhs)`。
消元主元均非零的一般三对角方程组可用
`num_dense_solve_tridiagonal(sub, diag, sup, rhs)`。两者的时间和空间开销均
随方程组大小线性增长。`diag` 和 `rhs` 的长度为 `n`，次对角线长度为
`max(n - 1, 0)`。各条对角线的量纲一致，解的量纲为右端量纲除以矩阵量纲。
一般求解器采用不选主元的 Thomas 算法。对于需要选主元的非奇异方程组，
可将完整矩阵传入 `num_dense_solve`。

### 测量网格

`num_meshgen_delaunay_2d_incremental(points)` 接收 3 到 4096 个互不相同的
二维点，返回逆时针三角形的零基顶点索引。坐标张量中各项须使用同一量纲。
`num_meshgen_triangle_areas(points, triangles)` 和
`num_meshgen_tetra_volumes(points, tetrahedra)` 按单元输入顺序返回面积和体积，
量纲分别为坐标量纲的平方和立方。越界顶点索引会报告定义域错误，退化单元的测量结果为零。

`num_meshgen_min_angle_deg` 以度数返回最小角。`num_meshgen_boundary_faces` 提取只与
一个四面体相邻的三角面并保持顶点绕序。正定向四面体的边界面朝外，内部共用面被移除。

### 精确组合计数与浮点数值入口

`num_harmonic(0)` 返回 0。`num_ln_binomial(n, k)` 在
`0 <= k <= n` 时直接计算对数，即使系数本身超出 binary64 范围也可用。
`num_multinomial(n, parts)` 要求各分量之和精确等于 `n`，并在结果溢出时报告定义域错误。
需要任意精度整数结果时，使用符号算子 `binomial` 和 `multinomial`。

`num_inverse_beta_reg(a, b, p)` 接收有限正形状参数和区间 `[0, 1]` 内的概率。
端点及形状参数等于 1 的情形用解析公式。一般情形在 binary64 上有界求解，
精度取决于正向概率计算。正向计算超出可表示范围时返回错误。

### 零相位滤波

`num_signal_filtfilt(b, a, samples)` 采用奇反射填充和稳态初始化，沿正反两个方向滤波。
系数向量 `b`、`a` 无量纲，返回值保留输入样本数和量纲。每侧填充长度取滤波器阶数的
三倍与样本数减一中的较小值。单样本也按两次滤波的稳态响应处理。系数向量为空、
分母首项为零、稳态无定义或存在非有限值时，返回定义域错误。

### 多项式根界与单位圆盘

`num_cauchy_root_bound`、`num_cauchy_sharp_root_bound` 和
`num_fujiwara_root_bound` 接收按幂次升序排列的有限复系数。
`num_enestrom_kakeya_bounds` 接收至少两个严格为正的实系数，并返回包含
`lower`、`upper` 半径的记录。多项式自变量无量纲，所有系数可以共用同一量纲，
该量纲在根界计算中抵消。通过值接口调用时，超出 binary64 范围的根界会返回定义域错误。

对于单位圆盘上的归一化单叶函数，`num_koebe_growth_bounds(r)` 和
`num_koebe_distortion_bounds(r)` 分别返回函数模和导数模的上下界。
`num_koebe_rotation_bound(r)` 返回导数连续辐角绝对值的上界，单位为弧度。
输入半径无量纲且满足 `0 <= r < 1`。

### 原根与精确多项式乘法

`num_number_theory_extended_gcd_u64(a, b)` 返回精确的 `{gcd, x, y}`
记录，其中 Bézout 系数带符号。`num_number_theory_mod_inv_u64` 接收至少为 2 的模数，
在输入与模数互素时返回最小非负逆元。Jacobi 与 Legendre 接口区分非法分母和合法的零值。

`num_number_theory_primitive_root` 返回循环剩余类乘法群的最小正生成元。
`num_finite_field_primitive_root_prime` 要求模数为素数，也接受 2。
`num_finite_field_poly_mul_ntt_mod(a, b, p)` 返回奇素数 `p` 下的系数卷积。
两个输入的系数排列方向须一致。非空结果的长度向上取到最近的二次幂后，必须整除
`p - 1`。模数有效时，空输入返回空张量。值接口使用通用的原生张量长度上限。

`num_finite_field_ntt_plan(modulus, length)` 验证奇素数和整除
`modulus - 1` 的二次幂长度，并返回仅含 `modulus` 和 `length` 的记录。
`num_finite_field_ntt_transform(values, plan, inverse)` 接受等长的普通
无符号 64 位整数张量。正变换使用由最小正生成元得到的根，逆变换在域中
除以长度。恢复的记录会在调用时重新验证，长度受原生值预算限制。

### 生存曲线与 Cox 回归

`num_survival_kaplan_meier(time, event, confidence)` 返回事件时间、风险
人数、事件次数、生存概率、Greenwood 方差和 log-log 置信界。若时间
带量纲，`num_survival_median(curve)` 保留该量纲。曲线未达到一半时
返回空选项。

`num_survival_fit_cox(time, event, covariates, ties, max_iterations, tolerance)`
接受无量纲、无截距的协变量。可用 `"breslow"` 或 `"efron"` 处理并列
事件；删失观测仍在对应时刻的风险集中。将拟合记录传给
`num_survival_predict_cox(fit, covariates)` 可获取相对危险度。预测时
必须与原拟合保持列数和列顺序一致。

### 矩形 CSR 运算

`num_csr_transpose` 保留条目量纲。`num_csr_matmul_sparse` 和
`num_csr_matmul_dense` 的结果条目量纲是两个输入的乘积。CSR 记录包含
`rows`、`cols`、`indptr`、`indices` 和 `values`，列下标从零开始并
在每行严格递增。

`num_csr_lsqr(matrix, rhs, damping, tolerance, max_iterations)`
求解矩形最小二乘问题。阻尼系数与矩阵条目同量纲，解的量纲为右端项除以
矩阵条目量纲。结果包括解、原系统相对残差、迭代次数和收敛状态。最优
最小二乘解可以具有非零残差。

### 稀疏序列二次规划

`num_sqp_solve` 使用带类型的求值和 Hessian 回调求解归一化非线性规划。
求值记录包含目标值、梯度、约束值和稀疏 Jacobian。Hessian 回调接受迭代
点和两组乘子，返回 Lagrangian Hessian 的上三角。稀疏记录由零起点的
`rows`、`columns` 和 `values` 张量组成。约束先列等式，再列形式为
`g(x) <= 0` 的不等式。

回调使用归一化变量、目标和约束。同一个连续访问的点只求值一次，两个
回调共享预算。结果包括迭代点、乘子、KKT 残差、收敛状态和调用次数。
还可用 `num_sqp_oracle_*` 单独检查各回调分量。
回调接口最多接受 128 个坐标、128 个约束以及每个稀疏记录 16,384 个三元组。

### 精确十进制交换

`num_mp_decimal_from_text(text, precision)` 将十进制字面量直接舍入为二进制值，
等距时取偶。`from_text_in` 可指定舍入模式。`num_mp_decimal_to_text(value, digits)`
按有效位数格式化精确二进制值，`to_text_in` 可指定格式化方向。有效位数范围为
1 至 10,000。十进制解析和包络接口接受至多 100,000 字节的文本，去掉尾零前
的十进制尺度须在 ±100,000 内。结果超过二进制指数范围时返回定义域错误。

`num_mp_decimal_enclose` 返回包含精确输入的二进制区间。
`num_mp_decimal_format_ball` 只在两个端点舍入结果相同时返回文本。
`num_mp_decimal_round_text` 直接舍入十进制文本，接受 ±40,000 内的十进制尺度。

### 精确二进制多精度值

`num_mp_real_round(mantissa, exponent, precision, mode)` 返回无量纲记录
`{mantissa, exponent, precision}`。三个字段都是精确整数，数值为
`mantissa * 2^exponent`。非零尾数的绝对值恰占 `precision` 位。
零值的指数为零，并保留请求精度。

精度范围为 1 至 32,768 位。构造器最多接收 65,536 位的输入尾数。
输入和规格化后的指数绝对值均不得超过 262,144。舍入模式为
`nearest`（等距取偶）、`zero`、`away`、`up` 和 `down`。

用 `num_mp_real_exp(value, precision)` 计算精确二进制输入的指数函数，
或用 `num_mp_ball_exp(interval, precision)` 计算闭二进制区间上的指数函数。
两者均返回无量纲二进制区间。中点采用请求精度，半径包含舍入和级数误差。
指数超出支持范围时返回定义域错误。

用 `num_mp_real_zero` 和 `num_mp_real_one` 创建指定精度的常量。
`num_mp_real_from_i64` 对有符号 64 位整数按等距取偶舍入。
`num_mp_real_from_bigint` 支持更宽的整数及显式舍入模式。
`num_mp_real_from_f64` 把有限 binary64 值舍入到请求精度，两个有符号零均变为
同一规格化零值。`num_mp_real_to_f64` 直接按等距取偶舍入至 binary64，保留次正规数
和下溢零的符号，溢出时返回定义域错误。`num_mp_real_sign` 从精确尾数返回 -1、0 或 1。

`num_mp_real_compare(left, right)` 可精确比较不同精度的值，返回 -1、0 或 1。
`num_mp_real_add`、`num_mp_real_sub`、`num_mp_real_mul` 和 `num_mp_real_div` 接收
`(left, right, precision, mode)`，对精确结果只舍入一次。除法拒绝零除数。
通过值接口再次传入结果时，保持记录字段为整数。

`num_mp_ball_add`、`num_mp_ball_sub`、`num_mp_ball_mul` 和 `num_mp_ball_div` 均接收
`(left, right, precision)`，依次返回两闭区间的和、差、积及商的包络。
`{midpoint, radius}` 的两个字段都是上述精确二进制记录，`radius` 非负。
该记录表示 `[midpoint - radius, midpoint + radius]`。请求精度决定中点的精度，
半径可以采用不同精度以保持包络。除数区间不得含零，包括端点为零的情况。
指数超出可表示范围时返回定义域错误。

复数矩形是无量纲的 `{real, imag}` 记录，两个字段均为闭二进制区间。
可用 `num_mp_complex_new` 组合两个区间，用 `num_mp_complex_point` 创建精确坐标，
或用 `num_mp_complex_from_real` 将实区间嵌入虚部为零的矩形。
加、减、乘、除、倒数、模平方和有符号整数幂沿用实区间算术的精度范围。
模平方返回二进制区间。除数包络含原点时，除法、倒数和负整数幂返回定义域错误。

### 时序统计、重建与状态空间模型

`num_timeseries_integrate(differences, initial)` 还原多阶一阶差分，
阶数由 `initial` 的长度决定。两组输入须共用量纲，输出保留该量纲。

`num_timeseries_rolling_moments(samples, window, ddof)` 按窗口起点顺序返回完整窗口的
`mean` 和 `variance` 向量。均值保留样本量纲，方差使用其平方。总体方差使用
`ddof=0`。样本方差使用 `ddof=1`，窗口宽度须至少为 2。正的窗口宽度大于序列
长度时返回两个空向量。样本须为有限数，且 `ddof` 须小于窗口宽度。

`num_timeseries_kalman_filter(y, f, h, q, r, x0, p0, observation_scale)`
返回滤波均值、预测均值、按时间排列的三阶滤波协方差张量及对数似然。
第一次预测使用给定的先验。状态分量共用一个量纲，观测分量共用另一个量纲，
观测矩阵负责在两者之间转换。正的 `observation_scale` 使用观测单位，
作为无量纲对数似然的参考量。协方差使用相应状态或观测量纲的平方。

`num_timeseries_kalman_smooth(y, f, h, q, r, x0, p0)` 按时间顺序返回每条观测的
平滑状态均值。反向求解每一步都需要正定的预测协方差。

### 精确二进制值与区间的定向求根

`num_mp_real_cbrt_in(value, precision, rounding)` 按指定方向舍入实立方根，
支持负输入。`num_mp_real_cbrt` 和 `num_mp_real_hypot` 分别返回实立方根
和两个精确二进制值的欧几里得范数的包络。
`num_mp_ball_sqrt` 返回非负闭区间上平方根的包络。

`num_mp_real_limits()` 返回支持的位精度和二进制指数范围。二进制值和区间
的读取算子保留精确整数字段。`num_mp_real_rounding_mode(rounding)` 返回
受支持舍入模式的标准拼写。

### 带量纲状态空间仿真

`num_control_c2d_zoh(a, b, c, d, dt)` 在零阶保持下离散化连续模型，
返回可移植记录 `{a, b, c, d, dt}`。采样周期 `dt` 须为正且带时间量纲。
对于 `n` 个状态、`m` 个输入和 `p` 个输出，四个矩阵的形状依次为
`[n,n]`、`[n,m]`、`[p,n]` 和 `[p,m]`。连续模型的 `a` 使用
逆时间量纲。离散模型的 `a` 无量纲，`b` 增加一个时间因子。

`num_control_simulate_discrete(model, x0, inputs)` 可直接使用离散化记录。
`num_control_simulate_continuous(a, b, c, d, x0, inputs, times)`
还要求一个严格递增、带时间量纲的网格，每个时刻对应一行输入。两个算子均按输入行
返回输出，首行使用给定初始状态，并计入直接馈通项 `d`。连续仿真将每行输入保持到下一采样时刻。

## 插件

插件通过 `Registry` 注册算子与重写规则。`SimplePlugin` 接收一段注册函数，向 `Registry::register_op` 提供 `OpDescriptor`，并可用 `register_rewrites` 设置规则和重写层级。每个注册算子需要签名或 `DynamicReason`。同名但形状不兼容的注册失败。为同一算子增加签名时使用描述符的 `.overload()`。[Rust 插件示例](KERNEL.md#plugins)

## 复数矩形函数的可靠包络

`num_mp_complex_abs(value, precision)` 返回包住复数矩形模长的闭二进制区间。
当虚部中点位于 ±2^20 弧度以内时，`num_mp_complex_exp` 返回包住复指数的矩形。
`num_mp_complex_ln` 返回复对数主值的矩形包络。输入接触非正实轴时返回定义域错误。
计算超出指数或精度上限时返回错误。

`num_mp_complex_sin`、`num_mp_complex_cos`、`num_mp_complex_tan` 及双曲变体
`num_mp_complex_sinh`、`num_mp_complex_cosh`、`num_mp_complex_tanh` 同样接受并
返回二进制矩形。三角函数对实部中点约化角度，双曲函数对虚部中点约化角度，范围为
±2^20 弧度。计算得到的分母矩形包含零时，正切返回定义域错误。纯虚输入的正切与纯实
输入的双曲正切调用实数双曲正切，可以处理较大的有限参数。

`num_mp_complex_sqrt` 包住主平方根。负实轴上的主平方根具有正虚部。
横跨支割线的输入矩形可能返回含有两种虚部符号的包络。
`num_mp_complex_pow` 接受底数和指数矩形。指数为精确实整数单点时使用有符号 64 位
整数幂，包括零次幂。其余情形通过可靠复对数和复指数计算主值，底数必须避开非正实轴，
并满足复指数的角度范围。

`num_mp_complex_{asin,acos,atan,asinh,acosh,atanh}` 返回无量纲二进制矩形上
复反函数主值的包络。对 [-1, 1] 内的实区间，反正弦和反余弦包含两个端点。
反双曲正弦包含对应的虚轴端点。反双曲余弦包含点 1，输入矩形与实轴上小于 1
的部分相交时返回定义域错误。反正切和反双曲正切遇到各自的虚轴和实轴支割线或
分支点时返回同类错误。计算中的包络触及奇点或超出指数范围时也返回定义域错误。

需要取得算术表达式的精确二进制中点和向外舍入半径时，可调用
`num_mp_expression_enclose(text, precision)`。精度参数以位为单位。
例如，`sqrt(81) + hypot(3, 4) + cbrt(-8)` 的准确结果是 12。
语法错误、无法认证的定义域和超出范围的指数会返回定义域错误。

`num_mp_ball_lipschitz(interval, mapped_midpoint, bound, precision)`
用调用方提供的非负 Lipschitz 界传播输入半径，并纳入映射中点的舍入误差。
调用方应提供输入中点的精确像，并确保给出的界在整个区间成立。

带类型的 `num_mp_complex_erf` 和 `num_mp_complex_erfc` 返回无量纲二进制
闭矩形的像的包络。泰勒级数对整个矩形使用统一的几何余项界。在
`Re(z) >= |Im(z)|` 的右侧扇形中，余误差函数可用渐近展开的首个遗漏项
控制误差，从而保留正实轴上的小尾值。输入过宽或级数预算不足时返回定义域错误。
参见 [DLMF §7.6](https://dlmf.nist.gov/7.6) 和 [DLMF §7.12](https://dlmf.nist.gov/7.12)。

`num_mp_complex_bessel_j(value, order, precision)` 返回复 Bessel J 主分支
在两个二进制矩形操作数上的包络。绝对值不超过 512 的精确整数阶允许负实
轴参数。一般阶的整个矩形必须满足实部非负，精度上限为 512 位，且参数
矩形必须避开非正实轴。在原点，零阶返回一，非零精确整数阶或实部严格为
正的阶返回零。其余参数的整个矩形必须满足
`abs(Re(value)) + abs(Im(value)) <= 256`。可根据输出半径判断输入不确定
度带来的影响。

复 gamma 系列与 Bessel J 共用带严格余项界的 Stirling 展开。
`num_mp_complex_gamma` 和 `num_mp_complex_digamma` 支持避开所有非正整数
极点的二进制矩形。`num_mp_complex_log_gamma` 沿正实轴取实值并解析延拓，
因此不接受包含非正实轴的矩形。三种接口的输入精度上限均为 1024 位，
递推预算耗尽时返回错误。可直接计算的正整数阶 gamma 先求精确阶乘，再
向外舍入。

`num_mp_complex_riemann_zeta(order, precision)` 返回 zeta 的外向包络，
阶矩形必须避开一点的极点。Euler–Maclaurin 路径要求实部至少为 -32，
模至多为 256。具有闭式的整数阶可直接求值。
`num_mp_complex_polylog(value, order, precision)` 要求完整参数矩形严格
位于单位圆盘内，阶的实部至少为 -256。两种接口的精度上限均为 512 位，
预算耗尽时返回错误。

## 实数三角函数的可靠包络

`num_mp_real_sin`、`num_mp_real_cos` 和 `num_mp_real_tan` 返回精确二进制
弧度角的三角函数包络。支持的角度绝对值上限为 2^20 弧度。
`num_mp_real_sin_cos` 只做一次范围缩减，同时得到正弦和余弦的包络。
对应的 `num_mp_ball_*` 算子接受中点位于支持范围内的闭角度区间。
计算出的余弦包络含零时，正切返回定义域错误。

## 实数双曲函数的可靠包络

`num_mp_real_sinh`、`num_mp_real_cosh` 和 `num_mp_real_tanh`
为精确实输入返回外向二进制包络。对应的 `num_mp_ball_*` 算子接受闭区间。
区间跨越零时，双曲余弦保留最小值 1；双曲正切在大有限输入下接近 ±1 时
仍返回有限包络。

## 主辐角的可靠包络

`num_mp_real_atan` 和 `num_mp_ball_atan` 分别返回实值或闭区间的反正切包络。
`num_mp_real_atan2(y, x, precision)` 返回 (-π, π] 内主辐角的包络，
原点超出定义域，负实轴上的点取 +π。矩形接口 `num_mp_ball_atan2` 拒绝含原点
的矩形，也拒绝接触负实轴但并非全在轴上的矩形。复数矩形接口
`num_mp_complex_arg` 拒绝触及非正实轴的矩形。

## 反正弦与反余弦的可靠包络

`num_mp_real_asin` 和 `num_mp_real_acos` 为 [-1, 1] 内的精确实值
返回主值包络，包括两个端点。对应的 `num_mp_ball_*` 变体接受完全位于
该定义域内的闭区间，区间有任何一部分超出 [-1, 1] 时返回定义域错误。
输入为 1 时，反余弦返回零半径的精确零。

## 实数幂的可靠包络

`num_mp_real_powi` 和 `num_mp_ball_powi` 接受有符号 64 位整数指数。
指数为零时返回 1，负指数拒绝含零的底数区间。分数幂接口
`num_mp_real_pow_ratio` 和 `num_mp_ball_pow_ratio` 接受有符号 64 位分子
和正的 32 位分母。负底数要求约分后的分母为奇数。

`num_mp_real_pow` 和 `num_mp_ball_pow` 接受二进制指数。非零底数采用
非整数指数时，底数须为正。指数区间包含多个值时，底数区间须严格为正。

## 实对数与二进制指数

`num_mp_real_ln`、`num_mp_real_log2`、`num_mp_real_log10` 和
`num_mp_real_exp2` 返回精确实输入的外向二进制包络。
对应的 `num_mp_ball_*` 算子接受闭二进制区间。对数要求输入严格为正，
区间的两个端点也须为正。输入为整数且结果可表示时，二进制指数返回精确的二次幂。

`num_mp_real_scale_decimal(value, exponent, precision)` 将二进制值乘以
十的有符号整数次幂，只对精确结果舍入一次。
`num_mp_real_scale_decimal_in` 接受指定的舍入方向。

## 多精度数学常数

`num_mp_real_const_named(name, precision)` 对数学常数按等距取偶返回经过
验证的二进制值。支持 `pi`、`tau`、`e`、`ln2`、`sqrt2` 和 `phi`。
可用 `num_mp_real_const_enclose` 保留外向包络。
`num_mp_real_const_pi_scaled(bits)` 返回 π × 2^bits 的精确整数下取整，
小数位数范围为 0 至 32,768。名称列表与多精度解析器支持的常数及表达式函数一致。

## 物理常数

`dim::constants` 保存 CODATA 2022 常数的数值、量纲和 SI 单位。按定义精确的常数与有标准不确定度的测量常数分别标记。`c` 精确为整数 `299792458`，`hbar` 以有理数乘 `π^-1` 保存 `h/(2π)`。宿主可从 `KernelManifest.constants` 生成常数模块。这里的 `G` 是万有引力常数。单位表的 `G` 是高斯，位于不同的命名空间。

## 共享库插件

`kernel::cdylib` 的 C ABI 输入 bincode-next 编码的 `Vec<Value>`，返回同格式编码的 `Result<Value, String>`。Rust 枚举内存布局不作为 ABI。共享库代码在宿主进程运行。启用 `plugin-cdylib` 后，`Features::plugin_kinds` 报告共享库加载能力。`tests/cdylib_plugin.rs` 执行真实 C 插件的端到端调用。[C 构建指南](BUILDING.zh.md)

## 构建能力

| `KernelManifest.features` 字段 | 能力 |
|---|---|
| `dynamic_call`、`parallel_symbolic`、`higher_order` | 动态调用、并行符号阶段和宿主回调 |
| `jit`、`jit_branch`、`jit_index` | 可选 JIT 及可降低的分支或索引 |
| `exact_units`、`affine_units`、`user_base_dimensions` | 精确单位、仿射单位与自定义基本量纲 |
| `structured_results` | 结构化方程解和积分结果可跨 ABI 交换 |
| `sde` | SDE 规则集合 |
| `plugin_abi`、`plugin_kinds` | 插件 ABI 以及构建支持的插件来源 |

`sde_wiener` 生成带种子的 Brownian 增量，形状为 `steps × m`，量纲为步长量纲的平方根。`sde_levy_areas` 生成在给定增量条件下的迭代积分，其量纲等于步长量纲。`sde_scheme_info` 给出算法的强阶、随机微积分约定、交换性、自适应策略和 Lévy 面积要求。原生 SDE 绑定可持有漂移和扩散程序、复用噪声样本并运行积分。`num_analysis_sde_ensemble` 使用共享求值预算计算多条轨迹。Rust 驱动还可直接接受本地闭包。[回调契约](CALLBACKS.zh.md) · [数值回调算子](NUMERICAL_CALLBACKS.zh.md)

## 可进入表达式图的值

清单的 `value_kinds` 以 `in_term_graph` 标记可物化为图节点的种类。`Int`、`Real`、`Float`、`Bool`、`String`、`Symbol`、`Unit`、`Quantity`、`Tensor`、`Term` 和 `Dyn` 可进入图，仍须满足具体算子的形状和量纲规则。虚部为零的 `Complex` 可按实数进入图。`Sparse`、`Record`、`Variant`、`Function` 和 `Opaque` 由动态求值器或动态结果接口处理。

## 结构化结果

`Value::Record` 携带指定字段及其 `FlattenSchema`。结构化值支持二进制序列化，并随值保留布局。

| 类型 | 字段 |
|---|---|
| `Solutions` | `variable`、`roots`、`multiplicities`、`conditions`、`complete` |
| `OdeSolutionView` | `t`、`y`、`flatten_order`、`events`、`nfev`、`rejected` 及量纲 |
| `SimResult` | `scenario`、`t`、`fields`、`schema` |

`Session::solve_structured` 将 `solve` 的列表结果提升为 `Solutions`。空 `multiplicities` 表示求解器未区分根的重数。只有穷尽解集时 `complete` 才为真。`OdeSolution::view` 为积分器的列命名，量纲数据随记录传输。

## 量纲与单位

`UnitDef` 保存名称、别名、以 `DimExpr` 表示的量纲、精确比例、可选仿射偏置和显示信息。`UnitRegistry` 提供枚举、查询、复合解析、注册与快照。冲突的名称或别名会被拒绝。`KernelManifest.units` 可供宿主直接生成单位表。

`DimExpr` 用名称引用基本量纲，这样另一进程也可读取同一表达式。`Config::dimension_term` 接受结构化表达式，可避免文本单位名的歧义。`DimExpr::to_term` 在接收域中注册此前未知的基本名字。[量纲与单位](DIMENSIONS.zh.md)

```rust
let config = Config::new().dimension_term("x", DimExpr::base("Currency"));
```

## 线程与句柄

每个 `Session` 限于一个线程。`Term<'s>` 借用会话。`TermRef` 是 `Copy + Send + Sync` 的会话 ID 与节点引用，只能由生成它的会话解析。`fork` 产生可移交给另一线程的独立图。`merge` 导入对方节点并转换量纲假设，遇到冲突返回错误。一个会话内的并行批量接口是 `Session::compute_all`。

Rust 宿主可用 `Session::with_native_object_scope` 在当前会话运行临时操作。该操作返回或发生栈展开时，作用域释放本次创建的原生对象。已有对象仍可使用，作用域内创建的句柄须在该操作内使用。

## 带量纲的 CAD 网格

`num_cad_mesh_new` 使用两个三列张量，分别存放带长度量纲的顶点坐标和从零开始的整数
三角形索引。拓扑查询重用缓存的边关联统计。有符号体积的量纲为 `Length³`，要求网格闭合
且面朝向一致。焊接与修复接受带量纲的容差，返回新网格和变更统计。

`num_cad_mesh_boolean` 对两个连通、面朝外的凸壳求并集、差集或交集。操作参数
分别为 `union`、`difference` 和 `intersection`。结果为新的带量纲网格，相离实体的
交集为空网格。不支持的几何或超出工作量限制时返回错误。Aevia 中可用
`Rssn_Cad.mesh_boolean` 调用。

STL 入口使用内存中的文本或字节张量。STL 自身没有单位，请提供一个坐标单位所代表的
长度。从 STL 读出的三角形可能各自保存顶点。查询拓扑前应先焊接，使用完毕后调用
`num_native_release` 释放句柄。

`num_cad_step_read_text` 读取有工作量上限的平面 STEP 实体，边使用 `LINE` 几何，
结果是网格句柄张量。通过 STEP 长度单位转换为米，未声明单位时按米解释。

## C ABI

```c
char *json = rssn_manifest_json();
uint64_t hash = rssn_manifest_hash();
```

两个函数查询同一份标准清单快照，其中包含首次调用时可见的进程级单位。JSON 编码可以复用，每次调用分别返回一份字符串，宿主用同一 ABI 的 `rssn_string_free` 释放。注册了插件的会话可用 Rust API 查询该会话的完整清单。[C 接口概览](../ARCHITECTURE.md#c-interface)

## JIT

算术、`sqrt`、`abs`、`floor`、`ceil`、比较、选择和整数次幂有原生降低。`Inst::Gather` 访问 `inputs[base + offset]`，先检查索引再计算安全地址，超界取 `NaN`。`OpDescriptor::jit` 可声明 `Intrinsic`、`Extern1`、`Extern2` 或 `Gather` 等降低。其余路径可调用标量求值器，或在解释器中运行。缓存键同时包含运算带结构、清单和目标特征。[执行后端](BACKENDS.md)

## 添加算子

对于可组合编解码器支持的数值函数，在 `codec-policy.tsv` 中登记 Rust 路径、ABI 名称、命名空间及量纲策略。运行 `bazel run //:refresh_contracts` 生成重载，宿主仍通过内核清单与值调用接口访问它。 [绑定生成与支持的契约结构](../tools/rssn-api-gen/README.md)

符号规则、回调与原生对象适配器按下列步骤直接注册。

带类型的原生对象适配器从同一参数列表生成签名与参数个数。运行时先检查参数个数，再按深度优先顺序检查各参数的深度、元素数量与张量形状，最后调用原生方法。

1. 用 `Installer::op` 或 `Registry::register_op` 注册。
2. 给出 `Signature` 或 `DynamicReason`，并声明量纲规则。
3. 填写 `.doc(..)`、`.category(..)` 和 `.since(..)`。这些字段随清单提供给宿主。
4. 运行 `bazel test //:contracts_test` 和该算子的测试。

## 预计算传递函数

`num_control_transfer_new` 接受按降幂排列的分子和分母系数，创建连续时间单输入、
单输出模型。需要阶跃响应时可用 `_new_proper` 检查分子次数是否不超过分母次数。模型可复用于
频率求值、稳定性分析以及串联、并联与负反馈，重复查询极点和零点会重用会话内缓存。

调用前选择统一的参考时间，按该参考量归一化系数、频率及采样时间。ABI 在此使用
无量纲数值。内部对数网格上没有相应交越频率时，裕量为 `None`。遇到奇点或数值
失败时返回错误。增益裕量以 dB 表示，相位裕量以度表示。用完可调用 `num_native_release`。

## 归一化 PID 与系统辨识

`num_control_pid_response` 使用已准备的传递函数和归一化时间网格。
`num_control_simulate_pid` 接受误差采样、正时间步长与可选输出范围。两者使用并联增益
`kp`、`ki`、`kd` 及可选的正滤波时间 `filter_time`。当积分项会使输出进一步越过限幅时，
积分保持不变。

`num_control_era` 根据归一化的输入采样和逐采样、逐输出排列的张量辨识离散单输入模型。
第一个输入采样须非零。`num_control_era_matrices` 返回 A、B、C，
`num_control_impulse_response` 返回逐输出脉冲响应，首行是零。模型无直接馈通，使用结束后
通过 `num_native_release` 释放。

## 预计算离散模型预测控制

`num_control_mpc_prepare` 读取归一化的离散 A、B 矩阵、正采样间隔，以及包含预测与控制时域、
对角权重和可选分量限值的配置。输入权重须为正。计划可重复用于 `num_control_mpc_closed_loop`，
以复用预测矩阵、Hessian 和约束行。每一步会按当前状态更新约束的右侧项。

原生接口限制矩阵尺寸、时域及轨迹计算量。无可行解或结果非有限时返回错误。
轨迹包含初态。使用结束后调用 `num_native_release`。

### 有限区间数值变换

`num_transform_*_checked` 返回数值、绝对误差估计、收敛状态和实际求值次数。
Mellin、逆 Mellin、Hilbert、Hankel 和 Stieltjes 的回调坐标无量纲，结果保留
回调振幅的量纲。物理坐标需要显式归一化，再按变换公式还原尺度。
自适应误差只覆盖给定的有限区间。Laplace 轮廓和 Gaver–Stehfest 则比较
两种分辨率，给出经验误差。平移参数应位于相关奇点右侧。

Mellin–Barnes 的有限垂直积分线段不能经过极点。Meijer G 的积分线还须分开两组分子极点。
误差估计不含截掉的积分尾部。Borel 多项式保留系数的共同量纲，
自变量须无量纲。

## 经验证的复数圆盘

复数圆盘的张量末轴为 `[实部, 虚部, 半径]`。特殊函数
`num_complex_disk_erf`、`num_complex_disk_gamma`、
`num_complex_disk_log_gamma` 和 `num_complex_disk_bessel_j` 接收无量纲坐标，
并返回覆盖整个输入圆盘的可选包络。碰到极点、跨越分支割线或超出有限计算范围时，
结果可能为 `None`。

标量求根认证 `num_complex_disk_certify_root` 使用一阶 Jet。系统求根认证
`num_complex_disk_certify_system` 接收最多 64 个圆盘所成的无量纲复数盒。
回调需在整个圆盘或盒上包络函数和导数或 Jacobian。标量认证需要调用两次，
系统认证需要调用三次。Krawczyk 检验成功时，返回包络内存在唯一根。

`num_complex_disk_certify_rouche` 检验圆周上是否处处满足
`|g(z)| < |f(z)|`。调用方提供解析函数，并由回调包络覆盖各段圆弧的整个
圆盘上的函数值。接口接收 4 到 4,096 段圆弧，每段至少预留两次回调。
返回 `false` 表示给定包络未能建立该不等式。

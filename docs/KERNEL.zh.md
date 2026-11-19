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

## 插件

插件通过 `Registry` 注册算子与重写规则。`SimplePlugin` 接收一段注册函数，向 `Registry::register_op` 提供 `OpDescriptor`，并可用 `register_rewrites` 设置规则和重写层级。每个注册算子需要签名或 `DynamicReason`。同名但形状不兼容的注册失败。为同一算子增加签名时使用描述符的 `.overload()`。[Rust 插件示例](KERNEL.md#plugins)

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

## C ABI

```c
char *json = rssn_manifest_json();
uint64_t hash = rssn_manifest_hash();
```

两个函数查询标准构建。宿主负责用同一 ABI 的 `rssn_string_free` 释放返回字符串。注册了插件的会话可用 Rust API 查询该会话的完整清单。[C 接口概览](../ARCHITECTURE.md#c-interface)

## JIT

算术、`sqrt`、`abs`、`floor`、`ceil`、比较、选择和整数次幂有原生降低。`Inst::Gather` 访问 `inputs[base + offset]`，先检查索引再计算安全地址，超界取 `NaN`。`OpDescriptor::jit` 可声明 `Intrinsic`、`Extern1`、`Extern2` 或 `Gather` 等降低。其余路径可调用标量求值器，或在解释器中运行。缓存键同时包含运算带结构、清单和目标特征。[执行后端](BACKENDS.md)

## 添加算子

对于可组合编解码器支持的数值函数，在 `codec-policy.tsv` 中登记 Rust 路径、ABI 名称、命名空间及量纲策略。运行 `cargo run --locked -p rssn-api-gen -- --codecs-write` 生成重载，宿主仍通过内核清单与值调用接口访问它。 [绑定生成与支持的契约结构](../tools/rssn-api-gen/README.md)

符号规则、回调与原生对象适配器按下列步骤直接注册。

1. 用 `Installer::op` 或 `Registry::register_op` 注册。
2. 给出 `Signature` 或 `DynamicReason`，并声明量纲规则。
3. 填写 `.doc(..)`、`.category(..)` 和 `.since(..)`。这些字段随清单提供给宿主。
4. 运行 `cargo run --locked --bin signature-audit` 和该算子的测试。

# RSSN 实现参考

[English](ARCHITECTURE.md) | [简体中文](ARCHITECTURE.zh-CN.md)

[系统架构手册](https://aevia.apich.org/architecture/03-rssn-engine?lang=zh)介绍了 Aevia 上下文中的 RSSN 与注册表设计；本文档记录内核内部实现细节。

**rssn** 围绕一个核心理念构建：**所有符号运算都是恒等变换**。`diff(f, x)`、`solve(eq, x)` 等算子在带有等价类（e-graph）的哈希驻留 DAG 上由规则集重写，再由提取器按目标*阶段*（闭式或数值）选出最优等价项。

非恒等计算（时间步进、离散场、粒子系统）位于 `sim`，与符号侧共享 `kernels` 中的纯数值算法。

```
 Rust 用户 ────────► api ──────────────────────────────────┐
 C / C++ / Python ► ffi（基于 api + sim 的 C ABI）─────────┤
 宿主、插件 ───────► kernel（动态调用 ABI、清单）──────────┤
                                                          │ 请求项 + 目标阶段
               ┌──────────────────────────────────────────▼───┐
               │ graph：哈希驻留 DAG + 等价类                  │  算子注册表、模式、
               │ 分层调度器 · 提取 · 量纲推导                   │  事实、数值区间球
               └───────┬──────────────────────────┬───────────┘
            重写 + 内核 │                          │ 数值证据
               ┌───────▼───────┐   ┌──────────────▼─┐   ┌──────────────┐
               │ rules/<domain>│──▶│ kernels/<area> │   │ sim/（数据，  │
               │ （恒等式）     │   │ （纯数值算法）  │   │  非符号项）   │
               └───────────────┘   └────────────────┘   └──────────────┘
                      io/（解析、LaTeX、Typst、排版、绘图）· backend/（编译、JIT）
                      dim/（单位、物理常数、单位制 — graph::dim 的门面）
```

## 模块结构

| 模块 | 职责 |
|---|---|
| `graph` | 与具体数学领域无关的等价图引擎：`store`（哈希驻留节点、并查集等价类、字面量折叠）、`op`（开放算子注册表，含 `HEAVY`/`OPAQUE_ON_APPLY`/`COMMUTATIVE` 标志、代价与绑定器）、`pattern`/`subst`、`rule`（`RuleSet`、`Installer`）、`schedule`（Reduce → Normalize → Explore 分层不动点调度，带预算与禁用机制）、`window`（局部树窗口搜索）、`extract`（`ClosedForm`、`SizeCost`）、`facts`（假设与符号/整数性分析）、`eval`（区间球求值）、`soundness` 及 `dim*`（量纲类型层）。 |
| `rules` | 每个数学领域对应一个 `RuleSet`；`rules::standard()` 安装全部标准规则集。 |
| `kernels` | 基于 `f64`、切片与 `ndarray` 的纯数值算法：求积、ODE/SDE/PDE（`pde.rs`、`pde/three.rs`）、特殊函数、线性代数（稠密/CSR/BSR、Krylov、特征值）、FFT、优化（`kernels::optimize`）、有限域、图论、同调、李代数、`convex`（LP/QP/SOCP/SDP）、`milp`、`glm`、`timeseries`、`survival`、`mcmc`、`control`、`meshgen`（2D/3D Delaunay）、`cad`（STL/STEP、BSP CSG）、`interpolate` 与 `mp`（任意精度）。Checked API 记录输入校验规则；底层构造器（如 `csr_from_triplets`）在越界时可能 panic（[docs/MATH_COVERAGE.zh.md](docs/MATH_COVERAGE.zh.md)）。 |
| `api` | `Session`、`Term<'s>`（重载运算符的 `Copy` 句柄）、`Config`、`Answer`、`ComputeError`、`api::cache` 及 `api::value`（`call_value`、`call_key`、`call_ho`、`fork`/`merge`、插件注册）。 |
| `kernel` | 嵌入表面：签名/类型语言、动态 `Value` ABI、`CallCtx`/`OpError`、插件与回调、`KernelManifest`、`bincode-next` 编解码及 `kernel::cdylib`（`plugin-cdylib`）（[docs/KERNEL.zh.md](docs/KERNEL.zh.md)）。 |
| `dim` | 公开量纲与单位层（`UnitDef`、`UnitRegistry`、CODATA 常数、自然/普朗克/原子/几何单位制），重导出 `graph::dim` 类型。 |
| `io` | 可往返文本、LaTeX、Typst、Unicode 排版、绘图、`.npy`（`npy`）、`io::frame`（`DataFrame`、CSV、JSON）、`io::parquet`（v1/v2 页、嵌套列表/映射、Snappy/gzip/Zstd/LZ4）、`io::hdf5`（superblock 0/2/3、分形堆、v2 B 树、分块过滤）与 `io::netcdf`（CDF-1/2/5）。编解码器将文件长度与偏移视为不可信输入，出错返回 `Err` 而非 panic。 |
| `backend` | `Backend` trait 将闭式编译为 `Interpreter`、SSA 磁带 IR 或 Cranelift JIT（`jit`）。`backend::device` 提供类型化内核 IR、CPU 设备、wgpu（`gpu.rs`，`ds.rs` 提供双单精度 `f64` 模拟）、WGSL/GLSL/OpenCL/CUDA/Metal 代码生成、启动融合（`plan.rs`）与数组程序（[docs/BACKENDS.zh.md](docs/BACKENDS.zh.md)）。 |
| `sim` | 时间步进、场与粒子：FDM/FEM/FVM/BEM、CFD、MD、SPH、谱方法、多重网格、`fem_triangles`、`gas_dynamics`、`integrators`、Ising 及 `sim::dimensional`（`[nx]`/`[nx, ny]`/`[nx, ny, nz]` 内核）。`sim::scenario::run(name, json)` 为统一入口，`sim::schema::catalog()` 描述所有仿真场景。 |
| `ffi` | 覆盖 `api`、`kernel::manifest` 与 `sim::scenario` 的精简 C ABI；由 cbindgen 生成 `rssn.h`/`rssn.hpp`（`bazel run //:refresh_headers`；[docs/BUILDING.zh.md](docs/BUILDING.zh.md)）。 |
| `prelude`、`constant` | `use rssn::prelude::*` 常用重导出；构建元数据与数学常数。 |

### 规则集

标准会话共享已安装的规则模板与当前单位注册表修订版的不可变清单，同时各自持有独立的项图、原生对象存储与答案缓存。`rules::standard()` 安装：

- **核心代数与分析**：`arith`、`elementary`、`calculus`、`poly`、`solve`、`ode`、`linalg`、`sparse_symbolic`（模行列式、Dixon 提升、Wiedemann）、`holonomic`（D-有限方程与恒等式）、`geometry`、`complex`、`number_theory`、`combinatorics`、`logic`、`special`、`stats`、`transforms`。
- **数值算子**（`src/rules/numeric/` 子模块）：`bridge.rs`、`decomp.rs`、`functions.rs`、`ode.rs`、`resample.rs`、`solvers.rs`、`special.rs`、`transform.rs`、`conic.rs`、`milp.rs`、`glm.rs`、`timeseries.rs`、`control.rs`、`survival.rs`、`mcmc.rs`、`mesh.rs`、`mp.rs`（[docs/MATH_COVERAGE.zh.md](docs/MATH_COVERAGE.zh.md)），以及生成的动态 ABI 族（`numerical_vector`、`numerical_stats`、`numerical_special`、`numerical_matrix`）、`numerical_callback`（[docs/NUMERICAL_CALLBACKS.zh.md](docs/NUMERICAL_CALLBACKS.zh.md)）与 `sde`。
- **应用数学与物理**：`variational`、`pde`、`functional`、`physics`、`physics_quantum`、`optimize`、`verify`。
- **离散与抽象代数**：`discrete`（图、有限群与置换群、同调、编码、密码）、`rewriting`（Knuth–Bendix）、`lie`、`lie_structure`、`geometric_algebra`。
- **单位与工具**：`units`、`radicals`、`cse`。

每个规则集包含**算子**、带守卫的**重写规则**、**定义**、**内核**（返回 `Equal`、`Pinned`、`Approx` 或 `Pass`）及算子**属性**（量纲规则、签名、后端降级、积分表）。

数论内核采用两层设计：`u64` 字长快速路径（`kernels::modular` Montgomery/Barrett/NTT）与同一 `Domain` 接口下的通用精确路径，二者相互进行差分测试。重型算法（Risch–Bronstein/Risch–Norman、Hermite 归约、Meijer G、Gosper、Zassenhaus/Kronecker 分解、Gröbner/F4/F5/FGLM/CAD、ODE/PDE 求解器、二次筛、Schreier–Sims/Todd–Coxeter、`kernels::permgroup`）均作为内核实现并将结果以等式写回等价图。

## 阶段

- **`Target::Symbolic`**：使用 `ClosedForm` 提取（不允许残留重算子）；若无闭式，返回最小尺寸项并标记 `reduced = false`。
- **`Target::Numeric { tolerance }`**：对最优项求值，失败时回退至内核生成的数值证据（区间球）。近似数值不会作为精确项混入等价图。

## 会话、隔离与并发

每次 `Session::compute` 克隆规则配置模板图、导入请求、在克隆上运行引擎，并仅将答案导回会话。`Session` 为单线程（内部使用 `RefCell`）；跨线程请使用 `Session::fork` / `Session::merge`，语言绑定使用 `TermRef`（`Copy + Send + Sync`）。

- **缓存（`api::cache`、`CachePolicy`）**：答案按请求节点、模板及影响结果的 `Config` 指纹记忆化（超时运行不缓存）。可选**引理（Lemmas）**缓存重型子请求的闭式。注册插件或单位会清空两张表。
- **批处理（`Session::compute_all`）**：入口处快照命中与引理，在工作线程克隆上并行执行未命中项（每个 Rayon worker 最多保留 2 个私有任务），按输入顺序导入答案。`execution` 遵循 `RAYON_NUM_THREADS`；`kernels::reuse` 池化有界 FFT 计划与临时缓冲区。
- **增量匹配与索引**：等价类记录 `Graph::touched` 纪元，深度 `d` 重写与窗口遍历跳过未变更子图（`RSSN_FULL_SEARCH=1` 可关闭）。`Engine` 构建一次 `TierIndex`，按需构建算子索引，并按图修订号记忆化 `Graph::facts`、`looks_reduced` 与失败窗口单元。
- **并行分层与性能剖析**：`Engine::apply_tier` 在 512 个候选节点以上并行搜索重写批次，并按确定性顺序串行提交。两阶段内核（`Kernel::two_phase`）并行执行 `PreparedWork` 后按节点顺序提交。设置 `RSSN_PROFILE_RULES=1` 开启逐规则计时统计（`graph::profile`、`examples/engine_timing.rs`）。
- **取消机制**：运行在达到不动点、预算上限或协作式 `rssn::cancellation` 检查点时终止。

## 嵌入：`kernel` 与 `dim`

宿主读取 `Session::manifest`（算子、签名、单位、常数与构建能力的确定性可哈希目录），并通过 `Session::call_value`、`call_key` 或 `call_ho` 调用算子。每个算子必须具备类型签名或 `DynamicReason`（`bazel run //:signature_audit`）。签名由测试轨迹（`signatures.witnessed.json`）与人工覆写（`signatures.overlay.json`）共同提供，`build.rs` 要求两者均存在。

`kernel` 与 `dim` 作为模块而非独立 crate 保留，因为二者共享表达式图、算子表与量纲竞技场，避免跨 crate 门面开销。详见 [docs/KERNEL.zh.md](docs/KERNEL.zh.md)、[docs/CALLBACKS.zh.md](docs/CALLBACKS.zh.md) 与 [docs/NUMERICAL_CALLBACKS.zh.md](docs/NUMERICAL_CALLBACKS.zh.md)。

## 诊断

当求值器拒绝参数时，引擎保留未归约调用项。`kernel::diag` 在不改变归约行为的前提下提供原因说明：

- 求值器调用 `diag::fail(kind, arg, ...)`（或 `rules::numeric::bridge` 辅助函数），仅在激活时向**线程局部收集器**写入 `Diag { op, arg, kind, message }`。
- `Session::apply_diag` 与 `Session::explain` 在调用线程的 `diag::capture` 下重放该调用。
- 表驱动数值算子（`rules::numeric::TABLES`）记录精确参数错误；重型算子通过 `diag::classify` 推导 `Arity`、`Shape`、`Symbolic` 或 `Unsupported`。`From<Diag> for OpError` 将诊断映射为 `OpError::HostDiagnostic`（`rssn.diag.<kind>`）。

## 任意精度

`graph::number::Number` 仅包含整数、有理数与 `f64`。任意精度作为值层能力位于 `kernels::mp`：

| 层 | 内容 |
|---|---|
| `kernels::mp::real` | `Real`（`mant * 2^exp`），在 `BigInt` 中精确计算并按 `Round` 模式单次舍入。 |
| `kernels::mp::funcs` | 初等与超越函数返回严谨的 `Ball` 区间包络，以及按需缓存的整数级数常数（`pi`、`e`、`ln 2`、`sqrt 2`、黄金比例）。 |
| `kernels::mp::decimal` | 精确十进制输入输出；仅当区间两端点舍入一致时输出数字。 |

`rules::numeric::mp`（`mp_eval`、`mp_ball`、`mp_to_float` 等）在项图中以**十进制字符串字面量**传递数值，使高精度值可直接参与哈希驻留与缓存而无需改动 `Value` 或引擎。`mp_eval` 采用 Ziv 策略提升精度（上限 `MAX_BITS` 约 9,800 位十进制数字，指数上限 `2^18`）。

## 量纲

`graph::dim` 为标量项推导精确有理数或有理函数指数的物理量纲（`x^n`、`x^(a+b)`、`x^(1/n)`），并通过合一求解未知符号与多态插件签名（`dim_signature`）。每个标准算子必须绑定 `DimRule`。`Advisory` 模式在 `Answer::diagnostics` 中报告冲突；`DimensionMode::Strict` 抛出 `ComputeError::Dimension`；`Config::bypass_dimension_checks` 可显式跳过检查（[docs/DIMENSIONS.zh.md](docs/DIMENSIONS.zh.md)）。

## 编译求值与 JIT

`backend::current()` 按 `backend::with_backend` → `backend::set_default` → `Interpreter` 顺序解析后端。

- **Tape IR（`backend::tape`）**：将项降级为 SSA 磁带（常量折叠、死代码消除、结构哈希），支持 `f64`、对偶数（`tape::autodiff` 反向微分与 Hessian-向量积）与外向舍入区间。`Tape::eval_columns` 按最多 256 行分块求值。
- **Cranelift JIT（`backend::jit`，`jit` 特性）**：将磁带编译为 `fn(inputs, params, out)` 并按结构哈希与清单/目标键缓存；`TieredBackend` 先解释执行并在后台编译。`unsafe` 仅限于 `abi` 与两处调用点。
- **LIR（`backend::lir`、`backend::lir_jit`）**：带控制流的类型化寄存器机，服务于原生行循环、`device::native` 及 Aevia `@jit` 函数（[docs/BACKENDS.zh.md](docs/BACKENDS.zh.md)）。
- **调用方**：数值积分、无穷级数、`odeint`、`nminimize`、`chebyshev_approx`/`rational_approx`、数值求解器与 `Term::compile`。

## C 接口

`src/ffi`（`term.rs`、`compute.rs`、`manifest.rs`）提供基于句柄与 JSON/`bincode-next` 的 C ABI（会话内 `uint32_t` 项句柄、`rssn_term_apply`、`rssn_compute*`、`rssn_manifest_*`、`rssn_sim_run`）。在 `ffi` 配置下，受保护入口会捕获展开的 panic 并返回 `RSSN_STATUS_PANIC`。详见 `examples/c/smoke.c` 与 [docs/BUILDING.zh.md](docs/BUILDING.zh.md)。

## 扩展 rssn

- **新增算子或领域**：通过 `Config::with` 或 `Session::with_rules` 添加 `RuleSet`，并附带 `DimRule` 与 `Signature`（或 `DynamicReason`）。
- **外部插件**：通过 `Session::register` 或 `plugin-cdylib` 动态库加载（[docs/KERNEL.zh.md](docs/KERNEL.zh.md)）。
- **积分表 / 降级 / 后端**：扩展 `rules::calculus::IntegralTable`、附加 `backend::Lowering` 或实现 `backend::Backend`。

## 代码检查策略

`src/lib.rs` 拒绝 `clippy::all`、`pedantic`、`nursery`、`unwrap_used`、`expect_used`、`indexing_slicing`、`arithmetic_side_effects`、`single_call_fn` 与 `missing_docs`，并对科学计算常见模式保留带说明的豁免（`suboptimal_flops`、`cast_precision_loss`、`many_single_char_names`、`similar_names`、`tuple_array_conversions`、`needless_range_loop`、`unreadable_literal`、`excessive_precision`、`while_float` 及基础类型转换 lint）。局部 `#[allow]` 须附原因注释。

## 测试

- **单元测试与正确性验证**：各模块单元测试、`graph::soundness`（`no_standard_rewrite_is_refuted`）、`graph::proptests` 与 `graph::dim` 覆盖率扫描。
- **集成测试（`tests/`）**：`tests/kernels.rs`、`tests/sim.rs`、`tests/dimensions*.rs` 与 `tests/coverage_*.rs`（`coverage_algebra`、`coverage_api_complex`、`coverage_numeric_a`、`coverage_numeric_b`）。`src/backend/device/tests.rs` 中的 GPU 测试需启用 `--features gpu-wgpu` 与真实适配器。
- **生成绑定与基准**：`bazel test //:generated_sources_check_test //:contracts_test`（[tools/rssn-api-gen/README.zh-CN.md](tools/rssn-api-gen/README.zh-CN.md)）、`bazel run //:benchmark -- engine`。快速过滤运行见 [README.zh-CN.md](README.zh-CN.md#测试)。

## 延伸阅读

| 文档 | 主题 |
|---|---|
| [README.zh-CN.md](README.zh-CN.md) | 特性与快速开始 |
| [docs/KERNEL.zh.md](docs/KERNEL.zh.md) | 嵌入内核：清单、动态 ABI、插件与单位 |
| [docs/DIMENSIONS.zh.md](docs/DIMENSIONS.zh.md) | 量纲代数、推导与嵌入 API |
| [docs/CALLBACKS.zh.md](docs/CALLBACKS.zh.md) · [docs/NUMERICAL_CALLBACKS.zh.md](docs/NUMERICAL_CALLBACKS.zh.md) | 宿主回调与数值回调算子 |
| [docs/BUILDING.zh.md](docs/BUILDING.zh.md) | Rust 与 C 产物、头文件与构建元数据 |
| [docs/MATH_COVERAGE.zh.md](docs/MATH_COVERAGE.zh.md) | 数值边界约定与覆盖范围 |
| [CONTRIBUTING.zh-CN.md](CONTRIBUTING.zh-CN.md) | 开发流程与检查项 |

# rssn 统一恒等算子架构重构计划

## Context

`/home/user/dev/rssn`（分支 `feature/egraph-symbolic-overhaul`，1220 个 .rs / 25 万行）处于半截迁移状态：

- `src/symbolic/egraph/`（8.9k 行）是建在旧 `Expr` 之上的 e-graph：`ENode` 把 156 变体的 `Expr` 又抄了一遍；变量是 `String`；`rules/oracles.rs`（4k 行）把 e-class 抽取回 `Expr` → 调旧算法 → 再塞回去；微分规则里靠变量名（`"y"|"z"|"psi"…`）猜函数依赖；`compute()` 用 `union(var, const)` 做绑定（不健全）；数值相态是抽取之后的 `if let` 兜底链。
- 旧 `Expr`/`DagManager` 全局单例、`simplify`/`simplify_dag`/`rewriting` 三套化简并存；`src/ffi_apis` 9 万行机械包装；`prelude.rs` 6k 行；`rssn.h/.hpp` 各 140 万字节。
- `ARCHITECTURE_MIGRATION_PLAN.md` 的方向（算子即恒等变换、相态、统一 `compute`、极简 FFI）可取，但它保留 oracle 桥接、封闭 `OpKind` 枚举、把数值结果当普通节点，这些不采用。

已与用户确认的决策：
1. 图内核在 rssn 内**全新设计**；rssn-advanced 本次不接入，只预留后端接口（未来做数值/模拟阶段的 JIT）。
2. 旧 `Expr` **彻底删除**，**全部域原生重写**，不留桥接层。
3. **全部域保留**；外围域（密码、纠错码、图形学、图算法、有限域、拓扑、分形…）改为“依赖图引擎的外层实用 API 包装”。
4. 数据结构采用“**环链 e-class + 树窗口**”。
5. 子代理：仅少量 Sonnet 5.5，用于非创造性工作；内核与核心算法由主线程亲自写。

目标终态：引擎 = **图计算**（规则集驱动的启发式 e-graph）+ **模拟**（纯数值步进等无法归为恒等算子的部分）。

---

## 1. 目标目录结构

```
src/
  graph/      图内核（本次核心创新，无任何域知识）
  rules/      域规则集：每个域 = 符号重写规则 + 规约内核(符号/精确/数值) 并排注册
  kernels/    不依赖图的纯数值例程（由 src/numerical 中 27 个已与 Expr 无关的文件整理而来）
  sim/        模拟：src/physics/* + numerical/physics_{cfd,fea,md}；非恒等算子
  toolkit/    外围域的实用包装 API（crypto、ecc、graphics、graph_algo、finite_field、topology、fractal…），内部全走 compute()
  api/        Session / Term / Config / compute / DSL 构造器
  io/         parser（输出到图）、latex/typst/pretty（从图读取）、plotting
  ffi/        极简 C-FFI（<1000 行）
  backend/    执行后端 trait + 解释器；JIT 占位（现有 src/jit 适配为一个后端，rssn-advanced 为未来后端）
  plugins/    插件 = 注册 Op + RuleSet
```

删除：`src/compute`、`src/symbolic/`（整体，逐域重写进 `rules/`）、`src/ffi_apis`、`src/ffi_blindings`、`src/prelude.rs`（重写为 <200 行）、`src/nightly`（并入 kernels 的 feature 门）、`rssn.h/.hpp`（重新生成）、`CODE_STASTICS.md`、`poc*.py`、全部 `*_ffi_test.rs` 与 `benches/compute_*`。旧实现仅通过 `git show <基线提交>:path` 参考，不保留 legacy 目录。

---

## 2. 图内核设计（`src/graph/`）

### 2.1 主体：开放算子的 hash-consed arena DAG
- `NodeId(u32)` / `ClassId(u32)` / `SymbolId(u32)` / `OpId(u32)`；`Graph` 为显式拥有的对象，**无全局单例**。
- 节点紧凑定长：`{ op: OpId, children: 子节点池中的 span, hash, ring_next: NodeId, uses_head, flags }`。子节点放连续 `Vec<NodeId>` 池，变长元数不另行分配。
- **Op 是开放注册表而非枚举**：`OpDescriptor { name, arity, flags(交换/结合/绑定子/线性/纯), sort 签名, 基础代价, eval 内核?, lower 钩子? }`。核心算子为常量 `OpId`，域与插件可动态注册 —— 这是 JIT/插件可扩展性的根基，替代 156 变体 `Expr` 与重复的 `ENode`。
- 叶子载荷走侧表驻留（`PayloadId`）：`Int(BigInt) / Rat / Float / Complex / Bool / Symbol / Tensor(Arc<ArrayD>) / Poly / Trajectory / Opaque(Arc<dyn Payload>)`。
- 绑定子（`Diff/Integral/Sum/Limit/Lambda/ForAll`）的约束变量是显式子节点；未知函数用 `Apply(f, args…)` 显式声明依赖，取代按名字猜测。

### 2.2 环链 e-class（动态链表 ①）
- e-class 成员不存 `Vec<ENode>`：每个节点带 `ring_next`，同类节点构成**侵入式循环链表**。`union` = 并查集合并 + 两环 O(1) 拼接，零拷贝；遍历某类 = 走环。
- use-list（父指针）同为侵入式链表（`uses_head` + use 单元池），合并时 O(1) 接链；延迟 `rebuild` 沿 use 链做同余闭包修复（egg 式 worklist，但不重建全 memo —— 现实现每次 rebuild 全表重哈希，要去掉）。
- 每类挂 **分析格**（`ClassData`）：常量值、自由变量集、sort/形状、**相态集合**、最优代价缓存、**数值见证** `approx: Option<Ball>`（见 2.5）。

### 2.3 树窗口（动态链表 ②，局部搜索/树优化）
- 从某根向下，把“单引用（use 计数 = 1）且未越过预算”的极大子区域切出为 **TreeWindow**：slab 支撑的双向链表树单元（`parent / first_child / next_sibling / prev`），可**就地破坏式改写**。窗口边界处的共享节点作为不透明叶（指回 `ClassId`）。
- 窗口内跑贪心 / 束搜索 / 有界回溯的局部重写（结合律重排、同类项合并、常量折叠、多项式规范化等会令全局 e-graph 爆炸的规则只在这里跑），不产生 e-node 膨胀。
- 提交：把窗口最优结果重新驻留进 DAG，并与原根类 `union`。全局饱和只负责**跨窗口的共享结构**与高阶算子。
- 即：DAG 负责共享与等价，链表树负责便宜的局部最优化；二者通过“切出—提交”闭环。

### 2.4 规则与规则集
- `Rewrite`：声明式模式 → 模式 + 守卫，`rules! { "sin2+cos2": sin(?x)^2 + cos(?x)^2 => 1; … }`，编译为按 `OpId` 索引的判别网匹配器（只在含该 op 的类上尝试，替代现在每条规则扫全图）。
- `Kernel`（规约内核）：过程式，“算法即恒等变换算子”的载体。签名 `fn reduce(&self, cx: &mut ReduceCx, node) -> Outcome`，声明 `from 相态 → to 相态`、触发 op、代价估计。RK45、Gauss–Kronrod、Risch、Buchberger、QR 都是 Kernel，和符号规则**同表注册、同一调度器**。
- `RuleSet { name, tier, deps, rules }`，`Registry` 负责组合；`Config::with(ruleset)` 选集。窗口规则与全局规则用标志区分。

### 2.5 数值/符号融合
- 相态：`Symbolic | Exact | Numeric{tol}`。目标相态是**调度目标**而非后处理：根类出现满足目标相态的成员即可提前停止。
- 近似值**不作为 e-node 并入等价类**（`∫…` 与 `0.333333` 不是恒等），而是写入类分析的 `approx: Ball{mid, rad}`（张量/轨迹同理为 witness）。这样同余闭包保持健全，同时数值内核的结果可被父节点的数值内核继续消费、可用于守卫判定（符号、非零、分支选择）与剪枝 —— 数值为符号搜索导航，符号为数值降低代价。
- 变量绑定是求值环境（传给数值内核 / 显式 `Subst` 算子），不再 `union(var, const)`。

### 2.6 启发式调度与抽取
- 分层：T0 高权算子“脱茧”（Diff/Integral/Solve/Ode…）→ T1 窗口内规范化 → T2 全局结构探索；每规则带 backoff（命中过多即禁用若干轮）；节点/类/时间预算；目标相态提前终止；数值内核仅在输入类已有数值见证或符号路径耗尽预算时触发。
- 抽取：`CostModel` trait（符号规模 / 数值求值代价 / 后端代价），DAG 感知（共享子式只计一次），可带相态约束。

### 2.7 后端与 JIT 扩展点（本次只做接口 + 解释器）
- `trait Backend { fn lower(&self, g: &Graph, roots: &[NodeId], sig: &Signature) -> Result<Box<dyn Compiled>> }`；`Compiled::call(&[f64]) / call_batch(列)`。
- `OpDescriptor.lower` 钩子让每个 op 自描述如何降级；`CostModel` 可由后端提供。
- 本次实现：树遍历解释器 + 把现有 `src/jit`（cranelift 指令级）包成 `feature = "jit"` 的一个后端适配。rssn-advanced 接入点记录在 `ARCHITECTURE.md` 的“未来规划”一节：数值相态内核与 `sim` 步进器通过 `Backend` 获取已编译闭包，并由调度器按调用次数动态决定是否 JIT。

---

## 3. 对外 API 与 FFI

```rust
let s = Session::new();
let x = s.sym("x");
let cfg = Config::new().with(rules::calculus()).target(Phase::Numeric { tol: 1e-9 }).bind(x, 3.0);
let r = s.compute(diff(sin(x) * cos(x), x), &cfg)?;   // Result<Term, ComputeError>
```
- `Term<'s>` 为 Copy 句柄，重载算术运算符；结果取值 `as_f64 / as_rational / as_tensor / to_latex`。
- 错误统一 `ComputeError`（不再 panic / 返回未求值表达式装作成功；未能规约时返回带 `Unreduced` 状态的结果）。
- C-FFI（`src/ffi/`，<1000 行，cbindgen 重新生成头文件）：`rssn_session_*`、`rssn_term_{sym,int,float,apply(op_name, args)}`、`rssn_config_*`、`rssn_compute`、`rssn_term_{to_f64,to_string,to_latex,tensor_data}`、`rssn_sim_run(name, json)`、错误码 + last_error。

---

## 4. 里程碑（每个里程碑结束时 `cargo test --all-features` 全绿再进入下一个）

**M0 基线**：把当前未提交的工作区做一次检查点提交（便于 `git show` 参考旧算法）；记录现有能通过的测试名单作为迁移对照。

**M1 图内核**（主线程）：`graph/` 全部 + 内核测试（proptest：hash-consing 唯一性、环完整性、同余不变量、窗口提交健全性、抽取最优性小例穷举）。

**M2 API + 基础规则集**（主线程）：`api/`、数值塔与常量折叠、代数/幂/指数对数/三角、微分、多项式规范化、解释器后端；双相态端到端测试。

**M3 切换**：删除 `compute`、`symbolic/core`、`symbolic/egraph`、`simplify*`、`rewriting`、`handles`、`ffi_apis`、`ffi_blindings`、旧 prelude 与对应测试/bench；`src/numerical` 中与 Expr 无关的文件整理进 `kernels/`，`physics` 进 `sim/`；`Cargo.toml` feature 与依赖清理（去掉 uuid、lazy_static、once_cell、双 rand 等不再需要者）。此时 crate 只含新架构，能编译。

**M4 全域原生重写**（分波次；每个域 = 规则集 + 内核 + 迁移并加强的测试，完成一个删一个旧参考）：
- 波 A（主线程，算法核心）：积分（表驱动 + 启发式 Risch + 数值求积内核）、极限/级数/求和、方程求解、ODE、PDE、多项式/因式分解/Gröbner/实根/根式、线性代数（符号 + faer 数值内核）、张量/微分几何。
- 波 B（主线程定模板，Sonnet 子代理每次 ≤2 个、按域独占目录执行）：变换表、特殊函数恒等式、数论、组合、逻辑、统计四模块、向量微积分、复分析、泛函/积分方程/变分、群论/李代数、几何代数、物理符号域（经典力学、电磁、量子、QFT、相对论、热力学、固体）、坐标/单位。
- 波 C（Sonnet 子代理）：`toolkit/` 外围包装（crypto、纠错码、图形学、图算法、有限域、拓扑、分形、CAD、proof、multi_valued），内部经 `compute()`；`io/` 的打印器与 parser 适配；`plugins/`。
- 子代理交付物由主线程逐一审阅并跑测试后才合入；不给子代理内核或调度器的修改权。

**M5 FFI 与周边**：`ffi/` + cbindgen + `build.rs`；`sim` 的统一入口；示例、`ARCHITECTURE.md` 重写（含 JIT/rssn-advanced 规划），删除 `ARCHITECTURE_MIGRATION_PLAN.md`；CI/脚本同步。

**M6 测试强化**：
- 每个算子的双相态测试（符号结果 + 数值结果互相校验）。
- **规则健全性自动检查**：每条声明式规则随机取值，数值求值两侧比对（规则注册即入测）。
- 由 M0 基线迁移来的全部域测试（`tests/` 从 ~230 个散文件重组为按域的模块目录）。
- 内核基准（饱和、窗口、抽取）替换旧 bench。

规模说明：这是数十万行级别的周转，会跨多个会话；M1–M3 是不可分的关键路径，M4 各域彼此独立、可按波次逐步推进，任何时刻主干都可编译可测试。

---

## 5. 关键文件

- 新建：`src/graph/{mod,id,op,payload,node,arena,class,uses,window,pattern,rule,registry,schedule,extract,analysis}.rs`、`src/api/*`、`src/backend/*`、`src/rules/<domain>/*`、`src/toolkit/*`、`src/ffi/*`
- 重写：`src/lib.rs`、`src/prelude.rs`、`Cargo.toml`、`build.rs`、`cbindgen.toml`、`ARCHITECTURE.md`、`README.md` 示例
- 可直接复用（与 Expr 无耦合）：`src/numerical/{matrix,sparse,special,interpolate,transforms,signal,optimize,stats,polynomial,real_roots,tensor,vector,…}.rs`、`src/physics/*`、`src/jit/*`、`src/constant.rs`
- 仅作算法参考后删除：`src/symbolic/*.rs`、`src/symbolic/egraph/rules/*.rs`（其中 `rules/algebra.rs` 的 `NumVal` 数值塔、`union_find.rs` 思路可借鉴）

## 6. 验证

- 每里程碑：`cargo check --all-features`、`cargo test --all-features`、`cargo clippy --all-features`（沿用 `lib.rs` 现有 deny 集，新代码不加 `#![allow(missing_docs)]` 之类豁免）。
- M1：proptest 不变量 + 小规模穷举对照朴素 e-graph 实现。
- M2 起：README 中的目标示例（`diff(sin x·cos x)` 符号/数值、定积分、ODE 轨迹）作为 doctest。
- M5：`cbindgen` 生成头文件，编一个 C 冒烟程序调用 `rssn_compute`；确认 `src/ffi` < 1000 行、头文件体积回到 KB 级。
- M6：规则健全性检查全量通过；`cargo bench` 内核基准可运行。

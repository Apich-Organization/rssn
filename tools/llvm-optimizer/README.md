# RSSN LLVM Optimizer

这是 RSSN release 可显式加载的 LLVM New Pass Manager 插件。验证工具链为 rustc 1.100.0（`0ef1bc386059596dada8fb65b50425aab2ea637a`）与 LLVM 23.1.5。插件在 O3 最终优化阶段运行，适用于 RSSN 的 fat LTO。Linux x86-64 也通过了 ThinLTO 和关闭 LTO 的加载验证。

## 优化范围

- 影响局部且能在 RSSN 源码表达的优化优先改 Rust。把求积对称性、特定算法的恒等式和局部分段条件放在拥有这些信息的源码层，再用同工具链、同优化配置复核生成结果。 [双指数求积的源码复用](WINDOWS-STUDY.md#源码复用双指数求积的对称节点)
- 发现局部问题后，用 IR、MIR 和机器码还原它的成因。按数据依赖、控制流和目标约束描述模式，在其他函数、模块或依赖中查找同类实例及其变体。跨调用活跃区间、目标调用约定和后端拆分宽载入等问题，优先考虑可复用的 LLVM 分析。仅在某个算法上成立的条件继续由源码表达。 [调用约定实例](ABI-STUDY.md) [相邻字段读取实例](ZSTD-STUDY.md)
- `libm`、`zstd` 等依赖的源码如果不在本仓库控制下，可在能获得依赖 IR 的构建阶段应用插件。候选按 IR 结构、属性和数据布局匹配。若只有一个依赖实例，也记录实际执行频率、适用的版本或目标，以及该构建阶段能看到的 IR。对只在二进制中可见的调用，插件只能改写可见的调用方。
- 把候选识别、合法性证明和收益估算放在报告和变换可以共用的分析中。报告区分 IR 代理量、MIR 的物理寄存器和 spill slot、最终指令以及运行时间。启用变换前，用真实热路径核对输出、机器码和运行 A/B。源码方案与插件可以叠加时，补做四组消融并记录插件的增量收益。 [Parquet 与 HDF5 的四组对照](ZSTD-STUDY.md#rssn-的-parquet-和-hdf5-读取入口)

## 构建与使用

在 RSSN 仓库根目录运行。Bazel 使用与 rustc 同版本的 LLVM SDK、本机 C++17 编译器和 CMake 构建宿主插件。测试还需要 Python 3.11 和 SDK 的 `opt`、`llc`。

```bash
rustc -vV
llvm-config-23 --version
bazel build --repo_env=RSSN_LLVM_CONFIG="$(command -v llvm-config-23)" //:llvm_optimizer
bazel test --repo_env=RSSN_LLVM_CONFIG="$(command -v llvm-config-23)" //:llvm_optimizer_test
bazel test --config=release-max --repo_env=RSSN_LLVM_CONFIG="$(command -v llvm-config-23)" //:release_pass_test
```

Bazel 按完整 LLVM 版本校验 rustc 与 SDK。动态库从宿主解析 LLVM 符号，避免在 rustc 进程里装载第二份 LLVM。`//:llvm_optimizer_tools_test` 跑独立 LLVM 工具的 IR、MIR 和目标对象测试，可以在无法直接加载插件的 rustc 宿主上使用。`//:llvm_optimizer_rustc_test` 核对 rustc 加载以及无 LTO、ThinLTO 和 fat LTO 的编译流程。旧 CTest 测试名分别转发给这两个 Bazel 目标。更换 rustc 或 SDK 后重新构建并运行测试。需要固定源码时，可以使用该 rustc 的 LLVM 子模块提交 `0a5579eb649d3634a9c30ac873300cc2944fa070` 生成 SDK。 [Rust 的 LLVM fork](https://github.com/rust-lang/llvm-project/tree/0a5579eb649d3634a9c30ac873300cc2944fa070)

### 宿主与编译目标

插件机器码取决于运行 `rustc` 或 `opt` 的宿主。各 pass 再根据模块目标和合法性条件决定是否执行。Linux 和 FreeBSD 等 ELF 宿主从进程解析 LLVM 符号，macOS 模块使用 Mach-O 动态查找。Windows 构建需要 `RSSN_LLVM_IMPORT_LIBRARY` 指向实际宿主的 LLVM C++ 导入库。若使用带 `LLVM_EXPORT_SYMBOLS_FOR_PLUGINS` 的 `opt.exe`，其导入库只服务该 `opt`，可用 `-DRSSN_TEST_RUSTC=OFF` 只运行 LLVM 测试。官方 Windows x64 rustc 1.100.0 的驱动 DLL 缺少所需 LLVM C++ 导出，原生 rustc 加载需要提供这些符号的构建。 [LLVM 的接口导出说明](https://www.llvm.org/docs/InterfaceExportAnnotations.html)

`rssn-pow-unit-base`、`rssn-fp-endpoint-guards` 和 `rssn-libm-unit-log` 使用 IEEE binary64 的 IR 条件，可跨目标运行。CTest 在 SDK 包含相应后端时，对 Linux x86/x64、ARM/ARM64，Windows x86/x64、ARM/ARM64，macOS x64/ARM64 和 FreeBSD x64/ARM64 验证 IR 并生成目标对象。调用约定、GPR/XMM 保存和内联汇编变换继续按各自的目标条件启用。交叉生成对象的验证与目标平台本机加载和运行分开记录。 [Appell 的跨目标验证](WINDOWS-STUDY.md#共享单位底数与跨目标验证)

Linux x86-64 已通过本机加载和完整 CTest。Linux ARM64 已用同版本 ARM64 SDK 构建插件，并在 QEMU 中用 ARM64 `opt` 加载、变换实际 Appell IR，以及用官方 ARM64 rustc 加载后生成 IR。Linux ARM32 插件已在 QEMU 的同架构 `opt` 和 LLD 中处理实际 Lerch IR。官方 ARMv7 rustc 1.100.0 静态链接 LLVM，驱动没有导出插件所需的 C++ 符号，比如 `IRBuilderDefaultInserter` 的虚表。这个宿主可让独立 `opt` 或支持 `--load-pass-plugin` 的同版本 LLD 加载插件，并用 `-DRSSN_TEST_RUSTC=OFF` 运行 LLVM 侧测试。Windows x64 的 PE 插件已在 Wine 中通过 CLANG64 `opt.exe` 加载并处理 RSSN 和 zstd 的 IR。FreeBSD 15.1 amd64 已验证同版本的独立 `opt`、`llc` 加载和实际 CPU lane 运行。Windows 其他架构、macOS 和 FreeBSD ARM64 已通过跨目标对象验证，宿主加载仍需验证。 [Windows PE 宿主实验](WINDOWS-STUDY.md#用-windows-opt-加载插件) [FreeBSD 独立 LLVM 宿主](WINDOWS-STUDY.md#freebsd-amd64-独立-llvm-宿主)

ELF 链接方式把 bitcode 交给 LLD，插件在链接器进程中运行。下面的命令用官方 ARMv7 rustc 和匹配的 ARM32 LLD 验证过。独立 LLD 处理 RSSN Lerch bitcode 时也命中了四条单位节点链，链接后的输出与基线逐位一致。

```bash
plugin="$PWD/bazel-bin/llvm_optimizer/RssnLLVMOptimizer.so"
RSSN_SPILL_MODE=off RSSN_LIBM_UNIT_LOG=transform \
  rustc -O -Clinker-plugin-lto -Clink-arg=-fuse-ld=lld \
  "-Clink-arg=-Wl,--load-pass-plugin=$plugin" -Clink-arg=-Wl,--lto-O3 main.rs
```

以一个最终链接的 RSSN release 目标为例。编码后的 Cargo rustflags 保留仓库原有的 `target-cpu=native`。指定 rustc 时，在 CMake 和 Cargo 中使用同一个可执行文件。

```bash
plugin="$PWD/bazel-bin/llvm_optimizer/RssnLLVMOptimizer.so"
flags=$(printf '%s\037%s' '-Ctarget-cpu=native' "-Zllvm-plugins=$plugin")

RSSN_SPILL_MODE=transform RUSTC_BOOTSTRAP=1 \
  CARGO_ENCODED_RUSTFLAGS="$flags" CARGO_TARGET_DIR="$PWD/target/spill-transform" \
  cargo build --release --example numeric_kernels

RSSN_SPILL_MODE=off RUSTC_BOOTSTRAP=1 \
  CARGO_ENCODED_RUSTFLAGS="$flags" CARGO_TARGET_DIR="$PWD/target/spill-off" \
  cargo build --release --example numeric_kernels
```

- `RSSN_SPILL_MODE=transform` 执行分析和改写。显式加载插件时采用这个模式。
- `RSSN_SPILL_MODE=report` 仅分析。rustc 加上 `-Cremark=rssn-spill-pressure` 和 `-Cllvm-args=-pass-remarks-analysis=rssn-spill-pressure` 可查看逐函数诊断。`-Cdebuginfo=1` 可补充源码位置。仅在比较的两组都使用相同调试级别时，才应比较它们的机器码和运行时间。
- `RSSN_SPILL_MODE=off` 禁用自动运行，适合在保持同一插件装载路径的情况下做基线。`RSSN_SPILL_FILTER` 是可选的 LLVM 符号名子串，两种运行模式都支持。Cargo 不会根据插件读取的这些环境变量切换缓存。切换模式、过滤器或插件版本时使用独立的目标目录或重新构建。
- `RSSN_SPILL_PLAN_WIDTH=1` 只允许单成员移动，默认值为 4。用同一个插件可以对联合规划做消融。合并和融合读取、跨块重物化仍按各自的合法性和收益模型运行。
- `opt` 可以显式运行 `function(rssn-spill-pressure)` 或 `function(rssn-spill-pressure-report)`。`-pass-remarks-analysis=rssn-spill-pressure` 输出诊断，`-verify-each` 检查改写后的 IR。Win64 报告按保存寄存器容量区分 GP、XMM 大小和原生宽 SIMD 值，分别定位显式调用与可能降为 libm 的 intrinsic，并单独报告标量 ABI 栈参数。 [Win64 跨调用诊断](WINDOWS-STUDY.md#win64-跨调用压力诊断) [栈参数与真实内核 A/B](ABI-STUDY.md#win64-栈参数与跨调用保存)
- `RSSN_SPILL_LIBM=transform` 启用统一条件下的向量 `exp`、`log` 分支恢复，默认关闭。`RSSN_SPILL_LIBM=report` 另外在 SysV 上报告标量 libm 调用附近的活跃区间移动机会。可用 `function(rssn-libm-pressure)` 和 `function(rssn-libm-pressure-report)` 单独运行。Windows DLL 取证、实际 RSSN 内核 A/B 和复现命令见 [WINDOWS-STUDY.md](WINDOWS-STUDY.md)。
- `RSSN_SPILL_LIBM_REUSE=transform` 为同一循环内的镜像 `exp` 及其常量运算添加按位校验的复用分支，默认关闭。`RSSN_SPILL_LIBM_REUSE=report` 只输出候选。显式入口为 `function(rssn-libm-reuse)` 和 `function(rssn-libm-reuse-report)`，支持 Linux x86-64 与 Win64 MSVC。该选项可独立于向量 libm 和其他压力优化运行。 [Appell F₁ 实验](WINDOWS-STUDY.md#appell-f₁-正负节点的标量复用)
- `RSSN_LIBM_UNIT_LOG=transform` 在 exp-sinh 的单位节点复用 `exp(C*log(x)-x)` 的运行时结果，默认关闭。已由前驱条件排除单位节点的链会跳过。`RSSN_LIBM_UNIT_LOG=report` 仅报告候选。显式入口为 `function(rssn-libm-unit-log)` 和 `function(rssn-libm-unit-log-report)`。 [小参数 Lerch 实验](WINDOWS-STUDY.md#复用-exp-sinh-单位节点的-logexp-结果)
- `RSSN_LIBM_TAIL_CUTOFF=transform` 为 exp-sinh 远端节点提前证明截断或零分子，支持带循环不变线性尺度和分数偏移的表达式，默认关闭。`RSSN_LIBM_TAIL_CUTOFF=report` 仅报告候选。显式入口为 `function(rssn-libm-tail-cutoff)` 和 `function(rssn-libm-tail-cutoff-report)`，支持 Linux x86-64 与 Win64 MSVC。 [Lerch 远端节点实验](WINDOWS-STUDY.md#提前截断-exp-sinh-远端节点) [负阶 Lerch 热循环](WINDOWS-STUDY.md#负阶-lerch-的偏移量与跨调用重读)
- `RSSN_POW_EXPONENT_PACK=transform` 为 Win64 MSVC 的相邻 `pow` 打包重复使用的循环不变指数，默认关闭。窗口可跨到同一循环内的单前驱后继块。`RSSN_POW_EXPONENT_PACK=report` 只输出候选。显式入口为 `function(rssn-pow-exponent-pack)` 和 `function(rssn-pow-exponent-pack-report)`。要求目标启用 AVX，Linux 可交叉构建。 [Win64 指数打包实验](WINDOWS-STUDY.md#win64-跨-pow-调用保留指数对)
- `RSSN_LOOP_EXIT_COMPRESSION=transform` 提前计算循环退出处的纯浮点表达式，在调用点用结果替代多个循环不变输入，默认关闭。`RSSN_LOOP_EXIT_COMPRESSION=report` 只输出候选。显式入口为 `function(rssn-loop-exit-compression)` 和 `function(rssn-loop-exit-compression-report)`，支持 x86-64 SysV 与 Win64 MSVC。 [Appell 系数实验](WINDOWS-STUDY.md#appell-退出表达式与跨调用活跃区间)
- `RSSN_LOOP_PREHEADER_PRESSURE=transform` 将循环前置块里的分支专用纯整数、地址和向量构造移到使用分支。它按对应 ABI 的保存预算比较调用点活跃区间，并限制目的块的估计执行频率。支持 x86-64 Win64、Linux SysV，以及 Linux、macOS、Windows 和 FreeBSD 上使用 64 位指针的普通 AArch64 调用。AArch64 模型区分可装入 `D8–D15` 的值和需要完整向量寄存器的值。默认关闭。`RSSN_LOOP_PREHEADER_PRESSURE=report` 只报告候选。显式入口为 `function(rssn-loop-preheader-pressure)` 和 `function(rssn-loop-preheader-pressure-report)`。 [AArch64 ABI 与真实 CPU lane 对照](WINDOWS-STUDY.md#aarch64-的低-64-位保存预算)
- 计算前置块的压力变化前，先在基线中模拟 CodeGenPrepare 对单用途整数分支比较的下沉。比较的操作数也按新位置计入活跃度。IR 改写只在移动其他定义会破坏原有支配关系时，才一并移动这类比较。 [分支比较消融与评分校准](WINDOWS-STUDY.md#分支比较消融与评分校准)
- `RSSN_FP_ENDPOINT_GUARDS=transform` 合并双指数求积中的端点与权重检查，并在证明安全时跳过无效节点的 `llvm.exp`，默认关闭。`RSSN_FP_ENDPOINT_GUARDS=report` 只输出候选。显式入口为 `function(rssn-fp-endpoint-guards)` 和 `function(rssn-fp-endpoint-guards-report)`。 [Appell 与 Lerch 端点实验](WINDOWS-STUDY.md#双指数求积的端点检查)
- `RSSN_FP_CLASS_INTERVAL=transform` 将循环中控制分支的 binary64 有限非零检查降低为一次无符号区间判断，默认关闭。`RSSN_FP_CLASS_INTERVAL=report` 只输出候选。显式入口为 `function(rssn-fp-class-interval)` 和 `function(rssn-fp-class-interval-report)`，支持 x86-64 SysV 与 Win64 MSVC。 [带状 LU 与求积对照](WINDOWS-STUDY.md#合并有限非零类别的区间检查)
- `RSSN_SCALAR_LOAD_PACK=transform` 合并高 GP 压力循环中的相邻窄整数读取，默认关闭。`RSSN_SCALAR_LOAD_PACK=report` 只输出候选。显式入口为 `function(rssn-scalar-load-pack)` 和 `function(rssn-scalar-load-pack-report)`。x86-64 Linux 使用 4 字节组。Win64 的常规组宽是 8 字节，能证明元素步长为 4 字节时也考虑合并整条记录。若字段值参与后续记录的寻址，额外检查各个字段读取点的平均 GP 压力。C 依赖可在各自的 Clang 23 构建动作中启用。 [zstd IR、机器码和 RSSN 数据 A/B](ZSTD-STUDY.md)
- `RSSN_SWITCH_LOAD_REUSE=transform` 合并多分支都执行的宽读取，并复用其低位字节替代同址窄读取。默认关闭，`report` 仅输出候选。显式入口为 `function(rssn-switch-load-reuse)` 和 `function(rssn-switch-load-reuse-report)`。 [zstd 双哈希实验](ZSTD-STUDY.md#doublefast-短哈希与长哈希共用读取)
- `RSSN_FP_REGISTER_BANK=transform` 在 Win64 MSVC AVX 目标上用通用寄存器携带部分跨 libm 调用的双精度标量，默认关闭。`RSSN_FP_REGISTER_BANK=report` 只输出候选。显式入口为 `function(rssn-fp-register-bank)` 和 `function(rssn-fp-register-bank-report)`。 [Appell 与 Lerch 的跨寄存器组实验](WINDOWS-STUDY.md#用非易失通用寄存器保存浮点位模式)
- `RSSN_POW_UNIT_BASE=transform` 在双精度节点满足舍入边界时跳过一组底数已为 `+1.0` 的 `llvm.pow`，默认关闭。`RSSN_POW_UNIT_BASE=report` 只输出候选。显式入口为 `function(rssn-pow-unit-base)` 和 `function(rssn-pow-unit-base-report)`。RSSN 的 Appell 积分在 Rust 源码中处理这个特例。 [Appell 的源码与 pass 对照](WINDOWS-STUDY.md#appell-单位底数的源码对照)
- `llc -run-pass=rssn-edge-restore` 在 x86-64 寄存器分配后合并多条退出边的相同恢复序列，也可将恢复延后至仅有一个前驱的后继块。`rssn-edge-restore-report` 只报告。两者通过 `-load-pass-plugin` 显式加载，输入为 `virtregrewriter` 之后、栈帧下标消除之前的 MIR。常规 rustc 入口继续运行 IR 管线。 [MIR 布局、代价和运行对照](WINDOWS-STUDY.md#退出边的布局和逐边代价)
- `llc -run-pass=rssn-machine-stack-traffic -pass-remarks-analysis=rssn-machine-stack-traffic` 使用同一阶段的 MIR 输出静态栈流量。报告读取机器内存操作数和帧对象，计入折入其他指令的 spill 重读，并区分 spill、不可变固定栈槽和常量池。`after calls` 统计同一块中至少一条调用之后的访问。`RSSN_SPILL_FILTER` 可选择函数。 [Tape 的 IR 和机器码对照](WINDOWS-STUDY.md#tape-的循环头去重与物理重读)
- 设置 `RSSN_MACHINE_STACK_DETAIL=1` 可沿机器 CFG 追踪最近的调用点。只有所有传播到读取位置的路径都对应同一个静态调用点时，才把读取记到该调用。分歧的前驱、异常入口和含调用的指令束保留为归属不明。每个调用报告 spill 重读宽度、已知字节数和跨块计数。该诊断统计控制流先后关系，无法单独证明重读由这个调用引起。 [Appell、Tape 与 CPU lane 的逐调用对照](WINDOWS-STUDY.md#逐调用核对跨块重读)
- 同一明细模式在 x86-64 报告由常量池指令产生的 spill 槽。它要求槽只有一条全宽写入指令、写入支配每次全宽读取，并在同一块中追踪写入源寄存器的最后定义。报告列出常量来源并单独统计直接重读，供比较栈读取与常量重物化的代价。 [常量重物化对照](WINDOWS-STUDY.md#常量池值的落栈与重物化)

插件使用优化器尾部扩展点处理非 LTO 和 ThinLTO 后链接，使用 Full LTO 的尾部扩展点处理 fat LTO。在 LTO 预链接阶段跳过改写，避免后续内联、向量化和公共子表达式消除放置决策。Bazel 的 FFI 发布目标默认使用全特性、fat LTO 和此宿主插件。插件内容及启用的模式共同参与 Rust 构建缓存键。 [rustc 的 PassBuilder 接入](https://github.com/rust-lang/rust/blob/0ef1bc386059596dada8fb65b50425aab2ea637a/compiler/rustc_llvm/llvm-wrapper/PassWrapper.cpp)

上述独立开关以仅装载插件时的默认值为准。Bazel 发布配置通过 `bazel/runner.py` 启用已验证机器码变化的变换，各 pass 仍按目标和 IR 合法性选择候选。Linux x86-64 启用 ABI 改写和两项 libm 变换。Win64 MSVC 启用两项 libm 变换、指数打包和跨寄存器组保存。两个 MIR 入口使用独立的寄存器分配后工具流程。

已经完成 LLVM 优化的模块可以用 `opt -load-pass-plugin="$plugin" -passes='rssn-release,verify'` 显式运行同组 IR pass。此入口与 rustc 尾部扩展点读取同一组模式、过滤器和参数，并使用同一顺序。输入应是优化完成的 IR 或 fat LTO bitcode。Windows MSVC、Linux、macOS 和 FreeBSD 可用 `--config=release-max --define=rssn_llvm_standalone=1` 让 Bazel 的 FFI 目标在独立 LLVM 工具中运行插件。完整 FFI 交叉编译需要目标 Rust 和系统工具链。macOS x64/ARM64 和 Windows ARM32 的无系统库发布探针覆盖了目标链接与 Cargo 缓存。 [独立 LLVM 链接验证](WINDOWS-STUDY.md#用独立-llvm-工具处理-fat-lto-bitcode)

## ABI 诊断和改写

`RSSN_SPILL_ABI=report` 在 SysV 报告寄存器与栈参数候选，并在 SysV 和 Win64 检测本地被调函数中的间接标量字段读取。`RSSN_SPILL_ABI=transform` 对通过跨过程收益模型的 SysV 本地调用改用 `x86_regcallcc`，默认为 `off`。两者独立于函数 pressure pass，LLVM pipeline 名分别为 `rssn-call-abi-report` 和 `rssn-call-abi`。Cargo A/B 也需按 ABI 模式隔离目标目录。 [实际 RSSN 反汇编及 ABI 分析](ABI-STUDY.md)

## SSA 活跃度与目标代价

反汇编里的栈读写可能来自源码对象、调用约定、保存寄存器或寄存器分配。分配器创建的 spill slot 到机器指令阶段才存在。IR pass 通过改变同时活跃的值、数据表示和定义位置，改变后端的分配输入。测量时要区分 IR 超限代理量、实际机器码中的 spill/reload，以及运行时间。[LLVM Greedy 分配器](https://github.com/rust-lang/llvm-project/blob/0a5579eb649d3634a9c30ac873300cc2944fa070/llvm/lib/CodeGen/RegAllocGreedy.cpp) · [InlineSpiller](https://github.com/rust-lang/llvm-project/blob/0a5579eb649d3634a9c30ac873300cc2944fa070/llvm/lib/CodeGen/InlineSpiller.cpp)

- 对被建模的 SSA 值，PHI 定义在所属块内杀死，操作数使用计入对应入边的前驱活出。块内反向求 `LiveBefore = Use ∪ (LiveAfter − Def)`。SysV 及同样使用 SysV 分配路径的 x86-64 `fastcc` 调用处按 `LiveAfter − Def` 统计跨调用集合，独立于调用参数和返回值。
- GP 和 SIMD 共享各自物理寄存器预算。类型占用来自 TTI 的 legalization 后寄存器数，同一个 `<4 x double>` 在 AVX 下计 1 个寄存器，在只有 SSE2 的目标上计 2 个。标量 FP 与固定向量共同占用 SIMD bank。后端没有提供可靠的通用 IR mask 到 x86 物理 bank 映射，`<N x i1>` 不参与该估计。
- 普通 SysV C ABI 调用点保守计 5 个可保留 GP 和 0 个 SIMD。每个 SIMD 值的 spill 代价依据该类型的 TTI store/load 成本加 16 字节粒度传输项。同占一个寄存器的一个标量和一个宽向量因此可以具有不同的 spill 代价。
- SysV 压力报告把跨调用的 GP 与 SIMD 占用分别按静态站点求和。GP 超限也在各站点单独计算后求和。这样可以看到延后 load 时新增的跨调用地址及其是否用尽保存预算。实际的寄存器保存与重读继续通过 MIR 核对。 [当前 Appell 与 Lerch 发布配置对照](WINDOWS-STUDY.md#当前源码的发布配置消融)
- 每个普通程序点统计接近预算的压力和超限量。超限代价为线性项加有界二次项，几个值同时挤过同一预算时边际成本上升。调用点用独立的估算代价并提高权重。候选计划需降低总分，并逐点限制各寄存器 bank 原有超限量。两地址指令、AVX mask、隐式 libcall、后端调度和实际 fold 仍可能改变物理分配。

## 联合移动与跨块重物化

局部搜索枚举定义移动以及与单用 GEP 一起移动两种方案。候选插入点考虑用 TTI 指令延迟和吞吐估算的依赖暴露，并在移动纯定义时保留较晚的插入点供比较。移动普通浮点加、减、乘只调整定义位置，保留原运算及其 flags。

对候选计划 `S`，将每个结果去掉的旧区间和每个输入新增的区间写成差分事件。一个输入由多个成员共享时，新活跃区间截止到最晚使用点，延长量只记一次。在调用点单独扫跨调用区间。多个成员插进同一间隙时按实际插入顺序检查瞬时 GEP 和结果。计划内的依赖与插入锚点冲突会让组合退出搜索。

候选按跨越的 ABI 调用或局部区域分桶，每桶最多 12 个候选，保留宽度 8 的束搜索，最多联合移动 4 个成员。搜索可以扩展暂时超限的部分组合，只有完整且逐点可行的计划才能落地。例如已有 5 个 GP 值跨越调用时，单独延后一个载入会引入第 6 个 GP 地址。同一个计划可以移走一个整数计算，并让两次载入共享唯一延长的基地址，释放两个 SIMD 值而没有增加跨调用的 GP 数。

跨块重物化处理至多三个同一循环层级的目标块，要求输入在插入点已有独立活跃区间。块频率限制合计动态执行次数，目标代码体积和吞吐代价限制复制，多目标复制前还要证明输入不含 undef 或 poison。CFG 分析按函数有界求解，局部窗口、候选数和改写次数也有上限。这里的搜索保留多变量相互作用，但无法覆盖任意 CFG 上的全局最优调度与分配。[Rematerialization](https://doi.org/10.1145/143095.143143) · [Register Spilling and Live-Range Splitting for SSA-Form Programs](https://publikationen.bibliothek.kit.edu/1000017591) · [Combinatorial Register Allocation and Instruction Scheduling](https://arxiv.org/abs/1804.02452)

## 三种读取合并

| 方案 | 合法性与收益 |
| --- | --- |
| 消除重复地址读取 | 将同一基础指针的 inbounds 常偏移正规化，能识别经不同 GEP 计算的同地址。比较类型和非调试元数据，验证中间无有效写入。收益扣除沿间隙延长的压力，并计入两个副本原先同时活跃时的压力下降。 |
| 融合相邻读取供现有向量使用 | 识别完整标量 insert 链或两个向量半部的顺序 concat。只接受无间隙、同一基础地址的 128/256 位合法目标类型。比较窄 load、装配与宽 load 的 TTI 吞吐、代码体积和依赖延迟。宽 load 写在旧向量定义点，检查延长基地址的压力。 |
| 给必须跨调用的标量建立向量载体 | 当至少两个相邻读无法越过可能修改输入的调用，在原读取区域内用一个 128 位向量保存各 lane。调用后，具有相同操作数位置、共享标量操作数及相同 FP flags 的加、减或乘合成向量操作，在旧使用点投影结果。成本含广播、向量算术、投影、瞬时压力与旧标量 folded reload 的机会成本。只有独立 extract 的载体实测会被 SelectionDAG 拆回标量形式，因此这一分支要求真正的向量消费者。 |

合并、融合和 load 移动都要求非 volatile、非原子、精确的访存范围。AA 的无修改结论之外，还禁止跨越 fence、对象生命周期和相关 intrinsic、inline asm。普通调用要通过 `nounwind`、`nofree`、`nosync` 和 ModRef 检查，或者由私有栈对象证明判定该调用无法接触其存储。融合只复用最低地址 load 的对齐保证，各 lane 都有 `!noundef` 时才把它附到宽 load。TBAA、alias scope 等不会从某一 lane 直接传播到整块宽访问。[LLVM IR 内存和函数属性](https://llvm.org/docs/LangRef.html)

测试包含相邻与非相邻地址、别名写入、同步与对象生命周期、不同目标 legalization、共同输入延长、联合搜索与单候选消融，以及实际 X86 后端的 spill/reload 数。

## 复现实验

```bash
python3 tools/llvm-optimizer/bench/run.py "$plugin" \
  --opt opt-23 --llc llc-23 --cpu 3
```

`--cpu` 选择当前 affinity mask 中的一个 CPU，可省略。脚本在独立目录编译两组，使用相同输入轮流运行 A/B 和 B/A，对每组输出计算覆盖整个数组的位级哈希，并保存每对耗时、区间、汇编和编译诊断。RustFFT 编译使用锁定的 crate 6.4.1 与 rustc fat LTO。两个回调内核都从同一份 O3 IR 生成两组机器码，唯一的 IR A/B 变量是这个 pass。一个回调可安全越过，另一个实际修改输入，必须保留调用之前的值。汇编统计读取 LLVM 的 `Spill` 和 `Reload` 注释并包含 folded 形式。RustFFT 的统计覆盖最后二进制中带 `7rustfft` 符号片段的整个函数集合，包含没有在定时窗口里运行的函数。

测量结果见 [BENCHMARKS.md](BENCHMARKS.md)。

# 构建 Rust 和 C 产物

[English](BUILDING.md)

发布构建、测试和生成操作的入口见 [Bazel 构建指南](BAZEL.zh.md)。

普通的 `cargo build` 和下游 Rust 构建使用 Rust `rlib`。`cargo check` 生成用于类型检查的 Rust 元数据。在 Linux 上，可用 Rust 工具链运行以下命令，刷新 cbindgen 头文件并构建 C 静态库和共享库。

```sh
bazel run //:refresh_headers
bazel build //:ffi_distribution
cc -I. examples/c/smoke.c bazel-bin/ffi_distribution/ffi/librssn.a -lm -lpthread -ldl -o smoke
./smoke
```

默认 FFI 包含 Cranelift JIT、NumPy `.npy` 读写和 wgpu。发布构建采用通用 CPU 参数，FFT 计划仍可在运行时选择 CPU 指令。构建启用栈展开，受保护的 C 入口可将 Rust panic 转为 `RssnStatus::Panic`。C 链接参数因平台而异。

先为宿主配置匹配的 LLVM SDK，再运行 `bazel run //:package_ffi`、`bash scripts/build_ffi.sh --package`，或在 Windows 运行 `scripts/build_ffi.ps1 -Package`。归档包含库、头文件、示例与双语指南，并附有 SHA-256 文件。正式发布用 `bazel run --config=release-max //:release_ffi`。它启用全特性、fat LTO 和 LLVM 插件，并打包 C++、Python、Julia、R 绑定及验证归档中的客户端。设置 `RSSN_FFI_PROFILE=dev` 可通过 Bazel 构建本地开发版 C 库并复制到 `target/debug`。`--llvm-overlay` 与 `-LlvmOverlay` 保留归档名的 `-llvm` 后缀。

独立头文件目标通过 cbindgen 生成头文件，只需编译生成工具。Cargo 的 `DEV=1` 与 `dev-headers` 组合仍可使用。

`signatures.overlay.json` 和 `signatures.witnessed.json` 都是必需的源码文件，源码归档与 Cargo 包也必须包含它们。构建会原样复制 overlay，缺少文件直接报错，签名无效则在构建内核清单时报错。手写覆盖应与使用它的代码一起更新。见证表根据实际调用生成。修改算子后，用 `RSSN_SORT_TRACE=sorts.tsv cargo test --lib` 收集新的测试调用，再运行 `cargo run --profile tools --bin signature-mine -- sorts.tsv` 更新。无见证的算子仍具有其直接声明、覆盖或推导的签名。[签名生成规则](KERNEL.zh.md#3a-见证签名)

构建元数据支持 `SOURCE_DATE_EPOCH`、`VERGEN_IDEMPOTENT` 和显式 `VERGEN_*` 覆盖。Linux 可从 `/etc/os-release`、`/proc` 与 `uname` 查询主机信息。相应来源或 Git 不可用时，缺失值显示为 `unknown`。CPU 数取自 Rust 可用并行度，可能反映容器配额。跨宿主复现元数据时同时设置 `VERGEN_IDEMPOTENT=1` 与 `SOURCE_DATE_EPOCH`，并按需要为归档来源设置 `VERGEN_*`。

## 数值基准

`polynomials` 基准组测量两种拉格朗日插值输入布局，并用相同采样点比较
逐点霍纳求值和批量内核。两个版本使用同一份基准源码、Cargo 构建目录、
功能选项、Rust 工具链和工作线程数。先在旧版运行第一条命令保存基线，
再在新版运行第二条命令。

```sh
RAYON_NUM_THREADS=4 cargo bench --locked --bench engine -- polynomials --save-baseline before
RAYON_NUM_THREADS=4 cargo bench --locked --bench engine -- polynomials --baseline before
```

使用 `fixed_prefix_scan` 基准组可比较单块扫描与每块 256 个元素的扫描。
两轮运行使用相同的工作线程数和基线选项。

`numerical_codec_dispatch` 测量同一会话中的重复数值调用。
`mesh_triangle_areas` 用相同的网格单元分别测量工作线程上限为 1 和 4 时的耗时，
计时范围不包含线程池初始化。对比时使用相同的基准源码和执行配置。

`special_combinatorics` 测量均衡和偏斜输入的对数二项式系数、二项式溢出路径和多项式分组。
`beta_quantiles` 分别测量内部概率、上尾概率和有解析解的极小尾概率。对比时沿用上面的基线方法。

`finite_field_ntt` 基准组测量同一素数的原根冷启动和缓存查询，以及两个输入长度下
已缓存变换计划的多项式卷积。冷启动在每次计时前清空当前进程的复用池。对比时固定工作线程数。

`kalman_history` 基准组用同一模型，在三种观测数量下测量前向滤波及平滑均值。
对比版本时使用相同的基准源码和工作线程数。

## 发布频道

可复用的 CI 工作流用启用全部特性的 Linux 动态库核对生成绑定，并让 Python、Julia、R 和 C++ 使用同一份动态库执行联验。版本发布和每日构建都依赖这项检查。

推送和 `Cargo.toml` 一致的版本标签会启动发布工作流。工作流测试标签对应的提交，并将带校验和的 C 归档上传至草稿。发布草稿后，GitHub 锁定标签与附件并生成发布证明。链接前核验归档的版本与 SHA-256。有已登录的 GitHub CLI 时，用 `gh release verify-asset TAG ARCHIVE --repo Apich-Organization/rssn` 验证下载的归档。

每日构建使用仓库默认分支和 UTC 标签 `nightly-YYYY-MM-DD`。每次构建发布为独立的不可变 GitHub 预发布。crates.io 任务仅适用于版本发布，并需要环境批准。启用 **Release** 工作流的 `dry_run` 可运行验证和构建，且不会发布。

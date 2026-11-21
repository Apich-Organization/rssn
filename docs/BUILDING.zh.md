# 构建 Rust 和 C 产物

[English](BUILDING.md)

普通的 `cargo build` 和下游 Rust 构建使用 Rust `rlib`。`cargo check` 生成用于类型检查的 Rust 元数据。需要 C 静态库或共享库时，使用独立命令构建。

在 Linux 上使用以下命令构建 C 库并重新生成 cbindgen 头文件。

```sh
DEV=1 bash scripts/build_ffi.sh --features dev-headers
cc -I. examples/c/smoke.c target/ffi/librssn.a -lm -lpthread -ldl -o smoke
./smoke
```

脚本调用 `cargo rustc --lib --crate-type staticlib,cdylib`。产物保存在对应的目标和 profile 目录中。可设置 `CARGO_TARGET_DIR`，也可传入 `--features jit`、`--target`、`--locked` 或 `--offline` 等 Cargo 参数。默认 `ffi` 配置继承发布优化并启用栈展开，受保护的 C 入口可将触发栈展开的 Rust panic 转为 `RssnStatus::Panic`。若需开发构建，设置 `RSSN_FFI_PROFILE=dev`，产物目录为 `debug`。C 链接参数因平台而异。

运行 `bash scripts/build_ffi.sh --package`，或在 Windows 运行 `scripts/build_ffi.ps1 -Package`，将本机库、头文件、示例、双语指南及校验值一起打包。

RSSN 源码仓库将 `target-cpu=native` 设为通用 Cargo 构建默认值。对于 `x86_64-pc-windows-msvc`，Cargo 选用该目标专用的栈大小链接选项。发布工作流为发行归档设置 `RUSTFLAGS="-C target-cpu=generic"`。自行构建供分发的归档或为其他机器构建时，可使用同一设置。FFT 计划仍可在运行时选择 CPU 指令。

使用已有头文件时可省略 `DEV=1`。重新生成头文件时，须同时设置 `DEV=1` 并启用 `dev-headers`，以免普通 Rust 构建引入 cbindgen。单独启用 Cargo 特性不会改变产物类型。从 `cargo build --release` 获得 C 库的既有宿主请改用上述脚本。

`signatures.overlay.json` 和 `signatures.witnessed.json` 都是必需的源码文件，源码归档与 Cargo 包也必须包含它们。构建会原样复制 overlay，缺少文件直接报错，签名无效则在构建内核清单时报错。手写覆盖应与使用它的代码一起更新。见证表根据实际调用生成。修改算子后，用 `RSSN_SORT_TRACE=sorts.tsv cargo test --lib` 收集新的测试调用，再运行 `cargo run --profile tools --bin signature-mine -- sorts.tsv` 更新。无见证的算子仍具有其直接声明、覆盖或推导的签名。[签名生成规则](KERNEL.zh.md#3a-见证签名)

构建元数据支持 `SOURCE_DATE_EPOCH`、`VERGEN_IDEMPOTENT` 和显式 `VERGEN_*` 覆盖。Linux 可从 `/etc/os-release`、`/proc` 与 `uname` 查询主机信息。相应来源或 Git 不可用时，缺失值显示为 `unknown`。CPU 数取自 Rust 可用并行度，可能反映容器配额。跨宿主复现元数据时同时设置 `VERGEN_IDEMPOTENT=1` 与 `SOURCE_DATE_EPOCH`，并按需要为归档来源设置 `VERGEN_*`。

## 发布频道

推送和 `Cargo.toml` 一致的版本标签会启动发布工作流。工作流测试标签对应的提交，并将带校验和的 C 归档上传至草稿。发布草稿后，GitHub 锁定标签与附件并生成发布证明。链接前核验归档的版本与 SHA-256。有已登录的 GitHub CLI 时，用 `gh release verify-asset TAG ARCHIVE --repo Apich-Organization/rssn` 验证下载的归档。

每日构建使用仓库默认分支和 UTC 标签 `nightly-YYYY-MM-DD`。每次构建发布为独立的不可变 GitHub 预发布。crates.io 任务仅适用于版本发布，并需要环境批准。启用 **Release** 工作流的 `dry_run` 可运行验证和构建，且不会发布。

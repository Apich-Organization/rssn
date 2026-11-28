# 参与 RSSN 贡献

[English](CONTRIBUTING.md) | [简体中文](CONTRIBUTING.zh-CN.md)

感谢你考虑为 **rssn** 做出贡献！无论是修复缺陷、实现新算法、优化性能还是完善文档与测试，我们都非常欢迎。如有新想法，建议先提交 Issue 讨论。

## 代码分布

建议先阅读 [ARCHITECTURE.zh-CN.md](ARCHITECTURE.zh-CN.md) 了解整体架构。各模块位置如下：

| 目标 | 对应目录 |
|---|---|
| 新增或修复符号恒等式、求解器或积分器 | `src/rules/<domain>`（单个 `RuleSet` 的重写、定义与内核） |
| 新增或修复纯数值算法 | `src/kernels`（不涉及符号项，供规则与用户共同调用） |
| 新增或修复仿真、场或粒子求解器 | `src/sim` |
| 修改引擎本身（存储、调度器、提取、事实推导） | `src/graph` |
| 修改量纲或单位 | `src/graph/dim*.rs`（推导）与 `src/dim`（单位表、常数、单位制）；设计见 [docs/DIMENSIONS.zh.md](docs/DIMENSIONS.zh.md) |
| 修改宿主可见表面（清单、签名、动态调用、插件） | `src/kernel` 与 `src/api/value.rs`；契约见 [docs/KERNEL.zh.md](docs/KERNEL.zh.md) |
| 修改 C 接口 | `src/ffi`；头文件生成见 [docs/BUILDING.zh.md](docs/BUILDING.zh.md) |
| 修改解析/输出、LaTeX、Typst 或绘图 | `src/io` |
| 修改编译求值或 JIT | `src/backend` |
| 共享 CPU 调度策略、临时缓冲区或不可变计划 | `src/execution.rs`、`src/kernels/par.rs` 与 `src/kernels/reuse.rs` |
| 生成数值绑定代码 | `tools/rssn-api-gen`（[README](tools/rssn-api-gen/README.zh-CN.md)） |

单元测试紧邻对应源码；跨模块集成测试位于 `tests/`。

## 开始之前

1. **安装 Rust**：本项目使用 Rust 2024 edition，需要 `Cargo.toml` 中 `rust-version` 指定的工具链（当前为 1.100）或更新版本，可通过 [rustup](https://rustup.rs/) 安装。
2. **安装 Bazel 并克隆仓库**：安装 Bazel 或 Bazelisk 后克隆你的 fork：
   ```bash
   git clone https://github.com/YOUR-USERNAME/rssn.git
   cd rssn
   ```
3. **构建与测试**：符号测试套件会运行大量嵌套引擎求值（测试配置已开启 `opt-level = 1`），首次构建耗时较长：
   ```bash
   bazel test //:rust_tests
   ```
   可使用 `bazel test //:rust_tests --test_arg=--test --test_arg=dimensions` 运行单个集成测试目标。`//:rust_all_features_test` 覆盖所有特性，`//:gpu_test` 需要可用的 wgpu 适配器。日常开发可用 `bazel run //:check -- --lib --tests` 快速检查编译。

## 开发流程

1. **从 `main` 创建特性分支**：
   ```bash
   git checkout -b feature/my-new-feature
   ```
2. **编写代码与测试**：新增重写规则会自动由 `rules::tests::no_standard_rewrite_is_refuted` 进行数值验证，同时请补充对应的恒等式测试。
3. **遵守代码规范**：
   - `src/lib.rs` 启用了严格的 Clippy 检查（`pedantic`、`nursery`、`unwrap_used`、`expect_used`、`indexing_slicing`、`arithmetic_side_effects`、`missing_docs` 等）。新代码必须通过检查；局部 `#[allow]` 须附说明注释。
   - 所有公开项须编写文档注释，格式遵循 `rustfmt.toml`。
4. **提交 PR 前运行检查**：
   ```bash
   bazel run //:fmt
   bazel test //:clippy_test //:rust_all_features_test
   ```
   若修改了对应区域，还需运行：
   ```bash
   bazel test //:api_generator_test //:contracts_test //:headers_check_test
   bazel test //:bindings_check_test //:bindings_test
   ```
5. **提交并推送**：使用清晰的提交信息（如 `feat(calculus): add Hermite reduction for algebraic integrands`）。
6. **向 `Apich-Organization/rssn` 的 `main` 分支发起 Pull Request**，说明变更内容与原因。

### 新增算子

对于编解码器支持的纯 Rust 数值函数，请将算法与测试放入对应的 `src/kernels` 或 `src/sim` 模块，沿用 `num_<family>_<operation>` 命名并在 `codec-policy.tsv` 中注册 Rust 路径与量纲契约（详见 [绑定生成文档](tools/rssn-api-gen/README.zh-CN.md)）。

对于符号规则、宿主回调或原生对象，请直接注册：

1. 在规则集安装函数中注册（`Installer::op`）。
2. 绑定量纲规则（`i.dim(op, rule)`），否则 `graph::dim` 覆盖率扫描会失败。
3. 声明 `Signature`（或无法静态类型化时声明 `DynamicReason`），并补充 `.doc(..)`、`.category(..)` 与 `.since(..)`。
4. 运行 `bazel test //:contracts_test` 与算子测试。

`signatures.overlay.json` 提供人工审核后的签名修正，`signatures.witnessed.json` 记录从测试语料提取的签名；若需刷新见证签名，运行 `bazel run //:refresh_witness`（[内核签名来源](docs/KERNEL.zh.md)）。

### 修改 C 接口

`rssn.h` 与 `rssn.hpp` 由 cbindgen 生成并签入仓库。修改 `src/ffi` 后请运行 `bazel run //:refresh_headers` 并一并提交更新后的头文件。

## CI/CD

`.github/workflows` 包含以下工作流：

| 工作流 | 触发条件 | 职责 |
|---|---|---|
| `ci.yml` | 推送至 `main`/`v0.3`、PR、手动触发、每周定时、发布门禁 | Nightly rustfmt、严格 Clippy、Linux/Windows 全特性测试、macOS 默认特性测试、契约与生成检查、C/C++ 客户端、Rust 示例与 rustdoc（Linux 含软件 Vulkan wgpu 计算测试）。必需分支检查为 `CI success`。 |
| `security.yml` | `main`/`v0.3`、PR、每周定时 | CodeQL、PR 依赖审查、RustSec 安全通告检查。 |
| `release.yml` | 版本标签、每日构建、手动触发 | 在发布提交上运行 CI 并构建多平台 C 归档；支持版本发布至 crates.io 与 `dry_run` 模式。 |

## 致谢

所有贡献者均会在发行说明与 GitHub 页面致谢。感谢每一份改进！

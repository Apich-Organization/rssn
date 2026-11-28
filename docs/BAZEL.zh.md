# 用 Bazel 构建 RSSN

[English](BAZEL.md)

如果安装的命令名是 `bazelisk`，可用它替换示例中的 `bazel`。Windows 上请将 Bazel 输出根目录放在 NTFS 卷，例如 `bazelisk --output_user_root=D:/bz build //:headers`。

使用 Bazelisk 启动 Bazel 7.7.1，并安装 Python 3.11 或更新版本及 `Cargo.toml` 指定的 Rust 工具链。Bazel 和 Cargo 共用 `Cargo.lock`。RSSN 模块也向 Aevia 提供源文件、工作流规则和 LLVM 插件。

| 命令 | 产物或检查 |
| --- | --- |
| `bazel build //:rustlib` | 开发用 Rust 库 |
| `bazel run //:check` | 工作区类型检查，省去测试可执行文件的链接 |
| `bazel test //:fmt_test` / `bazel run //:fmt` | 检查或格式化 Rust 源码，使用 nightly rustfmt |
| `bazel test //:docs_test` | 检查文档，再构建示例并执行功能示例 |
| `bazel run //:example -- quickstart` | 运行 Rust 示例，可用 `--profile` 选择配置 |
| `bazel run //:benchmark -- --help` | 查看计时用例和基准参数 |
| `bazel run //:benchmark -- engine -- polynomials` | 运行指定基准并向其传递参数 |
| `bazel run //:benchmark -- engine -- economy_qr` | 单线程测量四种固定尺寸的薄 QR |
| `bazel run //:benchmark_ci -- --help` | 运行 Criterion 冒烟检查或收集全部 RSSN 计时报表 |
| `bazel run //:llvm_benchmark -- --help` | 使用 Bazel 构建的优化器运行 LLVM 回调、lanes、ABI 和 RustFFT 基准 |
| `bazel run //:llvm_benchmark_build -- --help` | 按指定配置和工具链构建 LLVM 基准程序 |
| `bazel run //:llvm_benchmark_compare -- --help` | 以相同输入比较预构建的研究程序 |
| `bazel run //:llvm_c_benchmark_build -- zstd --help` | 构建解码、doubleFast、FSE 或 Win64 Appell 与 Lerch 研究程序 |
| `bazel run //:llvm_benchmark_ci -- --help` | 运行 x86-64 Linux LLVM 基准任务 |
| `bazel test //:test` | 默认 Rust 测试和生成契约检查 |
| `bazel test //:rust_minimal_test` | 关闭默认特性后运行 Rust 测试 |
| `bazel test //:integration_test` | Rust 全特性、契约、头文件、见证重放及四语言绑定测试 |
| `bazel test //:clippy_test` | 对全特性工作区的全部 Rust 目标执行严格 Clippy |
| `bazel test //:ffi_native_clients_test` | 运行原生 C、C++ 客户端测试，Linux 使用静态库，macOS 和 Windows 使用动态库 |
| `bazel test //:rust_all_features_test` | 全特性 Rust 测试，包括构建机可执行的 GPU 测试 |
| `bazel test --define=rssn_bounded_tests=1 //:rust_all_features_test` | 每批运行两个 RSSN 集成测试程序，完成后回收该批的可执行文件 |
| `bazel test //:rust_all_features_test --test_arg=--batches=4` | 每批运行四个集成测试程序并回收其缓存空间，仍执行全量测试 |
| `bazel test //:rust_ignored_test //:witness_test` | 耗时数值回归和已签入见证表的重放 |
| `bazel test //:bindings_check_test //:bindings_test` | 针对同一个全特性 cdylib 检查 API 以及 Python、Julia、R、C++ 绑定 |
| `bazel run //:bindings_from_library -- --library /path/to/librssn.so` | 校验生成的绑定并测试已有 cdylib |
| `bazel test //:generated_sources_check_test` | 检查契约、结构化规则、绑定和编解码器的生成结果 |
| `bazel test //:api_generator_test` | 测试源码模型与生成器，无需链接 RSSN |
| `bazel build //:manifest` | 在 `bazel-bin/manifest/` 生成内核清单 JSON |
| `bazel build //:headers` | 在 `bazel-bin/headers/` 生成 C 和 C++ 头文件 |
| `bazel run //:refresh_headers` | 重新生成检出的 `rssn.h` 和 `rssn.hpp`，无需编译 libRSSN |
| `bazel test //:headers_check_test` | 校验生成的 C 和 C++ 头文件 |
| `bazel run //:refresh_contracts` | 刷新本机接口契约、绑定和编解码器 |
| `bazel run //:refresh_bindings` | 基于新内核清单刷新语言绑定接口 |
| `bazel run //:refresh_bindings_from_manifest -- --manifest manifest.json` | 用已链接的 RSSN 产物清单刷新各语言接口 |
| `bazel run //:refresh_witness` | 采集测试轨迹并更新 `signatures.witnessed.json` |

绑定集成测试使用带 NumPy 的 Python、Julia、带开发头文件的 R、CMake 和 C++ 编译器。R 包安装在构建缓存中，所有语言都使用同一个指定的库文件。加上 `--config=release` 可测试优化后的 cdylib。只需全特性开发版动态库时，可运行 `bazel build //:ffi_dev`。

如果 R 环境提供独立的 C++ 工具链，可用 `--repo_env=CXX=/path/to/c++` 选择与待测库兼容的编译器。

使用 `bazel run //:benchmark -- NAME` 运行计时用例。

| 名称 | 测量内容 |
| --- | --- |
| `backends` | 参考解释器、Tape 和 JIT 的批量执行 |
| `engine-timing` | `tools` 配置下的端到端符号请求 |
| `extract` | 表达式提取和可达等价类遍历 |
| `eval-cost` | 首次及缓存后的求值代价计算 |
| `cp-als` | 秩为 3 的张量分解和单次迭代成本 |
| `groebner-elimination` | 分块消元和字典序 Gröbner 基 |
| `fft-plans` | 合数长度和 2 的幂长度的正向变换 |
| `parallel-kernels` | 数值内核单线程和多线程执行对照 |
| `parallel-sim` | 仿真内核单线程和多线程执行对照 |

计时示例默认使用 `tools` 配置。基于测试的计时用例默认使用 `test` 配置，并逐个执行。
可用 `RAYON_NUM_THREADS` 选择并行侧的线程数，例如
`RAYON_NUM_THREADS=4 bazel run //:benchmark -- parallel-kernels`。

比较两个工作树时，请为两边设置同一个绝对路径 `CRITERION_HOME`。
用 `bazel run //:benchmark -- engine -- economy_qr --save-baseline before` 保存第一版的结果，
再用 `bazel run //:benchmark -- engine -- economy_qr --baseline before` 比较第二版。
两次应使用同一台机器和同一套 Rust 工具链。

常规 Linux CI 以测试模式执行 Criterion。定时任务会收集 Criterion 结果和全部九项计时报表，
手动触发时也可以启用 `benchmark_reports`。Criterion 结果保存为工件，其他计时报表写入任务日志。

启用 `llvm_benchmarks` 可运行回调、lanes、ABI、RustFFT、zstd 解码、doubleFast 和 FSE
的配对用例。Linux CI 未指定 `RSSN_LLVM_CONFIG` 时，由 Bazel 安装固定版本的 LLVM SDK。该任务
还会构建 libm 调用计数探针、校验对照组输出，并上传 JSON 计时报表和构建日志。Win64 研究
及四组 zstd I/O 对照的输入见各自研究文档。

## FFI 分发

```sh
bazel build //:ffi_distribution
bazel run //:package_ffi -- --dist-dir dist
```

`//:ffi_distribution` 在一次 rustc 调用中生成静态库与动态库。产物位于 `bazel-bin/ffi_distribution/ffi/`。默认使用全部 Cargo 特性并加载 LLVM 插件。归档包含库、生成的头文件、示例和文档，并附有 SHA-256 校验文件。

开发时用 `bazel build //:ffi_distribution_dev` 在 `bazel-bin/ffi_distribution_dev/ffi/` 生成两种库，或用 `bazel run //:package_ffi_dev -- --dist-dir dist` 将开发版库和生成的头文件一起打包。默认归档名以 `-dev.tar.gz` 结尾。

分发构建默认使用 O3、栈展开、通用 CPU 参数、fat LTO 和单 CGU。i686 目标使用 rustc 提供的 SSE2 基线。Bazel 依赖 `//:release_pass_gate` 先运行发布版 LLVM 探针，再构建两种库。`//:release_ffi` 打包后解包并运行 C、C++ 和 Python 客户端测试，通过后才把归档写入分发目录。检查已有归档可用 `bazel run //:ffi_package_smoke -- --archive /path/to/rssn.tar.gz`。

`//:package_ffi_from_binaries` 接收 `--native-dir` 和 `--headers-dir`，对配套的已有产物执行相同的打包验证。

### LLVM 优化插件

LLVM 插件需要 CMake、C++ 编译器和与 `rustc -vV` 完整版本一致的宿主 LLVM SDK。Bazel 自动查找对应的 `llvm-config-<主版本>`，也接受显式路径。Windows 构建机通过 `--repo_env=RSSN_LLVM_IMPORT_LIBRARY=/path/to/host-llvm.lib` 提供加载进程的 LLVM C++ 导入库。若官方 rustc 缺少插件所需的 LLVM 导出，可通过 `--repo_env=RUSTC=/absolute/path/to/rustc` 指定兼容编译器，并用 `CARGO` 指定对应的 Cargo。CI 通过按 Rust 三元组索引的 JSON 映射读取运行器上已安装的工具链路径，变量为 `RSSN_RELEASE_LLVM_CONFIGS`、`RSSN_RELEASE_LLVM_IMPORT_LIBRARIES`、`RSSN_RELEASE_RUSTCS` 和 `RSSN_RELEASE_CARGOS`。原有共享路径 `RSSN_RUSTC` 与 `RSSN_CARGO` 仍可用。

Unix 发布 CI 使用 `--repo_env=RSSN_FETCH_LLVM_SDK=1`。LLVM 23.1.5 对应的 Linux x86-64、
Linux Arm64 和 macOS Arm64 SDK 由 Bazel 下载并校验。Intel macOS 和 Windows 需要提供匹配的
SDK 和所需工具链路径。可用 `RSSN_RELEASE_RUNNERS` JSON 映射按 Rust 三元组选择已配置的运行器。
显式的 `RSSN_LLVM_CONFIG` 优先于下载。

```sh
bazel test --repo_env=RSSN_LLVM_CONFIG=/path/to/llvm-config //:llvm_optimizer_test //:release_pass_test
bazel build --config=release --repo_env=RSSN_LLVM_CONFIG=/path/to/llvm-config //:rustlib
bazel build --repo_env=RSSN_LLVM_CONFIG=/path/to/llvm-config //:ffi //:ffi_static
bazel run --repo_env=RSSN_LLVM_CONFIG=/path/to/llvm-config //:release_ffi -- --dist-dir dist
```

`//:ffi` 构建动态库和适用平台的导入库。`//:ffi_static` 构建静态库。`//:package_ffi_portable` 可显式构建无插件的对照包。插件是 Rust action 的声明输入，经 `-Zllvm-plugins` 传给目标 rustc。单独输出的 Rust 库关闭 LTO，使插件直接在其 O3 流水线运行。各 pass 根据目标平台和函数选择是否应用。

为 Windows MSVC、Linux GNU 或 musl、FreeBSD 目标构建 FFI 分发包时，可用 `--define=rssn_llvm_standalone=1` 让 SDK 的 `opt` 在 LTO 后的 bitcode 上运行插件。Windows 构建机使用匹配 SDK 的 LLVM 导入库。Bazel 继续启用全特性、fat LTO 和单 CGU，通过 `llc` 生成替换对象，最后用 Cargo 原有链接参数调用 rustc。ELF 路线生成位置无关对象以构建共享库。给 `//:release_pass_test` 传入同一参数可检查共享库、静态库、适用平台的导入库以及 Cargo 依赖跟踪。

```sh
bazel test --config=release-max --define=rust_target=x86_64-pc-windows-msvc --define=rssn_llvm_standalone=1 //:release_pass_test
bazel run --config=release-max --define=rust_target=x86_64-pc-windows-msvc --define=rssn_llvm_standalone=1 //:release_ffi -- --dist-dir dist
```

### 交叉编译

通过 `--define=rust_target=<三元组>` 选择目标，并安装对应 Rust 标准库、原生链接器、C 工具链和目标 SDK。可选的插件在构建机运行。通过 `--repo_env=RSSN_TARGET_ENV=/absolute/path/to/toolchains.json` 提供按 Rust 三元组索引的环境变量文件，Bazel 将它作为配置输入。以下示例从 Linux 交叉编译 AArch64：

```json
{
  "aarch64-unknown-linux-gnu": {
    "CARGO_TARGET_AARCH64_UNKNOWN_LINUX_GNU_LINKER": "/usr/bin/aarch64-linux-gnu-gcc",
    "CARGO_TARGET_AARCH64_UNKNOWN_LINUX_GNU_RUNNER": "/usr/bin/qemu-aarch64 -L /usr/aarch64-linux-gnu",
    "CC_aarch64_unknown_linux_gnu": "/usr/bin/aarch64-linux-gnu-gcc",
    "AR_aarch64_unknown_linux_gnu": "/usr/bin/aarch64-linux-gnu-ar"
  }
}
```

```sh
bazel build --config=release-max --define=rust_target=aarch64-unknown-linux-gnu \
  --repo_env=RSSN_TARGET_ENV=/absolute/path/to/toolchains.json //:ffi_distribution
```

为 Windows MSVC、macOS、FreeBSD 或其他 Rust 目标使用对应的 Cargo `_LINKER` 和 `_RUNNER` 变量。目标工具链可提供 C/C++ 编译器、归档工具、sysroot、`PKG_CONFIG_*`、Windows 的 `INCLUDE` 和 `LIB` 或 macOS 的 `SDKROOT`。交叉发布探针可以只编译和链接。Rust 测试也可以通过配置的 runner 执行。绑定集成测试应在能原生加载对应库的目标系统运行。本机开发目标无需 LLVM SDK。

## 缓存与计时

Bazel 分别跟踪声明的项目输入、Cargo 清单与锁文件、LLVM 插件。Rust action 在稳定的本地源码镜像中委托 Cargo 调度 rustc。未变化的文件保留修改时间，删除的输入从镜像移除；因此源码编辑触发重新执行时仍能复用 Cargo 指纹、第三方构建结果和开发版 CGU 增量缓存。工具链配置是本地状态键的一部分。插件路径包含库文件和优化设置的摘要，任一项变化都将使依赖它的 Cargo 产物失效。分发构建使用通用 CPU 参数和确定性 RSSN 元数据。

本地状态默认位于 `.bazel-cache/`。Windows 可用 `--repo_env=BAZEL_SCI_CACHE=C:/b/r` 设置较短的绝对路径。可另用 `--disk_cache` 保留 Bazel action 产物。有现成编译器缓存时可传入 `--repo_env=RSSN_SCCACHE=/path/to/sccache`。原位更换工具链或 SDK 后，在计时前运行 `bazel sync --configure`。这些 action 使用构建机的本地工具链。

在同一机器、工具链、输入和目标上比较发布配置。测冷构建时采用独立缓存目录，测热构建时沿用缓存和原命令。`--define=rssn_timings=1` 将 Cargo 报告保存在输出中指示的位置。`--profile=/path/to/profile.gz` 记录 Bazel action 与调度耗时。`--define=rssn_jobs=N` 单独调整 Cargo 和 CMake 并行任务数。`--define=rssn_frontend_threads=4` 可比较 RSSN 的并行前端，包装器仅给 RSSN 的 rustc 调用附加 `-Zthreads`，并使用单独的 Cargo 缓存。

- `--config=release` 使用 ThinLTO 和 16 个 CGU。
- `--config=release-local` 关闭 LTO 并保留 O3。LLVM 目标仍使用优化插件。
- `--config=release-max` 使用 fat LTO 和单 CGU。受控实验可以用 `--define=rssn_lto` 和 `--define=rssn_codegen_units` 单独覆盖。

这些对比适用于链接产物。Bazel 发布配置独立于兼容 Cargo 发布配置。评估工作负载的数值一致性和运行时间后再调整分发配置。

# RSSN 语言绑定：Python、Julia、C++ 与 R

[English](README.md) | [简体中文](README.zh-CN.md)

基于同一套 C ABI 与一致性测试表构建的四种宿主语言完整封装。

| 目录 | 包名 | 入口 |
|---|---|---|
| [`python/`](python/) | `rssn`（ctypes） | `import rssn` |
| [`julia/`](julia/) | `RSSN`（ccall） | `using RSSN` |
| [`cpp/`](cpp/) | `rssn::client`（纯头文件） | `#include "rssn/client.hpp"` |
| [`r/`](r/) | `rssn`（`.Call`） | `library(rssn)` |
| [`api/`](api/) | 生成的 API 描述 | `rssn_api.json` |
| [`tools/`](tools/) | 代码生成器 | `extract_api.py`、`generate.py` |
| [`tests/`](tests/) | 跨语言一致性测试套件 | `cases.tsv` |
| [`examples/`](examples/) | 宿主回调与精确张量示例 | `python/callbacks.py`、`julia/callbacks.jl`、`cpp/callbacks.cpp`、`r/callbacks.R` |

## 覆盖范围

rssn 暴露精简的 C ABI（`rssn.h`，约 50 个入口）与由内核清单描述的完整算子表。各语言绑定完整覆盖以下两层：

- **C ABI 一一映射**：会话与配置生命周期、项构建与四种渲染（文本/LaTeX/Typst/Unicode）、单次与并行批量计算（`compute_all`）、动态值 ABI（`call_value`）、宿主回调（`Value::Function`）、量纲报告、内核清单与仿真调度。所有由 rssn 分配的字符串、字节缓冲区和张量均自动通过 `rssn_*_free` 释放，错误统一映射为携带 `RssnStatus` 的异常。
- **算子表逐算子生成**：`python/rssn/_ops.py`、`julia/src/ops.jl`、`cpp/include/rssn/ops.hpp` 与 `r/R/ops.R` 由内核清单生成，保留参数名与文档；语言关键字做后缀或前缀处理，原始拼写仍可通过 `apply` 与 `OPS` 表访问。`eval = dynamic` 算子支持 checked 形式（Python `checked=True`、Julia `checked=true`、R `checked = TRUE`、C++ `rssn::client::checked`）。
- **清单元数据**：82 个单位（精确有理缩放、仿射偏移、别名）、30 个物理常数与 29 个仿真场景。

## 快速开始

先构建共享库并设置 `RSSN_LIB` 指向库文件或所在目录（源码检出下各绑定会自动探测 `bazel-bin/ffi_distribution/ffi` 与 `bazel-bin/ffi_dev/ffi`）：

```sh
bazel build //:ffi_distribution                 # 全特性 LLVM 分发版
bazel build //:ffi_dev                          # 全特性开发共享库
export RSSN_LIB="$PWD/bazel-bin/ffi_distribution/ffi"
```

各语言最小示例见 [Python](python/README.zh-CN.md)、[Julia](julia/README.zh-CN.md)、[C++](cpp/README.zh-CN.md) 与 [R](r/README.zh-CN.md)。

## 执行后端、NumPy 与宽整数

- **执行后端**：通过 `available_backends`、`current_backend` 与 `set_default_backend` 查询或切换进程默认后端（`auto` 恢复自动选择），或使用 `open_backend` / `release_backend` 配合 `num_backend_compile`。详见 [后端指南](../docs/BACKENDS.zh.md)。
- **NumPy 数组**：原生 NumPy 编解码器支持 C/Fortran 序 `float32`/`float64` 读取与小端 `float64` 写出（通过 `encode_npy`/`decode_npy`，特性标志位于 `manifest.features.npy`）。
- **任意精度整数与精确矩阵**：Python `int` 与 Julia `BigInt` 超出 64 位有符号范围时自动走十进制文本通道；C++ 与 R 使用 `value_integer_decimal`。精确模运算与有理矩阵运算使用 `call_value` 与 `num_native_qmatrix_new`。

## 宿主回调

189 个高阶算子（如 `num_brent`、`num_bfgs`）接受函数参数。各语言绑定将原生函数封装为 C 回调（`rssn_callback_new_for` / `rssn_callback_new` / `rssn_callback_signatures`），自动完成参数与返回值的细胞转换。

## 测试与代码生成

```sh
bazel test //:bindings_check_test //:bindings_test
python bindings/tests/run_conformance.py
```

- 四个语言运行器共享 `tests/cases.tsv`；`run_readme_examples.py` 会直接提取并执行四个语言 README 中的代码块。
- 重新生成绑定：

```sh
bazel build //:manifest
bazel run //:refresh_bindings
bazel test //:bindings_check_test
```

## 已知限制

与 `src/ffi/mod.rs` 一致：所有注册算子均可通过 `apply`、`call_value` 或回调调用；无项形式的原生句柄、记录与回调不能 `to_term`；`bincode-next` 缓冲区对封装层不透明；清单反映标准规则集快照；`ffi` 构建配置下受保护 C 入口捕获 panic 并返回 `Panic` 状态。

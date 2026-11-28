# RSSN C++ 绑定

[English](README.md) | [简体中文](README.zh-CN.md)

基于 rssn C ABI 的纯头文件 RAII 封装，并为每个注册算子提供对应函数。仅需 C++17 标准库与 `rssn.h`。

```cpp
#include <iostream>
#include "rssn/client.hpp"
#include "rssn/ops.hpp"
#include "rssn/value.hpp"

int main() {
    rssn::client::Session session;
    const auto x = session.sym("x");

    const auto answer = session.compute(diff(sin(x) * cos(x), x));
    std::cout << answer.text << "\n";                     // cos(2*x)

    rssn::client::Config numeric;
    numeric.numeric(1e-12).bind("x", 0.0);
    std::cout << *session.compute(diff(sin(x) * cos(x), x), &numeric).value << "\n";
}
```

## 构建

需要包含 `rssn.h`（仓库根目录）与 C 动态/静态库：

```sh
bazel build //:ffi_dev
bazel test //:ffi_native_clients_test
```

CMake 支持通过 `RSSN_LIB` 传入共享库路径或目录，或通过 `-DRSSN_LIBRARY=/path/to/lib` 显式指定链接库（Linux/macOS/Windows 分别对应 `librssn.so`/`librssn.a`、`librssn.dylib`/`librssn.a`、`rssn.dll.lib`+`rssn.dll`/`rssn.lib`；MinGW 可用 `gendef` + `dlltool` 生成导入库）。

## API 概览

| 类型 | 职责 |
|---|---|
| `rssn::client::Session` | 拥有所有项；仅可移动；单线程 |
| `rssn::client::Term` | 轻量句柄；`.text()`、`.latex()`、`.typst()`、`.pretty()`、`.value()`、`.tensor()` |
| `rssn::client::Config` | 目标阶段、变量绑定、假设、预算与量纲模式（链式调用） |
| `rssn::client::Answer` | `.term`、`.text`、`.value`、`.error`、`.reduced`、`.iterations`、`.dimension`、`.diagnostic_codes`、`.latex()`、`.tensor()` |
| `rssn::client::diff(f, x)` | 生成的算子函数，返回**未求值**请求项 |
| `rssn::client::kOps` / `op_spec(name)` | 算子元数据表（名称、元数、`ffi` 可达性） |
| `rssn::client::Value` | 动态值（`rssn/value.hpp`）：`kind()`、`as_f64()`、`text()`、`json()`、`to_term()`、`tensor()`、`size()`、`at()` |
| `rssn::client::kUnits` / `kConstants` / `kScenarios` | 单位表、物理常数表与仿真场景目录 |
| `rssn::client::Session::run_scenario(name, params)` | 按名称运行仿真 |
| `rssn::client::Error` | 携带 `.status()` 的异常类 |

## 常见用法

- **符号与数值求值**：`Term` 重载算术运算符，数值与字符串可隐式转换为默认会话中的 `Term`；`session.compute("...")` 在当前会话解析。`session.compute_all({...})` 并行求值多个请求。
- **张量与 NumPy**：`answer.tensor(&shape)` 返回行主序 `double` 数组；`session.encode_npy` / `session.decode_npy` 读写 `.npy` 字节流。
- **宽整数、后端与动态算子**：超 64 位整数使用 `value_integer_decimal` 与 `Value::text()`；后端通过 `available_backends()`、`current_backend()`、`set_default_backend()` 及 `open_backend`/`release_backend` 管理；`checked::num_*` 以原生参数直接调用动态算子；`with_coherent_dimension` 附带 `DimExpr` JSON 量纲。
- **宿主回调与注意事项**：`Callback::for_argument` 封装 `std::function` 且生命周期必须覆盖调用；C++ 关键字加 `_op` 后缀（`and_op`、`not_op` 等）；减法/取负/除法对应引擎 `sub`/`neg`/`divide`。

## 测试

```sh
ctest --test-dir build/cpp --output-on-failure
```

包含 `tests/smoke.cpp`、`tests/workflow.cpp` 与跨语言一致性测试 `../tests/conformance.cpp`。

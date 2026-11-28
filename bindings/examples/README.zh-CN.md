# 语言绑定示例

[English](README.md) | [简体中文](README.zh-CN.md)

## 宿主回调

每种语言各有一个对应程序，演示向高阶算子（`num_brent`、`num_bfgs`）传入宿主语言函数并读取结果：

| 文件 | 运行命令 |
|---|---|
| [`python/callbacks.py`](python/callbacks.py) | `RSSN_LIB=… python bindings/examples/python/callbacks.py` |
| [`julia/callbacks.jl`](julia/callbacks.jl) | `RSSN_LIB=… julia --project=bindings/julia bindings/examples/julia/callbacks.jl` |
| [`cpp/callbacks.cpp`](cpp/callbacks.cpp) | `cmake --build build/cpp --target rssn_callbacks && ./build/cpp/rssn_callbacks` |
| [`r/callbacks.R`](r/callbacks.R) | `RSSN_LIB=… Rscript bindings/examples/r/callbacks.R` |

源码检出下先构建共享库并设置 `RSSN_LIB`（Linux `librssn.so`、macOS `librssn.dylib`、Windows `rssn.dll`）：

```sh
bazel build //:ffi_dev
export RSSN_LIB="$PWD/bazel-bin/ffi_dev/ffi/librssn.so"
```

示例依次演示：`num_brent` 匿名/闭包回调求根、生成门面调用、`num_bfgs` 目标函数与梯度双回调、显式 `Callback` 复用，以及回调异常向 `RssnStatus` 错误文本的传播。

## 精确有理张量

以下程序构建包含 65 位整数与精确分数 `2/3` 的 2×2 矩阵，基于参考值 `1` 归一化并计算精确行列式 `73786976294838206466`：

| 文件 | 运行命令 |
|---|---|
| [Python](python/exact_tensor.py) | `PYTHONPATH=bindings/python python bindings/examples/python/exact_tensor.py` |
| [Julia](julia/exact_tensor.jl) | `julia --project=bindings/julia bindings/examples/julia/exact_tensor.jl` |
| [R](r/exact_tensor.R) | `Rscript bindings/examples/r/exact_tensor.R` |
| [C++](cpp/exact_tensor.cpp) | `cmake --build build/cpp --target rssn_exact_tensor` 后运行可执行文件 |

运行前请设置 `RSSN_LIB` 指向已构建的共享库（R 示例还需已安装 `rssn` 包）。

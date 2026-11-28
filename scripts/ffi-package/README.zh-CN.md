# RSSN 原生库与 SDK

[English](README.md)

此归档含 RSSN 静态库、共享库、[C 头文件](rssn.h)、[C ABI 的 C++ 头文件](rssn.hpp)，以及 [C++、Python、Julia 和 R SDK 源码](bindings/README.md)。请选择与归档的操作系统、处理器架构和 ABI 相匹配的语言运行时和编译器。

预编译库包含 Cranelift JIT、NumPy 数组读写和 wgpu。新会话会自动选择 CPU 后端，运行环境无法生成机器码时仍可使用解释器。通过 `rssn_manifest_json()` 查看编译时启用的能力，并通过 `num_backend_current_name()` 查询当前后端。GPU 执行还需要运行环境提供设备适配器。[执行后端](docs/BACKENDS.zh.md)

## 编译示例

在包含本说明和头文件的解压目录运行下列命令。

### Linux

```sh
cc -I. examples/c/smoke.c ./librssn.a -lm -lpthread -ldl -o smoke
./smoke
```

### macOS

```sh
cc -I. examples/c/smoke.c ./librssn.dylib -Wl,-rpath,@loader_path -o smoke
./smoke
```

### Windows（MSVC）

打开适用于目标架构的 Visual Studio 开发人员命令提示符。

```bat
cl /nologo /utf-8 /I. examples\c\smoke.c /link rssn.dll.lib /OUT:smoke.exe
smoke.exe
```

`rssn.dll.lib` 是导入库。在其他目录部署程序时，请将 `rssn.dll` 与可执行文件放在一起。静态归档为 `rssn.lib`。静态链接还需加入 RSSN Rust 构建通过 `--print native-static-libs` 列出的目标系统库。

## 使用语言 SDK

在解压目录内使用 CMake 构建和测试 C++ SDK。

```sh
cmake -S bindings/cpp -B build-cpp
cmake --build build-cpp --config Release
ctest --test-dir build-cpp --build-config Release --output-on-failure
```

直接使用解压目录中的 Python 或 Julia SDK 时，它们会自动查找解压目录根部的共享库。在 Unix shell 中，可运行以下命令查看加载的库。

```sh
PYTHONPATH="$PWD/bindings/python" python3 -c 'import rssn; print(rssn.library_path())'
julia --project=bindings/julia -e 'using RSSN; println(libpath())'
```

通过 `python -m pip install ./bindings/python` 或 `R CMD INSTALL bindings/r` 安装后，请将 `RSSN_LIB` 设置为解压目录的绝对路径。在 Unix shell 中运行 `export RSSN_LIB="$PWD"`，在 PowerShell 中运行 `$env:RSSN_LIB = (Get-Location).Path`。在 PowerShell 中使用 `R.exe CMD INSTALL bindings/r` 安装 R 包。各语言的后端选择和 NumPy 示例见 [SDK 指南](bindings/README.md#execution-backends)。

## 相关资料

- [精确重载示例](examples/c/operator_contract.c) 从已加载的库取得签名键，
  再计算两个宽整数的最大公约数。
- [嵌入指南](docs/KERNEL.zh.md#c-abi)介绍 C ABI、句柄所有权和错误报告。
- [文档索引](docs/README.zh.md)列出数值 API、执行后端和量纲资料。
- 从[源码检出](@RSSN_SOURCE@)运行 Rust 示例，或按[构建指南](docs/BUILDING.zh.md)重新构建库。
- [变更日志](CHANGELOG.md)和 [Apache 2.0 许可证](LICENSE)。

校验下载时，将归档的 SHA-256 与配套 `.sha256` 文件的第一项比较。Linux 使用 `sha256sum -c ARCHIVE.tar.gz.sha256`，macOS 使用 `shasum -a 256 -c ARCHIVE.tar.gz.sha256`，Windows 使用 `Get-FileHash ARCHIVE.tar.gz -Algorithm SHA256`。将 `ARCHIVE` 替换为实际名称，并在保存下载文件的目录执行命令。

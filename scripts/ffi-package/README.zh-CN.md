# RSSN C 和 C++ 库

[English](README.md)

此归档含 RSSN 静态库、共享库、[C 头文件](rssn.h)、[C++ 头文件](rssn.hpp)和 [C 示例](examples/c/smoke.c)。请选择与归档的操作系统、处理器架构和 ABI 相匹配的编译器。

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

`rssn.dll.lib` 是导入库。在其他目录部署程序时，请将 `rssn.dll` 与可执行文件放在一起。需要静态链接时，使用 `rssn.lib`。

## 相关资料

- [嵌入指南](docs/KERNEL.zh.md#c-abi)介绍 C ABI、句柄所有权和错误报告。
- [文档索引](docs/README.zh.md)列出数值 API、执行后端和量纲资料。
- 从[源码检出](@RSSN_SOURCE@)运行 Rust 示例，或按[构建指南](docs/BUILDING.zh.md)重新构建库。
- [变更日志](CHANGELOG.md)和 [Apache 2.0 许可证](LICENSE)。

校验下载时，将归档的 SHA-256 与配套 `.sha256` 文件的第一项比较。Linux 使用 `sha256sum -c ARCHIVE.tar.gz.sha256`，macOS 使用 `shasum -a 256 -c ARCHIVE.tar.gz.sha256`，Windows 使用 `Get-FileHash ARCHIVE.tar.gz -Algorithm SHA256`。将 `ARCHIVE` 替换为实际名称，并在保存下载文件的目录执行命令。

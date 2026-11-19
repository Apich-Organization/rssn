# RSSN C and C++ libraries

[简体中文](README.zh-CN.md)

This archive contains the RSSN static and shared libraries, the [C header](rssn.h) and [C++ header](rssn.hpp), and a [C example](examples/c/smoke.c). Use a compiler that matches the archive's operating system, architecture, and ABI.

## Compile the example

Run the commands from the directory containing this README and the headers.

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

### Windows with MSVC

Open the Visual Studio Developer Command Prompt for the target architecture.

```bat
cl /nologo /utf-8 /I. examples\c\smoke.c /link rssn.dll.lib /OUT:smoke.exe
smoke.exe
```

`rssn.dll.lib` is the import library. Keep `rssn.dll` beside your executable when you use it outside this directory. Link `rssn.lib` if your program needs the static library.

## Further reading

- [Embedding guide](docs/KERNEL.md#the-c-abi) for the C ABI, handle ownership, and error reporting.
- [Documentation index](docs/README.md) for numerical APIs, backends, and dimensions.
- [Source checkout](@RSSN_SOURCE@) and the [build guide](docs/BUILDING.md) for rebuilding the libraries and running the Rust examples.
- [Changelog](CHANGELOG.md) and [Apache 2.0 license](LICENSE).

To verify the download, compare the archive's SHA-256 digest with the first field of its companion `.sha256` file. On Linux, run `sha256sum -c ARCHIVE.tar.gz.sha256`. On macOS, run `shasum -a 256 -c ARCHIVE.tar.gz.sha256`. On Windows, use `Get-FileHash ARCHIVE.tar.gz -Algorithm SHA256`. Substitute the archive's filename and run the check beside the two downloaded files.

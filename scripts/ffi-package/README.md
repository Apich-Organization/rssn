# RSSN native libraries and SDKs

[简体中文](README.zh-CN.md)

This archive contains the RSSN static and shared libraries, the [C header](rssn.h), the [C ABI's C++ header](rssn.hpp), and the [C++, Python, Julia and R SDK sources](bindings/README.md). Use language runtimes and compilers that match the archive's operating system, architecture and ABI.

The prebuilt libraries include Cranelift JIT, NumPy array I/O, and wgpu. Fresh sessions select a CPU backend automatically and can keep using the interpreter when native code generation is unavailable. Find compiled capabilities in `rssn_manifest_json()` and query the selected backend with `num_backend_current_name()`. GPU execution also needs a runtime adapter. [Execution backends](docs/BACKENDS.md)

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

`rssn.dll.lib` is the import library. Keep `rssn.dll` beside your executable when you use it outside this directory. The static archive is `rssn.lib`. Static linking also needs the native system libraries from the RSSN Rust build's `--print native-static-libs` output.

## Use the language SDKs

From the extracted directory, build and test the C++ SDK with CMake.

```sh
cmake -S bindings/cpp -B build-cpp
cmake --build build-cpp --config Release
ctest --test-dir build-cpp --build-config Release --output-on-failure
```

The in-place Python and Julia SDKs discover the shared library at the extracted
root. In a Unix shell, try them from that directory.

```sh
PYTHONPATH="$PWD/bindings/python" python3 -c 'import rssn; print(rssn.library_path())'
julia --project=bindings/julia -e 'using RSSN; println(libpath())'
```

For Python installed with `python -m pip install ./bindings/python` or R
installed with `R CMD INSTALL bindings/r`, set `RSSN_LIB` to the absolute path
of the extracted directory. Use `export RSSN_LIB="$PWD"` in a Unix shell or
`$env:RSSN_LIB = (Get-Location).Path` in PowerShell. Run `R.exe CMD INSTALL bindings/r`
from PowerShell. The [SDK guide](bindings/README.md#execution-backends) shows
backend selection and NumPy examples for each language.

## Further reading

- [Exact overload example](examples/c/operator_contract.c) reads the loaded
  library's signature key and uses it to compute a greatest common divisor of
  two wide integers.
- [Embedding guide](docs/KERNEL.md#the-c-abi) for the C ABI, handle ownership, and error reporting.
- [Documentation index](docs/README.md) for numerical APIs, backends, and dimensions.
- [Source checkout](@RSSN_SOURCE@) and the [build guide](docs/BUILDING.md) for rebuilding the libraries and running the Rust examples.
- [Changelog](CHANGELOG.md) and [Apache 2.0 license](LICENSE).

To verify the download, compare the archive's SHA-256 digest with the first field of its companion `.sha256` file. On Linux, run `sha256sum -c ARCHIVE.tar.gz.sha256`. On macOS, run `shasum -a 256 -c ARCHIVE.tar.gz.sha256`. On Windows, use `Get-FileHash ARCHIVE.tar.gz -Algorithm SHA256`. Substitute the archive's filename and run the check beside the two downloaded files.

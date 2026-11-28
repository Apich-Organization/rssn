# RSSN Bazel builds

[简体中文](BAZEL.zh.md)

If Bazelisk is installed as `bazelisk` without a `bazel` alias, substitute `bazelisk` in the examples. On Windows, put Bazel's output root on an NTFS volume, for example `bazelisk --output_user_root=D:/bz build //:headers`.

Use Bazel 7.7.1 through Bazelisk, Python 3.11 or later, and the Rust toolchain specified by `Cargo.toml`. Both Bazel and Cargo use `Cargo.lock`. The RSSN module also exports its source, workflow rules and LLVM pass for Aevia.

| Command | Output or check |
| --- | --- |
| `bazel build //:rustlib` | Development Rust library |
| `bazel run //:check` | Workspace type checking without linking test binaries |
| `bazel test //:fmt_test` / `bazel run //:fmt` | Check or format Rust sources using nightly rustfmt |
| `bazel test //:docs_test` | Check documentation, then build and run the functional examples |
| `bazel run //:example -- quickstart` | Run a Rust example with an optional `--profile` |
| `bazel run //:benchmark -- engine -- polynomials` | Run the named benchmark with its arguments |
| `bazel test //:test` | Default Rust suite and generated contracts |
| `bazel test //:rust_minimal_test` | Rust tests with default features disabled |
| `bazel test //:integration_test` | Rust features, contracts, headers, witness replay and all four bindings |
| `bazel test //:clippy_test` | Strict Clippy across the full-feature workspace and all Rust targets |
| `bazel test //:ffi_native_clients_test` | Run native C and C++ client tests with static linkage on Linux and shared linkage on macOS and Windows |
| `bazel test //:rust_all_features_test` | Full-feature Rust tests, including the GPU tests the host can run |
| `bazel test --define=rssn_bounded_tests=1 //:rust_all_features_test` | Run each RSSN integration binary separately and reclaim its generated executables before the next group |
| `bazel test //:rust_all_features_test --test_arg=--batches=4` | Run the full Rust suite in groups of four integration binaries and reclaim their cache space after each group |
| `bazel test //:rust_ignored_test //:witness_test` | Long numerical regressions and replay of the checked-in signatures |
| `bazel test //:bindings_check_test //:bindings_test` | API and Python, Julia, R and C++ tests against the same all-feature cdylib |
| `bazel run //:bindings_from_library -- --library /path/to/librssn.so` | Check generated bindings and exercise an existing cdylib |
| `bazel test //:generated_sources_check_test` | Check generated contracts, structured rules, bindings and codecs |
| `bazel test //:api_generator_test` | Test the source model and generator without linking RSSN |
| `bazel build //:manifest` | Kernel manifest JSON under `bazel-bin/manifest/` |
| `bazel build //:headers` | Generate C and C++ header artifacts under `bazel-bin/headers/` |
| `bazel run //:refresh_headers` | Regenerate `rssn.h` and `rssn.hpp` without compiling libRSSN |
| `bazel test //:headers_check_test` | Check the generated C and C++ headers |
| `bazel run //:refresh_contracts` | Refresh native contracts, bindings and codecs |
| `bazel run //:refresh_bindings` | Refresh the language binding facades with a fresh manifest |
| `bazel run //:refresh_bindings_from_manifest -- --manifest manifest.json` | Refresh the facades from an already linked RSSN build's manifest |
| `bazel run //:refresh_witness` | Collect a fresh test trace and mine `signatures.witnessed.json` |

The binding integration test needs Python with NumPy, Julia, R with development headers, and CMake with a C++ compiler. It installs the R package into the build state and passes the same explicit library file to every suite. Add `--config=release` to exercise the optimized cdylib. For a machine without the full language toolchain, `bazel build //:ffi_dev` supplies the all-feature development cdylib.

When R is installed in an environment with a separate C++ toolchain, select a
compiler compatible with the supplied library using `--repo_env=CXX=/path/to/c++`.

## FFI distribution

```sh
bazel build //:ffi_distribution
bazel run //:package_ffi -- --dist-dir dist
```

`//:ffi_distribution` builds both libraries in one rustc invocation and places
them under `bazel-bin/ffi_distribution/ffi/`. It enables every Cargo feature and
loads the LLVM optimizer. The archive contains the libraries, generated
headers, examples and documentation, with a companion SHA-256 file.

For development, `bazel build //:ffi_distribution_dev` builds both libraries
under `bazel-bin/ffi_distribution_dev/ffi/`. Run
`bazel run //:package_ffi_dev -- --dist-dir dist` to package these libraries
with generated headers. The default archive name ends in `-dev.tar.gz`.

The distribution profile defaults to O3, generic CPU flags, unwind semantics,
fat LTO and one codegen unit. The i686 targets retain Rust's SSE2-capable target
baseline. The Bazel dependency `//:release_pass_gate` runs
the release LLVM probe before building the two libraries. `//:release_ffi`
packages the libraries and runs the extracted C, C++ and Python clients before
publishing the archive to the distribution directory. To check an existing archive, run
`bazel run //:ffi_package_smoke -- --archive /path/to/rssn.tar.gz`.

`//:package_ffi_from_binaries` accepts `--native-dir` and `--headers-dir` for a matching set of existing artifacts and runs the same archive validation.

### LLVM optimizer

The LLVM pass requires CMake, a C++ compiler and the host LLVM SDK
matching the full version in `rustc -vV`. Bazel looks for the corresponding
`llvm-config-<major>` or accepts an explicit path. On a Windows build host set
`--repo_env=RSSN_LLVM_IMPORT_LIBRARY=/path/to/host-llvm.lib` to the import library
exporting the LLVM C++ API of the process that loads the plugin. CMake selects
the host module format and link options. Where the stock rustc lacks the LLVM
exports used by the plugin, set `--repo_env=RUSTC=/absolute/path/to/rustc` to a
compatible compiler and optionally set `CARGO` to its Cargo executable. CI accepts
SDK and compiler paths through the JSON maps `RSSN_RELEASE_LLVM_CONFIGS`,
`RSSN_RELEASE_LLVM_IMPORT_LIBRARIES`, `RSSN_RELEASE_RUSTCS` and `RSSN_RELEASE_CARGOS`,
keyed by the Rust target triple. The SDK and compiler must be installed on the
runner. CI also accepts the existing shared `RSSN_RUSTC` and `RSSN_CARGO` paths.

```sh
bazel test --repo_env=RSSN_LLVM_CONFIG=/path/to/llvm-config //:llvm_optimizer_test //:release_pass_test
bazel build --config=release --repo_env=RSSN_LLVM_CONFIG=/path/to/llvm-config //:rustlib
bazel build --repo_env=RSSN_LLVM_CONFIG=/path/to/llvm-config //:ffi //:ffi_static
bazel run --repo_env=RSSN_LLVM_CONFIG=/path/to/llvm-config //:release_ffi -- --dist-dir dist
```

`//:ffi` builds the dynamic library and its import library where applicable.
`//:ffi_static` builds the static archive. Use `//:package_ffi_portable`
for an explicit plugin-free build.
The host plugin is a Rust action input passed to target rustc
with `-Zllvm-plugins`. The standalone Rust library uses its O3 pipeline with
LTO off so the pass runs during the library build. The optional optimizer
enables portable power, endpoint-guard and unit-log transforms. The Linux
x86-64 build also uses the standard spill-pressure transform. The Windows x64 MSVC target also
enables the Win64-compatible libm passes. Individual passes decide if a
target and function are eligible.

For Windows MSVC, Linux GNU or musl, and FreeBSD targets, set
`--define=rssn_llvm_standalone=1` to run the LLVM plugin on post-LTO bitcode in
the SDK's `opt`. The plugin uses the matching SDK's LLVM library, including its
import library on a Windows build host. The distribution uses all features,
fat LTO and one codegen unit. Bazel replaces the final codegen object with the
output of `llc` and passes the original Cargo linker arguments back to rustc.
The ELF route emits a position-independent object for the shared library.
Pass the same define to `//:release_pass_test` to check the platform's shared
and static libraries, any import library, and Cargo's dependency tracking.

```sh
bazel test --config=release-max --define=rust_target=x86_64-pc-windows-msvc --define=rssn_llvm_standalone=1 //:release_pass_test
bazel run --config=release-max --define=rust_target=x86_64-pc-windows-msvc --define=rssn_llvm_standalone=1 //:release_ffi -- --dist-dir dist
```

### Cross-compilation

Set `--define=rust_target=<triple>` and install that Rust target plus its native linker, C toolchain and SDK. The optional plugin runs on the build host. A JSON file selected by `--repo_env=RSSN_TARGET_ENV=/absolute/path/to/toolchains.json` maps Rust triples to toolchain environment variables. Bazel reads the file as a configuration input. Here is a Linux-to-AArch64 example:

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

Choose the corresponding Cargo `_LINKER` and `_RUNNER` keys for Windows MSVC, macOS, FreeBSD, or other Rust triples. The target toolchain can supply its own C/C++ compiler, archiver, sysroot, `PKG_CONFIG_*`, Windows `INCLUDE` and `LIB`, or macOS `SDKROOT`. A linker is sufficient for the cross-build probe. Rust tests can execute through a configured runner. The Python, Julia, R and C++ binding integration suite loads the library on the target platform. Native development targets do not need the LLVM SDK.

## Cache and timing

Bazel tracks the declared first-party sources, Cargo manifests and locks, and the LLVM plugin independently. The Rust actions use Cargo to schedule compiler invocations in a stable local source mirror. Unchanged files keep their original modification times across builds, so an action invalidated by a source edit can reuse Cargo fingerprints, third-party compilation and rustc's development CGU cache. Deleted inputs are removed from the mirror. Toolchain configuration is part of the local state key. The release plugin path in rustflags includes a digest of the library and its optimization settings, so changing either invalidates the Cargo artifacts that used them. Builds for distribution use generic CPU flags and deterministic RSSN metadata.

The default local build state is `.bazel-cache/`. On Windows use a short absolute cache path, for example `--repo_env=BAZEL_SCI_CACHE=C:/b/r`. A separate `--disk_cache` can retain Bazel action outputs. For an existing compiler cache, set `--repo_env=RSSN_SCCACHE=/path/to/sccache`. Cargo continues to use incremental development builds where applicable. After replacing toolchain binaries or an SDK in place, run `bazel sync --configure` before measuring. The local actions depend on host toolchains and are not scheduled on remote executors.

Compare release settings on the same machine, toolchain, input and target. Keep distinct cache directories for cold measurements. Reuse the cache and exact invocation for warm measurements. Use `--define=rssn_timings=1` to save the Cargo timing report in the printed target directory. A Bazel `--profile=/path/to/profile.gz` measures the actions and Bazel overhead. Use `--define=rssn_jobs=N` to set the nested Cargo and CMake job limit independently of Bazel’s scheduling limit. To compare the RSSN crate’s parallel frontend, add `--define=rssn_frontend_threads=4`. The wrapper applies `-Zthreads` only to RSSN and keeps that experiment in a separate Cargo cache.

- `--config=release` uses ThinLTO and 16 codegen units.
- `--config=release-local` disables LTO and retains O3. The optimizer remains enabled on LLVM targets.
- `--config=release-max` uses fat LTO and one codegen unit. Override either setting for controlled experiments with `--define=rssn_lto` and `--define=rssn_codegen_units`.

These comparisons apply to linked artifacts. Bazel release configurations are separate from the legacy Cargo release profile. Evaluate numerical parity and runtime of the workload before choosing a different distribution profile.

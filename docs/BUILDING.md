# Building Rust and C artifacts

[中文](BUILDING.zh.md)

Bazel builds release archives and runs Rust, binding and generator workflows.
The [RSSN Bazel guide](BAZEL.md) lists the targets and toolchain settings.

Ordinary `cargo build` and downstream Rust builds use the Rust library (`rlib`).
`cargo check` emits Rust metadata for type checking. Build the C static and
shared libraries explicitly when preparing C artifacts.

On Linux, run these commands to refresh the checked-in cbindgen headers and
build the libraries using the Rust toolchain.

```sh
bazel run //:refresh_headers
bazel build //:ffi_distribution
cc -I. examples/c/smoke.c bazel-bin/ffi_distribution/ffi/librssn.a -lm -lpthread -ldl -o smoke
./smoke
```

The default FFI distribution includes the Cranelift JIT, NumPy `.npy` I/O and
wgpu. It uses a generic target CPU and unwind semantics so guarded C entries
can return `RssnStatus::Panic`. FFT planning can still select CPU instructions
at runtime. C linker flags vary by platform.

Run `bazel run //:package_ffi` or `bash scripts/build_ffi.sh --package`.
On Windows, use `scripts/build_ffi.ps1 -Package`. The archive contains the
libraries, generated headers, examples and bilingual guides, with a companion
SHA-256 file. The release uses all features, fat LTO and the LLVM optimizer.
Configure the host's matching LLVM SDK before building. Set `RSSN_FFI_PROFILE=dev`
to stage local development libraries under `target/debug` through Bazel.
Use `bazel run //:release_ffi` to run the packaged client tests before publishing.
The `--llvm-overlay` and `-LlvmOverlay` wrapper flags retain the `-llvm` archive suffix.

The checkout configures `target-cpu=native` for direct Cargo development builds.
The Bazel `//:refresh_headers` target invokes cbindgen independently of the
Rust library.

The normal build script uses only the Rust standard library. Header generation
requires the optional `dev-headers` feature as well as `DEV=1`; this keeps
cbindgen and its dependencies out of ordinary Rust builds. Setting `DEV=1`
without that feature reports an error so stale headers cannot be mistaken for a
successful regeneration. `signatures.overlay.json` is required, version-controlled
source and is copied into the build output unchanged. Git checkouts, source
archives, and Cargo packages must include it. Missing source fails the build;
invalid signatures fail manifest construction instead of silently changing the
typed API. Local signature edits are source changes and should be committed
alongside the code that uses them.

`signatures.witnessed.json` (the signatures mined from the test corpus, see
[KERNEL.md](KERNEL.md)) is required in exactly the same way: the build fails with a
message naming the file when it is missing, so two builds of the same sources
cannot produce different manifests. Collect a fresh test trace with
`RSSN_SORT_TRACE=sorts.tsv cargo test --lib`, then regenerate the table with
`cargo run --profile tools --bin signature-mine -- sorts.tsv`. Operators without
a witness keep their declared, derived, or overlay signatures.

Build metadata honors `SOURCE_DATE_EPOCH`, `VERGEN_IDEMPOTENT`, and explicit
`VERGEN_*` overrides. On Linux it obtains host details from `/etc/os-release`,
`/proc`, and `uname`; platforms or source archives without those facilities or
Git use `unknown` for unavailable details. The reported CPU count uses Rust's
available parallelism and may reflect a container quota. For reproducible
metadata across hosts, set `VERGEN_IDEMPOTENT=1` with `SOURCE_DATE_EPOCH`, and
supply explicit `VERGEN_*` values for any desired archive provenance.

## Numerical benchmarks

The `polynomials` benchmark group measures both Lagrange input layouts and
compares scalar Horner evaluation with the batch kernel on the same points.
Use the same benchmark source, Cargo target directory, feature flags, Rust
toolchain and worker count for both revisions. Save a baseline with the first
command on the earlier revision. Run the second on the revised kernels.

```sh
RAYON_NUM_THREADS=4 cargo bench --locked --bench engine -- polynomials --save-baseline before
RAYON_NUM_THREADS=4 cargo bench --locked --bench engine -- polynomials --baseline before
```

Use the `fixed_prefix_scan` group to compare one-block scans with fixed
256-element blocks. Use the same worker count and baseline flags.

The `numerical_codec_dispatch` group measures repeated value calls on an
already constructed session. The `mesh_triangle_areas` group uses the same
cells with execution capped at one or four workers. Mesh timings exclude pool
setup. Use the same benchmark source and execution settings on both revisions.

The `special_combinatorics` group covers balanced and unbalanced log
coefficients, binomial overflow and multinomial partitions. The
`beta_quantiles` group covers an interior probability, an upper-tail
probability and an analytic tiny-tail case. Use the same baseline protocol.

The `finite_field_ntt` group measures the same prime's cold generator
setup and cached lookup, and polynomial convolutions at two input sizes
with warm plans. The cold case clears the process-local reuse pool before
each timed iteration. Keep the worker count fixed for comparisons.

The `kalman_history` group measures forward filtering and smoothed means
for one fixed model at three observation counts. Compare revisions with
the same benchmark source and fixed worker count.

## Distribution channels

The reusable CI workflow verifies the generated language bindings against an all-features Linux build and runs Python, Julia, R and C++ against that same shared library. The release workflow requires that CI result, including for nightly builds.

Push a version tag that matches `Cargo.toml` to start the versioned release workflow. The workflow tests the tagged commit and uploads the checksummed C archives to a draft. Publishing the draft locks its tag and assets and creates a GitHub release attestation. Check the package version and SHA-256 file against the archive before linking. With authenticated GitHub CLI, verify the downloaded archive with `gh release verify-asset TAG ARCHIVE --repo Apich-Organization/rssn`.

Daily builds from the repository's default branch use UTC tags `nightly-YYYY-MM-DD`. Each is an immutable GitHub prerelease. The crates.io job applies only to versioned publications and requires environment approval. The **Release** workflow's `dry_run` option runs validation and builds without publishing.

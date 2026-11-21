# Building Rust and C artifacts

[中文](BUILDING.zh.md)

Ordinary `cargo build` and downstream Rust builds use the Rust library (`rlib`).
`cargo check` emits Rust metadata for type checking. Build the C static and
shared libraries explicitly when preparing C artifacts.

On Linux, build the C libraries and refresh the checked-in cbindgen headers with the following commands.

```sh
DEV=1 bash scripts/build_ffi.sh --features dev-headers
cc -I. examples/c/smoke.c target/ffi/librssn.a -lm -lpthread -ldl -o smoke
./smoke
```

The script invokes `cargo rustc --lib --crate-type staticlib,cdylib`. Static and shared libraries are placed under the selected target and profile
directory. The default `ffi` profile inherits the release optimizations and
enables unwinding, so guarded C entry points can report unwinding Rust panics
as `RssnStatus::Panic`. Set `RSSN_FFI_PROFILE=dev` for a development C build in
`debug`. The scripts accept `CARGO_TARGET_DIR` and Cargo options such as
`--features jit`, `--target`, `--locked`, and `--offline`. C linker flags
vary by platform.

Use `bash scripts/build_ffi.sh --package` or
`scripts/build_ffi.ps1 -Package` on Windows to assemble a native archive with
its headers, examples, bilingual guides, and checksum.

The checkout configures `target-cpu=native` as the general Cargo build default.
For `x86_64-pc-windows-msvc`, Cargo selects the target-specific stack linker flags.
The release workflow sets `RUSTFLAGS="-C target-cpu=generic"` for the
distributed archives. Use that setting when building a redistributable archive
or targeting another machine. FFT planning can still select CPU instructions
at runtime.

Migration: C consumers previously using `cargo build --release` to obtain
`librssn.a` or a shared library should use this explicit distribution command.
Rust callers continue using ordinary Cargo commands. Merely enabling a Cargo
feature does not change the crate type. Header generation remains opt-in via
`DEV=1`; use the existing headers when regeneration is unnecessary.

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

## Distribution channels

Push a version tag that matches `Cargo.toml` to start the versioned release workflow. The workflow tests the tagged commit and uploads the checksummed C archives to a draft. Publishing the draft locks its tag and assets and creates a GitHub release attestation. Check the package version and SHA-256 file against the archive before linking. With authenticated GitHub CLI, verify the downloaded archive with `gh release verify-asset TAG ARCHIVE --repo Apich-Organization/rssn`.

Daily builds from the repository's default branch use UTC tags `nightly-YYYY-MM-DD`. Each is an immutable GitHub prerelease. The crates.io job applies only to versioned publications and requires environment approval. The **Release** workflow's `dry_run` option runs validation and builds without publishing.

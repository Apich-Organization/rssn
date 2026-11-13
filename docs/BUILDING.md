# Building Rust and C artifacts

Ordinary `cargo build`, `cargo check`, and downstream Rust dependency builds
produce the Rust library (`rlib`). They do not also link a shared C library or
copy all dependencies into a static C archive. This keeps each Rust build from
paying the distribution-artifact cost.

To produce the C libraries and refresh the checked-in cbindgen headers:

```sh
DEV=1 bash scripts/build_ffi.sh --features dev-headers
cc -I. examples/c/smoke.c target/release/librssn.a -lm -lpthread -ldl -o smoke
./smoke
```

The script invokes `cargo rustc --lib --crate-type staticlib,cdylib`. The C ABI
and filenames are unchanged: static and shared libraries are placed under the
selected target/profile directory. `CARGO_TARGET_DIR` and Cargo arguments such
as `--features jit`, `--target`, `--locked`, or `--offline` are respected. The
default profile is `release`; set `RSSN_FFI_PROFILE=dev` for a development C
build (its output directory is `debug`). Platform-specific C linker flags may
be needed outside Linux.

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
cannot produce different manifests. Regenerate it with
`cargo run --profile tools --bin signature-mine` after adding operators; until
then new operators carry their derived (and overlay) signatures only.

Build metadata honors `SOURCE_DATE_EPOCH`, `VERGEN_IDEMPOTENT`, and explicit
`VERGEN_*` overrides. On Linux it obtains host details from `/etc/os-release`,
`/proc`, and `uname`; platforms or source archives without those facilities or
Git use `unknown` for unavailable details. The reported CPU count uses Rust's
available parallelism and may reflect a container quota. For reproducible
metadata across hosts, set `VERGEN_IDEMPOTENT=1` with `SOURCE_DATE_EPOCH`, and
supply explicit `VERGEN_*` values for any desired archive provenance.

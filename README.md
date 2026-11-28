# RSSN: Scientific Computing for Rust

[English](README.md) | [简体中文](README.zh-CN.md)

[![Crates.io](https://img.shields.io/crates/v/rssn.svg)](https://crates.io/crates/rssn)
[![Docs.rs](https://docs.rs/rssn/badge.svg)](https://docs.rs/rssn)
[![License](https://img.shields.io/crates/l/rssn)](LICENSE)
[![Discord Server](https://img.shields.io/discord/1459399539403522074.svg?label=Discord&logo=discord&color=blue)](https://discord.gg/D5e2czMTT9)
[![DOI](https://zenodo.org/badge/DOI/10.6084/m9.figshare.31044715.svg)](https://doi.org/10.6084/m9.figshare.31044715)

**rssn** is an open-source scientific computing library for Rust: a computer algebra system, numerical kernels, and physics simulations behind one unified API.

Every symbolic operation is an *identity transformation*. Requests such as `diff(f, x)` are rewritten by rule sets on a hash-consed **e-graph DAG** until an equal closed form or numeric witness is extracted. See [ARCHITECTURE.md](ARCHITECTURE.md).

---

## Key Features

- **Single entry point**: `Session::compute(term, &config)` returns a symbolic closed form or numeric value depending on the target phase.
- **Dimensions & units**: Scalar terms infer physical dimensions with exact rational/symbolic exponents, SI/affine units, and advisory or strict checking ([docs/DIMENSIONS.md](docs/DIMENSIONS.md)).
- **Open e-graph engine**: Hash-consed DAG, extensible operator registry, tiered rule scheduling, cost-based extraction, assumptions (`positive`, `integer`, `real`), and `cse` common-subexpression extraction.
- **Calculus & differential equations**: Derivatives, indefinite/definite/improper integrals (including Hermite reduction), limits, series, Meijer G, Gosper summation, integral transforms, ODEs (exact, Lie symmetries, systems), 1D/2D/3D PDEs, and calculus of variations.
- **Algebra & analysis**: Multivariate factorization, Gröbner bases (Buchberger/F4/F5), CAD quantifier elimination, symbolic/numeric linear algebra, number theory, combinatorics, special functions, complex/functional analysis, and statistics.
- **Numeric kernels & signal processing**: Adaptive/DE quadrature, stiff & non-stiff ODEs, root finding, dense/sparse/Krylov linear algebra, conic/MILP optimization, FFT/wavelets, graphs, and arbitrary-precision (`kernels::mp`) ([docs/MATH_COVERAGE.md](docs/MATH_COVERAGE.md)).
- **Physics & simulations**: Rule sets for classical, EM, relativity, quantum, QFT, solid state, and thermodynamics, plus FDM/FEM/FVM/BEM/CFD, MD, SPH, spectral, and multigrid solvers.
- **Compilation & devices**: Lower closed forms via `Backend` to an SSA tape (with reverse autodiff and intervals), Cranelift JIT (`jit`), or CPU/wgpu devices (`gpu-wgpu`, with `f64` double-single emulation) ([docs/BACKENDS.md](docs/BACKENDS.md)).
- **I/O & visualization**: Round-trippable text, LaTeX, Typst, Unicode pretty-printing, plots, `.npy` (`npy`), DataFrames, Parquet, HDF5, and NetCDF.
- **Embedding & C ABI**: Deterministic *Kernel Manifest*, dynamic ABI (`Session::call_value`), host callbacks ([docs/CALLBACKS.md](docs/CALLBACKS.md)), shared-library plugins (`plugin-cdylib`), and C/C++ headers (`rssn.h`/`rssn.hpp`) ([docs/KERNEL.md](docs/KERNEL.md)).

---

## Quick Start (3.0 Beta API)

Add RSSN from a local checkout of the `v0.3` branch:

```bash
git clone --branch v0.3 https://github.com/Apich-Organization/rssn.git ../rssn
cargo add --path ../rssn
```

Save as `src/main.rs` and run `cargo run`:

```rust
use rssn::api::*;

fn main() -> Result<(), ComputeError> {
    let s = Session::new();
    let x = s.sym("x");

    // Symbolic phase: a closed form.
    let d = s.compute(diff(sin(x) * cos(x), x), &Config::new())?;
    assert_eq!(d.term.to_string(), "cos(2*x)");

    // Numeric phase: evaluated at x = 0.0.
    let n = s.compute(
        diff(sin(x) * cos(x), x),
        &Config::new().numeric(1e-12).bind("x", 0.0),
    )?;
    assert_eq!(n.value, Some(1.0));

    // Requests can also be parsed from text.
    let sum = s.compute(s.parse("sum(1/(k*(k+1)), k, 1, n)")?, &Config::new())?;
    assert_eq!(sum.term.to_string(), "1 - 1/(n + 1)");

    Ok(())
}
```

Explore more in the [examples index](examples/README.md):
- `bazel run //:example -- --profile release quickstart`
- `bazel run //:example -- --profile release simulation`
- `examples/c/smoke.c` (C client)
- `bazel run //:benchmark -- engine` (benchmarks)

### Cargo Features

Default features are `jit`, `npy`, and `gpu-wgpu`. Use `--no-default-features` for a minimal Rust build, or `bazel build //:ffi_distribution //:headers` for C/C++ libraries ([docs/BUILDING.md](docs/BUILDING.md)).

| Feature | Enables |
|---|---|
| `jit` | Cranelift JIT backend (`rssn::backend::jit`) |
| `npy` | NumPy `.npy` import/export in `rssn::io` |
| `plugin-cdylib` | Dynamic shared-library plugins (runs native code in host process) |
| `dev-headers` | Regenerating `rssn.h`/`rssn.hpp` with cbindgen (requires `DEV=1`) |
| `gpu-wgpu` | wgpu GPU device (`rssn::backend::device`, `f32` and double-single `f64`; not in `full`) |
| `full` | `npy`, `jit`, and `plugin-cdylib` |

---

## Testing

```bash
bazel test //:rust_tests
bazel test //:rust_tests --test_arg=--test --test_arg=dimensions
bazel test //:rust_tests --test_arg=--test --test_arg=sim --test_arg=acoustics
bazel test //:gpu_test
bazel run //:check -- --lib --tests   # fast compile check (~2 min)
```

Full unit and `sim` suites are compute-intensive (~19–25 min); use `--test_arg` filters for focused runs. CI runs `//:rust_all_features_test`, `//:api_generator_test`, `//:contracts_test`, and `//:headers_check_test`. See [CONTRIBUTING.md](CONTRIBUTING.md).

---

## Documentation

| Document | Subject |
|---|---|
| [ARCHITECTURE.md](ARCHITECTURE.md) | Engine, rule sets, kernels, simulation, C ABI, and JIT |
| [docs/README.md](docs/README.md) | Full documentation index |
| [docs/KERNEL.md](docs/KERNEL.md) | Embedding kernel: manifest, dynamic ABI, plugins, units |
| [docs/DIMENSIONS.md](docs/DIMENSIONS.md) | Dimension and unit layer |
| [docs/CALLBACKS.md](docs/CALLBACKS.md) · [docs/NUMERICAL_CALLBACKS.md](docs/NUMERICAL_CALLBACKS.md) | Host callbacks and numerical callback operators |
| [docs/BUILDING.md](docs/BUILDING.md) | Building Rust and C artifacts |
| [docs/BACKENDS.md](docs/BACKENDS.md) | Execution backends: CPU, wgpu, double-single `f64`, launch fusion |
| [docs/MATH_COVERAGE.md](docs/MATH_COVERAGE.md) | Numerical conventions, coverage, and known gaps |
| [CHANGELOG.md](CHANGELOG.md) | Release notes |
| [API reference](https://docs.rs/rssn) | Generated rustdoc |

---

## Community & Policies

- **Contributing**: [CONTRIBUTING.md](CONTRIBUTING.md)
- **License**: Apache 2.0 ([LICENSE](LICENSE)) · [ATTRIBUTIONS.md](ATTRIBUTIONS.md)
- **Security**: [SECURITY.md](SECURITY.md)
- **Code of Conduct**: [CODE_OF_CONDUCT.md](CODE_OF_CONDUCT.md)

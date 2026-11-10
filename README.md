# rssn: A High-Performance Scientific Computing Library for Rust

[![Crates.io](https://img.shields.io/crates/v/rssn.svg)](https://crates.io/crates/rssn)
[![Docs.rs](https://docs.rs/rssn/badge.svg)](https://docs.rs/rssn)
[![License](https://img.shields.io/crates/l/rssn)](LICENSE)
[![Scc Count Badge Code](https://sloc.xyz/github/Apich-Organization/rssn/?category=code)](https://github.com/Apich-Organization/rssn/)
[![Scc Count Badge Blanks](https://sloc.xyz/github/Apich-Organization/rssn/?category=blanks)](https://github.com/Apich-Organization/rssn/)
[![Scc Count Badge Lines](https://sloc.xyz/github/Apich-Organization/rssn/?category=lines)](https://github.com/Apich-Organization/rssn/)
[![Scc Count Badge Comments](https://sloc.xyz/github/Apich-Organization/rssn/?category=comments)](https://github.com/Apich-Organization/rssn/)
[![Scc Count Badge Cocomo](https://sloc.xyz/github/Apich-Organization/rssn/?category=cocomo)](https://github.com/Apich-Organization/rssn/)
[![Scc Count Badge Effort](https://sloc.xyz/github/Apich-Organization/rssn/?category=effort)](https://github.com/Apich-Organization/rssn/)
[![Discord Server](https://img.shields.io/discord/1459399539403522074.svg?label=Discord&logo=discord&color=blue)](https://discord.gg/D5e2czMTT9)
[![DOI](https://zenodo.org/badge/DOI/10.6084/m9.figshare.31044715.svg)](https://doi.org/10.6084/m9.figshare.31044715)

**rssn** is an open-source scientific computing library for Rust: a computer algebra system, numerical kernels and physics simulations behind one small API.

Every symbolic operation is an *identity transformation*. Differentiating, integrating, solving and simplifying are all requests like `diff(f, x)`, rewritten by rule sets on a hash-consed **e-graph DAG** until an equal closed form (or a number) is found. See [ARCHITECTURE.md](ARCHITECTURE.md).

---

## Key Features

- **One entry point**: `Session::compute(term, &config)`. The answer is symbolic or numeric, depending on the target phase.
- **Dimensions**: every scalar term has an inferred dimension, every operator declares how dimensions combine, and answers report theirs. Exact rational and symbolic exponents (`x^n`, `x^(a+b)`, `x^(1/n)`), SI and affine units, user-defined base dimensions; advisory by default, strict hard errors on request, explicit bypass. See [docs/DIMENSIONS.md](docs/DIMENSIONS.md).
- **Open e-graph engine**: a hash-consed DAG with equivalence classes, an operator registry open to plugins, tiered rule scheduling, cost-based extraction, assumptions such as `positive`, `integer` and `real`, and a numeric soundness check of every rewrite.
- **Calculus**: derivatives (including of undetermined functions); indefinite, definite and improper integrals (including Hermite reduction in `Q(x)(√D)` for algebraic integrands); limits; univariate and multivariate series; Meijer G evaluation, reduction, differentiation and integration; Gosper summation, indefinite sums and products (`indefinite_sum`, `indefinite_product`); integral transforms.
- **Algebra**: rational normal form, multivariate factorisation, Gröbner bases (Buchberger, F4 and F5), polynomial, transcendental and inequality solving, systems of equations (with free parameters), real quantifier elimination by CAD, general trigonometric solution families, linear algebra (symbolic and numeric), number theory (including a quadratic sieve for large composites), combinatorics, logic.
- **Differential equations**: ODE classification and solution (linear, exact, Bernoulli, Riccati, Clairaut/Lagrange, integrating factors, Lie symmetries, systems with first integrals); PDE solution in 1D/2D/3D on boxes, disks, cylinders, balls, annuli, shells, half-spaces and exterior domains with Dirichlet/Neumann/Robin/periodic conditions (eigenfunction expansions, images, Kirchhoff/Poisson/Riemann formulas, linear systems); calculus of variations (Euler–Lagrange).
- **Analysis**: special functions (Γ, ζ, Bessel, erf, orthogonal polynomials, …), complex analysis, functional analysis and integral equations (Fredholm, Volterra), optimisation, statistics.
- **Physics**: formula rule sets for classical mechanics, electromagnetism, relativity, quantum mechanics, QFT, solid state and thermodynamics. Also simulation modules: FDM, FEM, FVM, BEM, CFD, molecular dynamics, SPH, spectral and multigrid methods, plus ready-made scenarios.
- **Output**: round-trippable text, LaTeX, Typst, Unicode pretty-printing, plots, `.npy` export.
- **Compilation**: compile closed forms into fast callable functions through the `Backend` trait: an SSA tape (with reverse gradients, Hessian-vector products and interval enclosures) and, with the `jit` feature, Cranelift machine code with caching and tiered compilation.
- **Expression optimisation**: `cse(list(...))` extracts shared subexpressions of several expressions into temporaries.
- **Numeric kernels**: adaptive and double-exponential quadrature, cubature, stiff and non-stiff ODE solvers (Dormand–Prince, Rosenbrock, Radau IIA, BDF), root finding (Brent, Jenkins–Traub, Aberth, homotopy), dense and sparse linear algebra, Krylov solvers, Chebyshev/AAA approximation, optimisation, robust statistics; connected to the symbolic side (`nminimize`, `chebyshev_approx`, compiled numeric `defint`/`odeint`).
- **Host callbacks**: concrete checked callback signatures and serial numerical operators (quadrature, Brent, nonlinear least squares, four ODE integrators) through the dynamic ABI, including physical dimensions and error propagation. See [callback embedding](docs/CALLBACKS.md) and [the numerical callback operators](docs/NUMERICAL_CALLBACKS.md).
- **Embedding kernel**: a self-describing, deterministic *Kernel Manifest* of every operator, typed signature, unit and physical constant; a dynamic call ABI (`Session::call_value`), plugins (in-process, or shared libraries with the `plugin-cdylib` feature) and a deterministic `bincode-next` wire format. See [the embedding kernel](docs/KERNEL.md).
- **C interface**: a small, panic-safe C ABI (`rssn.h` / `rssn.hpp`), generated by cbindgen. Build C distribution libraries explicitly with `bash scripts/build_ffi.sh`; see [building Rust and C artifacts](docs/BUILDING.md).

---

## Quick Start

```bash
cargo add rssn
```

```rust
use rssn::api::*;

let s = Session::new();
let x = s.sym("x");

// Symbolic phase: a closed form.
let d = s.compute(diff(sin(x) * cos(x), x), &Config::new())?;
assert_eq!(d.term.to_string(), "cos(2*x)");

// Numeric phase: the same request, evaluated.
let n = s.compute(diff(sin(x) * cos(x), x), &Config::new().numeric(1e-12).bind("x", 0.0))?;
assert_eq!(n.value, Some(1.0));

// Any request can be parsed from text.
let sum = s.compute(s.parse("sum(1/(k*(k+1)), k, 1, n)")?, &Config::new())?;
assert_eq!(sum.term.to_string(), "1 - 1/(n + 1)");
# Ok::<(), ComputeError>(())
```

More examples:

- `cargo run --release --example quickstart`: a tour of the API.
- `cargo run --release --example simulation`: simulations by name with JSON parameters.
- `examples/c/smoke.c`: a C client. Build instructions are in the file.

Benchmarks: `cargo bench --bench engine`.

### Cargo features

| Feature | Enables |
|---|---|
| `jit` | the Cranelift JIT backend (`rssn::backend::jit`) |
| `npy` | reading and writing `.npy` files in `rssn::io` |
| `plugin-cdylib` | loading plugins from shared libraries (runs arbitrary code in the host process) |
| `dev-headers` | regenerating `rssn.h`/`rssn.hpp` with cbindgen (together with `DEV=1`) |
| `full` | `npy`, `jit` and `plugin-cdylib` |

---

## Documentation

| Document | Subject |
|---|---|
| [ARCHITECTURE.md](ARCHITECTURE.md) | How the engine, rule sets, kernels, simulation, C ABI and JIT fit together |
| [docs/KERNEL.md](docs/KERNEL.md) | The embedding kernel: manifest, dynamic ABI, plugins, units |
| [docs/DIMENSIONS.md](docs/DIMENSIONS.md) | The dimension and unit layer |
| [docs/CALLBACKS.md](docs/CALLBACKS.md), [docs/NUMERICAL_CALLBACKS.md](docs/NUMERICAL_CALLBACKS.md) | Host callbacks |
| [docs/BUILDING.md](docs/BUILDING.md) | Rust and C artifacts |
| [API reference](https://docs.rs/rssn) | Generated rustdoc |

---

## Contributing

We welcome contributions of all kinds — bug fixes, performance optimizations, new algorithms, and documentation improvements.
See [CONTRIBUTING.md](CONTRIBUTING.md) for detailed guidelines.

---

## License

Licensed under the **Apache 2.0** License.
Please see [LICENSE](LICENSE) for more details.

---

## Attributions

Please see [ATTRIBUTIONS](ATTRIBUTIONS.md) for more details.

---

## Security

Please see [SECURITY](SECURITY.md) for more details.

---

## Code Of Conduct

Please see [CODE_OF_CONDUCT](CODE_OF_CONDUCT.md) for more details.

Reports of abuse are fully available in this project.

---

## Project Wiki

Please see the GitHub wiki page for more details.

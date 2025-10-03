# rssn: A Comprehensive Scientific Computing Library for Rust

[![Crates.io](https://img.shields.io/crates/v/rssn.svg)](https://crates.io/crates/rssn)
[![Docs.rs](https://docs.rs/rssn/badge.svg)](https://docs.rs/rssn)
[![License](https://img.shields.io/crates/l/rssn)](LICENSE)

**rssn** is an open-source scientific computing library for Rust, combining **symbolic computation**, **numerical methods**, and **physics simulations** in a single ecosystem.  
It is designed to provide a foundation for building a **next-generation CAS (Computer Algebra System)** and numerical toolkit in Rust.

---

## ✨ Features

The library is organized into five major components:

- **Symbolic**:  
  Computer algebra system foundations, differentiation & integration, group theory, Lie algebras, polynomial algebra, PDE/ODE solvers, Grobner bases, quantum mechanics operators, graph algorithms, and more.

- **Numerical**:  
  Linear algebra, optimization (Rastrigin, Rosenbrock, Sphere, Linear Regression), numerical integration, probability distributions, FFT, combinatorics, special functions, PDE solvers (heat, wave, Schrödinger 1D–3D), root finding, and statistical analysis.

- **Physics**:  
  Simulation modules covering FDM/FEM/FVM solvers, multigrid methods, molecular mechanics (SPH), electrodynamics (FDTD), Navier–Stokes fluid dynamics, relativity (geodesics, Schwarzschild), elasticity, quantum simulations, and more.

- **Output**:  
  Pretty-printing, LaTeX/Typst export, NumPy-compatible I/O, and plotting utilities (2D/3D surfaces, vector fields, parametric curves).

- **Plugins**:  
  Optional extensions (enabled with the `full` feature).

---

## 🚀 Quick Start

Add **rssn** to your Rust project:

```bash
cargo add rssn
````

Then start exploring:

```rust
use rssn::symbolic::calculus::differentiate;
use rssn::numerical::matrix::Matrix;

fn main() {
    let expr = "x^2 + sin(x)";
    let d = differentiate(expr, "x");
    println!("d/dx of {} = {}", expr, d);

    let m = Matrix::identity(3);
    println!("3x3 Identity Matrix: {:?}", m);
}
```

For more examples, see the [project repository](https://github.com/Apich-Organization/rssn).

---

## 📚 Documentation

* API Docs: [docs.rs/rssn](https://docs.rs/rssn)
* Project Website: [Apich-Organization.github.io/rssn](https://Apich-Organization.github.io/rssn)

---

## 🗺️ Roadmap

* **v0.1.0** — First public release
* **v0.2.0** — Stabilization release
* **v0.3.0** — Performance improvements & broader coverage
* **v0.4.0** — Optional FFI for HPC, start development of **rsst** scripting toolkit
* **v1.0.0** — API stabilization

---

## 🤝 Contributing

We welcome contributions of all kinds — bug fixes, performance optimizations, new algorithms, and documentation improvements.
See [CONTRIBUTING.md](CONTRIBUTING.md) for detailed guidelines.

---

## 💰 Sponsorship & Donations

Scientific computing requires heavy resources for CI/CD, benchmarking, and cloud testing.
You can support development via **GitHub Sponsors**.

Enterprise sponsors will receive:

* Priority support from the core maintainers
* Ability to request features
* Direct collaboration on integration needs

Excess donations will be redirected to upstream Rust ecosystem projects (e.g., rust-LLVM) or community initiatives.

Updates:
Due to temporary issues, GitHub Sponsors is currently unavailable. If you would like to make a donation, please use PayPal to donate to [@panayang338](https://www.paypal.me/panayang338).

---

## 👥 Maintainers & Contributors

* **Author**: [Pana Yang](https://github.com/panayang) (ORCID: 0009-0007-2600-0948, email: [Pana.Yang@hotmail.com](mailto:Pana.Yang@hotmail.com))
* **Consultants**:

  * X. Zhang (Algorithm & Informatics, [@RheaCherry](https://github.com/RheaCherry), [3248998213@qq.com](mailto:3248998213@qq.com))
  * Z. Wang (Mathematics)
  * Y. Li (Physics) ([xian1360685019@qq.com](mailto:xian1360685019@qq.com))
* **Additional contributors**: Owen Yang ([yangguangyong@gmail.com](mailto:yangguangyong@gmail.com))

---

## 📜 License

Licensed under the **Apache 2.0**.
See [LICENSE](LICENSE) for details.



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

**rssn** is an open-source scientific computing library for Rust, combining a high-performance **symbolic computation** engine with **numerical methods** support and **physics simulations** functionalities.

At its core, `rssn` utilizes a **Directed Acyclic Graph (DAG)** to represent mathematical expressions, ensuring that they are always in a canonical form. This allows for highly efficient memory use and computational speed.

---

## Key Features

- **Efficient DAG-based Expression Model**: Expressions are stored as a canonical DAG, ensuring that identical subexpressions are represented by a single node in memory for maximum efficiency.
- **Advanced Symbolic Algebra**: A powerful Computer Algebra System (CAS) that goes beyond simple simplification:
  - **Polynomial Algebra**: Includes **Gröbner basis** computation for solving polynomial systems.
  - **Simplification with Relations**: Can simplify expressions with respect to polynomial side-relations (e.g., simplifying `x^2` to `1 - y^2` given `x^2 + y^2 - 1 = 0`).
- **Symbolic Calculus**: Functions for differentiation, integration, limits, and series expansion.
- **Numerical Methods**: A rich collection of algorithms for numerical integration, optimization, and solving differential equations.
- **Versatile Output**: Render expressions as pretty-printed text, LaTeX, or Typst.
- **Stable FFI Interface**: A robust C-compatible foreign function interface (`cdylib`) is available for integration with other languages like Python, C++, and Fortran.

---

## Quick Start

Add **rssn** to your Rust project:

```bash
cargo add rssn
```

Then, perform a simple symbolic differentiation:

```rust
use rssn::symbolic::core::Expr;
use rssn::symbolic::calculus::differentiate;

// Define a variable 'x'
let x = Expr::new_variable("x");

// Define the expression: sin(x)
let expr = Expr::new_sin(x);

// Differentiate with respect to 'x'
let derivative = differentiate(&expr, "x");

// The result will be cos(x)
println!("The derivative of {} is: {}", expr, derivative);
```

For more advanced examples, such as simplification with relations, please see the [API documentation](https://docs.rs/rssn).

---

## Roadmap

* **v0.1.0** — Finalize core symbolic engine, expand simplification rules.
* **v0.2.0** — Stabilization release, expand numerical methods.
* **v0.3.0** — Performance improvements & broader algorithm coverage.
* **v0.4.0** — Building up Computation Ecosystem, expand physics simulations.
* **v0.5.0** — Optional FFI for HPC, start development of **rsst** scripting toolkit.
* **v1.0.0** — Full API stabilization.

---

## Contributing

We welcome contributions of all kinds — bug fixes, performance optimizations, new algorithms, and documentation improvements.
See [CONTRIBUTING.md](CONTRIBUTING.md) for detailed guidelines.

---

## Maintainers & Contributors

* **Author**: [Pana Yang](https://github.com/panayang) (ORCID: 0009-0007-2600-0948, email: [Pana.Yang@hotmail.com](mailto:Pana.Yang@hotmail.com))
* **Consultants**:
  * X. Zhang (Algorithm & Informatics, [@RheaCherry](https://github.com/RheaCherry), [3248998213@qq.com](mailto:3248998213@qq.com))
  * Z. Wang (Mathematics)
  * Y. Li (Physics) ([xian1360685019@qq.com](mailto:xian1360685019@qq.com))
* **Project Reviewer**: Z. Li
* **Outside Collaborator**: Chahat Patel ([@chahat-101](https://github.com/chahat-101))
* **Contributors**: Owen Yang ([yangguangyong@gmail.com](mailto:yangguangyong@gmail.com))

---

## License

Licensed under the **Apache 2.0** License.
Please see [LICENSE](LICENSE) for more details.

---

## Architecture

Please see [ARCHITECTURE](ARCHITECTURE.md) for more details.

---

## Code Stastics

Please see [CODE_STASTICS](CODE_STASTICS.md) for more details.

---

## Attributions

Please see [ATTRIBUTIONS](ATTRIBUTIONS.md) for more details.

---

## Security

Please see [SECURITY](SECURITY.md) for more details.

---

## Code Of Conduct

Please see [CODE_OF_CONDUCT](CODE_OF_CONDUCT.md) for more details.

Report of abuse are fully avalible in this project.

---

## Project Wiki

Please see the GitHub wiki Page for more details.

---

## A Note from the Author

As one of the primary author, I extend my deepest gratitude for your interest in this project.

I am a high school student in mainland China with an interest in the field of hep-th and computing science. Due to my demanding academic commitments, sometimes my time is limited, and my responses to issues and core pull requests which need my review may sometimes be delayed.

And also, as one of the mission of Apich, we will continue to test the edges of the current AI system assisted coding and development. Discussions on that is welcomed but only without hate. 

I sincerely appreciate all of your patience and understanding, and I welcome any contribution from the community.

--- Pana Yang



# rssn: A Comprehensive Scientific Computing Library for Rust

[![Crates.io](https://img.shields.io/crates/v/rssn.svg)](https://crates.io/crates/rssn)
[![Docs.rs](https://docs.rs/rssn/badge.svg)](https://docs.rs/rssn)
[![License](https://img.shields.io/crates/l/rssn)](LICENSE)

**rssn** is an open-source scientific computing library for Rust, combining **symbolic computation**, **numerical methods**, and **physics simulations** in a single ecosystem.  
It is designed to provide a foundation for building a **next-generation CAS (Computer Algebra System)** and numerical toolkit in Rust.

## Project Status and Engineering Focus

Due to recent community discussions, some of which included unprofessional language, we have decided to **isolate the primary development focus** and move all related architectural discussions to **GitHub Discussions**. We have taken formal steps to address the inappropriate behavior.

Effective immediately, the majority of our resources will be dedicated to the **Dynamic Library (`cdylib`) version** of the core.

### Why the Pivot to FFI?

Our primary commitment is to provide **maximum stability, reliability, and institutional adoption** in high-stakes scientific computing environments (Fortran, C++, Python).

* **Focus:** We are implementing a highly robust **Handle-JSON Hybrid FFI** interface.
* **Goal:** To securely expose the `rssn` core's symbolic analysis capabilities via a stable C interface, ensuring **absolute isolation** from the internal Rust implementation.
* **Commitment:** We continue to validate the core with **property-based testing (`proptest`)** to guarantee professional-grade accuracy and zero failures in complex scenarios.

**Our best response to any doubt is uncompromising engineering quality and reliability.** Thank you for your support as we focus on delivering this critical FFI layer.

## rssn FFI Usage Guide

### Core Concepts

The FFI is built around two core concepts:

1.  **Handles**: Rust objects (like symbolic expressions) are exposed to the C API as opaque pointers called "handles". You can think of a handle as a ticket that refers to an object living in Rust's memory. You can pass these handles back to other FFI functions to operate on the objects they represent.
    - A handle for an `Expr` object is of type `*mut Expr`.

2.  **JSON Serialization**: Complex data is passed across the FFI boundary using JSON strings. For example, to create a symbolic expression, you provide a JSON representation of that expression. Similarly, some functions may return a JSON string to represent a complex result or an error.

### Memory Management

**The caller is responsible for memory management.**

When you create an object via an FFI function (e.g., `expr_from_json`), you receive a handle (a pointer). When you are finished with this handle, you **must** call the corresponding `_free` function (e.g., `expr_free`) to release the memory. Failure to do so will result in memory leaks.

Similarly, when an FFI function returns a string (`*mut c_char`), you **must** call `free_string` to release its memory.

**General Rule:** If you receive a pointer from the library, you own it, and you must free it.

### Basic Workflow

1.  **Create an object**: Use a `_from_json` function to create an object from a JSON string. You will get a handle.
2.  **Operate on the object**: Pass the handle to other FFI functions (e.g., `expr_simplify`, `expr_to_string`).
3.  **Inspect the result**: If a function returns a string (like `expr_to_string` or `expr_to_json`), you can read it. Remember to free it afterwards. If a function returns a new handle, you now own that handle.
4.  **Clean up**: When you are done with a handle, call its `_free` function.

### FFI Health Check

Before diving into complex operations, it is a good practice to verify that the FFI interface is working correctly. The following function is provided for this purpose.

- `rssn_test_string_passing() -> *mut c_char`
  This function allocates a simple test string ("pong") and returns a pointer to it. It serves two purposes:
  1.  Confirms that you can successfully call a function in the `rssn` library.
  2.  Allows you to test the memory management of strings. You should call `free_string` on the returned pointer to ensure that allocation and deallocation are working correctly across the FFI boundary.

**Example Verification Flow:**
1. Call `rssn_test_string_passing()` and receive a pointer.
2. Check if the pointer is not null.
3. (Optional) Read the string to verify it is "pong".
4. Call `free_string()` on the pointer.

If all these steps complete without errors, your FFI setup is likely correct.

### Available Functions for `Expr`

Below is a summary of the available FFI functions for `Expr` objects.

1. Object Creation and Destruction

- `expr_from_json(json_ptr: *const c_char) -> *mut Expr`
  Creates an `Expr` object from a JSON string. Returns a handle to the new object. Returns a null pointer if the JSON is invalid.

- `expr_to_json(handle: *mut Expr) -> *mut c_char`
  Serializes the `Expr` object pointed to by the handle into a JSON string. The caller must free the returned string.

- `expr_free(handle: *mut Expr)`
  Frees the memory of the `Expr` object associated with the handle.

2. Expression Operations

- `expr_to_string(handle: *mut Expr) -> *mut c_char`
  Returns a human-readable string representation of the expression. The caller must free the returned string.

- `expr_simplify(handle: *mut Expr) -> *mut Expr`
  Simplifies the expression and returns a handle to a **new** simplified expression. The caller owns the new handle and must free it.

- `expr_unify_expression(handle: *mut Expr) -> *mut c_char`
  Attempts to unify the physical units within an expression. This function returns a JSON string representing a result object. The result object will have one of two fields:
    - `ok`: If successful, this field will contain the JSON representation of the new, unified `Expr`. You can pass this JSON to `expr_from_json` to get a handle to it.
    - `err`: If it fails, this field will contain a string with the error message.

### Utility Functions

- `free_string(s: *mut c_char)`
  Frees a string that was allocated and returned by the library.

## Example `Expr` JSON Format

The JSON format for an `Expr` directly mirrors the Rust enum definition. Here are a few examples:

**A simple constant `3.14`:**
```json
{ "Constant": 3.14 }
```

**A variable `x`:**
```json
{ "Variable": "x" }
```

**The expression `x + 2`:**
```json
{
  "Add": [
    { "Variable": "x" },
    { "Constant": 2.0 }
  ]
}
```

**The expression `sin(x^2)`:**
```json
{
  "Sin": {
    "Power": [
      { "Variable": "x" },
      { "Constant": 2.0 }
    ]
  }
}
```


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
use num_bigint::BigInt;
use rssn::symbolic::calculus::differentiate;
use rssn::symbolic::core::Expr;
use std::sync::Arc;

fn test_differentiate_x_squared_stack_overflow() {
    let x = Expr::Variable("x".to_string());
    let x2 = Expr::Mul(Arc::new(x.clone()), Arc::new(x.clone()));
    let d = differentiate(&x2, "x");

    // The derivative of x^2 is 2*x.
    // The simplification process might result in Constant(2.0) or BigInt(2).
    let two_const = Expr::Constant(2.0);
    let expected_const = Expr::Mul(Arc::new(two_const), Arc::new(x.clone()));

    let two_int = Expr::BigInt(BigInt::from(2));
    let expected_int = Expr::Mul(Arc::new(two_int), Arc::new(x.clone()));

    println!("Derivative: {:?}", d);
    println!("Expected (const): {:?}", expected_const);
    println!("Expected (int): {:?}", expected_int);

    assert!(d == expected_const || d == expected_int);
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

## Code statistics (SCC 2025-10-08 CST)

```
───────────────────────────────────────────────────────────────────────────────
Language                 Files     Lines   Blanks  Comments     Code Complexity
───────────────────────────────────────────────────────────────────────────────
Rust                       417     60046     5743     18990    35313       4546
Markdown                     4      1211      181         0     1030          0
TOML                         2       113       11         0      102          2
YAML                         2       304       52        17      235          0
Fortran Modern               1       167       34        33      100          6
License                      1        73       32         0       41          0
───────────────────────────────────────────────────────────────────────────────
Total                      427     61914     6053     19040    36821       4554
───────────────────────────────────────────────────────────────────────────────
Estimated Cost to Develop (organic) $1,191,230
Estimated Schedule Effort (organic) 14.70 months
Estimated People Required (organic) 7.20
───────────────────────────────────────────────────────────────────────────────
Processed 2292078 bytes, 2.292 megabytes (SI)
───────────────────────────────────────────────────────────────────────────────

───────────────────────────────────────────────────────────────────────────────
Language                 Files     Lines   Blanks  Comments     Code Complexity
───────────────────────────────────────────────────────────────────────────────
Rust                       417     60046     5743     18990    35313       4546
───────────────────────────────────────────────────────────────────────────────
src\ffi_api.rs                      3334      208       349     2777        375
src\symbolic\core.rs                2300       66       242     1992        136
src\symbolic\calculus.rs            2154      122       365     1667        331
~\symbolic\graph_algorithms.rs      1632      124       346     1162        288
src\symbolic\polynomial.rs          1305       79       514      712        130
src\symbolic\simplify.rs            1203       59       119     1025        197
src\symbolic\matrix.rs              1070       88       309      673        175
src\symbolic\solve.rs               1004       79       156      769        135
src\numerical\testing.rs             985       74       218      693        123
src\symbolic\ode.rs                  952       98       164      690        111
src\symbolic\pde.rs                  947       74       148      725        101
src\symbolic\coordinates.rs          820       56       199      565         60
~ymbolic\integral_equations.rs       754       51       437      266          7
src\numerical\physics.rs             752       61       361      330         65
~symbolic\computer_graphics.rs       686       35       172      479          7
src\symbolic\integration.rs          684       74       220      390         35
src\symbolic\combinatorics.rs        670       61       259      350         48
src\symbolic\number_theory.rs        661       65       180      416        108
src\numerical\optimize.rs            644       41       442      161         27
~symbolic\special_functions.rs       640       24       223      393         64
src\symbolic\logic.rs                578       34        91      453         88
~ymbolic\poly_factorization.rs       572       47       155      370         58
src\symbolic\transforms.rs           558       28       193      337         39
src\symbolic\tensor.rs               551       37       179      335         63
src\symbolic\topology.rs             545       42       207      296         55
~ic\error_correction_helper.rs       522       51       142      329         57
~symbolic\stats_probability.rs       520       35       172      313          0
src\numerical\matrix.rs              483       43       113      327         65
src\physics\physics_fvm.rs           455       51        83      321         38
src\numerical\stats.rs               452       58       156      238         12
src\symbolic\finite_field.rs         441       43        89      309         52
src\symbolic\series.rs               420       20       133      267         31
~olic\differential_geometry.rs       416       22       239      155         20
src\physics\physics_mtm.rs           364       42        42      280         49
src\physics\physics_fem.rs           363       44        55      264         54
src\physics\physics_sm.rs            356       48        52      256         12
src\output\pretty_print.rs           349       27        64      258         28
~c\symbolic\cas_foundations.rs       347       22        15      310         56
src\physics\physics_rkm.rs           328       32        65      231         14
src\physics\physics_em.rs            327       74       100      153         10
src\symbolic\convergence.rs          325       15       189      121         34
src\symbolic\vector.rs               322       20       124      178          3
src\symbolic\elementary.rs           306       15        68      223         11
~\symbolic\unit_unification.rs       291       27        37      227         25
~symbolic\geometric_algebra.rs       287       18        86      183         38
src\output\plotting.rs               285       33        15      237          9
~ic\lie_groups_and_algebras.rs       275       26        81      168          7
~c\symbolic\discrete_groups.rs       271       22       113      136         25
src\symbolic\proof.rs                268       22        87      159         33
src\numerical\interpolate.rs         263       31        61      171         41
src\numerical\finite_field.rs        261       24        86      151         28
src\symbolic\cad.rs                  260       27        21      212         35
src\symbolic\rewriting.rs            259       21        49      189         46
~numerical\error_correction.rs       258       33        60      165         27
src\symbolic\graph.rs                258       21       110      127         15
src\physics\physics_bem.rs           250       23        55      172         17
~\symbolic\error_correction.rs       248       37        62      149         27
src\physics\physics_mm.rs            243       28        22      193         21
src\symbolic\grobner.rs              243       29        45      169         19
~mbolic\functional_analysis.rs       230       19       110      101          2
~c\stats_information_theory.rs       229       18        90      121         12
src\numerical\sparse.rs              227       23        69      135         18
src\physics\physics_fdm.rs           225       24        48      153         26
src\symbolic\real_roots.rs           222       25        52      145         28
src\numerical\polynomial.rs          219       20        46      153         30
src\output\latex.rs                  208       14        12      182          7
src\physics\physics_cnm.rs           201       24        41      136         17
~sics_sim\linear_elasticity.rs       198       27        36      135         20
~umerical\computer_graphics.rs       197       16        75      106          7
src\plugins\manager.rs               193       26        25      142         15
src\numerical\coordinates.rs         190       11        68      111         16
~h_isomorphism_and_coloring.rs       189       17        54      118         29
~rc\symbolic\thermodynamics.rs       189       11       100       78          0
~c\symbolic\vector_calculus.rs       182       21        62       99          0
~\symbolic\complex_analysis.rs       172       18        72       82          8
src\symbolic\optimize.rs             171       17        54      100         18
~mbolic\classical_mechanics.rs       171       13        99       59          0
~\numerical\vector_calculus.rs       169       26        52       91          7
src\numerical\transforms.rs          168       26        57       85         19
~symbolic\quantum_mechanics.rs       164       15        78       71          0
src\symbolic\group_theory.rs         160       11        63       86         14
~cs_sim\schrodinger_quantum.rs       159       20        21      118         14
src\numerical\solve.rs               155       16        41       98         20
src\symbolic\relativity.rs           155       15        75       65          0
~cs_sim\geodesic_relativity.rs       151       21        39       91          3
src\symbolic\cryptography.rs         150       12        60       78         13
src\lib.rs                           150        7        92       51          1
src\symbolic\multi_valued.rs         149       20        67       62          0
src\numerical\tensor.rs              148       16        46       86         15
~umerical\geometric_algebra.rs       148        7        27      114          4
src\numerical\physics_cfd.rs         146       13        48       85         11
~\symbolic\graph_operations.rs       145       15        10      120         41
src\numerical\real_roots.rs          145       13        29      103         43
~\symbolic\stats_regression.rs       143       16        53       74          5
~cs_sim\navier_stokes_fluid.rs       143       19        20      104         17
src\symbolic\numeric.rs              136       10        31       95         32
src\symbolic\stats.rs                134        7        59       68          7
~bolic\quantum_field_theory.rs       133       11        68       54          2
~lic\calculus_of_variations.rs       131       14        88       29          0
src\numerical\graph.rs               128       15        45       68          7
src\numerical\topology.rs            127       14        29       84         19
~rc\numerical\number_theory.rs       127       11        44       72         33
~\symbolic\electromagnetism.rs       126       13        56       57          0
src\numerical\vector.rs              125        9        68       48         14
src\output\io.rs                     123       13        44       66          9
~sics_sim\gpe_superfluidity.rs       118       17        21       80          5
src\numerical\convergence.rs         117        7        41       69          9
~sics_sim\ising_statistical.rs       115       16        23       76         11
src\numerical\physics_fea.rs         114       10        46       58         10
src\numerical\physics_md.rs          113       16        41       56          4
src\numerical\integrate.rs           112        8        20       84         10
~s_sim\fdtd_electrodynamics.rs       109       16        21       72         16
~rc\numerical\combinatorics.rs       106        8        40       58         15
src\numerical\ode.rs                 102       12        20       70          4
~erical\functional_analysis.rs       100        5        50       45          4
~c\symbolic\stats_inference.rs        99        9        22       68          0
src\numerical\signal.rs               97       14        26       57          4
src\output\typst.rs                   96        4         2       90          3
~ical\differential_geometry.rs        95        9        26       60         12
src\symbolic\special.rs               91        8        61       22          0
src\symbolic\mod.rs                   89        3        19       67          0
src\symbolic\radicals.rs              87        7        21       59         10
~fractal_geometry_and_chaos.rs        83        7        31       45          5
~mbolic\solid_state_physics.rs        77        6        42       29          0
tests\numerical\calculus.rs           76       14        11       51          3
src\plugins\mod.rs                    75       13        28       34          2
~ests\numerical\interpolate.rs        73       12        13       48          7
~fractal_geometry_and_chaos.rs        68        8        35       25          1
src\numerical\calculus.rs             68       10        29       29          4
~numerical\complex_analysis.rs        64        3        16       45          1
src\numerical\multi_valued.rs         63        9        19       35          5
src\symbolic\handles.rs               57        8        16       33          0
src\numerical\series.rs               56        7        17       32          2
tests\symbolic\calculus.rs            54        7        25       22          0
src\numerical\elementary.rs           54        3        16       35          1
~cal\calculus_of_variations.rs        49        5        26       18          0
tests\mod.rs                          49        6        42        1          0
~s_sim\fdtd_electrodynamics.rs        48        6        25       17          0
~mbolic\solid_state_physics.rs        48        6        25       17          0
tests\numerical\elementary.rs         48        6        25       17          0
tests\output\pretty_print.rs          48        6        25       17          0
~ical\differential_geometry.rs        48        6        25       17          0
~numerical\error_correction.rs        48        6        25       17          0
tests\output\typst.rs                 48        6        25       17          0
~sts\numerical\finite_field.rs        48        6        25       17          0
~ts\numerical\combinatorics.rs        48        6        25       17          0
tests\physics\physics_bem.rs          48        6        25       17          0
~umerical\computer_graphics.rs        48        6        25       17          0
~numerical\complex_analysis.rs        48        6        25       17          0
~cal\calculus_of_variations.rs        48        6        25       17          0
tests\physics\mod.rs                  48        6        25       17          0
~ests\numerical\coordinates.rs        48        6        25       17          0
tests\physics\physics_fdm.rs          48        6        25       17          0
~ests\numerical\convergence.rs        48        6        25       17          0
~erical\functional_analysis.rs        48        6        25       17          0
tests\physics\physics_cnm.rs          48        6        25       17          0
tests\numerical\graph.rs              48        6        25       17          0
~umerical\geometric_algebra.rs        48        6        25       17          0
~ts\numerical\number_theory.rs        48        6        25       17          0
tests\numerical\matrix.rs             48        6        25       17          0
tests\output\plotting.rs              48        6        25       17          0
~fractal_geometry_and_chaos.rs        48        6        25       17          0
tests\physics\physics_em.rs           48        6        25       17          0
tests\numerical\optimize.rs           48        6        25       17          0
~sts\numerical\multi_valued.rs        48        6        25       17          0
tests\numerical\pde.rs                48        6        25       17          0
~ests\numerical\physics_cfd.rs        48        6        25       17          0
tests\numerical\physics.rs            48        6        25       17          0
~ests\numerical\physics_fea.rs        48        6        25       17          0
tests\numerical\polynomial.rs         48        6        25       17          0
tests\numerical\physics_md.rs         48        6        25       17          0
tests\numerical\integrate.rs          48        6        25       17          0
tests\numerical\special.rs            48        6        25       17          0
tests\numerical\solve.rs              48        6        25       17          0
tests\numerical\real_roots.rs         48        6        25       17          0
tests\numerical\tensor.rs             48        6        25       17          0
tests\physics\physics_fem.rs          48        6        25       17          0
tests\numerical\sparse.rs             48        6        25       17          0
tests\numerical\series.rs             48        6        25       17          0
tests\physics\physics_fvm.rs          48        6        25       17          0
~s\symbolic\vector_calculus.rs        48        6        25       17          0
tests\numerical\testing.rs            48        6        25       17          0
tests\symbolic\topology.rs            48        6        25       17          0
tests\output\mod.rs                   48        6        25       17          0
tests\output\latex.rs                 48        6        25       17          0
tests\symbolic\vector.rs              48        6        25       17          0
tests\numerical\signal.rs             48        6        25       17          0
tests\output\io.rs                    48        6        25       17          0
tests\numerical\stats.rs              48        6        25       17          0
tests\symbolic\transforms.rs          48        6        25       17          0
~symbolic\stats_probability.rs        48        6        25       17          0
~\numerical\vector_calculus.rs        48        6        25       17          0
tests\symbolic\stats.rs               48        6        25       17          0
tests\symbolic\tensor.rs              48        6        25       17          0
~ts\symbolic\thermodynamics.rs        48        6        25       17          0
~s\symbolic\stats_inference.rs        48        6        25       17          0
tests\symbolic\special.rs             48        6        25       17          0
~symbolic\special_functions.rs        48        6        25       17          0
~\symbolic\stats_regression.rs        48        6        25       17          0
~c\stats_information_theory.rs        48        6        25       17          0
tests\physics\physics_mtm.rs          48        6        25       17          0
tests\physics\physics_mm.rs           48        6        25       17          0
~symbolic\geometric_algebra.rs        48        6        25       17          0
tests\physics\physics_sm.rs           48        6        25       17          0
~mbolic\functional_analysis.rs        48        6        25       17          0
~\symbolic\graph_algorithms.rs        48        6        25       17          0
tests\numerical\topology.rs           48        6        25       17          0
tests\physics\physics_rkm.rs          48        6        25       17          0
~sics_sim\gpe_superfluidity.rs        48        6        25       17          0
~sics_sim\linear_elasticity.rs        48        6        25       17          0
tests\symbolic\simplify.rs            48        6        25       17          0
~cs_sim\geodesic_relativity.rs        48        6        25       17          0
~cs_sim\schrodinger_quantum.rs        48        6        25       17          0
tests\symbolic\solve.rs               48        6        25       17          0
~sics_sim\ising_statistical.rs        48        6        25       17          0
~ts\physics\physics_sim\mod.rs        48        6        25       17          0
~cs_sim\navier_stokes_fluid.rs        48        6        25       17          0
~sts\symbolic\combinatorics.rs        48        6        25       17          0
tests\symbolic\cad.rs                 48        6        25       17          0
~s\symbolic\cas_foundations.rs        48        6        25       17          0
~\symbolic\complex_analysis.rs        48        6        25       17          0
~mbolic\classical_mechanics.rs        48        6        25       17          0
tests\symbolic\series.rs              48        6        25       17          0
~symbolic\computer_graphics.rs        48        6        25       17          0
tests\symbolic\rewriting.rs           48        6        25       17          0
~lic\calculus_of_variations.rs        48        6        25       17          0
tests\symbolic\convergence.rs         48        6        25       17          0
~ests\symbolic\cryptography.rs        48        6        25       17          0
~olic\differential_geometry.rs        48        6        25       17          0
tests\symbolic\coordinates.rs         48        6        25       17          0
tests\numerical\transforms.rs         48        6        25       17          0
~symbolic\quantum_mechanics.rs        48        6        25       17          0
~\symbolic\electromagnetism.rs        48        6        25       17          0
~s\symbolic\discrete_groups.rs        48        6        25       17          0
~ests\symbolic\finite_field.rs        48        6        25       17          0
~fractal_geometry_and_chaos.rs        48        6        25       17          0
tests\symbolic\elementary.rs          48        6        25       17          0
~ic\error_correction_helper.rs        48        6        25       17          0
~\symbolic\error_correction.rs        48        6        25       17          0
tests\symbolic\relativity.rs          48        6        25       17          0
~ymbolic\poly_factorization.rs        48        6        25       17          0
tests\symbolic\real_roots.rs          48        6        25       17          0
tests\symbolic\radicals.rs            48        6        25       17          0
~ic\lie_groups_and_algebras.rs        48        6        25       17          0
tests\symbolic\optimize.rs            48        6        25       17          0
~bolic\quantum_field_theory.rs        48        6        25       17          0
tests\symbolic\proof.rs               48        6        25       17          0
tests\symbolic\matrix.rs              48        6        25       17          0
~ests\symbolic\multi_valued.rs        48        6        25       17          0
~h_isomorphism_and_coloring.rs        48        6        25       17          0
tests\symbolic\polynomial.rs          48        6        25       17          0
~sts\symbolic\number_theory.rs        48        6        25       17          0
tests\symbolic\numeric.rs             48        6        25       17          0
tests\symbolic\ode.rs                 48        6        25       17          0
tests\symbolic\pde.rs                 48        6        25       17          0
tests\numerical\vector.rs             48        6        25       17          0
tests\symbolic\integration.rs         48        6        25       17          0
tests\symbolic\graph.rs               48        6        25       17          0
~ymbolic\integral_equations.rs        48        6        25       17          0
tests\symbolic\logic.rs               48        6        25       17          0
~ests\symbolic\group_theory.rs        48        6        25       17          0
tests\symbolic\grobner.rs             48        6        25       17          0
~\symbolic\graph_operations.rs        48        6        25       17          0
tests\prelude.rs                      47        5        42        0          0
~ins\example_plugin\src\lib.rs        47        7         9       31          1
src\numerical\mod.rs                  46        1         4       41          0
src\physics\mod.rs                    46        1        34       11          0
~enches\symbolic\relativity.rs        43        8        23       12          0
~enches\symbolic\transforms.rs        43        8        23       12          0
benches\symbolic\matrix.rs            43        8        23       12          0
~ches\symbolic\multi_valued.rs        43        8        23       12          0
benches\symbolic\mod.rs               43        8        23       12          0
benches\symbolic\pde.rs               43        8        23       12          0
benches\symbolic\series.rs            43        8        23       12          0
~bolic\quantum_field_theory.rs        43        8        23       12          0
benches\symbolic\optimize.rs          43        8        23       12          0
~sics_sim\linear_elasticity.rs        43        8        23       12          0
benches\symbolic\rewriting.rs         43        8        23       12          0
~enches\symbolic\real_roots.rs        43        8        23       12          0
benches\symbolic\numeric.rs           43        8        23       12          0
benches\symbolic\proof.rs             43        8        23       12          0
benches\symbolic\ode.rs               43        8        23       12          0
~enches\symbolic\polynomial.rs        43        8        23       12          0
~symbolic\quantum_mechanics.rs        43        8        23       12          0
benches\symbolic\simplify.rs          43        8        23       12          0
~ymbolic\poly_factorization.rs        43        8        23       12          0
benches\symbolic\logic.rs             43        8        23       12          0
benches\symbolic\solve.rs             43        8        23       12          0
~mbolic\solid_state_physics.rs        43        8        23       12          0
~nches\symbolic\integration.rs        43        8        23       12          0
~ymbolic\integral_equations.rs        43        8        23       12          0
benches\symbolic\radicals.rs          43        8        23       12          0
~es\numerical\number_theory.rs        43        8        23       12          0
~ic\lie_groups_and_algebras.rs        43        8        23       12          0
~enches\physics\physics_cnm.rs        43        8        23       12          0
benches\symbolic\grobner.rs           43        8        23       12          0
~ches\symbolic\group_theory.rs        43        8        23       12          0
~nches\numerical\physics_md.rs        43        8        23       12          0
~symbolic\computer_graphics.rs        43        8        23       12          0
~\symbolic\graph_operations.rs        43        8        23       12          0
~h_isomorphism_and_coloring.rs        43        8        23       12          0
~\symbolic\error_correction.rs        43        8        23       12          0
benches\numerical\matrix.rs           43        8        23       12          0
~\symbolic\electromagnetism.rs        43        8        23       12          0
~symbolic\geometric_algebra.rs        43        8        23       12          0
~\symbolic\graph_algorithms.rs        43        8        23       12          0
benches\numerical\physics.rs          43        8        23       12          0
~mbolic\functional_analysis.rs        43        8        23       12          0
benches\symbolic\graph.rs             43        8        23       12          0
benches\numerical\optimize.rs         43        8        23       12          0
~fractal_geometry_and_chaos.rs        43        8        23       12          0
~ic\error_correction_helper.rs        43        8        23       12          0
benches\numerical\mod.rs              43        8        23       12          0
benches\numerical\pde.rs              43        8        23       12          0
~hes\numerical\multi_valued.rs        43        8        23       12          0
~ches\symbolic\finite_field.rs        43        8        23       12          0
~olic\differential_geometry.rs        43        8        23       12          0
~enches\symbolic\elementary.rs        43        8        23       12          0
~s\symbolic\discrete_groups.rs        43        8        23       12          0
benches\symbolic\special.rs           43        8        23       12          0
~ches\numerical\physics_fea.rs        43        8        23       12          0
~ches\symbolic\cryptography.rs        43        8        23       12          0
benches\symbolic\stats.rs             43        8        23       12          0
benches\symbolic\core.rs              43        8        23       12          0
~nches\symbolic\coordinates.rs        43        8        23       12          0
~nches\symbolic\convergence.rs        43        8        23       12          0
~c\stats_information_theory.rs        43        8        23       12          0
~symbolic\special_functions.rs        43        8        23       12          0
~symbolic\stats_probability.rs        43        8        23       12          0
~s\symbolic\stats_inference.rs        43        8        23       12          0
~\symbolic\complex_analysis.rs        43        8        23       12          0
~ches\numerical\interpolate.rs        43        8        23       12          0
~enches\physics\physics_rkm.rs        43        8        23       12          0
benches\numerical\ode.rs              43        8        23       12          0
~nches\numerical\polynomial.rs        43        8        23       12          0
benches\numerical\signal.rs           43        8        23       12          0
~\symbolic\stats_regression.rs        43        8        23       12          0
~sics_sim\gpe_superfluidity.rs        43        8        23       12          0
benches\symbolic\tensor.rs            43        8        23       12          0
benches\symbolic\topology.rs          43        8        23       12          0
~ches\numerical\physics_cfd.rs        43        8        23       12          0
benches\numerical\solve.rs            43        8        23       12          0
~hes\symbolic\number_theory.rs        43        8        23       12          0
~nches\numerical\real_roots.rs        43        8        23       12          0
benches\numerical\sparse.rs           43        8        23       12          0
benches\symbolic\vector.rs            43        8        23       12          0
~enches\numerical\integrate.rs        43        8        23       12          0
benches\numerical\series.rs           43        8        23       12          0
~hes\symbolic\combinatorics.rs        43        8        23       12          0
benches\numerical\special.rs          43        8        23       12          0
~hes\numerical\finite_field.rs        43        8        23       12          0
~s\symbolic\vector_calculus.rs        43        8        23       12          0
benches\numerical\tensor.rs           43        8        23       12          0
benches\numerical\topology.rs         43        8        23       12          0
benches\numerical\testing.rs          43        8        23       12          0
~cs_sim\schrodinger_quantum.rs        43        8        23       12          0
~nches\numerical\transforms.rs        43        8        23       12          0
benches\mod.rs                        43        8        23       12          0
~umerical\geometric_algebra.rs        43        8        23       12          0
benches\numerical\vector.rs           43        8        23       12          0
~\numerical\vector_calculus.rs        43        8        23       12          0
~fractal_geometry_and_chaos.rs        43        8        23       12          0
benches\numerical\graph.rs            43        8        23       12          0
~erical\functional_analysis.rs        43        8        23       12          0
~es\symbolic\thermodynamics.rs        43        8        23       12          0
benches\physics_sim.rs                43        8        23       12          0
~nches\numerical\elementary.rs        43        8        23       12          0
~numerical\error_correction.rs        43        8        23       12          0
~ical\differential_geometry.rs        43        8        23       12          0
~ches\numerical\coordinates.rs        43        8        23       12          0
~ches\numerical\convergence.rs        43        8        23       12          0
benches\output\io.rs                  43        8        23       12          0
benches\numerical\stats.rs            43        8        23       12          0
~umerical\computer_graphics.rs        43        8        23       12          0
~cs_sim\navier_stokes_fluid.rs        43        8        23       12          0
~numerical\complex_analysis.rs        43        8        23       12          0
benches\output\latex.rs               43        8        23       12          0
benches\output\mod.rs                 43        8        23       12          0
benches\output\typst.rs               43        8        23       12          0
benches\physics\physics_em.rs         43        8        23       12          0
~mbolic\classical_mechanics.rs        43        8        23       12          0
~enches\output\pretty_print.rs        43        8        23       12          0
benches\output\plotting.rs            43        8        23       12          0
benches\physics\mod.rs                43        8        23       12          0
~s\symbolic\cas_foundations.rs        43        8        23       12          0
benches\physics\physics_sm.rs         43        8        23       12          0
~cal\calculus_of_variations.rs        43        8        23       12          0
~es\numerical\combinatorics.rs        43        8        23       12          0
benches\symbolic\calculus.rs          43        8        23       12          0
~enches\physics\physics_bem.rs        43        8        23       12          0
~lic\calculus_of_variations.rs        43        8        23       12          0
~cs_sim\geodesic_relativity.rs        43        8        23       12          0
~es\physics\physics_sim\mod.rs        43        8        23       12          0
~enches\physics\physics_fem.rs        43        8        23       12          0
benches\plugins\mod.rs                43        8        23       12          0
~sics_sim\ising_statistical.rs        43        8        23       12          0
~enches\physics\physics_fvm.rs        43        8        23       12          0
benches\numerical\calculus.rs         43        8        23       12          0
benches\prelude.rs                    43        8        23       12          0
benches\physics\physics_mm.rs         43        8        23       12          0
~enches\physics\physics_mtm.rs        43        8        23       12          0
~enches\physics\physics_fdm.rs        43        8        23       12          0
~s_sim\fdtd_electrodynamics.rs        43        8        23       12          0
benches\symbolic\cad.rs               43        8        23       12          0
tests\symbolic\core.rs                36        6        12       18          1
src\numerical\special.rs              35        7         6       22          0
~\symbolic\unit_unification.rs        28        4         1       23          5
src\prelude.rs                        28        3        19        6          0
tests\regression_test.rs              25        5         2       18          3
src\numerical\pde.rs                  23        3        15        5          0
src\output\mod.rs                     22        1        16        5          0
tests\numerical\ode.rs                17        1         0       16          0
~rc\physics\physics_sim\mod.rs         7        0         0        7          0
tests\numerical\mod.rs                 2        0         0        2          0
tests\symbolic\mod.rs                  2        0         0        2          0
benches\rssn_benches.rs                2        0         1        1          0
tests\plugins\mod.rs                   0        0         0        0          0
───────────────────────────────────────────────────────────────────────────────
Markdown                     4      1211      181         0     1030          0
───────────────────────────────────────────────────────────────────────────────
README.md                            914       95         0      819          0
CONTRIBUTING.md                      203       54         0      149          0
~\ISSUE_TEMPLATE\bug_report.md        49       15         0       34          0
~E_TEMPLATE\feature_request.md        45       17         0       28          0
───────────────────────────────────────────────────────────────────────────────
TOML                         2       113       11         0      102          2
───────────────────────────────────────────────────────────────────────────────
Cargo.toml                           103        9         0       94          2
~ins\example_plugin\Cargo.toml        10        2         0        8          0
───────────────────────────────────────────────────────────────────────────────
YAML                         2       304       52        17      235          0
───────────────────────────────────────────────────────────────────────────────
.gitea\workflows\ci.yml              161       29         1      131          0
.github\workflows\ci.yml             143       23        16      104          0
───────────────────────────────────────────────────────────────────────────────
Fortran Modern               1       167       34        33      100          6
───────────────────────────────────────────────────────────────────────────────
examples\fortran\test.f90            167       34        33      100          6
───────────────────────────────────────────────────────────────────────────────
License                      1        73       32         0       41          0
───────────────────────────────────────────────────────────────────────────────
LICENSE                               73       32         0       41          0
───────────────────────────────────────────────────────────────────────────────
Total                      427     61914     6053     19040    36821       4554
───────────────────────────────────────────────────────────────────────────────
Estimated Cost to Develop (organic) $1,191,230
Estimated Schedule Effort (organic) 14.70 months
Estimated People Required (organic) 7.20
───────────────────────────────────────────────────────────────────────────────
Processed 2292078 bytes, 2.292 megabytes (SI)
───────────────────────────────────────────────────────────────────────────────
```

```
───────────────────────────────────────────────────────────────────────────────
Language                 Files     Lines   Blanks  Comments     Code Complexity
───────────────────────────────────────────────────────────────────────────────
Rust                       140     47649     3824     12531    31294       4526
───────────────────────────────────────────────────────────────────────────────
Total                      140     47649     3824     12531    31294       4526
───────────────────────────────────────────────────────────────────────────────
Estimated Cost to Develop (organic) $1,004,221
Estimated Schedule Effort (organic) 13.78 months
Estimated People Required (organic) 6.48
───────────────────────────────────────────────────────────────────────────────
Processed 1703929 bytes, 1.704 megabytes (SI)
───────────────────────────────────────────────────────────────────────────────

───────────────────────────────────────────────────────────────────────────────
Language                 Files     Lines   Blanks  Comments     Code Complexity
───────────────────────────────────────────────────────────────────────────────
Rust                       140     47649     3824     12531    31294       4526
───────────────────────────────────────────────────────────────────────────────
ffi_api.rs                          3334      208       349     2777        375
symbolic\core.rs                    2300       66       242     1992        136
symbolic\calculus.rs                2154      122       365     1667        331
symbolic\graph_algorithms.rs        1632      124       346     1162        288
symbolic\polynomial.rs              1305       79       514      712        130
symbolic\simplify.rs                1203       59       119     1025        197
symbolic\matrix.rs                  1070       88       309      673        175
symbolic\solve.rs                   1004       79       156      769        135
numerical\testing.rs                 985       74       218      693        123
symbolic\ode.rs                      952       98       164      690        111
symbolic\pde.rs                      947       74       148      725        101
symbolic\coordinates.rs              820       56       199      565         60
~ymbolic\integral_equations.rs       754       51       437      266          7
numerical\physics.rs                 752       61       361      330         65
symbolic\computer_graphics.rs        686       35       172      479          7
symbolic\integration.rs              684       74       220      390         35
symbolic\combinatorics.rs            670       61       259      350         48
symbolic\number_theory.rs            661       65       180      416        108
numerical\optimize.rs                644       41       442      161         27
symbolic\special_functions.rs        640       24       223      393         64
symbolic\logic.rs                    578       34        91      453         88
~ymbolic\poly_factorization.rs       572       47       155      370         58
symbolic\transforms.rs               558       28       193      337         39
symbolic\tensor.rs                   551       37       179      335         63
symbolic\topology.rs                 545       42       207      296         55
~ic\error_correction_helper.rs       522       51       142      329         57
symbolic\stats_probability.rs        520       35       172      313          0
numerical\matrix.rs                  483       43       113      327         65
physics\physics_fvm.rs               455       51        83      321         38
numerical\stats.rs                   452       58       156      238         12
symbolic\finite_field.rs             441       43        89      309         52
symbolic\series.rs                   420       20       133      267         31
~olic\differential_geometry.rs       416       22       239      155         20
physics\physics_mtm.rs               364       42        42      280         49
physics\physics_fem.rs               363       44        55      264         54
physics\physics_sm.rs                356       48        52      256         12
output\pretty_print.rs               349       27        64      258         28
symbolic\cas_foundations.rs          347       22        15      310         56
physics\physics_rkm.rs               328       32        65      231         14
physics\physics_em.rs                327       74       100      153         10
symbolic\convergence.rs              325       15       189      121         34
symbolic\vector.rs                   322       20       124      178          3
symbolic\elementary.rs               306       15        68      223         11
symbolic\unit_unification.rs         291       27        37      227         25
symbolic\geometric_algebra.rs        287       18        86      183         38
output\plotting.rs                   285       33        15      237          9
~ic\lie_groups_and_algebras.rs       275       26        81      168          7
symbolic\discrete_groups.rs          271       22       113      136         25
symbolic\proof.rs                    268       22        87      159         33
numerical\interpolate.rs             263       31        61      171         41
numerical\finite_field.rs            261       24        86      151         28
symbolic\cad.rs                      260       27        21      212         35
symbolic\rewriting.rs                259       21        49      189         46
numerical\error_correction.rs        258       33        60      165         27
symbolic\graph.rs                    258       21       110      127         15
physics\physics_bem.rs               250       23        55      172         17
symbolic\error_correction.rs         248       37        62      149         27
physics\physics_mm.rs                243       28        22      193         21
symbolic\grobner.rs                  243       29        45      169         19
~mbolic\functional_analysis.rs       230       19       110      101          2
~c\stats_information_theory.rs       229       18        90      121         12
numerical\sparse.rs                  227       23        69      135         18
physics\physics_fdm.rs               225       24        48      153         26
symbolic\real_roots.rs               222       25        52      145         28
numerical\polynomial.rs              219       20        46      153         30
output\latex.rs                      208       14        12      182          7
physics\physics_cnm.rs               201       24        41      136         17
~sics_sim\linear_elasticity.rs       198       27        36      135         20
~umerical\computer_graphics.rs       197       16        75      106          7
plugins\manager.rs                   193       26        25      142         15
numerical\coordinates.rs             190       11        68      111         16
symbolic\thermodynamics.rs           189       11       100       78          0
~h_isomorphism_and_coloring.rs       189       17        54      118         29
symbolic\vector_calculus.rs          182       21        62       99          0
symbolic\complex_analysis.rs         172       18        72       82          8
~mbolic\classical_mechanics.rs       171       13        99       59          0
symbolic\optimize.rs                 171       17        54      100         18
numerical\vector_calculus.rs         169       26        52       91          7
numerical\transforms.rs              168       26        57       85         19
symbolic\quantum_mechanics.rs        164       15        78       71          0
symbolic\group_theory.rs             160       11        63       86         14
~cs_sim\schrodinger_quantum.rs       159       20        21      118         14
symbolic\relativity.rs               155       15        75       65          0
numerical\solve.rs                   155       16        41       98         20
~cs_sim\geodesic_relativity.rs       151       21        39       91          3
lib.rs                               150        7        92       51          1
symbolic\cryptography.rs             150       12        60       78         13
symbolic\multi_valued.rs             149       20        67       62          0
~umerical\geometric_algebra.rs       148        7        27      114          4
numerical\tensor.rs                  148       16        46       86         15
numerical\physics_cfd.rs             146       13        48       85         11
numerical\real_roots.rs              145       13        29      103         43
symbolic\graph_operations.rs         145       15        10      120         41
~cs_sim\navier_stokes_fluid.rs       143       19        20      104         17
symbolic\stats_regression.rs         143       16        53       74          5
symbolic\numeric.rs                  136       10        31       95         32
symbolic\stats.rs                    134        7        59       68          7
~bolic\quantum_field_theory.rs       133       11        68       54          2
~lic\calculus_of_variations.rs       131       14        88       29          0
numerical\graph.rs                   128       15        45       68          7
numerical\topology.rs                127       14        29       84         19
numerical\number_theory.rs           127       11        44       72         33
symbolic\electromagnetism.rs         126       13        56       57          0
numerical\vector.rs                  125        9        68       48         14
output\io.rs                         123       13        44       66          9
~sics_sim\gpe_superfluidity.rs       118       17        21       80          5
numerical\convergence.rs             117        7        41       69          9
~sics_sim\ising_statistical.rs       115       16        23       76         11
numerical\physics_fea.rs             114       10        46       58         10
numerical\physics_md.rs              113       16        41       56          4
numerical\integrate.rs               112        8        20       84         10
~s_sim\fdtd_electrodynamics.rs       109       16        21       72         16
numerical\combinatorics.rs           106        8        40       58         15
numerical\ode.rs                     102       12        20       70          4
~erical\functional_analysis.rs       100        5        50       45          4
symbolic\stats_inference.rs           99        9        22       68          0
numerical\signal.rs                   97       14        26       57          4
output\typst.rs                       96        4         2       90          3
~ical\differential_geometry.rs        95        9        26       60         12
symbolic\special.rs                   91        8        61       22          0
symbolic\mod.rs                       89        3        19       67          0
symbolic\radicals.rs                  87        7        21       59         10
~fractal_geometry_and_chaos.rs        83        7        31       45          5
~mbolic\solid_state_physics.rs        77        6        42       29          0
plugins\mod.rs                        75       13        28       34          2
~fractal_geometry_and_chaos.rs        68        8        35       25          1
numerical\calculus.rs                 68       10        29       29          4
numerical\complex_analysis.rs         64        3        16       45          1
numerical\multi_valued.rs             63        9        19       35          5
symbolic\handles.rs                   57        8        16       33          0
numerical\series.rs                   56        7        17       32          2
numerical\elementary.rs               54        3        16       35          1
~cal\calculus_of_variations.rs        49        5        26       18          0
numerical\mod.rs                      46        1         4       41          0
physics\mod.rs                        46        1        34       11          0
numerical\special.rs                  35        7         6       22          0
prelude.rs                            28        3        19        6          0
numerical\pde.rs                      23        3        15        5          0
output\mod.rs                         22        1        16        5          0
physics\physics_sim\mod.rs             7        0         0        7          0
───────────────────────────────────────────────────────────────────────────────
Total                      140     47649     3824     12531    31294       4526
───────────────────────────────────────────────────────────────────────────────
Estimated Cost to Develop (organic) $1,004,221
Estimated Schedule Effort (organic) 13.78 months
Estimated People Required (organic) 6.48
───────────────────────────────────────────────────────────────────────────────
Processed 1703929 bytes, 1.704 megabytes (SI)
───────────────────────────────────────────────────────────────────────────────
```




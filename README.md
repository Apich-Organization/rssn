# rssn: A Comprehensive Scientific Computing Library for Rust

[![Crates.io](https://img.shields.io/crates/v/rssn.svg)](https://crates.io/crates/rssn)
[![Docs.rs](https://docs.rs/rssn/badge.svg)](https://docs.rs/rssn)
[![License](https://img.shields.io/crates/l/rssn)](LICENSE)

**rssn** is an open-source scientific computing library for Rust, combining **symbolic computation**, **numerical methods**, and **physics simulations** in a single ecosystem.  
It is designed to provide a foundation for building a **next-generation CAS (Computer Algebra System)** and numerical toolkit in Rust.

## Code statistics (SCC 2025-10-07 CST)

───────────────────────────────────────────────────────────────────────────────
Language                 Files     Lines   Blanks  Comments     Code Complexity
───────────────────────────────────────────────────────────────────────────────
Rust                       140     47394     3796     12471    31127       4499
───────────────────────────────────────────────────────────────────────────────
Total                      140     47394     3796     12471    31127       4499
───────────────────────────────────────────────────────────────────────────────
Estimated Cost to Develop (organic) $998,595
Estimated Schedule Effort (organic) 13.75 months
Estimated People Required (organic) 6.45
───────────────────────────────────────────────────────────────────────────────
Processed 1690889 bytes, 1.691 megabytes (SI)
───────────────────────────────────────────────────────────────────────────────

───────────────────────────────────────────────────────────────────────────────
Language                 Files     Lines   Blanks  Comments     Code Complexity
───────────────────────────────────────────────────────────────────────────────
Rust                       140     47394     3796     12471    31127       4499
───────────────────────────────────────────────────────────────────────────────
ffi_api.rs                          3214      196       323     2695        364
symbolic\core.rs                    2219       61       215     1943        123
symbolic\calculus.rs                2149      123       365     1661        328
symbolic\graph_algorithms.rs        1632      124       346     1162        288
symbolic\polynomial.rs              1305       79       514      712        130
symbolic\simplify.rs                1203       59       119     1025        197
symbolic\matrix.rs                  1070       88       309      673        175
symbolic\solve.rs                   1001       79       156      766        134
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
~ymbolic\poly_factorization.rs       565       47       153      365         58
symbolic\transforms.rs               558       28       193      337         39
symbolic\tensor.rs                   551       37       179      335         63
symbolic\topology.rs                 544       42       207      295         55
~ic\error_correction_helper.rs       535       51       142      342         62
symbolic\stats_probability.rs        520       35       172      313          0
numerical\matrix.rs                  482       43       113      326         64
physics\physics_fvm.rs               455       51        83      321         38
symbolic\finite_field.rs             453       40        89      324         57
numerical\stats.rs                   452       58       156      238         12
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
symbolic\cad.rs                      260       27        21      212         35
symbolic\rewriting.rs                259       21        49      189         46
symbolic\graph.rs                    258       21       110      127         15
numerical\error_correction.rs        258       33        60      165         27
physics\physics_bem.rs               250       23        55      172         17
symbolic\error_correction.rs         248       37        62      149         27
symbolic\grobner.rs                  243       29        45      169         19
physics\physics_mm.rs                243       28        22      193         21
numerical\finite_field.rs            239       19        91      129         24
~mbolic\functional_analysis.rs       230       19       110      101          2
~c\stats_information_theory.rs       229       18        90      121         12
physics\physics_fdm.rs               225       24        48      153         26
symbolic\real_roots.rs               222       25        52      145         28
numerical\sparse.rs                  221       23        69      129         16
numerical\polynomial.rs              212       19        46      147         30
output\latex.rs                      208       14        12      182          7
physics\physics_cnm.rs               201       24        41      136         17
~sics_sim\linear_elasticity.rs       198       27        36      135         20
~umerical\computer_graphics.rs       197       16        75      106          7
numerical\coordinates.rs             190       11        68      111         16
symbolic\thermodynamics.rs           189       11       100       78          0
~h_isomorphism_and_coloring.rs       189       17        54      118         29
symbolic\vector_calculus.rs          182       21        62       99          0
symbolic\complex_analysis.rs         172       18        72       82          8
~mbolic\classical_mechanics.rs       171       13        99       59          0
symbolic\optimize.rs                 171       17        54      100         18
numerical\vector_calculus.rs         169       26        52       91          7
numerical\transforms.rs              168       26        57       85         19
plugins\manager.rs                   166       25        16      125         13
symbolic\quantum_mechanics.rs        164       15        78       71          0
symbolic\group_theory.rs             160       11        63       86         14
~cs_sim\schrodinger_quantum.rs       159       20        21      118         14
numerical\solve.rs                   155       16        41       98         20
symbolic\relativity.rs               155       15        75       65          0
~cs_sim\geodesic_relativity.rs       151       21        39       91          3
symbolic\cryptography.rs             150       12        60       78         13
symbolic\multi_valued.rs             149       20        67       62          0
lib.rs                               149        7        91       51          1
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
numerical\number_theory.rs           127       11        44       72         33
numerical\topology.rs                127       14        29       84         19
symbolic\electromagnetism.rs         126       13        56       57          0
numerical\vector.rs                  125        9        68       48         14
output\io.rs                         123       13        44       66          9
numerical\physics_md.rs              119       16        41       62          5
~sics_sim\gpe_superfluidity.rs       118       17        21       80          5
numerical\convergence.rs             117        7        41       69          9
~sics_sim\ising_statistical.rs       115       16        23       76         11
numerical\physics_fea.rs             114       10        46       58         10
numerical\integrate.rs               112        8        20       84         10
~s_sim\fdtd_electrodynamics.rs       109       16        21       72         16
numerical\combinatorics.rs           106        8        40       58         15
~erical\functional_analysis.rs       100        5        50       45          4
symbolic\stats_inference.rs           99        9        22       68          0
numerical\ode.rs                      97       10        20       67          3
numerical\signal.rs                   97       14        26       57          4
output\typst.rs                       96        4         2       90          3
~ical\differential_geometry.rs        95        9        26       60         12
symbolic\special.rs                   91        8        61       22          0
symbolic\mod.rs                       89        3        19       67          0
symbolic\radicals.rs                  87        7        21       59         10
~fractal_geometry_and_chaos.rs        83        7        31       45          5
~mbolic\solid_state_physics.rs        77        6        42       29          0
plugins\mod.rs                        75       13        28       34          2
numerical\calculus.rs                 68       10        29       29          4
~fractal_geometry_and_chaos.rs        68        8        35       25          1
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
Total                      140     47394     3796     12471    31127       4499
───────────────────────────────────────────────────────────────────────────────
Estimated Cost to Develop (organic) $998,595
Estimated Schedule Effort (organic) 13.75 months
Estimated People Required (organic) 6.45
───────────────────────────────────────────────────────────────────────────────
Processed 1690889 bytes, 1.691 megabytes (SI)
───────────────────────────────────────────────────────────────────────────────

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

fn test_differentiate_x_squared_stack_overflow() {
    let x = Expr::Variable("x".to_string());
    let x2 = Expr::Mul(Box::new(x.clone()), Box::new(x.clone()));
    let d = differentiate(&x2, "x");

    // The derivative of x^2 is 2*x.
    // The simplification process might result in Constant(2.0) or BigInt(2).
    let two_const = Expr::Constant(2.0);
    let expected_const = Expr::Mul(Box::new(two_const), Box::new(x.clone()));

    let two_int = Expr::BigInt(BigInt::from(2));
    let expected_int = Expr::Mul(Box::new(two_int), Box::new(x.clone()));

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



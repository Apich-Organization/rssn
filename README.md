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



//! # RSSN: Rust Symbolic and Scientific Numerics
//!
//! `rssn` is a comprehensive library for symbolic mathematics and scientific computing in Rust.
//! It aims to provide a powerful and expressive toolkit for a wide range of mathematical tasks,
//! from symbolic algebra and calculus to advanced numerical simulation.
//!
//! ## Key Features
//!
//! - **Symbolic Computation**: A powerful Computer Algebra System (CAS) for manipulating
//!   mathematical expressions, performing calculus (derivatives, integrals, limits), and solving equations.
//! - **Numerical Methods**: A rich collection of algorithms for numerical integration, solving
//!   differential equations (ODEs and PDEs), optimization, and more.
//! - **Physics Simulation**: High-level tools and examples for simulating physical systems,
//!   including fluid dynamics, electromagnetism, and quantum mechanics.
//! - **Extensibility**: A plugin system (under development) to allow for easy extension of core functionality.
//! - **Versatile Output**: Render expressions as pretty-printed text, LaTeX, or plots.
//!
//! ## Crate Structure
//!
//! The `rssn` crate is organized into the following main modules:
//!
//! - **`symbolic`**: The core of the CAS. It defines the `Expr` tree and provides all
//!   functionality for symbolic manipulation.
//! - **`numerical`**: Contains implementations of various numerical algorithms, such as
//!   quadrature, root-finding, and interpolation.
//! - **`physics`**: Implements numerical methods specifically for physics simulations, such as
//!   the Finite Element Method (FEM), Finite Difference Method (FDM), and various time-stepping schemes.
//! - **`output`**: Provides tools for formatting and displaying expressions in different formats.
//! - **`plugins`**: A placeholder for a future plugin system to extend the library's capabilities.
//! - **`prelude`**: Re-exports the most common types and functions for convenient use.
//!
//! ## Example: Symbolic Differentiation
//!
//! ```rust
//! use rssn::symbolic::calculus::differentiate;
//! use rssn::symbolic::core::Expr;
//! use std::sync::Arc;
//!
//! // Create a symbolic variable 'x'
//! let x = Expr::Variable("x".to_string());
//!
//! // Define an expression: sin(x^2)
//! let expr = Expr::Sin(Arc::new(Expr::Power(Arc::new(x.clone()), Arc::new(Expr::Constant(2.0)))));
//!
//! // Differentiate the expression with respect to 'x'
//! let derivative = differentiate(&expr, "x");
//!
//! // The result will be: (cos(x^2) * (2 * x))
//! // Note: The actual output format may vary.
//! println!("The derivative is: {}", derivative);
//! ```
//!
//! This library is currently in active development. The API may change, and contributions
//! from the community are welcome.
// =========================================================================
// RUST LINT CONFIGURATION: rssn (Scientific Computing Library)
// =========================================================================

// -------------------------------------------------------------------------
// LEVEL 1: CRITICAL ERRORS (Deny)
// -------------------------------------------------------------------------
#![deny(
    // Rust Compiler Errors
    dead_code,
    unreachable_code,
    improper_ctypes_definitions,
    future_incompatible,
    nonstandard_style,
    rust_2018_idioms,
	clippy::perf,
    clippy::correctness,
    clippy::suspicious,
    clippy::unwrap_used,
// =========================================================================
// == LINT DISCUSSIONS AND SUPPRESSIONS (Performance vs. Safety) ===========
// =========================================================================

// clippy::expect_used:
// We allow `expect()` in specific controlled scenarios where an unrecoverable 
// logic error (e.g., failed quantity parsing where input is guaranteed clean)
// is assumed to be a bug, not a runtime failure. For unit testing and quick 
// prototyping, this is often preferred over verbose unwrap_or_else.
//
// ENGINEERING ASSUMPTION: Scientific computing environments often operate 
// under the assumption that user-facing input validation is handled **up-stream** // (at the API or parsing layer). Internal calculations proceed based on the 
// assumption that values (like Quantity components) are valid, non-exceptional 
// numbers, allowing us to favor performance in hot spots over repetitive 
// internal validity checks.

// clippy::indexing_slicing:
// We suppress this globally or per-file because our FVM and grid calculations 
// rely on flattened `Vec<f64>` and index arithmetic (e.g., `idx + width`) 
// for performance.
// 
// DISCUSSION: We all know that in many numerical computing cases, performance 
// and security must be balanced. Direct indexing is structurally safe
// in our validated internal loops, and using it allows the compiler to 
// perform **Bounds Check Elision (BCE)**, which is critical for performance 
// in hot loops. The long-term plan is to refactor these loops to use 
// Ghost Cells or clear range iteration to make the safety provable to Clippy 
// without suppression.

// clippy::arithmetic_side_effects:
// We suppress this because we overload arithmetic operators for custom types 
// (`SupportedQuantity`) where the underlying operations (e.g., `uom`'s types) 
// are **mathematically pure** (side-effect-free). 
// 
// PERFORMANCE CONCERN: The lint encourages explicit function calls over 
// operator overloading for complex types. However, using idiomatic operator 
// overloading (`*`, `/`, `+`, `-`) here keeps the code clean and allows the 
// compiler to better inline and optimize these fundamental algebraic steps 
// within a complex expression tree, which is vital for the performance and 
// readability of scientific calculations. We assert that these operations 
// are side-effect-free.

// clippy::missing_safety_doc:
// This lint is suppressed due to the project's current instability and the 
// high volume of FFI-exported functions (30,000+ lines). 
// 
// REASON FOR SUPPRESSION: All FFI export functions have been correctly marked 
// as `unsafe` (the structural fix). However, manually writing the required 
// `/// # Safety` documentation for every single one of these functions is 
// an insurmountable task for an early-stage project with limited resources. 
//
// TECHNICAL DEBT: This is considered a **critical technical debt item** that 
// must be addressed before the project hits its Beta milestone. A tool 
// or a specialized automated script will be required to audit and enforce 
// detailed safety contracts at that stage. For now, the focus remains on 
// functional correctness.
)]
// -------------------------------------------------------------------------
// LEVEL 2: STYLE WARNINGS (Warn)
// -------------------------------------------------------------------------
#![warn(
    warnings,
    unsafe_code,
    clippy::all,
    clippy::pedantic,
    //clippy::restriction,
    //clippy::nursery,
    clippy::dbg_macro,
    clippy::todo,
    clippy::implicit_clone,
    clippy::unnecessary_safety_comment,
    clippy::same_item_push
)]
// -------------------------------------------------------------------------
// LEVEL 3: ALLOW/IGNORABLE (Allow)
// -------------------------------------------------------------------------
#![allow(
    missing_docs,
    clippy::indexing_slicing,
    clippy::match_same_arms,
    clippy::comparison_chain,
    clippy::redundant_closure_for_method_calls,
    clippy::if_not_else,
    clippy::single_match_else,
    clippy::redundant_else,
    clippy::missing_safety_doc,
    clippy::single_call_fn,
    clippy::min_ident_chars,
    clippy::missing_docs_in_private_items,
    clippy::missing_errors_doc,
    clippy::missing_panics_doc,
    clippy::undocumented_unsafe_blocks,
    clippy::doc_markdown,
    unused_doc_comments,
    clippy::float_arithmetic,
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss,
    clippy::suboptimal_flops,
    clippy::manual_midpoint,
    clippy::non_std_lazy_statics,
    clippy::unreadable_literal,
    clippy::manual_let_else,
    clippy::manual_map,
    clippy::option_if_let_else,
    clippy::empty_line_after_doc_comments,
    clippy::many_single_char_names,
    clippy::module_name_repetitions,
    clippy::redundant_field_names,
    clippy::similar_names,
    clippy::redundant_pub_crate,
    clippy::too_many_lines,
    clippy::must_use_candidate,
    clippy::shadow_unrelated,
    clippy::use_self,
    clippy::str_to_string,
    clippy::uninlined_format_args
)]
#[cfg(feature = "ffi_api")]
pub mod ffi_apis;
#[cfg(feature = "ffi_blinding")]
pub mod ffi_blindings;
pub mod numerical;
#[cfg(feature = "output")]
pub mod output;
#[cfg(feature = "physics")]
pub mod physics;
#[cfg(feature = "plugins")]
pub mod plugins;
#[cfg(feature = "full")]
pub mod prelude;
pub mod symbolic;
use std::sync::Arc;
/// Checks if an `Arc` has exclusive ownership (strong count is 1).
///
/// This is useful for optimizations where you want to mutate the contained data
/// in-place, avoiding a clone. If this returns `true`, `Arc::get_mut` or
/// `Arc::try_unwrap` will succeed.
#[allow(clippy::inline_always)]
#[inline(always)]
pub fn is_exclusive<T>(arc: &Arc<T>) -> bool {
    Arc::strong_count(arc) == 1
}

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
#![deny(
    dead_code,
    unreachable_code,
    improper_ctypes_definitions,
    future_incompatible,
    nonstandard_style,
    rust_2018_idioms,
    clippy::perf,
    clippy::correctness,
    clippy::suspicious,
    clippy::unwrap_used
)]
#![warn(
    warnings,
    unsafe_code,
    clippy::all,
    clippy::pedantic,
    clippy::dbg_macro,
    clippy::todo,
    clippy::implicit_clone,
    clippy::unnecessary_safety_comment,
    clippy::same_item_push
)]
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

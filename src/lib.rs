//! # RSSN: Rust Symbolic and Scientific Numerics
//!
//! `rssn` is a high-performance, modern library for symbolic mathematics and scientific
//! computing in Rust. It leverages a Directed Acyclic Graph (DAG) based model for
//! efficient, canonical representation of mathematical expressions.
//!
//! ## Key Features
//!
//! - **Efficient DAG-based Expression Model**: Expressions are stored as a Directed Acyclic
//!   Graph, ensuring that identical subexpressions are represented by a single node in memory.
//!   This provides automatic canonicalization and significant performance gains.
//!
//! - **Advanced Symbolic Algebra**: A powerful Computer Algebra System (CAS) for manipulating
//!   expressions. It goes beyond simple simplification by supporting:
//!     - Polynomial algebra including **Gröbner basis** computation.
//!     - Simplification and normalization of expressions with respect to polynomial side-relations.
//!
//! - **Symbolic Calculus**: Functions for differentiation, integration, limits, and series expansion.
//!
//! - **Numerical Methods**: A rich collection of algorithms for numerical integration, solving
//!   differential equations (ODEs and PDEs), optimization, and more.
//!
//! - **Versatile Output**: Render expressions as pretty-printed text, LaTeX, or plots.
//!
//! ## Crate Structure
//!
//! The `rssn` crate is organized into the following main modules:
//!
//! - **`symbolic`**: The core of the CAS. It defines the `Expr` and `DagNode` representations and
//!   provides all functionality for symbolic manipulation, including the `simplify_dag` engine
//!   and advanced algebraic tools in `cas_foundations` and `grobner`.
//! - **`numerical`**: Contains implementations of various numerical algorithms.
//! - **`physics`**: Implements numerical methods specifically for physics simulations.
//! - **`output`**: Provides tools for formatting and displaying expressions.
//! - **`prelude`**: Re-exports the most common types and functions for convenient use.
//!
//! ## Example: Simplification with Relations
//!
//! `rssn` can simplify expressions within the context of an algebraic variety. For example,
//! simplifying `x^2` given the side-relation that `x^2 + y^2 - 1 = 0` (the unit circle).
//!
//! ```rust
//! 
//! use rssn::symbolic::cas_foundations::simplify_with_relations;
//! use rssn::symbolic::core::Expr;
//! use rssn::symbolic::grobner::MonomialOrder;
//!
//! // Define variables x and y
//! let x = Expr::new_variable("x");
//!
//! let y = Expr::new_variable("y");
//!
//! // Expression to simplify: 2*x^2
//! let two = Expr::new_bigint(2.into());
//!
//! let x_sq = Expr::new_pow(
//!     x.clone(),
//!     Expr::new_bigint(2.into()),
//! );
//!
//! let expr_to_simplify = Expr::new_mul(two, x_sq);
//!
//! // Define the side-relation: x^2 + y^2 - 1 = 0
//! let y_sq = Expr::new_pow(
//!     y.clone(),
//!     Expr::new_bigint(2.into()),
//! );
//!
//! let one = Expr::new_bigint(1.into());
//!
//! let relation = Expr::new_sub(
//!     Expr::new_add(x.clone(), y.clone()),
//!     one,
//! );
//!
//! // Simplify the expression with respect to the relation
//! let simplified_expr = simplify_with_relations(
//!     &expr_to_simplify,
//!     &[relation],
//!     &["x", "y"],
//!     MonomialOrder::Lexicographical,
//! );
//!
//! // The result will be 2 - 2y^2
//! // Note: The exact output format and canonical form may vary.
//! println!(
//!     "Original expression: {}",
//!     expr_to_simplify
//! );
//!
//! println!(
//!     "Simplified expression: {}",
//!     simplified_expr
//! );
//! ```
//!
//! This library is in active development. The API may change, and community
//! contributions are welcome.

// =========================================================================
// RUST LINT CONFIGURATION: rssn (Scientific Computing Library) -- version 1
// =========================================================================
//
// -------------------------------------------------------------------------
// LEVEL 1: CRITICAL ERRORS (Deny)
// -------------------------------------------------------------------------
// #![deny(
// Rust Compiler Errors
// dead_code,
// unreachable_code,
// improper_ctypes_definitions,
// future_incompatible,
// nonstandard_style,
// rust_2018_idioms,
// clippy::perf,
// clippy::correctness,
// clippy::suspicious,
// clippy::unwrap_used,
// =========================================================================
// == LINT DISCUSSIONS AND SUPPRESSIONS (Performance vs. Safety) ===========
// =========================================================================
//
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
//
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
//
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
//
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
// )]
// -------------------------------------------------------------------------
// LEVEL 2: STYLE WARNINGS (Warn)
// -------------------------------------------------------------------------
// #![warn(
// warnings,
// unsafe_code,
// clippy::all,
// clippy::pedantic,
// clippy::nursery,
// clippy::dbg_macro,
// clippy::todo,
// clippy::implicit_clone,
// clippy::unnecessary_safety_comment,
// clippy::same_item_push
// )]
// -------------------------------------------------------------------------
// LEVEL 3: ALLOW/IGNORABLE (Allow)
// -------------------------------------------------------------------------
// #![allow(
// missing_docs,  // Temporary: To be addressed gradually
// clippy::indexing_slicing,  // Performance-critical in numerical computations
// clippy::match_same_arms,  // Allowed for code clarity in some cases
// clippy::comparison_chain,  // Sometimes more readable than alternatives
// clippy::redundant_closure_for_method_calls,  // For clarity in some cases
// clippy::if_not_else,  // Sometimes the negative condition is clearer
// clippy::single_match_else,  // Sometimes clearer than refactoring
// clippy::redundant_else,  // Sometimes clearer to explicitly show all branches
// clippy::missing_safety_doc,  // FFI functions that need safety docs to be completed
// clippy::single_call_fn,  // Sometimes used for clarity or future expansion
// clippy::min_ident_chars,  // Mathematical notation often uses short variable names
// clippy::missing_docs_in_private_items,  // Private items temporarily missing docs
// clippy::missing_errors_doc,  // Documentation to be completed
// clippy::missing_panics_doc,  // Documentation to be completed
// clippy::undocumented_unsafe_blocks,  // Safety comments to be added
// clippy::doc_markdown,  // Technical terms sometimes need backticks
// unused_doc_comments,  // To be addressed in documentation improvements
// clippy::float_arithmetic,  // Required in numerical computing
// clippy::cast_possible_truncation,  // Sometimes necessary in numerical computing
// clippy::cast_precision_loss,  // Sometimes necessary in numerical computing
// clippy::cast_sign_loss,  // Sometimes necessary in numerical computing
// clippy::suboptimal_flops,  // Performance trade-offs in numerical methods
// clippy::manual_midpoint,  // Manual implementation for numerical stability
// clippy::non_std_lazy_statics,  // For performance in some cases
// clippy::unreadable_literal,  // Sometimes mathematical constants need specific formatting
// clippy::manual_let_else,  // Not yet stabilized feature
// clippy::manual_map,  // Sometimes explicit control flow is clearer
// clippy::option_if_let_else,  // Sometimes match is clearer
// clippy::empty_line_after_doc_comments,  // Formatting preference for some cases
// clippy::many_single_char_names,  // Mathematical notation often uses single letters
// clippy::module_name_repetitions,  // Sometimes module names naturally repeat
// clippy::redundant_field_names,  // Sometimes clearer to be explicit
// clippy::similar_names,  // Sometimes mathematical variables are naturally similar
// clippy::redundant_pub_crate,  // For API consistency
// clippy::too_many_lines,  // Complex functions that need refactoring over time
// clippy::must_use_candidate,  // To be addressed gradually
// clippy::shadow_unrelated,  // Sometimes appropriate for local variables
// clippy::use_self,  // For consistency in some cases
// clippy::str_to_string,  // Sometimes needed for API compatibility
// clippy::uninlined_format_args,  // Performance considerations in hot paths
// clippy::collapsible_if,
// clippy::single_match,
// clippy::needless_pass_by_value,
// clippy::needless_pass_by_ref_mut,
// clippy::used_underscore_binding,
// )]
//
// =========================================================================
// RUST LINT CONFIGURATION: rssn (Scientific Computing Library) -- version 2
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
    clippy::missing_safety_doc,
    clippy::same_item_push,
    clippy::implicit_clone,
    clippy::all,
    clippy::pedantic,
    clippy::nursery,
    clippy::single_call_fn,
)]
// -------------------------------------------------------------------------
// LEVEL 2: STYLE WARNINGS (Warn)
// -------------------------------------------------------------------------
#![warn(
    warnings,
    missing_docs,
    unsafe_code,
    // To avoid performance issues in hot paths
    clippy::expect_used,
    // To avoid simd optimization issues
    clippy::indexing_slicing,
    // To avoid simd optimization issues
    clippy::arithmetic_side_effects,
    // Precision loss is allowed here because we need to introduce bigint to all places otherwise --- that will require a lot of work and results in breaking change and loss of performance. So we will have to handle it later.
    // DEBT: Handle precision loss when possible (add more suites of code)
    // Possible Truncation warnned, due to CPU branch prediction and simd optimization programs, we will just warn this problems instead deny it.
    clippy::cast_precision_loss,
    clippy::cast_possible_wrap,
    clippy::cast_possible_truncation,
    clippy::dbg_macro,
    clippy::todo,
    clippy::unnecessary_safety_comment
)]
// -------------------------------------------------------------------------
// LEVEL 3: ALLOW/IGNORABLE (Allow)
// -------------------------------------------------------------------------
#![allow(
    clippy::restriction,
    clippy::inline_always,
    unused_doc_comments,
    clippy::many_single_char_names,
    clippy::similar_names,
    clippy::redundant_else,
    clippy::needless_continue,
    clippy::empty_line_after_doc_comments,
    clippy::empty_line_after_outer_attr,
    clippy::manual_let_else,
    // It is always reporting on normal math writings.
    clippy::doc_markdown,
    // We thinks do not collapsible if makes the code more extensible.
    clippy::collapsible_if,
    clippy::collapsible_match,
    clippy::collapsible_else_if
)]

/// Computation engine and task management.
pub mod compute;
/// System and physical constants.
pub mod constant;
#[cfg(feature = "ffi_api")]
pub mod ffi_apis;

// #[instability::unstable(feature = "experimental")]
// Disabled because it only works on nightly rust
#[cfg(feature = "ffi_blinding")]
/// FFI blinding and security utilities.
pub mod ffi_blindings;
/// Input parsing and handling.
pub mod input;
#[cfg(feature = "jit")]
/// Just-In-Time (JIT) compilation for expressions.
pub mod jit;
#[cfg(feature = "nightly")]
/// Features requiring nightly Rust.
pub mod nightly;
pub mod numerical;
#[cfg(feature = "output")]
pub mod output;
#[cfg(feature = "physics")]
pub mod physics;
#[cfg(feature = "plugins")]
/// Plugin system for extending functionality.
pub mod plugins;
pub mod prelude;
pub mod symbolic;

use std::sync::Arc;

/// Checks if an `Arc` has exclusive ownership (strong count is 1).
///
/// This is useful for optimizations where you want to mutate the contained data
/// in-place, avoiding a clone. If this returns `true`, `Arc::get_mut` or
/// `Arc::try_unwrap` will succeed.
///
/// # Arguments
/// * `arc` - A reference to the `Arc` to check for exclusive ownership
///
/// # Returns
/// * `bool` - True if the Arc has exclusive ownership, false otherwise
#[allow(clippy::inline_always)]
#[inline(always)]

pub fn is_exclusive<T>(
    arc: &Arc<T>
) -> bool {

    Arc::strong_count(arc) == 1
}

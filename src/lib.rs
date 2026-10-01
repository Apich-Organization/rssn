//! # rssn
//!
//! A scientific computing engine built on one idea: **almost every
//! numerical or symbolic operation is an identity transformation** — it
//! replaces an expression by an equal one. Differentiating, integrating,
//! solving, simplifying, factoring a matrix and evaluating to a float differ
//! only in which equal form is wanted.
//!
//! So instead of a catalogue of separate methods, rssn has two parts:
//!
//! * **Graph computation** — [`graph`] is a domain-free kernel: a
//!   hash-consed expression DAG whose nodes are linked into equivalence
//!   classes, tree windows for destructive local search, and a heuristic
//!   scheduler. [`rules`] supplies the mathematics as rule sets: declarative
//!   rewrites and procedural reduction kernels, symbolic and numeric side by
//!   side. [`api`] is the single entry point, [`Session::compute`].
//! * **Simulation** — [`sim`] holds what is not an identity transformation:
//!   time stepping, discretised fields, particle systems. [`kernels`] are
//!   the plain numeric routines both parts share, and [`backend`] turns a
//!   reduced term into an executable function for them.
//!
//! ```
//! use rssn::prelude::*;
//!
//! let s = Session::new();
//! let x = s.sym("x");
//!
//! // One request, two phases.
//! let request = diff(x.pow(3) * sin(x), x);
//! let closed_form = s.compute(request, &Config::new())?.term;
//! let at_two = s.compute(request, &Config::new().numeric(1e-10).bind("x", 2.0))?.value;
//!
//! assert_eq!(closed_form.eval(&[("x", 2.0)]), at_two);
//! # Ok::<(), ComputeError>(())
//! ```
//!
//! See `ARCHITECTURE.md` for the design.

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
    missing_docs,
    unsafe_code,
)]
// -------------------------------------------------------------------------
// LEVEL 2: STYLE WARNINGS (Warn)
// -------------------------------------------------------------------------
#![warn(
    warnings,
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
    // This is usually a sign of dead code --- but for development purposes, we will just warn it.
    clippy::used_underscore_binding,
    clippy::unnecessary_safety_comment
)]
// -------------------------------------------------------------------------
// LEVEL 3: ALLOW/IGNORABLE (Allow)
// -------------------------------------------------------------------------
#![allow(
    clippy::restriction,
    // `mul_add` is a slow libm call on targets without hardware FMA, and
    // fused rounding changes results across targets.
    clippy::suboptimal_flops,
    // Fires on `&[a, b]` built from separately computed bindings.
    clippy::tuple_array_conversions,
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
    clippy::collapsible_else_if,
    clippy::no_effect_underscore_binding,
    clippy::must_use_candidate
)]
#![doc(
    html_logo_url = "https://raw.githubusercontent.com/Apich-Organization/rssn/refs/heads/dev/doc/logo.png"
)]
#![doc(
    html_favicon_url = "https://raw.githubusercontent.com/Apich-Organization/rssn/refs/heads/dev/doc/favicon.ico"
)]

/// Unified mathematical compute engine.
pub mod api;
pub mod backend;
/// Build metadata and mathematical constants.
pub mod constant;
pub mod ffi;
pub mod graph;
pub mod io;
pub mod kernels;
pub mod prelude;
pub mod rules;
pub mod sim;

pub use api::Session;

/// The README's code blocks, compiled and run as doctests.
#[cfg(doctest)]
#[doc = include_str!("../README.md")]
pub struct ReadmeDoctests;

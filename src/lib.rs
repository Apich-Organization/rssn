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

// -------------------------------------------------------------------------
// LEVEL 1: CRITICAL ERRORS (Deny)
// -------------------------------------------------------------------------
#![deny(
    unreachable_code,
    improper_ctypes_definitions,
    future_incompatible,
    nonstandard_style,
    rust_2018_idioms,
    clippy::perf,
    clippy::correctness,
    clippy::suspicious,
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    clippy::missing_safety_doc,
    clippy::same_item_push,
    clippy::implicit_clone,
    clippy::all,
    clippy::pedantic,
    missing_docs,
    clippy::nursery,
    clippy::single_call_fn,
)]
// -------------------------------------------------------------------------
// LEVEL 2: STYLE WARNINGS (Warn)
// -------------------------------------------------------------------------
#![warn(
    dead_code,
    warnings,
    unsafe_code,
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
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap,
    clippy::empty_line_after_doc_comments
)]
// -------------------------------------------------------------------------
// Scientific-computing compromises (documented downgrades)
// -------------------------------------------------------------------------
#![allow(
    // `mul_add` is a slow libm call on targets without hardware FMA, and a
    // fused rounding changes results between targets; the kernels keep the
    // written evaluation order so results are reproducible everywhere.
    clippy::suboptimal_flops,
    // Counts and indices become `f64` throughout numerics (sample sizes,
    // node positions, orders); they stay far below 2^53.
    clippy::cast_precision_loss,
    // Formulas follow the notation of the literature (`a, b, c, x, y, z`,
    // `x0, x1`, `p, q, r, s`), which reads better than spelled-out names.
    clippy::many_single_char_names,
    clippy::similar_names,
    // Fires on `[a, b]` built from separately computed bindings.
    clippy::tuple_array_conversions
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

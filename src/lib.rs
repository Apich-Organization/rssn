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
//! use rssn::prelude::*;
//!
//! // Create a symbolic variable 'x'
//! let x = Expr::Variable("x".to_string());
//!
//! // Define an expression: sin(x^2)
//! let expr = Expr::Sin(Box::new(Expr::Power(Box::new(x.clone()), Box::new(Expr::Constant(2.0)))));
//!
//! // Differentiate the expression with respect to 'x'
//! let derivative = diff(&expr, "x");
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
    //warnings,                
    //unsafe_code,
    future_incompatible,
    nonstandard_style,
    rust_2018_idioms,
    dead_code,
    unreachable_code,

/*
Disabled during v0.1.x releases.
    clippy::all,              
    clippy::cargo,            
	
    clippy::pedantic,         

    clippy::unwrap_used,      
    clippy::expect_used,      
    clippy::indexing_slicing, 
    clippy::get_unwrap,       
    clippy::integer_arithmetic, 
    clippy::float_arithmetic, 
    clippy::default_trait_access, 
    clippy::unnecessary_safety_comment, 
    clippy::redundant_closure_for_method_call,
    clippy::same_item_push,   
*/
)]
// -------------------------------------------------------------------------
// LEVEL 2: STYLE WARNINGS (Warn)
// -------------------------------------------------------------------------
#![warn(
    //missing_docs,
	//Temporarily disabled due to project progress reasons.
    clippy::missing_errors_doc,
    clippy::missing_panics_doc,
    clippy::dbg_macro,
    clippy::todo,
    clippy::implicit_clone,
    clippy::str_to_string,
    clippy::string_add,
    clippy::undocumented_unsafe_blocks
)]
// -------------------------------------------------------------------------
// LEVEL 3: ALLOW/IGNORABLE (Allow)
// -------------------------------------------------------------------------
#![allow(
    non_snake_case,
    unused_variables,
    unused_imports,
    unused_doc_comments,
    clippy::module_name_repetitions,
    clippy::too_many_lines,
    clippy::must_use_candidate,
    clippy::similar_names,
    clippy::redundant_pub_crate,
    clippy::shadow_unrelated,
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss,
    clippy::use_self,
    clippy::redundant_field_names
)]

pub mod numerical;
#[cfg(feature = "full")]
pub mod output;
#[cfg(feature = "full")]
pub mod physics;
#[cfg(feature = "full")]
pub mod plugins;
pub mod prelude;
pub mod symbolic;

#[cfg(feature = "full")]
pub mod ffi_api;



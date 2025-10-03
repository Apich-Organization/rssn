//! # Crate Prelude
//!
//! This module re-exports the most commonly used types, traits, and functions from across the `rssn`
//! crate. The purpose of a prelude is to make it easier to use the library's key features without
//! having to import them one by one.
//!
//! By including `use rssn::prelude::*;` at the top of a file, you can gain convenient access to
//! the items exported here.
//!
//! ## Contents
//!
//! - **`Expr`**: The core symbolic expression enum and all of its variants.
//! - **Symbolic Operations**: Key functions like `diff` (for differentiation), `integrate`, `limit`,
//!   `simplify`, and `solve`.
//! - **Numerical Types**: Common data structures used in numerical computation, such as `CsMat` for
//!   sparse matrices.

// Re-export common types and functions for easier use across the crate.

pub use crate::symbolic::core::Expr;
pub use crate::symbolic::core::Expr::*;
// Common symbolic operations
pub use crate::symbolic::calculus::{differentiate as diff, integrate, limit};
pub use crate::symbolic::simplify::simplify;
pub use crate::symbolic::solve::solve;

// Common numerical types
pub use sprs::CsMat;

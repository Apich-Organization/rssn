//! FFI bindings for symbolic logic operations.
//!
//! This module exposes foreign function interface (FFI) bindings for propositional
//! and first-order logic, including formula construction, simplification, and
//! satisfiability checking.

pub mod bincode_api;
/// Handle-based FFI bindings for logic formulas using opaque `Expr` handles.
pub mod handle;
/// JSON-based FFI bindings for logic formulas using serialized expressions.
pub mod json;

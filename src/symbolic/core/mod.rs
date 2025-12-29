//! # Symbolic Expression Core Module
//!
//! This module provides the foundational data structures and operations for symbolic
//! mathematics in the RSSN library. It implements a hybrid AST/DAG (Abstract Syntax Tree /
//! Directed Acyclic Graph) representation system for mathematical expressions.
//!
//! ## Architecture Overview
//!
//! The core of this module is the [`Expr`] enum, which represents symbolic mathematical
//! expressions. The system is designed around a dual representation strategy:
//!
//! 1. **AST Representation**: Traditional tree-based expression structure (legacy)
//! 2. **DAG Representation**: Graph-based structure with shared subexpressions (modern)
//!
//! The DAG representation is managed by the [`DagManager`], which provides:
//! - Automatic deduplication of identical subexpressions
//! - Hash-based node lookup for O(1) retrieval
//! - Canonical normalization of expressions
//! - Memory-efficient storage through structural sharing
//!
//! ## Key Components
//!
//! ### Expression Types ([`Expr`])
//!
//! The [`Expr`] enum supports a comprehensive set of mathematical operations:
//!
//! - **Atomic Values**: Constants, variables, patterns, special constants (π, e, ∞)
//! - **Arithmetic**: Addition, subtraction, multiplication, division, power, negation
//! - **Trigonometric**: sin, cos, tan, sec, csc, cot and their inverses
//! - **Hyperbolic**: sinh, cosh, tanh, sech, csch, coth and their inverses
//! - **Special Functions**: Gamma, Beta, Bessel, Legendre, Laguerre, Hermite, etc.
//! - **Calculus**: Derivatives, integrals, limits, series, summations
//! - **Linear Algebra**: Matrices, vectors, transpose, inverse, matrix multiplication
//! - **Logic**: Boolean operations, predicates, quantifiers
//! - **Advanced**: ODEs, PDEs, distributions, complex analysis
//!
//! ### N-ary Operations
//!
//! The module includes efficient n-ary operation variants:
//!
//! - [`Expr::AddList`]: Sum of multiple terms in a single operation
//! - [`Expr::MulList`]: Product of multiple factors in a single operation
//!
//! These variants improve performance by reducing tree depth and enabling better
//! optimization opportunities during simplification.
//!
//! ### Dynamic Operations
//!
//! The module supports runtime-extensible operations through:
//!
//! - [`Expr::UnaryList`]: Custom unary operations
//! - [`Expr::BinaryList`]: Custom binary operations  
//! - [`Expr::NaryList`]: Custom n-ary operations
//!
//! These are registered in the [`DYNAMIC_OP_REGISTRY`] with properties like
//! associativity and commutativity, enabling plugin systems and domain-specific
//! extensions without modifying the core enum.
//!
//! ### DAG Management
//!
//! The [`DagManager`] provides centralized management of DAG nodes:
//!
//! ```rust
//! 
//! use rssn::symbolic::core::{Expr, DAG_MANAGER};
//!
//! // Create expressions using smart constructors
//! let x = Expr::new_variable("x");
//!
//! let two = Expr::new_constant(2.0);
//!
//! let expr = Expr::new_add(x, two);
//!
//! // The DAG_MANAGER automatically deduplicates identical subexpressions
//! ```
//!
//! ### Smart Constructors
//!
//! All operations have corresponding smart constructors (e.g., `new_add`, `new_mul`)
//! that automatically:
//! - Convert to DAG representation
//! - Normalize the expression
//! - Deduplicate subexpressions
//! - Apply basic simplifications
//!
//! ## AST to DAG Migration
//!
//! The module is undergoing a gradual migration from AST to DAG representation:
//!
//! - **Legacy AST forms** (e.g., `Expr::Add(Arc<Expr>, Arc<Expr>)`) remain for compatibility
//! - **Modern DAG forms** use `Expr::Dag(Arc<DagNode>)` wrapper
//! - **Smart constructors** automatically create DAG forms
//! - **Conversion utilities** (`to_dag()`, `to_ast()`) enable interoperability
//!
//! This hybrid approach ensures backward compatibility while enabling new optimizations.
//!
//! ## Expression Traversal
//!
//! The module provides multiple traversal methods:
//!
//! - [`Expr::pre_order_walk`]: Visit parent before children
//! - [`Expr::post_order_walk`]: Visit children before parent
//! - [`Expr::in_order_walk`]: Visit left child, parent, then right child
//!
//! ## Normalization and Canonicalization
//!
//! The [`Expr::normalize`] method provides canonical forms:
//! - Sorts commutative operation children
//! - Flattens nested associative operations
//! - Applies consistent ordering for hashing
//!
//! ## Examples
//!
//! ### Basic Expression Creation
//!
//! ```rust
//! 
//! use rssn::symbolic::core::Expr;
//!
//! // Using smart constructors (recommended)
//! let x = Expr::new_variable("x");
//!
//! let y = Expr::new_variable("y");
//!
//! let sum = Expr::new_add(x.clone(), y.clone());
//!
//! let product = Expr::new_mul(x, y);
//! ```
//!
//! ### N-ary Operations
//!
//! ```rust
//! 
//! use rssn::symbolic::core::Expr;
//!
//! // Efficient multi-term addition
//! let sum = Expr::AddList(vec![
//!     Expr::Variable("a".to_string()),
//!     Expr::Variable("b".to_string()),
//!     Expr::Variable("c".to_string()),
//!     Expr::Variable("d".to_string()),
//! ]);
//! ```
//!
//! ### Dynamic Operations
//!
//! ```rust
//! 
//! use std::sync::Arc;
//!
//! use rssn::symbolic::core::register_dynamic_op;
//! use rssn::symbolic::core::DynamicOpProperties;
//! use rssn::symbolic::core::Expr;
//!
//! // Register a custom operation
//! register_dynamic_op(
//!     "custom_func",
//!     DynamicOpProperties {
//!         name : "custom_func".to_string(),
//!         description : "My custom function".to_string(),
//!         is_associative : false,
//!         is_commutative : false,
//!     },
//! );
//!
//! // Use it
//! let expr = Expr::UnaryList(
//!     "custom_func".to_string(),
//!     Arc::new(Expr::Variable(
//!         "x".to_string(),
//!     )),
//! );
//! ```
//!
//! ## Performance Considerations
//!
//! - **DAG representation** reduces memory usage through structural sharing
//! - **Hash-based deduplication** provides O(1) lookup for common subexpressions
//! - **N-ary operations** reduce tree depth and improve cache locality
//! - **Lazy evaluation** defers expensive operations until needed
//!
//! ## Thread Safety
//!
//! - The [`DAG_MANAGER`] uses internal locking for thread-safe access
//! - The [`DYNAMIC_OP_REGISTRY`] uses `RwLock` for concurrent reads
//! - Individual [`Expr`] values are immutable and can be shared across threads
//!
//! ## See Also
//!
//! - [`simplify_dag`](crate::symbolic::simplify_dag) - Modern DAG-based simplification
//! - [`simplify`](crate::symbolic::simplify) - Legacy AST-based simplification (deprecated)
//! - [`calculus`](crate::symbolic::calculus) - Symbolic differentiation and integration
//! - [`elementary`](crate::symbolic::elementary) - Elementary function transformations

// Unavoidable for intermodule functionality issues.
#![allow(deprecated)]
// Unavoidable for intermodule functionality development issues.
#![allow(unused_imports)]

pub use api::*;

pub use dag_mgr::*;
pub use expr::*;
pub use expr_impl::*;


/// Public API and constructors for symbolic expressions.
pub mod api;
/// Abstract Syntax Tree (AST) specific implementations.
pub mod ast_impl;
/// DAG management and deduplication logic.
pub mod dag_mgr;
/// Core expression type definition.
pub mod expr;
/// Core expression implementation details.
pub mod expr_impl;
/// Traits and utilities for converting values to Expressions.
pub mod to_expr;

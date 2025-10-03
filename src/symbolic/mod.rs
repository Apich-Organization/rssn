//! # Symbolic Computation Module
//!
//! This module is the heart of the `rssn` computer algebra system (CAS). It provides the data structures
//! and algorithms for symbolic manipulation of mathematical expressions.
//!
//! ## Core Components
//!
//! - `core`: Defines the fundamental `Expr` enum, which represents the Abstract Syntax Tree (AST) of
//!   mathematical expressions.
//! - `calculus`: Contains functions for symbolic differentiation, integration, limits, and series expansion.
//! - `algebra`: Provides tools for expression simplification, expansion, and solving equations.
//! - `linear_algebra`: Implements symbolic operations on matrices and vectors.
//! - `logic`: Handles Boolean algebra, predicates, and logical quantifiers.
//!
//! ## Sub-modules
//!
//! The module is organized into various sub-modules, each dedicated to a specific area of mathematics,
//! such as `complex_analysis`, `number_theory`, `geometry`, and more. This modular structure allows
//! for clear separation of concerns and facilitates future expansion.

pub mod cad;
pub mod calculus;
pub mod calculus_of_variations;
pub mod cas_foundations;
pub mod classical_mechanics;
pub mod combinatorics;
pub mod complex_analysis;
pub mod computer_graphics;
pub mod convergence;
pub mod coordinates;
pub mod core;
pub mod cryptography;
pub mod differential_geometry;
pub mod discrete_groups;
pub mod electromagnetism;
pub mod elementary;
pub mod error_correction;
pub mod error_correction_helper;

pub mod finite_field;
pub mod fractal_geometry_and_chaos;
pub mod functional_analysis;
pub mod geometric_algebra;
pub mod graph;
pub mod graph_algorithms;
pub mod graph_isomorphism_and_coloring;
pub mod graph_operations;
pub mod grobner;
pub mod group_theory;
pub mod integral_equations;
pub mod integration;
pub mod lie_groups_and_algebras;

pub mod logic;
pub mod matrix;
pub mod multi_valued;
pub mod number_theory;
pub mod numeric;
pub mod ode;
pub mod optimize;
pub mod pde;
pub mod poly_factorization;
pub mod polynomial;
pub mod proof;
pub mod quantum_field_theory;
pub mod quantum_mechanics;
pub mod radicals;
pub mod real_roots;
pub mod relativity;
pub mod rewriting;
pub mod series;
pub mod simplify;
pub mod solid_state_physics;
pub mod solve;
pub mod special;
pub mod special_functions;
pub mod stats;
pub mod stats_inference;
pub mod stats_information_theory;
pub mod stats_probability;
pub mod stats_regression;
pub mod tensor;
pub mod thermodynamics;
pub mod topology;
pub mod transforms;
pub mod unit_unification;
pub mod vector;
pub mod vector_calculus;

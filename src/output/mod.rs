//! # Output Formatting Module
//!
//! This module provides functionalities for rendering and presenting mathematical expressions
//! in various formats. It allows the internal `Expr` representation to be converted into
//! human-readable or machine-processable outputs.
//!
//! ## Sub-modules
//!
//! - `pretty_print`: Generates formatted, indented, and colored text output suitable for
//!   display in a terminal.
//! - `latex`: Converts expressions into LaTeX strings, allowing for high-quality typesetting
//!   of mathematical formulas.
//! - `plotting`: Provides tools to generate plots and visualizations of functions and data.
//!   (Note: This might interface with external plotting libraries).
//! - `io`: Handles saving and loading of expressions and results to/from files.
//! - `typst`: Converts expressions into Typst code, a modern typesetting system.

pub mod io;
pub mod latex;
pub mod plotting;
pub mod pretty_print;
pub mod typst;

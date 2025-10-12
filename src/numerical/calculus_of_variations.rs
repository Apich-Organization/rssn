//! # Numerical Calculus of Variations
//!
//! This module provides numerical tools for the calculus of variations.
//! It focuses on evaluating the action of a functional for a given path,
//! which is a fundamental step in solving problems like finding geodesics
//! or optimizing paths in physics and engineering.

use crate::numerical::integrate::{quadrature, QuadratureMethod};
use crate::symbolic::calculus::differentiate;
use crate::symbolic::calculus::substitute;
use crate::symbolic::core::Expr;

/// Evaluates the action of a functional for a given path.
///
/// The functional is `S[y] = integral from a to b of L(t, y, y_dot) dt`.
/// This function numerically computes this integral by first substituting the `path`
/// and its derivative into the `lagrangian`, and then performing numerical quadrature.
///
/// # Arguments
/// * `lagrangian` - The Lagrangian expression `L`. It should be an expression in terms of
///   `t_var`, `path_var`, and `path_dot_var`.
/// * `path` - The actual path `y(t)` as an expression.
/// * `t_var` - The name of the independent variable (e.g., "t").
/// * `path_var` - The name of the path variable used in the Lagrangian (e.g., "y").
/// * `path_dot_var` - The name of the path's derivative used in the Lagrangian (e.g., "y_dot").
/// * `t_range` - The interval of integration `(a, b)`.
///
/// # Returns
/// A `Result` containing the numerical value of the action integral.
pub fn evaluate_action(
    lagrangian: &Expr,
    path: &Expr,
    t_var: &str,
    path_var: &str,
    path_dot_var: &str,
    t_range: (f64, f64),
) -> Result<f64, String> {
    // 1. Find the derivative of the given path.
    let path_dot = differentiate(path, t_var);

    // 2. Substitute the path and its derivative into the Lagrangian.
    let integrand_with_y = substitute(lagrangian, path_var, path);
    let integrand = substitute(&integrand_with_y, path_dot_var, &path_dot);

    // 3. Numerically integrate the resulting expression over the given range.
    quadrature(&integrand, t_var, t_range, 1000, QuadratureMethod::Simpson)
}

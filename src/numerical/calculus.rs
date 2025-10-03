//! # Numerical Calculus
//!
//! This module provides numerical calculus operations, primarily focusing on
//! finite difference methods for approximating derivatives. It includes functions
//! for computing the numerical gradient of multivariate functions.

use crate::numerical::elementary::eval_expr;
use crate::symbolic::core::Expr;
use std::collections::HashMap;

/// Computes the numerical gradient of a multivariate function at a given point.
///
/// The gradient is approximated using finite differences. For each variable, the partial
/// derivative is estimated by evaluating the function at `x + h` and `x - h`.
///
/// # Arguments
/// * `f` - The expression representing the function.
/// * `vars` - The variables of the function.
/// * `point` - The point at which to compute the gradient.
///
/// # Returns
/// A `Result` containing the gradient vector, or an error string.
pub fn gradient(f: &Expr, vars: &[&str], point: &[f64]) -> Result<Vec<f64>, String> {
    if vars.len() != point.len() {
        return Err("Number of variables must match number of point dimensions".to_string());
    }

    let mut grad = Vec::with_capacity(vars.len());
    let h = 1e-6; // A small step for finite differences

    for i in 0..vars.len() {
        let mut point_plus_h = point.to_vec();
        point_plus_h[i] += h;

        let mut point_minus_h = point.to_vec();
        point_minus_h[i] -= h;

        let f_plus_h = eval_at_point(f, vars, &point_plus_h)?;
        let f_minus_h = eval_at_point(f, vars, &point_minus_h)?;

        let partial_deriv = (f_plus_h - f_minus_h) / (2.0 * h);
        grad.push(partial_deriv);
    }

    Ok(grad)
}

/// Helper to evaluate a multivariate expression at a point.
///
/// This function substitutes the numerical values from `point` into the `vars`
/// of the `expr` and then numerically evaluates the resulting expression.
///
/// # Arguments
/// * `expr` - The expression to evaluate.
/// * `vars` - The variables of the expression.
/// * `point` - The numerical values for the variables.
///
/// # Returns
/// A `Result` containing the numerical value of the expression, or an error string.
pub(crate) fn eval_at_point(expr: &Expr, vars: &[&str], point: &[f64]) -> Result<f64, String> {
    let mut vars_map = HashMap::new();
    for (i, &var) in vars.iter().enumerate() {
        vars_map.insert(var.to_string(), point[i]);
    }
    eval_expr(expr, &vars_map)
}

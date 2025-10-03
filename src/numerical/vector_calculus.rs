//! # Numerical Vector Calculus
//!
//! This module provides numerical implementations of vector calculus operations.
//! It includes functions for computing the gradient of a scalar field, and the
//! divergence and curl of a vector field, using finite difference approximations.

use crate::numerical::elementary::eval_expr;
use crate::symbolic::core::Expr;
use std::collections::HashMap;

/// Computes the numerical gradient of a scalar field `f` at a given point.
///
/// The scalar field is represented by a symbolic expression. The gradient is approximated
/// using finite differences.
///
/// # Arguments
/// * `f` - The symbolic expression representing the scalar field.
/// * `vars` - A slice of string slices representing the independent variables.
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

/// Computes the numerical divergence of a vector field at a given point.
///
/// The vector field is represented by a closure `F(&[f64]) -> Result<Vec<f64>, String>`
/// that returns the vector components at a given point. Divergence is approximated
/// using finite differences.
/// Assumes Cartesian coordinates.
///
/// # Arguments
/// * `vector_field` - A closure representing the vector field `F`.
/// * `point` - The point at which to compute the divergence.
///
/// # Returns
/// A `Result` containing the numerical divergence, or an error string.
pub fn divergence<F>(vector_field: F, point: &[f64]) -> Result<f64, String>
where
    F: Fn(&[f64]) -> Result<Vec<f64>, String>,
{
    let dim = point.len();
    let h = 1e-6;
    let mut div = 0.0;

    for i in 0..dim {
        let mut point_plus_h = point.to_vec();
        point_plus_h[i] += h;

        let mut point_minus_h = point.to_vec();
        point_minus_h[i] -= h;

        let f_plus_h = vector_field(&point_plus_h)?;
        let f_minus_h = vector_field(&point_minus_h)?;

        // Partial derivative of the i-th component with respect to the i-th variable
        let partial_deriv = (f_plus_h[i] - f_minus_h[i]) / (2.0 * h);
        div += partial_deriv;
    }

    Ok(div)
}

/// Computes the numerical curl of a 3D vector field at a given point.
///
/// The vector field is represented by a closure `F(&[f64]) -> Result<Vec<f64>, String>`
/// that returns the vector components at a given point. Curl is approximated
/// using finite differences.
/// Assumes Cartesian coordinates.
///
/// # Arguments
/// * `vector_field` - A closure representing the vector field `F`.
/// * `point` - The point at which to compute the curl.
///
/// # Returns
/// A `Result` containing the numerical curl vector, or an error string.
pub fn curl<F>(vector_field: F, point: &[f64]) -> Result<Vec<f64>, String>
where
    F: Fn(&[f64]) -> Result<Vec<f64>, String>,
{
    if point.len() != 3 {
        return Err("Curl is only defined for 3D vector fields.".to_string());
    }
    let h = 1e-6;

    // Partial derivatives needed for curl components
    let mut p_plus_h = point.to_vec();
    let mut p_minus_h = point.to_vec();

    // dVz/dy
    p_plus_h[1] += h;
    p_minus_h[1] -= h;
    let d_vz_dy = (vector_field(&p_plus_h)?[2] - vector_field(&p_minus_h)?[2]) / (2.0 * h);
    p_plus_h[1] = point[1];
    p_minus_h[1] = point[1]; // reset

    // dVy/dz
    p_plus_h[2] += h;
    p_minus_h[2] -= h;
    let d_vy_dz = (vector_field(&p_plus_h)?[1] - vector_field(&p_minus_h)?[1]) / (2.0 * h);
    p_plus_h[2] = point[2];
    p_minus_h[2] = point[2];

    // dVx/dz
    p_plus_h[2] += h;
    p_minus_h[2] -= h;
    let d_vx_dz = (vector_field(&p_plus_h)?[0] - vector_field(&p_minus_h)?[0]) / (2.0 * h);
    p_plus_h[2] = point[2];
    p_minus_h[2] = point[2];

    // dVz/dx
    p_plus_h[0] += h;
    p_minus_h[0] -= h;
    let d_vz_dx = (vector_field(&p_plus_h)?[2] - vector_field(&p_minus_h)?[2]) / (2.0 * h);
    p_plus_h[0] = point[0];
    p_minus_h[0] = point[0];

    // dVy/dx
    p_plus_h[0] += h;
    p_minus_h[0] -= h;
    let d_vy_dx = (vector_field(&p_plus_h)?[1] - vector_field(&p_minus_h)?[1]) / (2.0 * h);
    p_plus_h[0] = point[0];
    p_minus_h[0] = point[0];

    // dVx/dy
    p_plus_h[1] += h;
    p_minus_h[1] -= h;
    let d_vx_dy = (vector_field(&p_plus_h)?[0] - vector_field(&p_minus_h)?[0]) / (2.0 * h);

    let curl_x = d_vz_dy - d_vy_dz;
    let curl_y = d_vx_dz - d_vz_dx;
    let curl_z = d_vy_dx - d_vx_dy;

    Ok(vec![curl_x, curl_y, curl_z])
}

/// Helper to evaluate a multivariate expression at a point.
pub(crate) fn eval_at_point(expr: &Expr, vars: &[&str], point: &[f64]) -> Result<f64, String> {
    let mut vars_map = HashMap::new();
    for (i, &var) in vars.iter().enumerate() {
        vars_map.insert(var.to_string(), point[i]);
    }
    eval_expr(expr, &vars_map)
}

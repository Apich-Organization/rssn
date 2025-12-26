//! # Numerical Vector Calculus
//!
//! This module provides numerical implementations of vector calculus operations.
//! It includes functions for computing the gradient, divergence, curl, and Laplacian.

use crate::numerical::calculus::{eval_at_point, gradient as scalar_gradient};
use crate::symbolic::core::Expr;

/// Computes the numerical gradient of a scalar field `f` at a given point.
///
/// This is a re-export of `crate::numerical::calculus::gradient`.
///
/// # Arguments
/// * `f` - The symbolic expression representing the scalar field.
/// * `vars` - A slice of string slices representing the independent variables.
/// * `point` - The point at which to compute the gradient.
pub fn gradient(f: &Expr, vars: &[&str], point: &[f64]) -> Result<Vec<f64>, String> {
    scalar_gradient(f, vars, point)
}

/// Computes the numerical divergence of a vector field at a given point.
///
/// The vector field is represented by a closure `F(&[f64]) -> Result<Vec<f64>, String>`.
///
/// # Arguments
/// * `vector_field` - A closure representing the vector field `F`.
/// * `point` - The point at which to compute the divergence.
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
        let partial_deriv = (f_plus_h[i] - f_minus_h[i]) / (2.0 * h);
        div += partial_deriv;
    }
    Ok(div)
}

/// Computes the numerical divergence of a vector field represented by symbolic expressions.
///
/// # Arguments
/// * `funcs` - A slice of symbolic expressions representing the vector field components.
/// * `vars` - A slice of string slices representing the independent variables.
/// * `point` - The point at which to compute the divergence.
///
/// # Example
/// ```rust
/// use rssn::symbolic::core::Expr;
/// use rssn::numerical::vector_calculus::divergence_expr;
///
/// let x = Expr::new_variable("x");
/// let y = Expr::new_variable("y");
/// // F = [x^2, y^2] -> div F = 2x + 2y
/// let f1 = Expr::new_pow(x.clone(), Expr::new_constant(2.0));
/// let f2 = Expr::new_pow(y.clone(), Expr::new_constant(2.0));
/// let div = divergence_expr(&[f1, f2], &["x", "y"], &[1.0, 2.0]).unwrap();
/// assert!((div - 6.0).abs() < 1e-5);
/// ```
pub fn divergence_expr(funcs: &[Expr], vars: &[&str], point: &[f64]) -> Result<f64, String> {
    if funcs.len() != vars.len() {
        return Err("Number of functions must match number of variables".to_string());
    }
    let vector_field = |p: &[f64]| -> Result<Vec<f64>, String> {
        let mut res = Vec::with_capacity(funcs.len());
        for f in funcs {
            res.push(eval_at_point(f, vars, p)?);
        }
        Ok(res)
    };
    divergence(vector_field, point)
}

/// Computes the numerical curl of a 3D vector field at a given point.
///
/// # Arguments
/// * `vector_field` - A closure representing the vector field `F`.
/// * `point` - The point at which to compute the curl.
pub fn curl<F>(vector_field: F, point: &[f64]) -> Result<Vec<f64>, String>
where
    F: Fn(&[f64]) -> Result<Vec<f64>, String>,
{
    if point.len() != 3 {
        return Err("Curl is only defined for 3D vector fields.".to_string());
    }
    let h = 1e-6;
    let mut p_plus_h = point.to_vec();
    let mut p_minus_h = point.to_vec();

    // x component: (dVz/dy - dVy/dz)
    p_plus_h[1] += h;
    p_minus_h[1] -= h;
    let dvz_dy = (vector_field(&p_plus_h)?[2] - vector_field(&p_minus_h)?[2]) / (2.0 * h);
    p_plus_h[1] = point[1];
    p_minus_h[1] = point[1];

    p_plus_h[2] += h;
    p_minus_h[2] -= h;
    let dvy_dz = (vector_field(&p_plus_h)?[1] - vector_field(&p_minus_h)?[1]) / (2.0 * h);
    p_plus_h[2] = point[2];
    p_minus_h[2] = point[2];

    // y component: (dVx/dz - dVz/dx)
    p_plus_h[2] += h;
    p_minus_h[2] -= h;
    let dvx_dz = (vector_field(&p_plus_h)?[0] - vector_field(&p_minus_h)?[0]) / (2.0 * h);
    p_plus_h[2] = point[2];
    p_minus_h[2] = point[2];

    p_plus_h[0] += h;
    p_minus_h[0] -= h;
    let dvz_dx = (vector_field(&p_plus_h)?[2] - vector_field(&p_minus_h)?[2]) / (2.0 * h);
    p_plus_h[0] = point[0];
    p_minus_h[0] = point[0];

    // z component: (dVy/dx - dVx/dy)
    p_plus_h[0] += h;
    p_minus_h[0] -= h;
    let dvy_dx = (vector_field(&p_plus_h)?[1] - vector_field(&p_minus_h)?[1]) / (2.0 * h);
    p_plus_h[0] = point[0];
    p_minus_h[0] = point[0];

    p_plus_h[1] += h;
    p_minus_h[1] -= h;
    let dvx_dy = (vector_field(&p_plus_h)?[0] - vector_field(&p_minus_h)?[0]) / (2.0 * h);

    Ok(vec![dvz_dy - dvy_dz, dvx_dz - dvz_dx, dvy_dx - dvx_dy])
}

/// Computes the numerical curl of a 3D vector field represented by symbolic expressions.
pub fn curl_expr(funcs: &[Expr], vars: &[&str], point: &[f64]) -> Result<Vec<f64>, String> {
    if funcs.len() != 3 || vars.len() != 3 {
        return Err("Curl is only defined for 3D vector fields.".to_string());
    }
    let vector_field = |p: &[f64]| -> Result<Vec<f64>, String> {
        let mut res = Vec::with_capacity(3);
        for f in funcs {
            res.push(eval_at_point(f, vars, p)?);
        }
        Ok(res)
    };
    curl(vector_field, point)
}

/// Computes the numerical Laplacian of a scalar function at a given point.
///
/// # Example
/// ```rust
/// use rssn::symbolic::core::Expr;
/// use rssn::numerical::vector_calculus::laplacian;
///
/// let x = Expr::new_variable("x");
/// let y = Expr::new_variable("y");
/// // f = x^2 + y^2 -> laplacian = 2 + 2 = 4
/// let f = Expr::new_add(Expr::new_pow(x, Expr::new_constant(2.0)), Expr::new_pow(y, Expr::new_constant(2.0)));
/// let lap = laplacian(&f, &["x", "y"], &[1.0, 1.0]).unwrap();
/// assert!((lap - 4.0).abs() < 1e-5);
/// ```
pub fn laplacian(f: &Expr, vars: &[&str], point: &[f64]) -> Result<f64, String> {
    let mut lap = 0.0;
    let h = 1e-4;
    for i in 0..vars.len() {
        let mut p_plus = point.to_vec();
        p_plus[i] += h;
        let mut p_minus = point.to_vec();
        p_minus[i] -= h;
        let f_plus = eval_at_point(f, vars, &p_plus)?;
        let f_center = eval_at_point(f, vars, point)?;
        let f_minus = eval_at_point(f, vars, &p_minus)?;
        lap += (f_plus - 2.0 * f_center + f_minus) / (h * h);
    }
    Ok(lap)
}

/// Computes the numerical directional derivative of a function at a given point.
pub fn directional_derivative(
    f: &Expr,
    vars: &[&str],
    point: &[f64],
    direction: &[f64],
) -> Result<f64, String> {
    let grad = gradient(f, vars, point)?;
    if grad.len() != direction.len() {
        return Err("Direction vector dimension must match point dimension".to_string());
    }
    let mut dot = 0.0;
    let mut mag_sq = 0.0;
    for i in 0..direction.len() {
        dot += grad[i] * direction[i];
        mag_sq += direction[i] * direction[i];
    }
    if mag_sq == 0.0 {
        return Err("Direction vector cannot be zero".to_string());
    }
    Ok(dot / mag_sq.sqrt())
}

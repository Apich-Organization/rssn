//! # Numerical Coordinate Transformations
//!
//! This module provides numerical tools for coordinate transformations.
//! It supports converting points between various coordinate systems (Cartesian, Cylindrical, Spherical)
//! and computing numerical Jacobians of these transformations.

use crate::numerical::calculus::gradient;
use crate::numerical::matrix::Matrix;
use crate::symbolic::coordinates::{self, CoordinateSystem};
use crate::symbolic::core::Expr;
use std::collections::HashMap;

/// Transforms a numerical point from one coordinate system to another.
///
/// This function leverages symbolic transformation rules and then numerically evaluates
/// the resulting expressions. For performance-critical applications, consider using
/// `transform_point_pure`.
///
/// # Arguments
/// * `point` - The point to transform as a slice of `f64` coordinates.
/// * `from` - The `CoordinateSystem` of the input point.
/// * `to` - The target `CoordinateSystem`.
///
/// # Returns
/// A `Result` containing a `Vec<f64>` of the transformed coordinates, or an error string.
pub fn transform_point(
    point: &[f64],
    from: CoordinateSystem,
    to: CoordinateSystem,
) -> Result<Vec<f64>, String> {
    // This implementation leverages the symbolic transformation and then evaluates the result.
    // A more optimized version would perform direct numerical calculations.
    let point_expr: Vec<Expr> = point.iter().map(|&v| Expr::Constant(v)).collect();
    let transformed_expr = coordinates::transform_point(&point_expr, from, to)?;

    let mut result = Vec::new();
    for expr in transformed_expr {
        result.push(crate::numerical::elementary::eval_expr(
            &expr,
            &HashMap::new(),
        )?);
    }
    Ok(result)
}

/// Computes the numerical Jacobian matrix of a coordinate transformation at a specific point.
///
/// The Jacobian matrix `J` contains the partial derivatives of the new coordinates with
/// respect to the old coordinates: `J_ij = ∂(new_coord_i) / ∂(old_coord_j)`.
/// This is computed using symbolic differentiation followed by numerical evaluation.
///
/// # Arguments
/// * `from` - The source `CoordinateSystem`.
/// * `to` - The target `CoordinateSystem`.
/// * `at_point` - The point at which to evaluate the Jacobian.
///
/// # Returns
/// A `Result` containing a `Matrix<f64>` representing the Jacobian matrix, or an error string.
pub fn numerical_jacobian(
    from: CoordinateSystem,
    to: CoordinateSystem,
    at_point: &[f64],
) -> Result<Matrix<f64>, String> {
    let (from_vars, _, rules) = coordinates::get_transform_rules(from, to)?;
    let mut jacobian_rows = Vec::new();

    for rule in &rules {
        let grad = gradient(
            rule,
            &from_vars.iter().map(|s| s.as_str()).collect::<Vec<_>>(),
            at_point,
        )?;
        jacobian_rows.push(grad);
    }

    let rows = jacobian_rows.len();
    let cols = if rows > 0 { jacobian_rows[0].len() } else { 0 };
    Ok(Matrix::new(rows, cols, jacobian_rows.concat()))
}

// =====================================================================================
// region: Pure Numerical Implementations
// =====================================================================================

/// Transforms a numerical point using direct `f64` calculations for high performance.
///
/// This function provides a more optimized approach for coordinate transformations
/// by directly applying numerical formulas without intermediate symbolic manipulation.
///
/// # Arguments
/// * `point` - The point to transform as a slice of `f64` coordinates.
/// * `from` - The `CoordinateSystem` of the input point.
/// * `to` - The target `CoordinateSystem`.
///
/// # Returns
/// A `Result` containing a `Vec<f64>` of the transformed coordinates, or an error string.
pub fn transform_point_pure(
    point: &[f64],
    from: CoordinateSystem,
    to: CoordinateSystem,
) -> Result<Vec<f64>, String> {
    if from == to {
        return Ok(point.to_vec());
    }
    let cartesian_point = to_cartesian_pure(point, from)?;
    from_cartesian_pure(&cartesian_point, to)
}

pub(crate) fn to_cartesian_pure(point: &[f64], from: CoordinateSystem) -> Result<Vec<f64>, String> {
    /// Converts a numerical point from a given coordinate system to Cartesian coordinates.
    ///
    /// This is a helper function for `transform_point_pure`.
    ///
    /// # Arguments
    /// * `point` - The point to convert.
    /// * `from` - The source `CoordinateSystem`.
    ///
    /// # Returns
    /// A `Result` containing a `Vec<f64>` of the Cartesian coordinates, or an error string.
    match from {
        CoordinateSystem::Cartesian => Ok(point.to_vec()),
        CoordinateSystem::Cylindrical => {
            if point.len() != 3 {
                return Err("Cylindrical point must have 3 components (r, theta, z)".to_string());
            }
            let r = point[0];
            let theta = point[1];
            let z = point[2];
            let x = r * theta.cos();
            let y = r * theta.sin();
            Ok(vec![x, y, z])
        }
        CoordinateSystem::Spherical => {
            if point.len() != 3 {
                return Err("Spherical point must have 3 components (rho, theta, phi)".to_string());
            }
            let rho = point[0];
            let theta = point[1];
            let phi = point[2];
            let x = rho * phi.sin() * theta.cos();
            let y = rho * phi.sin() * theta.sin();
            let z = rho * phi.cos();
            Ok(vec![x, y, z])
        }
    }
}

pub(crate) fn from_cartesian_pure(point: &[f64], to: CoordinateSystem) -> Result<Vec<f64>, String> {
    /// Converts a numerical point from Cartesian coordinates to a given target coordinate system.
    ///
    /// This is a helper function for `transform_point_pure`.
    ///
    /// # Arguments
    /// * `point` - The Cartesian point to convert.
    /// * `to` - The target `CoordinateSystem`.
    ///
    /// # Returns
    /// A `Result` containing a `Vec<f64>` of the transformed coordinates, or an error string.
    match to {
        CoordinateSystem::Cartesian => Ok(point.to_vec()),
        CoordinateSystem::Cylindrical => {
            if point.len() < 2 {
                return Err("Cartesian point must have at least 2 components (x, y)".to_string());
            }
            let x = point[0];
            let y = point[1];
            let r = (x.powi(2) + y.powi(2)).sqrt();
            let theta = y.atan2(x);
            let mut result = vec![r, theta];
            if point.len() > 2 {
                result.push(point[2]); // Preserve z component
            }
            Ok(result)
        }
        CoordinateSystem::Spherical => {
            if point.len() != 3 {
                return Err("Cartesian point must have 3 components (x, y, z)".to_string());
            }
            let x = point[0];
            let y = point[1];
            let z = point[2];
            let rho = (x.powi(2) + y.powi(2) + z.powi(2)).sqrt();
            let theta = y.atan2(x);
            let phi = (z / rho).acos();
            Ok(vec![rho, theta, phi])
        }
    }
}

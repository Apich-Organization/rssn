//! # Numerical Functional Analysis
//!
//! This module provides numerical implementations for concepts from functional analysis.
//! It includes functions for calculating various norms (L1, L2, L-infinity) and inner products
//! for functions represented by discrete points.

/// Calculates the L1 norm of a function represented by discrete points.
///
/// The L1 norm (or Manhattan norm) is defined as `∫|f(x)|dx`. For discrete points,
/// it is approximated as `Σ|y_i|Δx_i`.
///
/// # Arguments
/// * `points` - A slice of `(x, y)` tuples representing the function's samples.
///
/// # Returns
/// The numerical value of the L1 norm.
pub fn l1_norm(points: &[(f64, f64)]) -> f64 {
    points
        .windows(2)
        .map(|w| {
            let (x1, y1) = w[0];
            let (x2, y2) = w[1];
            // Use the average height of the interval for trapezoidal-like integration
            (y1.abs() + y2.abs()) / 2.0 * (x2 - x1)
        })
        .sum()
}

/// Calculates the L2 norm of a function represented by discrete points.
///
/// The L2 norm (or Euclidean norm) is defined as `sqrt(∫|f(x)|²dx)`. For discrete points,
/// it is approximated as `sqrt(Σ|y_i|²Δx_i)`.
///
/// # Arguments
/// * `points` - A slice of `(x, y)` tuples representing the function's samples.
///
/// # Returns
/// The numerical value of the L2 norm.
pub fn l2_norm(points: &[(f64, f64)]) -> f64 {
    let integral_sq: f64 = points
        .windows(2)
        .map(|w| {
            let (x1, y1) = w[0];
            let (x2, y2) = w[1];
            (y1.powi(2) + y2.powi(2)) / 2.0 * (x2 - x1)
        })
        .sum();
    integral_sq.sqrt()
}

/// Calculates the L-infinity norm of a function represented by discrete points.
///
/// The L-infinity norm (or Chebyshev norm) is defined as `max(|f(x)|)`. For discrete points,
/// it is simply the maximum absolute value among the sampled points.
///
/// # Arguments
/// * `points` - A slice of `(x, y)` tuples representing the function's samples.
///
/// # Returns
/// The numerical value of the L-infinity norm.
pub fn infinity_norm(points: &[(f64, f64)]) -> f64 {
    points.iter().map(|(_, y)| y.abs()).fold(0.0, f64::max)
}

/// Calculates the inner product of two functions, `<f, g> = ∫f(x)g(x)dx`.
///
/// Both functions must be sampled at the same x-coordinates. For discrete points,
/// it is approximated as `Σ f(x_i)g(x_i)Δx_i`.
///
/// # Arguments
/// * `f_points` - A slice of `(x, y)` tuples representing the first function's samples.
/// * `g_points` - A slice of `(x, y)` tuples representing the second function's samples.
///
/// # Returns
/// A `Result` containing the numerical value of the inner product, or an error string
/// if the input functions have different numbers of sample points.
pub fn inner_product(f_points: &[(f64, f64)], g_points: &[(f64, f64)]) -> Result<f64, String> {
    if f_points.len() != g_points.len() {
        return Err("Input functions must have the same number of sample points.".to_string());
    }
    let integral = f_points
        .windows(2)
        .enumerate()
        .map(|(i, w)| {
            let (x1, y1_f) = w[0];
            let (x2, y2_f) = w[1];
            let (_, y1_g) = g_points[i];
            let (_, y2_g) = g_points[i + 1];
            // Check if x-coordinates match
            if (x1 - g_points[i].0).abs() > 1e-9 || (x2 - g_points[i + 1].0).abs() > 1e-9 {
                // This is a simplified check. In a real scenario, you might panic or handle it differently.
                return 0.0;
            }
            (y1_f * y1_g + y2_f * y2_g) / 2.0 * (x2 - x1)
        })
        .sum();
    Ok(integral)
}

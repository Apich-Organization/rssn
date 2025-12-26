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
#[must_use]
pub fn l1_norm(points: &[(f64, f64)]) -> f64 {
    points
        .windows(2)
        .map(|w| {
            let (x1, y1) = w[0];
            let (x2, y2) = w[1];
            ((y1.abs() + y2.abs()) / 2.0) * (x2 - x1)
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
#[must_use]
pub fn l2_norm(points: &[(f64, f64)]) -> f64 {
    let integral_sq: f64 = points
        .windows(2)
        .map(|w| {
            let (x1, y1) = w[0];
            let (x2, y2) = w[1];
            (y2 * y2 + y1 * y1) / 2.0 * (x2 - x1)
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
            if (x1 - g_points[i].0).abs() > 1e-9 || (x2 - g_points[i + 1].0).abs() > 1e-9 {
                return 0.0;
            }
            (y1_f * y1_g + y2_f * y2_g) / 2.0 * (x2 - x1)
        })
        .sum();
    Ok(integral)
}

/// Computes the projection of function `f` onto function `g`.
/// Both functions must be sampled at the same x-coordinates.
pub fn project(
    f_points: &[(f64, f64)],
    g_points: &[(f64, f64)],
) -> Result<Vec<(f64, f64)>, String> {
    let ip_fg = inner_product(f_points, g_points)?;
    let ip_gg = inner_product(g_points, g_points)?;

    if ip_gg.abs() < 1e-15 {
        return Ok(g_points.iter().map(|&(x, _)| (x, 0.0)).collect());
    }

    let coeff = ip_fg / ip_gg;
    Ok(g_points.iter().map(|&(x, y)| (x, y * coeff)).collect())
}

/// Normalizes a function relative to its L2 norm.
pub fn normalize(points: &[(f64, f64)]) -> Vec<(f64, f64)> {
    let n = l2_norm(points);
    if n.abs() < 1e-15 {
        return points.to_vec();
    }
    points.iter().map(|&(x, y)| (x, y / n)).collect()
}

/// Performs the Gram-Schmidt process to orthogonalize a set of functions.
/// All functions must have the same sampling points.
pub fn gram_schmidt(basis: &[Vec<(f64, f64)>]) -> Result<Vec<Vec<(f64, f64)>>, String> {
    let mut orthogonal_basis: Vec<Vec<(f64, f64)>> = Vec::new();
    for i in 0..basis.len() {
        let mut v = basis[i].clone();
        for u in &orthogonal_basis {
            let proj = project(&basis[i], u)?;
            for (j, (_, py)) in proj.into_iter().enumerate() {
                v[j].1 -= py;
            }
        }
        orthogonal_basis.push(v);
    }
    Ok(orthogonal_basis)
}

/// Performs the Gram-Schmidt process to orthonormalize a set of functions.
pub fn gram_schmidt_orthonormal(basis: &[Vec<(f64, f64)>]) -> Result<Vec<Vec<(f64, f64)>>, String> {
    let orthogonal_basis = gram_schmidt(basis)?;
    Ok(orthogonal_basis
        .into_iter()
        .map(|v| normalize(&v))
        .collect())
}

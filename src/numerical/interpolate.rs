//! # Numerical Interpolation
//!
//! This module provides numerical interpolation methods for constructing new data points
//! within the range of a discrete set of known data points. It includes implementations
//! for Lagrange interpolation, cubic spline interpolation, and Bézier and B-spline curves.

use crate::numerical::polynomial::Polynomial;

/// Constructs a Lagrange interpolating polynomial that passes through a given set of points.
///
/// Lagrange interpolation is a method of finding a polynomial that takes on certain values
/// at specified points. The resulting polynomial is unique and passes exactly through all given points.
///
/// # Arguments
/// * `points` - A slice of `(x, y)` tuples.
///
/// # Returns
/// A `Result` containing a `Polynomial` struct representing the interpolating polynomial,
/// or an error string if duplicate x-coordinates are found.
pub fn lagrange_interpolation(points: &[(f64, f64)]) -> Result<Polynomial, String> {
    if points.is_empty() {
        return Ok(Polynomial { coeffs: vec![0.0] });
    }

    let mut total_poly = Polynomial { coeffs: vec![0.0] };

    for (j, (xj, yj)) in points.iter().enumerate() {
        let mut basis_poly = Polynomial { coeffs: vec![1.0] };
        for (i, (xi, _)) in points.iter().enumerate() {
            if i == j {
                continue;
            }
            // basis_poly *= (x - xi) / (xj - xi)
            let numerator = Polynomial {
                coeffs: vec![1.0, -xi],
            }; // Represents (x - xi)
            let denominator = xj - xi;
            if denominator.abs() < 1e-9 {
                return Err(format!("Duplicate x-coordinates found: {}", xj));
            }
            basis_poly = basis_poly * (numerator / denominator);
        }
        total_poly = total_poly + (basis_poly * *yj);
    }

    Ok(total_poly)
}

/// Creates a cubic spline interpolator for a given set of points.
///
/// Cubic spline interpolation constructs a piecewise cubic polynomial that passes
/// through each data point, ensuring smoothness (continuity of first and second derivatives)
/// at the interior points. The input points must be sorted by their x-coordinates.
///
/// # Arguments
/// * `points` - A slice of `(x, y)` tuples. Must be sorted by x.
///
/// # Returns
/// A `Result` containing a closure `Box<dyn Fn(f64) -> f64>` that can be used to evaluate
/// the spline at any point, or an error string if fewer than two points are provided.
pub fn cubic_spline_interpolation(
    points: &[(f64, f64)],
) -> Result<Box<dyn Fn(f64) -> f64>, String> {
    let n = points.len();
    if n < 2 {
        return Err("At least two points are required for spline interpolation.".to_string());
    }

    let mut h = vec![0.0; n - 1];
    for i in 0..(n - 1) {
        h[i] = points[i + 1].0 - points[i].0;
    }

    let mut alpha = vec![0.0; n - 1];
    for i in 1..(n - 1) {
        alpha[i] = (3.0 / h[i]) * (points[i + 1].1 - points[i].1)
            - (3.0 / h[i - 1]) * (points[i].1 - points[i - 1].1);
    }

    let mut l = vec![1.0; n];
    let mut mu = vec![0.0; n];
    let mut z = vec![0.0; n];

    for i in 1..(n - 1) {
        l[i] = 2.0 * (points[i + 1].0 - points[i - 1].0) - h[i - 1] * mu[i - 1];
        mu[i] = h[i] / l[i];
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
    }

    let mut c = vec![0.0; n];
    let mut b = vec![0.0; n - 1];
    let mut d = vec![0.0; n - 1];

    for j in (0..(n - 1)).rev() {
        c[j] = z[j] - mu[j] * c[j + 1];
        b[j] = (points[j + 1].1 - points[j].1) / h[j] - h[j] * (c[j + 1] + 2.0 * c[j]) / 3.0;
        d[j] = (c[j + 1] - c[j]) / (3.0 * h[j]);
    }

    let points_owned: Vec<_> = points.to_vec();
    let spline = move |x: f64| -> f64 {
        let i = match points_owned.binary_search_by(|(px, _)| px.partial_cmp(&x).unwrap()) {
            Ok(idx) => idx,
            Err(idx) => (idx - 1).max(0),
        };

        if i >= n - 1 {
            // Handle evaluation at the last point or beyond
            return points_owned[n - 1].1;
        }

        let dx = x - points_owned[i].0;
        points_owned[i].1 + b[i] * dx + c[i] * dx.powi(2) + d[i] * dx.powi(3)
    };

    Ok(Box::new(spline))
}

/// Evaluates a point on a Bézier curve defined by a set of control points at parameter `t`.
///
/// This function uses De Casteljau's algorithm, a numerically stable and efficient method
/// for evaluating Bézier curves. The curve passes through the first and last control points
/// and is tangent to the segments connecting the first two and last two control points.
///
/// # Arguments
/// * `control_points` - A slice of points, where each point is a slice of `f64` coordinates.
/// * `t` - The parameter, typically in the range `[0, 1]`.
///
/// # Returns
/// The coordinates of the point on the curve.
pub fn bezier_curve(control_points: &[Vec<f64>], t: f64) -> Vec<f64> {
    if control_points.is_empty() {
        return vec![];
    }
    if control_points.len() == 1 {
        return control_points[0].clone();
    }

    let mut new_points = Vec::with_capacity(control_points.len() - 1);
    for i in 0..(control_points.len() - 1) {
        let p1 = &control_points[i];
        let p2 = &control_points[i + 1];
        let new_point: Vec<f64> = p1
            .iter()
            .zip(p2.iter())
            .map(|(&c1, &c2)| (1.0 - t) * c1 + t * c2)
            .collect();
        new_points.push(new_point);
    }
    bezier_curve(&new_points, t)
}

/// Evaluates a point on a B-spline curve.
///
/// B-splines are generalizations of Bézier curves, offering more flexibility and local control.
/// They are defined by control points, a degree, and a knot vector. This function uses
/// the Cox-de Boor recursion formula to evaluate the curve at a given parameter `t`.
///
/// # Arguments
/// * `control_points` - The control points.
/// * `degree` - The degree of the spline, `p`.
/// * `knots` - The knot vector, `U`.
/// * `t` - The parameter.
///
/// # Returns
/// The coordinates of the point on the curve.
pub fn b_spline(
    control_points: &[Vec<f64>],
    degree: usize,
    knots: &[f64],
    t: f64,
) -> Option<Vec<f64>> {
    let n = control_points.len() - 1;
    let m = knots.len() - 1;

    if degree > n || m != n + degree + 1 {
        return None; // Invalid input
    }

    // Find the knot span `i` in which `t` lies.
    let i = find_knot_span(n, degree, t, knots);

    // Compute the non-zero basis functions.
    let basis_vals = basis_functions(i, t, degree, knots);

    // Compute the curve point.
    let mut point = vec![0.0; control_points[0].len()];
    for j in 0..=degree {
        let pt_idx = i - degree + j;
        let p = &control_points[pt_idx];
        for d in 0..point.len() {
            point[d] += basis_vals[j] * p[d];
        }
    }
    Some(point)
}

/// Finds the knot span for a given parameter t.
pub(crate) fn find_knot_span(n: usize, p: usize, t: f64, knots: &[f64]) -> usize {
    if t >= knots[n + 1] {
        return n;
    }
    if t < knots[p] {
        return p;
    }

    let mut low = p;
    let mut high = n + 1;
    let mut mid = (low + high) / 2;

    while t < knots[mid] || t >= knots[mid + 1] {
        if t < knots[mid] {
            high = mid;
        } else {
            low = mid;
        }
        mid = (low + high) / 2;
    }
    mid
}

/// Computes the B-spline basis functions using the Cox-de Boor formula.
pub(crate) fn basis_functions(i: usize, t: f64, p: usize, knots: &[f64]) -> Vec<f64> {
    let mut n = vec![0.0; p + 1];
    n[0] = 1.0;

    for j in 1..=p {
        let mut left = vec![0.0; j + 1];
        let mut right = vec![0.0; j + 1];

        for r in 0..j {
            let den1 = knots[i + r + 1] - knots[i + r - j + 1];
            if den1 != 0.0 {
                left[r] = (t - knots[i + r - j + 1]) / den1;
            }
            let den2 = knots[i + r + 2] - knots[i + r - j + 2];
            if den2 != 0.0 {
                right[r] = (knots[i + r + 2] - t) / den2;
            }
        }

        let mut saved = 0.0;
        for r in 0..j {
            let temp = n[r] / (knots[i + r + j] - knots[i + r]);
            n[r] = saved + (knots[i + r + j] - t) * temp;
            saved = t * temp;
        }
        n[j] = saved;
    }
    n
}

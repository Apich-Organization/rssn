//! # Complex Analysis
//!
//! This module implements advanced concepts from complex analysis, with a primary
//! focus on analytic continuation along a path. It provides tools for estimating
//! the radius of convergence of Taylor series and performing continuation steps.

use crate::symbolic::core::Expr;
use crate::symbolic::series::{self, calculate_taylor_coefficients, taylor_series};

/// Represents the analytic continuation of a function along a path.
/// It is stored as a chain of series expansions, each centered at a point on the path.
#[derive(Debug, Clone)]
pub struct PathContinuation {
    pub var: String,
    pub order: usize,
    /// A vector of (center, series_expression) tuples.
    pub pieces: Vec<(Expr, Expr)>,
}

impl PathContinuation {
    /// Creates a new analytic continuation starting with a Taylor series for `func` centered at `start_point`.
    ///
    /// # Arguments
    /// * `func` - The initial function to continue.
    /// * `var` - The variable of the function (e.g., "z").
    /// * `start_point` - The center for the first series expansion.
    /// * `order` - The order of the Taylor series expansions.
    ///
    /// # Returns
    /// A new `PathContinuation` instance.
    pub fn new(func: &Expr, var: &str, start_point: &Expr, order: usize) -> Self {
        // The order for coefficient calculation needs to be higher to get a good estimate for the radius.
        let initial_series = taylor_series(func, var, start_point, order);
        PathContinuation {
            var: var.to_string(),
            order,
            pieces: vec![(start_point.clone(), initial_series)],
        }
    }

    /// Continues the function along a given path.
    ///
    /// This method iterates through the `path_points`, using the end of the previous
    /// segment as the start for the next. Before each step, it verifies that the next point
    /// is within the estimated radius of convergence of the current series expansion.
    /// If a point is outside the radius of convergence, the continuation fails.
    ///
    /// # Arguments
    /// * `path_points` - A vector of `Expr` representing the points to continue through.
    ///
    /// # Panics
    /// * Panics if the `pieces` vector is empty (if `new` was not called first).
    /// * Panics if the next point on the path is outside the estimated radius of convergence
    ///   of the current series.
    /// * Panics if the radius of convergence cannot be estimated.
    pub fn continue_along_path(&mut self, path_points: &[Expr]) {
        for next_point in path_points {
            let (last_center, last_series) = self
                .pieces
                .last()
                .expect("PathContinuation must be initialized with `new` before continuing.");

            // --- Convergence Radius Check ---
            // Estimate the radius of convergence for the current series.
            // We need a slightly higher order for a stable estimation.
            let radius = estimate_radius_of_convergence(last_series, &self.var, last_center, self.order + 5)
                .expect("Failed to estimate the radius of convergence. The series may be trivial or coefficients non-numeric.");

            // Calculate the distance to the next point.
            let distance = complex_distance(last_center, next_point)
                .expect("Failed to calculate distance between complex points. Ensure points are valid complex numbers or real numbers.");

            // Check if the next point is within the radius.
            if distance >= radius {
                panic!(
                    "Analytic continuation failed: The next point {} is outside the estimated radius of convergence ({}) of the series centered at {}.",
                    next_point,
                    radius,
                    last_center
                );
            }
            // --- End of Check ---

            // Use the existing single-step continuation function to get the next series.
            let next_series = series::analytic_continuation(
                last_series,
                &self.var,
                last_center,
                next_point,
                self.order,
            );

            self.pieces.push((next_point.clone(), next_series));
        }
    }

    /// Returns the final expression (Taylor series) after continuation to the last point.
    ///
    /// # Returns
    /// An `Option<&Expr>` containing the final series expression, or `None` if the continuation
    /// path is empty.
    pub fn get_final_expression(&self) -> Option<&Expr> {
        self.pieces.last().map(|(_, series)| series)
    }
}

/// Estimates the radius of convergence for a Taylor series using the ratio test on its coefficients.
/// R = 1 / lim |c_{n+1}/c_n|
///
/// # Arguments
/// * `series_expr` - The Taylor series expression.
/// * `var` - The variable of the series.
/// * `center` - The center of the series expansion.
/// * `order` - The order of coefficients to calculate for the estimation.
///
/// # Returns
/// An `Option<f64>` containing the estimated radius. Returns `None` if the radius cannot be determined
/// (e.g., if coefficients are not convertible to f64, or if there are not enough non-zero coefficients).
pub(crate) fn estimate_radius_of_convergence(
    series_expr: &Expr,
    var: &str,
    center: &Expr,
    order: usize,
) -> Option<f64> {
    let coeffs = calculate_taylor_coefficients(series_expr, var, center, order);

    // Find the last two consecutive non-zero coefficients to apply the ratio test.
    for n in (1..coeffs.len()).rev() {
        let cn = &coeffs[n];
        let cn_minus_1 = &coeffs[n - 1];

        if let (Some(c_n_val), Some(c_n_minus_1_val)) = (cn.to_f64(), cn_minus_1.to_f64()) {
            // Avoid division by zero if a coefficient is zero.
            if c_n_val.abs() > f64::EPSILON && c_n_minus_1_val.abs() > f64::EPSILON {
                let limit_ratio = c_n_val.abs() / c_n_minus_1_val.abs();
                if limit_ratio < f64::EPSILON {
                    // If the limit of ratios is zero, the radius is infinite.
                    return Some(f64::INFINITY);
                } else {
                    return Some(1.0 / limit_ratio);
                }
            }
        }
    }

    // If no two consecutive non-zero coefficients were found, the radius is effectively infinite
    // for a polynomial, or the estimation failed.
    Some(f64::INFINITY)
}

/// Calculates the Euclidean distance between two points represented as `Expr`.
/// The points can be real (`Expr::Constant`) or complex (`Expr::Complex`).
///
/// # Arguments
/// * `p1` - The first point.
/// * `p2` - The second point.
///
/// # Returns
/// An `Option<f64>` for the distance. Returns `None` if the points are not valid numbers.
pub(crate) fn complex_distance(p1: &Expr, p2: &Expr) -> Option<f64> {
    let re1 = p1.re().to_f64().unwrap_or(0.0);
    let im1 = p1.im().to_f64().unwrap_or(0.0);
    let re2 = p2.re().to_f64().unwrap_or(0.0);
    let im2 = p2.im().to_f64().unwrap_or(0.0);

    let dx = re1 - re2;
    let dy = im1 - im2;

    Some((dx * dx + dy * dy).sqrt())
}

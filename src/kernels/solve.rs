//! # Numerical Equation Solvers
//!
//! This module provides numerical methods for solving linear and non-linear systems of equations.
//! It includes functions for solving linear systems using Gaussian elimination (via RREF)
//! and non-linear systems using Newton's method. Functions are plain closures, so
//! callers decide how they are evaluated (interpreted, compiled, hand-written).

use serde::Deserialize;
use serde::Serialize;

use crate::kernels::calculus::jacobian;
use crate::kernels::exact_eq;
use crate::kernels::matrix::Matrix;

/// Represents the solution to a system of linear equations.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum LinearSolution {
    /// A single unique solution vector.
    Unique(Vec<f64>),
    /// A parametric solution (particular + null space basis).
    Parametric {
        /// A particular solution vector.
        particular: Vec<f64>,
        /// A basis for the null space of the matrix.
        null_space_basis: Matrix<f64>,
    },
    /// The system is inconsistent and has no solution.
    NoSolution,
}

/// Solves a system of linear equations `Ax = b`.
///
/// This function constructs an augmented matrix `[A | b]`, computes its Reduced Row Echelon Form (RREF),
/// and then analyzes the RREF to determine the nature of the solution:
/// - **Unique Solution**: Returns a unique solution vector.
/// - **Parametric Solution**: Returns a particular solution and a basis for the null space.
/// - **No Solution**: Indicates an inconsistent system.
///
/// # Arguments
/// * `a` - The coefficient matrix `A`.
/// * `b` - The constant vector `b`.
///
/// # Returns
/// A `Result` containing a `LinearSolution` enum, or an error string.
///
/// # Errors
/// Returns an error if the matrix and vector dimensions are incompatible, or if RREF/null space computation fails.
pub fn solve_linear_system(
    a: &Matrix<f64>,
    b: &[f64],
) -> Result<LinearSolution, String> {
    let (rows, cols) = (a.rows(), a.cols());

    if rows != b.len() {
        return Err("Matrix and \
                    vector dimensions \
                    are incompatible.\
                    "
        .to_string());
    }

    let mut augmented_data = vec![0.0; rows * (cols + 1)];

    for i in 0..rows {
        for j in 0..cols {
            augmented_data[i * (cols + 1) + j] = *a.get(i, j);
        }

        augmented_data[i * (cols + 1) + cols] = b[i];
    }

    let mut augmented = Matrix::new(rows, cols + 1, augmented_data);

    let rank = augmented.rref()?;

    // Check for inconsistency: if any row has a leading 1 in the last column (the constant vector column)
    for i in 0..rank {
        let mut pivot_col = 0;

        while pivot_col < cols + 1 && augmented.get(i, pivot_col).abs() < 1e-9 {
            pivot_col += 1;
        }

        if pivot_col == cols {
            return Ok(LinearSolution::NoSolution);
        }
    }

    if rank < cols {
        let mut particular = vec![0.0; cols];

        #[warn(clippy::collection_is_never_read)]
        let mut lead = 0;

        for r in 0..rank {
            let mut i = lead;

            while i < cols && augmented.get(r, i).abs() < 1e-9 {
                i += 1;
            }

            if i < cols {
                particular[i] = *augmented.get(r, cols);

                lead = i + 1;
            }
        }

        let null_space = a.null_space()?;

        Ok(LinearSolution::Parametric {
            particular,
            null_space_basis: null_space,
        })
    } else {
        let mut solution = vec![0.0; cols];

        for (i, var) in solution.iter_mut().enumerate().take(rank) {
            *var = *augmented.get(i, cols);
        }

        Ok(LinearSolution::Unique(solution))
    }
}

/// Solves the square nonlinear system `f(x) = 0` by Newton's method with a
/// finite-difference Jacobian.
///
/// `f` writes the residuals for the point given as its first argument into
/// its second argument; both have the length of `start`.
///
/// # Errors
/// Returns an error when a Jacobian is singular, a residual is not finite,
/// or the iteration does not converge within `max_iter` steps.
pub fn solve_nonlinear_system(
    f: impl Fn(&[f64], &mut [f64]),
    start: &[f64],
    tolerance: f64,
    max_iter: usize,
) -> Result<Vec<f64>, String> {
    let n = start.len();
    let mut x = start.to_vec();
    let mut residual = vec![0.0; n];
    for _ in 0..max_iter {
        f(&x, &mut residual);
        if residual.iter().any(|r| !r.is_finite()) {
            return Err("Residual is not finite.".to_string());
        }
        let jac = Matrix::new(n, n, jacobian(&f, &x, n));
        let rhs: Vec<f64> = residual.iter().map(|r| -r).collect();
        let LinearSolution::Unique(delta) = solve_linear_system(&jac, &rhs)? else {
            return Err("Jacobian is singular; Newton's method failed.".to_string());
        };
        for (xi, di) in x.iter_mut().zip(&delta) {
            *xi += di;
        }
        if delta.iter().map(|d| d * d).sum::<f64>().sqrt() < tolerance {
            return Ok(x);
        }
    }
    Err("Newton's method did not converge.".to_string())
}

/// Finds a root of `f` by Newton's method, using `df` for the derivative.
///
/// # Errors
/// Returns an error when the derivative vanishes or the iteration does not
/// converge.
pub fn solve_root_newton(
    f: impl Fn(f64) -> f64,
    df: impl Fn(f64) -> f64,
    start: f64,
    tolerance: f64,
    max_iter: usize,
) -> Result<f64, String> {
    let mut x = start;
    for _ in 0..max_iter {
        let fx = f(x);
        if !fx.is_finite() {
            return Err("Function value is not finite.".to_string());
        }
        if fx.abs() < tolerance {
            return Ok(x);
        }
        let slope = df(x);
        if slope.abs() < 1e-14 {
            return Err("Derivative too close to zero in Newton's method.".to_string());
        }
        let delta = fx / slope;
        x -= delta;
        if delta.abs() < tolerance {
            return Ok(x);
        }
    }
    if f(x).abs() < tolerance * 100.0 {
        Ok(x)
    } else {
        Err("Newton's method did not converge.".to_string())
    }
}

/// Finds a root of `f` in `interval` by bisection.
///
/// Endpoints may be in either order. A finite endpoint or midpoint whose
/// residual is zero or has magnitude less than `tolerance` is accepted.
/// Otherwise the absolute interval width is compared with the tolerance.
/// Zero tolerance disables approximate stopping. The two endpoint evaluations
/// are followed by at most `max_iter` midpoint evaluations. At the budget or
/// when no representable point remains inside the bracket, returns its midpoint.
///
/// # Errors
/// Returns an error for non-finite endpoints or function values, a negative
/// or non-finite tolerance, or endpoints with the same sign outside the tolerance.
pub fn solve_root_bisection(
    f: impl Fn(f64) -> f64,
    interval: (f64, f64),
    tolerance: f64,
    max_iter: usize,
) -> Result<f64, String> {
    let (a, b) = interval;
    if !a.is_finite() || !b.is_finite() {
        return Err("Root interval endpoints must be finite.".to_string());
    }
    check_root_tolerance(tolerance)?;
    let (fa, fb) = (f(a), f(b));
    if !fa.is_finite() || !fb.is_finite() {
        return Err("Function value is not finite.".to_string());
    }
    // With zero tolerance, a strict comparison admits only signed zeros.
    let residual_tolerance = tolerance.max(f64::from_bits(1));
    if fa.abs() < residual_tolerance {
        return Ok(a);
    }
    if fb.abs() < residual_tolerance {
        return Ok(b);
    }
    let sign_at_a = fa.is_sign_negative();
    if sign_at_a == fb.is_sign_negative() {
        return Err(format!("Interval [{a}, {b}] does not bracket a root"));
    }
    bisect_sampled(f, interval, sign_at_a, max_iter, |value, width, _| {
        value.abs() < residual_tolerance || width < tolerance
    })
    .ok_or_else(|| "Function value is not finite.".to_string())
}

/// Bisects finite endpoints with already evaluated opposite, nonzero signs.
///
/// The stopping predicate receives the midpoint residual, absolute bracket
/// width and midpoint coordinate. A non-finite midpoint residual gives `None`.
pub(super) fn bisect_sampled(
    mut f: impl FnMut(f64) -> f64,
    (mut a, mut b): (f64, f64),
    sign_at_a: bool,
    max_iter: usize,
    converged: impl Fn(f64, f64, f64) -> bool,
) -> Option<f64> {
    // Subsequent endpoints stay within the original interval.
    let half_max = f64::MAX * 0.5;
    let direct_midpoint = a.abs() <= half_max && b.abs() <= half_max;
    #[allow(
        clippy::manual_midpoint,
        reason = "both endpoints stay within the prechecked half range"
    )]
    let midpoint = |a, b| {
        if direct_midpoint {
            (a + b) * 0.5
        } else {
            f64::midpoint(a, b)
        }
    };
    for _ in 0..max_iter {
        let mid = midpoint(a, b);
        if exact_eq(mid, a) || exact_eq(mid, b) {
            return Some(mid);
        }
        let fm = f(mid);
        if !fm.is_finite() {
            return None;
        }
        if converged(fm, (b - a).abs(), mid) {
            return Some(mid);
        }
        if fm.is_sign_negative() == sign_at_a {
            a = mid;
        } else {
            b = mid;
        }
    }
    Some(midpoint(a, b))
}

fn check_root_tolerance(tolerance: f64) -> Result<(), String> {
    if tolerance.is_finite() && tolerance >= 0.0 {
        Ok(())
    } else {
        Err("Root tolerance must be finite and nonnegative.".to_string())
    }
}

/// Finds a root of `f` near `guess`: Newton's method with a central
/// difference derivative from several starting points, then bisection over
/// growing brackets. A non-finite sample discards the current attempt.
///
/// # Errors
/// Returns an error for a non-finite guess or a negative or non-finite
/// tolerance, or when every strategy fails.
pub fn solve_root(
    f: impl Fn(f64) -> f64,
    guess: f64,
    tolerance: f64,
    max_iter: usize,
) -> Result<f64, String> {
    if !guess.is_finite() {
        return Err("Root guess must be finite.".to_string());
    }
    check_root_tolerance(tolerance)?;
    let finite_sample = |x: f64| {
        if x.is_finite() {
            f(x)
        } else {
            f64::NAN
        }
    };
    let df = |x: f64| {
        let h = 1e-6 * x.abs().max(1.0);
        let slope = (finite_sample(x + h) - finite_sample(x - h)) / (2.0 * h);
        // A non-finite derivative cannot establish a converged zero step.
        if slope.is_finite() {
            slope
        } else {
            f64::NAN
        }
    };
    for start in [guess, 0.0, -1.0, 2.0, 0.5, -0.5, 10.0] {
        if let Ok(root) = solve_root_newton(finite_sample, df, start, tolerance, max_iter)
            && root.is_finite()
        {
            return Ok(root);
        }
    }
    for span in [2.0, 5.0, 10.0, 50.0, 100.0] {
        if let Ok(root) =
            solve_root_bisection(&f, (guess - span, guess + span), tolerance, max_iter)
        {
            return Ok(root);
        }
    }
    Err("Failed to find numerical root.".to_string())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn newton_system() {
        // x^2 + y^2 = 4, x*y = 1
        let f = |x: &[f64], out: &mut [f64]| {
            out[0] = x[0] * x[0] + x[1] * x[1] - 4.0;
            out[1] = x[0] * x[1] - 1.0;
        };
        let root = solve_nonlinear_system(f, &[2.0, 0.5], 1e-12, 50).unwrap_or_default();
        assert!((root[0] * root[0] + root[1] * root[1] - 4.0).abs() < 1e-9);
        assert!((root[0] * root[1] - 1.0).abs() < 1e-9);
    }

    #[test]
    fn scalar_roots() {
        let f = |x: f64| x * x * x - 2.0 * x - 5.0;
        let root = solve_root(f, 2.0, 1e-12, 100).unwrap_or(f64::NAN);
        assert!(f(root).abs() < 1e-9);
        let root = solve_root_bisection(f64::cos, (1.0, 2.0), 1e-12, 200).unwrap_or(f64::NAN);
        assert!((root - std::f64::consts::FRAC_PI_2).abs() < 1e-9);
        assert!(solve_root_bisection(|x| x * x + 1.0, (-1.0, 1.0), 1e-12, 50).is_err());
        assert!(solve_root(|x| x * x + 1.0, 0.3, 1e-12, 50).is_err());
    }
}

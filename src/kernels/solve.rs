//! # Numerical Equation Solvers
//!
//! This module provides numerical methods for solving linear and non-linear systems of equations.
//! It includes functions for solving linear systems using Gaussian elimination (via RREF)
//! and non-linear systems using Newton's method. Functions are plain closures, so
//! callers decide how they are evaluated (interpreted, compiled, hand-written).

use serde::Deserialize;
use serde::Serialize;

use crate::kernels::calculus::jacobian;
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
        let mut _pivot_cols = Vec::new();

        let mut lead = 0;

        for r in 0..rank {
            let mut i = lead;

            while i < cols && augmented.get(r, i).abs() < 1e-9 {
                i += 1;
            }

            if i < cols {
                _pivot_cols.push(i);

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
/// # Errors
/// Returns an error when the end points do not bracket a sign change.
pub fn solve_root_bisection(
    f: impl Fn(f64) -> f64,
    interval: (f64, f64),
    tolerance: f64,
    max_iter: usize,
) -> Result<f64, String> {
    let (mut a, mut b) = interval;
    let (mut fa, fb) = (f(a), f(b));
    if fa.abs() < tolerance {
        return Ok(a);
    }
    if fb.abs() < tolerance {
        return Ok(b);
    }
    if !(fa * fb < 0.0) {
        return Err(format!("Interval [{a}, {b}] does not bracket a root"));
    }
    for _ in 0..max_iter {
        let mid = f64::midpoint(a, b);
        let fm = f(mid);
        if fm.abs() < tolerance || (b - a).abs() < tolerance {
            return Ok(mid);
        }
        if fa * fm < 0.0 {
            b = mid;
        } else {
            a = mid;
            fa = fm;
        }
    }
    Ok(f64::midpoint(a, b))
}

/// Finds a root of `f` near `guess`: Newton's method with a central
/// difference derivative from several starting points, then bisection over
/// growing brackets.
///
/// # Errors
/// Returns an error when every strategy fails.
pub fn solve_root(
    f: impl Fn(f64) -> f64,
    guess: f64,
    tolerance: f64,
    max_iter: usize,
) -> Result<f64, String> {
    let df = |x: f64| {
        let h = 1e-6 * x.abs().max(1.0);
        (f(x + h) - f(x - h)) / (2.0 * h)
    };
    for start in [guess, 0.0, -1.0, 2.0, 0.5, -0.5, 10.0] {
        if let Ok(root) = solve_root_newton(&f, df, start, tolerance, max_iter) {
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

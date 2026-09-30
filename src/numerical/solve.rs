//! # Numerical Equation Solvers
//!
//! This module provides numerical methods for solving linear and non-linear systems of equations.
//! It includes functions for solving linear systems using Gaussian elimination (via RREF)
//! and non-linear systems using Newton's method.

use std::collections::HashMap;

use serde::Deserialize;
use serde::Serialize;

use crate::numerical::calculus::gradient;
use crate::numerical::elementary::eval_expr;
use crate::numerical::matrix::Matrix;
use crate::symbolic::core::Expr;

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

/// Solves a system of non-linear equations `F(X) = 0` using Newton's method.
///
/// Newton's method for systems of equations is an iterative algorithm that approximates
/// the solution by repeatedly solving a linear system involving the Jacobian matrix.
///
/// # Arguments
/// * `funcs` - A slice of expressions representing the functions in the system.
/// * `vars` - The variables to solve for.
/// * `start_point` - An initial guess for the solution vector.
/// * `tolerance` - The desired precision of the solution.
/// * `max_iter` - The maximum number of iterations.
///
/// # Returns
/// A `Result` containing the solution vector, or an error string.
///
/// # Errors
/// Returns an error if symbolic expression evaluation fails, if the Jacobian is singular, or if the method fails to converge.
pub fn solve_nonlinear_system(
    funcs: &[Expr],
    vars: &[&str],
    start_point: &[f64],
    tolerance: f64,
    max_iter: usize,
) -> Result<Vec<f64>, String> {
    let mut x_n = start_point.to_vec();

    for _ in 0..max_iter {
        let mut f_at_x = Vec::new();

        for func in funcs {
            let mut vars_map = HashMap::new();

            for (i, &var) in vars.iter().enumerate() {
                vars_map.insert(var.to_string(), x_n[i]);
            }

            f_at_x.push(eval_expr(func, &vars_map)?);
        }

        let mut jacobian_rows = Vec::new();

        for func in funcs {
            jacobian_rows.push(gradient(func, vars, &x_n)?);
        }

        let jacobian = Matrix::new(funcs.len(), vars.len(), jacobian_rows.concat());

        let neg_f: Vec<f64> = f_at_x.iter().map(|v| -v).collect();

        let delta_x = match solve_linear_system(&jacobian, &neg_f)? {
            | LinearSolution::Unique(sol) => sol,
            | _ => return Err("Jacobian is singular; Newton's method failed.".to_string()),
        };

        for i in 0..x_n.len() {
            x_n[i] += delta_x[i];
        }

        let norm_delta = delta_x.iter().map(|v| v * v).sum::<f64>().sqrt();

        if norm_delta < tolerance {
            return Ok(x_n);
        }
    }

    Err("Newton's method did not \
         converge."
        .to_string())
}

/// Solves for a single root of `func = 0` with respect to `var` using the Newton-Raphson method.
///
/// # Arguments
/// * `func` - An expression representing the function whose root is sought.
/// * `var` - The variable name to solve for.
/// * `start_point` - The initial guess.
/// * `tolerance` - The desired error tolerance.
/// * `max_iter` - The maximum number of iterations.
///
/// # Errors
/// Returns an error if evaluation fails, derivative is near zero, or convergence is not achieved.
pub fn solve_root_newton(
    func: &Expr,
    var: &str,
    start_point: f64,
    tolerance: f64,
    max_iter: usize,
) -> Result<f64, String> {
    let mut x_n = start_point;
    let mut vars_map = HashMap::new();

    for _ in 0..max_iter {
        vars_map.insert(var.to_string(), x_n);
        let f_val = eval_expr(func, &vars_map)?;

        if f_val.abs() < tolerance {
            return Ok(x_n);
        }

        // Numerical derivative via central differences
        let h = 1e-6;
        vars_map.insert(var.to_string(), x_n + h);
        let f_plus = eval_expr(func, &vars_map)?;
        vars_map.insert(var.to_string(), x_n - h);
        let f_minus = eval_expr(func, &vars_map)?;
        let df = (f_plus - f_minus) / (2.0 * h);

        if df.abs() < 1e-12 {
            return Err("Derivative too close to zero in Newton's method.".to_string());
        }

        let delta = f_val / df;
        x_n -= delta;

        if delta.abs() < tolerance {
            return Ok(x_n);
        }
    }

    vars_map.insert(var.to_string(), x_n);
    let final_val = eval_expr(func, &vars_map)?;
    if final_val.abs() < tolerance * 100.0 {
        Ok(x_n)
    } else {
        Err("Newton's method did not converge.".to_string())
    }
}

/// Solves for a single root of `func = 0` with respect to `var` using the bisection method.
///
/// # Arguments
/// * `func` - An expression representing the function.
/// * `var` - The variable name.
/// * `interval` - `(a, b)` bracketing interval such that `f(a) * f(b) <= 0`.
/// * `tolerance` - Convergence tolerance.
/// * `max_iter` - Maximum iterations.
pub fn solve_root_bisection(
    func: &Expr,
    var: &str,
    interval: (f64, f64),
    tolerance: f64,
    max_iter: usize,
) -> Result<f64, String> {
    let (mut a, mut b) = interval;
    let mut vars_map = HashMap::new();

    vars_map.insert(var.to_string(), a);
    let fa = eval_expr(func, &vars_map)?;
    vars_map.insert(var.to_string(), b);
    let fb = eval_expr(func, &vars_map)?;

    if fa.abs() < tolerance {
        return Ok(a);
    }
    if fb.abs() < tolerance {
        return Ok(b);
    }

    if fa * fb > 0.0 {
        return Err(format!("Interval [{a}, {b}] does not bracket a root"));
    }

    for _ in 0..max_iter {
        let mid = f64::midpoint(a, b);
        vars_map.insert(var.to_string(), mid);
        let f_mid = eval_expr(func, &vars_map)?;

        if f_mid.abs() < tolerance || (b - a).abs() < tolerance {
            return Ok(mid);
        }

        vars_map.insert(var.to_string(), a);
        let fa = eval_expr(func, &vars_map)?;
        if fa * f_mid < 0.0 {
            b = mid;
        } else {
            a = mid;
        }
    }

    Ok(f64::midpoint(a, b))
}

/// Universal 1D root finder combining Newton-Raphson and Bisection fallbacks.
pub fn solve_root(
    func: &Expr,
    var: &str,
    start_guess: Option<f64>,
    tolerance: f64,
    max_iter: usize,
) -> Result<f64, String> {
    let guess = start_guess.unwrap_or(1.0);

    // 1. Try Newton's method from guess
    if let Ok(root) = solve_root_newton(func, var, guess, tolerance, max_iter) {
        return Ok(root);
    }

    // 2. Try alternative starting points
    for alt_guess in [0.0, -1.0, 2.0, 0.5, -0.5, 10.0] {
        if (alt_guess - guess).abs() > 1e-5 {
            if let Ok(root) = solve_root_newton(func, var, alt_guess, tolerance, max_iter) {
                return Ok(root);
            }
        }
    }

    // 3. Try Bisection search over brackets around guess
    for span in [2.0, 5.0, 10.0, 50.0, 100.0] {
        let interval = (guess - span, guess + span);
        if let Ok(root) = solve_root_bisection(func, var, interval, tolerance, max_iter) {
            return Ok(root);
        }
    }

    Err("Failed to find numerical root.".to_string())
}


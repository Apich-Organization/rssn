//! # Numerical Equation Solvers
//!
//! This module provides numerical methods for solving linear and non-linear systems of equations.
//! It includes functions for solving linear systems using Gaussian elimination (via RREF)
//! and non-linear systems using Newton's method.

use crate::numerical::calculus::gradient;
use crate::numerical::elementary::eval_expr;
use crate::numerical::matrix::Matrix;
use crate::prelude::Expr;
use std::collections::HashMap;

/// Represents the solution to a system of linear equations.
pub enum LinearSolution {
    Unique(Vec<f64>),
    Parametric {
        particular: Vec<f64>,
        null_space_basis: Matrix<f64>,
    },
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
pub fn solve_linear_system(a: &Matrix<f64>, b: &Vec<f64>) -> Result<LinearSolution, String> {
    let (rows, cols) = (a.rows(), a.cols());
    if rows != b.len() {
        return Err("Matrix and vector dimensions are incompatible.".to_string());
    }

    let mut augmented_data = vec![0.0; rows * (cols + 1)];
    for i in 0..rows {
        for j in 0..cols {
            augmented_data[i * (cols + 1) + j] = *a.get(i, j);
        }
        augmented_data[i * (cols + 1) + cols] = b[i];
    }
    let mut augmented = Matrix::new(rows, cols + 1, augmented_data);

    let rank = augmented.rref();

    for i in rank..rows {
        if augmented.get(i, cols).abs() > 1e-9 {
            return Ok(LinearSolution::NoSolution);
        }
    }

    if rank < cols {
        // Infinite solutions
        let mut particular = vec![0.0; cols];
        let mut pivot_cols = Vec::new();
        let mut lead = 0;
        for r in 0..rank {
            let mut i = lead;
            while i < cols && augmented.get(r, i).abs() < 1e-9 {
                i += 1;
            }
            if i < cols {
                pivot_cols.push(i);
                particular[i] = *augmented.get(r, cols);
                lead = i + 1;
            }
        }

        let null_space = a.null_space();
        Ok(LinearSolution::Parametric {
            particular,
            null_space_basis: null_space,
        })
    } else {
        // Unique solution
        let mut solution = vec![0.0; cols];
        for i in 0..rank {
            solution[i] = *augmented.get(i, cols);
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
pub fn solve_nonlinear_system(
    funcs: &[Expr],
    vars: &[&str],
    start_point: &[f64],
    tolerance: f64,
    max_iter: usize,
) -> Result<Vec<f64>, String> {
    let mut x_n = start_point.to_vec();

    for _ in 0..max_iter {
        // 1. Evaluate F(X_n)
        let mut f_at_x = Vec::new();
        for func in funcs {
            let mut vars_map = HashMap::new();
            for (i, &var) in vars.iter().enumerate() {
                vars_map.insert(var.to_string(), x_n[i]);
            }
            f_at_x.push(eval_expr(func, &vars_map)?);
        }

        // 2. Compute Jacobian matrix J(X_n)
        let mut jacobian_rows = Vec::new();
        for func in funcs {
            jacobian_rows.push(gradient(func, vars, &x_n)?);
        }
        let jacobian = Matrix::new(funcs.len(), vars.len(), jacobian_rows.concat());

        // 3. Solve the linear system J * delta_x = -F
        let neg_f: Vec<f64> = f_at_x.iter().map(|v| -v).collect();
        let delta_x = match solve_linear_system(&jacobian, &neg_f)? {
            LinearSolution::Unique(sol) => sol,
            _ => return Err("Jacobian is singular; Newton's method failed.".to_string()),
        };

        // 4. Update X_{n+1} = X_n + delta_x
        for i in 0..x_n.len() {
            x_n[i] += delta_x[i];
        }

        // 5. Check for convergence
        let norm_delta = delta_x.iter().map(|v| v * v).sum::<f64>().sqrt();
        if norm_delta < tolerance {
            return Ok(x_n);
        }
    }

    Err("Newton's method did not converge.".to_string())
}

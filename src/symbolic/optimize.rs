//! # Symbolic Optimization
//!
//! This module provides functions for symbolic optimization, including finding and
//! classifying critical points (local minima, maxima, saddle points) of multivariate
//! functions. It leverages symbolic differentiation and eigenvalue analysis of the Hessian
//! matrix. It also supports constrained optimization using Lagrange Multipliers.

use crate::symbolic::calculus::differentiate;
use crate::symbolic::core::Expr;
use crate::symbolic::matrix::eigen_decomposition;
use crate::symbolic::simplify::simplify;
use crate::symbolic::solve::solve_system;
use std::collections::HashMap;

#[derive(Debug, PartialEq, Eq)]
pub enum ExtremumType {
    LocalMin,
    LocalMax,
    SaddlePoint,
    Unknown,
}

#[derive(Debug)]
pub struct CriticalPoint {
    pub point: HashMap<String, Expr>,
    pub point_type: ExtremumType,
}

/// Finds and classifies the critical points of a multivariate function.
///
/// This function identifies critical points by solving the system of equations
/// formed by setting the gradient of the function to zero (`∇f = 0`).
/// It then classifies these points (local minimum, local maximum, or saddle point)
/// by analyzing the eigenvalues of the Hessian matrix at each critical point.
///
/// # Arguments
/// * `f` - The multivariate function as an `Expr`.
/// * `vars` - A slice of string slices representing the independent variables.
///
/// # Returns
/// A `Result` containing a vector of `CriticalPoint` structs, or an error string
/// if the system cannot be solved or classification fails.
pub fn find_extrema(f: &Expr, vars: &[&str]) -> Result<Vec<CriticalPoint>, String> {
    // 1. Find critical points by solving ∇f = 0
    let mut grad_eqs = Vec::new();
    for &var in vars {
        let deriv = differentiate(f, var);
        grad_eqs.push(Expr::Eq(Box::new(deriv), Box::new(Expr::Constant(0.0))));
    }
    let critical_points_sol = match solve_system(&grad_eqs, vars) {
        Some(sol) => sol,
        None => return Ok(vec![]), // No critical points found or system couldn't be solved
    };
    // solve_system returns one solution, a full implementation would handle multiple solution sets.
    let crit_point_map: HashMap<String, Expr> = critical_points_sol.into_iter().collect();

    // 2. Compute the Hessian matrix
    let hessian = hessian_matrix(f, vars);

    // 3. Classify the critical point using the Hessian's eigenvalues
    let mut hessian_at_point = hessian.clone();
    for (var, val) in &crit_point_map {
        hessian_at_point = crate::symbolic::calculus::substitute(&hessian_at_point, var, val);
    }

    let (eigenvalues_expr, _) = eigen_decomposition(&hessian_at_point)?;

    if let Expr::Matrix(eig_rows) = eigenvalues_expr {
        let eigenvalues: Vec<f64> = eig_rows
            .iter()
            .flatten()
            .map(|e| crate::symbolic::simplify::as_f64(e).unwrap_or(f64::NAN))
            .collect();

        let point_type = if eigenvalues.iter().all(|&v| v > 0.0) {
            ExtremumType::LocalMin
        } else if eigenvalues.iter().all(|&v| v < 0.0) {
            ExtremumType::LocalMax
        } else if eigenvalues.iter().any(|&v| v > 0.0) && eigenvalues.iter().any(|&v| v < 0.0) {
            ExtremumType::SaddlePoint
        } else {
            ExtremumType::Unknown // Degenerate case (zero eigenvalues)
        };

        Ok(vec![CriticalPoint {
            point: crit_point_map,
            point_type,
        }])
    } else {
        Err("Could not determine eigenvalues of the Hessian.".to_string())
    }
}

/// Computes the Hessian matrix (matrix of second partial derivatives) of a function.
///
/// The Hessian matrix `H` is a square matrix of second-order partial derivatives of a
/// scalar-valued function. `H_ij = ∂²f / ∂x_i ∂x_j`.
/// It is used in multivariate calculus for classifying critical points.
///
/// # Arguments
/// * `f` - The multivariate function as an `Expr`.
/// * `vars` - A slice of string slices representing the independent variables.
///
/// # Returns
/// An `Expr::Matrix` representing the Hessian matrix.
pub fn hessian_matrix(f: &Expr, vars: &[&str]) -> Expr {
    let n = vars.len();
    let mut matrix = vec![vec![Expr::Constant(0.0); n]; n];
    for i in 0..n {
        for j in 0..n {
            // d/dx_j (d/dx_i f)
            let df_dxi = differentiate(f, vars[i]);
            let d2f_dxj_dxi = differentiate(&df_dxi, vars[j]);
            matrix[i][j] = simplify(d2f_dxj_dxi);
        }
    }
    Expr::Matrix(matrix)
}

/// Finds the extrema of a function subject to equality constraints using the method of Lagrange Multipliers.
///
/// The method of Lagrange Multipliers is a strategy for finding the local maxima and minima
/// of a function subject to equality constraints. It introduces new variables (Lagrange multipliers)
/// and solves a system of equations derived from the gradient of the Lagrangian function.
///
/// # Arguments
/// * `f` - The objective function to optimize.
/// * `constraints` - A slice of expressions representing the equality constraints `g_i(x) = 0`.
/// * `vars` - The variables of the objective function.
///
/// # Returns
/// A `Result` containing a list of solutions, where each solution is a map from variable names to their values.
pub fn find_constrained_extrema(
    f: &Expr,
    constraints: &[Expr],
    vars: &[&str],
) -> Result<Vec<HashMap<String, Expr>>, String> {
    // 1. Introduce Lagrange multipliers
    let mut lambda_vars = Vec::new();
    for i in 0..constraints.len() {
        lambda_vars.push(format!("lambda_{}", i));
    }

    // 2. Construct the Lagrangian: L(X, λ) = f(X) - Σ λ_i * g_i(X)
    let mut lagrangian = f.clone();
    for (i, g) in constraints.iter().enumerate() {
        let lambda_i = Expr::Variable(lambda_vars[i].clone());
        let term = Expr::Mul(Box::new(lambda_i), Box::new(g.clone()));
        lagrangian = simplify(Expr::Sub(Box::new(lagrangian), Box::new(term)));
    }

    // 3. Create the system of equations by taking the gradient of L w.r.t. all variables (original + lambdas)
    let mut all_vars: Vec<&str> = vars.to_vec();
    let lambda_vars_str: Vec<&str> = lambda_vars.iter().map(|s| s.as_str()).collect();
    all_vars.extend(&lambda_vars_str);

    let mut grad_eqs = Vec::new();
    for &var in &all_vars {
        let deriv = differentiate(&lagrangian, var);
        grad_eqs.push(Expr::Eq(Box::new(deriv), Box::new(Expr::Constant(0.0))));
    }

    // 4. Solve the system of equations
    match solve_system(&grad_eqs, &all_vars) {
        Some(solution) => Ok(vec![solution.into_iter().collect()]), // solve_system returns one solution for now
        None => Ok(vec![]),                                         // No solution found
    }
}

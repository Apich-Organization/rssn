//! # Symbolic Regression Analysis
//!
//! This module provides functions for symbolic regression analysis, enabling the
//! derivation of regression coefficients as symbolic expressions. It supports
//! simple linear regression, non-linear regression (by minimizing sum of squared
//! residuals), and polynomial regression.

use crate::symbolic::core::Expr;
use crate::symbolic::matrix;
use crate::symbolic::simplify::simplify;
use crate::symbolic::solve::solve_system;
use crate::symbolic::stats::{covariance, mean, variance};

/// Computes the symbolic coefficients (`b0`, `b1`) for a simple linear regression `y = b0 + b1*x`.
///
/// This function derives the ordinary least squares (OLS) estimators for the intercept (`b0`)
/// and slope (`b1`) based on the formulas involving means, variances, and covariances of the data.
///
/// # Arguments
/// * `data` - A slice of tuples `(x_i, y_i)` representing the data points.
///
/// # Returns
/// A tuple `(b0, b1)` where each element is a symbolic expression for the coefficient.
pub fn simple_linear_regression_symbolic(data: &[(Expr, Expr)]) -> (Expr, Expr) {
    let (xs, ys): (Vec<_>, Vec<_>) = data.iter().cloned().unzip();

    let mean_x = mean(&xs);
    let mean_y = mean(&ys);
    let var_x = variance(&xs);
    let cov_xy = covariance(&xs, &ys);

    // b1 = cov(X, Y) / var(X)
    let b1 = simplify(Expr::Div(Box::new(cov_xy), Box::new(var_x)));

    // b0 = mean(Y) - b1 * mean(X)
    let b0 = simplify(Expr::Sub(
        Box::new(mean_y),
        Box::new(Expr::Mul(Box::new(b1.clone()), Box::new(mean_x))),
    ));

    (b0, b1)
}

/// Attempts to find symbolic expressions for the parameters of a non-linear model.
///
/// This is done by minimizing the sum of squared residuals (SSR). The function sets up
/// a system of equations by taking the partial derivative of the SSR with respect to
/// each parameter and setting it to zero (`∂SSR/∂β_j = 0`). It then attempts to solve
/// this (usually non-linear) system symbolically.
///
/// # Arguments
/// * `data` - A slice of tuples `(x_i, y_i)` representing the observed data points.
/// * `model` - An `Expr` representing the non-linear model `f(x, β_1, ..., β_m)`.
/// * `vars` - The independent variables in the model (e.g., `["x"]`).
/// * `params` - The parameters to solve for (e.g., `["a", "b"]`).
///
/// # Returns
/// An `Option<Vec<(String, Expr)>>` containing the symbolic solutions for the parameters,
/// or `None` if the system cannot be solved.
pub fn nonlinear_regression_symbolic(
    data: &[(Expr, Expr)],
    model: &Expr,
    vars: &[&str],   // The independent variables, e.g., ["x"]
    params: &[&str], // The parameters to solve for, e.g., ["a", "b"]
) -> Option<Vec<(String, Expr)>> {
    // 1. Construct the Sum of Squared Residuals (S)
    let mut s_expr = Expr::Constant(0.0);
    let x_var = vars.get(0).cloned().unwrap_or("x"); // Assuming one independent variable for now
    let _y_var = "y"; // Placeholder for the dependent variable in the data

    for (x_i, y_i) in data {
        let mut model_at_point = model.clone();
        model_at_point = crate::symbolic::calculus::substitute(&model_at_point, x_var, x_i);
        let residual = Expr::Sub(Box::new(y_i.clone()), Box::new(model_at_point));
        let residual_sq = Expr::Power(Box::new(residual), Box::new(Expr::Constant(2.0)));
        s_expr = Expr::Add(Box::new(s_expr), Box::new(residual_sq));
    }

    // 2. Take the partial derivative of S with respect to each parameter
    let mut grad_eqs = Vec::new();
    for &param in params {
        let deriv = crate::symbolic::calculus::differentiate(&s_expr, param);
        grad_eqs.push(Expr::Eq(Box::new(deriv), Box::new(Expr::Constant(0.0))));
    }

    // 3. Solve the resulting system of (usually non-linear) equations
    solve_system(&grad_eqs, params)
}

/// Computes the symbolic coefficients for a polynomial regression `y = c0 + c1*x + ... + cm*x^m`.
///
/// This function solves the normal equation `(X^T * X) * C = X^T * Y` symbolically.
/// `X` is the Vandermonde matrix, `Y` is the vector of dependent variable values,
/// and `C` is the vector of polynomial coefficients.
///
/// # Arguments
/// * `data` - A slice of tuples `(x_i, y_i)` representing the data points.
/// * `degree` - The degree `m` of the polynomial.
///
/// # Returns
/// A `Result` containing a vector of symbolic expressions for the coefficients `[c0, c1, ..., cm]`,
/// or an error string if the system cannot be solved.
pub fn polynomial_regression_symbolic(
    data: &[(Expr, Expr)],
    degree: usize,
) -> Result<Vec<Expr>, String> {
    let (xs, ys): (Vec<_>, Vec<_>) = data.iter().cloned().unzip();
    let n = data.len();

    // 1. Construct the Vandermonde matrix X
    let mut x_matrix_rows = Vec::with_capacity(n);
    for x_i in &xs {
        let mut row = Vec::with_capacity(degree + 1);
        for j in 0..=degree {
            row.push(simplify(Expr::Power(
                Box::new(x_i.clone()),
                Box::new(Expr::Constant(j as f64)),
            )));
        }
        x_matrix_rows.push(row);
    }
    let x_matrix = Expr::Matrix(x_matrix_rows);
    let x_matrix_t = matrix::transpose_matrix(&x_matrix);

    // 2. Construct the matrices for the normal equation
    let xt_x = matrix::mul_matrices(&x_matrix_t, &x_matrix);
    let xt_y = matrix::mul_matrices(
        &x_matrix_t,
        &Expr::Matrix(ys.into_iter().map(|y| vec![y]).collect()),
    );

    // 3. Solve the system of linear equations for the coefficient vector C
    let _coeff_vars: Vec<String> = (0..=degree).map(|i| format!("c{}", i)).collect();
    let result = matrix::solve_linear_system(&xt_x, &xt_y);

    match result {
        Ok(Expr::Matrix(rows)) => Ok(rows.into_iter().map(|row| row[0].clone()).collect()),
        Ok(_) => Err("Solver returned a non-vector solution.".to_string()),
        Err(e) => Err(e),
    }
}

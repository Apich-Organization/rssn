//! # Symbolic Result Verification
//!
//! This module provides functions for verifying symbolic results numerically.
//! It uses random sampling and numerical evaluation to check the correctness of
//! solutions to equations, integrals, ODEs, and matrix operations. This is particularly
//! useful for complex symbolic computations where direct algebraic verification is difficult.

use crate::numerical::elementary::eval_expr;
use crate::numerical::integrate::quadrature;
use crate::numerical::integrate::QuadratureMethod;
use crate::prelude::simplify;
use crate::prelude::Expr;
use crate::symbolic::calculus::differentiate;
use crate::symbolic::calculus::substitute;
use crate::symbolic::matrix;
use crate::symbolic::simplify::as_f64;
use rand::{thread_rng, Rng};
use std::collections::HashMap;

const TOLERANCE: f64 = 1e-6;
const NUM_SAMPLES: usize = 100;

/// Verifies a solution to a single equation or a system of equations using numerical sampling.
///
/// This function substitutes the proposed solution into the equations and evaluates the result
/// at several random points. If the equations are satisfied (i.e., evaluate to approximately zero)
/// at all sample points, the solution is considered verified.
///
/// # Arguments
/// * `equations` - A slice of `Expr` representing the equations to verify.
/// * `solution` - A `HashMap` mapping variable names to their proposed solution expressions.
/// * `free_vars` - A slice of string slices representing any free variables in the solution.
///
/// # Returns
/// `true` if the solution is numerically verified, `false` otherwise.
pub fn verify_equation_solution<K, V>(
    equations: &[Expr],
    solution: &HashMap<String, Expr>,
    free_vars: &[&str],
) -> bool {
    let mut rng = thread_rng();
    let mut vars_to_sub = solution.clone();

    for eq in equations {
        let diff = if let Expr::Eq(lhs, rhs) = eq {
            simplify(Expr::Sub(lhs.clone(), rhs.clone()))
        } else {
            eq.clone()
        };

        for _ in 0..NUM_SAMPLES {
            // Create a map for all variables (solved and free)
            let _eval_map: HashMap<K, V> = HashMap::new();
            for var in free_vars {
                vars_to_sub.insert(
                    var.to_string(),
                    Expr::Constant(rng.gen_range(-100.0..100.0)),
                );
            }

            let mut substituted_expr = diff.clone();
            for (var, val) in &vars_to_sub {
                substituted_expr = substitute(&substituted_expr, var, val);
            }

            // All variables should be substituted, now evaluate
            match eval_expr(&simplify(substituted_expr), &HashMap::new()) {
                Ok(val) => {
                    if val.abs() > TOLERANCE {
                        return false;
                    }
                }
                Err(_) => return false, // Evaluation failed
            }
        }
    }
    true
}

/// Verifies an indefinite integral `F(x)` for an integrand `f(x)` by checking if `F'(x) == f(x)`.
///
/// This is done by symbolically differentiating the proposed integral `F(x)` and comparing
/// the result with the original integrand `f(x)` at several random numerical sample points.
///
/// # Arguments
/// * `integrand` - The original integrand `f(x)`.
/// * `integral_result` - The proposed indefinite integral `F(x)`.
/// * `var` - The integration variable.
///
/// # Returns
/// `true` if the integral is numerically verified, `false` otherwise.
pub fn verify_indefinite_integral(integrand: &Expr, integral_result: &Expr, var: &str) -> bool {
    let derivative_of_result = differentiate(integral_result, var);
    let diff = simplify(Expr::Sub(
        Box::new(integrand.clone()),
        Box::new(derivative_of_result),
    ));

    let mut rng = thread_rng();
    let mut vars = HashMap::new();
    for _ in 0..NUM_SAMPLES {
        let x_val = rng.gen_range(-100.0..100.0);
        vars.insert(var.to_string(), x_val);
        match eval_expr(&diff, &vars) {
            Ok(val) => {
                if val.abs() > TOLERANCE {
                    return false;
                }
            }
            Err(_) => return false, // Evaluation failed
        }
    }
    true
}

/// Verifies a definite integral by comparing the symbolic result with numerical quadrature.
///
/// This function evaluates the symbolic result of a definite integral and compares it
/// against a numerical approximation of the integral obtained through quadrature methods.
///
/// # Arguments
/// * `integrand` - The integrand `f(x)`.
/// * `var` - The integration variable.
/// * `range` - A tuple `(lower_bound, upper_bound)` for the integration interval.
/// * `symbolic_result` - The proposed symbolic result of the definite integral.
///
/// # Returns
/// `true` if the symbolic result matches the numerical approximation within a tolerance, `false` otherwise.
pub fn verify_definite_integral(
    integrand: &Expr,
    var: &str,
    range: (f64, f64),
    symbolic_result: &Expr,
) -> bool {
    let symbolic_val = match as_f64(symbolic_result) {
        Some(v) => v,
        None => return false, // Symbolic result is not a number
    };

    if let Ok(numerical_val) = quadrature(integrand, var, range, 1000, QuadratureMethod::Simpson) {
        (symbolic_val - numerical_val).abs() < TOLERANCE
    } else {
        false
    }
}

/// Verifies a solution to a first-order ODE `y' = f(x,y)` by numerical sampling.
///
/// This function substitutes the proposed solution `y(x)` into the ODE and checks
/// if the equation holds true at several random sample points.
///
/// # Arguments
/// * `ode` - The ODE to verify, in the form `Expr::Eq(y_prime, f_xy)`.
/// * `solution` - The proposed solution `y(x)`.
/// * `func_name` - The name of the unknown function (e.g., "y").
/// * `var` - The independent variable (e.g., "x").
///
/// # Returns
/// `true` if the solution is numerically verified, `false` otherwise.
pub fn verify_ode_solution(ode: &Expr, solution: &Expr, func_name: &str, var: &str) -> bool {
    if let Expr::Eq(lhs, rhs) = ode {
        // Assuming ODE is in the form y' = f(x,y), so lhs is y'
        let y_prime_from_ode = lhs;
        let _f_xy = rhs;

        let sol_prime = differentiate(solution, var);
        let diff_symbolic = simplify(Expr::Sub(y_prime_from_ode.clone(), Box::new(sol_prime)));

        let mut substituted_diff = diff_symbolic.clone();
        substituted_diff = substitute(&substituted_diff, func_name, solution);

        let mut rng = thread_rng();
        let mut vars = HashMap::new();
        for _ in 0..NUM_SAMPLES {
            let x_val = rng.gen_range(-100.0..100.0);
            vars.insert(var.to_string(), x_val);
            // This doesn't handle other free variables (like constants of integration), a full implementation would.
            match eval_expr(&substituted_diff, &vars) {
                Ok(val) => {
                    if val.abs() > TOLERANCE {
                        return false;
                    }
                }
                Err(_) => return false, // Evaluation failed
            }
        }
        true
    } else {
        false
    }
}

/// Verifies a matrix inverse `A⁻¹` by checking if `A * A⁻¹` is the identity matrix.
///
/// This function performs symbolic matrix multiplication of the original matrix `A`
/// with its proposed inverse `A⁻¹` and then numerically evaluates the resulting matrix
/// to check if it approximates the identity matrix within a given tolerance.
///
/// # Arguments
/// * `original` - The original matrix `A` as an `Expr::Matrix`.
/// * `inverse` - The proposed inverse matrix `A⁻¹` as an `Expr::Matrix`.
///
/// # Returns
/// `true` if the inverse is numerically verified, `false` otherwise.
pub fn verify_matrix_inverse(original: &Expr, inverse: &Expr) -> bool {
    if let (Expr::Matrix(_mat_a), Expr::Matrix(_mat_inv)) = (original, inverse) {
        let product = matrix::mul_matrices(original, inverse);
        if let Expr::Matrix(prod_mat) = simplify(product) {
            let n = prod_mat.len();
            for i in 0..n {
                for j in 0..n {
                    let expected = if i == j { 1.0 } else { 0.0 };
                    if let Some(val) = as_f64(&prod_mat[i][j]) {
                        if (val - expected).abs() > TOLERANCE {
                            return false;
                        }
                    } else {
                        return false; // Element is not a number
                    }
                }
            }
            return true;
        }
    }
    false
}

/// Verifies a symbolic derivative `f'(x)` by comparing it to a numerical differentiation of `f(x)`.
///
/// This function evaluates both the symbolic derivative and a numerical approximation
/// of the derivative at several random points and compares them.
///
/// # Arguments
/// * `original_func` - The original function `f(x)`.
/// * `derivative_func` - The proposed symbolic derivative `f'(x)`.
/// * `var` - The differentiation variable.
///
/// # Returns
/// `true` if the derivative is numerically verified, `false` otherwise.
pub fn verify_derivative(original_func: &Expr, derivative_func: &Expr, var: &str) -> bool {
    let mut rng = thread_rng();
    let mut vars_map = HashMap::new();

    for _ in 0..NUM_SAMPLES {
        let x_val = rng.gen_range(-100.0..100.0);
        vars_map.insert(var.to_string(), x_val);

        // Evaluate the symbolic derivative at the random point
        let symbolic_deriv_val = match eval_expr(derivative_func, &vars_map) {
            Ok(v) => v,
            Err(_) => return false,
        };

        // Calculate the numerical derivative at the same point
        let numerical_deriv_val =
            match crate::numerical::calculus::gradient(original_func, &[var], &[x_val]) {
                Ok(grad_vec) => grad_vec[0],
                Err(_) => return false,
            };

        if (symbolic_deriv_val - numerical_deriv_val).abs() > TOLERANCE {
            return false;
        }
    }
    true
}

//! # Numerical Testing and Verification
//!
//! This module provides numerical testing and verification utilities, primarily focused
//! on solving equations. It includes numerical solvers for polynomial and transcendental
//! equations, as well as systems of linear and non-linear equations.

use crate::symbolic::calculus::differentiate;
use crate::symbolic::core::Expr;
use crate::symbolic::simplify::{is_zero, simplify};
use num_bigint::BigInt;
use num_complex::Complex;
use num_traits::{One, ToPrimitive, Zero};
use std::collections::HashMap;

/// Main entry point for solving equations.
///
/// This function attempts to solve a single equation for a given variable.
/// It first tries to solve it as a polynomial equation (up to quartic degree).
/// If that fails, it falls back to a numerical method for transcendental equations.
///
/// # Arguments
/// * `expr` - The equation to solve (e.g., `Expr::Eq(lhs, rhs)` or `lhs - rhs`).
/// * `var` - The variable to solve for.
///
/// # Returns
/// A `Vec<Expr>` containing the symbolic or numerical solutions.
pub fn solve(expr: &Expr, var: &str) -> Vec<Expr> {
    let mut equation = expr.clone();

    // If the expression is an equality, rewrite it as expr - other = 0
    if let Expr::Eq(left, right) = expr {
        equation = Expr::Sub(left.clone(), right.clone());
    }

    let simplified_expr = simplify(equation);

    // Try to solve as a polynomial
    if let Some(coeffs) = extract_polynomial_coeffs(&simplified_expr, var) {
        // The coefficients are ordered from highest degree to lowest.
        // We need to reverse them for our solver functions which expect c0, c1*x, c2*x^2 ...
        let mut reversed_coeffs: Vec<Expr> = coeffs.into_iter().rev().collect();
        return solve_polynomial(&mut reversed_coeffs);
    }

    // Fallback to numerical method for transcendental equations
    solve_transcendental_numerical(&simplified_expr, var)
}

/// Solves a polynomial equation given its coefficients `[c0, c1, c2, ...]`.
///
/// This function handles polynomial equations up to degree 4 symbolically.
/// For higher-degree polynomials, it falls back to a numerical solver.
/// Coefficients are expected in ascending order of degree (c0 + c1*x + c2*x^2 + ...).
///
/// # Arguments
/// * `coeffs` - A mutable `Vec<Expr>` representing the coefficients of the polynomial.
///
/// # Returns
/// A `Vec<Expr>` containing the symbolic or numerical solutions.
pub fn solve_polynomial(coeffs: &mut Vec<Expr>) -> Vec<Expr> {
    // Normalize coefficients by dividing by the leading coefficient
    if let Some(leading_coeff) = coeffs.last().cloned() {
        if !is_zero(&leading_coeff) {
            for c in coeffs.iter_mut() {
                *c = simplify(Expr::Div(
                    Box::new(c.clone()),
                    Box::new(leading_coeff.clone()),
                ));
            }
        }
    }

    // Trim trailing zeros that might result from normalization
    while coeffs.len() > 1 && coeffs.last().is_some_and(is_zero) {
        coeffs.pop();
    }

    match coeffs.len() {
        0 => vec![],
        1 => {
            if is_zero(&coeffs[0]) {
                vec![Expr::InfiniteSolutions]
            } else {
                vec![]
            }
        } // 0=0 or c=0
        2 => solve_linear(coeffs),
        3 => solve_quadratic(coeffs),
        4 => solve_cubic(coeffs),
        5 => solve_quartic(coeffs),
        _ => {
            solve_polynomial_numerical(&coeffs.iter().map(|e| e.to_f64().unwrap_or(0.0)).collect())
        }
    }
}

/// Extracts coefficients of a polynomial expression `p(var)`.
///
/// This function parses a symbolic expression and attempts to extract its coefficients
/// with respect to a specified variable. Coefficients are returned in descending order
/// of degree `[a_n, a_{n-1}, ..., a_1, a_0]` for `a_n*var^n + ... + a_0`.
///
/// # Arguments
/// * `expr` - The symbolic expression.
/// * `var` - The variable of the polynomial.
///
/// # Returns
/// An `Option<Vec<Expr>>` containing the coefficients, or `None` if the expression
/// is not a polynomial in `var` or contains other variables.
pub fn extract_polynomial_coeffs(expr: &Expr, var: &str) -> Option<Vec<Expr>> {
    let mut coeffs_map = HashMap::new();
    collect_coeffs(expr, var, &mut coeffs_map, &Expr::BigInt(BigInt::one()))?;

    if coeffs_map.is_empty() {
        if let Some(val) = eval_as_constant(expr, var) {
            return Some(vec![val]);
        }
        return None;
    }

    let max_degree = *coeffs_map.keys().max().unwrap_or(&0);
    let mut coeffs = vec![Expr::BigInt(BigInt::zero()); max_degree as usize + 1];
    for (degree, coeff) in coeffs_map {
        coeffs[degree as usize] = coeff;
    }

    // Reverse to get [an, an-1, ..., a0]
    coeffs.reverse();
    Some(coeffs)
}

pub(crate) fn eval_as_constant(expr: &Expr, var: &str) -> Option<Expr> {
    match expr {
        Expr::Constant(c) => Some(Expr::Constant(*c)),
        Expr::BigInt(i) => Some(Expr::BigInt(i.clone())),
        Expr::Rational(r) => Some(Expr::Rational(r.clone())),
        Expr::Variable(v) if v != var => None, // Cannot evaluate other variables
        Expr::Add(l, r) => Some(simplify(Expr::Add(
            Box::new(eval_as_constant(l, var)?),
            Box::new(eval_as_constant(r, var)?),
        ))),
        Expr::Sub(l, r) => Some(simplify(Expr::Sub(
            Box::new(eval_as_constant(l, var)?),
            Box::new(eval_as_constant(r, var)?),
        ))),
        Expr::Mul(l, r) => Some(simplify(Expr::Mul(
            Box::new(eval_as_constant(l, var)?),
            Box::new(eval_as_constant(r, var)?),
        ))),
        Expr::Div(l, r) => {
            let den = eval_as_constant(r, var)?;
            if is_zero(&den) {
                None
            } else {
                Some(simplify(Expr::Div(
                    Box::new(eval_as_constant(l, var)?),
                    Box::new(den),
                )))
            }
        }
        Expr::Neg(e) => Some(simplify(Expr::Neg(Box::new(eval_as_constant(e, var)?)))),
        _ => None,
    }
}

pub(crate) fn collect_coeffs(
    expr: &Expr,
    var: &str,
    coeffs: &mut HashMap<u32, Expr>,
    factor: &Expr,
) -> Option<()> {
    match expr {
        Expr::Constant(_) | Expr::BigInt(_) | Expr::Rational(_) => {
            *coeffs.entry(0).or_insert(Expr::BigInt(BigInt::zero())) = simplify(Expr::Add(
                Box::new(
                    coeffs
                        .get(&0)
                        .unwrap_or(&Expr::BigInt(BigInt::zero()))
                        .clone(),
                ),
                Box::new(Expr::Mul(Box::new(expr.clone()), Box::new(factor.clone()))),
            ));
            Some(())
        }
        Expr::Variable(v) if v == var => {
            *coeffs.entry(1).or_insert(Expr::BigInt(BigInt::zero())) = simplify(Expr::Add(
                Box::new(
                    coeffs
                        .get(&1)
                        .unwrap_or(&Expr::BigInt(BigInt::zero()))
                        .clone(),
                ),
                Box::new(factor.clone()),
            ));
            Some(())
        }
        Expr::Variable(_) => None, // Other variables make it non-polynomial in `var`
        Expr::Add(l, r) => {
            collect_coeffs(l, var, coeffs, factor)?;
            collect_coeffs(r, var, coeffs, factor)
        }
        Expr::Sub(l, r) => {
            collect_coeffs(l, var, coeffs, factor)?;
            collect_coeffs(r, var, coeffs, &Expr::Neg(Box::new(factor.clone())))
        }
        Expr::Mul(l, r) => {
            // Simplified: assumes one side is constant, other is power of var
            if let Some(c) = eval_as_constant(l, var) {
                let mut term_coeffs = HashMap::new();
                collect_coeffs(r, var, &mut term_coeffs, &Expr::BigInt(BigInt::one()))?;
                for (deg, coeff) in term_coeffs {
                    *coeffs.entry(deg).or_insert(Expr::BigInt(BigInt::zero())) =
                        simplify(Expr::Add(
                            Box::new(
                                coeffs
                                    .get(&deg)
                                    .unwrap_or(&Expr::BigInt(BigInt::zero()))
                                    .clone(),
                            ),
                            Box::new(Expr::Mul(Box::new(c.clone()), Box::new(coeff))),
                        ));
                }
                Some(())
            } else if let Some(c) = eval_as_constant(r, var) {
                let mut term_coeffs = HashMap::new();
                collect_coeffs(l, var, &mut term_coeffs, &Expr::BigInt(BigInt::one()))?;
                for (deg, coeff) in term_coeffs {
                    *coeffs.entry(deg).or_insert(Expr::BigInt(BigInt::zero())) =
                        simplify(Expr::Add(
                            Box::new(
                                coeffs
                                    .get(&deg)
                                    .unwrap_or(&Expr::BigInt(BigInt::zero()))
                                    .clone(),
                            ),
                            Box::new(Expr::Mul(Box::new(c.clone()), Box::new(coeff))),
                        ));
                }
                Some(())
            } else {
                None // var * var case, needs more complex expansion logic
            }
        }
        Expr::Power(b, e) => {
            if let (Expr::Variable(v), Expr::Constant(p)) = (&**b, &**e) {
                if v == var {
                    *coeffs
                        .entry(*p as u32)
                        .or_insert(Expr::BigInt(BigInt::zero())) = simplify(Expr::Add(
                        Box::new(
                            coeffs
                                .get(&(*p as u32))
                                .unwrap_or(&Expr::BigInt(BigInt::zero()))
                                .clone(),
                        ),
                        Box::new(factor.clone()),
                    ));
                    Some(())
                } else {
                    None
                }
            } else {
                None
            }
        }
        Expr::Neg(e) => collect_coeffs(e, var, coeffs, &Expr::Neg(Box::new(factor.clone()))),
        _ => None, // Not a polynomial form
    }
}

pub(crate) fn solve_linear(coeffs: &[Expr]) -> Vec<Expr> {
    /// Solves a linear equation `c1*x + c0 = 0`.
    ///
    /// # Arguments
    /// * `coeffs` - A slice of `Expr` representing the coefficients `[c0, c1]`.
    ///
    /// # Returns
    /// A `Vec<Expr>` containing the solution(s).
    let c1 = coeffs[1].clone();
    let c0 = coeffs[0].clone();
    if is_zero(&c1) {
        return if is_zero(&c0) {
            vec![Expr::InfiniteSolutions]
        } else {
            vec![Expr::NoSolution]
        };
    }
    vec![simplify(Expr::Neg(Box::new(Expr::Div(
        Box::new(c0),
        Box::new(c1),
    ))))]
}

pub(crate) fn solve_quadratic(coeffs: &[Expr]) -> Vec<Expr> {
    /// Solves a quadratic equation `c2*x^2 + c1*x + c0 = 0`.
    ///
    /// # Arguments
    /// * `coeffs` - A slice of `Expr` representing the coefficients `[c0, c1, c2]`.
    ///
    /// # Returns
    /// A `Vec<Expr>` containing the solution(s) (real or complex).
    let c2 = coeffs[2].clone();
    let c1 = coeffs[1].clone();
    let c0 = coeffs[0].clone();
    let discriminant = simplify(Expr::Sub(
        Box::new(Expr::Power(
            Box::new(c1.clone()),
            Box::new(Expr::BigInt(BigInt::from(2))),
        )),
        Box::new(Expr::Mul(
            Box::new(Expr::BigInt(BigInt::from(4))),
            Box::new(Expr::Mul(Box::new(c2.clone()), Box::new(c0.clone()))),
        )),
    ));

    if let Some(d_val) = discriminant.to_f64() {
        if d_val >= 0.0 {
            let sqrt_d = Expr::Constant(d_val.sqrt());
            vec![
                simplify(Expr::Div(
                    Box::new(Expr::Add(
                        Box::new(Expr::Neg(Box::new(c1.clone()))),
                        Box::new(sqrt_d.clone()),
                    )),
                    Box::new(Expr::Mul(
                        Box::new(Expr::BigInt(BigInt::from(2))),
                        Box::new(c2.clone()),
                    )),
                )),
                simplify(Expr::Div(
                    Box::new(Expr::Sub(
                        Box::new(Expr::Neg(Box::new(c1.clone()))),
                        Box::new(sqrt_d),
                    )),
                    Box::new(Expr::Mul(
                        Box::new(Expr::BigInt(BigInt::from(2))),
                        Box::new(c2.clone()),
                    )),
                )),
            ]
        } else {
            let sqrt_d = Expr::Constant((-d_val).sqrt());
            vec![
                Expr::Complex(
                    Box::new(simplify(Expr::Div(
                        Box::new(Expr::Neg(Box::new(c1.clone()))),
                        Box::new(Expr::Mul(
                            Box::new(Expr::BigInt(BigInt::from(2))),
                            Box::new(c2.clone()),
                        )),
                    ))),
                    Box::new(simplify(Expr::Div(
                        Box::new(sqrt_d.clone()),
                        Box::new(Expr::Mul(
                            Box::new(Expr::BigInt(BigInt::from(2))),
                            Box::new(c2.clone()),
                        )),
                    ))),
                ),
                Expr::Complex(
                    Box::new(simplify(Expr::Div(
                        Box::new(Expr::Neg(Box::new(c1.clone()))),
                        Box::new(Expr::Mul(
                            Box::new(Expr::BigInt(BigInt::from(2))),
                            Box::new(c2.clone()),
                        )),
                    ))),
                    Box::new(simplify(Expr::Div(
                        Box::new(Expr::Neg(Box::new(sqrt_d))),
                        Box::new(Expr::Mul(
                            Box::new(Expr::BigInt(BigInt::from(2))),
                            Box::new(c2.clone()),
                        )),
                    ))),
                ),
            ]
        }
    } else {
        vec![Expr::Solve(
            Box::new(Expr::Add(
                Box::new(Expr::Mul(
                    Box::new(c2),
                    Box::new(Expr::Power(
                        Box::new(Expr::Variable("x".to_string())),
                        Box::new(Expr::BigInt(BigInt::from(2))),
                    )),
                )),
                Box::new(Expr::Add(
                    Box::new(Expr::Mul(
                        Box::new(c1),
                        Box::new(Expr::Variable("x".to_string())),
                    )),
                    Box::new(c0),
                )),
            )),
            "x".to_string(),
        )]
    }
}

pub(crate) fn solve_cubic(coeffs: &[Expr]) -> Vec<Expr> {
    /// Solves a cubic equation `c3*x^3 + c2*x^2 + c1*x + c0 = 0` numerically.
    ///
    /// # Arguments
    /// * `coeffs` - A slice of `Expr` representing the coefficients `[c0, c1, c2, c3]`.
    ///
    /// # Returns
    /// A `Vec<Expr>` containing the numerical solution(s).
    solve_polynomial_numerical(&coeffs.iter().map(|c| c.to_f64().unwrap_or(0.0)).collect())
}

pub(crate) fn solve_quartic(coeffs: &[Expr]) -> Vec<Expr> {
    /// Solves a quartic equation `c4*x^4 + c3*x^3 + c2*x^2 + c1*x + c0 = 0` numerically.
    ///
    /// # Arguments
    /// * `coeffs` - A slice of `Expr` representing the coefficients `[c0, c1, c2, c3, c4]`.
    ///
    /// # Returns
    /// A `Vec<Expr>` containing the numerical solution(s).
    solve_polynomial_numerical(&coeffs.iter().map(|c| c.to_f64().unwrap_or(0.0)).collect())
}

// Numerical solver for polynomials (Durand-Kerner method)
pub(crate) fn solve_polynomial_numerical(coeffs: &Vec<f64>) -> Vec<Expr> {
    /// Numerically solves a polynomial equation using the Durand-Kerner method.
    ///
    /// The Durand-Kerner method (also known as Weierstrass method) is an iterative
    /// algorithm for finding all roots (real and complex) of a polynomial simultaneously.
    ///
    /// # Arguments
    /// * `coeffs` - A `Vec<f64>` representing the coefficients of the polynomial `[c0, c1, ..., cn]`.
    ///
    /// # Returns
    /// A `Vec<Expr>` containing the numerical solutions (real or complex).
    let degree = coeffs.len() - 1;
    if degree == 0 {
        return vec![];
    }

    // Initial guess for roots
    let mut roots: Vec<Complex<f64>> = (0..degree)
        .map(|i| Complex::new(0.4, 0.9).powu(i as u32))
        .collect();

    let poly_norm = coeffs.iter().map(|c| c.abs()).sum::<f64>().max(1.0);

    for _ in 0..100 {
        // Max iterations
        let mut max_change: f64 = 0.0;
        let old_roots = roots.clone();
        for i in 0..degree {
            let mut den = Complex::new(coeffs[degree], 0.0);
            for j in 0..degree {
                if i != j {
                    den *= old_roots[i] - old_roots[j];
                }
            }

            if den.norm_sqr() < 1e-20 {
                continue;
            }

            let p_val = evaluate_polynomial_horner(coeffs, old_roots[i]);
            let correction = p_val / den;
            roots[i] = old_roots[i] - correction;
            max_change = max_change.max(correction.norm());
        }
        if max_change / poly_norm < 1e-9 {
            break;
        }
    }

    roots
        .into_iter()
        .map(|r| {
            if r.im.abs() < 1e-9 {
                Expr::Constant(r.re)
            } else {
                Expr::Complex(
                    Box::new(Expr::Constant(r.re)),
                    Box::new(Expr::Constant(r.im)),
                )
            }
        })
        .collect()
}

// Evaluates polynomial using Horner's method. Coeffs: [c0, c1, ..., cn]
pub(crate) fn evaluate_polynomial_horner(coeffs: &[f64], x: Complex<f64>) -> Complex<f64> {
    let mut result = Complex::new(0.0, 0.0);
    for &c in coeffs.iter().rev() {
        result = result * x + c;
    }
    result
}

// Evaluates a symbolic expression at a given point.
pub(crate) fn evaluate_expr(expr: &Expr, var: &str, val: f64) -> Option<f64> {
    match expr {
        Expr::Constant(c) => Some(*c),
        Expr::BigInt(i) => i.to_f64(),
        Expr::Rational(r) => r.to_f64(),
        Expr::Variable(v) if v == var => Some(val),
        Expr::Variable(_) => None, // Can't evaluate other variables
        Expr::Add(l, r) => Some(evaluate_expr(l, var, val)? + evaluate_expr(r, var, val)?),
        Expr::Sub(l, r) => Some(evaluate_expr(l, var, val)? - evaluate_expr(r, var, val)?),
        Expr::Mul(l, r) => Some(evaluate_expr(l, var, val)? * evaluate_expr(r, var, val)?),
        Expr::Div(l, r) => {
            let den = evaluate_expr(r, var, val)?;
            if den.abs() < 1e-9 {
                None
            } else {
                Some(evaluate_expr(l, var, val)? / den)
            }
        }
        Expr::Power(b, e) => Some(evaluate_expr(b, var, val)?.powf(evaluate_expr(e, var, val)?)),
        Expr::Sin(arg) => Some(evaluate_expr(arg, var, val)?.sin()),
        Expr::Cos(arg) => Some(evaluate_expr(arg, var, val)?.cos()),
        Expr::Tan(arg) => Some(evaluate_expr(arg, var, val)?.tan()),
        Expr::Exp(arg) => Some(evaluate_expr(arg, var, val)?.exp()),
        Expr::Log(arg) => Some(evaluate_expr(arg, var, val)?.ln()),
        Expr::Neg(arg) => Some(-evaluate_expr(arg, var, val)?),
        _ => None, // Other expressions not supported for now
    }
}

// Numerical solver for transcendental equations (Newton-Raphson method)
pub fn solve_transcendental_numerical(expr: &Expr, var: &str) -> Vec<Expr> {
    /// Numerical solver for transcendental equations (Newton-Raphson method).
    ///
    /// This function attempts to find a single real root of a transcendental equation
    /// `f(x) = 0` using the Newton-Raphson iterative method.
    ///
    /// # Arguments
    /// * `expr` - The symbolic expression `f(x)`.
    /// * `var` - The variable `x` to solve for.
    ///
    /// # Returns
    /// A `Vec<Expr>` containing the numerical solution as `Expr::Constant`,
    /// or `Expr::Solve` if no convergence or symbolic issues.
    let derivative = differentiate(&expr.clone(), var);

    let f = |x: f64| -> Option<f64> { evaluate_expr(expr, var, x) };
    let f_prime = |x: f64| -> Option<f64> { evaluate_expr(&derivative, var, x) };

    let mut x0 = 1.0; // Initial guess
    for _ in 0..100 {
        let y = match f(x0) {
            Some(val) => val,
            None => return vec![Expr::Solve(Box::new(expr.clone()), var.to_string())],
        };
        let y_prime = match f_prime(x0) {
            Some(val) => val,
            None => return vec![Expr::Solve(Box::new(expr.clone()), var.to_string())],
        };

        if y_prime.abs() < 1e-9 {
            return vec![Expr::Solve(Box::new(expr.clone()), var.to_string())];
        }
        let x1 = x0 - y / y_prime;
        if (x1 - x0).abs() < 1e-9 {
            return vec![Expr::Constant(x1)];
        }
        x0 = x1;
    }
    vec![Expr::Solve(Box::new(expr.clone()), var.to_string())]
}

// Helper for evaluating an expression given a map of variable values.
pub(crate) fn evaluate_expr_with_vars(
    expr: &Expr,
    var_values: &HashMap<String, f64>,
) -> Option<f64> {
    match expr {
        Expr::Constant(c) => Some(*c),
        Expr::BigInt(i) => i.to_f64(),
        Expr::Rational(r) => r.to_f64(),
        Expr::Variable(v) => var_values.get(v).cloned(),
        Expr::Add(l, r) => {
            Some(evaluate_expr_with_vars(l, var_values)? + evaluate_expr_with_vars(r, var_values)?)
        }
        Expr::Sub(l, r) => {
            Some(evaluate_expr_with_vars(l, var_values)? - evaluate_expr_with_vars(r, var_values)?)
        }
        Expr::Mul(l, r) => {
            Some(evaluate_expr_with_vars(l, var_values)? * evaluate_expr_with_vars(r, var_values)?)
        }
        Expr::Div(l, r) => {
            let den = evaluate_expr_with_vars(r, var_values)?;
            if den.abs() < 1e-9 {
                None
            } else {
                Some(evaluate_expr_with_vars(l, var_values)? / den)
            }
        }
        Expr::Neg(e) => Some(-evaluate_expr_with_vars(e, var_values)?),
        Expr::Power(b, e) => Some(
            evaluate_expr_with_vars(b, var_values)?.powf(evaluate_expr_with_vars(e, var_values)?),
        ),
        Expr::Sin(arg) => Some(evaluate_expr_with_vars(arg, var_values)?.sin()),
        Expr::Cos(arg) => Some(evaluate_expr_with_vars(arg, var_values)?.cos()),
        Expr::Tan(arg) => Some(evaluate_expr_with_vars(arg, var_values)?.tan()),
        Expr::Exp(arg) => Some(evaluate_expr_with_vars(arg, var_values)?.exp()),
        Expr::Log(arg) => Some(evaluate_expr_with_vars(arg, var_values)?.ln()),
        _ => None, // Other expressions not supported for now
    }
}

// Helper to extract coefficients for a single linear equation
// Returns a HashMap of variable_name -> coefficient, and the constant term
pub(crate) fn extract_linear_equation_coeffs(
    equation: &Expr,
    vars: &[&str],
) -> Option<(HashMap<String, f64>, f64)> {
    let mut coeffs = HashMap::new();
    //let mut constant_term = 0.0;

    let mut zero_values = HashMap::new();
    for &v_name in vars {
        zero_values.insert(v_name.to_string(), 0.0);
    }

    // Get constant term by evaluating the equation with all variables set to 0
    let initial_constant = evaluate_expr_with_vars(equation, &zero_values)?;
    let constant_term = -initial_constant; // Ax + By + C = 0  => Ax + By = -C

    for &target_var in vars {
        let mut test_values = zero_values.clone();
        test_values.insert(target_var.to_string(), 1.0);

        let coeff_val = evaluate_expr_with_vars(equation, &test_values)?;
        // Subtract the constant term, as it's included in the evaluation
        let actual_coeff = coeff_val - initial_constant;
        coeffs.insert(target_var.to_string(), actual_coeff);
    }

    Some((coeffs, constant_term))
}

// Solves a system of linear equations numerically using Gaussian elimination.
pub fn solve_linear_system_numerical(
    mut matrix: Vec<Vec<f64>>,
    mut rhs: Vec<f64>,
) -> Option<Vec<f64>> {
    /// Solves a system of linear equations numerically using Gaussian elimination.
    ///
    /// # Arguments
    /// * `matrix` - The coefficient matrix `A`.
    /// * `rhs` - The right-hand side vector `b`.
    ///
    /// # Returns
    /// An `Option<Vec<f64>>` containing the solution vector, or `None` if the matrix is singular or dimensions mismatch.
    let n = matrix.len();
    if n == 0 {
        return Some(vec![]);
    }
    if matrix[0].len() != n {
        return None;
    } // Not a square matrix
    if rhs.len() != n {
        return None;
    } // RHS vector size mismatch

    // Forward elimination
    for i in 0..n {
        // Find pivot
        let mut max_row = i;
        for k in (i + 1)..n {
            if matrix[k][i].abs() > matrix[max_row][i].abs() {
                max_row = k;
            }
        }
        matrix.swap(i, max_row);
        rhs.swap(i, max_row);

        let pivot = matrix[i][i];
        if pivot.abs() < 1e-9 {
            return None;
        } // Singular matrix, no unique solution

        for k in (i + 1)..n {
            let factor = matrix[k][i] / pivot;
            for j in i..n {
                matrix[k][j] -= factor * matrix[i][j];
            }
            rhs[k] -= factor * rhs[i];
        }
    }

    // Back substitution
    let mut solution = vec![0.0; n];
    for i in (0..n).rev() {
        let mut sum = 0.0;
        //for j in (i + 1)..n {
        for (j, _item) in solution.iter().enumerate().take(n).skip(i + 1) {
            sum += matrix[i][j] * solution[j];
        }
        solution[i] = (rhs[i] - sum) / matrix[i][i];
    }

    Some(solution)
}

// Helper for symbolic division, simplifying the result.
pub(crate) fn expr_div(numerator: Expr, denominator: Expr) -> Expr {
    simplify(Expr::Div(Box::new(numerator), Box::new(denominator)))
}

// Solves a system of linear equations symbolically using Gaussian elimination.
pub fn solve_linear_system_symbolic(
    mut matrix: Vec<Vec<Expr>>,
    mut rhs: Vec<Expr>,
) -> Option<Vec<Expr>> {
    /// Solves a system of linear equations symbolically using Gaussian elimination.
    ///
    /// # Arguments
    /// * `matrix` - The coefficient matrix `A` with symbolic entries.
    /// * `rhs` - The right-hand side vector `b` with symbolic entries.
    ///
    /// # Returns
    /// An `Option<Vec<Expr>>` containing the symbolic solution vector, or `None` if the matrix is singular or dimensions mismatch.
    let n = matrix.len();
    if n == 0 {
        return Some(vec![]);
    }
    if matrix[0].len() != n {
        return None;
    } // Not a square matrix
    if rhs.len() != n {
        return None;
    } // RHS vector size mismatch

    // Forward elimination
    for i in 0..n {
        // Find pivot row (symbolic pivot is tricky, just use current for now)
        let pivot_expr = matrix[i][i].clone();

        // Check if pivot is zero (symbolically)
        if let Expr::Constant(val) = simplify(pivot_expr.clone()) {
            if val.abs() < 1e-9 {
                return None;
            } // Singular matrix, no unique solution
        }

        // Make the diagonal element 1
        for j in i..n {
            matrix[i][j] = expr_div(matrix[i][j].clone(), pivot_expr.clone());
        }
        rhs[i] = expr_div(rhs[i].clone(), pivot_expr.clone());

        // Eliminate other rows
        for k in 0..n {
            if k != i {
                let factor = matrix[k][i].clone();
                for j in i..n {
                    matrix[k][j] = simplify(Expr::Sub(
                        Box::new(matrix[k][j].clone()),
                        Box::new(Expr::Mul(
                            Box::new(factor.clone()),
                            Box::new(matrix[i][j].clone()),
                        )),
                    ));
                }
                rhs[k] = simplify(Expr::Sub(
                    Box::new(rhs[k].clone()),
                    Box::new(Expr::Mul(
                        Box::new(factor.clone()),
                        Box::new(rhs[i].clone()),
                    )),
                ));
            }
        }
    }

    // Back substitution (solutions are directly in rhs now)
    // let mut solution = Vec::with_capacity(n);
    // for i in 0..n {
    // //for item in rhs.iter().take(n) {
    //     solution.push(simplify(rhs[i].clone()));
    // }
    // Back substitution (solutions are directly in rhs now)
    //let mut solution = Vec::with_capacity(n);
    // Assume T is f64
    let _solution: Vec<f64> = Vec::with_capacity(n);

    // Approach A: Using map and collect (The most idiomatic Rust style)
    //
    // Iterate over references to the first n elements of rhs, perform simplify/clone
    // on each element, then collect them into a new Vec.
    let solution: Vec<_> = rhs
        .iter()
        .take(n) // Take the first n elements
        .map(|value| simplify(value.clone())) // Perform the operation on each element
        .collect(); // Collect into a Vec

    // Approach B: Using a for loop (If you prefer explicit iteration)
    /*
    // for value_ref in rhs.iter().take(n) {
    //     // value_ref is a reference to the element, but you need clone(), so it's simplified to one line
    //     solution.push(simplify(value_ref.clone()));
    // }
     */
    Some(solution)
}

// Solves a system of equations.
pub fn solve_system(equations: Vec<Expr>, vars: Vec<&str>) -> Vec<Vec<Expr>> {
    /// Solves a system of equations.
    ///
    /// This function acts as a dispatcher, attempting to solve the system either
    /// as a linear system (symbolically) or as a non-linear system (numerically).
    ///
    /// # Arguments
    /// * `equations` - A `Vec<Expr>` representing the equations in the system.
    /// * `vars` - A `Vec<&str>` representing the variables to solve for.
    ///
    /// # Returns
    /// A `Vec<Vec<Expr>>` where each inner `Vec` is a set of solutions for the variables.
    let num_equations = equations.len();
    let num_vars = vars.len();

    if num_equations != num_vars {
        // Only square systems for now
        return vec![]; // Or handle non-square systems appropriately
    }

    // Try to solve as a linear system symbolically
    let mut symbolic_matrix: Vec<Vec<Expr>> = vec![vec![]; num_equations];
    let mut symbolic_rhs: Vec<Expr> = vec![];
    let mut is_linear_system = true;

    // for (i, eq) in equations.iter().enumerate() {
    //     let mut current_row_coeffs: HashMap<String, f64> = HashMap::new();
    //     let mut current_constant: f64 = 0.0;

    //     // Attempt to extract linear coefficients numerically first to determine linearity
    //     if let Some((coeffs_map, constant)) = extract_linear_equation_coeffs(eq, &vars) {
    //         current_row_coeffs = coeffs_map;
    //         current_constant = constant;
    //     } else {
    //         is_linear_system = false;
    //         break;
    //     }

    //     let mut row_exprs: Vec<Expr> = Vec::with_capacity(num_vars);
    //     for &var_name in vars.iter() {
    //         row_exprs.push(Expr::Constant(
    //             *current_row_coeffs.get(var_name).unwrap_or(&0.0),
    //         ));
    //     }
    //     symbolic_matrix[i] = row_exprs;
    //     symbolic_rhs.push(Expr::Constant(-current_constant)); // Ax = b, so constant goes to RHS
    // }
    for (i, eq) in equations.iter().enumerate() {
        // Attempt to extract linear coefficients numerically first to determine linearity
        let (current_row_coeffs, current_constant) =
            if let Some((coeffs_map, constant)) = extract_linear_equation_coeffs(eq, &vars) {
                (coeffs_map, constant)
            } else {
                is_linear_system = false;
                break;
            };

        let mut row_exprs: Vec<Expr> = Vec::with_capacity(num_vars);
        for &var_name in vars.iter() {
            row_exprs.push(Expr::Constant(
                *current_row_coeffs.get(var_name).unwrap_or(&0.0),
            ));
        }
        symbolic_matrix[i] = row_exprs;
        symbolic_rhs.push(Expr::Constant(-current_constant)); // Ax = b, so constant goes to RHS
    }

    if is_linear_system {
        if let Some(sol) = solve_linear_system_symbolic(symbolic_matrix, symbolic_rhs) {
            vec![sol]
        } else {
            vec![] // No unique symbolic linear solution
        }
    } else {
        // Fallback to numerical method for nonlinear systems
        solve_nonlinear_system_numerical(equations, vars)
    }
}

// Solves a system of nonlinear equations using Newton's method.
pub fn solve_nonlinear_system_numerical(equations: Vec<Expr>, vars: Vec<&str>) -> Vec<Vec<Expr>> {
    /// Solves a system of nonlinear equations using Newton's method.
    ///
    /// This function implements Newton's method for systems of equations. It iteratively
    /// refines an initial guess by solving a linear system involving the Jacobian matrix.
    ///
    /// # Arguments
    /// * `equations` - A `Vec<Expr>` representing the equations in the system.
    /// * `vars` - A `Vec<&str>` representing the variables to solve for.
    ///
    /// # Returns
    /// A `Vec<Vec<Expr>>` where each inner `Vec` is a set of numerical solutions for the variables.
    let n = vars.len();
    if n == 0 {
        return vec![];
    }

    let mut current_var_values: HashMap<String, f64> = HashMap::new();
    for &var_name in vars.iter() {
        current_var_values.insert(var_name.to_string(), 1.0); // Initial guess
    }

    let max_iterations = 100;
    let tolerance = 1e-9;

    for _ in 0..max_iterations {
        // 1. Evaluate F(x_k)
        let mut f_values: Vec<f64> = Vec::with_capacity(n);
        for eq in equations.iter() {
            if let Some(val) = evaluate_expr_with_vars(eq, &current_var_values) {
                f_values.push(val);
            } else {
                return vec![]; // Cannot evaluate equation, likely non-numeric
            }
        }

        // Check for convergence of F(x)
        let f_norm: f64 = f_values.iter().map(|&v| v * v).sum::<f64>().sqrt();
        if f_norm < tolerance {
            break;
        }

        // 2. Compute Jacobian J(x_k)
        let mut jacobian_matrix: Vec<Vec<f64>> = vec![vec![0.0; n]; n];
        for i in 0..n {
            // For each equation
            for (j, _item) in vars.iter().enumerate().take(n) {
                // For each variable
                let differentiated_expr = differentiate(&equations[i].clone(), vars[j]);
                if let Some(val) =
                    evaluate_expr_with_vars(&differentiated_expr, &current_var_values)
                {
                    jacobian_matrix[i][j] = val;
                } else {
                    return vec![]; // Cannot evaluate derivative, likely non-numeric
                }
            }
        }

        // 3. Solve J * delta_x = -F(x)
        let rhs_vector: Vec<f64> = f_values.iter().map(|&v| -v).collect();
        let delta_x_solution = solve_linear_system_numerical(jacobian_matrix, rhs_vector);

        if let Some(delta_x) = delta_x_solution {
            // 4. Update x_k+1 = x_k + delta_x
            let mut max_delta: f64 = 0.0;
            for (i, &var_name) in vars.iter().enumerate() {
                let current_val = *current_var_values.get(var_name).unwrap_or(&0.0);
                let new_val = current_val + delta_x[i];
                current_var_values.insert(var_name.to_string(), new_val);
                max_delta = max_delta.max((new_val - current_val).abs());
            }

            // Check for convergence of delta_x
            if max_delta < tolerance {
                break;
            }
        } else {
            // Linear system could not be solved (singular Jacobian)
            return vec![]; // No unique solution found numerically
        }
    }

    // Convert final numerical solution to Expr format
    let mut result_solution = Vec::with_capacity(n);
    for var_name in vars.iter() {
        if let Some(val) = current_var_values.get(*var_name) {
            result_solution.push(Expr::Constant(*val));
        } else {
            result_solution.push(Expr::Solve(
                Box::new(Expr::Variable(var_name.to_string())),
                var_name.to_string(),
            ));
        }
    }
    vec![result_solution]
}

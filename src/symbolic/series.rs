//! # Symbolic Series Expansions
//!
//! This module provides functions for symbolic series expansions, including Taylor,
//! Laurent, and Fourier series. These tools are fundamental for approximating functions,
//! analyzing their local and global behavior, and solving differential equations.

use crate::symbolic::calculus::{
    definite_integrate, differentiate, evaluate_at_point, factorial, substitute,
};
use crate::symbolic::core::Expr;
use crate::symbolic::simplify::simplify;
use num_bigint::BigInt;
use num_traits::{One, Zero};

/// Computes the Taylor series expansion of an expression around a given center.
///
/// The Taylor series provides a polynomial approximation of a function around a point.
/// It is defined as: `f(x) = Σ_{n=0 to ∞} [f^(n)(a) / n!] * (x - a)^n`.
///
/// # Arguments
/// * `expr` - The expression `f(x)` to expand.
/// * `var` - The variable `x` to expand with respect to.
/// * `center` - The point `a` around which to expand the series.
/// * `order` - The maximum order `N` of the series to compute.
///
/// # Returns
/// An `Expr` representing the truncated Taylor series.
pub fn taylor_series(expr: &Expr, var: &str, center: &Expr, order: usize) -> Expr {
    let coeffs = calculate_taylor_coefficients(expr, var, center, order);
    let mut series_sum = Expr::BigInt(BigInt::zero());

    for (n, coeff) in coeffs.iter().enumerate() {
        let power_term = Expr::Power(
            Box::new(Expr::Sub(
                Box::new(Expr::Variable(var.to_string())),
                Box::new(center.clone()),
            )),
            Box::new(Expr::BigInt(BigInt::from(n))),
        );
        series_sum = simplify(Expr::Add(
            Box::new(series_sum),
            Box::new(Expr::Mul(Box::new(coeff.clone()), Box::new(power_term))),
        ));
    }
    series_sum
}

/// Calculates the coefficients of the Taylor series for a given expression.
/// `c_n = f^(n)(center) / n!`
///
/// # Arguments
/// * `expr` - The expression to expand.
/// * `var` - The variable to expand around.
/// * `center` - The point at which to center the series.
/// * `order` - The order of the series.
///
/// # Returns
/// A vector of `Expr` representing the coefficients `[c_0, c_1, ..., c_order]`.
pub fn calculate_taylor_coefficients(
    expr: &Expr,
    var: &str,
    center: &Expr,
    order: usize,
) -> Vec<Expr> {
    let mut coeffs = Vec::with_capacity(order + 1);
    let mut current_derivative = expr.clone();

    for n in 0..=order {
        let evaluated_derivative = evaluate_at_point(&current_derivative, var, center);
        let n_factorial = factorial(n);
        let term_coefficient = simplify(Expr::Div(
            Box::new(evaluated_derivative),
            Box::new(Expr::Constant(n_factorial)),
        ));
        coeffs.push(term_coefficient);

        if n < order {
            current_derivative = differentiate(&current_derivative, var);
        }
    }
    coeffs
}

/// Computes the Laurent series expansion of an expression around a given center.
///
/// The Laurent series is a generalization of the Taylor series, allowing for terms
/// with negative powers of `(z - c)`. It is particularly useful for analyzing
/// functions with singularities.
///
/// # Arguments
/// * `expr` - The expression `f(z)` to expand.
/// * `var` - The variable `z` to expand with respect to.
/// * `center` - The point `c` around which to expand the series.
/// * `order` - The maximum order of the series (both positive and negative powers).
///
/// # Returns
/// An `Expr` representing the truncated Laurent series.
pub fn laurent_series(expr: &Expr, var: &str, center: &Expr, order: usize) -> Expr {
    let mut k = 0;
    let mut g_z = expr.clone();
    let _help = g_z;
    loop {
        let term = Expr::Power(
            Box::new(Expr::Sub(
                Box::new(Expr::Variable(var.to_string())),
                Box::new(center.clone()),
            )),
            Box::new(Expr::BigInt(BigInt::from(k))),
        );
        let test_expr = simplify(Expr::Mul(Box::new(expr.clone()), Box::new(term)));
        let val_at_center = simplify(evaluate_at_point(&test_expr, var, center));
        if let Expr::Constant(c) = val_at_center {
            if c.is_finite() && c.abs() > 1e-9 {
                g_z = test_expr;
                break;
            }
        }
        k += 1;
        if k > order + 5 {
            return Expr::Series(
                Box::new(expr.clone()),
                var.to_string(),
                Box::new(center.clone()),
                Box::new(Expr::BigInt(BigInt::from(order))),
            );
        }
    }
    let taylor_part = taylor_series(&g_z, var, center, order);
    let divisor = Expr::Power(
        Box::new(Expr::Sub(
            Box::new(Expr::Variable(var.to_string())),
            Box::new(center.clone()),
        )),
        Box::new(Expr::BigInt(BigInt::from(k))),
    );
    simplify(Expr::Div(Box::new(taylor_part), Box::new(divisor)))
}

/// Computes the Fourier series expansion of a periodic expression.
///
/// The Fourier series decomposes a periodic function into a sum of sines and cosines.
/// It is defined as: `f(x) = a_0/2 + Σ_{n=1 to ∞} [a_n cos(nπx/L) + b_n sin(nπx/L)]`.
///
/// # Arguments
/// * `expr` - The periodic expression `f(x)` to expand.
/// * `var` - The variable `x` to expand with respect to.
/// * `period` - The period `T` of the function.
/// * `order` - The maximum order `N` of the series to compute.
///
/// # Returns
/// An `Expr` representing the truncated Fourier series.
pub fn fourier_series(expr: &Expr, var: &str, period: &Expr, order: usize) -> Expr {
    let l = simplify(Expr::Div(
        Box::new(period.clone()),
        Box::new(Expr::BigInt(BigInt::from(2))),
    ));
    let neg_l = simplify(Expr::Neg(Box::new(l.clone())));
    let a0_integrand = expr.clone();
    let a0_integral = definite_integrate(&a0_integrand, var, &neg_l, &l);
    let a0 = simplify(Expr::Div(Box::new(a0_integral), Box::new(l.clone())));
    let mut series_sum = simplify(Expr::Div(
        Box::new(a0),
        Box::new(Expr::BigInt(BigInt::from(2))),
    ));
    for n in 1..=order {
        let n_f64 = n as f64;
        let n_pi_x_over_l = Expr::Div(
            Box::new(Expr::Mul(
                Box::new(Expr::Constant(n_f64 * std::f64::consts::PI)),
                Box::new(Expr::Variable(var.to_string())),
            )),
            Box::new(l.clone()),
        );
        let an_integrand = Expr::Mul(
            Box::new(expr.clone()),
            Box::new(Expr::Cos(Box::new(n_pi_x_over_l.clone()))),
        );
        let an_integral = definite_integrate(&an_integrand, var, &neg_l, &l);
        let an = simplify(Expr::Div(Box::new(an_integral), Box::new(l.clone())));
        let an_term = Expr::Mul(
            Box::new(an),
            Box::new(Expr::Cos(Box::new(n_pi_x_over_l.clone()))),
        );
        series_sum = simplify(Expr::Add(Box::new(series_sum), Box::new(an_term)));
        let bn_integrand = Expr::Mul(
            Box::new(expr.clone()),
            Box::new(Expr::Sin(Box::new(n_pi_x_over_l.clone()))),
        );
        let bn_integral = definite_integrate(&bn_integrand, var, &neg_l, &l);
        let bn = simplify(Expr::Div(Box::new(bn_integral), Box::new(l.clone())));
        let bn_term = Expr::Mul(
            Box::new(bn),
            Box::new(Expr::Sin(Box::new(n_pi_x_over_l.clone()))),
        );
        series_sum = simplify(Expr::Add(Box::new(series_sum), Box::new(bn_term)));
    }
    series_sum
}

/// Computes the symbolic summation of an expression over a given range.
///
/// This function attempts to evaluate finite sums directly. For infinite sums
/// or sums with symbolic bounds, it returns a symbolic `Expr::Summation`.
/// It includes basic rules for arithmetic series and geometric series.
///
/// # Arguments
/// * `expr` - The expression to sum.
/// * `var` - The summation variable.
/// * `lower_bound` - The lower bound of the summation.
/// * `upper_bound` - The upper bound of the summation.
///
/// # Returns
/// An `Expr` representing the sum.
pub fn summation(expr: &Expr, var: &str, lower_bound: &Expr, upper_bound: &Expr) -> Expr {
    if let (Expr::Constant(lower), Expr::Variable(upper_name)) = (lower_bound, upper_bound) {
        // Arithmetic series: sum(a + d*n, n, 0, N) = (N+1)/2 * (2a + d*N)
        if let Expr::Add(a, d_n) = expr {
            if let Expr::Mul(d, n_var) = &**d_n {
                if let Expr::Variable(n_name) = &**n_var {
                    if n_name == var && *lower == 0.0 {
                        let n = Box::new(Expr::Variable(upper_name.clone()));
                        let term1 = Expr::Div(
                            Box::new(Expr::Add(n.clone(), Box::new(Expr::BigInt(BigInt::one())))),
                            Box::new(Expr::BigInt(BigInt::from(2))),
                        );
                        let term2 = Expr::Add(
                            Box::new(Expr::Mul(
                                Box::new(Expr::BigInt(BigInt::from(2))),
                                a.clone(),
                            )),
                            Box::new(Expr::Mul(d.clone(), n)),
                        );
                        return simplify(Expr::Mul(Box::new(term1), Box::new(term2)));
                    }
                }
            }
        }
    }
    if let (Expr::Constant(0.0), Expr::Infinity) = (lower_bound, upper_bound) {
        if let Expr::Power(base, exp) = expr {
            if let Expr::Variable(exp_var_name) = &**exp {
                if exp_var_name == var {
                    return Expr::Div(
                        Box::new(Expr::BigInt(BigInt::one())),
                        Box::new(Expr::Sub(
                            Box::new(Expr::BigInt(BigInt::one())),
                            base.clone(),
                        )),
                    );
                }
            }
        }
    }
    if let (Some(lower_val), Some(upper_val)) = (lower_bound.to_f64(), upper_bound.to_f64()) {
        let mut sum = Expr::BigInt(BigInt::zero());
        for i in lower_val as i64..=upper_val as i64 {
            sum = simplify(Expr::Add(
                Box::new(sum),
                Box::new(evaluate_at_point(expr, var, &Expr::BigInt(BigInt::from(i)))),
            ));
        }
        return sum;
    }
    Expr::Summation(
        Box::new(expr.clone()),
        var.to_string(),
        Box::new(lower_bound.clone()),
        Box::new(upper_bound.clone()),
    )
}

/// Computes the symbolic product of an expression over a given range.
///
/// This function attempts to evaluate finite products directly. For products
/// with symbolic bounds, it returns a symbolic `Expr::Product`.
///
/// # Arguments
/// * `expr` - The expression to multiply.
/// * `var` - The product variable.
/// * `lower_bound` - The lower bound of the product.
/// * `upper_bound` - The upper bound of the product.
///
/// # Returns
/// An `Expr` representing the product.
pub fn product(expr: &Expr, var: &str, lower_bound: &Expr, upper_bound: &Expr) -> Expr {
    if let (Some(lower_val), Some(upper_val)) = (lower_bound.to_f64(), upper_bound.to_f64()) {
        let mut prod = Expr::BigInt(BigInt::one());
        for i in lower_val as i64..=upper_val as i64 {
            prod = simplify(Expr::Mul(
                Box::new(prod),
                Box::new(evaluate_at_point(expr, var, &Expr::BigInt(BigInt::from(i)))),
            ));
        }
        prod
    } else {
        Expr::Product(
            Box::new(expr.clone()),
            var.to_string(),
            Box::new(lower_bound.clone()),
            Box::new(upper_bound.clone()),
        )
    }
}

/// Analyzes the convergence of a series using the Ratio Test.
///
/// The Ratio Test states that for a series `Σ a_n`, if `L = lim (n→∞) |a_{n+1}/a_n|` exists,
/// then the series converges absolutely if `L < 1`, diverges if `L > 1`, and the test is
/// inconclusive if `L = 1`.
///
/// # Arguments
/// * `series_expr` - The series expression, typically `Expr::Summation`.
/// * `var` - The index variable of the series.
///
/// # Returns
/// An `Expr` representing the convergence condition (e.g., `L < 1`).
pub fn analyze_convergence(series_expr: &Expr, var: &str) -> Expr {
    if let Expr::Summation(term, index_var, _, _) = series_expr {
        if index_var == var {
            let an = term;
            let an_plus_1 = evaluate_at_point(
                an,
                var,
                &Expr::Add(
                    Box::new(Expr::Variable(var.to_string())),
                    Box::new(Expr::BigInt(BigInt::one())),
                ),
            );
            let ratio = simplify(Expr::Abs(Box::new(Expr::Div(
                Box::new(an_plus_1),
                Box::new(*an.clone()),
            ))));
            // We need a limit function here. For now, let's assume a simple case.
            // Placeholder for Limit[ratio, var -> Infinity]
            let limit_expr =
                Expr::Limit(Box::new(ratio), var.to_string(), Box::new(Expr::Infinity));
            return Expr::Lt(Box::new(limit_expr), Box::new(Expr::BigInt(BigInt::one())));
        }
    }
    Expr::ConvergenceAnalysis(Box::new(series_expr.clone()), var.to_string())
}

/// Computes the asymptotic expansion of an expression around a given point (e.g., infinity).
///
/// An asymptotic expansion is a series that approximates a function as its argument
/// approaches a particular value (often infinity). It is not necessarily convergent,
/// but provides a good approximation for large arguments.
///
/// # Arguments
/// * `expr` - The expression to expand.
/// * `var` - The variable to expand with respect to.
/// * `point` - The point around which to expand (e.g., `Expr::Infinity`).
/// * `order` - The maximum order of the expansion.
///
/// # Returns
/// An `Expr` representing the asymptotic expansion.
pub fn asymptotic_expansion(expr: &Expr, var: &str, point: &Expr, order: usize) -> Expr {
    // For now, we only implement expansion at infinity for rational functions.
    if !matches!(point, Expr::Infinity) {
        return expr.clone(); // Fallback for other points
    }

    if let Expr::Div(_p, _q) = expr {
        // Substitute x = 1/y
        let y = Expr::Variable("y".to_string());
        let one_over_y = Expr::Div(Box::new(Expr::Constant(1.0)), Box::new(y.clone()));
        let substituted_expr = substitute(expr, var, &one_over_y);

        // The result is a rational function in y. We need to simplify it to a single P(y)/Q(y) form.
        // This simplification step is non-trivial and is assumed to be handled by a robust simplify function.
        let simplified_expr_in_y = simplify(substituted_expr);

        // Compute Taylor series around y = 0
        let taylor_series_in_y =
            taylor_series(&simplified_expr_in_y, "y", &Expr::Constant(0.0), order);

        // Substitute y = 1/x back
        let one_over_x = Expr::Div(
            Box::new(Expr::Constant(1.0)),
            Box::new(Expr::Variable(var.to_string())),
        );
        let final_series = substitute(&taylor_series_in_y, "y", &one_over_x);

        return simplify(final_series);
    }
    let _ = Box::new(expr.clone());
    var.to_string();
    let _ = Box::new(point.clone());
    let _ = Box::new(Expr::BigInt(BigInt::from(order)));

    expr.clone()
}

/// Performs analytic continuation of a function represented by a power series.
///
/// Analytic continuation extends the domain of an analytic function initially
/// defined by a power series in a smaller region. This is achieved by re-expanding
/// the series around a new center within the function's analytic domain.
///
/// # Arguments
/// * `expr` - The original expression (or its power series representation).
/// * `var` - The variable of the function.
/// * `original_center` - The center of the original power series.
/// * `new_center` - The new center for the analytic continuation.
/// * `order` - The order of the new series expansion.
///
/// # Returns
/// An `Expr` representing the analytically continued series.
pub fn analytic_continuation(
    expr: &Expr,
    var: &str,
    original_center: &Expr,
    new_center: &Expr,
    order: usize,
) -> Expr {
    let series_representation = taylor_series(expr, var, original_center, order + 5); // Get a higher order series
    taylor_series(&series_representation, var, new_center, order)
}

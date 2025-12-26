//! # Symbolic Series Expansions
//!
//! This module provides functions for symbolic series expansions, including Taylor,
//! Laurent, and Fourier series. These tools are fundamental for approximating functions,
//! analyzing their local and global behavior, and solving differential equations.
//!
//! ## Examples
//!
//! ### Taylor Series
//! ```
//! 
//! use rssn::symbolic::core::Expr;
//! use rssn::symbolic::series::taylor_series;
//!
//! // Taylor series of e^x around 0 to order 3
//! let x = Expr::new_variable("x");
//!
//! let expr = Expr::new_exp(x);
//!
//! let series = taylor_series(
//!     &expr,
//!     "x",
//!     &Expr::new_constant(0.0),
//!     3,
//! );
//! // Result: 1 + x + x^2/2 + x^3/6
//! ```
//!
//! ### Summation
//! ```
//! 
//! use rssn::symbolic::core::Expr;
//! use rssn::symbolic::series::summation;
//!
//! // Sum of i from 1 to 5
//! let i = Expr::new_variable("i");
//!
//! let sum = summation(
//!     &i,
//!     "i",
//!     &Expr::new_constant(1.0),
//!     &Expr::new_constant(5.0),
//! );
//! // Result: 15
//! ```

use std::sync::Arc;

use num_bigint::BigInt;
use num_traits::One;
use num_traits::Zero;

use crate::symbolic::calculus::definite_integrate;
use crate::symbolic::calculus::differentiate;
use crate::symbolic::calculus::evaluate_at_point;
use crate::symbolic::calculus::factorial;
use crate::symbolic::calculus::substitute;
use crate::symbolic::core::Expr;
use crate::symbolic::simplify_dag::simplify;

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
#[must_use]

pub fn taylor_series(
    expr : &Expr,
    var : &str,
    center : &Expr,
    order : usize,
) -> Expr {

    let coeffs =
        calculate_taylor_coefficients(
            expr,
            var,
            center,
            order,
        );

    let mut series_sum =
        Expr::BigInt(BigInt::zero());

    for (n, coeff) in coeffs
        .iter()
        .enumerate()
    {

        let power_term = Expr::new_pow(
            Expr::new_sub(
                Expr::Variable(
                    var.to_string(),
                ),
                center.clone(),
            ),
            Expr::BigInt(BigInt::from(
                n,
            )),
        );

        series_sum =
            simplify(&Expr::new_add(
                series_sum,
                Expr::new_mul(
                    coeff.clone(),
                    power_term,
                ),
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
#[must_use]

pub fn calculate_taylor_coefficients(
    expr : &Expr,
    var : &str,
    center : &Expr,
    order : usize,
) -> Vec<Expr> {

    let mut coeffs =
        Vec::with_capacity(order + 1);

    let mut current_derivative =
        expr.clone();

    for n in 0 ..= order {

        let evaluated_derivative =
            evaluate_at_point(
                &current_derivative,
                var,
                center,
            );

        let n_factorial = factorial(n);

        let term_coefficient =
            simplify(&Expr::new_div(
                evaluated_derivative,
                Expr::Constant(
                    n_factorial,
                ),
            ));

        coeffs.push(term_coefficient);

        if n < order {

            current_derivative =
                differentiate(
                    &current_derivative,
                    var,
                );
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
#[must_use]

pub fn laurent_series(
    expr : &Expr,
    var : &str,
    center : &Expr,
    order : usize,
) -> Expr {

    let mut k = 0;

    let mut g_z = expr.clone();

    let _help = g_z;

    loop {

        let term = Expr::new_pow(
            Expr::new_sub(
                Expr::Variable(
                    var.to_string(),
                ),
                center.clone(),
            ),
            Expr::BigInt(BigInt::from(
                k,
            )),
        );

        let test_expr =
            simplify(&Expr::new_mul(
                expr.clone(),
                term,
            ));

        let val_at_center = simplify(
            &evaluate_at_point(
                &test_expr,
                var,
                center,
            ),
        );

        if let Expr::Constant(c) =
            val_at_center
        {

            if c.is_finite()
                && c.abs() > 1e-9
            {

                g_z = test_expr;

                break;
            }
        }

        k += 1;

        if k > order + 5 {

            return Expr::Series(
                Arc::new(expr.clone()),
                var.to_string(),
                Arc::new(
                    center.clone(),
                ),
                Arc::new(Expr::BigInt(
                    BigInt::from(order),
                )),
            );
        }
    }

    let taylor_part = taylor_series(
        &g_z,
        var,
        center,
        order,
    );

    let divisor = Expr::new_pow(
        Expr::new_sub(
            Expr::Variable(
                var.to_string(),
            ),
            center.clone(),
        ),
        Expr::BigInt(BigInt::from(k)),
    );

    simplify(&Expr::new_div(
        taylor_part,
        divisor,
    ))
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
#[must_use]

pub fn fourier_series(
    expr : &Expr,
    var : &str,
    period : &Expr,
    order : usize,
) -> Expr {

    let l = simplify(&Expr::new_div(
        period.clone(),
        Expr::BigInt(BigInt::from(2)),
    ));

    let neg_l = simplify(
        &Expr::new_neg(l.clone()),
    );

    let a0_integrand = expr.clone();

    let a0_integral =
        definite_integrate(
            &a0_integrand,
            var,
            &neg_l,
            &l,
        );

    let a0 = simplify(&Expr::new_div(
        a0_integral,
        l.clone(),
    ));

    let mut series_sum =
        simplify(&Expr::new_div(
            a0,
            Expr::BigInt(BigInt::from(
                2,
            )),
        ));

    for n in 1 ..= order {

        let n_f64 = n as f64;

        let n_pi_x_over_l = Expr::new_div(
            Expr::new_mul(
                Expr::Constant(n_f64 * std::f64::consts::PI),
                Expr::Variable(var.to_string()),
            ),
            l.clone(),
        );

        let an_integrand =
            Expr::new_mul(
                expr.clone(),
                Expr::new_cos(
                    n_pi_x_over_l
                        .clone(),
                ),
            );

        let an_integral =
            definite_integrate(
                &an_integrand,
                var,
                &neg_l,
                &l,
            );

        let an =
            simplify(&Expr::new_div(
                an_integral,
                l.clone(),
            ));

        let an_term = Expr::new_mul(
            an,
            Expr::new_cos(
                n_pi_x_over_l.clone(),
            ),
        );

        series_sum =
            simplify(&Expr::new_add(
                series_sum,
                an_term,
            ));

        let bn_integrand =
            Expr::new_mul(
                expr.clone(),
                Expr::new_sin(
                    n_pi_x_over_l
                        .clone(),
                ),
            );

        let bn_integral =
            definite_integrate(
                &bn_integrand,
                var,
                &neg_l,
                &l,
            );

        let bn =
            simplify(&Expr::new_div(
                bn_integral,
                l.clone(),
            ));

        let bn_term = Expr::new_mul(
            bn,
            Expr::new_sin(
                n_pi_x_over_l.clone(),
            ),
        );

        series_sum =
            simplify(&Expr::new_add(
                series_sum,
                bn_term,
            ));
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
#[must_use]

pub fn summation(
    expr : &Expr,
    var : &str,
    lower_bound : &Expr,
    upper_bound : &Expr,
) -> Expr {

    if let (
        Expr::Constant(lower),
        Expr::Variable(upper_name),
    ) = (
        lower_bound,
        upper_bound,
    ) {

        if let Expr::Add(a, d_n) = expr
        {

            if let Expr::Mul(d, n_var) =
                &**d_n
            {

                if let Expr::Variable(
                    n_name,
                ) = &**n_var
                {

                    if n_name == var
                        && *lower == 0.0
                    {

                        let n = Arc::new(Expr::Variable(
                            upper_name.clone(),
                        ));

                        let term1 = Expr::new_div(
                            Expr::new_add(
                                n.clone(),
                                Expr::BigInt(BigInt::one()),
                            ),
                            Expr::BigInt(BigInt::from(2)),
                        );

                        let term2 = Expr::new_add(
                            Expr::new_mul(
                                Expr::BigInt(BigInt::from(2)),
                                a.clone(),
                            ),
                            Expr::new_mul(d.clone(), n),
                        );

                        return simplify(&Expr::new_mul(
                            term1, term2,
                        ));
                    }
                }
            }
        }
    }

    if matches!(
        (
            lower_bound,
            upper_bound
        ),
        (
            Expr::Constant(0.0),
            Expr::Infinity
        )
    ) {

        if let Expr::Power(base, exp) =
            expr
        {

            if let Expr::Variable(
                exp_var_name,
            ) = &**exp
            {

                if exp_var_name == var {

                    return Expr::new_div(
                        Expr::BigInt(BigInt::one()),
                        Expr::new_sub(
                            Expr::BigInt(BigInt::one()),
                            base.clone(),
                        ),
                    );
                }
            }
        }
    }

    if let (
        Some(lower_val),
        Some(upper_val),
    ) = (
        lower_bound.to_f64(),
        upper_bound.to_f64(),
    ) {

        let mut sum = Expr::BigInt(
            BigInt::zero(),
        );

        for i in lower_val as i64
            ..= upper_val as i64
        {

            sum = simplify(&Expr::new_add(
                sum,
                evaluate_at_point(
                    expr,
                    var,
                    &Expr::BigInt(BigInt::from(i)),
                ),
            ));
        }

        return sum;
    }

    Expr::Summation(
        Arc::new(expr.clone()),
        var.to_string(),
        Arc::new(lower_bound.clone()),
        Arc::new(upper_bound.clone()),
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
#[must_use]

pub fn product(
    expr : &Expr,
    var : &str,
    lower_bound : &Expr,
    upper_bound : &Expr,
) -> Expr {

    if let (
        Some(lower_val),
        Some(upper_val),
    ) = (
        lower_bound.to_f64(),
        upper_bound.to_f64(),
    ) {

        let mut prod =
            Expr::BigInt(BigInt::one());

        for i in lower_val as i64
            ..= upper_val as i64
        {

            prod = simplify(&Expr::new_mul(
                prod,
                evaluate_at_point(
                    expr,
                    var,
                    &Expr::BigInt(BigInt::from(i)),
                ),
            ));
        }

        prod
    } else {

        Expr::Product(
            Arc::new(expr.clone()),
            var.to_string(),
            Arc::new(
                lower_bound.clone(),
            ),
            Arc::new(
                upper_bound.clone(),
            ),
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
#[must_use]

pub fn analyze_convergence(
    series_expr : &Expr,
    var : &str,
) -> Expr {

    if let Expr::Summation(
        term,
        index_var,
        _,
        _,
    ) = series_expr
    {

        if index_var == var {

            let an = term;

            let an_plus_1 = evaluate_at_point(
                an,
                var,
                &Expr::new_add(
                    Expr::Variable(var.to_string()),
                    Expr::BigInt(BigInt::one()),
                ),
            );

            let ratio = simplify(
                &Expr::new_abs(
                    Expr::new_div(
                        an_plus_1,
                        an.as_ref()
                            .clone(),
                    ),
                ),
            );

            let limit_expr =
                Expr::Limit(
                    Arc::new(ratio),
                    var.to_string(),
                    Arc::new(
                        Expr::Infinity,
                    ),
                );

            return Expr::Lt(
                Arc::new(limit_expr),
                Arc::new(Expr::BigInt(
                    BigInt::one(),
                )),
            );
        }
    }

    Expr::ConvergenceAnalysis(
        Arc::new(series_expr.clone()),
        var.to_string(),
    )
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
#[must_use]

pub fn asymptotic_expansion(
    expr : &Expr,
    var : &str,
    point : &Expr,
    order : usize,
) -> Expr {

    if !matches!(
        point,
        Expr::Infinity
    ) {

        return expr.clone();
    }

    if let Expr::Div(_p, _q) = expr {

        let y = Expr::Variable(
            "y".to_string(),
        );

        let one_over_y = Expr::new_div(
            Expr::Constant(1.0),
            y,
        );

        let substituted_expr =
            substitute(
                expr,
                var,
                &one_over_y,
            );

        let simplified_expr_in_y =
            simplify(&substituted_expr);

        let taylor_series_in_y =
            taylor_series(
                &simplified_expr_in_y,
                "y",
                &Expr::Constant(0.0),
                order,
            );

        let one_over_x = Expr::new_div(
            Expr::Constant(1.0),
            Expr::Variable(
                var.to_string(),
            ),
        );

        let final_series = substitute(
            &taylor_series_in_y,
            "y",
            &one_over_x,
        );

        return simplify(&final_series);
    }

    let _ = Arc::new(expr.clone());

    var.to_string();

    let _ = Arc::new(point.clone());

    let _ = Arc::new(Expr::BigInt(
        BigInt::from(order),
    ));

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
#[must_use]

pub fn analytic_continuation(
    expr : &Expr,
    var : &str,
    original_center : &Expr,
    new_center : &Expr,
    order : usize,
) -> Expr {

    let series_representation =
        taylor_series(
            expr,
            var,
            original_center,
            order + 5,
        );

    taylor_series(
        &series_representation,
        var,
        new_center,
        order,
    )
}

//! # Convergence Analysis of Series
//!
//! This module provides functions to analyze the convergence of infinite series.
//! It implements several standard convergence tests, including the p-series test,
//! term test, alternating series test, ratio test, root test, and integral test.
use crate::symbolic::calculus::{differentiate, improper_integral, limit, substitute};
use crate::symbolic::core::Expr;
use crate::symbolic::elementary::infinity;
use crate::symbolic::simplify_dag::simplify;
use crate::symbolic::simplify::is_zero;
use num_bigint::BigInt;
use num_traits::One;
/// Represents the result of a convergence test.
#[derive(Debug, PartialEq, Eq)]
pub enum ConvergenceResult {
    /// The series is determined to converge.
    Converges,
    /// The series is determined to diverge.
    Diverges,
    /// The convergence could not be determined with the available tests.
    Inconclusive,
}
/// Checks if a function is eventually positive for large values of the variable.
///
/// This is a heuristic check that evaluates the function at a large number (n=1000)
/// and checks if the result is positive. This is not a formal proof but works for many common functions.
///
/// # Arguments
/// * `f_n` - The expression for the function `f(n)`.
/// * `n` - The variable name.
///
/// # Returns
/// `true` if the function is likely positive for large n, `false` otherwise.
pub(crate) fn is_positive(f_n: &Expr, n: &str) -> bool {
    let large_n = Expr::Constant(1000.0);
    let val_at_large_n = simplify(&substitute(f_n, n, &large_n));
    if let Some(v) = val_at_large_n.to_f64() {
        v > 0.0
    } else {
        false
    }
}
/// Checks if a function is eventually monotonically decreasing for large values of the variable.
///
/// This is a heuristic check that evaluates the derivative of the function at a large number (n=1000)
/// and checks if the result is negative. This indicates that the function is decreasing.
///
/// # Arguments
/// * `f_n` - The expression for the function `f(n)`.
/// * `n` - The variable name.
///
/// # Returns
/// `true` if the function is likely decreasing for large n, `false` otherwise.
pub(crate) fn is_eventually_decreasing(f_n: &Expr, n: &str) -> bool {
    let derivative = differentiate(f_n, n);
    let large_n = Expr::Constant(1000.0);
    let deriv_at_large_n = simplify(&substitute(&derivative, n, &large_n));
    if let Some(v) = deriv_at_large_n.to_f64() {
        v <= 0.0
    } else {
        false
    }
}
/// Analyzes the convergence of a series given its general term `a_n`.
///
/// It applies a sequence of standard convergence tests:
/// 1.  **p-series Test**: For series of the form `1/n^p`.
/// 2.  **Term Test (Test for Divergence)**: If `lim(n->inf) a_n != 0`, the series diverges.
/// 3.  **Alternating Series Test**: For alternating series `sum((-1)^n * b_n)`.
/// 4.  **Ratio Test**: For series `sum(a_n)`, examines `lim(n->inf) |a_{n+1}/a_n|`.
/// 5.  **Root Test**: For series `sum(a_n)`, examines `lim(n->inf) |a_n|^(1/n)`.
/// 6.  **Integral Test**: If `f(x) = a_n` is positive, continuous, and decreasing.
///
/// # Arguments
/// * `a_n` - The general term of the series as an `Expr`.
/// * `n` - The name of the index variable (e.g., "n").
///
/// # Returns
/// A `ConvergenceResult` enum indicating whether the series converges, diverges, or if the test is inconclusive.
pub fn analyze_convergence(a_n: &Expr, n: &str) -> ConvergenceResult {
    if let Expr::Div(one, power) = a_n {
        if let Expr::BigInt(b) = &**one {
            if b.is_one() {
                if let Expr::Power(var, p) = &**power {
                    if let Expr::Variable(name) = &**var {
                        if name == n {
                            if let Some(p_val) = simplify(&p.as_ref().clone()).to_f64() {
                                return if p_val > 1.0 {
                                    ConvergenceResult::Converges
                                } else {
                                    ConvergenceResult::Diverges
                                };
                            }
                        }
                    }
                }
            }
        }
    }
    let term_limit = limit(a_n, n, &infinity());
    let simplified_limit = simplify(&term_limit.clone());
    if !is_zero(&simplified_limit)
        && (simplified_limit.to_f64().is_some() || matches!(simplified_limit, Expr::Infinity)) {
            return ConvergenceResult::Diverges;
        }
    let mut is_alternating = false;
    let mut b_n = a_n.clone();
    if let Expr::Mul(factor1, factor2) = a_n {
        if let Expr::Power(neg_one, _) = &**factor1 {
            if let Expr::BigInt(base) = &**neg_one {
                if base == &BigInt::from(-1) {
                    is_alternating = true;
                    b_n = factor2.as_ref().clone();
                }
            }
        }
    }
    if is_alternating
        && is_eventually_decreasing(&b_n, n) {
            return ConvergenceResult::Converges;
        }
    let n_plus_1 = Expr::new_add(Expr::Variable(n.to_string()), Expr::BigInt(BigInt::one()));
    let a_n_plus_1 = substitute(a_n, n, &n_plus_1);
    let ratio = simplify(&Expr::new_abs(Expr::new_div(a_n_plus_1, a_n.clone())));
    let ratio_limit = limit(&ratio, n, &infinity());
    if let Some(l) = simplify(&ratio_limit).to_f64() {
        if l < 1.0 {
            return ConvergenceResult::Converges;
        }
        if l > 1.0 {
            return ConvergenceResult::Diverges;
        }
    }
    let root_expr = simplify(&Expr::new_pow(
        Expr::new_abs(a_n.clone()),
        Expr::new_div(Expr::BigInt(BigInt::one()), Expr::Variable(n.to_string())),
    ));
    let root_limit = limit(&root_expr, n, &infinity());
    if let Some(l) = simplify(&root_limit).to_f64() {
        if l < 1.0 {
            return ConvergenceResult::Converges;
        }
        if l > 1.0 {
            return ConvergenceResult::Diverges;
        }
    }
    if is_positive(a_n, n) && is_eventually_decreasing(a_n, n) {
        let integral_result = improper_integral(a_n, n);
        if matches!(integral_result, Expr::Infinity) {
            return ConvergenceResult::Diverges;
        }
        if !matches!(integral_result, Expr::Integral { .. }) {
            return ConvergenceResult::Converges;
        }
    }
    ConvergenceResult::Inconclusive
}

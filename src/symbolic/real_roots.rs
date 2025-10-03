//! # Real Root Isolation and Counting
//!
//! This module provides algorithms for finding and isolating real roots of polynomials.
//! It includes an implementation of Sturm's theorem to count the number of distinct
//! real roots in an interval and to generate isolating intervals for these roots.

use crate::symbolic::core::{Expr, SparsePolynomial};
use crate::symbolic::polynomial::{differentiate_poly, gcd};
use crate::symbolic::simplify::as_f64;
use num_traits::ToPrimitive;
use std::collections::HashMap;
use std::ops::Neg;

/// Generates the Sturm sequence for a given polynomial.
///
/// The Sturm sequence is a sequence of polynomials derived from the original polynomial
/// and its derivative. It is used to determine the number of distinct real roots
/// of a polynomial in a given interval.
/// The input polynomial is first made square-free.
///
/// # Arguments
/// * `poly` - The input polynomial as a `SparsePolynomial`.
/// * `var` - The variable of the polynomial.
///
/// # Returns
/// A `Vec<SparsePolynomial>` representing the Sturm sequence.
pub fn sturm_sequence(poly: &SparsePolynomial, var: &str) -> Vec<SparsePolynomial> {
    if poly.terms.is_empty() {
        return vec![];
    }

    // Make the polynomial square-free
    let p_prime = differentiate_poly(poly, var);
    let common_divisor = gcd(poly.clone(), p_prime.clone(), var);
    let p0 = poly.clone().long_division(common_divisor, var).0;

    let mut seq = Vec::new();
    seq.push(p0.clone());

    let p1 = differentiate_poly(&p0, var);
    if p1.terms.is_empty() {
        return seq;
    }
    seq.push(p1);

    let mut i = 1;
    while !seq[i].terms.is_empty() && seq[i].degree(var) > 0 {
        let p_prev = &seq[i - 1];
        let p_curr = &seq[i];
        let (_, remainder) = p_prev.clone().long_division(p_curr.clone(), var);

        if remainder.terms.is_empty() {
            break;
        }
        seq.push(remainder.neg());
        i += 1;
    }

    seq
}

/// Counts the number of sign changes in the Sturm sequence at a given point.
pub(crate) fn count_sign_changes(sequence: &[SparsePolynomial], point: f64, var: &str) -> usize {
    let mut changes = 0;
    let mut last_sign: Option<i8> = None;
    let mut vars = HashMap::new();
    vars.insert(var.to_string(), point);

    for poly in sequence {
        let val = poly.eval(&vars);

        let sign = if val > 1e-9 {
            Some(1)
        } else if val < -1e-9 {
            Some(-1)
        } else {
            None // Zero
        };

        if let Some(s) = sign {
            if let Some(ls) = last_sign {
                if s != ls {
                    changes += 1;
                }
            }
            last_sign = Some(s);
        }
    }
    changes
}

/// Counts the number of distinct real roots of a polynomial in an interval `(a, b]`.
///
/// This function uses Sturm's theorem, which states that the number of distinct real roots
/// of a polynomial in an interval `(a, b]` is equal to the difference in the number of
/// sign changes of the Sturm sequence evaluated at `a` and `b`.
///
/// # Arguments
/// * `poly` - The input polynomial as a `SparsePolynomial`.
/// * `var` - The variable of the polynomial.
/// * `a` - The lower bound of the interval.
/// * `b` - The upper bound of the interval.
///
/// # Returns
/// A `Result` containing the number of distinct real roots as a `usize`,
/// or an error string if evaluation fails.
pub fn count_real_roots_in_interval(
    poly: &SparsePolynomial,
    var: &str,
    a: f64,
    b: f64,
) -> Result<usize, String> {
    let seq = sturm_sequence(poly, var);
    let changes_a = count_sign_changes(&seq, a, var);
    let changes_b = count_sign_changes(&seq, b, var);
    Ok(changes_a.saturating_sub(changes_b))
}

/// Finds isolating intervals for all distinct real roots of a polynomial.
///
/// This function uses a bisection method combined with Sturm's theorem to recursively
/// narrow down intervals until each interval contains exactly one distinct real root.
///
/// # Arguments
/// * `poly` - The input polynomial as a `SparsePolynomial`.
/// * `var` - The variable of the polynomial.
/// * `precision` - The desired maximum width of the isolating intervals.
///
/// # Returns
/// A `Result` containing a `Vec<(f64, f64)>` of tuples, where each tuple `(a, b)`
/// represents an interval `[a, b]` containing exactly one root. Returns an error string
/// if root bounding or counting fails.
pub fn isolate_real_roots(
    poly: &SparsePolynomial,
    var: &str,
    precision: f64,
) -> Result<Vec<(f64, f64)>, String> {
    let sq_free = poly
        .clone()
        .long_division(gcd(poly.clone(), differentiate_poly(poly, var), var), var)
        .0;
    let seq = sturm_sequence(&sq_free, var);

    let bound = root_bound(&sq_free, var)?;
    let mut roots = Vec::new();
    let mut stack = vec![(-bound, bound)];

    while let Some((a, b)) = stack.pop() {
        if b - a < precision {
            continue;
        }

        let changes_a = count_sign_changes(&seq, a, var);
        let changes_b = count_sign_changes(&seq, b, var);
        let num_roots = changes_a.saturating_sub(changes_b);

        if num_roots == 1 {
            // We found an isolating interval. We can shrink it further.
            let mut low = a;
            let mut high = b;
            while high - low > precision {
                let mid = (low + high) / 2.0;
                if count_sign_changes(&seq, low, var) - count_sign_changes(&seq, mid, var) > 0 {
                    high = mid;
                } else {
                    low = mid;
                }
            }
            roots.push((low, high));
        } else if num_roots > 1 {
            let mid = (a + b) / 2.0;
            stack.push((a, mid));
            stack.push((mid, b));
        }
    }

    roots.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    Ok(roots)
}

/// Computes an upper bound for the absolute value of the real roots of a polynomial (Cauchy's bound).
pub(crate) fn root_bound(poly: &SparsePolynomial, var: &str) -> Result<f64, String> {
    let coeffs = poly.get_coeffs_as_vec(var);
    if coeffs.is_empty() {
        return Ok(1.0);
    }

    let leading_coeff_expr = coeffs.first().unwrap();
    let lc = as_f64(leading_coeff_expr).ok_or("Leading coefficient is not numerical.")?;
    if lc == 0.0 {
        return Err("Leading coefficient cannot be zero.".to_string());
    }

    let max_coeff = coeffs
        .iter()
        .skip(1)
        .map(|c| as_f64(c).unwrap_or(0.0).abs())
        .fold(0.0, f64::max);

    Ok(1.0 + max_coeff / lc.abs())
}

// Helper needed for eval
pub fn eval_expr(expr: &Expr, vars: &HashMap<String, f64>) -> f64 {
    match expr {
        Expr::Constant(c) => *c,
        Expr::BigInt(i) => i.to_f64().unwrap_or(0.0),
        Expr::Variable(v) => *vars.get(v).unwrap_or(&0.0),
        Expr::Add(a, b) => eval_expr(a, vars) + eval_expr(b, vars),
        Expr::Sub(a, b) => eval_expr(a, vars) - eval_expr(b, vars),
        Expr::Mul(a, b) => eval_expr(a, vars) * eval_expr(b, vars),
        Expr::Div(a, b) => eval_expr(a, vars) / eval_expr(b, vars),
        Expr::Power(b, e) => eval_expr(b, vars).powf(eval_expr(e, vars)),
        Expr::Neg(a) => -eval_expr(a, vars),
        _ => 0.0, // Fallback for non-numerical expressions
    }
}

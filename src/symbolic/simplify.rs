//! # Symbolic Expression Simplification
//!
//! This module provides algebraic simplification and analysis utilities,
//! powered by the unified E-Graph rewrite engine.

use std::collections::BTreeMap;

use num_bigint::BigInt;
use num_traits::{One, ToPrimitive, Zero};

use crate::symbolic::core::Expr;

pub use crate::symbolic::simplify_dag::{pattern_match, substitute_patterns};

/// The main simplification function.
/// Simplifies the expression using the heuristic E-Graph engine.
#[must_use]
pub fn simplify(expr: Expr) -> Expr {
    crate::symbolic::egraph::simplify(&expr)
}

/// Heuristic simplification function.
/// Delegates directly to the multi-tier heuristic E-Graph engine.
#[must_use]
pub fn heuristic_simplify(expr: Expr) -> Expr {
    crate::symbolic::egraph::simplify(&expr)
}

/// Checks if a symbolic expression is equivalent to zero.
#[inline]
#[must_use]
pub fn is_zero(expr: &Expr) -> bool {
    match expr {
        Expr::Constant(val) if *val == 0.0 => true,
        Expr::BigInt(val) if val.is_zero() => true,
        Expr::Rational(val) if val.is_zero() => true,
        _ => false,
    }
}

/// Checks if a symbolic expression is infinite.
#[inline]
#[must_use]
pub fn is_infinite(expr: &Expr) -> bool {
    match expr {
        Expr::Infinity | Expr::NegativeInfinity => true,
        Expr::Constant(val) => val.is_infinite(),
        _ => false,
    }
}

/// Checks if a symbolic expression is equivalent to the constant value one.
#[inline]
#[must_use]
pub fn is_one(expr: &Expr) -> bool {
    match expr {
        Expr::Constant(val) if (*val - 1.0).abs() < f64::EPSILON => true,
        Expr::BigInt(val) if val.is_one() => true,
        Expr::Rational(val) if val.is_one() => true,
        _ => false,
    }
}

/// Attempts to convert a symbolic expression to its numerical floating-point representation.
#[inline]
#[must_use]
pub fn as_f64(expr: &Expr) -> Option<f64> {
    match expr {
        Expr::Constant(val) => Some(*val),
        Expr::BigInt(val) => val.to_f64(),
        Expr::Rational(val) => val.to_f64(),
        _ => None,
    }
}

/// Checks if an expression is a pure numeric type.
#[must_use]
pub const fn is_numeric(expr: &Expr) -> bool {
    matches!(
        expr,
        Expr::Constant(_) | Expr::BigInt(_) | Expr::Rational(_)
    )
}

/// A `RewriteRule` consists of a pattern and a replacement expression.
#[derive(Debug, Clone)]
pub struct RewriteRule {
    /// The name of the rewrite rule.
    pub name: &'static str,
    /// The pattern to match against.
    pub pattern: Expr,
    /// The replacement expression.
    pub replacement: Expr,
}

/// Returns the name of a rewrite rule.
#[must_use]
pub fn get_name(rule: &RewriteRule) -> String {
    rule.name.to_string()
}

pub(crate) fn complexity(expr: &Expr) -> usize {
    match expr {
        Expr::BigInt(_) => 1,
        Expr::Rational(_) => 2,
        Expr::Constant(_) => 3,
        Expr::Variable(_) | Expr::Pattern(_) => 5,
        Expr::Add(a, b) | Expr::Sub(a, b) | Expr::Mul(a, b) | Expr::Div(a, b) => {
            complexity(a) + complexity(b) + 1
        }
        Expr::Power(a, b) => complexity(a) + complexity(b) + 2,
        Expr::Sin(a) | Expr::Cos(a) | Expr::Tan(a) | Expr::Exp(a) | Expr::Log(a) | Expr::Neg(a) => {
            complexity(a) + 2
        }
        _ => 10,
    }
}

pub(crate) fn fold_constants(expr: Expr) -> Expr {
    let expr = match expr {
        Expr::Add(a, b) => {
            Expr::new_add(
                fold_constants(a.as_ref().clone()),
                fold_constants(b.as_ref().clone()),
            )
        }
        Expr::Sub(a, b) => {
            Expr::new_sub(
                fold_constants(a.as_ref().clone()),
                fold_constants(b.as_ref().clone()),
            )
        }
        Expr::Mul(a, b) => {
            Expr::new_mul(
                fold_constants(a.as_ref().clone()),
                fold_constants(b.as_ref().clone()),
            )
        }
        Expr::Div(a, b) => {
            Expr::new_div(
                fold_constants(a.as_ref().clone()),
                fold_constants(b.as_ref().clone()),
            )
        }
        Expr::Power(base, exp) => {
            Expr::new_pow(
                fold_constants((*base).clone()),
                fold_constants((*exp).clone()),
            )
        }
        Expr::Neg(arg) => Expr::new_neg(fold_constants((*arg).clone())),
        _ => expr,
    };

    match expr {
        Expr::Add(a, b) => {
            if let (Some(va), Some(vb)) = (as_f64(&a), as_f64(&b)) {
                Expr::Constant(va + vb)
            } else {
                Expr::new_add(a, b)
            }
        }
        Expr::Sub(a, b) => {
            if let (Some(va), Some(vb)) = (as_f64(&a), as_f64(&b)) {
                Expr::Constant(va - vb)
            } else {
                Expr::new_sub(a, b)
            }
        }
        Expr::Mul(a, b) => {
            if let (Some(va), Some(vb)) = (as_f64(&a), as_f64(&b)) {
                Expr::Constant(va * vb)
            } else {
                Expr::new_mul(a, b)
            }
        }
        Expr::Div(a, b) => {
            if let (Some(va), Some(vb)) = (as_f64(&a), as_f64(&b)) {
                if vb == 0.0 {
                    Expr::new_div(a, b)
                } else {
                    Expr::Constant(va / vb)
                }
            } else {
                Expr::new_div(a, b)
            }
        }
        Expr::Power(b, e) => {
            if let (Some(vb), Some(ve)) = (as_f64(&b), as_f64(&e)) {
                Expr::Constant(vb.powf(ve))
            } else {
                Expr::new_pow(b, e)
            }
        }
        Expr::Neg(arg) => as_f64(&arg).map_or_else(|| Expr::new_neg(arg), |v| Expr::Constant(-v)),
        _ => expr,
    }
}

pub(crate) fn collect_terms_recursive(
    expr: &Expr,
    coeff: &Expr,
    terms: &mut BTreeMap<Expr, Expr>,
) {
    let mut stack = vec![(expr.clone(), coeff.clone())];

    while let Some((current_expr, current_coeff)) = stack.pop() {
        match &current_expr {
            Expr::Add(a, b) => {
                stack.push((a.as_ref().clone(), current_coeff.clone()));
                stack.push((b.as_ref().clone(), current_coeff));
            }
            Expr::AddList(terms_list) => {
                for term in terms_list {
                    stack.push((term.clone(), current_coeff.clone()));
                }
            }
            Expr::Sub(a, b) => {
                stack.push((a.as_ref().clone(), current_coeff.clone()));
                stack.push((
                    b.as_ref().clone(),
                    fold_constants(Expr::new_neg(current_coeff)),
                ));
            }
            Expr::Mul(a, b) => {
                if is_numeric(a) {
                    stack.push((
                        b.as_ref().clone(),
                        fold_constants(Expr::new_mul(current_coeff, a.as_ref().clone())),
                    ));
                } else if is_numeric(b) {
                    stack.push((
                        a.as_ref().clone(),
                        fold_constants(Expr::new_mul(current_coeff, b.as_ref().clone())),
                    ));
                } else {
                    let base = current_expr;
                    let entry = terms
                        .entry(base)
                        .or_insert_with(|| Expr::BigInt(BigInt::zero()));
                    *entry = fold_constants(Expr::new_add(entry.clone(), current_coeff));
                }
            }
            Expr::MulList(factors) => {
                let mut numeric_part = Expr::BigInt(BigInt::one());
                let mut non_numeric_parts = Vec::new();

                for factor in factors {
                    if is_numeric(factor) {
                        numeric_part = fold_constants(Expr::new_mul(numeric_part, factor.clone()));
                    } else {
                        non_numeric_parts.push(factor.clone());
                    }
                }

                if non_numeric_parts.is_empty() {
                    let base = Expr::BigInt(BigInt::one());
                    let new_coeff = fold_constants(Expr::new_mul(current_coeff, numeric_part));
                    let entry = terms
                        .entry(base)
                        .or_insert_with(|| Expr::BigInt(BigInt::zero()));
                    *entry = fold_constants(Expr::new_add(entry.clone(), new_coeff));
                } else {
                    let base = if non_numeric_parts.len() == 1 {
                        non_numeric_parts[0].clone()
                    } else {
                        Expr::MulList(non_numeric_parts)
                    };
                    let new_coeff = fold_constants(Expr::new_mul(current_coeff, numeric_part));
                    let entry = terms
                        .entry(base)
                        .or_insert_with(|| Expr::BigInt(BigInt::zero()));
                    *entry = fold_constants(Expr::new_add(entry.clone(), new_coeff));
                }
            }
            _ => {
                let base = current_expr;
                let entry = terms
                    .entry(base)
                    .or_insert_with(|| Expr::BigInt(BigInt::zero()));
                *entry = fold_constants(Expr::new_add(entry.clone(), current_coeff));
            }
        }
    }
}

/// Collects terms from an expression and orders them by complexity.
#[must_use]
pub fn collect_and_order_terms(expr: &Expr) -> (Expr, Vec<(Expr, Expr)>) {
    let mut terms = BTreeMap::new();
    collect_terms_recursive(expr, &Expr::BigInt(BigInt::one()), &mut terms);

    let mut sorted_terms: Vec<(Expr, Expr)> = terms.into_iter().collect();
    sorted_terms.sort_by(|(b1, _), (b2, _)| complexity(b2).cmp(&complexity(b1)));

    let constant_term = sorted_terms
        .iter()
        .position(|(b, _)| is_one(b))
        .map_or_else(
            || Expr::BigInt(BigInt::zero()),
            |pos| {
                let (_, c) = sorted_terms.remove(pos);
                c
            },
        );

    (constant_term, sorted_terms)
}

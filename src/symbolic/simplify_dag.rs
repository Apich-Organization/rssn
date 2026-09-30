//! # Symbolic Expression Simplification
//!
//! This module provides the primary entry point for symbolic expression simplification
//! and structural pattern matching, backed by the Heuristic E-Graph engine.

use std::collections::HashMap;
use std::sync::Arc;

use super::core::Expr;

/// The main entry point for symbolic expression simplification.
///
/// Simplifies the input expression into canonical, cost-minimal form
/// using the multi-tier heuristic E-Graph rewrite engine.
///
/// # Arguments
/// * `expr` - The expression to simplify.
///
/// # Returns
/// A new, simplified `Expr`.
#[must_use]
pub fn simplify(expr: &Expr) -> Expr {
    crate::symbolic::egraph::simplify(expr)
}

/// Recursively attempts to match an expression against a pattern.
///
/// Handles `Expr::Dag` nodes by unwrapping them to their underlying structure for matching.
#[must_use]
pub fn pattern_match(
    expr: &Expr,
    pattern: &Expr,
) -> Option<HashMap<String, Expr>> {
    let mut assignments = HashMap::new();
    if pattern_match_recursive(expr, pattern, &mut assignments) {
        Some(assignments)
    } else {
        None
    }
}

pub(crate) fn pattern_match_recursive(
    expr: &Expr,
    pattern: &Expr,
    assignments: &mut HashMap<String, Expr>,
) -> bool {
    let expr_unwrapped = match expr {
        Expr::Dag(node) => node.to_expr().unwrap_or_else(|_| expr.clone()),
        _ => expr.clone(),
    };

    let pattern_unwrapped = match pattern {
        Expr::Dag(node) => node.to_expr().unwrap_or_else(|_| pattern.clone()),
        _ => pattern.clone(),
    };

    match (&expr_unwrapped, &pattern_unwrapped) {
        (_, Expr::Pattern(name)) => {
            if let Some(existing) = assignments.get(name) {
                return existing == expr;
            }
            assignments.insert(name.clone(), expr.clone());
            true
        }
        (Expr::Add(e1, e2), Expr::Add(p1, p2)) | (Expr::Mul(e1, e2), Expr::Mul(p1, p2)) => {
            let original_assignments = assignments.clone();
            if pattern_match_recursive(e1, p1, assignments)
                && pattern_match_recursive(e2, p2, assignments)
            {
                return true;
            }
            *assignments = original_assignments;
            pattern_match_recursive(e1, p2, assignments)
                && pattern_match_recursive(e2, p1, assignments)
        }
        (Expr::Sub(e1, e2), Expr::Sub(p1, p2))
        | (Expr::Div(e1, e2), Expr::Div(p1, p2))
        | (Expr::Power(e1, e2), Expr::Power(p1, p2)) => {
            pattern_match_recursive(e1, p1, assignments)
                && pattern_match_recursive(e2, p2, assignments)
        }
        (Expr::Sin(e), Expr::Sin(p))
        | (Expr::Cos(e), Expr::Cos(p))
        | (Expr::Tan(e), Expr::Tan(p))
        | (Expr::Exp(e), Expr::Exp(p))
        | (Expr::Log(e), Expr::Log(p))
        | (Expr::Neg(e), Expr::Neg(p))
        | (Expr::Abs(e), Expr::Abs(p))
        | (Expr::Sqrt(e), Expr::Sqrt(p)) => pattern_match_recursive(e, p, assignments),

        (Expr::NaryList(s1, v1), Expr::NaryList(s2, v2)) => {
            if s1 != s2 || v1.len() != v2.len() {
                return false;
            }
            let original_assignments = assignments.clone();
            for (e, p) in v1.iter().zip(v2.iter()) {
                if !pattern_match_recursive(e, p, assignments) {
                    *assignments = original_assignments;
                    return false;
                }
            }
            true
        }
        (Expr::UnaryList(s1, e1), Expr::UnaryList(s2, p1)) => {
            s1 == s2 && pattern_match_recursive(e1, p1, assignments)
        }
        (Expr::BinaryList(s1, e1a, e1b), Expr::BinaryList(s2, p1a, p1b)) => {
            s1 == s2
                && pattern_match_recursive(e1a, p1a, assignments)
                && pattern_match_recursive(e1b, p1b, assignments)
        }
        _ => expr_unwrapped == pattern_unwrapped,
    }
}

/// Substitutes pattern variables in a template expression with assigned expressions.
#[must_use]
pub fn substitute_patterns<S: std::hash::BuildHasher>(
    template: &Expr,
    assignments: &HashMap<String, Expr, S>,
) -> Expr {
    let template_unwrapped = match template {
        Expr::Dag(node) => node.to_expr().unwrap_or_else(|_| template.clone()),
        _ => template.clone(),
    };

    match template_unwrapped {
        Expr::Pattern(name) => {
            assignments
                .get(&name)
                .cloned()
                .unwrap_or_else(|| template.clone())
        }
        Expr::Add(a, b) => {
            Expr::new_add(
                substitute_patterns(&a, assignments),
                substitute_patterns(&b, assignments),
            )
        }
        Expr::Sub(a, b) => {
            Expr::new_sub(
                substitute_patterns(&a, assignments),
                substitute_patterns(&b, assignments),
            )
        }
        Expr::Mul(a, b) => {
            Expr::new_mul(
                substitute_patterns(&a, assignments),
                substitute_patterns(&b, assignments),
            )
        }
        Expr::Div(a, b) => {
            Expr::new_div(
                substitute_patterns(&a, assignments),
                substitute_patterns(&b, assignments),
            )
        }
        Expr::Power(b, e) => {
            Expr::new_pow(
                substitute_patterns(&b, assignments),
                substitute_patterns(&e, assignments),
            )
        }
        Expr::Sin(a) => Expr::new_sin(substitute_patterns(&a, assignments)),
        Expr::Cos(a) => Expr::new_cos(substitute_patterns(&a, assignments)),
        Expr::Tan(a) => Expr::new_tan(substitute_patterns(&a, assignments)),
        Expr::Exp(a) => Expr::new_exp(substitute_patterns(&a, assignments)),
        Expr::Log(a) => Expr::new_log(substitute_patterns(&a, assignments)),
        Expr::Neg(a) => Expr::new_neg(substitute_patterns(&a, assignments)),
        Expr::Abs(a) => Expr::new_abs(substitute_patterns(&a, assignments)),
        Expr::Sqrt(a) => Expr::new_sqrt(substitute_patterns(&a, assignments)),

        Expr::NaryList(s, v) => {
            Expr::NaryList(
                s,
                v.iter()
                    .map(|e| substitute_patterns(e, assignments))
                    .collect(),
            )
        }
        Expr::UnaryList(s, e) => {
            Expr::UnaryList(s, Arc::new(substitute_patterns(&e, assignments)))
        }
        Expr::BinaryList(s, a, b) => {
            Expr::BinaryList(
                s,
                Arc::new(substitute_patterns(&a, assignments)),
                Arc::new(substitute_patterns(&b, assignments)),
            )
        }
        _ => template.clone(),
    }
}

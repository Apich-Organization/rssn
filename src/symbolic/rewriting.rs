//! # Term Rewriting Systems
//!
//! This module provides tools for working with term rewriting systems (TRS).
//! It includes structures for `RewriteRule`s, functions to apply these rules
//! to expressions, and an implementation of the Knuth-Bendix completion algorithm
//! to attempt to convert a set of equations into a confluent and Noetherian
//! term rewriting system.

use crate::symbolic::calculus::substitute;
use crate::symbolic::core::Expr;
use crate::symbolic::polynomial::contains_var;
use crate::symbolic::simplify::{pattern_match, substitute_patterns};
use std::collections::HashMap;

/// Represents a rewrite rule, e.g., `lhs -> rhs`.
#[derive(Debug, Clone)]
pub struct RewriteRule {
    pub lhs: Expr,
    pub rhs: Expr,
}

/// A simple term ordering based on the complexity (number of nodes) of an expression.
/// Returns true if e1 > e2.
pub(crate) fn is_greater(e1: &Expr, e2: &Expr) -> bool {
    complexity(e1) > complexity(e2)
}

/// Applies a set of rewrite rules repeatedly to an expression until a fixed point is reached.
///
/// This function computes the "normal form" of the expression with respect to the given
/// set of rewrite rules. It iteratively applies rules until no further changes can be made.
///
/// # Arguments
/// * `expr` - The expression to transform.
/// * `rules` - A slice of `RewriteRule`s to apply.
///
/// # Returns
/// A new `Expr` representing the normal form of the expression.
pub fn apply_rules_to_normal_form(expr: &Expr, rules: &[RewriteRule]) -> Expr {
    let mut current_expr = expr.clone();
    let mut changed = true;

    while changed {
        changed = false;
        let (next_expr, applied) = apply_rules_once(&current_expr, rules);
        if applied {
            current_expr = next_expr;
            changed = true;
        }
    }
    current_expr
}

/// Applies the first applicable rule to the expression tree in a pre-order traversal.
pub(crate) fn apply_rules_once(expr: &Expr, rules: &[RewriteRule]) -> (Expr, bool) {
    // Try to apply a rule at the root first
    for rule in rules {
        if let Some(assignments) = pattern_match(expr, &rule.lhs) {
            return (substitute_patterns(&rule.rhs, &assignments), true);
        }
    }

    // If no rule applies at the root, try children
    match expr {
        Expr::Add(a, b) => {
            let (na, ca) = apply_rules_once(a, rules);
            if ca {
                return (Expr::Add(Box::new(na), b.clone()), true);
            }
            let (nb, cb) = apply_rules_once(b, rules);
            if cb {
                return (Expr::Add(a.clone(), Box::new(nb)), true);
            }
        }
        Expr::Mul(a, b) => {
            let (na, ca) = apply_rules_once(a, rules);
            if ca {
                return (Expr::Mul(Box::new(na), b.clone()), true);
            }
            let (nb, cb) = apply_rules_once(b, rules);
            if cb {
                return (Expr::Mul(a.clone(), Box::new(nb)), true);
            }
        }
        // ... other expression types would follow ...
        _ => {}
    }

    (expr.clone(), false) // No rule applied
}

/// Attempts to produce a complete term-rewriting system from a set of equations
/// using the Knuth-Bendix completion algorithm.
///
/// The Knuth-Bendix algorithm takes a set of equations and tries to convert them
/// into a confluent and Noetherian (terminating) set of rewrite rules. This is done
/// by generating and resolving "critical pairs" (overlaps between rules).
///
/// # Arguments
/// * `equations` - A slice of `Expr::Eq` representing the initial equations.
///
/// # Returns
/// A `Result` containing a `Vec<RewriteRule>` if the completion is successful,
/// or an error string if the input is invalid or the algorithm fails to complete.
pub fn knuth_bendix(equations: &[Expr]) -> Result<Vec<RewriteRule>, String> {
    let mut rules: Vec<RewriteRule> = Vec::new();
    for eq in equations {
        if let Expr::Eq(lhs, rhs) = eq {
            if is_greater(lhs, rhs) {
                rules.push(RewriteRule {
                    lhs: *lhs.clone(),
                    rhs: *rhs.clone(),
                });
            } else if is_greater(rhs, lhs) {
                rules.push(RewriteRule {
                    lhs: *rhs.clone(),
                    rhs: *lhs.clone(),
                });
            }
        } else {
            return Err("Input must be a list of equations (Expr::Eq).".to_string());
        }
    }

    let mut i = 0;
    while i < rules.len() {
        let mut j = 0;
        while j <= i {
            let (rule1, rule2) = (&rules[i].clone(), &rules[j].clone());
            let critical_pairs = find_critical_pairs(rule1, rule2);

            for (t1, t2) in critical_pairs {
                let n1 = apply_rules_to_normal_form(&t1, &rules);
                let n2 = apply_rules_to_normal_form(&t2, &rules);

                if n1 != n2 {
                    // Add new rule
                    let new_rule = if is_greater(&n1, &n2) {
                        RewriteRule { lhs: n1, rhs: n2 }
                    } else {
                        RewriteRule { lhs: n2, rhs: n1 }
                    };
                    // Check if the new rule is not trivial and not already present
                    if new_rule.lhs != new_rule.rhs && !rules.iter().any(|r| r.lhs == new_rule.lhs)
                    {
                        rules.push(new_rule);
                        // Restart the process since the rule set has changed
                        i = 0;
                        j = 0;
                    }
                }
            }
            j += 1;
        }
        i += 1;
    }

    Ok(rules)
}

/// Finds critical pairs between two rewrite rules.
pub(crate) fn find_critical_pairs(r1: &RewriteRule, r2: &RewriteRule) -> Vec<(Expr, Expr)> {
    let mut pairs = Vec::new();
    let mut sub_expressions = Vec::new();
    r1.lhs.pre_order_walk(&mut |sub_expr| {
        sub_expressions.push(sub_expr.clone());
    });

    for sub_expr in &sub_expressions {
        if let Some(subst) = unify(sub_expr, &r2.lhs) {
            // Critical pair 1: Apply rule 2 to the sub-expression of rule 1's LHS
            let t1 = substitute(&r1.lhs, &sub_expr.to_string(), &r2.rhs);
            let t1_subst = substitute_patterns(&t1, &subst);

            // Critical pair 2: Apply rule 1 at the top level
            let t2 = substitute_patterns(&r1.rhs, &subst);

            if t1_subst != t2 {
                pairs.push((t1_subst, t2));
            }
        }
    }
    pairs
}

/// Unifies two expressions, finding a substitution that makes them equal.
/// Returns a map of substitutions if successful.
pub(crate) fn unify(e1: &Expr, e2: &Expr) -> Option<HashMap<String, Expr>> {
    let mut subst = HashMap::new();
    if unify_recursive(e1, e2, &mut subst) {
        Some(subst)
    } else {
        None
    }
}

pub(crate) fn unify_recursive(e1: &Expr, e2: &Expr, subst: &mut HashMap<String, Expr>) -> bool {
    match (e1, e2) {
        (Expr::Pattern(p), _) => {
            if let Some(val) = subst.get(p) {
                return val == e2;
            }
            // Occurs check: ensure the variable is not contained in the expression to avoid infinite loops.
            if contains_var(e2, p) {
                return false;
            }
            subst.insert(p.clone(), e2.clone());
            true
        }
        (_, Expr::Pattern(p)) => {
            if let Some(val) = subst.get(p) {
                return val == e1;
            }
            if contains_var(e1, p) {
                return false;
            }
            subst.insert(p.clone(), e1.clone());
            true
        }
        (Expr::Add(a1, b1), Expr::Add(a2, b2)) | (Expr::Mul(a1, b1), Expr::Mul(a2, b2)) => {
            // Commutative unification
            let original_subst = subst.clone();
            if unify_recursive(a1, a2, subst) && unify_recursive(b1, b2, subst) {
                true
            } else {
                *subst = original_subst;
                unify_recursive(a1, b2, subst) && unify_recursive(b1, a2, subst)
            }
        }
        (Expr::Sub(a1, b1), Expr::Sub(a2, b2))
        | (Expr::Div(a1, b1), Expr::Div(a2, b2))
        | (Expr::Power(a1, b1), Expr::Power(a2, b2)) => {
            unify_recursive(a1, a2, subst) && unify_recursive(b1, b2, subst)
        }
        (Expr::Sin(a1), Expr::Sin(a2))
        | (Expr::Cos(a1), Expr::Cos(a2))
        | (Expr::Tan(a1), Expr::Tan(a2))
        | (Expr::Log(a1), Expr::Log(a2))
        | (Expr::Exp(a1), Expr::Exp(a2))
        | (Expr::Neg(a1), Expr::Neg(a2)) => unify_recursive(a1, a2, subst),
        _ => e1 == e2,
    }
}

/// Calculates a simple complexity measure for an expression.
pub(crate) fn complexity(expr: &Expr) -> usize {
    match expr {
        Expr::Add(a, b) | Expr::Mul(a, b) | Expr::Sub(a, b) | Expr::Div(a, b) => {
            complexity(a) + complexity(b) + 1
        }
        Expr::Power(b, e) => complexity(b) + complexity(e) + 2,
        Expr::Sin(a) | Expr::Cos(a) | Expr::Tan(a) | Expr::Log(a) | Expr::Exp(a) | Expr::Neg(a) => {
            complexity(a) + 1
        }
        _ => 1, // Variables, Constants
    }
}

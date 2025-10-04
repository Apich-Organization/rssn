//! # Symbolic Expression Simplification
//!
//! This module provides functions for symbolic expression simplification.
//! It includes a core `simplify` function that applies deterministic algebraic rules,
//! and a `heuristic_simplify` function that uses pattern matching and rewrite rules
//! to find simpler forms of expressions. It also contains utilities for term collection
//! and rational expression simplification.

use crate::symbolic::calculus::substitute;
use crate::symbolic::core::Expr;
use num_bigint::BigInt;
use num_traits::{One, ToPrimitive, Zero};
use std::collections::{BTreeMap, HashMap};

// =====================================================================================
// region: Helper Functions
// =====================================================================================

pub fn is_zero(expr: &Expr) -> bool {
    // Checks if an expression is numerically equal to zero.
    //
    // Handles `Expr::Constant`, `Expr::BigInt`, and `Expr::Rational` variants.
    matches!(expr, Expr::Constant(val) if *val == 0.0)
        || matches!(expr, Expr::BigInt(val) if val.is_zero())
        || matches!(expr, Expr::Rational(val) if val.is_zero())
}

pub fn is_one(expr: &Expr) -> bool {
    // Checks if an expression is numerically equal to one.
    //
    // Handles `Expr::Constant`, `Expr::BigInt`, and `Expr::Rational` variants.
    matches!(expr, Expr::Constant(val) if *val == 1.0)
        || matches!(expr, Expr::BigInt(val) if val.is_one())
        || matches!(expr, Expr::Rational(val) if val.is_one())
}

pub fn as_f64(expr: &Expr) -> Option<f64> {
    // Attempts to convert an expression to an `f64` value.
    //
    // Handles `Expr::Constant`, `Expr::BigInt`, and `Expr::Rational` variants.
    //
    // # Returns
    // `Some(f64)` if the conversion is successful, `None` otherwise.
    match expr {
        Expr::Constant(val) => Some(*val),
        Expr::BigInt(val) => val.to_f64(),
        Expr::Rational(val) => val.to_f64(),
        _ => None,
    }
}

// endregion

// =====================================================================================
// region: Core Simplification Logic
// =====================================================================================

/// The main simplification function.
/// It recursively simplifies an expression tree by applying deterministic algebraic rules.
///
/// This function performs a deep simplification, traversing the expression tree
/// and applying various algebraic identities and arithmetic evaluations.
/// It also includes a step for simplifying rational expressions by canceling common factors.
///
/// # Arguments
/// * `expr` - The expression to simplify.
///
/// # Returns
/// A new, simplified `Expr`.
pub fn simplify(expr: Expr) -> Expr {
    let simplified_expr = match expr {
        Expr::Add(a, b) => Expr::Add(Box::new(simplify(*a)), Box::new(simplify(*b))),
        Expr::Sub(a, b) => Expr::Sub(Box::new(simplify(*a)), Box::new(simplify(*b))),
        Expr::Mul(a, b) => Expr::Mul(Box::new(simplify(*a)), Box::new(simplify(*b))),
        Expr::Div(a, b) => Expr::Div(Box::new(simplify(*a)), Box::new(simplify(*b))),
        Expr::Power(b, e) => Expr::Power(Box::new(simplify(*b)), Box::new(simplify(*e))),
        Expr::Sin(arg) => Expr::Sin(Box::new(simplify(*arg))),
        Expr::Cos(arg) => Expr::Cos(Box::new(simplify(*arg))),
        Expr::Tan(arg) => Expr::Tan(Box::new(simplify(*arg))),
        Expr::Exp(arg) => Expr::Exp(Box::new(simplify(*arg))),
        Expr::Log(arg) => Expr::Log(Box::new(simplify(*arg))),
        Expr::Neg(arg) => Expr::Neg(Box::new(simplify(*arg))),
        Expr::Sum {
            body,
            var,
            from,
            to,
        } => Expr::Sum {
            body: Box::new(simplify(*body)),
            var: Box::new(simplify(*var)),
            from: Box::new(simplify(*from)),
            to: Box::new(simplify(*to)),
        },
        _ => expr,
    };
    let simplified_expr = apply_rules(simplified_expr);
    simplify_rational_expression(&simplified_expr)
}

/// Applies a set of deterministic simplification rules to an expression.
#[allow(clippy::unnecessary_to_owned)]
pub(crate) fn apply_rules(expr: Expr) -> Expr {
    match expr {
        Expr::Add(a, b) => {
            if let (Some(va), Some(vb)) = (as_f64(&a), as_f64(&b)) {
                return Expr::Constant(va + vb);
            }

            let original_expr = Expr::Add(a, b);
            let (mut new_constant, mut terms) = collect_and_order_terms(&original_expr);

            // Check for sin(x)^2 + cos(x)^2 = 1 identity
            let mut changed = true;
            while changed {
                changed = false;
                let mut i = 0;
                while i < terms.len() {
                    let mut j = i + 1;
                    let mut found_match = false;
                    while j < terms.len() {
                        let (base1, coeff1) = &terms[i];
                        let (base2, coeff2) = &terms[j];
                        let mut matched = false;

                        if coeff1 == coeff2 {
                            if let (Expr::Power(b1, e1), Expr::Power(b2, e2)) = (base1, base2) {
                                let two = Expr::BigInt(BigInt::from(2));
                                let two_f = Expr::Constant(2.0);
                                //if (*e1 == Box::new(two.clone()) || **e1 == two_f)
                                if (**e1 == two || **e1 == two_f) && (**e2 == two || **e2 == two_f)
                                {
                                    if let (Expr::Sin(arg1), Expr::Cos(arg2)) = (&**b1, &**b2) {
                                        if arg1 == arg2 {
                                            matched = true;
                                        }
                                    } else if let (Expr::Cos(arg1), Expr::Sin(arg2)) =
                                        (&**b1, &**b2)
                                    {
                                        if arg1 == arg2 {
                                            matched = true;
                                        }
                                    }
                                }
                            }
                        }

                        if matched {
                            new_constant = simplify(Expr::Add(
                                Box::new(new_constant.clone()),
                                Box::new(coeff1.clone()),
                            ));
                            terms.remove(j); // Remove higher index first
                            terms.remove(i);
                            found_match = true;
                            break;
                        }
                        j += 1;
                    }
                    if found_match {
                        changed = true;
                        break;
                    }
                    i += 1;
                }
            }

            let mut term_iter = terms.into_iter().filter(|(_, coeff)| !is_zero(coeff));
            let mut result_expr = match term_iter.next() {
                Some((base, coeff)) => {
                    let first_term = if is_one(&coeff) {
                        base
                    } else {
                        Expr::Mul(Box::new(coeff), Box::new(base))
                    };
                    if !is_zero(&new_constant) {
                        Expr::Add(Box::new(new_constant), Box::new(first_term))
                    } else {
                        first_term
                    }
                }
                None => new_constant,
            };

            for (base, coeff) in term_iter {
                let term = if is_one(&coeff) {
                    base
                } else {
                    Expr::Mul(Box::new(coeff), Box::new(base))
                };
                result_expr = Expr::Add(Box::new(result_expr), Box::new(term));
            }
            result_expr
        }
        Expr::Sub(a, b) => {
            if let (Some(va), Some(vb)) = (as_f64(&a), as_f64(&b)) {
                return Expr::Constant(va - vb);
            }
            if is_zero(&b) {
                return *a;
            }
            if a == b {
                return Expr::BigInt(BigInt::zero());
            }
            // 1 - cos(x)^2 -> sin(x)^2
            if is_one(&a) {
                if let Expr::Power(base, exp) = &*b {
                    let two = Expr::BigInt(BigInt::from(2));
                    let two_f = Expr::Constant(2.0);
                    if **exp == two || **exp == two_f {
                        if let Expr::Cos(arg) = &**base {
                            return simplify(Expr::Power(
                                Box::new(Expr::Sin(arg.clone())),
                                Box::new(Expr::Constant(2.0)),
                            ));
                        }
                        if let Expr::Sin(arg) = &**base {
                            return simplify(Expr::Power(
                                Box::new(Expr::Cos(arg.clone())),
                                Box::new(Expr::Constant(2.0)),
                            ));
                        }
                    }
                }
            }
            Expr::Sub(a, b)
        }
        Expr::Mul(a, b) => {
            if let (Some(va), Some(vb)) = (as_f64(&a), as_f64(&b)) {
                return Expr::Constant(va * vb);
            }
            if is_zero(&a) || is_zero(&b) {
                return Expr::BigInt(BigInt::zero());
            }
            if is_one(&a) {
                return *b;
            }
            if is_one(&b) {
                return *a;
            }

            // exp(a) * exp(b) -> exp(a+b)
            if let (Expr::Exp(a_inner), Expr::Exp(b_inner)) = (&*a, &*b) {
                return simplify(Expr::Exp(Box::new(Expr::Add(
                    a_inner.clone(),
                    b_inner.clone(),
                ))));
            }
            // x^a * x^b -> x^(a+b)
            if let (Expr::Power(base1, exp1), Expr::Power(base2, exp2)) = (&*a, &*b) {
                if base1 == base2 {
                    return simplify(Expr::Power(
                        base1.clone(),
                        Box::new(Expr::Add(exp1.clone(), exp2.clone())),
                    ));
                }
            }

            if let Expr::Add(b_inner, c_inner) = *b {
                // Distribute a*(b+c)
                return simplify(Expr::Add(
                    Box::new(Expr::Mul(a.clone(), b_inner)),
                    Box::new(Expr::Mul(a, c_inner)),
                ));
            }
            Expr::Mul(a, b)
        }
        Expr::Div(a, b) => {
            if let (Some(va), Some(vb)) = (as_f64(&a), as_f64(&b)) {
                if vb != 0.0 {
                    return Expr::Constant(va / vb);
                }
            }
            if is_zero(&a) {
                return Expr::BigInt(BigInt::zero());
            }
            if is_one(&b) {
                return *a;
            }
            if a == b {
                return Expr::BigInt(BigInt::one());
            }
            Expr::Div(a, b)
        }
        Expr::Power(b, e) => {
            if let (Some(vb), Some(ve)) = (as_f64(&b), as_f64(&e)) {
                return Expr::Constant(vb.powf(ve));
            }
            if is_zero(&e) {
                return Expr::BigInt(BigInt::one());
            }
            if is_one(&e) {
                return *b;
            }
            if is_zero(&b) {
                return Expr::BigInt(BigInt::zero());
            }
            if is_one(&b) {
                return Expr::BigInt(BigInt::one());
            }
            if let Expr::Power(inner_b, inner_e) = *b {
                // (b^e1)^e2 = b^(e1*e2)
                return simplify(Expr::Power(inner_b, Box::new(Expr::Mul(inner_e, e))));
            }
            // (exp(x))^y -> exp(x*y)
            if let Expr::Exp(base_inner) = *b {
                return simplify(Expr::Exp(Box::new(Expr::Mul(base_inner, e))));
            }
            Expr::Power(b, e)
        }
        Expr::Sqrt(arg) => {
            let simplified_arg = simplify(*arg);
            // Try to denest the square root
            let denested = crate::symbolic::radicals::denest_sqrt(&Expr::Sqrt(Box::new(
                simplified_arg.clone(),
            )));
            if let Expr::Sqrt(_) = denested {
                // Denesting failed, apply other rules
                if let Expr::Power(ref b, ref e) = simplified_arg {
                    if let Some(val) = as_f64(e) {
                        return simplify(Expr::Power(
                            b.clone(),
                            Box::new(Expr::Constant(val / 2.0)),
                        ));
                    }
                }
                Expr::Sqrt(Box::new(simplified_arg))
            } else {
                denested
            }
        }
        Expr::Neg(arg) => {
            if let Expr::Neg(inner_arg) = *arg {
                return *inner_arg;
            } // --x = x
            if let Some(v) = as_f64(&arg) {
                return Expr::Constant(-v);
            }
            Expr::Neg(arg)
        }
        Expr::Log(arg) => {
            // Rule: Log(Complex(re, im)) -> Complex(ln(sqrt(re^2+im^2)), atan2(im, re))
            if let Expr::Complex(re, im) = &*arg {
                let magnitude_sq = Expr::Add(
                    Box::new(Expr::Power(re.clone(), Box::new(Expr::Constant(2.0)))),
                    Box::new(Expr::Power(im.clone(), Box::new(Expr::Constant(2.0)))),
                );
                let magnitude = Expr::Sqrt(Box::new(magnitude_sq));

                let real_part = Expr::Log(Box::new(magnitude));
                let imag_part = Expr::Atan2(im.clone(), re.clone());

                return simplify(Expr::Complex(Box::new(real_part), Box::new(imag_part)));
            }

            if let Expr::E = *arg {
                return Expr::BigInt(BigInt::one());
            } // ln(e) = 1
            if let Expr::Exp(inner) = *arg {
                return *inner;
            } // log(exp(x)) = x
            if is_one(&arg) {
                return Expr::BigInt(BigInt::zero());
            } // log(1) = 0
            if let Expr::Power(base, exp) = *arg {
                // log(x^y) = y*log(x)
                return simplify(Expr::Mul(exp, Box::new(Expr::Log(base))));
            }
            Expr::Log(arg)
        }
        Expr::Exp(arg) => {
            if let Expr::Log(inner) = *arg {
                return *inner;
            } // exp(log(x)) = x
            if is_zero(&arg) {
                return Expr::BigInt(BigInt::one());
            } // exp(0) = 1
            Expr::Exp(arg)
        }
        Expr::Sin(arg) => {
            if let Expr::Pi = *arg {
                return Expr::BigInt(BigInt::zero());
            } // sin(pi) = 0
            if let Expr::Neg(inner_arg) = *arg {
                // sin(-x) = -sin(x)
                return simplify(Expr::Neg(Box::new(Expr::Sin(inner_arg))));
            }
            Expr::Sin(arg)
        }
        Expr::Cos(arg) => {
            if let Expr::Pi = *arg {
                return Expr::Neg(Box::new(Expr::BigInt(BigInt::one())));
            } // cos(pi) = -1
            if let Expr::Neg(inner_arg) = *arg {
                // cos(-x) = cos(x)
                return simplify(Expr::Cos(inner_arg));
            }
            Expr::Cos(arg)
        }
        Expr::Tan(arg) => {
            if let Expr::Pi = *arg {
                return Expr::BigInt(BigInt::zero());
            } // tan(pi) = 0
            if let Expr::Neg(inner_arg) = *arg {
                // tan(-x) = -tan(x)
                return simplify(Expr::Neg(Box::new(Expr::Tan(inner_arg))));
            }
            Expr::Tan(arg)
        }
        Expr::Sum {
            body,
            var,
            from,
            to,
        } => {
            if let (Some(start), Some(end)) = (as_f64(&from), as_f64(&to)) {
                let mut total = Expr::Constant(0.0);
                for i in (start.round() as i64)..=(end.round() as i64) {
                    let i_expr = Expr::Constant(i as f64);
                    if let Expr::Variable(ref v) = *var {
                        let term = substitute(&body, v, &i_expr);
                        total = simplify(Expr::Add(Box::new(total), Box::new(term)));
                    } else {
                        return Expr::Sum {
                            body,
                            var,
                            from,
                            to,
                        }; // Cannot expand with non-variable iterator
                    }
                }
                total
            } else {
                Expr::Sum {
                    body,
                    var,
                    from,
                    to,
                } // Cannot simplify symbolic bounds
            }
        }
        _ => expr,
    }
}

// endregion

// =====================================================================================
// region: Heuristic Simplification
// =====================================================================================

pub struct RewriteRule {
    name: &'static str,
    pattern: Expr,
    replacement: Expr,
}

pub fn get_name(rule: &RewriteRule) -> String {
    println!("{}", rule.name);
    rule.name.to_string()
}
pub(crate) fn get_default_rules() -> Vec<RewriteRule> {
    vec![
        // --- Factoring and Distribution ---
        // a*b + a*c -> a*(b+c)
        RewriteRule {
            name: "factor_common_term",
            pattern: Expr::Add(
                Box::new(Expr::Mul(
                    Box::new(Expr::Pattern("a".to_string())),
                    Box::new(Expr::Pattern("b".to_string())),
                )),
                Box::new(Expr::Mul(
                    Box::new(Expr::Pattern("a".to_string())),
                    Box::new(Expr::Pattern("c".to_string())),
                )),
            ),
            replacement: Expr::Mul(
                Box::new(Expr::Pattern("a".to_string())),
                Box::new(Expr::Add(
                    Box::new(Expr::Pattern("b".to_string())),
                    Box::new(Expr::Pattern("c".to_string())),
                )),
            ),
        },
        // a*(b+c) -> a*b + a*c
        RewriteRule {
            name: "distribute_mul_add",
            pattern: Expr::Mul(
                Box::new(Expr::Pattern("a".to_string())),
                Box::new(Expr::Add(
                    Box::new(Expr::Pattern("b".to_string())),
                    Box::new(Expr::Pattern("c".to_string())),
                )),
            ),
            replacement: Expr::Add(
                Box::new(Expr::Mul(
                    Box::new(Expr::Pattern("a".to_string())),
                    Box::new(Expr::Pattern("b".to_string())),
                )),
                Box::new(Expr::Mul(
                    Box::new(Expr::Pattern("a".to_string())),
                    Box::new(Expr::Pattern("c".to_string())),
                )),
            ),
        },
        // --- Trigonometric Identities ---
        // tan(x) -> sin(x)/cos(x)
        RewriteRule {
            name: "tan_to_sin_cos",
            pattern: Expr::Tan(Box::new(Expr::Pattern("x".to_string()))),
            replacement: Expr::Div(
                Box::new(Expr::Sin(Box::new(Expr::Pattern("x".to_string())))),
                Box::new(Expr::Cos(Box::new(Expr::Pattern("x".to_string())))),
            ),
        },
        // sin(x)/cos(x) -> tan(x)
        RewriteRule {
            name: "sin_cos_to_tan",
            pattern: Expr::Div(
                Box::new(Expr::Sin(Box::new(Expr::Pattern("x".to_string())))),
                Box::new(Expr::Cos(Box::new(Expr::Pattern("x".to_string())))),
            ),
            replacement: Expr::Tan(Box::new(Expr::Pattern("x".to_string()))),
        },
        // 2*sin(x)*cos(x) -> sin(2*x)
        RewriteRule {
            name: "double_angle_sin",
            pattern: Expr::Mul(
                Box::new(Expr::BigInt(BigInt::from(2))),
                Box::new(Expr::Mul(
                    Box::new(Expr::Sin(Box::new(Expr::Pattern("x".to_string())))),
                    Box::new(Expr::Cos(Box::new(Expr::Pattern("x".to_string())))),
                )),
            ),
            replacement: Expr::Sin(Box::new(Expr::Mul(
                Box::new(Expr::BigInt(BigInt::from(2))),
                Box::new(Expr::Pattern("x".to_string())),
            ))),
        },
        // cos(x)^2 - sin(x)^2 -> cos(2*x)
        RewriteRule {
            name: "double_angle_cos_1",
            pattern: Expr::Sub(
                Box::new(Expr::Power(
                    Box::new(Expr::Cos(Box::new(Expr::Pattern("x".to_string())))),
                    Box::new(Expr::BigInt(BigInt::from(2))),
                )),
                Box::new(Expr::Power(
                    Box::new(Expr::Sin(Box::new(Expr::Pattern("x".to_string())))),
                    Box::new(Expr::BigInt(BigInt::from(2))),
                )),
            ),
            replacement: Expr::Cos(Box::new(Expr::Mul(
                Box::new(Expr::BigInt(BigInt::from(2))),
                Box::new(Expr::Pattern("x".to_string())),
            ))),
        },
        // 2*cos(x)^2 - 1 -> cos(2*x)
        RewriteRule {
            name: "double_angle_cos_2",
            pattern: Expr::Sub(
                Box::new(Expr::Mul(
                    Box::new(Expr::BigInt(BigInt::from(2))),
                    Box::new(Expr::Power(
                        Box::new(Expr::Cos(Box::new(Expr::Pattern("x".to_string())))),
                        Box::new(Expr::BigInt(BigInt::from(2))),
                    )),
                )),
                Box::new(Expr::BigInt(BigInt::from(1))),
            ),
            replacement: Expr::Cos(Box::new(Expr::Mul(
                Box::new(Expr::BigInt(BigInt::from(2))),
                Box::new(Expr::Pattern("x".to_string())),
            ))),
        },
        // 1 - 2*sin(x)^2 -> cos(2*x)
        RewriteRule {
            name: "double_angle_cos_3",
            pattern: Expr::Sub(
                Box::new(Expr::BigInt(BigInt::from(1))),
                Box::new(Expr::Mul(
                    Box::new(Expr::BigInt(BigInt::from(2))),
                    Box::new(Expr::Power(
                        Box::new(Expr::Sin(Box::new(Expr::Pattern("x".to_string())))),
                        Box::new(Expr::BigInt(BigInt::from(2))),
                    )),
                )),
            ),
            replacement: Expr::Cos(Box::new(Expr::Mul(
                Box::new(Expr::BigInt(BigInt::from(2))),
                Box::new(Expr::Pattern("x".to_string())),
            ))),
        },
    ]
}

pub fn substitute_patterns(template: &Expr, assignments: &HashMap<String, Expr>) -> Expr {
    match template {
        Expr::Pattern(name) => assignments
            .get(name)
            .cloned()
            .unwrap_or_else(|| template.clone()),
        Expr::Add(a, b) => Expr::Add(
            Box::new(substitute_patterns(a, assignments)),
            Box::new(substitute_patterns(b, assignments)),
        ),
        Expr::Sub(a, b) => Expr::Sub(
            Box::new(substitute_patterns(a, assignments)),
            Box::new(substitute_patterns(b, assignments)),
        ),
        Expr::Mul(a, b) => Expr::Mul(
            Box::new(substitute_patterns(a, assignments)),
            Box::new(substitute_patterns(b, assignments)),
        ),
        Expr::Div(a, b) => Expr::Div(
            Box::new(substitute_patterns(a, assignments)),
            Box::new(substitute_patterns(b, assignments)),
        ),
        Expr::Power(b, e) => Expr::Power(
            Box::new(substitute_patterns(b, assignments)),
            Box::new(substitute_patterns(e, assignments)),
        ),
        Expr::Sin(arg) => Expr::Sin(Box::new(substitute_patterns(arg, assignments))),
        Expr::Cos(arg) => Expr::Cos(Box::new(substitute_patterns(arg, assignments))),
        Expr::Tan(arg) => Expr::Tan(Box::new(substitute_patterns(arg, assignments))),
        Expr::Exp(arg) => Expr::Exp(Box::new(substitute_patterns(arg, assignments))),
        Expr::Log(arg) => Expr::Log(Box::new(substitute_patterns(arg, assignments))),
        Expr::Neg(arg) => Expr::Neg(Box::new(substitute_patterns(arg, assignments))),
        _ => template.clone(),
    }
}

pub(crate) fn apply_rules_recursively(expr: &Expr, rules: &[RewriteRule]) -> (Expr, bool) {
    let mut current_expr = expr.clone();
    let mut changed = false;

    // Apply to children first (post-order traversal)
    let simplified_children = match &current_expr {
        Expr::Add(a, b) => {
            let (na, ca) = apply_rules_recursively(a, rules);
            let (nb, cb) = apply_rules_recursively(b, rules);
            if ca || cb {
                Some(Expr::Add(Box::new(na), Box::new(nb)))
            } else {
                None
            }
        }
        Expr::Sub(a, b) => {
            let (na, ca) = apply_rules_recursively(a, rules);
            let (nb, cb) = apply_rules_recursively(b, rules);
            if ca || cb {
                Some(Expr::Sub(Box::new(na), Box::new(nb)))
            } else {
                None
            }
        }
        Expr::Mul(a, b) => {
            let (na, ca) = apply_rules_recursively(a, rules);
            let (nb, cb) = apply_rules_recursively(b, rules);
            if ca || cb {
                Some(Expr::Mul(Box::new(na), Box::new(nb)))
            } else {
                None
            }
        }
        Expr::Div(a, b) => {
            let (na, ca) = apply_rules_recursively(a, rules);
            let (nb, cb) = apply_rules_recursively(b, rules);
            if ca || cb {
                Some(Expr::Div(Box::new(na), Box::new(nb)))
            } else {
                None
            }
        }
        Expr::Power(b, e) => {
            let (nb, cb) = apply_rules_recursively(b, rules);
            let (ne, ce) = apply_rules_recursively(e, rules);
            if cb || ce {
                Some(Expr::Power(Box::new(nb), Box::new(ne)))
            } else {
                None
            }
        }
        Expr::Sin(arg) => {
            let (narg, carg) = apply_rules_recursively(arg, rules);
            if carg {
                Some(Expr::Sin(Box::new(narg)))
            } else {
                None
            }
        }
        Expr::Cos(arg) => {
            let (narg, carg) = apply_rules_recursively(arg, rules);
            if carg {
                Some(Expr::Cos(Box::new(narg)))
            } else {
                None
            }
        }
        Expr::Tan(arg) => {
            let (narg, carg) = apply_rules_recursively(arg, rules);
            if carg {
                Some(Expr::Tan(Box::new(narg)))
            } else {
                None
            }
        }
        Expr::Exp(arg) => {
            let (narg, carg) = apply_rules_recursively(arg, rules);
            if carg {
                Some(Expr::Exp(Box::new(narg)))
            } else {
                None
            }
        }
        Expr::Log(arg) => {
            let (narg, carg) = apply_rules_recursively(arg, rules);
            if carg {
                Some(Expr::Log(Box::new(narg)))
            } else {
                None
            }
        }
        Expr::Neg(arg) => {
            let (narg, carg) = apply_rules_recursively(arg, rules);
            if carg {
                Some(Expr::Neg(Box::new(narg)))
            } else {
                None
            }
        }
        _ => None,
    };

    if let Some(new_expr) = simplified_children {
        current_expr = new_expr;
        changed = true;
    }

    // Apply rules at the current node
    for rule in rules {
        if let Some(assignments) = pattern_match(&current_expr, &rule.pattern) {
            let new_expr = substitute_patterns(&rule.replacement, &assignments);
            let simplified_new_expr = simplify(new_expr); // Simplify the result of the rewrite

            if complexity(&simplified_new_expr) < complexity(&current_expr) {
                current_expr = simplified_new_expr;
                changed = true;
            }
        }
    }

    (current_expr, changed)
}

/// Applies a set of heuristic transformations to find a simpler form of an expression.
///
/// This function uses pattern matching and rewrite rules to transform the expression.
/// It iteratively applies rules until a fixed point is reached or a maximum number
/// of iterations is exceeded. After each pass of rule application, it performs a
/// deterministic simplification using `simplify`.
///
/// # Arguments
/// * `expr` - The expression to heuristically simplify.
///
/// # Returns
/// A new, heuristically simplified `Expr`.
pub fn heuristic_simplify(expr: Expr) -> Expr {
    let mut current_expr = expr;
    let rules = get_default_rules();
    const MAX_ITERATIONS: usize = 10;

    for _ in 0..MAX_ITERATIONS {
        let (next_expr, changed) = apply_rules_recursively(&current_expr, &rules);
        current_expr = simplify(next_expr); // Apply deterministic simplification after each pass
        if !changed {
            break; // Fixed point reached
        }
    }
    current_expr
}

// endregion

// =====================================================================================
// region: Utility Implementations
// =====================================================================================

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
            complexity(a) + 3
        }
        _ => 100,
    }
}

/// Attempts to match an expression against a pattern.
///
/// If a match is found, it returns a `HashMap` containing the assignments
/// for the pattern variables. Pattern variables are represented by `Expr::Pattern(name)`.
///
/// # Arguments
/// * `expr` - The expression to match.
/// * `pattern` - The pattern to match against.
///
/// # Returns
/// `Some(HashMap<String, Expr>)` with variable assignments if a match is found,
/// `None` otherwise.
pub fn pattern_match(expr: &Expr, pattern: &Expr) -> Option<HashMap<String, Expr>> {
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
    match (expr, pattern) {
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
        | (Expr::Neg(e), Expr::Neg(p)) => pattern_match_recursive(e, p, assignments),
        _ => expr == pattern,
    }
}

pub fn collect_and_order_terms(expr: &Expr) -> (Expr, Vec<(Expr, Expr)>) {
    /// Collects terms from an expression and orders them by complexity.
    ///
    /// This function is useful for canonicalizing expressions, especially sums and differences.
    /// It extracts a constant term and a vector of `(base, coefficient)` pairs for other terms.
    /// Terms are ordered heuristically by their complexity.
    ///
    /// # Arguments
    /// * `expr` - The expression to collect terms from.
    ///
    /// # Returns
    /// A tuple `(constant_term, terms)` where `constant_term` is an `Expr` and `terms` is a
    /// `Vec<(Expr, Expr)>` of `(base, coefficient)` pairs.
    let mut terms = BTreeMap::new();
    collect_terms_recursive(expr, &Expr::BigInt(BigInt::one()), &mut terms);
    let mut sorted_terms: Vec<(Expr, Expr)> = terms.into_iter().collect();
    sorted_terms.sort_by(|(b1, _), (b2, _)| complexity(b2).cmp(&complexity(b1)));
    let constant_term = if let Some(pos) = sorted_terms.iter().position(|(b, _)| is_one(b)) {
        let (_, c) = sorted_terms.remove(pos);
        c
    } else {
        Expr::BigInt(BigInt::zero())
    };
    (constant_term, sorted_terms)
}

fn fold_constants(expr: Expr) -> Expr {
    let expr = match expr {
        Expr::Add(a, b) => Expr::Add(Box::new(fold_constants(*a)), Box::new(fold_constants(*b))),
        Expr::Sub(a, b) => Expr::Sub(Box::new(fold_constants(*a)), Box::new(fold_constants(*b))),
        Expr::Mul(a, b) => Expr::Mul(Box::new(fold_constants(*a)), Box::new(fold_constants(*b))),
        Expr::Div(a, b) => Expr::Div(Box::new(fold_constants(*a)), Box::new(fold_constants(*b))),
        Expr::Power(base, exp) => Expr::Power(
            Box::new(fold_constants(*base)),
            Box::new(fold_constants(*exp)),
        ),
        Expr::Neg(arg) => Expr::Neg(Box::new(fold_constants(*arg))),
        _ => expr,
    };

    match expr {
        Expr::Add(a, b) => {
            if let (Some(va), Some(vb)) = (as_f64(&a), as_f64(&b)) {
                Expr::Constant(va + vb)
            } else {
                Expr::Add(a, b)
            }
        }
        Expr::Sub(a, b) => {
            if let (Some(va), Some(vb)) = (as_f64(&a), as_f64(&b)) {
                Expr::Constant(va - vb)
            } else {
                Expr::Sub(a, b)
            }
        }
        Expr::Mul(a, b) => {
            if let (Some(va), Some(vb)) = (as_f64(&a), as_f64(&b)) {
                Expr::Constant(va * vb)
            } else {
                Expr::Mul(a, b)
            }
        }
        Expr::Div(a, b) => {
            if let (Some(va), Some(vb)) = (as_f64(&a), as_f64(&b)) {
                if vb != 0.0 {
                    Expr::Constant(va / vb)
                } else {
                    Expr::Div(a, b)
                }
            } else {
                Expr::Div(a, b)
            }
        }
        Expr::Power(b, e) => {
            if let (Some(vb), Some(ve)) = (as_f64(&b), as_f64(&e)) {
                Expr::Constant(vb.powf(ve))
            } else {
                Expr::Power(b, e)
            }
        }
        Expr::Neg(arg) => {
            if let Some(v) = as_f64(&arg) {
                Expr::Constant(-v)
            } else {
                Expr::Neg(arg)
            }
        }
        _ => expr,
    }
}

pub(crate) fn collect_terms_recursive(expr: &Expr, coeff: &Expr, terms: &mut BTreeMap<Expr, Expr>) {
    match expr {
        Expr::Add(a, b) => {
            collect_terms_recursive(a, coeff, terms);
            collect_terms_recursive(b, coeff, terms);
        }
        Expr::Sub(a, b) => {
            collect_terms_recursive(a, coeff, terms);
            collect_terms_recursive(
                b,
                &fold_constants(Expr::Neg(Box::new(coeff.clone()))),
                terms,
            );
        }
        Expr::Mul(a, b) => {
            if as_f64(a).is_some() || !a.to_string().contains('x') {
                // Heuristic to identify constant-like factors
                collect_terms_recursive(
                    b,
                    &fold_constants(Expr::Mul(Box::new(coeff.clone()), Box::new(*a.clone()))),
                    terms,
                );
            } else if as_f64(b).is_some() || !b.to_string().contains('x') {
                collect_terms_recursive(
                    a,
                    &fold_constants(Expr::Mul(Box::new(coeff.clone()), Box::new(*b.clone()))),
                    terms,
                );
            } else {
                let base = expr.clone();
                let entry = terms
                    .entry(base)
                    .or_insert_with(|| Expr::BigInt(BigInt::zero()));
                *entry =
                    fold_constants(Expr::Add(Box::new(entry.clone()), Box::new(coeff.clone())));
            }
        }
        _ => {
            let base = expr.clone();
            let entry = terms
                .entry(base)
                .or_insert_with(|| Expr::BigInt(BigInt::zero()));
            *entry = fold_constants(Expr::Add(Box::new(entry.clone()), Box::new(coeff.clone())));
        }
    }
}

pub(crate) fn as_rational(expr: &Expr) -> (Expr, Expr) {
    if let Expr::Div(num, den) = expr {
        (num.as_ref().clone(), den.as_ref().clone())
    } else {
        (expr.clone(), Expr::Constant(1.0))
    }
}

pub(crate) fn simplify_rational_expression(expr: &Expr) -> Expr {
    if let Expr::Add(a, b) | Expr::Sub(a, b) | Expr::Mul(a, b) | Expr::Div(a, b) = expr {
        let (num1, den1) = as_rational(a);
        let (num2, den2) = as_rational(b);

        let (new_num_expr, new_den_expr) = match expr {
            Expr::Add(_, _) => (
                apply_rules(Expr::Add(
                    Box::new(Expr::Mul(Box::new(num1), Box::new(den2.clone()))),
                    Box::new(Expr::Mul(Box::new(num2), Box::new(den1.clone()))),
                )),
                apply_rules(Expr::Mul(Box::new(den1), Box::new(den2))),
            ),
            Expr::Sub(_, _) => (
                apply_rules(Expr::Sub(
                    Box::new(Expr::Mul(Box::new(num1), Box::new(den2.clone()))),
                    Box::new(Expr::Mul(Box::new(num2), Box::new(den1.clone()))),
                )),
                apply_rules(Expr::Mul(Box::new(den1), Box::new(den2))),
            ),
            Expr::Mul(_, _) => (
                apply_rules(Expr::Mul(Box::new(num1), Box::new(num2))),
                apply_rules(Expr::Mul(Box::new(den1), Box::new(den2))),
            ),
            Expr::Div(_, _) => (
                apply_rules(Expr::Mul(Box::new(num1), Box::new(den2.clone()))),
                apply_rules(Expr::Mul(Box::new(den1), Box::new(num2))),
            ),
            _ => unreachable!(),
        };

        if is_one(&new_den_expr) {
            return new_num_expr;
        }
        if is_zero(&new_num_expr) {
            return Expr::Constant(0.0);
        }

        // Cancel common factors using GCD
        // This assumes a single variable "x" for now. A robust solution would find all variables.
        let var = "x";
        let p_num = crate::symbolic::polynomial::expr_to_sparse_poly(&new_num_expr, &[var]);
        let p_den = crate::symbolic::polynomial::expr_to_sparse_poly(&new_den_expr, &[var]);
        let common_divisor = crate::symbolic::polynomial::gcd(p_num.clone(), p_den.clone(), var);

        if common_divisor.degree(var) > 0 {
            let final_num_poly = p_num.long_division(common_divisor.clone(), var).0;
            let final_den_poly = p_den.long_division(common_divisor, var).0;
            let final_num = crate::symbolic::polynomial::sparse_poly_to_expr(&final_num_poly);
            let final_den = crate::symbolic::polynomial::sparse_poly_to_expr(&final_den_poly);
            if is_one(&final_den) {
                return final_num;
            }
            return Expr::Div(Box::new(final_num), Box::new(final_den));
        }

        return Expr::Div(Box::new(new_num_expr), Box::new(new_den_expr));
    }
    expr.clone()
}

// endregion

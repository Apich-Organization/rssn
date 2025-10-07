//! # Symbolic Expression Simplification
//!
//! This module provides functions for symbolic expression simplification.
//! It includes a core `simplify` function that applies deterministic algebraic rules,
//! and a `heuristic_simplify` function that uses pattern matching and rewrite rules
//! to find simpler forms of expressions. It also contains utilities for term collection
//! and rational expression simplification.

use std::sync::Arc;

use crate::symbolic::calculus::substitute;
use crate::symbolic::core::Expr;
use num_bigint::BigInt;
use num_traits::{One, ToPrimitive, Zero};
use std::collections::{BTreeMap, HashMap};

// =====================================================================================
// region: Helper Functions
// =====================================================================================

#[inline]
pub fn is_zero(expr: &Expr) -> bool {
    // Checks if an expression is numerically equal to zero.
    //
    // Handles `Expr::Constant`, `Expr::BigInt`, and `Expr::Rational` variants.
    matches!(expr, Expr::Constant(val) if *val == 0.0)
        || matches!(expr, Expr::BigInt(val) if val.is_zero())
        || matches!(expr, Expr::Rational(val) if val.is_zero())
}

#[inline]
pub fn is_one(expr: &Expr) -> bool {
    // Checks if an expression is numerically equal to one.
    //
    // Handles `Expr::Constant`, `Expr::BigInt`, and `Expr::Rational` variants.
    matches!(expr, Expr::Constant(val) if *val == 1.0)
        || matches!(expr, Expr::BigInt(val) if val.is_one())
        || matches!(expr, Expr::Rational(val) if val.is_one())
}

#[inline]
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
    let mut cache = HashMap::new();
    simplify_with_cache(&expr, &mut cache)
}

/// The main simplification function with caching.
/// It recursively simplifies an expression tree by applying deterministic algebraic rules.
pub(crate) fn simplify_with_cache(expr: &Expr, cache: &mut HashMap<Expr, Expr>) -> Expr {
    if let Some(cached_result) = cache.get(expr) {
        return cached_result.clone();
    }

    let result = {
        let simplified_children_expr = match expr {
            Expr::Add(a, b) => Expr::Add(
                Arc::new(simplify_with_cache(a, cache)),
                Arc::new(simplify_with_cache(b, cache)),
            ),
            Expr::Sub(a, b) => Expr::Sub(
                Arc::new(simplify_with_cache(a, cache)),
                Arc::new(simplify_with_cache(b, cache)),
            ),
            Expr::Mul(a, b) => Expr::Mul(
                Arc::new(simplify_with_cache(a, cache)),
                Arc::new(simplify_with_cache(b, cache)),
            ),
            Expr::Div(a, b) => Expr::Div(
                Arc::new(simplify_with_cache(a, cache)),
                Arc::new(simplify_with_cache(b, cache)),
            ),
            Expr::Power(b, e) => Expr::Power(
                Arc::new(simplify_with_cache(b, cache)),
                Arc::new(simplify_with_cache(e, cache)),
            ),
            Expr::Sin(arg) => Expr::Sin(Arc::new(simplify_with_cache(arg, cache))),
            Expr::Cos(arg) => Expr::Cos(Arc::new(simplify_with_cache(arg, cache))),
            Expr::Tan(arg) => Expr::Tan(Arc::new(simplify_with_cache(arg, cache))),
            Expr::Exp(arg) => Expr::Exp(Arc::new(simplify_with_cache(arg, cache))),
            Expr::Log(arg) => Expr::Log(Arc::new(simplify_with_cache(arg, cache))),
            Expr::Neg(arg) => Expr::Neg(Arc::new(simplify_with_cache(arg, cache))),
            Expr::Sum {
                body,
                var,
                from,
                to,
            } => Expr::Sum {
                body: Arc::new(simplify_with_cache(body, cache)),
                var: Arc::new(simplify_with_cache(var, cache)),
                from: Arc::new(simplify_with_cache(from, cache)),
                to: Arc::new(simplify_with_cache(to, cache)),
            },
            _ => expr.clone(),
        };

        let simplified_expr = apply_rules(simplified_children_expr);
        simplify_rational_expression(&simplified_expr)
    };

    cache.insert(expr.clone(), result.clone());
    result
}

/// Applies a set of deterministic simplification rules to an expression.
#[allow(clippy::unnecessary_to_owned)]
pub(crate) fn apply_rules(expr: Expr) -> Expr {
    match expr {
        Expr::Add(a, b) => {
            let result_expr = match simplify_add((*a).clone(), (*b).clone()) {
                Ok(value) => value,
                Err(value) => return value,
            };
            result_expr
        }
        Expr::Sub(a, b) => {
            if let Some(value) = simplify_sub(&a, &b) {
                return value;
            }
            Expr::Sub(a, b)
        }
        Expr::Mul(a, b) => {
            if let Some(value) = simplify_mul(&a, &b) {
                return value;
            }
            Expr::Mul(a, b)
        }
        Expr::Div(a, b) => {
            if let Some(value) = simplify_div(&a, &b) {
                return value;
            }
            Expr::Div(a, b)
        }
        Expr::Power(b, e) => {
            if let Some(value) = simplify_power(&b, &e) {
                return value;
            }
            Expr::Power(b, e)
        }
        Expr::Sqrt(arg) => simplify_sqrt((*arg).clone()),
        Expr::Neg(mut arg) => {
            if matches!(*arg, Expr::Neg(_)) && crate::is_exclusive(&arg) {
                let temp_arg = arg; // Move out of the original `arg` binding
                match Arc::try_unwrap(temp_arg) {
                    Ok(Expr::Neg(inner_arc)) => {
                        return Arc::try_unwrap(inner_arc).unwrap_or_else(|a| (*a).clone());
                    }
                    Ok(other) => {
                        arg = Arc::new(other); // Re-assign if type was unexpected
                    }
                    Err(reclaimed_arg) => {
                        arg = reclaimed_arg; // Re-assign if optimization failed
                    }
                }
            }

            // Fallback path: `arg` is guaranteed to be valid here.
            if let Expr::Neg(ref inner_arg) = *arg {
                return inner_arg.as_ref().clone();
            }
            if let Some(v) = as_f64(&arg) {
                return Expr::Constant(-v);
            }
            Expr::Neg(arg)
        }
        Expr::Log(arg) => {
            if let Some(value) = simplify_log(&arg) {
                return value;
            }
            Expr::Log(arg)
        }
        Expr::Exp(arg) => {
            if let Expr::Log(ref inner) = *arg {
                return inner.as_ref().clone();
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
            if let Expr::Neg(ref inner_arg) = *arg {
                // sin(-x) = -sin(x)
                return simplify(Expr::Neg(Arc::new(Expr::Sin(inner_arg.clone()))));
            }
            Expr::Sin(arg)
        }
        Expr::Cos(arg) => {
            if let Expr::Pi = *arg {
                return Expr::Neg(Arc::new(Expr::BigInt(BigInt::one())));
            } // cos(pi) = -1
            if let Expr::Neg(ref inner_arg) = *arg {
                // cos(-x) = cos(x)
                return simplify(Expr::Cos(inner_arg.clone()));
            }
            Expr::Cos(arg)
        }
        Expr::Tan(arg) => {
            if let Expr::Pi = *arg {
                return Expr::BigInt(BigInt::zero());
            } // tan(pi) = 0
            if let Expr::Neg(ref inner_arg) = *arg {
                // tan(-x) = -tan(x)
                return simplify(Expr::Neg(Arc::new(Expr::Tan(inner_arg.clone()))));
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
                        total = simplify(Expr::Add(Arc::new(total), Arc::new(term)));
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

#[inline]
pub(crate) fn simplify_log(arg: &Expr) -> Option<Expr> {
    if let Expr::Complex(re, im) = &arg {
        let magnitude_sq = Expr::Add(
            Arc::new(Expr::Power(re.clone(), Arc::new(Expr::Constant(2.0)))),
            Arc::new(Expr::Power(im.clone(), Arc::new(Expr::Constant(2.0)))),
        );
        let magnitude = Expr::Sqrt(Arc::new(magnitude_sq));

        let real_part = Expr::Log(Arc::new(magnitude));
        let imag_part = Expr::Atan2(im.clone(), re.clone());

        return Some(simplify(Expr::Complex(
            Arc::new(real_part),
            Arc::new(imag_part),
        )));
    }
    if let Expr::E = arg {
        return Some(Expr::BigInt(BigInt::one()));
    }
    if let Expr::Exp(inner) = arg {
        return Some(inner.as_ref().clone());
    }
    if is_one(arg) {
        return Some(Expr::BigInt(BigInt::zero()));
    }
    // Rule: Log(Complex(re, im)) -> Complex(ln(sqrt(re^2+im^2)), atan2(im, re))

    // ln(e) = 1
    // log(exp(x)) = x
    // log(1) = 0
    if let Expr::Power(base, exp) = arg {
        // log(x^y) = y*log(x)
        return Some(simplify(Expr::Mul(
            exp.clone(),
            Arc::new(Expr::Log(base.clone())),
        )));
    }
    None
}

#[inline]
pub(crate) fn simplify_sqrt(arg: Expr) -> Expr {
    let simplified_arg = simplify(arg);
    // Try to denest the square root
    let denested =
        crate::symbolic::radicals::denest_sqrt(&Expr::Sqrt(Arc::new(simplified_arg.clone())));
    if let Expr::Sqrt(_) = denested {
        // Denesting failed, apply other rules
        if let Expr::Power(ref b, ref e) = simplified_arg {
            if let Some(val) = as_f64(e) {
                return simplify(Expr::Power(b.clone(), Arc::new(Expr::Constant(val / 2.0))));
            }
        }
        Expr::Sqrt(Arc::new(simplified_arg))
    } else {
        denested
    }
}

#[inline]
pub(crate) fn simplify_power(b: &Expr, e: &Expr) -> Option<Expr> {
    if let (Some(vb), Some(ve)) = (as_f64(b), as_f64(e)) {
        return Some(Expr::Constant(vb.powf(ve)));
    }
    if is_zero(e) {
        return Some(Expr::BigInt(BigInt::one()));
    }
    if is_one(e) {
        return Some(b.clone());
    }
    if is_zero(b) {
        return Some(Expr::BigInt(BigInt::zero()));
    }
    if is_one(b) {
        return Some(Expr::BigInt(BigInt::one()));
    }
    if let Expr::Power(inner_b, inner_e) = b {
        // (b^e1)^e2 = b^(e1*e2)
        return Some(simplify(Expr::Power(
            inner_b.clone(),
            Arc::new(Expr::Mul(inner_e.clone(), Arc::new(e.clone()))),
        )));
    }
    // (exp(x))^y -> exp(x*y)
    if let Expr::Exp(base_inner) = b {
        return Some(simplify(Expr::Exp(Arc::new(Expr::Mul(
            base_inner.clone(),
            Arc::new(e.clone()),
        )))));
    }
    None
}

#[inline]
pub(crate) fn simplify_div(a: &Expr, b: &Expr) -> Option<Expr> {
    if let (Some(va), Some(vb)) = (as_f64(a), as_f64(b)) {
        if vb != 0.0 {
            return Some(Expr::Constant(va / vb));
        }
    }
    if is_zero(a) {
        return Some(Expr::BigInt(BigInt::zero()));
    }
    if is_one(b) {
        return Some(a.clone());
    }
    if *a == *b {
        return Some(Expr::BigInt(BigInt::one()));
    }
    None
}

#[inline]
pub(crate) fn simplify_mul(a: &Expr, b: &Expr) -> Option<Expr> {
    if let (Some(va), Some(vb)) = (as_f64(a), as_f64(b)) {
        return Some(Expr::Constant(va * vb));
    }
    if is_zero(a) || is_zero(b) {
        return Some(Expr::BigInt(BigInt::zero()));
    }
    if is_one(a) {
        return Some(b.clone());
    }
    if is_one(b) {
        return Some(a.clone());
    }
    if let (Expr::Exp(a_inner), Expr::Exp(b_inner)) = (&a, &b) {
        return Some(simplify(Expr::Exp(Arc::new(Expr::Add(
            a_inner.clone(),
            b_inner.clone(),
        )))));
    }
    if let (Expr::Power(base1, exp1), Expr::Power(base2, exp2)) = (&a, &b) {
        if base1 == base2 {
            return Some(simplify(Expr::Power(
                base1.clone(),
                Arc::new(Expr::Add(exp1.clone(), exp2.clone())),
            )));
        }
    }

    // exp(a) * exp(b) -> exp(a+b)
    // x^a * x^b -> x^(a+b)

    if let Expr::Add(b_inner, c_inner) = b {
        // Distribute a*(b+c)
        return Some(simplify(Expr::Add(
            Arc::new(Expr::Mul(Arc::new(a.clone()), b_inner.clone())),
            Arc::new(Expr::Mul(Arc::new(a.clone()), c_inner.clone())),
        )));
    }
    None
}

#[inline]
#[allow(unused_allocation)]
pub(crate) fn simplify_sub(a: &Expr, b: &Expr) -> Option<Expr> {
    if let (Some(va), Some(vb)) = (as_f64(a), as_f64(b)) {
        return Some(Expr::Constant(va - vb));
    }
    if is_zero(b) {
        return Some(a.clone());
    }
    if *a == *b {
        return Some(Expr::BigInt(BigInt::zero()));
    }
    // 1 - cos(x)^2 -> sin(x)^2
    if is_one(a) {
        if let Expr::Power(base, exp) = &*b {
            let two = Expr::BigInt(BigInt::from(2));
            let two_f = Expr::Constant(2.0);
            if *exp == Arc::new(two) || *exp == Arc::new(two_f) {
                if let Expr::Cos(arg) = &**base {
                    return Some(simplify(Expr::Power(
                        Arc::new(Expr::Sin(arg.clone())),
                        Arc::new(Expr::Constant(2.0)),
                    )));
                }
                if let Expr::Sin(arg) = &**base {
                    return Some(simplify(Expr::Power(
                        Arc::new(Expr::Cos(arg.clone())),
                        Arc::new(Expr::Constant(2.0)),
                    )));
                }
            }
        }
    }
    None
}

#[inline]
pub(crate) fn simplify_add(a: Expr, b: Expr) -> Result<Expr, Expr> {
    if let (Some(va), Some(vb)) = (as_f64(&a), as_f64(&b)) {
        return Err(Expr::Constant(va + vb));
    }
    let original_expr = Expr::Add(Arc::new(a), Arc::new(b));
    let (mut new_constant, mut terms) = collect_and_order_terms(&original_expr);
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
                        //if (*e1 == Arc::new(two.clone()) || **e1 == two_f)
                        if (**e1 == two || **e1 == two_f) && (**e2 == two || **e2 == two_f) {
                            if let (Expr::Sin(arg1), Expr::Cos(arg2)) = (&**b1, &**b2) {
                                if arg1 == arg2 {
                                    matched = true;
                                }
                            } else if let (Expr::Cos(arg1), Expr::Sin(arg2)) = (&**b1, &**b2) {
                                if arg1 == arg2 {
                                    matched = true;
                                }
                            }
                        }
                    }
                }

                if matched {
                    new_constant = simplify(Expr::Add(
                        Arc::new(new_constant.clone()),
                        Arc::new(coeff1.clone()),
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
                Expr::Mul(Arc::new(coeff), Arc::new(base))
            };
            if !is_zero(&new_constant) {
                Expr::Add(Arc::new(new_constant), Arc::new(first_term))
            } else {
                first_term
            }
        }
        none => new_constant,
    };
    for (base, coeff) in term_iter {
        let term = if is_one(&coeff) {
            base
        } else {
            Expr::Mul(Arc::new(coeff), Arc::new(base))
        };
        result_expr = Expr::Add(Arc::new(result_expr), Arc::new(term));
    }
    Ok(result_expr)
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
                Arc::new(Expr::Mul(
                    Arc::new(Expr::Pattern("a".to_string())),
                    Arc::new(Expr::Pattern("b".to_string())),
                )),
                Arc::new(Expr::Mul(
                    Arc::new(Expr::Pattern("a".to_string())),
                    Arc::new(Expr::Pattern("c".to_string())),
                )),
            ),
            replacement: Expr::Mul(
                Arc::new(Expr::Pattern("a".to_string())),
                Arc::new(Expr::Add(
                    Arc::new(Expr::Pattern("b".to_string())),
                    Arc::new(Expr::Pattern("c".to_string())),
                )),
            ),
        },
        // a*(b+c) -> a*b + a*c
        RewriteRule {
            name: "distribute_mul_add",
            pattern: Expr::Mul(
                Arc::new(Expr::Pattern("a".to_string())),
                Arc::new(Expr::Add(
                    Arc::new(Expr::Pattern("b".to_string())),
                    Arc::new(Expr::Pattern("c".to_string())),
                )),
            ),
            replacement: Expr::Add(
                Arc::new(Expr::Mul(
                    Arc::new(Expr::Pattern("a".to_string())),
                    Arc::new(Expr::Pattern("b".to_string())),
                )),
                Arc::new(Expr::Mul(
                    Arc::new(Expr::Pattern("a".to_string())),
                    Arc::new(Expr::Pattern("c".to_string())),
                )),
            ),
        },
        // --- Trigonometric Identities ---
        // tan(x) -> sin(x)/cos(x)
        RewriteRule {
            name: "tan_to_sin_cos",
            pattern: Expr::Tan(Arc::new(Expr::Pattern("x".to_string()))),
            replacement: Expr::Div(
                Arc::new(Expr::Sin(Arc::new(Expr::Pattern("x".to_string())))),
                Arc::new(Expr::Cos(Arc::new(Expr::Pattern("x".to_string())))),
            ),
        },
        // sin(x)/cos(x) -> tan(x)
        RewriteRule {
            name: "sin_cos_to_tan",
            pattern: Expr::Div(
                Arc::new(Expr::Sin(Arc::new(Expr::Pattern("x".to_string())))),
                Arc::new(Expr::Cos(Arc::new(Expr::Pattern("x".to_string())))),
            ),
            replacement: Expr::Tan(Arc::new(Expr::Pattern("x".to_string()))),
        },
        // 2*sin(x)*cos(x) -> sin(2*x)
        RewriteRule {
            name: "double_angle_sin",
            pattern: Expr::Mul(
                Arc::new(Expr::BigInt(BigInt::from(2))),
                Arc::new(Expr::Mul(
                    Arc::new(Expr::Sin(Arc::new(Expr::Pattern("x".to_string())))),
                    Arc::new(Expr::Cos(Arc::new(Expr::Pattern("x".to_string())))),
                )),
            ),
            replacement: Expr::Sin(Arc::new(Expr::Mul(
                Arc::new(Expr::BigInt(BigInt::from(2))),
                Arc::new(Expr::Pattern("x".to_string())),
            ))),
        },
        // cos(x)^2 - sin(x)^2 -> cos(2*x)
        RewriteRule {
            name: "double_angle_cos_1",
            pattern: Expr::Sub(
                Arc::new(Expr::Power(
                    Arc::new(Expr::Cos(Arc::new(Expr::Pattern("x".to_string())))),
                    Arc::new(Expr::BigInt(BigInt::from(2))),
                )),
                Arc::new(Expr::Power(
                    Arc::new(Expr::Sin(Arc::new(Expr::Pattern("x".to_string())))),
                    Arc::new(Expr::BigInt(BigInt::from(2))),
                )),
            ),
            replacement: Expr::Cos(Arc::new(Expr::Mul(
                Arc::new(Expr::BigInt(BigInt::from(2))),
                Arc::new(Expr::Pattern("x".to_string())),
            ))),
        },
        // 2*cos(x)^2 - 1 -> cos(2*x)
        RewriteRule {
            name: "double_angle_cos_2",
            pattern: Expr::Sub(
                Arc::new(Expr::Mul(
                    Arc::new(Expr::BigInt(BigInt::from(2))),
                    Arc::new(Expr::Power(
                        Arc::new(Expr::Cos(Arc::new(Expr::Pattern("x".to_string())))),
                        Arc::new(Expr::BigInt(BigInt::from(2))),
                    )),
                )),
                Arc::new(Expr::BigInt(BigInt::from(1))),
            ),
            replacement: Expr::Cos(Arc::new(Expr::Mul(
                Arc::new(Expr::BigInt(BigInt::from(2))),
                Arc::new(Expr::Pattern("x".to_string())),
            ))),
        },
        // 1 - 2*sin(x)^2 -> cos(2*x)
        RewriteRule {
            name: "double_angle_cos_3",
            pattern: Expr::Sub(
                Arc::new(Expr::BigInt(BigInt::from(1))),
                Arc::new(Expr::Mul(
                    Arc::new(Expr::BigInt(BigInt::from(2))),
                    Arc::new(Expr::Power(
                        Arc::new(Expr::Sin(Arc::new(Expr::Pattern("x".to_string())))),
                        Arc::new(Expr::BigInt(BigInt::from(2))),
                    )),
                )),
            ),
            replacement: Expr::Cos(Arc::new(Expr::Mul(
                Arc::new(Expr::BigInt(BigInt::from(2))),
                Arc::new(Expr::Pattern("x".to_string())),
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
            Arc::new(substitute_patterns(a, assignments)),
            Arc::new(substitute_patterns(b, assignments)),
        ),
        Expr::Sub(a, b) => Expr::Sub(
            Arc::new(substitute_patterns(a, assignments)),
            Arc::new(substitute_patterns(b, assignments)),
        ),
        Expr::Mul(a, b) => Expr::Mul(
            Arc::new(substitute_patterns(a, assignments)),
            Arc::new(substitute_patterns(b, assignments)),
        ),
        Expr::Div(a, b) => Expr::Div(
            Arc::new(substitute_patterns(a, assignments)),
            Arc::new(substitute_patterns(b, assignments)),
        ),
        Expr::Power(b, e) => Expr::Power(
            Arc::new(substitute_patterns(b, assignments)),
            Arc::new(substitute_patterns(e, assignments)),
        ),
        Expr::Sin(arg) => Expr::Sin(Arc::new(substitute_patterns(arg, assignments))),
        Expr::Cos(arg) => Expr::Cos(Arc::new(substitute_patterns(arg, assignments))),
        Expr::Tan(arg) => Expr::Tan(Arc::new(substitute_patterns(arg, assignments))),
        Expr::Exp(arg) => Expr::Exp(Arc::new(substitute_patterns(arg, assignments))),
        Expr::Log(arg) => Expr::Log(Arc::new(substitute_patterns(arg, assignments))),
        Expr::Neg(arg) => Expr::Neg(Arc::new(substitute_patterns(arg, assignments))),
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
                Some(Expr::Add(Arc::new(na), Arc::new(nb)))
            } else {
                None
            }
        }
        Expr::Sub(a, b) => {
            let (na, ca) = apply_rules_recursively(a, rules);
            let (nb, cb) = apply_rules_recursively(b, rules);
            if ca || cb {
                Some(Expr::Sub(Arc::new(na), Arc::new(nb)))
            } else {
                None
            }
        }
        Expr::Mul(a, b) => {
            let (na, ca) = apply_rules_recursively(a, rules);
            let (nb, cb) = apply_rules_recursively(b, rules);
            if ca || cb {
                Some(Expr::Mul(Arc::new(na), Arc::new(nb)))
            } else {
                None
            }
        }
        Expr::Div(a, b) => {
            let (na, ca) = apply_rules_recursively(a, rules);
            let (nb, cb) = apply_rules_recursively(b, rules);
            if ca || cb {
                Some(Expr::Div(Arc::new(na), Arc::new(nb)))
            } else {
                None
            }
        }
        Expr::Power(b, e) => {
            let (nb, cb) = apply_rules_recursively(b, rules);
            let (ne, ce) = apply_rules_recursively(e, rules);
            if cb || ce {
                Some(Expr::Power(Arc::new(nb), Arc::new(ne)))
            } else {
                None
            }
        }
        Expr::Sin(arg) => {
            let (narg, carg) = apply_rules_recursively(arg, rules);
            if carg {
                Some(Expr::Sin(Arc::new(narg)))
            } else {
                None
            }
        }
        Expr::Cos(arg) => {
            let (narg, carg) = apply_rules_recursively(arg, rules);
            if carg {
                Some(Expr::Cos(Arc::new(narg)))
            } else {
                None
            }
        }
        Expr::Tan(arg) => {
            let (narg, carg) = apply_rules_recursively(arg, rules);
            if carg {
                Some(Expr::Tan(Arc::new(narg)))
            } else {
                None
            }
        }
        Expr::Exp(arg) => {
            let (narg, carg) = apply_rules_recursively(arg, rules);
            if carg {
                Some(Expr::Exp(Arc::new(narg)))
            } else {
                None
            }
        }
        Expr::Log(arg) => {
            let (narg, carg) = apply_rules_recursively(arg, rules);
            if carg {
                Some(Expr::Log(Arc::new(narg)))
            } else {
                None
            }
        }
        Expr::Neg(arg) => {
            let (narg, carg) = apply_rules_recursively(arg, rules);
            if carg {
                Some(Expr::Neg(Arc::new(narg)))
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
#[inline]
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
        Expr::Add(a, b) => Expr::Add(
            Arc::new(fold_constants(a.as_ref().clone())),
            Arc::new(fold_constants(b.as_ref().clone())),
        ),
        Expr::Sub(a, b) => Expr::Sub(
            Arc::new(fold_constants(a.as_ref().clone())),
            Arc::new(fold_constants(b.as_ref().clone())),
        ),
        Expr::Mul(a, b) => Expr::Mul(
            Arc::new(fold_constants(a.as_ref().clone())),
            Arc::new(fold_constants(b.as_ref().clone())),
        ),
        Expr::Div(a, b) => Expr::Div(
            Arc::new(fold_constants(a.as_ref().clone())),
            Arc::new(fold_constants(b.as_ref().clone())),
        ),
        Expr::Power(base, exp) => Expr::Power(
            Arc::new(fold_constants((*base).clone())),
            Arc::new(fold_constants((*exp).clone())),
        ),
        Expr::Neg(arg) => Expr::Neg(Arc::new(fold_constants((*arg).clone()))),
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
                &fold_constants(Expr::Neg(Arc::new(coeff.clone()))),
                terms,
            );
        }
        Expr::Mul(a, b) => {
            if as_f64(a).is_some() || !a.to_string().contains('x') {
                // Heuristic to identify constant-like factors
                collect_terms_recursive(
                    b,
                    &fold_constants(Expr::Mul(
                        Arc::new(coeff.clone()),
                        Arc::new(a.as_ref().clone()),
                    )),
                    terms,
                );
            } else if as_f64(b).is_some() || !b.to_string().contains('x') {
                collect_terms_recursive(
                    a,
                    &fold_constants(Expr::Mul(
                        Arc::new(coeff.clone()),
                        Arc::new(b.as_ref().clone()),
                    )),
                    terms,
                );
            } else {
                let base = expr.clone();
                let entry = terms
                    .entry(base)
                    .or_insert_with(|| Expr::BigInt(BigInt::zero()));
                *entry =
                    fold_constants(Expr::Add(Arc::new(entry.clone()), Arc::new(coeff.clone())));
            }
        }
        _ => {
            let base = expr.clone();
            let entry = terms
                .entry(base)
                .or_insert_with(|| Expr::BigInt(BigInt::zero()));
            *entry = fold_constants(Expr::Add(Arc::new(entry.clone()), Arc::new(coeff.clone())));
        }
    }
}

#[inline]
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
                    Arc::new(Expr::Mul(Arc::new(num1), Arc::new(den2.clone()))),
                    Arc::new(Expr::Mul(Arc::new(num2), Arc::new(den1.clone()))),
                )),
                apply_rules(Expr::Mul(Arc::new(den1), Arc::new(den2))),
            ),
            Expr::Sub(_, _) => (
                apply_rules(Expr::Sub(
                    Arc::new(Expr::Mul(Arc::new(num1), Arc::new(den2.clone()))),
                    Arc::new(Expr::Mul(Arc::new(num2), Arc::new(den1.clone()))),
                )),
                apply_rules(Expr::Mul(Arc::new(den1), Arc::new(den2))),
            ),
            Expr::Mul(_, _) => (
                apply_rules(Expr::Mul(Arc::new(num1), Arc::new(num2))),
                apply_rules(Expr::Mul(Arc::new(den1), Arc::new(den2))),
            ),
            Expr::Div(_, _) => (
                apply_rules(Expr::Mul(Arc::new(num1), Arc::new(den2.clone()))),
                apply_rules(Expr::Mul(Arc::new(den1), Arc::new(num2))),
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
            return Expr::Div(Arc::new(final_num), Arc::new(final_den));
        }

        return Expr::Div(Arc::new(new_num_expr), Arc::new(new_den_expr));
    }
    expr.clone()
}

// endregion

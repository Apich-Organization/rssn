//! # Symbolic Calculus Engine
//!
//! This module provides a comprehensive suite of symbolic calculus operations.
//! It includes functions for differentiation, indefinite and definite integration,
//! path integrals, limits, and series expansions. The integration capabilities
//! are supported by a multi-strategy approach including rule-based integration,
//! u-substitution, integration by parts, and more.
use crate::symbolic::core::DagOp;
use crate::symbolic::core::{Expr, PathType};
use crate::symbolic::polynomial::{is_polynomial, leading_coefficient, polynomial_degree};
use crate::symbolic::simplify::is_zero;
use crate::symbolic::simplify::normalize_constants_deep;
use crate::symbolic::simplify::simplify;
use crate::symbolic::solve::solve;
use num_bigint::BigInt;
use num_traits::{One, Zero};
use std::collections::HashMap;
use std::sync::Arc;
const ERROR_MARGIN: f64 = 1e-9;
/// Recursively substitutes all occurrences of a variable in an expression with a replacement expression.
///
/// This function traverses the expression tree and replaces every instance of the specified
/// variable `var` with the provided `replacement` expression.
///
/// # Arguments
/// * `expr` - The expression in which to perform substitutions.
/// * `var` - The name of the variable to be replaced.
/// * `replacement` - The expression to substitute in place of the variable.
///
/// # Returns
/// A new `Expr` with all occurrences of `var` replaced by `replacement`.
pub fn substitute(expr: &Expr, var: &str, replacement: &Expr) -> Expr {
    let mut stack = vec![expr.clone()];
    let mut cache: HashMap<Expr, Expr> = HashMap::new();

    while let Some(current_expr) = stack.last().cloned() {
        if cache.contains_key(&current_expr) {
            stack.pop();
            continue;
        }

        let mut children_pending = false;
        match &current_expr {
            Expr::Dag(node) => {
                return substitute(&node.to_expr().unwrap(), var, replacement);
            }
            Expr::Add(a, b)
            | Expr::Sub(a, b)
            | Expr::Mul(a, b)
            | Expr::Div(a, b)
            | Expr::Power(a, b) => {
                if !cache.contains_key(a.as_ref()) {
                    stack.push(a.as_ref().clone());
                    children_pending = true;
                }
                if !cache.contains_key(b.as_ref()) {
                    stack.push(b.as_ref().clone());
                    children_pending = true;
                }
            }
            Expr::Sin(arg)
            | Expr::Cos(arg)
            | Expr::Tan(arg)
            | Expr::Exp(arg)
            | Expr::Log(arg)
            | Expr::Neg(arg) => {
                if !cache.contains_key(arg.as_ref()) {
                    stack.push(arg.as_ref().clone());
                    children_pending = true;
                }
            }
            Expr::Integral {
                integrand,
                var: _,
                lower_bound,
                upper_bound,
            } => {
                if !cache.contains_key(integrand.as_ref()) {
                    stack.push(integrand.as_ref().clone());
                    children_pending = true;
                }
                if !cache.contains_key(lower_bound.as_ref()) {
                    stack.push(lower_bound.as_ref().clone());
                    children_pending = true;
                }
                if !cache.contains_key(upper_bound.as_ref()) {
                    stack.push(upper_bound.as_ref().clone());
                    children_pending = true;
                }
            }
            Expr::Sum {
                body,
                var: sum_var,
                from,
                to,
            } => {
                if let Expr::Variable(v) = sum_var.as_ref() {
                    if v != var {
                        if !cache.contains_key(body.as_ref()) {
                            stack.push(body.as_ref().clone());
                            children_pending = true;
                        }
                    }
                } else if !cache.contains_key(body.as_ref()) {
                    stack.push(body.as_ref().clone());
                    children_pending = true;
                }
                if !cache.contains_key(from.as_ref()) {
                    stack.push(from.as_ref().clone());
                    children_pending = true;
                }
                if !cache.contains_key(to.as_ref()) {
                    stack.push(to.as_ref().clone());
                    children_pending = true;
                }
            }
            _ => { /* Terminals or expressions with no children to substitute into */ }
        }

        if children_pending {
            continue;
        }

        let processed_expr = stack.pop().unwrap();
        let result = match &processed_expr {
            Expr::Variable(name) if name == var => replacement.clone(),
            Expr::Add(a, b) => Expr::new_add(cache[a.as_ref()].clone(), cache[b.as_ref()].clone()),
            Expr::Sub(a, b) => Expr::new_sub(cache[a.as_ref()].clone(), cache[b.as_ref()].clone()),
            Expr::Mul(a, b) => Expr::new_mul(cache[a.as_ref()].clone(), cache[b.as_ref()].clone()),
            Expr::Div(a, b) => Expr::new_div(cache[a.as_ref()].clone(), cache[b.as_ref()].clone()),
            Expr::Power(base, exp) => {
                Expr::new_pow(cache[base.as_ref()].clone(), cache[exp.as_ref()].clone())
            }
            Expr::Sin(arg) => Expr::new_sin(cache[arg.as_ref()].clone()),
            Expr::Cos(arg) => Expr::new_cos(cache[arg.as_ref()].clone()),
            Expr::Tan(arg) => Expr::new_tan(cache[arg.as_ref()].clone()),
            Expr::Exp(arg) => Expr::new_exp(cache[arg.as_ref()].clone()),
            Expr::Log(arg) => Expr::new_log(cache[arg.as_ref()].clone()),
            Expr::Neg(arg) => Expr::new_neg(cache[arg.as_ref()].clone()),
            Expr::Integral {
                integrand,
                var: int_var,
                lower_bound,
                upper_bound,
            } => Expr::Integral {
                integrand: Arc::new(cache[integrand.as_ref()].clone()),
                var: int_var.clone(),
                lower_bound: Arc::new(cache[lower_bound.as_ref()].clone()),
                upper_bound: Arc::new(cache[upper_bound.as_ref()].clone()),
            },
            Expr::Sum {
                body,
                var: sum_var,
                from,
                to,
            } => {
                let new_from = cache[from.as_ref()].clone();
                let new_to = cache[to.as_ref()].clone();
                let new_body = if let Expr::Variable(v) = &**sum_var {
                    if v == var {
                        body.clone()
                    } else {
                        Arc::new(cache[body.as_ref()].clone())
                    }
                } else {
                    Arc::new(cache[body.as_ref()].clone())
                };
                Expr::Sum {
                    body: new_body,
                    var: sum_var.clone(),
                    from: Arc::new(new_from),
                    to: Arc::new(new_to),
                }
            }
            _ => processed_expr.clone(),
        };
        cache.insert(processed_expr, result);
    }

    cache.get(expr).cloned().unwrap_or_else(|| expr.clone())
}
pub(crate) fn get_real_imag_parts(expr: &Expr) -> (Expr, Expr) {
    match simplify(expr.clone()) {
        Expr::Complex(re, im) => ((*re).clone(), (*im).clone()),
        other => (other, Expr::BigInt(BigInt::zero())),
    }
}
/// Symbolically differentiates an expression with respect to a variable.
///
/// This function implements standard differentiation rules, including the product rule, quotient rule,
/// and chain rule for nested expressions. It covers a wide range of mathematical functions
/// (polynomials, trigonometric, exponential, logarithmic, hyperbolic, inverse trigonometric, etc.).
///
/// # Arguments
/// * `expr` - The expression to differentiate.
/// * `var` - The variable with respect to which to differentiate.
///
/// # Returns
/// A new `Expr` representing the symbolic derivative.
pub fn differentiate(expr: &Expr, var: &str) -> Expr {
    let mut stack = vec![expr.clone()];
    let mut cache: HashMap<Expr, Expr> = HashMap::new();

    while let Some(current_expr) = stack.last().cloned() {
        if cache.contains_key(&current_expr) {
            stack.pop();
            continue;
        }

        let mut children_pending = false;
        match &current_expr {
            Expr::Dag(node) => {
                return differentiate(&node.to_expr().unwrap(), var);
            }
            Expr::Add(a, b)
            | Expr::Sub(a, b)
            | Expr::Mul(a, b)
            | Expr::Div(a, b)
            | Expr::Power(a, b) => {
                if !cache.contains_key(a.as_ref()) {
                    stack.push(a.as_ref().clone());
                    children_pending = true;
                }
                if !cache.contains_key(b.as_ref()) {
                    stack.push(b.as_ref().clone());
                    children_pending = true;
                }
            }
            Expr::Sin(arg)
            | Expr::Cos(arg)
            | Expr::Tan(arg)
            | Expr::Sec(arg)
            | Expr::Csc(arg)
            | Expr::Cot(arg)
            | Expr::Sinh(arg)
            | Expr::Cosh(arg)
            | Expr::Tanh(arg)
            | Expr::Exp(arg)
            | Expr::Log(arg)
            | Expr::ArcCot(arg)
            | Expr::ArcSec(arg)
            | Expr::ArcCsc(arg)
            | Expr::Coth(arg)
            | Expr::Sech(arg)
            | Expr::Csch(arg)
            | Expr::ArcSinh(arg)
            | Expr::ArcCosh(arg)
            | Expr::ArcTanh(arg)
            | Expr::ArcCoth(arg)
            | Expr::ArcSech(arg)
            | Expr::ArcCsch(arg)
            | Expr::Neg(arg) => {
                if !cache.contains_key(arg.as_ref()) {
                    stack.push(arg.as_ref().clone());
                    children_pending = true;
                }
            }
            Expr::Integral { integrand, .. } => {
                if !cache.contains_key(integrand.as_ref()) {
                    stack.push(integrand.as_ref().clone());
                    children_pending = true;
                }
            }
            Expr::Sum { body, .. } => {
                if !cache.contains_key(body.as_ref()) {
                    stack.push(body.as_ref().clone());
                    children_pending = true;
                }
            }
            _ => { /* Terminals */ }
        }

        if children_pending {
            continue;
        }

        let processed_expr = stack.pop().unwrap();
        let result = match &processed_expr {
            Expr::Constant(_) | Expr::BigInt(_) | Expr::Rational(_) | Expr::Pi | Expr::E => {
                Expr::BigInt(BigInt::zero())
            }
            Expr::Variable(name) if name == var => Expr::BigInt(BigInt::one()),
            Expr::Variable(_) => Expr::BigInt(BigInt::zero()),
            Expr::Add(a, b) => simplify(Expr::new_add(
                cache[a.as_ref()].clone(),
                cache[b.as_ref()].clone(),
            )),
            Expr::Sub(a, b) => simplify(Expr::new_sub(
                cache[a.as_ref()].clone(),
                cache[b.as_ref()].clone(),
            )),
            Expr::Mul(a, b) => simplify(Expr::new_add(
                Expr::new_mul(cache[a.as_ref()].clone(), b.as_ref().clone()),
                Expr::new_mul(a.as_ref().clone(), cache[b.as_ref()].clone()),
            )),
            Expr::Div(a, b) => simplify(Expr::new_div(
                Expr::new_sub(
                    Expr::new_mul(cache[a.as_ref()].clone(), b.as_ref().clone()),
                    Expr::new_mul(a.as_ref().clone(), cache[b.as_ref()].clone()),
                ),
                Expr::new_pow(b.as_ref().clone(), Expr::BigInt(BigInt::from(2))),
            )),
            Expr::Power(base, exp) => {
                let d_base = cache[base.as_ref()].clone();
                let d_exp = cache[exp.as_ref()].clone();
                let term1_val = Expr::new_log(base.as_ref().clone());
                let term1 = Expr::new_mul(d_exp, term1_val);
                let term1_arc = Arc::new(term1);
                let div_val = Expr::new_div(d_base, base.as_ref().clone());
                let term2_val = Expr::new_mul(exp.as_ref().clone(), div_val);
                let term2_arc = Arc::new(term2_val);
                let combined_term = Arc::new(Expr::new_add(term1_arc, term2_arc));
                simplify(Expr::new_mul(
                    Expr::new_pow(base.as_ref().clone(), exp.as_ref().clone()),
                    combined_term,
                ))
            }
            Expr::Sin(arg) => simplify(Expr::new_mul(
                Expr::new_cos(arg.as_ref().clone()),
                cache[arg.as_ref()].clone(),
            )),
            Expr::Cos(arg) => simplify(Expr::new_mul(
                Expr::new_neg(Expr::new_sin(arg.as_ref().clone())),
                cache[arg.as_ref()].clone(),
            )),
            Expr::Tan(arg) => simplify(Expr::new_mul(
                Expr::new_pow(
                    Expr::new_sec(arg.as_ref().clone()),
                    Expr::BigInt(BigInt::from(2)),
                ),
                cache[arg.as_ref()].clone(),
            )),
            Expr::Sec(arg) => simplify(Expr::new_mul(
                Expr::new_sec(arg.as_ref().clone()),
                Expr::new_mul(
                    Expr::new_tan(arg.as_ref().clone()),
                    cache[arg.as_ref()].clone(),
                ),
            )),
            Expr::Csc(arg) => simplify(Expr::new_mul(
                Expr::new_neg(Expr::new_csc(arg.as_ref().clone())),
                Expr::new_mul(
                    Expr::new_cot(arg.as_ref().clone()),
                    cache[arg.as_ref()].clone(),
                ),
            )),
            Expr::Cot(arg) => simplify(Expr::new_mul(
                Expr::new_neg(Expr::new_pow(
                    Expr::new_csc(arg.as_ref().clone()),
                    Expr::BigInt(BigInt::from(2)),
                )),
                cache[arg.as_ref()].clone(),
            )),
            Expr::Sinh(arg) => simplify(Expr::new_mul(
                Expr::new_cosh(arg.as_ref().clone()),
                cache[arg.as_ref()].clone(),
            )),
            Expr::Cosh(arg) => simplify(Expr::new_mul(
                Expr::new_sinh(arg.as_ref().clone()),
                cache[arg.as_ref()].clone(),
            )),
            Expr::Tanh(arg) => simplify(Expr::new_mul(
                Expr::new_pow(
                    Expr::new_sech(arg.as_ref().clone()),
                    Expr::BigInt(BigInt::from(2)),
                ),
                cache[arg.as_ref()].clone(),
            )),
            Expr::Exp(arg) => simplify(Expr::new_mul(
                Expr::new_exp(arg.as_ref().clone()),
                cache[arg.as_ref()].clone(),
            )),
            Expr::Log(arg) => simplify(Expr::new_div(
                cache[arg.as_ref()].clone(),
                arg.as_ref().clone(),
            )),
            Expr::ArcCot(arg) => simplify(Expr::new_neg(Expr::new_div(
                cache[arg.as_ref()].clone(),
                Expr::new_add(
                    Expr::BigInt(BigInt::one()),
                    Expr::new_pow(arg.as_ref().clone(), Expr::BigInt(BigInt::from(2))),
                ),
            ))),
            Expr::ArcSec(arg) => simplify(Expr::new_div(
                cache[arg.as_ref()].clone(),
                Expr::new_mul(
                    Expr::new_abs(arg.as_ref().clone()),
                    Expr::new_sqrt(Expr::new_sub(
                        Expr::new_pow(arg.as_ref().clone(), Expr::BigInt(BigInt::from(2))),
                        Expr::BigInt(BigInt::one()),
                    )),
                ),
            )),
            Expr::ArcCsc(arg) => simplify(Expr::new_neg(Expr::new_div(
                cache[arg.as_ref()].clone(),
                Expr::new_mul(
                    Expr::new_abs(arg.as_ref().clone()),
                    Expr::new_sqrt(Expr::new_sub(
                        Expr::new_pow(arg.as_ref().clone(), Expr::BigInt(BigInt::from(2))),
                        Expr::BigInt(BigInt::one()),
                    )),
                ),
            ))),
            Expr::Coth(arg) => simplify(Expr::new_neg(Expr::new_pow(
                Expr::new_csch(arg.as_ref().clone()),
                Expr::BigInt(BigInt::from(2)),
            ))),
            Expr::Sech(arg) => simplify(Expr::new_neg(Expr::new_mul(
                Expr::new_sech(arg.as_ref().clone()),
                Expr::new_tanh(arg.as_ref().clone()),
            ))),
            Expr::Csch(arg) => simplify(Expr::new_neg(Expr::new_mul(
                Expr::new_csch(arg.as_ref().clone()),
                Expr::new_coth(arg.as_ref().clone()),
            ))),
            Expr::ArcSinh(arg) => simplify(Expr::new_div(
                cache[arg.as_ref()].clone(),
                Expr::new_sqrt(Expr::new_add(
                    Expr::new_pow(arg.as_ref().clone(), Expr::BigInt(BigInt::from(2))),
                    Expr::BigInt(BigInt::one()),
                )),
            )),
            Expr::ArcCosh(arg) => simplify(Expr::new_div(
                cache[arg.as_ref()].clone(),
                Expr::new_sqrt(Expr::new_sub(
                    Expr::new_pow(arg.as_ref().clone(), Expr::BigInt(BigInt::from(2))),
                    Expr::BigInt(BigInt::one()),
                )),
            )),
            Expr::ArcTanh(arg) => simplify(Expr::new_div(
                cache[arg.as_ref()].clone(),
                Expr::new_sub(
                    Expr::BigInt(BigInt::one()),
                    Expr::new_pow(arg.as_ref().clone(), Expr::BigInt(BigInt::from(2))),
                ),
            )),
            Expr::ArcCoth(arg) => simplify(Expr::new_div(
                cache[arg.as_ref()].clone(),
                Expr::new_sub(
                    Expr::BigInt(BigInt::one()),
                    Expr::new_pow(arg.as_ref().clone(), Expr::BigInt(BigInt::from(2))),
                ),
            )),
            Expr::ArcSech(arg) => simplify(Expr::new_neg(Expr::new_div(
                cache[arg.as_ref()].clone(),
                Expr::new_mul(
                    arg.as_ref().clone(),
                    Expr::new_sqrt(Expr::new_sub(
                        Expr::BigInt(BigInt::one()),
                        Expr::new_pow(arg.as_ref().clone(), Expr::BigInt(BigInt::from(2))),
                    )),
                ),
            ))),
            Expr::ArcCsch(arg) => simplify(Expr::new_neg(Expr::new_div(
                cache[arg.as_ref()].clone(),
                Expr::new_mul(
                    Expr::new_abs(arg.as_ref().clone()),
                    Expr::new_sqrt(Expr::new_add(
                        Expr::BigInt(BigInt::one()),
                        Expr::new_pow(arg.as_ref().clone(), Expr::BigInt(BigInt::from(2))),
                    )),
                ),
            ))),
            Expr::Integral {
                integrand,
                var: int_var,
                ..
            } => {
                if *int_var.clone() == Expr::Variable(var.to_string()) {
                    integrand.as_ref().clone()
                } else {
                    Expr::Derivative(Arc::new(processed_expr.clone()), var.to_string())
                }
            }
            Expr::Sum {
                body,
                var: sum_var,
                from,
                to,
            } => {
                let diff_body = cache[body.as_ref()].clone();
                Expr::Sum {
                    body: Arc::new(diff_body),
                    var: sum_var.clone(),
                    from: from.clone(),
                    to: to.clone(),
                }
            }
            _ => Expr::Derivative(Arc::new(processed_expr.clone()), var.to_string()),
        };
        cache.insert(processed_expr, result);
    }

    cache.get(expr).cloned().unwrap_or_else(|| expr.clone())
}
/// Performs symbolic integration of an expression with respect to a variable.
///
/// This function acts as a dispatcher, attempting a series of integration strategies in order:
/// 1.  **Rule-based Integration**: Applies a comprehensive list of basic integration rules.
/// 2.  **U-Substitution**: Attempts to find a suitable u-substitution.
/// 3.  **Integration by Parts**: Uses the LIATE heuristic and tabular integration for applicable cases.
/// 4.  **Partial Fractions**: Decomposes rational functions, including handling of repeated roots and improper fractions (via long division).
/// 5.  **Trigonometric Substitution**: Handles integrals involving `sqrt(a^2-x^2)`, `sqrt(a^2+x^2)`, and `sqrt(x^2-a^2)`.
/// 6.  **Tangent Half-Angle Substitution**: For rational functions of trigonometric expressions.
///
/// If all strategies fail, it returns an unevaluated `Integral` expression.
///
/// # Arguments
/// * `expr` - The expression to integrate.
/// * `var` - The variable of integration.
/// * `lower_bound` - Optional: The lower bound for definite integration. If `Some`, `upper_bound` must also be `Some`.
/// * `upper_bound` - Optional: The upper bound for definite integration. If `Some`, `lower_bound` must also be `Some`.
///
/// # Returns
/// An `Expr` representing the symbolic integral.
pub fn integrate(
    expr: &Expr,
    var: &str,
    lower_bound: Option<&Expr>,
    upper_bound: Option<&Expr>,
) -> Expr {
    if let (Some(lower), Some(upper)) = (lower_bound, upper_bound) {
        return definite_integrate(expr, var, lower, upper);
    }
    let simplified_expr = simplify(expr.clone());
    if let Some(result) = integrate_by_rules(&simplified_expr, var) {
        return simplify(result);
    }
    if let Some(result) = u_substitution(&simplified_expr, var) {
        return simplify(result);
    }
    if let Some(result) = integrate_by_parts_master(&simplified_expr, var, 0) {
        return simplify(result);
    }
    if let Some(result) = integrate_by_partial_fractions(&simplified_expr, var) {
        return simplify(result);
    }
    if let Some(result) = trig_substitution(&simplified_expr, var) {
        return simplify(result);
    }
    if let Some(result) = tangent_half_angle_substitution(&simplified_expr, var) {
        return simplify(result);
    }
    let basic_result = integrate_basic(&simplified_expr, var);
    if let Expr::Integral { .. } = basic_result {
        Expr::Integral {
            integrand: Arc::new(expr.clone()),
            var: Arc::new(Expr::Variable(var.to_string())),
            lower_bound: Arc::new(Expr::Variable("a".to_string())),
            upper_bound: Arc::new(Expr::Variable("b".to_string())),
        }
    } else {
        simplify(basic_result)
    }
}
pub(crate) fn integrate_basic(expr: &Expr, var: &str) -> Expr {
    match expr {
        Expr::Dag(node) => {
            return integrate_basic(&node.to_expr().unwrap(), var);
        }
        Expr::Constant(c) => Expr::new_mul(Expr::Constant(*c), Expr::Variable(var.to_string())),
        Expr::BigInt(i) => Expr::new_mul(Expr::BigInt(i.clone()), Expr::Variable(var.to_string())),
        Expr::Rational(r) => {
            Expr::new_mul(Expr::Rational(r.clone()), Expr::Variable(var.to_string()))
        }
        Expr::Variable(name) if name == var => Expr::new_div(
            Expr::new_pow(
                Expr::Variable(var.to_string()),
                Expr::BigInt(BigInt::from(2)),
            ),
            Expr::BigInt(BigInt::from(2)),
        ),
        Expr::Add(a, b) => simplify(Expr::new_add(
            integrate(a, var, None, None),
            integrate(b, var, None, None),
        )),
        Expr::Sub(a, b) => simplify(Expr::new_sub(
            integrate(a, var, None, None),
            integrate(b, var, None, None),
        )),
        Expr::Power(base, exp) => {
            if let (Expr::Variable(name), Expr::Constant(n)) = (&**base, &**exp) {
                if name == var {
                    if (*n + 1.0).abs() < 1e-9 {
                        return Expr::new_log(Expr::new_abs(Expr::Variable(var.to_string())));
                    }
                    return Expr::new_div(
                        Expr::new_pow(Expr::Variable(var.to_string()), Expr::Constant(n + 1.0)),
                        Expr::Constant(n + 1.0),
                    );
                }
            }
            Expr::Integral {
                integrand: Arc::new(expr.clone()),
                var: Arc::new(Expr::Variable(var.to_string())),
                lower_bound: Arc::new(Expr::Variable("a".to_string())),
                upper_bound: Arc::new(Expr::Variable("b".to_string())),
            }
        }
        Expr::Exp(arg) => {
            if let Expr::Variable(name) = &**arg {
                if name == var {
                    return Expr::new_exp(Expr::Variable(var.to_string()));
                }
            }
            Expr::Integral {
                integrand: Arc::new(expr.clone()),
                var: Arc::new(Expr::Variable(var.to_string())),
                lower_bound: Arc::new(Expr::Variable("a".to_string())),
                upper_bound: Arc::new(Expr::Variable("b".to_string())),
            }
        }
        _ => Expr::Integral {
            integrand: Arc::new(expr.clone()),
            var: Arc::new(Expr::Variable(var.to_string())),
            lower_bound: Arc::new(Expr::Variable("a".to_string())),
            upper_bound: Arc::new(Expr::Variable("b".to_string())),
        },
    }
}
pub(crate) fn get_liate_type(expr: &Expr) -> i32 {
    match expr {
        Expr::Log(_) | Expr::LogBase(_, _) => 1,
        Expr::ArcSin(_) | Expr::ArcCos(_) | Expr::ArcTan(_) => 2,
        Expr::Variable(_) | Expr::Constant(_) | Expr::Power(_, _) => 3,
        Expr::Sin(_) | Expr::Cos(_) | Expr::Tan(_) => 4,
        Expr::Exp(_) => 5,
        _ => 6,
    }
}
pub(crate) fn integrate_by_parts(expr: &Expr, var: &str, depth: u32) -> Option<Expr> {
    if depth > 5 {
        return None;
    }
    if let Expr::Mul(f, g) = expr {
        let (u, dv) = if get_liate_type(f) <= get_liate_type(g) {
            (f, g)
        } else {
            (g, f)
        };
        let du_dx = differentiate(u, var);
        let v = integrate(dv, var, None, None);
        if let Expr::Integral { .. } = v {
            return None;
        }
        let uv = Expr::new_mul(u.as_ref().clone(), v.clone());
        let v_du = Expr::new_mul(v, du_dx);
        let integral_v_du = integrate(&v_du, var, None, None);
        return Some(simplify(Expr::new_sub(uv, integral_v_du)));
    }
    None
}
pub(crate) fn substitute_expr(expr: &Expr, to_replace: &Expr, replacement: &Expr) -> Expr {
    let mut stack = vec![expr.clone()];
    let mut cache = std::collections::HashMap::new();

    while let Some(current_expr) = stack.last().cloned() {
        if &current_expr == to_replace {
            cache.insert(current_expr, replacement.clone());
            stack.pop();
            continue;
        }

        if cache.contains_key(&current_expr) {
            stack.pop();
            continue;
        }

        let mut children_pending = false;
        for child in current_expr.children() {
            if !cache.contains_key(&child) {
                stack.push(child);
                children_pending = true;
            }
        }

        if children_pending {
            continue;
        }

        let processed_expr = stack.pop().unwrap();
        let op = processed_expr.op();
        let children = processed_expr.children();

        let result = match op {
            _ if processed_expr == *to_replace => replacement.clone(),
            DagOp::Add => Expr::new_add(cache[&children[0]].clone(), cache[&children[1]].clone()),
            DagOp::Sub => Expr::new_sub(cache[&children[0]].clone(), cache[&children[1]].clone()),
            DagOp::Mul => Expr::new_mul(cache[&children[0]].clone(), cache[&children[1]].clone()),
            DagOp::Div => Expr::new_div(cache[&children[0]].clone(), cache[&children[1]].clone()),
            DagOp::Power => Expr::new_pow(cache[&children[0]].clone(), cache[&children[1]].clone()),
            DagOp::Sin => Expr::new_sin(cache[&children[0]].clone()),
            _ => processed_expr.clone(),
        };
        cache.insert(processed_expr, result);
    }

    cache.get(expr).cloned().unwrap_or_else(|| expr.clone())
}
pub(crate) fn contains_var(expr: &Expr, var: &str) -> bool {
    let mut found = false;
    expr.pre_order_walk(&mut |e| {
        if let Expr::Variable(name) = e {
            if name == var {
                found = true;
            }
        }
    });
    found
}
pub(crate) fn get_u_candidates(expr: &Expr, candidates: &mut Vec<Expr>) {
    let mut stack = vec![expr.clone()];
    let mut visited = std::collections::HashSet::new();

    while let Some(current_expr) = stack.pop() {
        if !visited.insert(current_expr.clone()) {
            continue;
        }

        match &current_expr {
            Expr::Dag(node) => {
                return get_u_candidates(&node.to_expr().unwrap(), candidates);
            }
            Expr::Add(a, b) | Expr::Sub(a, b) | Expr::Mul(a, b) | Expr::Div(a, b) => {
                stack.push(a.as_ref().clone());
                stack.push(b.as_ref().clone());
            }
            Expr::Power(b, e) => {
                candidates.push(b.as_ref().clone());
                stack.push(b.as_ref().clone());
                stack.push(e.as_ref().clone());
            }
            Expr::Log(a)
            | Expr::Exp(a)
            | Expr::Sin(a)
            | Expr::Cos(a)
            | Expr::Tan(a)
            | Expr::Sec(a)
            | Expr::Csc(a)
            | Expr::Cot(a)
            | Expr::Sinh(a)
            | Expr::Cosh(a)
            | Expr::Tanh(a)
            | Expr::Sqrt(a) => {
                candidates.push(a.as_ref().clone());
                stack.push(a.as_ref().clone());
            }
            _ => {}
        }
    }
}
pub(crate) fn u_substitution(expr: &Expr, var: &str) -> Option<Expr> {
    if let Expr::Div(num, den) = expr {
        let den_prime = differentiate(den, var);
        let c = simplify(Expr::new_div(num.as_ref().clone(), den_prime));
        if !contains_var(&c, var) {
            let log_den = Expr::new_log(Expr::new_abs(den.clone()));
            return Some(simplify(Expr::new_mul(c, log_den)));
        }
    }
    let mut candidates = Vec::new();
    get_u_candidates(expr, &mut candidates);
    candidates.push(expr.clone());
    for u in candidates {
        if let Expr::Variable(_) = u {
            continue;
        }
        let du_dx = differentiate(&u, var);
        if is_zero(&du_dx) {
            continue;
        }
        let new_integrand_x = simplify(Expr::new_div(expr.clone(), du_dx));
        let temp_var = "t";
        let temp_expr = Expr::Variable(temp_var.to_string());
        let substituted = substitute_expr(&new_integrand_x, &u, &temp_expr);
        if !contains_var(&substituted, var) {
            let integral_in_t = integrate(&substituted, temp_var, None, None);
            if !matches!(integral_in_t, Expr::Integral { .. }) {
                return Some(substitute(&integral_in_t, temp_var, &u));
            }
        }
    }
    None
}
pub(crate) fn handle_trig_sub_sum(
    a_sq: &Expr,
    x_sq: &Expr,
    expr: &Expr,
    var: &str,
) -> Option<Expr> {
    if let (Expr::Constant(a_val), Expr::Power(x, two)) = (a_sq, x_sq) {
        if let (Expr::Variable(v), Expr::Constant(2.0)) = (&**x, two.as_ref().clone()) {
            if v == var && *a_val > 0.0 {
                let a = Expr::Constant(a_val.sqrt());
                let theta = Expr::Variable("theta".to_string());
                let x_sub = Expr::new_mul(a.clone(), Expr::new_tan(theta.clone()));
                let dx_dtheta = differentiate(&x_sub, "theta");
                let new_integrand =
                    simplify(Expr::new_mul(substitute(expr, var, &x_sub), dx_dtheta));
                let integral_theta = integrate(&new_integrand, "theta", None, None);
                let theta_sub = Expr::new_arctan(Expr::new_div(Expr::Variable(var.to_string()), a));
                return Some(substitute(&integral_theta, "theta", &theta_sub));
            }
        }
    }
    None
}
pub(crate) fn trig_substitution(expr: &Expr, var: &str) -> Option<Expr> {
    if let Expr::Sqrt(arg) = expr {
        if let Expr::Sub(a_sq, x_sq) = &**arg {
            if let (Expr::Constant(a_val), Expr::Power(x, two)) = (&**a_sq, &**x_sq) {
                if let (Expr::Variable(v), Expr::Constant(2.0)) = (&**x, two.as_ref().clone()) {
                    if v == var && *a_val > 0.0 {
                        let a = Expr::Constant(a_val.sqrt());
                        let theta = Expr::Variable("theta".to_string());
                        let x_sub = Expr::new_mul(a.clone(), Expr::new_sin(theta.clone()));
                        let dx_dtheta = differentiate(&x_sub, "theta");
                        let new_integrand =
                            simplify(Expr::new_mul(substitute(expr, var, &x_sub), dx_dtheta));
                        let integral_theta = integrate(&new_integrand, "theta", None, None);
                        let theta_sub =
                            Expr::new_arcsin(Expr::new_div(Expr::Variable(var.to_string()), a));
                        return Some(substitute(&integral_theta, "theta", &theta_sub));
                    }
                }
            }
        }
        if let Expr::Add(part1, part2) = &**arg {
            if let Some(result) = handle_trig_sub_sum(part1, part2, expr, var) {
                return Some(result);
            }
            if let Some(result) = handle_trig_sub_sum(part2, part1, expr, var) {
                return Some(result);
            }
        }
        if let Expr::Sub(x_sq, a_sq) = &**arg {
            if let (Expr::Power(x, two), Expr::Constant(a_val)) = (&**x_sq, &**a_sq) {
                if let (Expr::Variable(v), Expr::Constant(2.0)) = (&**x, two.as_ref().clone()) {
                    if v == var && *a_val > 0.0 {
                        let a = Expr::Constant(a_val.sqrt());
                        let theta = Expr::Variable("theta".to_string());
                        let x_sub = Expr::new_mul(a.clone(), Expr::new_sec(theta.clone()));
                        let dx_dtheta = differentiate(&x_sub, "theta");
                        let new_integrand =
                            simplify(Expr::new_mul(substitute(expr, var, &x_sub), dx_dtheta));
                        let integral_theta = integrate(&new_integrand, "theta", None, None);
                        let theta_sub =
                            Expr::new_arcsec(Expr::new_div(Expr::Variable(var.to_string()), a));
                        return Some(substitute(&integral_theta, "theta", &theta_sub));
                    }
                }
            }
        }
    }
    None
}
/// Evaluates an expression at a given point by substituting the variable with a value.
///
/// This is a wrapper around the `substitute` function, specifically for evaluating
/// an expression at a numerical or symbolic point.
///
/// # Arguments
/// * `expr` - The expression to evaluate.
/// * `var` - The variable to substitute.
/// * `value` - The value to substitute for the variable.
///
/// # Returns
/// A new `Expr` with the variable substituted by the given value.
pub fn evaluate_at_point(expr: &Expr, var: &str, value: &Expr) -> Expr {
    substitute(expr, var, value)
}
/// Computes the definite integral of an expression with respect to a variable from a lower to an upper bound.
///
/// It first finds the antiderivative (indefinite integral) using the `integrate` function.
/// Then, it applies the fundamental theorem of calculus, evaluating the antiderivative
/// at the upper and lower bounds and subtracting the results.
///
/// # Arguments
/// * `expr` - The expression to integrate.
/// * `var` - The variable of integration.
/// * `lower_bound` - The lower bound for definite integration.
/// * `upper_bound` - The upper bound for definite integration.
///
/// # Returns
/// An `Expr` representing the value of the definite integral.
/// If the indefinite integral cannot be found, it returns an unevaluated `Integral` expression.
pub fn definite_integrate(expr: &Expr, var: &str, lower_bound: &Expr, upper_bound: &Expr) -> Expr {
    let antiderivative = integrate(expr, var, None, None);
    if let Expr::Integral { .. } = antiderivative {
        return antiderivative;
    }
    let upper_eval = evaluate_at_point(&antiderivative, var, upper_bound);
    let lower_eval = evaluate_at_point(&antiderivative, var, lower_bound);
    simplify(Expr::new_sub(upper_eval, lower_eval))
}
/// Checks if a complex function `f(z)` is analytic by verifying the Cauchy-Riemann equations.
///
/// An analytic function is a function that is locally given by a convergent power series.
/// For a complex function `f(z) = u(x, y) + i*v(x, y)`, where `z = x + i*y`,
/// the Cauchy-Riemann equations state that `du/dx = dv/dy` and `du/dy = -dv/dx`.
///
/// # Arguments
/// * `expr` - The complex function `f(z)` as an `Expr`.
/// * `var` - The complex variable `z` (e.g., "z").
///
/// # Returns
/// `true` if the Cauchy-Riemann equations are satisfied (and thus the function is analytic),
/// `false` otherwise.
pub fn check_analytic(expr: &Expr, var: &str) -> bool {
    let z_replacement = Expr::new_complex(
        Expr::Variable("x".to_string()),
        Expr::Variable("y".to_string()),
    );
    let f_xy = substitute(expr, var, &z_replacement);
    let (u, v) = get_real_imag_parts(&f_xy);
    let du_dx = differentiate(&u, "x");
    let du_dy = differentiate(&u, "y");
    let dv_dx = differentiate(&v, "x");
    let dv_dy = differentiate(&v, "y");
    let cr1 = simplify(Expr::new_sub(du_dx, dv_dy));
    let cr2 = simplify(Expr::new_add(du_dy, dv_dx));
    is_zero(&cr1) && is_zero(&cr2)
}
/// Finds the poles of a rational expression by solving for the roots of the denominator.
///
/// A pole of a complex function is a point where the function's value becomes infinite.
/// For a rational function `P(z)/Q(z)`, poles occur at the roots of the denominator `Q(z)`.
///
/// # Arguments
/// * `expr` - The rational expression.
/// * `var` - The variable of the expression.
///
/// # Returns
/// A `Vec<Expr>` containing the symbolic expressions for the poles.
pub fn find_poles(expr: &Expr, var: &str) -> Vec<Expr> {
    if let Expr::Div(_, den) = expr {
        return solve(den, var);
    }
    vec![]
}
pub(crate) fn find_pole_order(expr: &Expr, var: &str, pole: &Expr) -> usize {
    let mut order = 1;
    loop {
        let term = Expr::new_pow(
            Expr::new_sub(Expr::Variable(var.to_string()), pole.clone()),
            Expr::BigInt(BigInt::from(order)),
        );
        let new_expr = simplify(Expr::new_mul(expr.clone(), term));
        let val_at_pole = simplify(evaluate_at_point(&new_expr, var, pole));
        if let Expr::Constant(c) = val_at_pole {
            if c.is_finite() && c.abs() > 1e-9 {
                return order;
            }
        }
        order += 1;
        if order > 10 {
            return 1;
        }
    }
}
/// Calculates the residue of a complex function at a given pole.
///
/// The residue is a complex number that describes the behavior of a function
/// around an isolated singularity (pole). It is crucial for evaluating contour
/// integrals using the Residue Theorem.
///
/// This function handles both simple poles (order 1) and poles of higher order `m > 1`.
/// - For a simple pole `c` of `f(z) = g(z)/h(z)`, `Res(f, c) = g(c) / h'(c)`.
/// - For a pole of order `m`, `Res(f, c) = 1/((m-1)!) * lim_{z->c} d^(m-1)/dz^(m-1) [(z-c)^m * f(z)]`.
///
/// # Arguments
/// * `expr` - The complex function.
/// * `var` - The complex variable.
/// * `pole` - The `Expr` representing the pole at which to calculate the residue.
///
/// # Returns
/// An `Expr` representing the calculated residue.
pub fn calculate_residue(expr: &Expr, var: &str, pole: &Expr) -> Expr {
    if let Expr::Div(num, den) = expr {
        let den_prime = differentiate(den, var);
        let num_at_pole = evaluate_at_point(num, var, pole);
        let den_prime_at_pole = evaluate_at_point(&den_prime, var, pole);
        if !is_zero(&simplify(den_prime_at_pole.clone())) {
            return simplify(Expr::new_div(num_at_pole, den_prime_at_pole));
        }
    }
    let m = find_pole_order(expr, var, pole);
    let m_minus_1_factorial = factorial(m - 1);
    let term = Expr::new_pow(
        Expr::new_sub(Expr::Variable(var.to_string()), pole.clone()),
        Expr::BigInt(BigInt::from(m)),
    );
    let g_z = simplify(Expr::new_mul(expr.clone(), term));
    let mut g_m_minus_1 = g_z;
    for _ in 0..(m - 1) {
        g_m_minus_1 = differentiate(&g_m_minus_1, var);
    }
    let limit = evaluate_at_point(&g_m_minus_1, var, pole);
    simplify(Expr::new_div(limit, Expr::Constant(m_minus_1_factorial)))
}
///
/// This function is used in complex analysis, particularly with the Residue Theorem,
/// to determine which poles of a function lie within a given integration path.
///
/// # Arguments
/// * `point` - The complex point to check, as an `Expr::Complex`.
/// * `contour` - The contour, currently supporting `Expr::Path(PathType::Circle, center, radius)`.
///
/// # Returns
/// `true` if the point is strictly inside the contour, `false` otherwise.
pub fn is_inside_contour(point: &Expr, contour: &Expr) -> bool {
    if let (Expr::Path(PathType::Circle, center, radius), Expr::Complex(re, im)) = (contour, point)
    {
        if let (Expr::Complex(center_re, center_im), Expr::Constant(r)) = (&**center, &**radius) {
            let dist_sq = Expr::new_add(
                Expr::new_pow(
                    Expr::new_sub(re.clone(), center_re.clone()),
                    Expr::BigInt(BigInt::from(2)),
                ),
                Expr::new_pow(
                    Expr::new_sub(im.clone(), center_im.clone()),
                    Expr::BigInt(BigInt::from(2)),
                ),
            );
            if let Expr::Constant(d2) = simplify(dist_sq) {
                return d2 < r * r;
            }
        }
    }
    false
}
/// Computes a path integral of a complex function over a given contour.
///
/// This function implements different strategies based on the type of contour:
/// - **Circular Contours**: Uses the Residue Theorem. If the function is analytic
///   inside the contour, the integral is zero. Otherwise, it sums the residues
///   of poles inside the contour and multiplies by `2 * pi * i`.
/// - **Line Segment Contours**: Parameterizes the line segment and converts the
///   path integral into a definite integral over a real variable.
/// - **Rectangular Contours**: Decomposes the rectangle into four line segments
///   and sums the integrals over each segment.
///
/// # Arguments
/// * `expr` - The complex function to integrate.
/// * `var` - The complex variable of integration.
/// * `contour` - An `Expr::Path` defining the integration contour.
///
/// # Returns
/// An `Expr` representing the value of the path integral.
pub fn path_integrate(expr: &Expr, var: &str, contour: &Expr) -> Expr {
    match contour {
        Expr::Path(path_type, param1, param2) => match path_type {
            PathType::Circle => {
                if check_analytic(expr, var) {
                    return Expr::BigInt(BigInt::zero());
                }
                let poles = find_poles(expr, var);
                let mut sum_of_residues = Expr::BigInt(BigInt::zero());
                for pole in poles {
                    if is_inside_contour(&pole, contour) {
                        let residue = calculate_residue(expr, var, &pole);
                        sum_of_residues = Expr::new_add(sum_of_residues, residue);
                    }
                }
                let two_pi_i = Expr::new_mul(
                    Expr::Constant(2.0 * std::f64::consts::PI),
                    Expr::new_complex(Expr::BigInt(BigInt::zero()), Expr::BigInt(BigInt::one())),
                );
                simplify(Expr::new_mul(two_pi_i, sum_of_residues))
            }
            PathType::Line => {
                let (z0, z1) = (&**param1, &**param2);
                let dz_dt = simplify(Expr::new_sub(z1.clone(), z0.clone()));
                let t_var = Expr::Variable("t".to_string());
                let z_t = simplify(Expr::new_add(
                    z0.clone(),
                    Expr::new_mul(t_var.clone(), dz_dt.clone()),
                ));
                let integrand_t = simplify(Expr::new_mul(substitute(expr, var, &z_t), dz_dt));
                definite_integrate(
                    &integrand_t,
                    "t",
                    &Expr::BigInt(BigInt::zero()),
                    &Expr::BigInt(BigInt::one()),
                )
            }
            PathType::Rectangle => {
                let (z_bl, z_tr) = (&**param1, &**param2);
                let z_br = Expr::new_complex(z_tr.re(), z_bl.im());
                let z_tl = Expr::new_complex(z_bl.re(), z_tr.im());
                let i1 = path_integrate(
                    expr,
                    var,
                    &Expr::Path(
                        PathType::Line,
                        Arc::new(z_bl.clone()),
                        Arc::new(z_br.clone()),
                    ),
                );
                let i2 = path_integrate(
                    expr,
                    var,
                    &Expr::Path(PathType::Line, Arc::new(z_br), Arc::new(z_tr.clone())),
                );
                let i3 = path_integrate(
                    expr,
                    var,
                    &Expr::Path(
                        PathType::Line,
                        Arc::new(z_tr.clone()),
                        Arc::new(z_tl.clone()),
                    ),
                );
                let i4 = path_integrate(
                    expr,
                    var,
                    &Expr::Path(PathType::Line, Arc::new(z_tl), Arc::new(z_bl.clone())),
                );
                simplify(Expr::new_add(i1, Expr::new_add(i2, Expr::new_add(i3, i4))))
            }
        },
        _ => Expr::Integral {
            integrand: Arc::new(expr.clone()),
            var: Arc::new(Expr::Variable(var.to_string())),
            lower_bound: Arc::new(Expr::Variable("C_lower".to_string())),
            upper_bound: Arc::new(Expr::Variable("C_upper".to_string())),
        },
    }
}
/// Computes the factorial of a non-negative integer `n`.
///
/// The factorial of `n` (denoted as `n!`) is the product of all positive integers
/// less than or equal to `n`. `0!` is defined as 1.
///
/// # Arguments
/// * `n` - The non-negative integer.
///
/// # Returns
/// The factorial of `n` as an `f64`.
pub fn factorial(n: usize) -> f64 {
    if n > 170 {
        return f64::INFINITY;
    }
    if n == 0 {
        1.0
    } else {
        (1..=n).map(|i| i as f64).product::<f64>()
    }
}
impl From<f64> for Expr {
    fn from(val: f64) -> Self {
        Expr::Constant(val)
    }
}
/// Calculates an improper integral from -infinity to +infinity using the residue theorem.
///
/// This function is designed for integrands `f(z)` that satisfy the following conditions:
/// 1. `f(z)` is analytic in the upper half-plane except for a finite number of poles.
/// 2. `f(z)` has no poles on the real axis.
/// 3. The integral of `f(z)` over a semi-circular arc in the upper half-plane vanishes as the radius tends to infinity.
///    A common condition for this is that `|f(z)|` behaves like `1/|z|^k` where `k >= 2` as `|z| -> infinity`.
///
/// The integral is calculated as `2 * pi * i * (sum of residues in the upper half-plane)`.
///
/// # Arguments
/// * `expr` - The expression to integrate.
/// * `var` - The variable of integration.
///
/// # Returns
/// An `Expr` representing the value of the improper integral.
pub fn improper_integral(expr: &Expr, var: &str) -> Expr {
    pub(crate) fn get_imag_part(expr: &Expr) -> Option<f64> {
        match simplify(expr.clone()) {
            Expr::Complex(_, im_part) => {
                if let Expr::Constant(val) = *im_part {
                    Some(val)
                } else if is_zero(&im_part) {
                    Some(0.0)
                } else {
                    None
                }
            }
            Expr::Constant(_val) => Some(0.0),
            expr if is_zero(&expr) => Some(0.0),
            _ => None,
        }
    }
    let poles = find_poles(expr, var);
    let mut sum_of_residues_in_uhp = Expr::BigInt(BigInt::zero());
    for pole in poles {
        if let Some(im_val) = get_imag_part(&pole) {
            if im_val > 1e-9 {
                let residue = calculate_residue(expr, var, &pole);
                sum_of_residues_in_uhp =
                    simplify(Expr::new_add(sum_of_residues_in_uhp.clone(), residue));
            }
        }
    }
    let two_pi_i = Expr::new_mul(
        Expr::Constant(2.0 * std::f64::consts::PI),
        Expr::new_complex(Expr::BigInt(BigInt::zero()), Expr::BigInt(BigInt::one())),
    );
    simplify(Expr::new_mul(two_pi_i, sum_of_residues_in_uhp))
}
pub(crate) fn integrate_by_rules(expr: &Expr, var: &str) -> Option<Expr> {
    match expr {
        Expr::Dag(node) => {
            return integrate_by_rules(&node.to_expr().unwrap(), var);
        }
        Expr::Constant(c) => Some(Expr::new_mul(
            Expr::Constant(*c),
            Expr::Variable(var.to_string()),
        )),
        Expr::BigInt(i) => Some(Expr::new_mul(
            Expr::BigInt(i.clone()),
            Expr::Variable(var.to_string()),
        )),
        Expr::Rational(r) => Some(Expr::new_mul(
            Expr::Rational(r.clone()),
            Expr::Variable(var.to_string()),
        )),
        Expr::Variable(name) if name == var => Some(Expr::new_div(
            Expr::new_pow(
                Expr::Variable(var.to_string()),
                Expr::BigInt(BigInt::from(2)),
            ),
            Expr::BigInt(BigInt::from(2)),
        )),
        Expr::Exp(arg) => {
            if let Expr::Variable(name) = &**arg {
                if name == var {
                    return Some(Expr::new_exp(Expr::Variable(var.to_string())));
                }
            }
            if let Expr::Mul(a, x) = &**arg {
                if let (Expr::Constant(coeff), Expr::Variable(v)) = (&**a, &**x) {
                    if v == var {
                        return Some(Expr::new_div(expr.clone(), Expr::Constant(*coeff)));
                    }
                }
                if let (Expr::Variable(v), Expr::Constant(coeff)) = (&**x, &**a) {
                    if v == var {
                        return Some(Expr::new_div(expr.clone(), Expr::Constant(*coeff)));
                    }
                }
            }
            None
        }
        Expr::Log(arg) => {
            if let Expr::Variable(name) = &**arg {
                if name == var {
                    let x = Expr::Variable(var.to_string());
                    return Some(Expr::new_sub(Expr::new_mul(x.clone(), expr.clone()), x));
                }
            }
            None
        }
        Expr::LogBase(base, arg) => {
            if let Expr::Variable(name) = &**arg {
                if name == var && !contains_var(base, var) {
                    let ln_x = Expr::new_log(arg.clone());
                    let ln_b = Expr::new_log(base.clone());
                    let new_expr = Expr::new_div(ln_x, ln_b);
                    return integrate(&new_expr, var, None, None).into();
                }
            }
            None
        }
        Expr::Div(num, den) => {
            if let (Expr::BigInt(one), Expr::Variable(name)) = (&**num, &**den) {
                if one.is_one() && name == var {
                    return Some(Expr::new_log(Expr::new_abs(Expr::Variable(
                        var.to_string(),
                    ))));
                }
            }
            if let Expr::BigInt(one) = &**num {
                if one.is_one() {
                    if let Expr::Add(part1, part2) = &**den {
                        let (a_sq_box, x_sq_box) = if let Expr::Power(_, _) = &**part1 {
                            (part2, part1)
                        } else {
                            (part1, part2)
                        };
                        if let (Expr::Constant(a_val), Expr::Power(x, two)) =
                            (&**a_sq_box, &**x_sq_box)
                        {
                            if let (Expr::Variable(v), Expr::Constant(val)) = (&**x, &**two) {
                                if v == var && (*val - 2.0).abs() < ERROR_MARGIN {
                                    let a = Expr::Constant(a_val.sqrt());
                                    return Some(Expr::new_mul(
                                        Expr::new_div(Expr::BigInt(BigInt::one()), a.clone()),
                                        Expr::new_arctan(Expr::new_div(
                                            Expr::Variable(var.to_string()),
                                            a,
                                        )),
                                    ));
                                }
                            }
                        }
                    }
                }
            }
            if let (Expr::BigInt(one), Expr::Sqrt(sqrt_arg)) = (&**num, &**den) {
                if one.is_one() {
                    if let Expr::Sub(a_sq, x_sq) = &**sqrt_arg {
                        if let (Expr::Constant(a_val), Expr::Power(x, two)) = (&**a_sq, &**x_sq) {
                            if let (Expr::Variable(v), Expr::Constant(val)) = (&**x, &**two) {
                                if v == var && (*val - 2.0).abs() < ERROR_MARGIN {
                                    let a = Expr::Constant(a_val.sqrt());
                                    return Some(Expr::new_arcsin(Expr::new_div(
                                        Expr::Variable(var.to_string()),
                                        a,
                                    )));
                                }
                            }
                        }
                    }
                }
            }
            None
        }
        Expr::Sin(arg) => {
            if let Expr::Variable(name) = &**arg {
                if name == var {
                    return Some(Expr::new_neg(Expr::new_cos(Expr::Variable(
                        var.to_string(),
                    ))));
                }
            }
            None
        }
        Expr::Cos(arg) => {
            if let Expr::Variable(name) = &**arg {
                if name == var {
                    return Some(Expr::new_sin(Expr::Variable(var.to_string())));
                }
            }
            None
        }
        Expr::Tan(arg) => {
            if let Expr::Variable(name) = &**arg {
                if name == var {
                    return Some(Expr::new_log(Expr::new_abs(Expr::new_sec(Expr::Variable(
                        var.to_string(),
                    )))));
                }
            }
            None
        }
        Expr::Sec(arg) => {
            if let Expr::Variable(name) = &**arg {
                if name == var {
                    return Some(Expr::new_log(Expr::new_abs(Expr::new_add(
                        Expr::new_sec(Expr::Variable(var.to_string())),
                        Expr::new_tan(Expr::Variable(var.to_string())),
                    ))));
                }
            }
            None
        }
        Expr::Csc(arg) => {
            if let Expr::Variable(name) = &**arg {
                if name == var {
                    return Some(Expr::new_log(Expr::new_abs(Expr::new_sub(
                        Expr::new_csc(Expr::Variable(var.to_string())),
                        Expr::new_cot(Expr::Variable(var.to_string())),
                    ))));
                }
            }
            None
        }
        Expr::Cot(arg) => {
            if let Expr::Variable(name) = &**arg {
                if name == var {
                    return Some(Expr::new_log(Expr::new_abs(Expr::new_sin(Expr::Variable(
                        var.to_string(),
                    )))));
                }
            }
            None
        }
        Expr::Power(base, exp)
            if matches!(&**base, Expr::Sec(_))
                && matches!(&** exp, Expr::BigInt(b) if * b == BigInt::from(2)) =>
        {
            if let Expr::Sec(arg) = &**base {
                if let Expr::Variable(name) = &**arg {
                    if name == var {
                        return Some(Expr::new_tan(Expr::Variable(var.to_string())));
                    }
                }
            }
            None
        }
        Expr::Power(base, exp)
            if matches!(&**base, Expr::Csc(_))
                && matches!(&** exp, Expr::BigInt(b) if * b == BigInt::from(2)) =>
        {
            if let Expr::Csc(arg) = &**base {
                if let Expr::Variable(name) = &**arg {
                    if name == var {
                        return Some(Expr::new_neg(Expr::new_cot(Expr::Variable(
                            var.to_string(),
                        ))));
                    }
                }
            }
            None
        }
        Expr::Mul(a, b)
            if (matches!(&**a, Expr::Sec(_)) && matches!(&**b, Expr::Tan(_)))
                || (matches!(&**b, Expr::Sec(_)) && matches!(&**a, Expr::Tan(_))) =>
        {
            let (sec_part, tan_part) = if let Expr::Sec(_) = &**a {
                (a, b)
            } else {
                (b, a)
            };
            if let (Expr::Sec(arg1), Expr::Tan(arg2)) = (&**sec_part, &**tan_part) {
                if arg1 == arg2 {
                    if let Expr::Variable(name) = &**arg1 {
                        if name == var {
                            return Some(Expr::new_sec(Expr::Variable(var.to_string())));
                        }
                    }
                }
            }
            None
        }
        Expr::Mul(a, b)
            if (matches!(&**a, Expr::Csc(_)) && matches!(&**b, Expr::Cot(_)))
                || (matches!(&**b, Expr::Csc(_)) && matches!(&**a, Expr::Cot(_))) =>
        {
            let (csc_part, cot_part) = if let Expr::Csc(_) = &**a {
                (a, b)
            } else {
                (b, a)
            };
            if let (Expr::Csc(arg1), Expr::Cot(arg2)) = (&**csc_part, &**cot_part) {
                if arg1 == arg2 {
                    if let Expr::Variable(name) = &**arg1 {
                        if name == var {
                            return Some(Expr::new_neg(Expr::new_csc(Expr::Variable(
                                var.to_string(),
                            ))));
                        }
                    }
                }
            }
            None
        }
        Expr::ArcTan(arg) => {
            if let Expr::Variable(name) = &**arg {
                if name == var {
                    let x = Expr::Variable(var.to_string());
                    let term1 = Expr::new_mul(x.clone(), expr.clone());
                    let term2 = Expr::new_mul(
                        Expr::Constant(0.5),
                        Expr::new_log(Expr::new_add(
                            Expr::BigInt(BigInt::one()),
                            Expr::new_pow(x, Expr::BigInt(BigInt::from(2))),
                        )),
                    );
                    return Some(Expr::new_sub(term1, term2));
                }
            }
            None
        }
        Expr::ArcSin(arg) => {
            if let Expr::Variable(name) = &**arg {
                if name == var {
                    let x = Expr::Variable(var.to_string());
                    let term1 = Expr::new_mul(x.clone(), expr.clone());
                    let term2 = Expr::new_sqrt(Expr::new_sub(
                        Expr::BigInt(BigInt::one()),
                        Expr::new_pow(x, Expr::BigInt(BigInt::from(2))),
                    ));
                    return Some(Expr::new_add(term1, term2));
                }
            }
            None
        }
        Expr::Sinh(arg) => {
            if let Expr::Variable(name) = &**arg {
                if name == var {
                    return Some(Expr::new_cosh(Expr::Variable(var.to_string())));
                }
            }
            None
        }
        Expr::Cosh(arg) => {
            if let Expr::Variable(name) = &**arg {
                if name == var {
                    return Some(Expr::new_sinh(Expr::Variable(var.to_string())));
                }
            }
            None
        }
        Expr::Tanh(arg) => {
            if let Expr::Variable(name) = &**arg {
                if name == var {
                    return Some(Expr::new_log(Expr::new_cosh(Expr::Variable(
                        var.to_string(),
                    ))));
                }
            }
            None
        }
        Expr::Power(base, exp)
            if matches!(&**base, Expr::Sech(_))
                && matches!(&** exp, Expr::BigInt(b) if * b == BigInt::from(2)) =>
        {
            if let Expr::Sech(arg) = &**base {
                if let Expr::Variable(name) = &**arg {
                    if name == var {
                        return Some(Expr::new_tanh(Expr::Variable(var.to_string())));
                    }
                }
            }
            None
        }
        Expr::Power(base, exp) => {
            if let (Expr::Variable(name), Expr::Constant(n)) = (&**base, &**exp) {
                if name == var {
                    if (*n + 1.0).abs() < 1e-9 {
                        return Some(Expr::new_log(Expr::new_abs(Expr::Variable(
                            var.to_string(),
                        ))));
                    }
                    return Some(Expr::new_div(
                        Expr::new_pow(Expr::Variable(var.to_string()), Expr::Constant(n + 1.0)),
                        Expr::Constant(n + 1.0),
                    ));
                }
            }
            None
        }
        _ => None,
    }
}
pub(crate) fn integrate_by_parts_tabular(expr: &Expr, var: &str) -> Option<Expr> {
    if let Expr::Mul(part1, part2) = expr {
        let (poly_part, other_part) = if is_polynomial(part1, var) {
            (part1, part2)
        } else if is_polynomial(part2, var) {
            (part2, part1)
        } else {
            return None;
        };
        let mut derivatives = vec![poly_part.clone()];
        while let Some(last_deriv) = derivatives.last() {
            if is_zero(&simplify((**last_deriv).clone())) {
                break;
            }
            derivatives.push(Arc::new(differentiate(last_deriv, var)));
        }
        derivatives.pop();
        let mut integrals = vec![other_part.clone()];
        for _ in 0..derivatives.len() {
            if let Some(last_integral) = integrals.last() {
                let next_integral = integrate(last_integral, var, None, None);
                if let Expr::Integral { .. } = next_integral {
                    return None;
                }
                integrals.push(Arc::new(simplify(next_integral)));
            } else {
                return None;
            }
        }
        if derivatives.len() >= integrals.len() {
            return None;
        }
        let mut total = Expr::BigInt(BigInt::zero());
        let mut sign = 1;
        for i in 0..derivatives.len() {
            let term = Expr::new_mul(
                derivatives[i].as_ref().clone(),
                integrals[i + 1].as_ref().clone(),
            );
            if sign == 1 {
                total = Expr::new_add(total, term);
            } else {
                total = Expr::new_sub(total, term);
            }
            sign *= -1;
        }
        return Some(simplify(total));
    }
    None
}
pub(crate) fn integrate_by_parts_master(expr: &Expr, var: &str, depth: u32) -> Option<Expr> {
    if depth == 0 {
        if let Some(result) = integrate_by_parts_tabular(expr, var) {
            return Some(result);
        }
    }
    integrate_by_parts(expr, var, depth)
}
pub(crate) fn find_roots_with_multiplicity(expr: &Expr, var: &str) -> Vec<(Expr, usize)> {
    let unique_poles = solve(expr, var);
    let mut roots_with_multiplicity = Vec::new();
    let mut processed_poles = std::collections::HashSet::new();
    for pole in unique_poles {
        if !processed_poles.insert(pole.clone()) {
            continue;
        }
        let mut m = 1;
        let mut current_deriv = expr.clone();
        while m < 10 {
            let next_deriv = differentiate(&current_deriv, var);
            let val_at_pole = simplify(evaluate_at_point(&next_deriv, var, &pole));
            if !is_zero(&val_at_pole) {
                break;
            }
            m += 1;
            current_deriv = next_deriv;
        }
        roots_with_multiplicity.push((pole, m));
    }
    roots_with_multiplicity
}
pub(crate) fn integrate_by_partial_fractions(expr: &Expr, var: &str) -> Option<Expr> {
    if let Expr::Div(num, den) = expr {
        use crate::symbolic::polynomial::{polynomial_degree, polynomial_long_division_coeffs};
        let num_deg = polynomial_degree(num, var);
        let den_deg = polynomial_degree(den, var);
        if num_deg >= 0 && den_deg >= 0 && num_deg >= den_deg {
            if let Ok((quotient, remainder)) = polynomial_long_division_coeffs(num, den, var) {
                let integral_of_quotient = integrate(&quotient, var, None, None);
                let integral_of_remainder = if is_zero(&remainder) {
                    Expr::BigInt(BigInt::zero())
                } else {
                    let remainder_fraction = Expr::new_div(remainder, den.as_ref().clone());
                    integrate(&remainder_fraction, var, None, None)
                };
                return Some(simplify(Expr::new_add(
                    integral_of_quotient,
                    integral_of_remainder,
                )));
            } else {
                return None;
            }
        }
        let roots = find_roots_with_multiplicity(den, var);
        if roots.is_empty() {
            return None;
        }
        let mut total_integral = Expr::BigInt(BigInt::zero());
        for (root, m) in roots {
            let term_to_multiply = Expr::new_pow(
                Expr::new_sub(Expr::Variable(var.to_string()), root.clone()),
                Expr::BigInt(BigInt::from(m)),
            );
            let g_z = simplify(Expr::new_mul(expr.clone(), term_to_multiply));
            for k in 0..m {
                let mut deriv_g = g_z.clone();
                for _ in 0..k {
                    deriv_g = differentiate(&deriv_g, var);
                }
                let val_at_root = evaluate_at_point(&deriv_g, var, &root);
                let k_factorial = Expr::Constant(factorial(k));
                let coefficient = simplify(Expr::new_div(val_at_root, k_factorial));
                let j = m - k;
                let integral_term = if j == 1 {
                    let log_arg = Expr::new_abs(simplify(Expr::new_sub(
                        Expr::Variable(var.to_string()),
                        root.clone(),
                    )));
                    Expr::new_mul(coefficient, Expr::new_log(log_arg))
                } else {
                    let j_i64 = match i64::try_from(j) {
                        Ok(val) => val,
                        Err(_) => {
                            eprintln!(
                                "Warning: usize value {} is too large to fit in i64 during partial fraction integration. Returning None.",
                                j
                            );
                            return None;
                        }
                    };
                    let new_power = 1_i64 - j_i64;
                    let new_denom = Expr::Constant(new_power as f64);
                    let integrated_power_term = Expr::new_pow(
                        Expr::new_sub(Expr::Variable(var.to_string()), root.clone()),
                        Expr::Constant(new_power as f64),
                    );
                    Expr::new_mul(coefficient, Expr::new_div(integrated_power_term, new_denom))
                };
                total_integral = simplify(Expr::new_add(total_integral, integral_term));
            }
        }
        return Some(total_integral);
    }
    None
}
pub(crate) fn contains_trig_function(expr: &Expr) -> bool {
    let mut stack = vec![expr.clone()];
    let mut visited = std::collections::HashSet::new();

    while let Some(current_expr) = stack.pop() {
        if !visited.insert(current_expr.clone()) {
            continue;
        }

        match &current_expr {
            Expr::Sin(_)
            | Expr::Cos(_)
            | Expr::Tan(_)
            | Expr::Sec(_)
            | Expr::Csc(_)
            | Expr::Cot(_) => {
                return true;
            }
            Expr::Add(a, b) | Expr::Sub(a, b) | Expr::Mul(a, b) | Expr::Div(a, b) => {
                stack.push(a.as_ref().clone());
                stack.push(b.as_ref().clone());
            }
            Expr::Power(base, exp) => {
                stack.push(base.as_ref().clone());
                stack.push(exp.as_ref().clone());
            }
            Expr::Log(arg) | Expr::Abs(arg) | Expr::Neg(arg) | Expr::Exp(arg) => {
                stack.push(arg.as_ref().clone());
            }
            Expr::Complex(re, im) => {
                stack.push(re.as_ref().clone());
                stack.push(im.as_ref().clone());
            }
            _ => {}
        }
    }

    false
}
pub(crate) fn tangent_half_angle_substitution(expr: &Expr, var: &str) -> Option<Expr> {
    if !contains_trig_function(expr) {
        return None;
    }
    let t = Expr::Variable("t".to_string());
    let t_squared = Expr::new_pow(t.clone(), Expr::BigInt(BigInt::from(2)));
    let one_plus_t_squared = Expr::new_add(Expr::BigInt(BigInt::one()), t_squared.clone());
    let sin_x_sub = Expr::new_div(
        Expr::new_mul(Expr::BigInt(BigInt::from(2)), t.clone()),
        one_plus_t_squared.clone(),
    );
    let cos_x_sub = Expr::new_div(
        Expr::new_sub(Expr::BigInt(BigInt::one()), t_squared.clone()),
        one_plus_t_squared.clone(),
    );
    let tan_x_sub = simplify(Expr::new_div(sin_x_sub.clone(), cos_x_sub.clone()));
    let dx_sub = Expr::new_div(Expr::BigInt(BigInt::from(2)), one_plus_t_squared.clone());
    let mut sub_expr = expr.clone();
    let x = Expr::Variable(var.to_string());
    sub_expr = substitute_expr(&sub_expr, &Expr::new_sin(x.clone()), &sin_x_sub);
    sub_expr = substitute_expr(&sub_expr, &Expr::new_cos(x.clone()), &cos_x_sub);
    sub_expr = substitute_expr(&sub_expr, &Expr::new_tan(x.clone()), &tan_x_sub);
    let new_integrand = simplify(Expr::new_mul(sub_expr, dx_sub));
    let integral_in_t = integrate(&new_integrand, "t", None, None);
    if let Expr::Integral { .. } = integral_in_t {
        return None;
    }
    let t_sub_back = Expr::new_tan(Expr::new_div(
        Expr::Variable(var.to_string()),
        Expr::BigInt(BigInt::from(2)),
    ));
    Some(substitute(&integral_in_t, "t", &t_sub_back))
}
/// Computes the limit of an expression as a variable approaches a certain value.
///
/// This is the public entry point, which calls the internal recursive implementation.
/// It handles various limit cases, including direct substitution, indeterminate forms
/// (using L'Hopital's Rule), and limits of rational functions at infinity.
///
/// # Arguments
/// * `expr` - The expression for which to compute the limit.
/// * `var` - The variable that is approaching a value.
/// * `to` - The value that the variable is approaching (e.g., a constant, `Expr::Infinity`, etc.).
///
/// # Returns
/// An `Expr` representing the computed limit.
pub fn limit(expr: &Expr, var: &str, to: &Expr) -> Expr {
    limit_internal(expr, var, to, 0)
}
/// Internal implementation of the limit function with a depth counter to prevent infinite recursion.
///
/// This function applies several strategies:
/// 1.  Checks for base cases (e.g., limit of `e^x` as `x -> oo`).
/// 2.  Attempts direct substitution.
/// 3.  If substitution results in an indeterminate form, it applies transformations or L'Hopital's Rule.
/// 4.  Falls back to specialized logic for rational functions at infinity.
/// 5.  If all else fails, returns an unevaluated `Limit` expression.
pub fn limit_internal(expr: &Expr, var: &str, to: &Expr, depth: u32) -> Expr {
    if depth > 7 {
        return Expr::Limit(
            Arc::new(expr.clone()),
            var.to_string(),
            Arc::new(to.clone()),
        );
    }
    let expr = &simplify(expr.clone());
    match to {
        Expr::Infinity => match expr {
            Expr::Exp(arg) if **arg == Expr::Variable(var.to_string()) => {
                return Expr::Infinity;
            }
            Expr::Log(arg) if **arg == Expr::Variable(var.to_string()) => {
                return Expr::Infinity;
            }
            Expr::ArcTan(arg) if **arg == Expr::Variable(var.to_string()) => {
                return Expr::Constant(std::f64::consts::PI / 2.0);
            }
            Expr::Variable(v) if v == var => return Expr::Infinity,
            _ => {}
        },
        Expr::NegativeInfinity => match expr {
            Expr::Exp(arg) if **arg == Expr::Variable(var.to_string()) => {
                return Expr::BigInt(BigInt::zero());
            }
            Expr::ArcTan(arg) if **arg == Expr::Variable(var.to_string()) => {
                return Expr::Constant(-std::f64::consts::PI / 2.0);
            }
            Expr::Variable(v) if v == var => return Expr::NegativeInfinity,
            _ => {}
        },
        _ => {}
    }
    if !contains_var(expr, var) {
        return expr.clone();
    }
    let val_at_point = simplify(evaluate_at_point(expr, var, to));
    if !matches!(val_at_point, Expr::Infinity | Expr::NegativeInfinity)
        && !contains_var(&val_at_point, var)
    {
        return val_at_point;
    }
    match expr {
        Expr::Div(num, den) => {
            let num_limit = limit_internal(num, var, to, depth + 1);
            let den_limit = limit_internal(den, var, to, depth + 1);
            let is_num_zero = is_zero(&num_limit);
            let is_den_zero = is_zero(&den_limit);
            let is_num_inf = matches!(num_limit, Expr::Infinity | Expr::NegativeInfinity);
            let is_den_inf = matches!(den_limit, Expr::Infinity | Expr::NegativeInfinity);
            if (is_num_zero && is_den_zero) || (is_num_inf && is_den_inf) {
                let d_num = differentiate(num, var);
                let d_den = differentiate(den, var);
                if is_zero(&d_den) {
                    return Expr::Infinity;
                }
                return limit_internal(&Expr::new_div(d_num, d_den), var, to, depth + 1);
            }
        }
        Expr::Mul(a, b) => {
            let a_limit = limit_internal(a, var, to, depth + 1);
            let b_limit = limit_internal(b, var, to, depth + 1);
            if is_zero(&a_limit) && matches!(b_limit, Expr::Infinity | Expr::NegativeInfinity) {
                let new_expr = Expr::new_div(
                    a.clone(),
                    Expr::new_div(Expr::BigInt(BigInt::one()), b.clone()),
                );
                return limit_internal(&new_expr, var, to, depth + 1);
            } else if is_zero(&b_limit)
                && matches!(a_limit, Expr::Infinity | Expr::NegativeInfinity)
            {
                let new_expr = Expr::new_div(
                    b.clone(),
                    Expr::new_div(Expr::BigInt(BigInt::one()), a.clone()),
                );
                return limit_internal(&new_expr, var, to, depth + 1);
            }
        }
        Expr::Power(base, exp) => {
            let base_limit = limit_internal(base, var, to, depth + 1);
            let exp_limit = limit_internal(exp, var, to, depth + 1);
            let is_base_one = is_zero(&simplify(Expr::new_sub(
                base_limit.clone(),
                Expr::BigInt(BigInt::one()),
            )));
            let is_base_zero = is_zero(&base_limit);
            let is_base_inf = matches!(base_limit, Expr::Infinity | Expr::NegativeInfinity);
            let is_exp_inf = matches!(exp_limit, Expr::Infinity | Expr::NegativeInfinity);
            let is_exp_zero = is_zero(&exp_limit);
            if (is_base_one && is_exp_inf)
                || (is_base_zero && is_exp_zero)
                || (is_base_inf && is_exp_zero)
            {
                let log_expr = Expr::new_mul(exp.clone(), Expr::new_log(base.clone()));
                let log_limit = limit_internal(&log_expr, var, to, depth + 1);
                if !contains_var(&log_limit, var) {
                    return Expr::new_exp(log_limit);
                }
            }
        }
        _ => {}
    }
    if let Expr::Infinity | Expr::NegativeInfinity = to {
        if let Expr::Div(num, den) = expr {
            if is_polynomial(num, var) && is_polynomial(den, var) {
                let deg_num = polynomial_degree(num, var);
                let deg_den = polynomial_degree(den, var);
                if deg_num < deg_den {
                    return Expr::BigInt(BigInt::zero());
                } else if deg_num > deg_den {
                    return if matches!(to, Expr::NegativeInfinity) {
                        Expr::NegativeInfinity
                    } else {
                        Expr::Infinity
                    };
                } else {
                    let lead_num = leading_coefficient(num, var);
                    let lead_den = leading_coefficient(den, var);
                    return simplify(Expr::new_div(lead_num, lead_den));
                }
            }
        }
    }
    if !contains_var(&val_at_point, var) {
        val_at_point
    } else {
        Expr::Limit(
            Arc::new(expr.clone()),
            var.to_string(),
            Arc::new(to.clone()),
        )
    }
}

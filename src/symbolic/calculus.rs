//! # Symbolic Calculus Engine
//!
//! This module provides a comprehensive suite of symbolic calculus operations.
//! It includes functions for differentiation, indefinite and definite integration,
//! path integrals, limits, and series expansions. The integration capabilities
//! are supported by a multi-strategy approach including rule-based integration,
//! u-substitution, integration by parts, and more.

use crate::symbolic::core::{Expr, PathType};
use crate::symbolic::polynomial::{is_polynomial, leading_coefficient, polynomial_degree};
use crate::symbolic::simplify::is_zero;
use crate::symbolic::simplify::simplify;
use crate::symbolic::solve::solve;
use num_bigint::BigInt;
use num_traits::{One, Zero};

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
    match expr {
        Expr::Variable(name) if name == var => replacement.clone(),
        Expr::Add(a, b) => Expr::Add(
            Box::new(substitute(a, var, replacement)),
            Box::new(substitute(b, var, replacement)),
        ),
        Expr::Sub(a, b) => Expr::Sub(
            Box::new(substitute(a, var, replacement)),
            Box::new(substitute(b, var, replacement)),
        ),
        Expr::Mul(a, b) => Expr::Mul(
            Box::new(substitute(a, var, replacement)),
            Box::new(substitute(b, var, replacement)),
        ),
        Expr::Div(a, b) => Expr::Div(
            Box::new(substitute(a, var, replacement)),
            Box::new(substitute(b, var, replacement)),
        ),
        Expr::Power(base, exp) => Expr::Power(
            Box::new(substitute(base, var, replacement)),
            Box::new(substitute(exp, var, replacement)),
        ),
        Expr::Sin(arg) => Expr::Sin(Box::new(substitute(arg, var, replacement))),
        Expr::Cos(arg) => Expr::Cos(Box::new(substitute(arg, var, replacement))),
        Expr::Tan(arg) => Expr::Tan(Box::new(substitute(arg, var, replacement))),
        Expr::Exp(arg) => Expr::Exp(Box::new(substitute(arg, var, replacement))),
        Expr::Log(arg) => Expr::Log(Box::new(substitute(arg, var, replacement))),
        Expr::Integral {
            integrand,
            var: int_var,
            lower_bound,
            upper_bound,
        } => Expr::Integral {
            integrand: Box::new(substitute(integrand, var, replacement)),
            var: int_var.clone(),
            lower_bound: Box::new(substitute(lower_bound, var, replacement)),
            upper_bound: Box::new(substitute(upper_bound, var, replacement)),
        },
        Expr::Sum {
            body,
            var: sum_var,
            from,
            to,
        } => {
            let new_from = substitute(from, var, replacement);
            let new_to = substitute(to, var, replacement);
            // If the substitution variable is the same as the summation index, don't substitute into the body.
            let new_body = if let Expr::Variable(v) = &**sum_var {
                if v == var {
                    body.clone()
                } else {
                    Box::new(substitute(body, var, replacement))
                }
            } else {
                Box::new(substitute(body, var, replacement))
            };
            Expr::Sum {
                body: Box::new(*new_body),
                var: sum_var.clone(),
                from: Box::new(new_from),
                to: Box::new(new_to),
            }
        }
        _ => expr.clone(),
    }
}

pub(crate) fn get_real_imag_parts(expr: &Expr) -> (Expr, Expr) {
    match simplify(expr.clone()) {
        Expr::Complex(re, im) => (*re, *im),
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
    match expr {
        Expr::Constant(_) | Expr::BigInt(_) | Expr::Rational(_) | Expr::Pi | Expr::E => {
            Expr::BigInt(BigInt::zero())
        }
        Expr::Variable(name) if name == var => Expr::BigInt(BigInt::one()),
        Expr::Variable(_) => Expr::BigInt(BigInt::zero()),
        Expr::Add(a, b) => simplify(Expr::Add(
            Box::new(differentiate(a, var)),
            Box::new(differentiate(b, var)),
        )),
        Expr::Sub(a, b) => simplify(Expr::Sub(
            Box::new(differentiate(a, var)),
            Box::new(differentiate(b, var)),
        )),
        Expr::Mul(a, b) => simplify(Expr::Add(
            Box::new(Expr::Mul(Box::new(differentiate(a, var)), b.clone())),
            Box::new(Expr::Mul(a.clone(), Box::new(differentiate(b, var)))),
        )),
        Expr::Div(a, b) => simplify(Expr::Div(
            Box::new(Expr::Sub(
                Box::new(Expr::Mul(Box::new(differentiate(a, var)), b.clone())),
                Box::new(Expr::Mul(a.clone(), Box::new(differentiate(b, var)))),
            )),
            Box::new(Expr::Power(
                b.clone(),
                Box::new(Expr::BigInt(BigInt::from(2))),
            )),
        )),
        Expr::Power(base, exp) => {
            let d_base = differentiate(base, var);
            let d_exp = differentiate(exp, var);
            // Use chain rule: d/dx(f(x)^g(x)) = f(x)^g(x) * [g'(x)ln(f(x)) + g(x) * f'(x)/f(x)]
            let term1 = Expr::Mul(Box::new(d_exp), Box::new(Expr::Log(base.clone())));
            let term2 = Expr::Mul(
                Box::new(*exp.clone()),
                Box::new(Expr::Div(Box::new(d_base), base.clone())),
            );
            simplify(Expr::Mul(
                Box::new(expr.clone()),
                Box::new(Expr::Add(Box::new(term1), Box::new(term2))),
            ))
        }
        Expr::Sin(arg) => simplify(Expr::Mul(
            Box::new(Expr::Cos(arg.clone())),
            Box::new(differentiate(arg, var)),
        )),
        Expr::Cos(arg) => simplify(Expr::Mul(
            Box::new(Expr::Neg(Box::new(Expr::Sin(arg.clone())))),
            Box::new(differentiate(arg, var)),
        )),
        Expr::Tan(arg) => simplify(Expr::Mul(
            Box::new(Expr::Power(
                Box::new(Expr::Sec(arg.clone())),
                Box::new(Expr::BigInt(BigInt::from(2))),
            )),
            Box::new(differentiate(arg, var)),
        )),

        Expr::Sec(arg) => simplify(Expr::Mul(
            Box::new(Expr::Sec(arg.clone())),
            Box::new(Expr::Mul(
                Box::new(Expr::Tan(arg.clone())),
                Box::new(differentiate(arg, var)),
            )),
        )),
        Expr::Csc(arg) => simplify(Expr::Mul(
            Box::new(Expr::Neg(Box::new(Expr::Csc(arg.clone())))),
            Box::new(Expr::Mul(
                Box::new(Expr::Cot(arg.clone())),
                Box::new(differentiate(arg, var)),
            )),
        )),
        Expr::Cot(arg) => simplify(Expr::Mul(
            Box::new(Expr::Neg(Box::new(Expr::Power(
                Box::new(Expr::Csc(arg.clone())),
                Box::new(Expr::BigInt(BigInt::from(2))),
            )))),
            Box::new(differentiate(arg, var)),
        )),
        Expr::Sinh(arg) => simplify(Expr::Mul(
            Box::new(Expr::Cosh(arg.clone())),
            Box::new(differentiate(arg, var)),
        )),
        Expr::Cosh(arg) => simplify(Expr::Mul(
            Box::new(Expr::Sinh(arg.clone())),
            Box::new(differentiate(arg, var)),
        )),
        Expr::Tanh(arg) => simplify(Expr::Mul(
            Box::new(Expr::Power(
                Box::new(Expr::Sech(arg.clone())),
                Box::new(Expr::BigInt(BigInt::from(2))),
            )),
            Box::new(differentiate(arg, var)),
        )),
        Expr::Exp(arg) => simplify(Expr::Mul(
            Box::new(Expr::Exp(arg.clone())),
            Box::new(differentiate(arg, var)),
        )),
        Expr::Log(arg) => simplify(Expr::Div(Box::new(differentiate(arg, var)), arg.clone())),
        Expr::ArcCot(arg) => simplify(Expr::Neg(Box::new(Expr::Div(
            Box::new(differentiate(arg, var)),
            Box::new(Expr::Add(
                Box::new(Expr::BigInt(BigInt::one())),
                Box::new(Expr::Power(
                    arg.clone(),
                    Box::new(Expr::BigInt(BigInt::from(2))),
                )),
            )),
        )))),
        Expr::ArcSec(arg) => simplify(Expr::Div(
            Box::new(differentiate(arg, var)),
            Box::new(Expr::Mul(
                Box::new(Expr::Abs(arg.clone())),
                Box::new(Expr::Sqrt(Box::new(Expr::Sub(
                    Box::new(Expr::Power(
                        arg.clone(),
                        Box::new(Expr::BigInt(BigInt::from(2))),
                    )),
                    Box::new(Expr::BigInt(BigInt::one())),
                )))),
            )),
        )),
        Expr::ArcCsc(arg) => simplify(Expr::Neg(Box::new(Expr::Div(
            Box::new(differentiate(arg, var)),
            Box::new(Expr::Mul(
                Box::new(Expr::Abs(arg.clone())),
                Box::new(Expr::Sqrt(Box::new(Expr::Sub(
                    Box::new(Expr::Power(
                        arg.clone(),
                        Box::new(Expr::BigInt(BigInt::from(2))),
                    )),
                    Box::new(Expr::BigInt(BigInt::one())),
                )))),
            )),
        )))),
        Expr::Coth(arg) => simplify(Expr::Neg(Box::new(Expr::Power(
            Box::new(Expr::Csch(arg.clone())),
            Box::new(Expr::BigInt(BigInt::from(2))),
        )))),
        Expr::Sech(arg) => simplify(Expr::Neg(Box::new(Expr::Mul(
            Box::new(Expr::Sech(arg.clone())),
            Box::new(Expr::Tanh(arg.clone())),
        )))),
        Expr::Csch(arg) => simplify(Expr::Neg(Box::new(Expr::Mul(
            Box::new(Expr::Csch(arg.clone())),
            Box::new(Expr::Coth(arg.clone())),
        )))),
        Expr::ArcSinh(arg) => simplify(Expr::Div(
            Box::new(differentiate(arg, var)),
            Box::new(Expr::Sqrt(Box::new(Expr::Add(
                Box::new(Expr::Power(
                    arg.clone(),
                    Box::new(Expr::BigInt(BigInt::from(2))),
                )),
                Box::new(Expr::BigInt(BigInt::one())),
            )))),
        )),
        Expr::ArcCosh(arg) => simplify(Expr::Div(
            Box::new(differentiate(arg, var)),
            Box::new(Expr::Sqrt(Box::new(Expr::Sub(
                Box::new(Expr::Power(
                    arg.clone(),
                    Box::new(Expr::BigInt(BigInt::from(2))),
                )),
                Box::new(Expr::BigInt(BigInt::one())),
            )))),
        )),
        Expr::ArcTanh(arg) => simplify(Expr::Div(
            Box::new(differentiate(arg, var)),
            Box::new(Expr::Sub(
                Box::new(Expr::BigInt(BigInt::one())),
                Box::new(Expr::Power(
                    arg.clone(),
                    Box::new(Expr::BigInt(BigInt::from(2))),
                )),
            )),
        )),
        Expr::ArcCoth(arg) => simplify(Expr::Div(
            Box::new(differentiate(arg, var)),
            Box::new(Expr::Sub(
                Box::new(Expr::BigInt(BigInt::one())),
                Box::new(Expr::Power(
                    arg.clone(),
                    Box::new(Expr::BigInt(BigInt::from(2))),
                )),
            )),
        )),
        Expr::ArcSech(arg) => simplify(Expr::Neg(Box::new(Expr::Div(
            Box::new(differentiate(arg, var)),
            Box::new(Expr::Mul(
                arg.clone(),
                Box::new(Expr::Sqrt(Box::new(Expr::Sub(
                    Box::new(Expr::BigInt(BigInt::one())),
                    Box::new(Expr::Power(
                        arg.clone(),
                        Box::new(Expr::BigInt(BigInt::from(2))),
                    )),
                )))),
            )),
        )))),
        Expr::ArcCsch(arg) => simplify(Expr::Neg(Box::new(Expr::Div(
            Box::new(differentiate(arg, var)),
            Box::new(Expr::Mul(
                Box::new(Expr::Abs(arg.clone())),
                Box::new(Expr::Sqrt(Box::new(Expr::Add(
                    Box::new(Expr::BigInt(BigInt::one())),
                    Box::new(Expr::Power(
                        arg.clone(),
                        Box::new(Expr::BigInt(BigInt::from(2))),
                    )),
                )))),
            )),
        )))),
        Expr::Integral {
            integrand,
            var: int_var,
            ..
        } => {
            if *int_var.clone() == Expr::Variable(var.to_string()) {
                *integrand.clone()
            } else {
                Expr::Derivative(Box::new(expr.clone()), var.to_string())
            }
        }
        Expr::Sum {
            body,
            var: sum_var,
            from,
            to,
        } => {
            // Differentiate the body of the sum
            let diff_body = differentiate(body, var);
            // Return a new sum with the differentiated body
            Expr::Sum {
                body: Box::new(diff_body),
                var: sum_var.clone(),
                from: from.clone(),
                to: to.clone(),
            }
        }
        _ => todo!(),
    }
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
    // If bounds are provided, perform definite integration.
    if let (Some(lower), Some(upper)) = (lower_bound, upper_bound) {
        return definite_integrate(expr, var, lower, upper);
    }

    // Otherwise, perform indefinite integration.
    let simplified_expr = simplify(expr.clone());

    // Strategy 1: Rule-based matching
    if let Some(result) = integrate_by_rules(&simplified_expr, var) {
        return simplify(result);
    }

    // Strategy 2: U-Substitution
    if let Some(result) = u_substitution(&simplified_expr, var) {
        return simplify(result);
    }

    // Strategy 3: Integration by Parts
    if let Some(result) = integrate_by_parts_master(&simplified_expr, var, 0) {
        return simplify(result);
    }

    // Strategy 4: Partial Fractions
    if let Some(result) = integrate_by_partial_fractions(&simplified_expr, var) {
        return simplify(result);
    }

    // Strategy 5: Trigonometric Substitution
    if let Some(result) = trig_substitution(&simplified_expr, var) {
        return simplify(result);
    }

    // Strategy 6: Tangent Half-Angle Substitution
    if let Some(result) = tangent_half_angle_substitution(&simplified_expr, var) {
        return simplify(result);
    }

    // Fallback to basic integration patterns or return unevaluated integral
    let basic_result = integrate_basic(&simplified_expr, var);
    if let Expr::Integral { .. } = basic_result {
        // Return unevaluated integral if basic integration also fails
        Expr::Integral {
            integrand: Box::new(expr.clone()),
            var: Box::new(Expr::Variable(var.to_string())),
            lower_bound: Box::new(Expr::Variable("a".to_string())),
            upper_bound: Box::new(Expr::Variable("b".to_string())),
        }
    } else {
        simplify(basic_result)
    }
}

pub(crate) fn integrate_basic(expr: &Expr, var: &str) -> Expr {
    match expr {
        Expr::Constant(c) => Expr::Mul(
            Box::new(Expr::Constant(*c)),
            Box::new(Expr::Variable(var.to_string())),
        ),
        Expr::BigInt(i) => Expr::Mul(
            Box::new(Expr::BigInt(i.clone())),
            Box::new(Expr::Variable(var.to_string())),
        ),
        Expr::Rational(r) => Expr::Mul(
            Box::new(Expr::Rational(r.clone())),
            Box::new(Expr::Variable(var.to_string())),
        ),
        Expr::Variable(name) if name == var => Expr::Div(
            Box::new(Expr::Power(
                Box::new(Expr::Variable(var.to_string())),
                Box::new(Expr::BigInt(BigInt::from(2))),
            )),
            Box::new(Expr::BigInt(BigInt::from(2))),
        ),
        Expr::Add(a, b) => simplify(Expr::Add(
            Box::new(integrate(a, var, None, None)),
            Box::new(integrate(b, var, None, None)),
        )),
        Expr::Sub(a, b) => simplify(Expr::Sub(
            Box::new(integrate(a, var, None, None)),
            Box::new(integrate(b, var, None, None)),
        )),
        Expr::Power(base, exp) => {
            if let (Expr::Variable(name), Expr::Constant(n)) = (&**base, &**exp) {
                if name == var {
                    if (*n + 1.0).abs() < 1e-9 {
                        return Expr::Log(Box::new(Expr::Abs(Box::new(Expr::Variable(
                            var.to_string(),
                        )))));
                    }
                    return Expr::Div(
                        Box::new(Expr::Power(
                            Box::new(Expr::Variable(var.to_string())),
                            Box::new(Expr::Constant(n + 1.0)),
                        )),
                        Box::new(Expr::Constant(n + 1.0)),
                    );
                }
            }
            Expr::Integral {
                integrand: Box::new(expr.clone()),
                var: Box::new(Expr::Variable(var.to_string())),
                lower_bound: Box::new(Expr::Variable("a".to_string())),
                upper_bound: Box::new(Expr::Variable("b".to_string())),
            }
        }
        Expr::Exp(arg) => {
            if let Expr::Variable(name) = &**arg {
                if name == var {
                    return Expr::Exp(Box::new(Expr::Variable(var.to_string())));
                }
            }
            Expr::Integral {
                integrand: Box::new(expr.clone()),
                var: Box::new(Expr::Variable(var.to_string())),
                lower_bound: Box::new(Expr::Variable("a".to_string())),
                upper_bound: Box::new(Expr::Variable("b".to_string())),
            }
        }
        _ => Expr::Integral {
            integrand: Box::new(expr.clone()),
            var: Box::new(Expr::Variable(var.to_string())),
            lower_bound: Box::new(Expr::Variable("a".to_string())),
            upper_bound: Box::new(Expr::Variable("b".to_string())),
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
    if depth > 2 {
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
        let uv = Expr::Mul(Box::new(*u.clone()), Box::new(v.clone()));
        let v_du = Expr::Mul(Box::new(v), Box::new(du_dx));
        let integral_v_du = integrate(&v_du, var, None, None);
        return Some(simplify(Expr::Sub(Box::new(uv), Box::new(integral_v_du))));
    }
    None
}

// Helper function to substitute an expression with another expression.
pub(crate) fn substitute_expr(expr: &Expr, to_replace: &Expr, replacement: &Expr) -> Expr {
    if expr == to_replace {
        return replacement.clone();
    }
    // Recursively traverse the expression tree.
    match expr {
        Expr::Add(a, b) => Expr::Add(
            Box::new(substitute_expr(a, to_replace, replacement)),
            Box::new(substitute_expr(b, to_replace, replacement)),
        ),
        Expr::Sub(a, b) => Expr::Sub(
            Box::new(substitute_expr(a, to_replace, replacement)),
            Box::new(substitute_expr(b, to_replace, replacement)),
        ),
        Expr::Mul(a, b) => Expr::Mul(
            Box::new(substitute_expr(a, to_replace, replacement)),
            Box::new(substitute_expr(b, to_replace, replacement)),
        ),
        Expr::Div(a, b) => Expr::Div(
            Box::new(substitute_expr(a, to_replace, replacement)),
            Box::new(substitute_expr(b, to_replace, replacement)),
        ),
        Expr::Power(base, exp) => Expr::Power(
            Box::new(substitute_expr(base, to_replace, replacement)),
            Box::new(substitute_expr(exp, to_replace, replacement)),
        ),
        Expr::Sin(arg) => Expr::Sin(Box::new(substitute_expr(arg, to_replace, replacement))),
        Expr::Cos(arg) => Expr::Cos(Box::new(substitute_expr(arg, to_replace, replacement))),
        Expr::Tan(arg) => Expr::Tan(Box::new(substitute_expr(arg, to_replace, replacement))),
        Expr::Sec(arg) => Expr::Sec(Box::new(substitute_expr(arg, to_replace, replacement))),
        Expr::Csc(arg) => Expr::Csc(Box::new(substitute_expr(arg, to_replace, replacement))),
        Expr::Cot(arg) => Expr::Cot(Box::new(substitute_expr(arg, to_replace, replacement))),
        Expr::Sinh(arg) => Expr::Sinh(Box::new(substitute_expr(arg, to_replace, replacement))),
        Expr::Cosh(arg) => Expr::Cosh(Box::new(substitute_expr(arg, to_replace, replacement))),
        Expr::Tanh(arg) => Expr::Tanh(Box::new(substitute_expr(arg, to_replace, replacement))),
        Expr::Exp(arg) => Expr::Exp(Box::new(substitute_expr(arg, to_replace, replacement))),
        Expr::Log(arg) => Expr::Log(Box::new(substitute_expr(arg, to_replace, replacement))),
        Expr::Complex(re, im) => Expr::Complex(
            Box::new(substitute_expr(re, to_replace, replacement)),
            Box::new(substitute_expr(im, to_replace, replacement)),
        ),
        Expr::Abs(arg) => Expr::Abs(Box::new(substitute_expr(arg, to_replace, replacement))),
        Expr::Neg(arg) => Expr::Neg(Box::new(substitute_expr(arg, to_replace, replacement))),
        Expr::ArcSin(arg) => Expr::ArcSin(Box::new(substitute_expr(arg, to_replace, replacement))),
        Expr::ArcCos(arg) => Expr::ArcCos(Box::new(substitute_expr(arg, to_replace, replacement))),
        Expr::ArcTan(arg) => Expr::ArcTan(Box::new(substitute_expr(arg, to_replace, replacement))),
        // Base cases: if no match, the expression does not contain `to_replace` in its branches.
        _ => expr.clone(),
    }
}

// Helper to check if an expression contains a specific variable.
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

// Gathers potential candidates for u-substitution.
pub(crate) fn get_u_candidates(expr: &Expr, candidates: &mut Vec<Expr>) {
    match expr {
        Expr::Add(a, b) | Expr::Sub(a, b) | Expr::Mul(a, b) | Expr::Div(a, b) => {
            get_u_candidates(a, candidates);
            get_u_candidates(b, candidates);
        }
        Expr::Power(b, e) => {
            candidates.push(*b.clone());
            get_u_candidates(b, candidates);
            get_u_candidates(e, candidates);
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
            candidates.push(*a.clone());
            get_u_candidates(a, candidates);
        }
        _ => {}
    }
}

// The main u-substitution strategy.
pub(crate) fn u_substitution(expr: &Expr, var: &str) -> Option<Expr> {
    // Strategy 1: Direct check for the g'(x)/g(x) pattern.
    if let Expr::Div(num, den) = expr {
        let den_prime = differentiate(den, var);
        let c = simplify(Expr::Div(Box::new(*num.clone()), Box::new(den_prime)));
        // If num is a constant multiple of den_prime
        if !contains_var(&c, var) {
            let log_den = Expr::Log(Box::new(Expr::Abs(den.clone())));
            return Some(simplify(Expr::Mul(Box::new(c), Box::new(log_den))));
        }
    }

    // Strategy 2: General substitution based on candidates.
    let mut candidates = Vec::new();
    get_u_candidates(expr, &mut candidates);
    // Also consider the expression itself as a candidate
    candidates.push(expr.clone());

    for u in candidates {
        if let Expr::Variable(_) = u {
            continue; // Skip simple variables
        }

        let du_dx = differentiate(&u, var);
        if is_zero(&du_dx) {
            continue;
        }

        // Try to rewrite the integrand in terms of u.
        // new_integrand = expr / (du/dx)
        let new_integrand_x = simplify(Expr::Div(Box::new(expr.clone()), Box::new(du_dx)));

        // Substitute x with a function of u. This is the hard part.
        // For now, we do a simpler check: substitute u with a temp variable `t`
        // in the new_integrand_x. If the original var `x` disappears, we have a success.
        let temp_var = "t";
        let temp_expr = Expr::Variable(temp_var.to_string());
        let substituted = substitute_expr(&new_integrand_x, &u, &temp_expr);

        if !contains_var(&substituted, var) {
            // We have a valid substitution. The new integrand is `substituted`.
            let integral_in_t = integrate(&substituted, temp_var, None, None);

            // If integration was successful (not an unevaluated integral)
            if !matches!(integral_in_t, Expr::Integral { .. }) {
                // Substitute back: t -> u
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
        if let (Expr::Variable(v), Expr::Constant(2.0)) = (&**x, *two.clone()) {
            if v == var && *a_val > 0.0 {
                let a = Expr::Constant(a_val.sqrt());
                let theta = Expr::Variable("theta".to_string());
                // Substitution: x = a*tan(theta), dx = a*sec^2(theta)d(theta)
                let x_sub = Expr::Mul(
                    Box::new(a.clone()),
                    Box::new(Expr::Tan(Box::new(theta.clone()))),
                );
                let dx_dtheta = differentiate(&x_sub, "theta");

                let new_integrand = simplify(Expr::Mul(
                    Box::new(substitute(expr, var, &x_sub)),
                    Box::new(dx_dtheta),
                ));
                let integral_theta = integrate(&new_integrand, "theta", None, None);

                // Substitute back: theta = atan(x/a)
                let theta_sub = Expr::ArcTan(Box::new(Expr::Div(
                    Box::new(Expr::Variable(var.to_string())),
                    Box::new(a),
                )));
                return Some(substitute(&integral_theta, "theta", &theta_sub));
            }
        }
    }
    None
}

pub(crate) fn trig_substitution(expr: &Expr, var: &str) -> Option<Expr> {
    // This function handles integrals involving expressions of the form:
    // sqrt(a^2 - x^2), sqrt(a^2 + x^2), and sqrt(x^2 - a^2).
    if let Expr::Sqrt(arg) = expr {
        // Case 1: sqrt(a^2 - x^2)  =>  x = a*sin(theta)
        if let Expr::Sub(a_sq, x_sq) = &**arg {
            if let (Expr::Constant(a_val), Expr::Power(x, two)) = (&**a_sq, &**x_sq) {
                if let (Expr::Variable(v), Expr::Constant(2.0)) = (&**x, *two.clone()) {
                    if v == var && *a_val > 0.0 {
                        let a = Expr::Constant(a_val.sqrt());
                        let theta = Expr::Variable("theta".to_string());
                        // Substitution: x = a*sin(theta), dx = a*cos(theta)d(theta)
                        let x_sub = Expr::Mul(
                            Box::new(a.clone()),
                            Box::new(Expr::Sin(Box::new(theta.clone()))),
                        );
                        let dx_dtheta = differentiate(&x_sub, "theta");

                        let new_integrand = simplify(Expr::Mul(
                            Box::new(substitute(expr, var, &x_sub)),
                            Box::new(dx_dtheta),
                        ));
                        let integral_theta = integrate(&new_integrand, "theta", None, None);

                        // Substitute back: theta = asin(x/a)
                        let theta_sub = Expr::ArcSin(Box::new(Expr::Div(
                            Box::new(Expr::Variable(var.to_string())),
                            Box::new(a),
                        )));
                        return Some(substitute(&integral_theta, "theta", &theta_sub));
                    }
                }
            }
        }

        // Case 2: sqrt(a^2 + x^2)  =>  x = a*tan(theta)
        if let Expr::Add(part1, part2) = &**arg {
            // The order could be x^2 + a^2, so we check both combinations.
            if let Some(result) = handle_trig_sub_sum(part1, part2, expr, var) {
                return Some(result);
            }
            if let Some(result) = handle_trig_sub_sum(part2, part1, expr, var) {
                return Some(result);
            }
        }

        // Case 3: sqrt(x^2 - a^2)  =>  x = a*sec(theta)
        if let Expr::Sub(x_sq, a_sq) = &**arg {
            if let (Expr::Power(x, two), Expr::Constant(a_val)) = (&**x_sq, &**a_sq) {
                if let (Expr::Variable(v), Expr::Constant(2.0)) = (&**x, *two.clone()) {
                    if v == var && *a_val > 0.0 {
                        let a = Expr::Constant(a_val.sqrt());
                        let theta = Expr::Variable("theta".to_string());
                        // Substitution: x = a*sec(theta), dx = a*sec(theta)tan(theta)d(theta)
                        let x_sub = Expr::Mul(
                            Box::new(a.clone()),
                            Box::new(Expr::Sec(Box::new(theta.clone()))),
                        );
                        let dx_dtheta = differentiate(&x_sub, "theta");

                        let new_integrand = simplify(Expr::Mul(
                            Box::new(substitute(expr, var, &x_sub)),
                            Box::new(dx_dtheta),
                        ));
                        let integral_theta = integrate(&new_integrand, "theta", None, None);

                        // Substitute back: theta = asec(x/a)
                        let theta_sub = Expr::ArcSec(Box::new(Expr::Div(
                            Box::new(Expr::Variable(var.to_string())),
                            Box::new(a),
                        )));
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
    } // Integration failed
    let upper_eval = evaluate_at_point(&antiderivative, var, upper_bound);
    let lower_eval = evaluate_at_point(&antiderivative, var, lower_bound);
    simplify(Expr::Sub(Box::new(upper_eval), Box::new(lower_eval)))
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
    let z_replacement = Expr::Complex(
        Box::new(Expr::Variable("x".to_string())),
        Box::new(Expr::Variable("y".to_string())),
    );
    let f_xy = substitute(expr, var, &z_replacement);
    let (u, v) = get_real_imag_parts(&f_xy);
    let du_dx = differentiate(&u, "x");
    let du_dy = differentiate(&u, "y");
    let dv_dx = differentiate(&v, "x");
    let dv_dy = differentiate(&v, "y");
    let cr1 = simplify(Expr::Sub(Box::new(du_dx), Box::new(dv_dy)));
    let cr2 = simplify(Expr::Add(Box::new(du_dy), Box::new(dv_dx)));
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

// Determines the order of a pole for a given expression.
pub(crate) fn find_pole_order(expr: &Expr, var: &str, pole: &Expr) -> usize {
    let mut order = 1;
    loop {
        let term = Expr::Power(
            Box::new(Expr::Sub(
                Box::new(Expr::Variable(var.to_string())),
                Box::new(pole.clone()),
            )),
            Box::new(Expr::BigInt(BigInt::from(order))),
        );
        let new_expr = simplify(Expr::Mul(Box::new(expr.clone()), Box::new(term)));
        let val_at_pole = simplify(evaluate_at_point(&new_expr, var, pole));
        if let Expr::Constant(c) = val_at_pole {
            if c.is_finite() && c.abs() > 1e-9 {
                return order;
            }
        }
        order += 1;
        if order > 10 {
            // Safety break
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
    // Formula for simple pole: Res(f, c) = g(c) / h'(c) where f = g/h
    if let Expr::Div(num, den) = expr {
        let den_prime = differentiate(den, var);
        let num_at_pole = evaluate_at_point(num, var, pole);
        let den_prime_at_pole = evaluate_at_point(&den_prime, var, pole);
        if !is_zero(&simplify(den_prime_at_pole.clone())) {
            return simplify(Expr::Div(
                Box::new(num_at_pole),
                Box::new(den_prime_at_pole),
            ));
        }
    }
    // Fallback for poles of order m > 1
    let m = find_pole_order(expr, var, pole);
    let m_minus_1_factorial = factorial(m - 1);
    let term = Expr::Power(
        Box::new(Expr::Sub(
            Box::new(Expr::Variable(var.to_string())),
            Box::new(pole.clone()),
        )),
        Box::new(Expr::BigInt(BigInt::from(m))),
    );
    let g_z = simplify(Expr::Mul(Box::new(expr.clone()), Box::new(term)));
    let mut g_m_minus_1 = g_z;
    for _ in 0..(m - 1) {
        g_m_minus_1 = differentiate(&g_m_minus_1, var);
    }
    let limit = evaluate_at_point(&g_m_minus_1, var, pole);
    simplify(Expr::Div(
        Box::new(limit),
        Box::new(Expr::Constant(m_minus_1_factorial)),
    ))
}

// /// Checks if a given complex point is inside a circular contour.
// pub fn is_inside_contour(point: &Expr, contour: &Expr) -> bool {
//     if let (Expr::Path(path_type, center, radius), Expr::Complex(re, im)) = (contour, point) {
//         if let PathType::Circle = path_type {
//             if let (Expr::Complex(center_re, center_im), Expr::Constant(r)) = (&**center, &**radius)
//             {
//                 let dist_sq = Expr::Add(
//                     Box::new(Expr::Power(
//                         Box::new(Expr::Sub(re.clone(), center_re.clone())),
//                         Box::new(Expr::BigInt(BigInt::from(2))),
//                     )),
//                     Box::new(Expr::Power(
//                         Box::new(Expr::Sub(im.clone(), center_im.clone())),
//                         Box::new(Expr::BigInt(BigInt::from(2))),
//                     )),
//                 );
//                 if let Expr::Constant(d2) = simplify(dist_sq) {
//                     return d2 < r * r;
//                 }
//             }
//         }
//     }
//     false
// }/// Checks if a given complex point is inside a specified circular contour.
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
    // Merged: Replaced 'path_type' with 'PathType::Circle'
    if let (Expr::Path(PathType::Circle, center, radius), Expr::Complex(re, im)) = (contour, point)
    {
        // This was the original inner 'if let' block
        if let (Expr::Complex(center_re, center_im), Expr::Constant(r)) = (&**center, &**radius) {
            let dist_sq = Expr::Add(
                Box::new(Expr::Power(
                    Box::new(Expr::Sub(re.clone(), center_re.clone())),
                    Box::new(Expr::BigInt(BigInt::from(2))),
                )),
                Box::new(Expr::Power(
                    Box::new(Expr::Sub(im.clone(), center_im.clone())),
                    Box::new(Expr::BigInt(BigInt::from(2))),
                )),
            );
            if let Expr::Constant(d2) = simplify(dist_sq) {
                return d2 < r * r;
            }
        }
    }
    // Note: You must ensure 'BigInt' is in scope (e.g., via 'use' statement)
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
                        sum_of_residues = Expr::Add(Box::new(sum_of_residues), Box::new(residue));
                    }
                }
                let two_pi_i = Expr::Mul(
                    Box::new(Expr::Constant(2.0 * std::f64::consts::PI)),
                    Box::new(Expr::Complex(
                        Box::new(Expr::BigInt(BigInt::zero())),
                        Box::new(Expr::BigInt(BigInt::one())),
                    )),
                );
                simplify(Expr::Mul(Box::new(two_pi_i), Box::new(sum_of_residues)))
            }
            PathType::Line => {
                let (z0, z1) = (&**param1, &**param2);
                let dz_dt = simplify(Expr::Sub(Box::new(z1.clone()), Box::new(z0.clone())));
                let t_var = Expr::Variable("t".to_string());
                let z_t = simplify(Expr::Add(
                    Box::new(z0.clone()),
                    Box::new(Expr::Mul(Box::new(t_var.clone()), Box::new(dz_dt.clone()))),
                ));
                let integrand_t = simplify(Expr::Mul(
                    Box::new(substitute(expr, var, &z_t)),
                    Box::new(dz_dt),
                ));
                definite_integrate(
                    &integrand_t,
                    "t",
                    &Expr::BigInt(BigInt::zero()),
                    &Expr::BigInt(BigInt::one()),
                )
            }
            PathType::Rectangle => {
                let (z_bl, z_tr) = (&**param1, &**param2);
                let z_br = Expr::Complex(Box::new(z_tr.re()), Box::new(z_bl.im()));
                let z_tl = Expr::Complex(Box::new(z_bl.re()), Box::new(z_tr.im()));
                let i1 = path_integrate(
                    expr,
                    var,
                    &Expr::Path(
                        PathType::Line,
                        Box::new(z_bl.clone()),
                        Box::new(z_br.clone()),
                    ),
                );
                let i2 = path_integrate(
                    expr,
                    var,
                    &Expr::Path(PathType::Line, Box::new(z_br), Box::new(z_tr.clone())),
                );
                let i3 = path_integrate(
                    expr,
                    var,
                    &Expr::Path(
                        PathType::Line,
                        Box::new(z_tr.clone()),
                        Box::new(z_tl.clone()),
                    ),
                );
                let i4 = path_integrate(
                    expr,
                    var,
                    &Expr::Path(PathType::Line, Box::new(z_tl), Box::new(z_bl.clone())),
                );
                simplify(Expr::Add(
                    Box::new(i1),
                    Box::new(Expr::Add(
                        Box::new(i2),
                        Box::new(Expr::Add(Box::new(i3), Box::new(i4))),
                    )),
                ))
            }
        },
        _ => Expr::Integral {
            integrand: Box::new(expr.clone()),
            var: Box::new(Expr::Variable(var.to_string())),
            lower_bound: Box::new(Expr::Variable("C_lower".to_string())),
            upper_bound: Box::new(Expr::Variable("C_upper".to_string())),
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
    // Helper to extract the imaginary part of a simplified expression.
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
            Expr::Constant(_val) => Some(0.0), // Real constants have zero imaginary part.
            expr if is_zero(&expr) => Some(0.0),
            _ => None, // Cannot determine imaginary part for other types.
        }
    }

    let poles = find_poles(expr, var);
    let mut sum_of_residues_in_uhp = Expr::BigInt(BigInt::zero());

    for pole in poles {
        // Check if the pole is in the upper half-plane.
        if let Some(im_val) = get_imag_part(&pole) {
            if im_val > 1e-9 {
                // Use a small epsilon to avoid floating point issues near the real axis.
                let residue = calculate_residue(expr, var, &pole);
                // If residue calculation fails, it might return a non-simplified expression.
                // For now, we add it directly. A more robust implementation might handle errors.
                sum_of_residues_in_uhp = simplify(Expr::Add(
                    Box::new(sum_of_residues_in_uhp.clone()),
                    Box::new(residue),
                ));
            }
        }
    }

    // The result is 2 * pi * i * sum_of_residues
    let two_pi_i = Expr::Mul(
        Box::new(Expr::Constant(2.0 * std::f64::consts::PI)),
        Box::new(Expr::Complex(
            Box::new(Expr::BigInt(BigInt::zero())),
            Box::new(Expr::BigInt(BigInt::one())),
        )),
    );

    simplify(Expr::Mul(
        Box::new(two_pi_i),
        Box::new(sum_of_residues_in_uhp),
    ))
}

// A more comprehensive rule-based integrator
pub(crate) fn integrate_by_rules(expr: &Expr, var: &str) -> Option<Expr> {
    match expr {
        // Basic rules
        Expr::Constant(c) => Some(Expr::Mul(
            Box::new(Expr::Constant(*c)),
            Box::new(Expr::Variable(var.to_string())),
        )),
        Expr::BigInt(i) => Some(Expr::Mul(
            Box::new(Expr::BigInt(i.clone())),
            Box::new(Expr::Variable(var.to_string())),
        )),
        Expr::Rational(r) => Some(Expr::Mul(
            Box::new(Expr::Rational(r.clone())),
            Box::new(Expr::Variable(var.to_string())),
        )),
        Expr::Variable(name) if name == var => Some(Expr::Div(
            Box::new(Expr::Power(
                Box::new(Expr::Variable(var.to_string())),
                Box::new(Expr::BigInt(BigInt::from(2))),
            )),
            Box::new(Expr::BigInt(BigInt::from(2))),
        )),

        // Exponential
        Expr::Exp(arg) => {
            if let Expr::Variable(name) = &**arg {
                if name == var {
                    return Some(Expr::Exp(Box::new(Expr::Variable(var.to_string()))));
                }
            }
            // Handle e^(ax)
            if let Expr::Mul(a, x) = &**arg {
                if let (Expr::Constant(coeff), Expr::Variable(v)) = (&**a, &**x) {
                    if v == var {
                        return Some(Expr::Div(
                            Box::new(expr.clone()),
                            Box::new(Expr::Constant(*coeff)),
                        ));
                    }
                }
                if let (Expr::Variable(v), Expr::Constant(coeff)) = (&**x, &**a) {
                    if v == var {
                        return Some(Expr::Div(
                            Box::new(expr.clone()),
                            Box::new(Expr::Constant(*coeff)),
                        ));
                    }
                }
            }
            None
        }

        // Logarithms
        Expr::Log(arg) => {
            if let Expr::Variable(name) = &**arg {
                if name == var {
                    let x = Expr::Variable(var.to_string());
                    // integral of ln(x) is x*ln(x) - x
                    return Some(Expr::Sub(
                        Box::new(Expr::Mul(Box::new(x.clone()), Box::new(expr.clone()))),
                        Box::new(x),
                    ));
                }
            }
            None
        }
        Expr::LogBase(base, arg) => {
            if let Expr::Variable(name) = &**arg {
                if name == var && !contains_var(base, var) {
                    // Convert log_b(x) to ln(x)/ln(b)
                    let ln_x = Expr::Log(arg.clone());
                    let ln_b = Expr::Log(base.clone());
                    let new_expr = Expr::Div(Box::new(ln_x), Box::new(ln_b));
                    // The integral is (1/ln(b)) * (x*ln(x) - x)
                    return integrate(&new_expr, var, None, None).into();
                }
            }
            None
        }

        // Division rules (e.g., 1/x, 1/(a^2+x^2))
        Expr::Div(num, den) => {
            // Rule: 1/x -> ln|x|
            if let (Expr::BigInt(one), Expr::Variable(name)) = (&**num, &**den) {
                if one.is_one() && name == var {
                    return Some(Expr::Log(Box::new(Expr::Abs(Box::new(Expr::Variable(
                        var.to_string(),
                    ))))));
                }
            }
            // Rule: 1/(a^2 + x^2) -> (1/a)atan(x/a)
            if let Expr::BigInt(one) = &**num {
                if one.is_one() {
                    if let Expr::Add(part1, part2) = &**den {
                        // Handle a^2 + x^2 and x^2 + a^2 by identifying which part is the power
                        let (a_sq_box, x_sq_box) = if let Expr::Power(_, _) = &**part1 {
                            (part2, part1)
                        } else {
                            (part1, part2)
                        };

                        if let (Expr::Constant(a_val), Expr::Power(x, two)) =
                            (&**a_sq_box, &**x_sq_box)
                        {
                            if let (Expr::Variable(v), Expr::Constant(val)) = (&**x, &**two) {
                                if v == var && *val == 2.0 {
                                    let a = Expr::Constant(a_val.sqrt());
                                    return Some(Expr::Mul(
                                        Box::new(Expr::Div(
                                            Box::new(Expr::BigInt(BigInt::one())),
                                            Box::new(a.clone()),
                                        )),
                                        Box::new(Expr::ArcTan(Box::new(Expr::Div(
                                            Box::new(Expr::Variable(var.to_string())),
                                            Box::new(a),
                                        )))),
                                    ));
                                }
                            }
                        }
                    }
                }
            }
            // Rule: 1/sqrt(a^2 - x^2) -> asin(x/a)
            if let (Expr::BigInt(one), Expr::Sqrt(sqrt_arg)) = (&**num, &**den) {
                if one.is_one() {
                    if let Expr::Sub(a_sq, x_sq) = &**sqrt_arg {
                        if let (Expr::Constant(a_val), Expr::Power(x, two)) = (&**a_sq, &**x_sq) {
                            if let (Expr::Variable(v), Expr::Constant(val)) = (&**x, &**two) {
                                if v == var && *val == 2.0 {
                                    let a = Expr::Constant(a_val.sqrt());
                                    return Some(Expr::ArcSin(Box::new(Expr::Div(
                                        Box::new(Expr::Variable(var.to_string())),
                                        Box::new(a),
                                    ))));
                                }
                            }
                        }
                    }
                }
            }
            None
        }

        // Trigonometric functions
        Expr::Sin(arg) => {
            if let Expr::Variable(name) = &**arg {
                if name == var {
                    return Some(Expr::Neg(Box::new(Expr::Cos(Box::new(Expr::Variable(
                        var.to_string(),
                    ))))));
                }
            }
            None
        }
        Expr::Cos(arg) => {
            if let Expr::Variable(name) = &**arg {
                if name == var {
                    return Some(Expr::Sin(Box::new(Expr::Variable(var.to_string()))));
                }
            }
            None
        }
        Expr::Tan(arg) => {
            if let Expr::Variable(name) = &**arg {
                if name == var {
                    return Some(Expr::Log(Box::new(Expr::Abs(Box::new(Expr::Sec(
                        Box::new(Expr::Variable(var.to_string())),
                    ))))));
                }
            }
            None
        }
        Expr::Sec(arg) => {
            if let Expr::Variable(name) = &**arg {
                if name == var {
                    return Some(Expr::Log(Box::new(Expr::Abs(Box::new(Expr::Add(
                        Box::new(Expr::Sec(Box::new(Expr::Variable(var.to_string())))),
                        Box::new(Expr::Tan(Box::new(Expr::Variable(var.to_string())))),
                    ))))));
                }
            }
            None
        }
        Expr::Csc(arg) => {
            if let Expr::Variable(name) = &**arg {
                if name == var {
                    return Some(Expr::Log(Box::new(Expr::Abs(Box::new(Expr::Sub(
                        Box::new(Expr::Csc(Box::new(Expr::Variable(var.to_string())))),
                        Box::new(Expr::Cot(Box::new(Expr::Variable(var.to_string())))),
                    ))))));
                }
            }
            None
        }
        Expr::Cot(arg) => {
            if let Expr::Variable(name) = &**arg {
                if name == var {
                    return Some(Expr::Log(Box::new(Expr::Abs(Box::new(Expr::Sin(
                        Box::new(Expr::Variable(var.to_string())),
                    ))))));
                }
            }
            None
        }
        // sec^2(x) -> tan(x)
        Expr::Power(base, exp)
            if matches!(&**base, Expr::Sec(_))
                && matches!(&**exp, Expr::BigInt(b) if *b == BigInt::from(2)) =>
        {
            if let Expr::Sec(arg) = &**base {
                if let Expr::Variable(name) = &**arg {
                    if name == var {
                        return Some(Expr::Tan(Box::new(Expr::Variable(var.to_string()))));
                    }
                }
            }
            None
        }
        // csc^2(x) -> -cot(x)
        Expr::Power(base, exp)
            if matches!(&**base, Expr::Csc(_))
                && matches!(&**exp, Expr::BigInt(b) if *b == BigInt::from(2)) =>
        {
            if let Expr::Csc(arg) = &**base {
                if let Expr::Variable(name) = &**arg {
                    if name == var {
                        return Some(Expr::Neg(Box::new(Expr::Cot(Box::new(Expr::Variable(
                            var.to_string(),
                        ))))));
                    }
                }
            }
            None
        }
        // sec(x)tan(x) -> sec(x)
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
                            return Some(Expr::Sec(Box::new(Expr::Variable(var.to_string()))));
                        }
                    }
                }
            }
            None
        }
        // csc(x)cot(x) -> -csc(x)
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
                            return Some(Expr::Neg(Box::new(Expr::Csc(Box::new(Expr::Variable(
                                var.to_string(),
                            ))))));
                        }
                    }
                }
            }
            None
        }

        // Inverse Trigonometric Functions
        Expr::ArcTan(arg) => {
            if let Expr::Variable(name) = &**arg {
                if name == var {
                    let x = Expr::Variable(var.to_string());
                    // integral of atan(x) is x*atan(x) - 1/2*ln(1+x^2)
                    let term1 = Expr::Mul(Box::new(x.clone()), Box::new(expr.clone()));
                    let term2 = Expr::Mul(
                        Box::new(Expr::Constant(0.5)),
                        Box::new(Expr::Log(Box::new(Expr::Add(
                            Box::new(Expr::BigInt(BigInt::one())),
                            Box::new(Expr::Power(
                                Box::new(x),
                                Box::new(Expr::BigInt(BigInt::from(2))),
                            )),
                        )))),
                    );
                    return Some(Expr::Sub(Box::new(term1), Box::new(term2)));
                }
            }
            None
        }
        Expr::ArcSin(arg) => {
            if let Expr::Variable(name) = &**arg {
                if name == var {
                    let x = Expr::Variable(var.to_string());
                    // integral of asin(x) is x*asin(x) + sqrt(1-x^2)
                    let term1 = Expr::Mul(Box::new(x.clone()), Box::new(expr.clone()));
                    let term2 = Expr::Sqrt(Box::new(Expr::Sub(
                        Box::new(Expr::BigInt(BigInt::one())),
                        Box::new(Expr::Power(
                            Box::new(x),
                            Box::new(Expr::BigInt(BigInt::from(2))),
                        )),
                    )));
                    return Some(Expr::Add(Box::new(term1), Box::new(term2)));
                }
            }
            None
        }

        // Hyperbolic functions
        Expr::Sinh(arg) => {
            if let Expr::Variable(name) = &**arg {
                if name == var {
                    return Some(Expr::Cosh(Box::new(Expr::Variable(var.to_string()))));
                }
            }
            None
        }
        Expr::Cosh(arg) => {
            if let Expr::Variable(name) = &**arg {
                if name == var {
                    return Some(Expr::Sinh(Box::new(Expr::Variable(var.to_string()))));
                }
            }
            None
        }
        Expr::Tanh(arg) => {
            if let Expr::Variable(name) = &**arg {
                if name == var {
                    return Some(Expr::Log(Box::new(Expr::Cosh(Box::new(Expr::Variable(
                        var.to_string(),
                    ))))));
                }
            }
            None
        }
        // sech^2(x) -> tanh(x)
        Expr::Power(base, exp)
            if matches!(&**base, Expr::Sech(_))
                && matches!(&**exp, Expr::BigInt(b) if *b == BigInt::from(2)) =>
        {
            if let Expr::Sech(arg) = &**base {
                if let Expr::Variable(name) = &**arg {
                    if name == var {
                        return Some(Expr::Tanh(Box::new(Expr::Variable(var.to_string()))));
                    }
                }
            }
            None
        }
        // Power rule
        Expr::Power(base, exp) => {
            if let (Expr::Variable(name), Expr::Constant(n)) = (&**base, &**exp) {
                if name == var {
                    if (*n + 1.0).abs() < 1e-9 {
                        // n == -1
                        return Some(Expr::Log(Box::new(Expr::Abs(Box::new(Expr::Variable(
                            var.to_string(),
                        ))))));
                    }
                    return Some(Expr::Div(
                        Box::new(Expr::Power(
                            Box::new(Expr::Variable(var.to_string())),
                            Box::new(Expr::Constant(n + 1.0)),
                        )),
                        Box::new(Expr::Constant(n + 1.0)),
                    ));
                }
            }
            None
        }

        // Add more rules here...
        _ => None,
    }
}

// Tabular integration by parts, for integrands like polynomial * other_function.
pub(crate) fn integrate_by_parts_tabular(expr: &Expr, var: &str) -> Option<Expr> {
    if let Expr::Mul(part1, part2) = expr {
        // Determine which part is the polynomial.
        let (poly_part, other_part) = if is_polynomial(part1, var) {
            (part1, part2)
        } else if is_polynomial(part2, var) {
            (part2, part1)
        } else {
            return None;
        };

        // Create the derivatives column (D) by differentiating the polynomial until it's zero.
        let mut derivatives = vec![poly_part.clone()];
        while !is_zero(&simplify(*derivatives.last().unwrap().clone())) {
            derivatives.push(Box::new(differentiate(derivatives.last().unwrap(), var)));
        }
        derivatives.pop(); // Remove the final zero.

        // Create the integrals column (I) by integrating the other part.
        let mut integrals = vec![other_part.clone()];
        for _ in 0..derivatives.len() {
            let next_integral = integrate(integrals.last().unwrap(), var, None, None);
            // If we can't integrate the other_part at any step, tabular method fails.
            if let Expr::Integral { .. } = next_integral {
                return None;
            }
            integrals.push(Box::new(simplify(next_integral)));
        }

        // If we don't have enough integrals, something went wrong.
        if derivatives.len() >= integrals.len() {
            return None;
        }

        // Sum the products of the diagonal terms with alternating signs.
        let mut total = Expr::BigInt(BigInt::zero());
        let mut sign = 1;
        for i in 0..derivatives.len() {
            let term = Expr::Mul(
                Box::new(*derivatives[i].clone()),
                Box::new(*integrals[i + 1].clone()),
            );
            if sign == 1 {
                total = Expr::Add(Box::new(total), Box::new(term));
            } else {
                total = Expr::Sub(Box::new(total), Box::new(term));
            }
            sign *= -1;
        }

        return Some(simplify(total));
    }
    None
}

// Master function for integration by parts.
pub(crate) fn integrate_by_parts_master(expr: &Expr, var: &str, depth: u32) -> Option<Expr> {
    // First, try the powerful tabular method for applicable cases.
    if depth == 0 {
        // Only try tabular method on the first call to avoid infinite loops
        if let Some(result) = integrate_by_parts_tabular(expr, var) {
            return Some(result);
        }
    }

    // If tabular method is not applicable, fall back to the standard single-step IBP.
    integrate_by_parts(expr, var, depth)
}

// Finds all unique roots of an expression and their multiplicities.
pub(crate) fn find_roots_with_multiplicity(expr: &Expr, var: &str) -> Vec<(Expr, usize)> {
    let unique_poles = solve(expr, var);
    let mut roots_with_multiplicity = Vec::new();
    let mut processed_poles = std::collections::HashSet::new();

    for pole in unique_poles {
        // Using a HashSet to ensure we process each unique pole only once,
        // as the solver might return duplicates.
        if !processed_poles.insert(pole.clone()) {
            continue;
        }

        let mut m = 1;
        let mut current_deriv = expr.clone();
        while m < 10 {
            // Safety break to avoid infinite loops in complex cases
            let next_deriv = differentiate(&current_deriv, var);
            let val_at_pole = simplify(evaluate_at_point(&next_deriv, var, &pole));
            if !is_zero(&val_at_pole) {
                break; // Found the first non-zero derivative, so multiplicity is m.
            }
            m += 1;
            current_deriv = next_deriv;
        }
        roots_with_multiplicity.push((pole, m));
    }

    roots_with_multiplicity
}

// Integration by partial fractions, now with support for repeated roots and long division.
pub(crate) fn integrate_by_partial_fractions(expr: &Expr, var: &str) -> Option<Expr> {
    if let Expr::Div(num, den) = expr {
        // Use functions from the polynomial module.
        // We prefer the _coeffs version for long division as it's generally more robust.
        use crate::symbolic::polynomial::{polynomial_degree, polynomial_long_division_coeffs};

        let num_deg = polynomial_degree(num, var);
        let den_deg = polynomial_degree(den, var);

        // Step 1: Perform polynomial long division if the fraction is improper.
        if num_deg >= 0 && den_deg >= 0 && num_deg >= den_deg {
            let (quotient, remainder) = polynomial_long_division_coeffs(num, den, var);

            // The integral of the quotient (which is a simple polynomial).
            let integral_of_quotient = integrate(&quotient, var, None, None);

            // The integral of the remainder fraction, which is now a proper fraction.
            let integral_of_remainder = if is_zero(&remainder) {
                Expr::BigInt(BigInt::zero())
            } else {
                let remainder_fraction = Expr::Div(Box::new(remainder), Box::new(*den.clone()));
                integrate(&remainder_fraction, var, None, None)
            };

            // The total integral is the sum of the two parts.
            return Some(simplify(Expr::Add(
                Box::new(integral_of_quotient),
                Box::new(integral_of_remainder),
            )));
        }

        // Step 2: If it's a proper fraction, proceed with partial fraction decomposition.
        let roots = find_roots_with_multiplicity(den, var);

        if roots.is_empty() {
            return None;
        }

        let mut total_integral = Expr::BigInt(BigInt::zero());

        for (root, m) in roots {
            let term_to_multiply = Expr::Power(
                Box::new(Expr::Sub(
                    Box::new(Expr::Variable(var.to_string())),
                    Box::new(root.clone()),
                )),
                Box::new(Expr::BigInt(BigInt::from(m))),
            );

            let g_z = simplify(Expr::Mul(
                Box::new(expr.clone()),
                Box::new(term_to_multiply),
            ));

            for k in 0..m {
                let mut deriv_g = g_z.clone();
                for _ in 0..k {
                    deriv_g = differentiate(&deriv_g, var);
                }

                let val_at_root = evaluate_at_point(&deriv_g, var, &root);
                let k_factorial = Expr::Constant(factorial(k));
                let coefficient = simplify(Expr::Div(Box::new(val_at_root), Box::new(k_factorial)));

                let j = m - k;
                let integral_term = if j == 1 {
                    let log_arg = Expr::Abs(Box::new(simplify(Expr::Sub(
                        Box::new(Expr::Variable(var.to_string())),
                        Box::new(root.clone()),
                    ))));
                    Expr::Mul(
                        Box::new(coefficient),
                        Box::new(Expr::Log(Box::new(log_arg))),
                    )
                } else {
                    let new_power = 1 - (j as i32);
                    let new_denom = Expr::Constant(new_power as f64);
                    let integrated_power_term = Expr::Power(
                        Box::new(Expr::Sub(
                            Box::new(Expr::Variable(var.to_string())),
                            Box::new(root.clone()),
                        )),
                        Box::new(Expr::Constant(new_power as f64)),
                    );
                    Expr::Mul(
                        Box::new(coefficient),
                        Box::new(Expr::Div(
                            Box::new(integrated_power_term),
                            Box::new(new_denom),
                        )),
                    )
                };
                total_integral =
                    simplify(Expr::Add(Box::new(total_integral), Box::new(integral_term)));
            }
        }
        return Some(total_integral);
    }

    None
}

// Helper to check for trig functions.
pub(crate) fn contains_trig_function(expr: &Expr) -> bool {
    match expr {
        Expr::Sin(_) | Expr::Cos(_) | Expr::Tan(_) | Expr::Sec(_) | Expr::Csc(_) | Expr::Cot(_) => {
            true
        }
        Expr::Add(a, b) | Expr::Sub(a, b) | Expr::Mul(a, b) | Expr::Div(a, b) => {
            contains_trig_function(a) || contains_trig_function(b)
        }
        Expr::Power(base, exp) => contains_trig_function(base) || contains_trig_function(exp),
        Expr::Log(arg) | Expr::Abs(arg) | Expr::Neg(arg) | Expr::Exp(arg) => {
            contains_trig_function(arg)
        }
        Expr::Complex(re, im) => contains_trig_function(re) || contains_trig_function(im),
        _ => false,
    }
}

// Tangent half-angle substitution (Weierstrass substitution)
pub(crate) fn tangent_half_angle_substitution(expr: &Expr, var: &str) -> Option<Expr> {
    // This substitution is used for rational functions of trigonometric functions.
    // We first check if the expression contains any trig function, and if so, try to apply it.
    if !contains_trig_function(expr) {
        return None;
    }

    // Let t = tan(x/2)
    let t = Expr::Variable("t".to_string());
    let t_squared = Expr::Power(Box::new(t.clone()), Box::new(Expr::BigInt(BigInt::from(2))));
    let one_plus_t_squared = Expr::Add(
        Box::new(Expr::BigInt(BigInt::one())),
        Box::new(t_squared.clone()),
    );

    // Define substitution rules
    let sin_x_sub = Expr::Div(
        Box::new(Expr::Mul(
            Box::new(Expr::BigInt(BigInt::from(2))),
            Box::new(t.clone()),
        )),
        Box::new(one_plus_t_squared.clone()),
    );
    let cos_x_sub = Expr::Div(
        Box::new(Expr::Sub(
            Box::new(Expr::BigInt(BigInt::one())),
            Box::new(t_squared.clone()),
        )),
        Box::new(one_plus_t_squared.clone()),
    );
    let tan_x_sub = simplify(Expr::Div(
        Box::new(sin_x_sub.clone()),
        Box::new(cos_x_sub.clone()),
    ));
    let dx_sub = Expr::Div(
        Box::new(Expr::BigInt(BigInt::from(2))),
        Box::new(one_plus_t_squared.clone()),
    );

    // Substitute all trig functions in the original expression
    let mut sub_expr = expr.clone();
    let x = Expr::Variable(var.to_string());
    sub_expr = substitute_expr(&sub_expr, &Expr::Sin(Box::new(x.clone())), &sin_x_sub);
    sub_expr = substitute_expr(&sub_expr, &Expr::Cos(Box::new(x.clone())), &cos_x_sub);
    sub_expr = substitute_expr(&sub_expr, &Expr::Tan(Box::new(x.clone())), &tan_x_sub);
    // A full implementation would also substitute sec, csc, cot based on sin and cos.

    // The new integrand is the substituted expression multiplied by dx/dt.
    let new_integrand = simplify(Expr::Mul(Box::new(sub_expr), Box::new(dx_sub)));

    // The result of this substitution should be a rational function of t.
    // We can now call the main integrate function to solve it.
    let integral_in_t = integrate(&new_integrand, "t", None, None);

    // If integration fails (returns an unevaluated Integral), this method was not applicable.
    if let Expr::Integral { .. } = integral_in_t {
        return None;
    }

    // Substitute back t = tan(var/2)
    let t_sub_back = Expr::Tan(Box::new(Expr::Div(
        Box::new(Expr::Variable(var.to_string())),
        Box::new(Expr::BigInt(BigInt::from(2))),
    )));
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
pub(crate) fn limit_internal(expr: &Expr, var: &str, to: &Expr, depth: u32) -> Expr {
    // Safety break for deep recursion, which can happen with L'Hopital's rule.
    if depth > 7 {
        return Expr::Limit(
            Box::new(expr.clone()),
            var.to_string(),
            Box::new(to.clone()),
        );
    }

    // First, try simplifying the expression.
    let expr = &simplify(expr.clone());

    // Strategy 1: Handle limits at +/- infinity for base cases.
    match to {
        Expr::Infinity => match expr {
            Expr::Exp(arg) if **arg == Expr::Variable(var.to_string()) => return Expr::Infinity,
            Expr::Log(arg) if **arg == Expr::Variable(var.to_string()) => return Expr::Infinity,
            Expr::ArcTan(arg) if **arg == Expr::Variable(var.to_string()) => {
                return Expr::Constant(std::f64::consts::PI / 2.0)
            }
            Expr::Variable(v) if v == var => return Expr::Infinity,
            _ => {}
        },
        Expr::NegativeInfinity => match expr {
            Expr::Exp(arg) if **arg == Expr::Variable(var.to_string()) => {
                return Expr::BigInt(BigInt::zero())
            }
            Expr::ArcTan(arg) if **arg == Expr::Variable(var.to_string()) => {
                return Expr::Constant(-std::f64::consts::PI / 2.0)
            }
            Expr::Variable(v) if v == var => return Expr::NegativeInfinity,
            _ => {}
        },
        _ => {}
    }

    // Strategy 2: Direct substitution.
    // If the expression does not contain the variable, the limit is the expression itself.
    if !contains_var(expr, var) {
        return expr.clone();
    }
    let val_at_point = simplify(evaluate_at_point(expr, var, to));
    // If substitution results in a concrete value (not infinity), return it.
    if !matches!(val_at_point, Expr::Infinity | Expr::NegativeInfinity)
        && !contains_var(&val_at_point, var)
    {
        return val_at_point;
    }

    // Strategy 3: Check for indeterminate forms and apply transformations or L'Hopital's Rule.
    match expr {
        Expr::Div(num, den) => {
            let num_limit = limit_internal(num, var, to, depth + 1);
            let den_limit = limit_internal(den, var, to, depth + 1);

            let is_num_zero = is_zero(&num_limit);
            let is_den_zero = is_zero(&den_limit);
            let is_num_inf = matches!(num_limit, Expr::Infinity | Expr::NegativeInfinity);
            let is_den_inf = matches!(den_limit, Expr::Infinity | Expr::NegativeInfinity);

            // L'Hopital's Rule for 0/0 or inf/inf
            if (is_num_zero && is_den_zero) || (is_num_inf && is_den_inf) {
                let d_num = differentiate(num, var);
                let d_den = differentiate(den, var);
                // If the denominator derivative is zero, L'Hopital's rule is inconclusive or leads to infinity.
                if is_zero(&d_den) {
                    return Expr::Infinity; // Or undefined, but Infinity is a common case.
                }
                return limit_internal(
                    &Expr::Div(Box::new(d_num), Box::new(d_den)),
                    var,
                    to,
                    depth + 1,
                );
            }
        }

        Expr::Mul(a, b) => {
            // Handles 0 * inf
            let a_limit = limit_internal(a, var, to, depth + 1);
            let b_limit = limit_internal(b, var, to, depth + 1);
            if is_zero(&a_limit) && matches!(b_limit, Expr::Infinity | Expr::NegativeInfinity) {
                // Rewrite a*b as a / (1/b)
                let new_expr = Expr::Div(
                    a.clone(),
                    Box::new(Expr::Div(Box::new(Expr::BigInt(BigInt::one())), b.clone())),
                );
                return limit_internal(&new_expr, var, to, depth + 1);
            } else if is_zero(&b_limit)
                && matches!(a_limit, Expr::Infinity | Expr::NegativeInfinity)
            {
                // Rewrite a*b as b / (1/a)
                let new_expr = Expr::Div(
                    b.clone(),
                    Box::new(Expr::Div(Box::new(Expr::BigInt(BigInt::one())), a.clone())),
                );
                return limit_internal(&new_expr, var, to, depth + 1);
            }
        }

        Expr::Power(base, exp) => {
            let base_limit = limit_internal(base, var, to, depth + 1);
            let exp_limit = limit_internal(exp, var, to, depth + 1);

            // Check for 1^inf, 0^0, inf^0
            let is_base_one = is_zero(&simplify(Expr::Sub(
                Box::new(base_limit.clone()),
                Box::new(Expr::BigInt(BigInt::one())),
            )));
            let is_base_zero = is_zero(&base_limit);
            let is_base_inf = matches!(base_limit, Expr::Infinity | Expr::NegativeInfinity);
            let is_exp_inf = matches!(exp_limit, Expr::Infinity | Expr::NegativeInfinity);
            let is_exp_zero = is_zero(&exp_limit);

            if (is_base_one && is_exp_inf)
                || (is_base_zero && is_exp_zero)
                || (is_base_inf && is_exp_zero)
            {
                // Transform y = f(x)^g(x) into exp(g(x) * ln(f(x)))
                let log_expr = Expr::Mul(exp.clone(), Box::new(Expr::Log(base.clone())));
                let log_limit = limit_internal(&log_expr, var, to, depth + 1);

                if !contains_var(&log_limit, var) {
                    return Expr::Exp(Box::new(log_limit));
                }
            }
        }
        _ => {}
    }

    // Strategy 4: Fallback for rational functions at infinity (from original implementation)
    if let Expr::Infinity | Expr::NegativeInfinity = to {
        if let Expr::Div(num, den) = expr {
            if is_polynomial(num, var) && is_polynomial(den, var) {
                let deg_num = polynomial_degree(num, var);
                let deg_den = polynomial_degree(den, var);

                if deg_num < deg_den {
                    return Expr::BigInt(BigInt::zero());
                } else if deg_num > deg_den {
                    // A more advanced version could check signs of leading coefficients
                    return if matches!(to, Expr::NegativeInfinity) {
                        Expr::NegativeInfinity
                    } else {
                        Expr::Infinity
                    };
                } else {
                    let lead_num = leading_coefficient(num, var);
                    let lead_den = leading_coefficient(den, var);
                    return simplify(Expr::Div(Box::new(lead_num), Box::new(lead_den)));
                }
            }
        }
    }

    // If no rule applies, return the direct substitution result or the unevaluated limit.
    if !contains_var(&val_at_point, var) {
        val_at_point
    } else {
        Expr::Limit(
            Box::new(expr.clone()),
            var.to_string(),
            Box::new(to.clone()),
        )
    }
}

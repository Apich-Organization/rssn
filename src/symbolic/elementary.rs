//! # Elementary Functions and Expression Manipulation
//!
//! This module provides constructor functions for creating elementary mathematical expressions
//! (like trigonometric, exponential, and power functions) and tools for manipulating
//! these expressions, such as `expand` using algebraic and trigonometric identities.
use crate::symbolic::core::Expr;
use crate::symbolic::simplify::heuristic_simplify;
use num_bigint::BigInt;
use num_traits::{One, ToPrimitive, Zero};
use std::sync::Arc;
/// Creates a sine expression: `sin(expr)`.
pub fn sin(expr: Expr) -> Expr {
    Expr::new_sin(expr)
}
/// Creates a cosine expression: `cos(expr)`.
pub fn cos(expr: Expr) -> Expr {
    Expr::new_cos(expr)
}
/// Creates a tangent expression: `tan(expr)`.
pub fn tan(expr: Expr) -> Expr {
    Expr::new_tan(expr)
}
/// Creates a hyperbolic sine expression: `sinh(expr)`.
pub fn sinh(expr: Expr) -> Expr {
    Expr::new_sinh(expr)
}
/// Creates a hyperbolic cosine expression: `cosh(expr)`.
pub fn cosh(expr: Expr) -> Expr {
    Expr::new_cosh(expr)
}
/// Creates a hyperbolic tangent expression: `tanh(expr)`.
pub fn tanh(expr: Expr) -> Expr {
    Expr::new_tanh(expr)
}
/// Creates a natural logarithm expression: `ln(expr)`.
pub fn ln(expr: Expr) -> Expr {
    Expr::new_log(expr)
}
/// Creates an exponential expression: `e^(expr)`.
pub fn exp(expr: Expr) -> Expr {
    Expr::new_exp(expr)
}
/// Creates a square root expression: `sqrt(expr)`.
pub fn sqrt(expr: Expr) -> Expr {
    Expr::new_sqrt(expr)
}
/// Creates a power expression: `base^exp`.
pub fn pow(base: Expr, exp: Expr) -> Expr {
    Expr::new_pow(base, exp)
}
/// Returns the symbolic representation of positive infinity.
pub const fn infinity() -> Expr {
    Expr::Infinity
}
/// Returns the symbolic representation of negative infinity.
pub const fn negative_infinity() -> Expr {
    Expr::NegativeInfinity
}
/// Creates a logarithm expression with a specified base: `log_base(expr)`.
pub fn log_base(base: Expr, expr: Expr) -> Expr {
    Expr::new_log_base(base, expr)
}
/// Creates a cotangent expression: `cot(expr)`.
pub fn cot(expr: Expr) -> Expr {
    Expr::new_cot(expr)
}
/// Creates a secant expression: `sec(expr)`.
pub fn sec(expr: Expr) -> Expr {
    Expr::new_sec(expr)
}
/// Creates a cosecant expression: `csc(expr)`.
pub fn csc(expr: Expr) -> Expr {
    Expr::new_csc(expr)
}
/// Creates an inverse cotangent expression: `acot(expr)`.
pub fn acot(expr: Expr) -> Expr {
    Expr::new_arccot(expr)
}
/// Creates an inverse secant expression: `asec(expr)`.
pub fn asec(expr: Expr) -> Expr {
    Expr::new_arcsec(expr)
}
/// Creates an inverse cosecant expression: `acsc(expr)`.
pub fn acsc(expr: Expr) -> Expr {
    Expr::new_arccsc(expr)
}
/// Creates a hyperbolic cotangent expression: `coth(expr)`.
pub fn coth(expr: Expr) -> Expr {
    Expr::new_coth(expr)
}
/// Creates a hyperbolic secant expression: `sech(expr)`.
pub fn sech(expr: Expr) -> Expr {
    Expr::new_sech(expr)
}
/// Creates a hyperbolic cosecant expression: `csch(expr)`.
pub fn csch(expr: Expr) -> Expr {
    Expr::new_csch(expr)
}
/// Creates an inverse hyperbolic sine expression: `asinh(expr)`.
pub fn asinh(expr: Expr) -> Expr {
    Expr::new_arcsinh(expr)
}
/// Creates an inverse hyperbolic cosine expression: `acosh(expr)`.
pub fn acosh(expr: Expr) -> Expr {
    Expr::new_arccosh(expr)
}
/// Creates an inverse hyperbolic tangent expression: `atanh(expr)`.
pub fn atanh(expr: Expr) -> Expr {
    Expr::new_arctanh(expr)
}
/// Creates an inverse hyperbolic cotangent expression: `acoth(expr)`.
pub fn acoth(expr: Expr) -> Expr {
    Expr::new_arccoth(expr)
}
/// Creates an inverse hyperbolic secant expression: `asech(expr)`.
pub fn asech(expr: Expr) -> Expr {
    Expr::new_arcsech(expr)
}
/// Creates an inverse hyperbolic cosecant expression: `acsch(expr)`.
pub fn acsch(expr: Expr) -> Expr {
    Expr::new_arccsch(expr)
}
/// Creates a 2-argument inverse tangent expression: `atan2(y, x)`.
pub fn atan2(y: Expr, x: Expr) -> Expr {
    Expr::new_atan2(y, x)
}
/// Returns the symbolic representation of Pi.
pub const fn pi() -> Expr {
    Expr::Pi
}
/// Returns the symbolic representation of Euler's number (e).
pub const fn e() -> Expr {
    Expr::E
}
/// Expands a symbolic expression by applying distributive, power, and trigonometric identities.
///
/// This is often the reverse of simplification and can reveal the underlying structure of an expression.
/// The expansion is applied recursively to all parts of the expression.
///
/// # Arguments
/// * `expr` - The expression to expand.
///
/// # Returns
/// A new, expanded `Expr`.
pub fn expand(expr: Expr) -> Expr {
    let expanded_expr = match expr {
        Expr::Dag(node) => {
            return expand(node.to_expr().expect("Expand"));
        }
        Expr::Add(a, b) => Expr::new_add(expand((*a).clone()), expand((*b).clone())),
        Expr::Sub(a, b) => Expr::new_sub(expand((*a).clone()), expand((*b).clone())),
        Expr::Mul(a, b) => expand_mul((*a).clone(), (*b).clone()),
        Expr::Div(a, b) => Expr::new_div(expand((*a).clone()), expand((*b).clone())),
        Expr::Power(b, e) => expand_power(&b, &e),
        Expr::Log(arg) => expand_log(&arg),
        Expr::Sin(arg) => expand_sin(&arg),
        Expr::Cos(arg) => expand_cos(&arg),
        _ => expr,
    };
    heuristic_simplify(expanded_expr)
}
/// Expands multiplication over addition: `a*(b+c) -> a*b + a*c`.
pub(crate) fn expand_mul(a: Expr, b: Expr) -> Expr {
    let a_exp = expand(a);
    let b_exp = expand(b);
    match (a_exp, b_exp) {
        (l, Expr::Add(m, n)) => Expr::new_add(
            expand(Expr::new_mul(l.clone(), m)),
            expand(Expr::new_mul(l, n)),
        ),
        (Expr::Add(m, n), r) => Expr::new_add(
            expand(Expr::new_mul(m, r.clone())),
            expand(Expr::new_mul(n, r)),
        ),
        (l, r) => Expr::new_mul(l, r),
    }
}
/// Expands powers, e.g., `(a*b)^c -> a^c * b^c` and `(a+b)^n -> a^n + ...` (binomial expansion).
pub(crate) fn expand_power(base: &Arc<Expr>, exp: &Arc<Expr>) -> Expr {
    let b_exp = expand(base.as_ref().clone());
    let e_exp = expand(exp.as_ref().clone());
    match (b_exp, e_exp) {
        (Expr::Mul(f, g), e) => Expr::new_mul(
            expand(Expr::new_pow(f, e.clone())),
            expand(Expr::new_pow(g, e)),
        ),
        (Expr::Add(a, b), Expr::BigInt(n)) => {
            if let Some(n_usize) = n.to_usize() {
                let mut sum = Expr::BigInt(BigInt::zero());
                for k in 0..=n_usize {
                    let bin_coeff = Expr::BigInt(binomial_coefficient(n_usize, k));
                    let term1 = Expr::new_pow(a.clone(), Expr::BigInt(BigInt::from(k)));
                    let term2 = Expr::new_pow(b.clone(), Expr::BigInt(BigInt::from(n_usize - k)));
                    let term = Expr::new_mul(bin_coeff, Expr::new_mul(term1, term2));
                    sum = Expr::new_add(sum, expand(term));
                }
                sum
            } else {
                Expr::new_pow(Expr::new_add(a, b), Expr::BigInt(n))
            }
        }
        (b, e) => Expr::new_pow(b, e),
    }
}
/// Expands logarithms using identities like `log(a*b) -> log(a) + log(b)`.
pub(crate) fn expand_log(arg: &Arc<Expr>) -> Expr {
    let arg_exp = expand(arg.as_ref().clone());
    match arg_exp {
        Expr::Mul(a, b) => Expr::new_add(expand(Expr::new_log(a)), expand(Expr::new_log(b))),
        Expr::Div(a, b) => Expr::new_sub(expand(Expr::new_log(a)), expand(Expr::new_log(b))),
        Expr::Power(b, e) => Expr::new_mul(e, expand(Expr::new_log(b))),
        a => Expr::new_log(a),
    }
}
/// Expands `sin` using sum-angle identities, e.g., `sin(a+b)`.
pub(crate) fn expand_sin(arg: &Arc<Expr>) -> Expr {
    let arg_exp = expand(arg.as_ref().clone());
    match arg_exp {
        Expr::Add(a, b) => Expr::new_add(
            Expr::new_mul(sin(a.as_ref().clone()), cos(b.as_ref().clone())),
            Expr::new_mul(cos(a.as_ref().clone()), sin(b.as_ref().clone())),
        ),
        a => Expr::new_sin(a),
    }
}
/// Expands `cos` using sum-angle identities, e.g., `cos(a+b)`.
pub(crate) fn expand_cos(arg: &Arc<Expr>) -> Expr {
    let arg_exp = expand(arg.as_ref().clone());
    match arg_exp {
        Expr::Add(a, b) => Expr::new_sub(
            Expr::new_mul(cos(a.as_ref().clone()), cos(b.as_ref().clone())),
            Expr::new_mul(sin(a.as_ref().clone()), sin(b.as_ref().clone())),
        ),
        a => Expr::new_cos(a),
    }
}
/// Helper to compute binomial coefficients C(n, k) = n! / (k! * (n-k)!).
pub(crate) fn binomial_coefficient(n: usize, k: usize) -> BigInt {
    if k > n {
        return BigInt::zero();
    }
    let mut res = BigInt::one();
    for i in 0..k {
        res = (res * (n - i)) / (i + 1);
    }
    res
}

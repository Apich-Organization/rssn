//! # Elementary Functions and Expression Manipulation
//!
//! This module provides constructor functions for creating elementary mathematical expressions
//! (like trigonometric, exponential, and power functions) and tools for manipulating
//! these expressions, such as `expand` using algebraic and trigonometric identities.

use std::sync::Arc;

use crate::symbolic::core::Expr;
use crate::symbolic::simplify::heuristic_simplify;
use num_bigint::BigInt;
use num_traits::{One, ToPrimitive, Zero};

// =====================================================================================
// region: Expression Constructors
// =====================================================================================

/// Creates a sine expression: `sin(expr)`.
pub fn sin(expr: Expr) -> Expr {
    Expr::Sin(Arc::new(expr))
}
/// Creates a cosine expression: `cos(expr)`.
pub fn cos(expr: Expr) -> Expr {
    Expr::Cos(Arc::new(expr))
}
/// Creates a tangent expression: `tan(expr)`.
pub fn tan(expr: Expr) -> Expr {
    Expr::Tan(Arc::new(expr))
}
/// Creates a hyperbolic sine expression: `sinh(expr)`.
pub fn sinh(expr: Expr) -> Expr {
    Expr::Sinh(Arc::new(expr))
}
/// Creates a hyperbolic cosine expression: `cosh(expr)`.
pub fn cosh(expr: Expr) -> Expr {
    Expr::Cosh(Arc::new(expr))
}
/// Creates a hyperbolic tangent expression: `tanh(expr)`.
pub fn tanh(expr: Expr) -> Expr {
    Expr::Tanh(Arc::new(expr))
}
/// Creates a natural logarithm expression: `ln(expr)`.
pub fn ln(expr: Expr) -> Expr {
    Expr::Log(Arc::new(expr))
}
/// Creates an exponential expression: `e^(expr)`.
pub fn exp(expr: Expr) -> Expr {
    Expr::Exp(Arc::new(expr))
}
/// Creates a square root expression: `sqrt(expr)`.
pub fn sqrt(expr: Expr) -> Expr {
    Expr::Sqrt(Arc::new(expr))
}
/// Creates a power expression: `base^exp`.
pub fn pow(base: Expr, exp: Expr) -> Expr {
    Expr::Power(Arc::new(base), Arc::new(exp))
}
/// Returns the symbolic representation of positive infinity.
pub fn infinity() -> Expr {
    Expr::Infinity
}
/// Returns the symbolic representation of negative infinity.
pub fn negative_infinity() -> Expr {
    Expr::NegativeInfinity
}
/// Creates a logarithm expression with a specified base: `log_base(expr)`.
pub fn log_base(base: Expr, expr: Expr) -> Expr {
    Expr::LogBase(Arc::new(base), Arc::new(expr))
}
/// Creates a cotangent expression: `cot(expr)`.
pub fn cot(expr: Expr) -> Expr {
    Expr::Cot(Arc::new(expr))
}
/// Creates a secant expression: `sec(expr)`.
pub fn sec(expr: Expr) -> Expr {
    Expr::Sec(Arc::new(expr))
}
/// Creates a cosecant expression: `csc(expr)`.
pub fn csc(expr: Expr) -> Expr {
    Expr::Csc(Arc::new(expr))
}
/// Creates an inverse cotangent expression: `acot(expr)`.
pub fn acot(expr: Expr) -> Expr {
    Expr::ArcCot(Arc::new(expr))
}
/// Creates an inverse secant expression: `asec(expr)`.
pub fn asec(expr: Expr) -> Expr {
    Expr::ArcSec(Arc::new(expr))
}
/// Creates an inverse cosecant expression: `acsc(expr)`.
pub fn acsc(expr: Expr) -> Expr {
    Expr::ArcCsc(Arc::new(expr))
}
/// Creates a hyperbolic cotangent expression: `coth(expr)`.
pub fn coth(expr: Expr) -> Expr {
    Expr::Coth(Arc::new(expr))
}
/// Creates a hyperbolic secant expression: `sech(expr)`.
pub fn sech(expr: Expr) -> Expr {
    Expr::Sech(Arc::new(expr))
}
/// Creates a hyperbolic cosecant expression: `csch(expr)`.
pub fn csch(expr: Expr) -> Expr {
    Expr::Csch(Arc::new(expr))
}
/// Creates an inverse hyperbolic sine expression: `asinh(expr)`.
pub fn asinh(expr: Expr) -> Expr {
    Expr::ArcSinh(Arc::new(expr))
}
/// Creates an inverse hyperbolic cosine expression: `acosh(expr)`.
pub fn acosh(expr: Expr) -> Expr {
    Expr::ArcCosh(Arc::new(expr))
}
/// Creates an inverse hyperbolic tangent expression: `atanh(expr)`.
pub fn atanh(expr: Expr) -> Expr {
    Expr::ArcTanh(Arc::new(expr))
}
/// Creates an inverse hyperbolic cotangent expression: `acoth(expr)`.
pub fn acoth(expr: Expr) -> Expr {
    Expr::ArcCoth(Arc::new(expr))
}
/// Creates an inverse hyperbolic secant expression: `asech(expr)`.
pub fn asech(expr: Expr) -> Expr {
    Expr::ArcSech(Arc::new(expr))
}
/// Creates an inverse hyperbolic cosecant expression: `acsch(expr)`.
pub fn acsch(expr: Expr) -> Expr {
    Expr::ArcCsch(Arc::new(expr))
}
/// Creates a 2-argument inverse tangent expression: `atan2(y, x)`.
pub fn atan2(y: Expr, x: Expr) -> Expr {
    Expr::Atan2(Arc::new(y), Arc::new(x))
}

/// Returns the symbolic representation of Pi.
pub fn pi() -> Expr {
    Expr::Pi
}

/// Returns the symbolic representation of Euler's number (e).
pub fn e() -> Expr {
    Expr::E
}

// endregion

// =====================================================================================
// region: Expression Expansion
// =====================================================================================

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
        Expr::Add(a, b) => Expr::Add(
            Arc::new(expand((*a).clone())),
            Arc::new(expand((*b).clone())),
        ),
        Expr::Sub(a, b) => Expr::Sub(
            Arc::new(expand((*a).clone())),
            Arc::new(expand((*b).clone())),
        ),
        Expr::Mul(a, b) => expand_mul((*a).clone(), (*b).clone()),
        Expr::Div(a, b) => Expr::Div(
            Arc::new(expand((*a).clone())),
            Arc::new(expand((*b).clone())),
        ),
        Expr::Power(b, e) => expand_power(&b, &e),
        Expr::Log(arg) => expand_log(&arg),
        Expr::Sin(arg) => expand_sin(&arg),
        Expr::Cos(arg) => expand_cos(&arg),
        _ => expr,
    };
    // After expansion, a light simplification can clean up the result (e.g., 1*x -> x)
    heuristic_simplify(expanded_expr)
}

/// Expands multiplication over addition: `a*(b+c) -> a*b + a*c`.
pub(crate) fn expand_mul(a: Expr, b: Expr) -> Expr {
    let a_exp = expand(a);
    let b_exp = expand(b);
    match (a_exp, b_exp) {
        (l, Expr::Add(m, n)) => Expr::Add(
            Arc::new(expand(Expr::Mul(Arc::new(l.clone()), m))),
            Arc::new(expand(Expr::Mul(Arc::new(l), n))),
        ),
        (Expr::Add(m, n), r) => Expr::Add(
            Arc::new(expand(Expr::Mul(m, Arc::new(r.clone())))),
            Arc::new(expand(Expr::Mul(n, Arc::new(r)))),
        ),
        (l, r) => Expr::Mul(Arc::new(l), Arc::new(r)),
    }
}
/// Expands powers, e.g., `(a*b)^c -> a^c * b^c` and `(a+b)^n -> a^n + ...` (binomial expansion).
pub(crate) fn expand_power(base: &Arc<Expr>, exp: &Arc<Expr>) -> Expr {
    let b_exp = expand(base.as_ref().clone());
    let e_exp = expand(exp.as_ref().clone());
    match (b_exp, e_exp) {
        // (a*b)^c -> a^c * b^c
        (Expr::Mul(f, g), e) => Expr::Mul(
            Arc::new(expand(Expr::Power(f, Arc::new(e.clone())))),
            Arc::new(expand(Expr::Power(g, Arc::new(e)))),
        ),
        // (a+b)^n for integer n
        (Expr::Add(a, b), Expr::BigInt(n)) => {
            if let Some(n_usize) = n.to_usize() {
                let mut sum = Expr::BigInt(BigInt::zero());
                for k in 0..=n_usize {
                    let bin_coeff = Expr::BigInt(binomial_coefficient(n_usize, k));
                    let term1 = Expr::Power(a.clone(), Arc::new(Expr::BigInt(BigInt::from(k))));
                    let term2 =
                        Expr::Power(b.clone(), Arc::new(Expr::BigInt(BigInt::from(n_usize - k))));
                    let term = Expr::Mul(
                        Arc::new(bin_coeff),
                        Arc::new(Expr::Mul(Arc::new(term1), Arc::new(term2))),
                    );
                    sum = Expr::Add(Arc::new(sum), Arc::new(expand(term)));
                }
                sum
            } else {
                Expr::Power(Arc::new(Expr::Add(a, b)), Arc::new(Expr::BigInt(n)))
            }
        }
        (b, e) => Expr::Power(Arc::new(b), Arc::new(e)),
    }
}

/// Expands logarithms using identities like `log(a*b) -> log(a) + log(b)`.
pub(crate) fn expand_log(arg: &Arc<Expr>) -> Expr {
    let arg_exp = expand(arg.as_ref().clone());
    match arg_exp {
        // log(a*b) -> log(a) + log(b)
        Expr::Mul(a, b) => Expr::Add(
            Arc::new(expand(Expr::Log(a))),
            Arc::new(expand(Expr::Log(b))),
        ),
        // log(a/b) -> log(a) - log(b)
        Expr::Div(a, b) => Expr::Sub(
            Arc::new(expand(Expr::Log(a))),
            Arc::new(expand(Expr::Log(b))),
        ),
        // log(a^b) -> b*log(a)
        Expr::Power(b, e) => Expr::Mul(e, Arc::new(expand(Expr::Log(b)))),
        a => Expr::Log(Arc::new(a)),
    }
}

/// Expands `sin` using sum-angle identities, e.g., `sin(a+b)`.
pub(crate) fn expand_sin(arg: &Arc<Expr>) -> Expr {
    let arg_exp = expand(arg.as_ref().clone());
    match arg_exp {
        // sin(a+b) -> sin(a)cos(b) + cos(a)sin(b)
        Expr::Add(a, b) => Expr::Add(
            Arc::new(Expr::Mul(
                Arc::new(sin(a.as_ref().clone())),
                Arc::new(cos(b.as_ref().clone())),
            )),
            Arc::new(Expr::Mul(
                Arc::new(cos(a.as_ref().clone())),
                Arc::new(sin(b.as_ref().clone())),
            )),
        ),
        a => Expr::Sin(Arc::new(a)),
    }
}

/// Expands `cos` using sum-angle identities, e.g., `cos(a+b)`.
pub(crate) fn expand_cos(arg: &Arc<Expr>) -> Expr {
    let arg_exp = expand(arg.as_ref().clone());
    match arg_exp {
        // cos(a+b) -> cos(a)cos(b) - sin(a)sin(b)
        Expr::Add(a, b) => Expr::Sub(
            Arc::new(Expr::Mul(
                Arc::new(cos(a.as_ref().clone())),
                Arc::new(cos(b.as_ref().clone())),
            )),
            Arc::new(Expr::Mul(
                Arc::new(sin(a.as_ref().clone())),
                Arc::new(sin(b.as_ref().clone())),
            )),
        ),
        a => Expr::Cos(Arc::new(a)),
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

// endregion

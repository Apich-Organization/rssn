//! # Elementary Functions and Expression Manipulation
//!
//! This module provides constructor functions for creating elementary mathematical expressions
//! (like trigonometric, exponential, and power functions) and tools for manipulating
//! these expressions, such as `expand` using algebraic and trigonometric identities.

use crate::symbolic::core::Expr;
use crate::symbolic::simplify::heuristic_simplify;
use num_bigint::BigInt;
use num_traits::{One, ToPrimitive, Zero};

// =====================================================================================
// region: Expression Constructors
// =====================================================================================

/// Creates a sine expression: `sin(expr)`.
pub fn sin(expr: Expr) -> Expr {
    Expr::Sin(Box::new(expr))
}
/// Creates a cosine expression: `cos(expr)`.
pub fn cos(expr: Expr) -> Expr {
    Expr::Cos(Box::new(expr))
}
/// Creates a tangent expression: `tan(expr)`.
pub fn tan(expr: Expr) -> Expr {
    Expr::Tan(Box::new(expr))
}
/// Creates a hyperbolic sine expression: `sinh(expr)`.
pub fn sinh(expr: Expr) -> Expr {
    Expr::Sinh(Box::new(expr))
}
/// Creates a hyperbolic cosine expression: `cosh(expr)`.
pub fn cosh(expr: Expr) -> Expr {
    Expr::Cosh(Box::new(expr))
}
/// Creates a hyperbolic tangent expression: `tanh(expr)`.
pub fn tanh(expr: Expr) -> Expr {
    Expr::Tanh(Box::new(expr))
}
/// Creates a natural logarithm expression: `ln(expr)`.
pub fn ln(expr: Expr) -> Expr {
    Expr::Log(Box::new(expr))
}
/// Creates an exponential expression: `e^(expr)`.
pub fn exp(expr: Expr) -> Expr {
    Expr::Exp(Box::new(expr))
}
/// Creates a square root expression: `sqrt(expr)`.
pub fn sqrt(expr: Expr) -> Expr {
    Expr::Sqrt(Box::new(expr))
}
/// Creates a power expression: `base^exp`.
pub fn pow(base: Expr, exp: Expr) -> Expr {
    Expr::Power(Box::new(base), Box::new(exp))
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
    Expr::LogBase(Box::new(base), Box::new(expr))
}
/// Creates a cotangent expression: `cot(expr)`.
pub fn cot(expr: Expr) -> Expr {
    Expr::Cot(Box::new(expr))
}
/// Creates a secant expression: `sec(expr)`.
pub fn sec(expr: Expr) -> Expr {
    Expr::Sec(Box::new(expr))
}
/// Creates a cosecant expression: `csc(expr)`.
pub fn csc(expr: Expr) -> Expr {
    Expr::Csc(Box::new(expr))
}
/// Creates an inverse cotangent expression: `acot(expr)`.
pub fn acot(expr: Expr) -> Expr {
    Expr::ArcCot(Box::new(expr))
}
/// Creates an inverse secant expression: `asec(expr)`.
pub fn asec(expr: Expr) -> Expr {
    Expr::ArcSec(Box::new(expr))
}
/// Creates an inverse cosecant expression: `acsc(expr)`.
pub fn acsc(expr: Expr) -> Expr {
    Expr::ArcCsc(Box::new(expr))
}
/// Creates a hyperbolic cotangent expression: `coth(expr)`.
pub fn coth(expr: Expr) -> Expr {
    Expr::Coth(Box::new(expr))
}
/// Creates a hyperbolic secant expression: `sech(expr)`.
pub fn sech(expr: Expr) -> Expr {
    Expr::Sech(Box::new(expr))
}
/// Creates a hyperbolic cosecant expression: `csch(expr)`.
pub fn csch(expr: Expr) -> Expr {
    Expr::Csch(Box::new(expr))
}
/// Creates an inverse hyperbolic sine expression: `asinh(expr)`.
pub fn asinh(expr: Expr) -> Expr {
    Expr::ArcSinh(Box::new(expr))
}
/// Creates an inverse hyperbolic cosine expression: `acosh(expr)`.
pub fn acosh(expr: Expr) -> Expr {
    Expr::ArcCosh(Box::new(expr))
}
/// Creates an inverse hyperbolic tangent expression: `atanh(expr)`.
pub fn atanh(expr: Expr) -> Expr {
    Expr::ArcTanh(Box::new(expr))
}
/// Creates an inverse hyperbolic cotangent expression: `acoth(expr)`.
pub fn acoth(expr: Expr) -> Expr {
    Expr::ArcCoth(Box::new(expr))
}
/// Creates an inverse hyperbolic secant expression: `asech(expr)`.
pub fn asech(expr: Expr) -> Expr {
    Expr::ArcSech(Box::new(expr))
}
/// Creates an inverse hyperbolic cosecant expression: `acsch(expr)`.
pub fn acsch(expr: Expr) -> Expr {
    Expr::ArcCsch(Box::new(expr))
}
/// Creates a 2-argument inverse tangent expression: `atan2(y, x)`.
pub fn atan2(y: Expr, x: Expr) -> Expr {
    Expr::Atan2(Box::new(y), Box::new(x))
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
        Expr::Add(a, b) => Expr::Add(Box::new(expand(*a)), Box::new(expand(*b))),
        Expr::Sub(a, b) => Expr::Sub(Box::new(expand(*a)), Box::new(expand(*b))),
        Expr::Mul(a, b) => expand_mul(*a, *b),
        Expr::Div(a, b) => Expr::Div(Box::new(expand(*a)), Box::new(expand(*b))),
        Expr::Power(b, e) => expand_power(b, e),
        Expr::Log(arg) => expand_log(arg),
        Expr::Sin(arg) => expand_sin(arg),
        Expr::Cos(arg) => expand_cos(arg),
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
            Box::new(expand(Expr::Mul(Box::new(l.clone()), m))),
            Box::new(expand(Expr::Mul(Box::new(l), n))),
        ),
        (Expr::Add(m, n), r) => Expr::Add(
            Box::new(expand(Expr::Mul(m, Box::new(r.clone())))),
            Box::new(expand(Expr::Mul(n, Box::new(r)))),
        ),
        (l, r) => Expr::Mul(Box::new(l), Box::new(r)),
    }
}
/// Expands powers, e.g., `(a*b)^c -> a^c * b^c` and `(a+b)^n -> a^n + ...` (binomial expansion).
pub(crate) fn expand_power(base: Box<Expr>, exp: Box<Expr>) -> Expr {
    let b_exp = expand(*base);
    let e_exp = expand(*exp);
    match (b_exp, e_exp) {
        // (a*b)^c -> a^c * b^c
        (Expr::Mul(f, g), e) => Expr::Mul(
            Box::new(expand(Expr::Power(f, Box::new(e.clone())))),
            Box::new(expand(Expr::Power(g, Box::new(e)))),
        ),
        // (a+b)^n for integer n
        (Expr::Add(a, b), Expr::BigInt(n)) => {
            if let Some(n_usize) = n.to_usize() {
                let mut sum = Expr::BigInt(BigInt::zero());
                for k in 0..=n_usize {
                    let bin_coeff = Expr::BigInt(binomial_coefficient(n_usize, k));
                    let term1 = Expr::Power(a.clone(), Box::new(Expr::BigInt(BigInt::from(k))));
                    let term2 =
                        Expr::Power(b.clone(), Box::new(Expr::BigInt(BigInt::from(n_usize - k))));
                    let term = Expr::Mul(
                        Box::new(bin_coeff),
                        Box::new(Expr::Mul(Box::new(term1), Box::new(term2))),
                    );
                    sum = Expr::Add(Box::new(sum), Box::new(expand(term)));
                }
                sum
            } else {
                Expr::Power(Box::new(Expr::Add(a, b)), Box::new(Expr::BigInt(n)))
            }
        }
        (b, e) => Expr::Power(Box::new(b), Box::new(e)),
    }
}

/// Expands logarithms using identities like `log(a*b) -> log(a) + log(b)`.
pub(crate) fn expand_log(arg: Box<Expr>) -> Expr {
    let arg_exp = expand(*arg);
    match arg_exp {
        // log(a*b) -> log(a) + log(b)
        Expr::Mul(a, b) => Expr::Add(
            Box::new(expand(Expr::Log(a))),
            Box::new(expand(Expr::Log(b))),
        ),
        // log(a/b) -> log(a) - log(b)
        Expr::Div(a, b) => Expr::Sub(
            Box::new(expand(Expr::Log(a))),
            Box::new(expand(Expr::Log(b))),
        ),
        // log(a^b) -> b*log(a)
        Expr::Power(b, e) => Expr::Mul(e, Box::new(expand(Expr::Log(b)))),
        a => Expr::Log(Box::new(a)),
    }
}

/// Expands `sin` using sum-angle identities, e.g., `sin(a+b)`.
pub(crate) fn expand_sin(arg: Box<Expr>) -> Expr {
    let arg_exp = expand(*arg);
    match arg_exp {
        // sin(a+b) -> sin(a)cos(b) + cos(a)sin(b)
        Expr::Add(a, b) => Expr::Add(
            Box::new(Expr::Mul(
                Box::new(sin(*a.clone())),
                Box::new(cos(*b.clone())),
            )),
            Box::new(Expr::Mul(Box::new(cos(*a)), Box::new(sin(*b)))),
        ),
        a => Expr::Sin(Box::new(a)),
    }
}

/// Expands `cos` using sum-angle identities, e.g., `cos(a+b)`.
pub(crate) fn expand_cos(arg: Box<Expr>) -> Expr {
    let arg_exp = expand(*arg);
    match arg_exp {
        // cos(a+b) -> cos(a)cos(b) - sin(a)sin(b)
        Expr::Add(a, b) => Expr::Sub(
            Box::new(Expr::Mul(
                Box::new(cos(*a.clone())),
                Box::new(cos(*b.clone())),
            )),
            Box::new(Expr::Mul(Box::new(sin(*a)), Box::new(sin(*b)))),
        ),
        a => Expr::Cos(Box::new(a)),
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

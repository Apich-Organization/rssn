//! # Multi-valued Function Handling in Complex Analysis
//!
//! This module provides structures and functions for representing and operating on
//! multi-valued complex functions, such as `log(z)` or `sqrt(z)`. It aims to track
//! the different branches of these functions, rather than just returning the principal value.
//! This is achieved by introducing a symbolic integer `k` to represent the branch number.

use std::sync::Arc;

use crate::prelude::simplify;
use crate::symbolic::core::Expr;

/// Returns the principal argument of a complex expression `z`.
/// `Arg(z)` is the angle in radians in the interval (-pi, pi].
pub(crate) fn arg(z: &Expr) -> Expr {
    // This is a symbolic representation. The simplification engine would handle
    // the actual computation, e.g., atan2(im, re).
    Expr::Apply(
        Arc::new(Expr::Variable("Arg".to_string())),
        Arc::new(z.clone()),
    )
}

/// Returns the absolute value (magnitude) of a complex expression `z`.
pub(crate) fn abs(z: &Expr) -> Expr {
    Expr::Abs(Arc::new(z.clone()))
}

/// Computes the general multi-valued logarithm of a complex expression `z`.
///
/// The formula is `log(z) = ln|z| + i * (Arg(z) + 2*pi*k)`,
/// where `k` is an integer representing the branch number.
///
/// # Arguments
/// * `z` - The complex expression.
/// * `k` - A symbolic expression representing an arbitrary integer (e.g., `Expr::Variable("k")`).
///
/// # Returns
/// An `Expr` representing the multi-valued logarithm.
pub fn general_log(z: &Expr, k: &Expr) -> Expr {
    let pi = Expr::Pi;
    let i = Expr::Complex(Arc::new(Expr::Constant(0.0)), Arc::new(Expr::Constant(1.0)));

    let term_2_pi_k = Expr::Mul(
        Arc::new(Expr::Constant(2.0)),
        Arc::new(Expr::Mul(Arc::new(pi), Arc::new(k.clone()))),
    );

    let full_arg = Expr::Add(Arc::new(arg(z)), Arc::new(term_2_pi_k));

    let result = Expr::Add(
        Arc::new(Expr::Log(Arc::new(abs(z)))),
        Arc::new(Expr::Mul(Arc::new(i), Arc::new(full_arg))),
    );

    simplify(result)
}

/// Computes the general multi-valued power `z^w`.
///
/// The formula is `z^w = exp(w * log(z))`, where `log(z)` is the multi-valued logarithm.
///
/// # Arguments
/// * `z` - The base, a complex expression.
/// * `w` - The exponent, a complex expression.
/// * `k` - A symbolic expression representing an arbitrary integer for the logarithm's branch.
///
/// # Returns
/// An `Expr` representing the multi-valued power.
pub fn general_power(z: &Expr, w: &Expr, k: &Expr) -> Expr {
    let log_z = general_log(z, k);
    simplify(Expr::Exp(Arc::new(Expr::Mul(
        Arc::new(w.clone()),
        Arc::new(log_z),
    ))))
}

/// Computes the general multi-valued arcsin of a complex expression `z`.
///
/// The formula is `k*pi + (-1)^k * asin(z)`,
/// where `asin(z)` is the principal value and `k` is an integer.
///
/// # Arguments
/// * `z` - The complex expression.
/// * `k` - A symbolic expression representing an arbitrary integer.
///
/// # Returns
/// An `Expr` representing the multi-valued arcsin.
pub fn general_arcsin(z: &Expr, k: &Expr) -> Expr {
    let pi = Expr::Pi;
    let principal_arcsin = Expr::ArcSin(Arc::new(z.clone()));

    let term1 = Expr::Mul(Arc::new(k.clone()), Arc::new(pi));
    let term2 = Expr::Mul(
        Arc::new(Expr::Power(
            Arc::new(Expr::Constant(-1.0)),
            Arc::new(k.clone()),
        )),
        Arc::new(principal_arcsin),
    );

    simplify(Expr::Add(Arc::new(term1), Arc::new(term2)))
}

/// Computes the general multi-valued arccos of a complex expression `z`.
///
/// The formula is `2*k*pi +/- acos(z)`,
/// where `acos(z)` is the principal value, `k` is an integer, and `s` determines the sign.
///
/// # Arguments
/// * `z` - The complex expression.
/// * `k` - A symbolic expression representing an arbitrary integer.
/// * `s` - A symbolic expression representing the sign (+1 or -1).
///
/// # Returns
/// An `Expr` representing the multi-valued arccos.
pub fn general_arccos(z: &Expr, k: &Expr, s: &Expr) -> Expr {
    let pi = Expr::Pi;
    let principal_arccos = Expr::ArcCos(Arc::new(z.clone()));

    let term1 = Expr::Mul(
        Arc::new(Expr::Constant(2.0)),
        Arc::new(Expr::Mul(Arc::new(k.clone()), Arc::new(pi))),
    );

    let term2 = Expr::Mul(Arc::new(s.clone()), Arc::new(principal_arccos));

    simplify(Expr::Add(Arc::new(term1), Arc::new(term2)))
}

/// Computes the general multi-valued arctan of a complex expression `z`.
///
/// The formula is `k*pi + atan(z)`,
/// where `atan(z)` is the principal value and `k` is an integer.
///
/// # Arguments
/// * `z` - The complex expression.
/// * `k` - A symbolic expression representing an arbitrary integer.
///
/// # Returns
/// An `Expr` representing the multi-valued arctan.
pub fn general_arctan(z: &Expr, k: &Expr) -> Expr {
    let pi = Expr::Pi;
    let principal_arctan = Expr::ArcTan(Arc::new(z.clone()));

    let term1 = Expr::Mul(Arc::new(k.clone()), Arc::new(pi));

    simplify(Expr::Add(Arc::new(term1), Arc::new(principal_arctan)))
}

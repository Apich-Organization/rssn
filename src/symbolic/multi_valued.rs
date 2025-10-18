//! # Multi-valued Function Handling in Complex Analysis
//!
//! This module provides structures and functions for representing and operating on
//! multi-valued complex functions, such as `log(z)` or `sqrt(z)`. It aims to track
//! the different branches of these functions, rather than just returning the principal value.
//! This is achieved by introducing a symbolic integer `k` to represent the branch number.
use crate::symbolic::core::Expr;
use crate::symbolic::simplify::simplify;
/// Returns the principal argument of a complex expression `z`.
/// `Arg(z)` is the angle in radians in the interval (-pi, pi].
pub(crate) fn arg(z: &Expr) -> Expr {
    Expr::new_apply(Expr::Variable("Arg".to_string()), z.clone())
}
/// Returns the absolute value (magnitude) of a complex expression `z`.
pub(crate) fn abs(z: &Expr) -> Expr {
    Expr::new_abs(z.clone())
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
    let i = Expr::new_complex(Expr::Constant(0.0), Expr::Constant(1.0));
    let term_2_pi_k = Expr::new_mul(Expr::Constant(2.0), Expr::new_mul(pi, k.clone()));
    let full_arg = Expr::new_add(arg(z), term_2_pi_k);
    let result = Expr::new_add(Expr::new_log(abs(z)), Expr::new_mul(i, full_arg));
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
    simplify(Expr::new_exp(Expr::new_mul(w.clone(), log_z)))
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
    let principal_arcsin = Expr::new_arcsin(z.clone());
    let term1 = Expr::new_mul(k.clone(), pi);
    let term2 = Expr::new_mul(
        Expr::new_pow(Expr::Constant(-1.0), k.clone()),
        principal_arcsin,
    );
    simplify(Expr::new_add(term1, term2))
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
    let principal_arccos = Expr::new_arccos(z.clone());
    let term1 = Expr::new_mul(Expr::Constant(2.0), Expr::new_mul(k.clone(), pi));
    let term2 = Expr::new_mul(s.clone(), principal_arccos);
    simplify(Expr::new_add(term1, term2))
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
    let principal_arctan = Expr::new_arctan(z.clone());
    let term1 = Expr::new_mul(k.clone(), pi);
    simplify(Expr::new_add(term1, principal_arctan))
}

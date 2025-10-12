//! # Numerical Special Functions
//!
//! This module provides numerical implementations of various special functions
//! commonly encountered in mathematics, physics, and engineering. These functions
//! are typically provided by external libraries (e.g., `statrs`) and are wrapped
//! here for convenience.

use statrs::function::beta::{beta, ln_beta};
use statrs::function::erf::{erf, erfc};
use statrs::function::gamma::{gamma, ln_gamma};

/// Computes the gamma function, `Γ(x)`.
///
/// The gamma function is an extension of the factorial function to real and complex numbers.
/// For positive integers `n`, `Γ(n) = (n-1)!`.
///
/// # Arguments
/// * `x` - The input value.
///
/// # Returns
/// The numerical value of `Γ(x)`.
pub fn gamma_numerical(x: f64) -> f64 {
    gamma(x)
}

/// Computes the natural logarithm of the gamma function, `ln(Γ(x))`.
///
/// This function is often used to avoid overflow/underflow issues when `Γ(x)` itself is very large or very small.
///
/// # Arguments
/// * `x` - The input value.
///
/// # Returns
/// The numerical value of `ln(Γ(x))`.
pub fn ln_gamma_numerical(x: f64) -> f64 {
    ln_gamma(x)
}

/// Computes the beta function, `B(a, b)`.
///
/// The beta function is closely related to the gamma function: `B(a, b) = Γ(a)Γ(b) / Γ(a+b)`.
/// It appears in probability theory and mathematical physics.
///
/// # Arguments
/// * `a` - The first input value.
/// * `b` - The second input value.
///
/// # Returns
/// The numerical value of `B(a, b)`.
pub fn beta_numerical(a: f64, b: f64) -> f64 {
    beta(a, b)
}

/// Computes the natural logarithm of the beta function, `ln(B(a, b))`.
///
/// # Arguments
/// * `a` - The first input value.
/// * `b` - The second input value.
///
/// # Returns
/// The numerical value of `ln(B(a, b))`.
pub fn ln_beta_numerical(a: f64, b: f64) -> f64 {
    ln_beta(a, b)
}

/// Computes the error function, `erf(x)`.
///
/// The error function is a special function of sigmoid shape that arises in probability,
/// statistics, and partial differential equations.
///
/// # Arguments
/// * `x` - The input value.
///
/// # Returns
/// The numerical value of `erf(x)`.
pub fn erf_numerical(x: f64) -> f64 {
    erf(x)
}

/// Computes the complementary error function, `erfc(x) = 1 - erf(x)`.
///
/// # Arguments
/// * `x` - The input value.
///
/// # Returns
/// The numerical value of `erfc(x)`.
pub fn erfc_numerical(x: f64) -> f64 {
    erfc(x)
}

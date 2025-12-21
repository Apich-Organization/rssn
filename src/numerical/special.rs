use statrs::function::beta::{beta, ln_beta};
use statrs::function::erf::{erf, erfc};
use statrs::function::gamma::{gamma, ln_gamma};
/// Computes the gamma function, Γ(x).
#[must_use]
pub fn gamma_numerical(x: f64) -> f64 {
    gamma(x)
}
/// Computes the natural logarithm of the gamma function, ln(Γ(x)).
#[must_use]
pub fn ln_gamma_numerical(x: f64) -> f64 {
    ln_gamma(x)
}
/// Computes the beta function, B(a, b).
#[must_use]
pub fn beta_numerical(a: f64, b: f64) -> f64 {
    beta(a, b)
}
/// Computes the natural logarithm of the beta function, ln(B(a, b)).
#[must_use]
pub fn ln_beta_numerical(a: f64, b: f64) -> f64 {
    ln_beta(a, b)
}
/// Computes the error function, erf(x).
#[must_use]
pub fn erf_numerical(x: f64) -> f64 {
    erf(x)
}
/// Computes the complementary error function, erfc(x) = 1 - erf(x).
#[must_use]
pub fn erfc_numerical(x: f64) -> f64 {
    erfc(x)
}

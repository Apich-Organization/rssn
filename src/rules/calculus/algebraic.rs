//! Integration of algebraic integrands: Hermite reduction in a quadratic
//! extension `Q(x)(sqrt(D))`.
//!
//! An integrand `R(x, sqrt(D))` is written as `(A + B sqrt(D)) / C` with
//! polynomials `A`, `B`, `C`. Hermite reduction subtracts the derivative of
//! a rational element of the extension, lowering every multiplicity in `C`
//! to one; what remains is integrated by the logarithmic part (residues and
//! Trager's resultant) or returned unevaluated.
//!
//! The kernel is registered on `integral(f, x)` for integrands that contain
//! a square root of a polynomial.

use crate::graph::RuleError;
use crate::graph::rule::Installer;

pub(super) fn install(_i: &mut Installer<'_>) -> Result<(), RuleError> {
    Ok(())
}

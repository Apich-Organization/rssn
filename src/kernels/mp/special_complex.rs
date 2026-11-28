//! Arbitrary-precision complex special functions on [`MpComplex`].

use super::ball::Ball;
use super::complex::MpComplex;
use super::real::Real;
use super::real::Round;

#[path = "bessel.rs"]
mod bessel;

#[path = "dirichlet.rs"]
mod dirichlet;

#[path = "error_functions.rs"]
mod error_functions;

#[path = "gamma_family.rs"]
mod gamma_family;

fn is_zero(value: &Ball) -> bool {
    value.mid.is_zero() && value.rad.is_zero()
}

fn divide_real(
    value: &MpComplex,
    divisor: &Ball,
    prec: u32,
) -> Result<MpComplex, String> {
    Ok(MpComplex::new(
        value.re.div(divisor, prec)?,
        value.im.div(divisor, prec)?,
    ))
}

fn add_radius(
    value: &mut Ball,
    tail: &Real,
    prec: u32,
) -> Result<(), String> {
    value.rad = value
        .rad
        .add(tail, prec, Round::TowardPositive)
        .or_else(|_| value.rad.add(tail, 1, Round::TowardPositive))?;
    Ok(())
}

fn add_error(
    value: &mut MpComplex,
    error: &Real,
    real_axis: bool,
    imaginary_axis: bool,
    prec: u32,
) -> Result<(), String> {
    if !imaginary_axis {
        add_radius(&mut value.re, error, prec)?;
    }
    if !real_axis {
        add_radius(&mut value.im, error, prec)?;
    }
    Ok(())
}

impl MpComplex {
    /// Certified complex error function `erf(z)` over a closed rectangle.
    ///
    /// # Errors
    /// Fails when the series cannot certify the tail within its work limit
    /// or a required exponent is outside the dyadic range.
    pub fn erf(
        &self,
        prec: u32,
    ) -> Result<Self, String> {
        error_functions::evaluate(self, false, prec)
    }

    /// Certified complex complementary error function `erfc(z)`.
    ///
    /// # Errors
    /// Fails when the series cannot certify the tail within its work limit
    /// or a required exponent is outside the dyadic range.
    pub fn erfc(
        &self,
        prec: u32,
    ) -> Result<Self, String> {
        error_functions::evaluate(self, true, prec)
    }

    /// Certified analytic log-gamma, real on the positive real axis.
    ///
    /// # Errors
    /// Fails if the rectangle meets the non-positive real axis or the
    /// requested precision, recurrence shift or exponents exceed the limits.
    pub fn ln_gamma(
        &self,
        prec: u32,
    ) -> Result<Self, String> {
        gamma_family::log_gamma(self, prec)
    }

    /// Certified complex gamma on a closed rectangle away from its poles.
    ///
    /// # Errors
    /// Fails for a rectangle meeting a non-positive integer pole, or when
    /// precision, recurrence shift or exponents exceed the limits.
    pub fn gamma(
        &self,
        prec: u32,
    ) -> Result<Self, String> {
        gamma_family::gamma(self, prec)
    }

    /// Certified complex logarithmic derivative of gamma.
    ///
    /// # Errors
    /// Fails for a rectangle meeting a non-positive integer pole, or when
    /// precision, recurrence shift or exponents exceed the limits.
    pub fn digamma(
        &self,
        prec: u32,
    ) -> Result<Self, String> {
        gamma_family::digamma(self, prec)
    }

    /// Certified principal Bessel J for closed dyadic argument and order rectangles.
    ///
    /// Exact integer orders through magnitude 512 admit the whole argument
    /// plane with abs(Re(z)) + abs(Im(z)) <= 256. General orders require a
    /// nonnegative real part, use the principal argument branch, and support
    /// up to 512 requested bits.
    ///
    /// # Errors
    /// Fails when an operand meets an excluded branch or exceeds the work
    /// limits, or when an intermediate is outside the dyadic exponent range.
    pub fn bessel_j(
        &self,
        nu: &Self,
        prec: u32,
    ) -> Result<Self, String> {
        bessel::evaluate(self, nu, prec)
    }

    /// Certified Riemann zeta enclosure for closed dyadic rectangles.
    ///
    /// # Errors
    /// Fails at a rectangle containing the pole at one, out-of-budget
    /// orders or precision, or exponents outside the dyadic range.
    pub fn riemann_zeta(
        &self,
        prec: u32,
    ) -> Result<Self, String> {
        dirichlet::zeta(self, prec)
    }

    /// Certified complex polylogarithm for arguments inside the open unit disk.
    ///
    /// # Errors
    /// Fails if the rectangle meets the unit circle or exceeds the order,
    /// precision, exponent or series work limits.
    pub fn polylog(
        &self,
        s: &Self,
        prec: u32,
    ) -> Result<Self, String> {
        dirichlet::polylog(self, s, prec)
    }
}

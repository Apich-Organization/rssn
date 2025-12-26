//! Handle-based FFI API for numerical special functions.
//!
//! This module provides C-compatible FFI functions for various special mathematical
//! functions including gamma, beta, error functions, Bessel functions, and combinatorial
//! functions.

use std::os::raw::c_double;

use crate::symbolic::special;

// ============================================================================
// Gamma and Related Functions
// ============================================================================

/// Computes the gamma function Γ(x).
#[no_mangle]

pub extern "C" fn rssn_gamma_numerical(
    x : c_double
) -> c_double {

    special::gamma_numerical(x)
}

/// Computes the natural logarithm of the gamma function ln(Γ(x)).
#[no_mangle]

pub extern "C" fn rssn_ln_gamma_numerical(
    x : c_double
) -> c_double {

    special::ln_gamma_numerical(x)
}

/// Computes the digamma function ψ(x) = d/dx ln(Γ(x)).
#[no_mangle]

pub extern "C" fn rssn_digamma_numerical(
    x : c_double
) -> c_double {

    special::digamma_numerical(x)
}

/// Computes the beta function B(a, b).
#[no_mangle]

pub extern "C" fn rssn_beta_numerical(
    a : c_double,
    b : c_double,
) -> c_double {

    special::beta_numerical(a, b)
}

/// Computes the natural logarithm of the beta function ln(B(a, b)).
#[no_mangle]

pub extern "C" fn rssn_ln_beta_numerical(
    a : c_double,
    b : c_double,
) -> c_double {

    special::ln_beta_numerical(a, b)
}

/// Computes the regularized incomplete beta function Iₓ(a, b).
#[no_mangle]

pub extern "C" fn rssn_regularized_incomplete_beta(
    a : c_double,
    b : c_double,
    x : c_double,
) -> c_double {

    special::regularized_incomplete_beta(
        a, b, x,
    )
}

/// Computes the regularized lower incomplete gamma function P(a, x).
#[no_mangle]

pub extern "C" fn rssn_regularized_gamma_p(
    a : c_double,
    x : c_double,
) -> c_double {

    special::regularized_gamma_p(a, x)
}

/// Computes the regularized upper incomplete gamma function Q(a, x).
#[no_mangle]

pub extern "C" fn rssn_regularized_gamma_q(
    a : c_double,
    x : c_double,
) -> c_double {

    special::regularized_gamma_q(a, x)
}

// ============================================================================
// Error Functions
// ============================================================================

/// Computes the error function erf(x).
#[no_mangle]

pub extern "C" fn rssn_erf_numerical(
    x : c_double
) -> c_double {

    special::erf_numerical(x)
}

/// Computes the complementary error function erfc(x).
#[no_mangle]

pub extern "C" fn rssn_erfc_numerical(
    x : c_double
) -> c_double {

    special::erfc_numerical(x)
}

/// Computes the inverse error function erf⁻¹(x).
#[no_mangle]

pub extern "C" fn rssn_inverse_erf(
    x : c_double
) -> c_double {

    special::inverse_erf(x)
}

/// Computes the inverse complementary error function erfc⁻¹(x).
#[no_mangle]

pub extern "C" fn rssn_inverse_erfc(
    x : c_double
) -> c_double {

    special::inverse_erfc(x)
}

// ============================================================================
// Combinatorial Functions
// ============================================================================

/// Computes the factorial n!.
#[no_mangle]

pub extern "C" fn rssn_factorial(
    n : u64
) -> u64 {

    special::factorial(n)
}

/// Computes the double factorial n!!.
#[no_mangle]

pub extern "C" fn rssn_double_factorial(
    n : u64
) -> u64 {

    special::double_factorial(n)
}

/// Computes the binomial coefficient C(n, k).
#[no_mangle]

pub extern "C" fn rssn_binomial(
    n : u64,
    k : u64,
) -> u64 {

    special::binomial(n, k)
}

/// Computes the rising factorial (Pochhammer symbol) (x)ₙ.
#[no_mangle]

pub extern "C" fn rssn_rising_factorial(
    x : c_double,
    n : u32,
) -> c_double {

    special::rising_factorial(x, n)
}

/// Computes the falling factorial (x)₍ₙ₎.
#[no_mangle]

pub extern "C" fn rssn_falling_factorial(
    x : c_double,
    n : u32,
) -> c_double {

    special::falling_factorial(x, n)
}

/// Computes the natural logarithm of the factorial ln(n!).
#[no_mangle]

pub extern "C" fn rssn_ln_factorial(
    n : u64
) -> c_double {

    special::ln_factorial(n)
}

// ============================================================================
// Bessel Functions
// ============================================================================

/// Computes the Bessel function of the first kind J₀(x).
#[no_mangle]

pub extern "C" fn rssn_bessel_j0(
    x : c_double
) -> c_double {

    special::bessel_j0(x)
}

/// Computes the Bessel function of the first kind J₁(x).
#[no_mangle]

pub extern "C" fn rssn_bessel_j1(
    x : c_double
) -> c_double {

    special::bessel_j1(x)
}

/// Computes the Bessel function of the second kind Y₀(x).
#[no_mangle]

pub extern "C" fn rssn_bessel_y0(
    x : c_double
) -> c_double {

    special::bessel_y0(x)
}

/// Computes the Bessel function of the second kind Y₁(x).
#[no_mangle]

pub extern "C" fn rssn_bessel_y1(
    x : c_double
) -> c_double {

    special::bessel_y1(x)
}

/// Computes the modified Bessel function of the first kind I₀(x).
#[no_mangle]

pub extern "C" fn rssn_bessel_i0(
    x : c_double
) -> c_double {

    special::bessel_i0(x)
}

/// Computes the modified Bessel function of the first kind I₁(x).
#[no_mangle]

pub extern "C" fn rssn_bessel_i1(
    x : c_double
) -> c_double {

    special::bessel_i1(x)
}

/// Computes the modified Bessel function of the second kind K₀(x).
#[no_mangle]

pub extern "C" fn rssn_bessel_k0(
    x : c_double
) -> c_double {

    special::bessel_k0(x)
}

/// Computes the modified Bessel function of the second kind K₁(x).
#[no_mangle]

pub extern "C" fn rssn_bessel_k1(
    x : c_double
) -> c_double {

    special::bessel_k1(x)
}

// ============================================================================
// Other Special Functions
// ============================================================================

/// Computes the normalized sinc function sin(πx)/(πx).
#[no_mangle]

pub extern "C" fn rssn_sinc(
    x : c_double
) -> c_double {

    special::sinc(x)
}

/// Computes the Riemann zeta function ζ(s).
#[no_mangle]

pub extern "C" fn rssn_zeta_numerical(
    s : c_double
) -> c_double {

    special::zeta(s)
}

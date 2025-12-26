//! Handle-based FFI API for numerical special functions.

use crate::numerical::special;

// Gamma functions
#[no_mangle]

pub extern "C" fn rssn_num_special_gamma(x: f64) -> f64 {

    special::gamma_numerical(x)
}

#[no_mangle]

pub extern "C" fn rssn_num_special_ln_gamma(x: f64) -> f64 {

    special::ln_gamma_numerical(x)
}

#[no_mangle]

pub extern "C" fn rssn_num_special_digamma(x: f64) -> f64 {

    special::digamma_numerical(x)
}

#[no_mangle]

pub extern "C" fn rssn_num_special_lower_incomplete_gamma(
    s: f64,
    x: f64,
) -> f64 {

    special::lower_incomplete_gamma(s, x)
}

#[no_mangle]

pub extern "C" fn rssn_num_special_upper_incomplete_gamma(
    s: f64,
    x: f64,
) -> f64 {

    special::upper_incomplete_gamma(s, x)
}

// Beta functions
#[no_mangle]

pub extern "C" fn rssn_num_special_beta(
    a: f64,
    b: f64,
) -> f64 {

    special::beta_numerical(a, b)
}

#[no_mangle]

pub extern "C" fn rssn_num_special_ln_beta(
    a: f64,
    b: f64,
) -> f64 {

    special::ln_beta_numerical(a, b)
}

#[no_mangle]

pub extern "C" fn rssn_num_special_regularized_beta(
    x: f64,
    a: f64,
    b: f64,
) -> f64 {

    special::regularized_beta(x, a, b)
}

// Error functions
#[no_mangle]

pub extern "C" fn rssn_num_special_erf(x: f64) -> f64 {

    special::erf_numerical(x)
}

#[no_mangle]

pub extern "C" fn rssn_num_special_erfc(x: f64) -> f64 {

    special::erfc_numerical(x)
}

#[no_mangle]

pub extern "C" fn rssn_num_special_inverse_erf(x: f64) -> f64 {

    special::inverse_erf_numerical(x)
}

// Bessel functions
#[no_mangle]

pub extern "C" fn rssn_num_special_bessel_j0(x: f64) -> f64 {

    special::bessel_j0(x)
}

#[no_mangle]

pub extern "C" fn rssn_num_special_bessel_j1(x: f64) -> f64 {

    special::bessel_j1(x)
}

#[no_mangle]

pub extern "C" fn rssn_num_special_bessel_y0(x: f64) -> f64 {

    special::bessel_y0(x)
}

#[no_mangle]

pub extern "C" fn rssn_num_special_bessel_y1(x: f64) -> f64 {

    special::bessel_y1(x)
}

#[no_mangle]

pub extern "C" fn rssn_num_special_bessel_i0(x: f64) -> f64 {

    special::bessel_i0(x)
}

#[no_mangle]

pub extern "C" fn rssn_num_special_bessel_i1(x: f64) -> f64 {

    special::bessel_i1(x)
}

// Orthogonal polynomials
#[no_mangle]

pub extern "C" fn rssn_num_special_legendre_p(
    n: u32,
    x: f64,
) -> f64 {

    special::legendre_p(n, x)
}

#[no_mangle]

pub extern "C" fn rssn_num_special_chebyshev_t(
    n: u32,
    x: f64,
) -> f64 {

    special::chebyshev_t(n, x)
}

#[no_mangle]

pub extern "C" fn rssn_num_special_chebyshev_u(
    n: u32,
    x: f64,
) -> f64 {

    special::chebyshev_u(n, x)
}

#[no_mangle]

pub extern "C" fn rssn_num_special_hermite_h(
    n: u32,
    x: f64,
) -> f64 {

    special::hermite_h(n, x)
}

#[no_mangle]

pub extern "C" fn rssn_num_special_laguerre_l(
    n: u32,
    x: f64,
) -> f64 {

    special::laguerre_l(n, x)
}

// Other special functions
#[no_mangle]

pub extern "C" fn rssn_num_special_factorial(n: u64) -> f64 {

    special::factorial(n)
}

#[no_mangle]

pub extern "C" fn rssn_num_special_double_factorial(n: u64) -> f64 {

    special::double_factorial(n)
}

#[no_mangle]

pub extern "C" fn rssn_num_special_binomial(
    n: u64,
    k: u64,
) -> f64 {

    special::binomial(n, k)
}

#[no_mangle]

pub extern "C" fn rssn_num_special_zeta(s: f64) -> f64 {

    special::riemann_zeta(s)
}

#[no_mangle]

pub extern "C" fn rssn_num_special_sinc(x: f64) -> f64 {

    special::sinc(x)
}

#[no_mangle]

pub extern "C" fn rssn_num_special_sigmoid(x: f64) -> f64 {

    special::sigmoid(x)
}

#[no_mangle]

pub extern "C" fn rssn_num_special_softplus(x: f64) -> f64 {

    special::softplus(x)
}

#[no_mangle]

pub extern "C" fn rssn_num_special_logit(p: f64) -> f64 {

    special::logit(p)
}

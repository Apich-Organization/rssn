//! JSON-based FFI API for symbolic special functions.
//!
//! This module provides C-compatible FFI functions using JSON strings
//! for language interoperability.

use std::os::raw::c_char;

use crate::ffi_apis::common::*;
use crate::symbolic::core::Expr;
use crate::symbolic::special_functions;

// ============================================================================
// Gamma and Related Functions
// ============================================================================

/// Computes the symbolic Gamma function Γ(z) via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_gamma(
    arg_json: *const c_char
) -> *mut c_char {

    let arg: Option<Expr> =
        from_json_string(arg_json);

    if let Some(a) = arg {

        to_json_string(
            &special_functions::gamma(
                a,
            ),
        )
    } else {

        std::ptr::null_mut()
    }
}

/// Computes the symbolic log-gamma function ln(Γ(z)) via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_ln_gamma(
    arg_json: *const c_char
) -> *mut c_char {

    let arg: Option<Expr> =
        from_json_string(arg_json);

    if let Some(a) = arg {

        to_json_string(&special_functions::ln_gamma(a))
    } else {

        std::ptr::null_mut()
    }
}

/// Computes the symbolic Beta function B(a, b) via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_beta(
    a_json: *const c_char,
    b_json: *const c_char,
) -> *mut c_char {

    let a: Option<Expr> =
        from_json_string(a_json);

    let b: Option<Expr> =
        from_json_string(b_json);

    if let (Some(val_a), Some(val_b)) =
        (a, b)
    {

        to_json_string(
            &special_functions::beta(
                val_a, val_b,
            ),
        )
    } else {

        std::ptr::null_mut()
    }
}

/// Computes the symbolic Digamma function ψ(z) via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_digamma(
    arg_json: *const c_char
) -> *mut c_char {

    let arg: Option<Expr> =
        from_json_string(arg_json);

    if let Some(a) = arg {

        to_json_string(
            &special_functions::digamma(
                a,
            ),
        )
    } else {

        std::ptr::null_mut()
    }
}

/// Computes the symbolic Polygamma function ψ⁽ⁿ⁾(z) via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_polygamma(
    n_json: *const c_char,
    z_json: *const c_char,
) -> *mut c_char {

    let n: Option<Expr> =
        from_json_string(n_json);

    let z: Option<Expr> =
        from_json_string(z_json);

    if let (Some(n), Some(z)) = (n, z) {

        to_json_string(&special_functions::polygamma(n, z))
    } else {

        std::ptr::null_mut()
    }
}

// ============================================================================
// Error Functions
// ============================================================================

/// Computes the symbolic error function erf(z) via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_erf(
    arg_json: *const c_char
) -> *mut c_char {

    let arg: Option<Expr> =
        from_json_string(arg_json);

    if let Some(a) = arg {

        to_json_string(
            &special_functions::erf(a),
        )
    } else {

        std::ptr::null_mut()
    }
}

/// Computes the symbolic complementary error function erfc(z) via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_erfc(
    arg_json: *const c_char
) -> *mut c_char {

    let arg: Option<Expr> =
        from_json_string(arg_json);

    if let Some(a) = arg {

        to_json_string(
            &special_functions::erfc(a),
        )
    } else {

        std::ptr::null_mut()
    }
}

/// Computes the symbolic imaginary error function erfi(z) via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_erfi(
    arg_json: *const c_char
) -> *mut c_char {

    let arg: Option<Expr> =
        from_json_string(arg_json);

    if let Some(a) = arg {

        to_json_string(
            &special_functions::erfi(a),
        )
    } else {

        std::ptr::null_mut()
    }
}

// ============================================================================
// Riemann Zeta Function
// ============================================================================

/// Computes the symbolic Riemann zeta function ζ(s) via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_zeta(
    arg_json: *const c_char
) -> *mut c_char {

    let arg: Option<Expr> =
        from_json_string(arg_json);

    if let Some(a) = arg {

        to_json_string(
            &special_functions::zeta(a),
        )
    } else {

        std::ptr::null_mut()
    }
}

// ============================================================================
// Bessel Functions
// ============================================================================

/// Computes the symbolic Bessel function J_n(x) via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_bessel_j(
    order_json: *const c_char,
    arg_json: *const c_char,
) -> *mut c_char {

    let order: Option<Expr> =
        from_json_string(order_json);

    let arg: Option<Expr> =
        from_json_string(arg_json);

    if let (Some(o), Some(a)) =
        (order, arg)
    {

        to_json_string(&special_functions::bessel_j(o, a))
    } else {

        std::ptr::null_mut()
    }
}

/// Computes the symbolic Bessel function Y_n(x) via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_bessel_y(
    order_json: *const c_char,
    arg_json: *const c_char,
) -> *mut c_char {

    let order: Option<Expr> =
        from_json_string(order_json);

    let arg: Option<Expr> =
        from_json_string(arg_json);

    if let (Some(o), Some(a)) =
        (order, arg)
    {

        to_json_string(&special_functions::bessel_y(o, a))
    } else {

        std::ptr::null_mut()
    }
}

/// Computes the symbolic modified Bessel function I_n(x) via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_bessel_i(
    order_json: *const c_char,
    arg_json: *const c_char,
) -> *mut c_char {

    let order: Option<Expr> =
        from_json_string(order_json);

    let arg: Option<Expr> =
        from_json_string(arg_json);

    if let (Some(o), Some(a)) =
        (order, arg)
    {

        to_json_string(&special_functions::bessel_i(o, a))
    } else {

        std::ptr::null_mut()
    }
}

/// Computes the symbolic modified Bessel function K_n(x) via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_bessel_k(
    order_json: *const c_char,
    arg_json: *const c_char,
) -> *mut c_char {

    let order: Option<Expr> =
        from_json_string(order_json);

    let arg: Option<Expr> =
        from_json_string(arg_json);

    if let (Some(o), Some(a)) =
        (order, arg)
    {

        to_json_string(&special_functions::bessel_k(o, a))
    } else {

        std::ptr::null_mut()
    }
}

// ============================================================================
// Orthogonal Polynomials
// ============================================================================

/// Computes the symbolic Legendre polynomial P_n(x) via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_legendre_p(
    degree_json: *const c_char,
    arg_json: *const c_char,
) -> *mut c_char {

    let degree: Option<Expr> =
        from_json_string(degree_json);

    let arg: Option<Expr> =
        from_json_string(arg_json);

    if let (Some(d), Some(a)) =
        (degree, arg)
    {

        to_json_string(&special_functions::legendre_p(d, a))
    } else {

        std::ptr::null_mut()
    }
}

/// Computes the symbolic Laguerre polynomial L_n(x) via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_laguerre_l(
    degree_json: *const c_char,
    arg_json: *const c_char,
) -> *mut c_char {

    let degree: Option<Expr> =
        from_json_string(degree_json);

    let arg: Option<Expr> =
        from_json_string(arg_json);

    if let (Some(d), Some(a)) =
        (degree, arg)
    {

        to_json_string(&special_functions::laguerre_l(d, a))
    } else {

        std::ptr::null_mut()
    }
}

/// Computes the symbolic Generalized Laguerre polynomial L_n^α(x) via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_generalized_laguerre(
    n_json: *const c_char,
    alpha_json: *const c_char,
    x_json: *const c_char,
) -> *mut c_char {

    let n: Option<Expr> =
        from_json_string(n_json);

    let alpha: Option<Expr> =
        from_json_string(alpha_json);

    let x: Option<Expr> =
        from_json_string(x_json);

    if let (
        Some(n),
        Some(alpha),
        Some(x),
    ) = (n, alpha, x)
    {

        to_json_string(&special_functions::generalized_laguerre(n, alpha, x))
    } else {

        std::ptr::null_mut()
    }
}

/// Computes the symbolic Hermite polynomial H_n(x) via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_hermite_h(
    degree_json: *const c_char,
    arg_json: *const c_char,
) -> *mut c_char {

    let degree: Option<Expr> =
        from_json_string(degree_json);

    let arg: Option<Expr> =
        from_json_string(arg_json);

    if let (Some(d), Some(a)) =
        (degree, arg)
    {

        to_json_string(&special_functions::hermite_h(d, a))
    } else {

        std::ptr::null_mut()
    }
}

/// Computes the symbolic Chebyshev polynomial T_n(x) via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_chebyshev_t(
    n_json: *const c_char,
    x_json: *const c_char,
) -> *mut c_char {

    let n: Option<Expr> =
        from_json_string(n_json);

    let x: Option<Expr> =
        from_json_string(x_json);

    if let (Some(n), Some(x)) = (n, x) {

        to_json_string(&special_functions::chebyshev_t(n, x))
    } else {

        std::ptr::null_mut()
    }
}

/// Computes the symbolic Chebyshev polynomial U_n(x) via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_chebyshev_u(
    n_json: *const c_char,
    x_json: *const c_char,
) -> *mut c_char {

    let n: Option<Expr> =
        from_json_string(n_json);

    let x: Option<Expr> =
        from_json_string(x_json);

    if let (Some(n), Some(x)) = (n, x) {

        to_json_string(&special_functions::chebyshev_u(n, x))
    } else {

        std::ptr::null_mut()
    }
}

// ============================================================================
// Differential Equations
// ============================================================================

/// Constructs Bessel's differential equation via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_bessel_differential_equation(
    y_json: *const c_char,
    x_json: *const c_char,
    n_json: *const c_char,
) -> *mut c_char {

    let y: Option<Expr> =
        from_json_string(y_json);

    let x: Option<Expr> =
        from_json_string(x_json);

    let n: Option<Expr> =
        from_json_string(n_json);

    if let (Some(y), Some(x), Some(n)) =
        (y, x, n)
    {

        to_json_string(&special_functions::bessel_differential_equation(&y, &x, &n))
    } else {

        std::ptr::null_mut()
    }
}

/// Constructs Legendre's differential equation via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_legendre_differential_equation(
    y_json: *const c_char,
    x_json: *const c_char,
    n_json: *const c_char,
) -> *mut c_char {

    let y: Option<Expr> =
        from_json_string(y_json);

    let x: Option<Expr> =
        from_json_string(x_json);

    let n: Option<Expr> =
        from_json_string(n_json);

    if let (Some(y), Some(x), Some(n)) =
        (y, x, n)
    {

        to_json_string(&special_functions::legendre_differential_equation(&y, &x, &n))
    } else {

        std::ptr::null_mut()
    }
}

/// Constructs Laguerre's differential equation via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_laguerre_differential_equation(
    y_json: *const c_char,
    x_json: *const c_char,
    n_json: *const c_char,
) -> *mut c_char {

    let y: Option<Expr> =
        from_json_string(y_json);

    let x: Option<Expr> =
        from_json_string(x_json);

    let n: Option<Expr> =
        from_json_string(n_json);

    if let (Some(y), Some(x), Some(n)) =
        (y, x, n)
    {

        to_json_string(&special_functions::laguerre_differential_equation(&y, &x, &n))
    } else {

        std::ptr::null_mut()
    }
}

/// Constructs Hermite's differential equation via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_hermite_differential_equation(
    y_json: *const c_char,
    x_json: *const c_char,
    n_json: *const c_char,
) -> *mut c_char {

    let y: Option<Expr> =
        from_json_string(y_json);

    let x: Option<Expr> =
        from_json_string(x_json);

    let n: Option<Expr> =
        from_json_string(n_json);

    if let (Some(y), Some(x), Some(n)) =
        (y, x, n)
    {

        to_json_string(&special_functions::hermite_differential_equation(&y, &x, &n))
    } else {

        std::ptr::null_mut()
    }
}

/// Constructs Chebyshev's differential equation via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_chebyshev_differential_equation(
    y_json: *const c_char,
    x_json: *const c_char,
    n_json: *const c_char,
) -> *mut c_char {

    let y: Option<Expr> =
        from_json_string(y_json);

    let x: Option<Expr> =
        from_json_string(x_json);

    let n: Option<Expr> =
        from_json_string(n_json);

    if let (Some(y), Some(x), Some(n)) =
        (y, x, n)
    {

        to_json_string(&special_functions::chebyshev_differential_equation(&y, &x, &n))
    } else {

        std::ptr::null_mut()
    }
}

// ============================================================================
// Rodrigues Formulas
// ============================================================================

/// Constructs Rodrigues' formula for Legendre polynomials via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_legendre_rodrigues_formula(
    n_json: *const c_char,
    x_json: *const c_char,
) -> *mut c_char {

    let n: Option<Expr> =
        from_json_string(n_json);

    let x: Option<Expr> =
        from_json_string(x_json);

    if let (Some(n), Some(x)) = (n, x) {

        to_json_string(&special_functions::legendre_rodrigues_formula(&n, &x))
    } else {

        std::ptr::null_mut()
    }
}

/// Constructs Rodrigues' formula for Hermite polynomials via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_hermite_rodrigues_formula(
    n_json: *const c_char,
    x_json: *const c_char,
) -> *mut c_char {

    let n: Option<Expr> =
        from_json_string(n_json);

    let x: Option<Expr> =
        from_json_string(x_json);

    if let (Some(n), Some(x)) = (n, x) {

        to_json_string(&special_functions::hermite_rodrigues_formula(&n, &x))
    } else {

        std::ptr::null_mut()
    }
}

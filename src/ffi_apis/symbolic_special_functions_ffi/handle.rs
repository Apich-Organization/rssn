//! Handle-based FFI API for symbolic special functions.
//!
//! This module provides C-compatible FFI functions using raw pointers
//! to Expr objects for maximum performance.

use crate::symbolic::core::Expr;
use crate::symbolic::special_functions;

// ============================================================================
// Gamma and Related Functions
// ============================================================================

/// Computes the symbolic Gamma function Γ(z).
#[no_mangle]

pub unsafe extern "C" fn rssn_gamma(arg: *const Expr) -> *mut Expr {

    if arg.is_null() {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(special_functions::gamma((*arg).clone())))
}

/// Computes the symbolic log-gamma function ln(Γ(z)).
#[no_mangle]

pub unsafe extern "C" fn rssn_ln_gamma(arg: *const Expr) -> *mut Expr {

    if arg.is_null() {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(special_functions::ln_gamma((*arg).clone())))
}

/// Computes the symbolic Beta function B(a, b).
#[no_mangle]

pub unsafe extern "C" fn rssn_beta(
    a: *const Expr,
    b: *const Expr,
) -> *mut Expr {

    if a.is_null() || b.is_null() {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(special_functions::beta(
        (*a).clone(),
        (*b).clone(),
    )))
}

/// Computes the symbolic Digamma function ψ(z).
#[no_mangle]

pub unsafe extern "C" fn rssn_digamma(arg: *const Expr) -> *mut Expr {

    if arg.is_null() {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(special_functions::digamma((*arg).clone())))
}

/// Computes the symbolic Polygamma function ψ⁽ⁿ⁾(z).
#[no_mangle]

pub unsafe extern "C" fn rssn_polygamma(
    n: *const Expr,
    z: *const Expr,
) -> *mut Expr {

    if n.is_null() || z.is_null() {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(special_functions::polygamma(
        (*n).clone(),
        (*z).clone(),
    )))
}

// ============================================================================
// Error Functions
// ============================================================================

/// Computes the symbolic error function erf(z).
#[no_mangle]

pub unsafe extern "C" fn rssn_erf(arg: *const Expr) -> *mut Expr {

    if arg.is_null() {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(special_functions::erf((*arg).clone())))
}

/// Computes the symbolic complementary error function erfc(z).
#[no_mangle]

pub unsafe extern "C" fn rssn_erfc(arg: *const Expr) -> *mut Expr {

    if arg.is_null() {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(special_functions::erfc((*arg).clone())))
}

/// Computes the symbolic imaginary error function erfi(z).
#[no_mangle]

pub unsafe extern "C" fn rssn_erfi(arg: *const Expr) -> *mut Expr {

    if arg.is_null() {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(special_functions::erfi((*arg).clone())))
}

// ============================================================================
// Riemann Zeta Function
// ============================================================================

/// Computes the symbolic Riemann zeta function ζ(s).
#[no_mangle]

pub unsafe extern "C" fn rssn_zeta(arg: *const Expr) -> *mut Expr {

    if arg.is_null() {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(special_functions::zeta((*arg).clone())))
}

// ============================================================================
// Bessel Functions
// ============================================================================

/// Computes the symbolic Bessel function of the first kind J_n(x).
#[no_mangle]

pub unsafe extern "C" fn rssn_bessel_j(
    order: *const Expr,
    arg: *const Expr,
) -> *mut Expr {

    if order.is_null() || arg.is_null() {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(special_functions::bessel_j(
        (*order).clone(),
        (*arg).clone(),
    )))
}

/// Computes the symbolic Bessel function of the second kind Y_n(x).
#[no_mangle]

pub unsafe extern "C" fn rssn_bessel_y(
    order: *const Expr,
    arg: *const Expr,
) -> *mut Expr {

    if order.is_null() || arg.is_null() {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(special_functions::bessel_y(
        (*order).clone(),
        (*arg).clone(),
    )))
}

/// Computes the symbolic modified Bessel function of the first kind I_n(x).
#[no_mangle]

pub unsafe extern "C" fn rssn_bessel_i(
    order: *const Expr,
    arg: *const Expr,
) -> *mut Expr {

    if order.is_null() || arg.is_null() {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(special_functions::bessel_i(
        (*order).clone(),
        (*arg).clone(),
    )))
}

/// Computes the symbolic modified Bessel function of the second kind K_n(x).
#[no_mangle]

pub unsafe extern "C" fn rssn_bessel_k(
    order: *const Expr,
    arg: *const Expr,
) -> *mut Expr {

    if order.is_null() || arg.is_null() {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(special_functions::bessel_k(
        (*order).clone(),
        (*arg).clone(),
    )))
}

// ============================================================================
// Orthogonal Polynomials
// ============================================================================

/// Computes the symbolic Legendre polynomial P_n(x).
#[no_mangle]

pub unsafe extern "C" fn rssn_legendre_p(
    degree: *const Expr,
    arg: *const Expr,
) -> *mut Expr {

    if degree.is_null() || arg.is_null() {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(special_functions::legendre_p(
        (*degree).clone(),
        (*arg).clone(),
    )))
}

/// Computes the symbolic Laguerre polynomial L_n(x).
#[no_mangle]

pub unsafe extern "C" fn rssn_laguerre_l(
    degree: *const Expr,
    arg: *const Expr,
) -> *mut Expr {

    if degree.is_null() || arg.is_null() {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(special_functions::laguerre_l(
        (*degree).clone(),
        (*arg).clone(),
    )))
}

/// Computes the symbolic Generalized Laguerre polynomial L_n^α(x).
#[no_mangle]

pub unsafe extern "C" fn rssn_generalized_laguerre(
    n: *const Expr,
    alpha: *const Expr,
    x: *const Expr,
) -> *mut Expr {

    if n.is_null() || alpha.is_null() || x.is_null() {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(special_functions::generalized_laguerre(
        (*n).clone(),
        (*alpha).clone(),
        (*x).clone(),
    )))
}

/// Computes the symbolic Hermite polynomial H_n(x).
#[no_mangle]

pub unsafe extern "C" fn rssn_hermite_h(
    degree: *const Expr,
    arg: *const Expr,
) -> *mut Expr {

    if degree.is_null() || arg.is_null() {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(special_functions::hermite_h(
        (*degree).clone(),
        (*arg).clone(),
    )))
}

/// Computes the symbolic Chebyshev polynomial of the first kind T_n(x).
#[no_mangle]

pub unsafe extern "C" fn rssn_chebyshev_t(
    n: *const Expr,
    x: *const Expr,
) -> *mut Expr {

    if n.is_null() || x.is_null() {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(special_functions::chebyshev_t(
        (*n).clone(),
        (*x).clone(),
    )))
}

/// Computes the symbolic Chebyshev polynomial of the second kind U_n(x).
#[no_mangle]

pub unsafe extern "C" fn rssn_chebyshev_u(
    n: *const Expr,
    x: *const Expr,
) -> *mut Expr {

    if n.is_null() || x.is_null() {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(special_functions::chebyshev_u(
        (*n).clone(),
        (*x).clone(),
    )))
}

// ============================================================================
// Differential Equations
// ============================================================================

/// Constructs Bessel's differential equation: x²y'' + xy' + (x² - n²)y = 0.
#[no_mangle]

pub unsafe extern "C" fn rssn_bessel_differential_equation(
    y: *const Expr,
    x: *const Expr,
    n: *const Expr,
) -> *mut Expr {

    if y.is_null() || x.is_null() || n.is_null() {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(special_functions::bessel_differential_equation(
        &*y, &*x, &*n,
    )))
}

/// Constructs Legendre's differential equation: (1-x²)y'' - 2xy' + n(n+1)y = 0.
#[no_mangle]

pub unsafe extern "C" fn rssn_legendre_differential_equation(
    y: *const Expr,
    x: *const Expr,
    n: *const Expr,
) -> *mut Expr {

    if y.is_null() || x.is_null() || n.is_null() {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(special_functions::legendre_differential_equation(
        &*y, &*x, &*n,
    )))
}

/// Constructs Laguerre's differential equation: xy'' + (1-x)y' + ny = 0.
#[no_mangle]

pub unsafe extern "C" fn rssn_laguerre_differential_equation(
    y: *const Expr,
    x: *const Expr,
    n: *const Expr,
) -> *mut Expr {

    if y.is_null() || x.is_null() || n.is_null() {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(special_functions::laguerre_differential_equation(
        &*y, &*x, &*n,
    )))
}

/// Constructs Hermite's differential equation: y'' - 2xy' + 2ny = 0.
#[no_mangle]

pub unsafe extern "C" fn rssn_hermite_differential_equation(
    y: *const Expr,
    x: *const Expr,
    n: *const Expr,
) -> *mut Expr {

    if y.is_null() || x.is_null() || n.is_null() {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(special_functions::hermite_differential_equation(
        &*y, &*x, &*n,
    )))
}

/// Constructs Chebyshev's differential equation: (1-x²)y'' - xy' + n²y = 0.
#[no_mangle]

pub unsafe extern "C" fn rssn_chebyshev_differential_equation(
    y: *const Expr,
    x: *const Expr,
    n: *const Expr,
) -> *mut Expr {

    if y.is_null() || x.is_null() || n.is_null() {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(
        special_functions::chebyshev_differential_equation(&*y, &*x, &*n),
    ))
}

// ============================================================================
// Rodrigues Formulas
// ============================================================================

/// Constructs Rodrigues' formula for Legendre polynomials.
#[no_mangle]

pub unsafe extern "C" fn rssn_legendre_rodrigues_formula(
    n: *const Expr,
    x: *const Expr,
) -> *mut Expr {

    if n.is_null() || x.is_null() {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(special_functions::legendre_rodrigues_formula(
        &*n, &*x,
    )))
}

/// Constructs Rodrigues' formula for Hermite polynomials.
#[no_mangle]

pub unsafe extern "C" fn rssn_hermite_rodrigues_formula(
    n: *const Expr,
    x: *const Expr,
) -> *mut Expr {

    if n.is_null() || x.is_null() {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(special_functions::hermite_rodrigues_formula(
        &*n, &*x,
    )))
}

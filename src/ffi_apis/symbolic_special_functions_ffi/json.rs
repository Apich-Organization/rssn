//! JSON-based FFI API for symbolic special functions.
//!
//! This module provides C-compatible FFI functions using JSON strings
//! for language interoperability.

use std::os::raw::c_char;

use crate::ffi_apis::common::from_json_string;
use crate::ffi_apis::common::to_json_string;
use crate::symbolic::core::Expr;
use crate::symbolic::special_functions;

// ============================================================================
// Gamma and Related Functions
// ============================================================================

/// Computes the symbolic Gamma function Γ(z) via JSON interface.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_gamma(
    arg_json: *const c_char
) -> *mut c_char {

    let arg: Option<Expr> =
        from_json_string(arg_json);

    if let Some(a) = arg {

        to_json_string(
            &special_functions::gamma(
                &a,
            ),
        )
    } else {

        std::ptr::null_mut()
    }
}

/// Computes the symbolic log-gamma function ln(Γ(z)) via JSON interface.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

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
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_beta(
    a_json: *const c_char,
    b_json: *const c_char,
) -> *mut c_char {

    let a: Option<Expr> =
        from_json_string(a_json);

    let b: Option<Expr> =
        from_json_string(b_json);

    match (a, b)
    { (Some(val_a), Some(val_b)) => {

        to_json_string(
            &special_functions::beta(
                val_a, val_b,
            ),
        )
    } _ => {

        std::ptr::null_mut()
    }}
}

/// Computes the symbolic Digamma function ψ(z) via JSON interface.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_digamma(
    arg_json: *const c_char
) -> *mut c_char {

    let arg: Option<Expr> =
        from_json_string(arg_json);

    if let Some(a) = arg {

        to_json_string(
            &special_functions::digamma(
                &a,
            ),
        )
    } else {

        std::ptr::null_mut()
    }
}

/// Computes the symbolic Polygamma function ψ⁽ⁿ⁾(z) via JSON interface.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_polygamma(
    n_json: *const c_char,
    z_json: *const c_char,
) -> *mut c_char {

    let n: Option<Expr> =
        from_json_string(n_json);

    let z: Option<Expr> =
        from_json_string(z_json);

    match (n, z) { (Some(n), Some(z)) => {

        to_json_string(&special_functions::polygamma(&n, &z))
    } _ => {

        std::ptr::null_mut()
    }}
}

// ============================================================================
// Error Functions
// ============================================================================

/// Computes the symbolic error function erf(z) via JSON interface.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_erf(
    arg_json: *const c_char
) -> *mut c_char {

    let arg: Option<Expr> =
        from_json_string(arg_json);

    if let Some(a) = arg {

        to_json_string(
            &special_functions::erf(&a),
        )
    } else {

        std::ptr::null_mut()
    }
}

/// Computes the symbolic complementary error function erfc(z) via JSON interface.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

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
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

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
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_zeta(
    arg_json: *const c_char
) -> *mut c_char {

    let arg: Option<Expr> =
        from_json_string(arg_json);

    if let Some(a) = arg {

        to_json_string(
            &special_functions::zeta(&a),
        )
    } else {

        std::ptr::null_mut()
    }
}

// ============================================================================
// Bessel Functions
// ============================================================================

/// Computes the symbolic Bessel function `J_n(x)` via JSON interface.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_bessel_j(
    order_json: *const c_char,
    arg_json: *const c_char,
) -> *mut c_char {

    let order: Option<Expr> =
        from_json_string(order_json);

    let arg: Option<Expr> =
        from_json_string(arg_json);

    match (order, arg)
    { (Some(o), Some(a)) => {

        to_json_string(&special_functions::bessel_j(&o, &a))
    } _ => {

        std::ptr::null_mut()
    }}
}

/// Computes the symbolic Bessel function `Y_n(x)` via JSON interface.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_bessel_y(
    order_json: *const c_char,
    arg_json: *const c_char,
) -> *mut c_char {

    let order: Option<Expr> =
        from_json_string(order_json);

    let arg: Option<Expr> =
        from_json_string(arg_json);

    match (order, arg)
    { (Some(o), Some(a)) => {

        to_json_string(&special_functions::bessel_y(&o, &a))
    } _ => {

        std::ptr::null_mut()
    }}
}

/// Computes the symbolic modified Bessel function `I_n(x)` via JSON interface.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_bessel_i(
    order_json: *const c_char,
    arg_json: *const c_char,
) -> *mut c_char {

    let order: Option<Expr> =
        from_json_string(order_json);

    let arg: Option<Expr> =
        from_json_string(arg_json);

    match (order, arg)
    { (Some(o), Some(a)) => {

        to_json_string(&special_functions::bessel_i(&o, &a))
    } _ => {

        std::ptr::null_mut()
    }}
}

/// Computes the symbolic modified Bessel function `K_n(x)` via JSON interface.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_bessel_k(
    order_json: *const c_char,
    arg_json: *const c_char,
) -> *mut c_char {

    let order: Option<Expr> =
        from_json_string(order_json);

    let arg: Option<Expr> =
        from_json_string(arg_json);

    match (order, arg)
    { (Some(o), Some(a)) => {

        to_json_string(&special_functions::bessel_k(&o, &a))
    } _ => {

        std::ptr::null_mut()
    }}
}

// ============================================================================
// Orthogonal Polynomials
// ============================================================================

/// Computes the symbolic Legendre polynomial `P_n(x)` via JSON interface.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_legendre_p(
    degree_json: *const c_char,
    arg_json: *const c_char,
) -> *mut c_char {

    let degree: Option<Expr> =
        from_json_string(degree_json);

    let arg: Option<Expr> =
        from_json_string(arg_json);

    match (degree, arg)
    { (Some(d), Some(a)) => {

        to_json_string(&special_functions::legendre_p(&d, a))
    } _ => {

        std::ptr::null_mut()
    }}
}

/// Computes the symbolic Laguerre polynomial `L_n(x)` via JSON interface.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_laguerre_l(
    degree_json: *const c_char,
    arg_json: *const c_char,
) -> *mut c_char {

    let degree: Option<Expr> =
        from_json_string(degree_json);

    let arg: Option<Expr> =
        from_json_string(arg_json);

    match (degree, arg)
    { (Some(d), Some(a)) => {

        to_json_string(&special_functions::laguerre_l(&d, a))
    } _ => {

        std::ptr::null_mut()
    }}
}

/// Computes the symbolic Generalized Laguerre polynomial `L_n^α(x)` via JSON interface.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

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

    match (n, alpha, x)
    { (
        Some(n),
        Some(alpha),
        Some(x),
    ) => {

        to_json_string(&special_functions::generalized_laguerre(&n, &alpha, &x))
    } _ => {

        std::ptr::null_mut()
    }}
}

/// Computes the symbolic Hermite polynomial `H_n(x)` via JSON interface.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_hermite_h(
    degree_json: *const c_char,
    arg_json: *const c_char,
) -> *mut c_char {

    let degree: Option<Expr> =
        from_json_string(degree_json);

    let arg: Option<Expr> =
        from_json_string(arg_json);

    match (degree, arg)
    { (Some(d), Some(a)) => {

        to_json_string(&special_functions::hermite_h(&d, a))
    } _ => {

        std::ptr::null_mut()
    }}
}

/// Computes the symbolic Chebyshev polynomial `T_n(x)` via JSON interface.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_chebyshev_t(
    n_json: *const c_char,
    x_json: *const c_char,
) -> *mut c_char {

    let n: Option<Expr> =
        from_json_string(n_json);

    let x: Option<Expr> =
        from_json_string(x_json);

    match (n, x) { (Some(n), Some(x)) => {

        to_json_string(&special_functions::chebyshev_t(&n, &x))
    } _ => {

        std::ptr::null_mut()
    }}
}

/// Computes the symbolic Chebyshev polynomial `U_n(x)` via JSON interface.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_chebyshev_u(
    n_json: *const c_char,
    x_json: *const c_char,
) -> *mut c_char {

    let n: Option<Expr> =
        from_json_string(n_json);

    let x: Option<Expr> =
        from_json_string(x_json);

    match (n, x) { (Some(n), Some(x)) => {

        to_json_string(&special_functions::chebyshev_u(&n, &x))
    } _ => {

        std::ptr::null_mut()
    }}
}

// ============================================================================
// Differential Equations
// ============================================================================

/// Constructs Bessel's differential equation via JSON interface.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

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

    match (y, x, n)
    { (Some(y), Some(x), Some(n)) => {

        to_json_string(&special_functions::bessel_differential_equation(&y, &x, &n))
    } _ => {

        std::ptr::null_mut()
    }}
}

/// Constructs Legendre's differential equation via JSON interface.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

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

    match (y, x, n)
    { (Some(y), Some(x), Some(n)) => {

        to_json_string(&special_functions::legendre_differential_equation(&y, &x, &n))
    } _ => {

        std::ptr::null_mut()
    }}
}

/// Constructs Laguerre's differential equation via JSON interface.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

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

    match (y, x, n)
    { (Some(y), Some(x), Some(n)) => {

        to_json_string(&special_functions::laguerre_differential_equation(&y, &x, &n))
    } _ => {

        std::ptr::null_mut()
    }}
}

/// Constructs Hermite's differential equation via JSON interface.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

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

    match (y, x, n)
    { (Some(y), Some(x), Some(n)) => {

        to_json_string(&special_functions::hermite_differential_equation(&y, &x, &n))
    } _ => {

        std::ptr::null_mut()
    }}
}

/// Constructs Chebyshev's differential equation via JSON interface.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

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

    match (y, x, n)
    { (Some(y), Some(x), Some(n)) => {

        to_json_string(&special_functions::chebyshev_differential_equation(&y, &x, &n))
    } _ => {

        std::ptr::null_mut()
    }}
}

// ============================================================================
// Rodrigues Formulas
// ============================================================================

/// Constructs Rodrigues' formula for Legendre polynomials via JSON interface.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_legendre_rodrigues_formula(
    n_json: *const c_char,
    x_json: *const c_char,
) -> *mut c_char {

    let n: Option<Expr> =
        from_json_string(n_json);

    let x: Option<Expr> =
        from_json_string(x_json);

    match (n, x) { (Some(n), Some(x)) => {

        to_json_string(&special_functions::legendre_rodrigues_formula(&n, &x))
    } _ => {

        std::ptr::null_mut()
    }}
}

/// Constructs Rodrigues' formula for Hermite polynomials via JSON interface.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_hermite_rodrigues_formula(
    n_json: *const c_char,
    x_json: *const c_char,
) -> *mut c_char {

    let n: Option<Expr> =
        from_json_string(n_json);

    let x: Option<Expr> =
        from_json_string(x_json);

    match (n, x) { (Some(n), Some(x)) => {

        to_json_string(&special_functions::hermite_rodrigues_formula(&n, &x))
    } _ => {

        std::ptr::null_mut()
    }}
}

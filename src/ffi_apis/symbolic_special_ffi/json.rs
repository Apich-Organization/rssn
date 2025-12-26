//! JSON-based FFI API for numerical special functions.
//!
//! This module provides JSON string-based FFI functions for various special mathematical
//! functions, enabling language-agnostic integration.

use crate::ffi_apis::common::*;
use crate::symbolic::special;
use std::os::raw::c_char;

// ============================================================================
// Gamma and Related Functions
// ============================================================================

/// Computes the gamma function Γ(x) via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_gamma_numerical(x_json: *const c_char) -> *mut c_char {

    let x: Option<f64> = from_json_string(x_json);

    if let Some(val) = x {

        to_json_string(&special::gamma_numerical(val))
    } else {

        std::ptr::null_mut()
    }
}

/// Computes ln(Γ(x)) via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_ln_gamma_numerical(x_json: *const c_char) -> *mut c_char {

    let x: Option<f64> = from_json_string(x_json);

    if let Some(val) = x {

        to_json_string(&special::ln_gamma_numerical(val))
    } else {

        std::ptr::null_mut()
    }
}

/// Computes the digamma function ψ(x) via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_digamma_numerical(x_json: *const c_char) -> *mut c_char {

    let x: Option<f64> = from_json_string(x_json);

    if let Some(val) = x {

        to_json_string(&special::digamma_numerical(val))
    } else {

        std::ptr::null_mut()
    }
}

/// Computes B(a, b) via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_beta_numerical(
    a_json: *const c_char,
    b_json: *const c_char,
) -> *mut c_char {

    let a: Option<f64> = from_json_string(a_json);

    let b: Option<f64> = from_json_string(b_json);

    if let (Some(val_a), Some(val_b)) = (a, b) {

        to_json_string(&special::beta_numerical(val_a, val_b))
    } else {

        std::ptr::null_mut()
    }
}

/// Computes ln(B(a, b)) via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_ln_beta_numerical(
    a_json: *const c_char,
    b_json: *const c_char,
) -> *mut c_char {

    let a: Option<f64> = from_json_string(a_json);

    let b: Option<f64> = from_json_string(b_json);

    if let (Some(val_a), Some(val_b)) = (a, b) {

        to_json_string(&special::ln_beta_numerical(val_a, val_b))
    } else {

        std::ptr::null_mut()
    }
}

/// Computes the regularized incomplete beta Iₓ(a, b) via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_regularized_incomplete_beta(
    a_json: *const c_char,
    b_json: *const c_char,
    x_json: *const c_char,
) -> *mut c_char {

    let a: Option<f64> = from_json_string(a_json);

    let b: Option<f64> = from_json_string(b_json);

    let x: Option<f64> = from_json_string(x_json);

    if let (Some(va), Some(vb), Some(vx)) = (a, b, x) {

        to_json_string(&special::regularized_incomplete_beta(va, vb, vx))
    } else {

        std::ptr::null_mut()
    }
}

/// Computes P(a, x) via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_regularized_gamma_p(
    a_json: *const c_char,
    x_json: *const c_char,
) -> *mut c_char {

    let a: Option<f64> = from_json_string(a_json);

    let x: Option<f64> = from_json_string(x_json);

    if let (Some(va), Some(vx)) = (a, x) {

        to_json_string(&special::regularized_gamma_p(va, vx))
    } else {

        std::ptr::null_mut()
    }
}

/// Computes Q(a, x) via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_regularized_gamma_q(
    a_json: *const c_char,
    x_json: *const c_char,
) -> *mut c_char {

    let a: Option<f64> = from_json_string(a_json);

    let x: Option<f64> = from_json_string(x_json);

    if let (Some(va), Some(vx)) = (a, x) {

        to_json_string(&special::regularized_gamma_q(va, vx))
    } else {

        std::ptr::null_mut()
    }
}

// ============================================================================
// Error Functions
// ============================================================================

/// Computes erf(x) via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_erf_numerical(x_json: *const c_char) -> *mut c_char {

    let x: Option<f64> = from_json_string(x_json);

    if let Some(val) = x {

        to_json_string(&special::erf_numerical(val))
    } else {

        std::ptr::null_mut()
    }
}

/// Computes erfc(x) via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_erfc_numerical(x_json: *const c_char) -> *mut c_char {

    let x: Option<f64> = from_json_string(x_json);

    if let Some(val) = x {

        to_json_string(&special::erfc_numerical(val))
    } else {

        std::ptr::null_mut()
    }
}

/// Computes erf⁻¹(x) via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_inverse_erf(x_json: *const c_char) -> *mut c_char {

    let x: Option<f64> = from_json_string(x_json);

    if let Some(val) = x {

        to_json_string(&special::inverse_erf(val))
    } else {

        std::ptr::null_mut()
    }
}

/// Computes erfc⁻¹(x) via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_inverse_erfc(x_json: *const c_char) -> *mut c_char {

    let x: Option<f64> = from_json_string(x_json);

    if let Some(val) = x {

        to_json_string(&special::inverse_erfc(val))
    } else {

        std::ptr::null_mut()
    }
}

// ============================================================================
// Combinatorial Functions
// ============================================================================

/// Computes n! via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_factorial(n_json: *const c_char) -> *mut c_char {

    let n: Option<u64> = from_json_string(n_json);

    if let Some(val) = n {

        to_json_string(&special::factorial(
            val,
        ))
    } else {

        std::ptr::null_mut()
    }
}

/// Computes n!! via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_double_factorial(n_json: *const c_char) -> *mut c_char {

    let n: Option<u64> = from_json_string(n_json);

    if let Some(val) = n {

        to_json_string(&special::double_factorial(val))
    } else {

        std::ptr::null_mut()
    }
}

/// Computes C(n, k) via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_binomial(
    n_json: *const c_char,
    k_json: *const c_char,
) -> *mut c_char {

    let n: Option<u64> = from_json_string(n_json);

    let k: Option<u64> = from_json_string(k_json);

    if let (Some(vn), Some(vk)) = (n, k) {

        to_json_string(&special::binomial(
            vn, vk,
        ))
    } else {

        std::ptr::null_mut()
    }
}

/// Computes the rising factorial (x)ₙ via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_rising_factorial(
    x_json: *const c_char,
    n_json: *const c_char,
) -> *mut c_char {

    let x: Option<f64> = from_json_string(x_json);

    let n: Option<u32> = from_json_string(n_json);

    if let (Some(vx), Some(vn)) = (x, n) {

        to_json_string(&special::rising_factorial(vx, vn))
    } else {

        std::ptr::null_mut()
    }
}

/// Computes the falling factorial (x)₍ₙ₎ via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_falling_factorial(
    x_json: *const c_char,
    n_json: *const c_char,
) -> *mut c_char {

    let x: Option<f64> = from_json_string(x_json);

    let n: Option<u32> = from_json_string(n_json);

    if let (Some(vx), Some(vn)) = (x, n) {

        to_json_string(&special::falling_factorial(vx, vn))
    } else {

        std::ptr::null_mut()
    }
}

/// Computes ln(n!) via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_ln_factorial(n_json: *const c_char) -> *mut c_char {

    let n: Option<u64> = from_json_string(n_json);

    if let Some(val) = n {

        to_json_string(&special::ln_factorial(val))
    } else {

        std::ptr::null_mut()
    }
}

// ============================================================================
// Bessel Functions
// ============================================================================

/// Computes J₀(x) via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_bessel_j0(x_json: *const c_char) -> *mut c_char {

    let x: Option<f64> = from_json_string(x_json);

    if let Some(val) = x {

        to_json_string(&special::bessel_j0(
            val,
        ))
    } else {

        std::ptr::null_mut()
    }
}

/// Computes J₁(x) via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_bessel_j1(x_json: *const c_char) -> *mut c_char {

    let x: Option<f64> = from_json_string(x_json);

    if let Some(val) = x {

        to_json_string(&special::bessel_j1(
            val,
        ))
    } else {

        std::ptr::null_mut()
    }
}

/// Computes Y₀(x) via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_bessel_y0(x_json: *const c_char) -> *mut c_char {

    let x: Option<f64> = from_json_string(x_json);

    if let Some(val) = x {

        to_json_string(&special::bessel_y0(
            val,
        ))
    } else {

        std::ptr::null_mut()
    }
}

/// Computes Y₁(x) via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_bessel_y1(x_json: *const c_char) -> *mut c_char {

    let x: Option<f64> = from_json_string(x_json);

    if let Some(val) = x {

        to_json_string(&special::bessel_y1(
            val,
        ))
    } else {

        std::ptr::null_mut()
    }
}

/// Computes I₀(x) via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_bessel_i0(x_json: *const c_char) -> *mut c_char {

    let x: Option<f64> = from_json_string(x_json);

    if let Some(val) = x {

        to_json_string(&special::bessel_i0(
            val,
        ))
    } else {

        std::ptr::null_mut()
    }
}

/// Computes I₁(x) via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_bessel_i1(x_json: *const c_char) -> *mut c_char {

    let x: Option<f64> = from_json_string(x_json);

    if let Some(val) = x {

        to_json_string(&special::bessel_i1(
            val,
        ))
    } else {

        std::ptr::null_mut()
    }
}

/// Computes K₀(x) via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_bessel_k0(x_json: *const c_char) -> *mut c_char {

    let x: Option<f64> = from_json_string(x_json);

    if let Some(val) = x {

        to_json_string(&special::bessel_k0(
            val,
        ))
    } else {

        std::ptr::null_mut()
    }
}

/// Computes K₁(x) via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_bessel_k1(x_json: *const c_char) -> *mut c_char {

    let x: Option<f64> = from_json_string(x_json);

    if let Some(val) = x {

        to_json_string(&special::bessel_k1(
            val,
        ))
    } else {

        std::ptr::null_mut()
    }
}

// ============================================================================
// Other Special Functions
// ============================================================================

/// Computes sinc(x) via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_sinc(x_json: *const c_char) -> *mut c_char {

    let x: Option<f64> = from_json_string(x_json);

    if let Some(val) = x {

        to_json_string(&special::sinc(val))
    } else {

        std::ptr::null_mut()
    }
}

/// Computes ζ(s) via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_zeta_numerical(s_json: *const c_char) -> *mut c_char {

    let s: Option<f64> = from_json_string(s_json);

    if let Some(val) = s {

        to_json_string(&special::zeta(val))
    } else {

        std::ptr::null_mut()
    }
}

//! JSON-based FFI API for symbolic elementary functions.

use std::os::raw::c_char;

use crate::ffi_apis::common::from_json_string;
use crate::ffi_apis::common::to_c_string;
use crate::ffi_apis::common::to_json_string;
use crate::symbolic::core::Expr;
use crate::symbolic::elementary;

/// Creates a sine expression from JSON: sin(expr).
///
/// # Arguments
/// * `json_expr` - JSON-serialized Expr
///
/// # Returns
/// JSON-serialized Expr or null on error
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_sin_json(
    json_expr: *const c_char
) -> *mut c_char {

    let expr : Expr = match from_json_string(json_expr) {
        | Some(e) => e,
        | None => return std::ptr::null_mut(),
    };

    to_json_string(&elementary::sin(
        expr,
    ))
}

/// Creates a cosine expression from JSON: cos(expr).
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_cos_json(
    json_expr: *const c_char
) -> *mut c_char {

    let expr : Expr = match from_json_string(json_expr) {
        | Some(e) => e,
        | None => return std::ptr::null_mut(),
    };

    to_json_string(&elementary::cos(
        expr,
    ))
}

/// Creates a tangent expression from JSON: tan(expr).
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_tan_json(
    json_expr: *const c_char
) -> *mut c_char {

    let expr : Expr = match from_json_string(json_expr) {
        | Some(e) => e,
        | None => return std::ptr::null_mut(),
    };

    to_json_string(&elementary::tan(
        expr,
    ))
}

/// Creates an exponential expression from JSON: e^(expr).
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_exp_json(
    json_expr: *const c_char
) -> *mut c_char {

    let expr : Expr = match from_json_string(json_expr) {
        | Some(e) => e,
        | None => return std::ptr::null_mut(),
    };

    to_json_string(&elementary::exp(
        expr,
    ))
}

/// Creates a natural logarithm expression from JSON: ln(expr).
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_ln_json(
    json_expr: *const c_char
) -> *mut c_char {

    let expr : Expr = match from_json_string(json_expr) {
        | Some(e) => e,
        | None => return std::ptr::null_mut(),
    };

    to_json_string(&elementary::ln(
        expr,
    ))
}

/// Creates a square root expression from JSON: sqrt(expr).
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_sqrt_json(
    json_expr: *const c_char
) -> *mut c_char {

    let expr : Expr = match from_json_string(json_expr) {
        | Some(e) => e,
        | None => return std::ptr::null_mut(),
    };

    to_json_string(&elementary::sqrt(
        expr,
    ))
}

/// Creates a power expression from JSON: base^exp.
///
/// # Arguments
/// * `json_base` - JSON-serialized base Expr
/// * `json_exp` - JSON-serialized exponent Expr
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_pow_json(
    json_base: *const c_char,
    json_exp: *const c_char,
) -> *mut c_char {

    let base : Expr = match from_json_string(json_base) {
        | Some(e) => e,
        | None => return std::ptr::null_mut(),
    };

    let exp : Expr = match from_json_string(json_exp) {
        | Some(e) => e,
        | None => return std::ptr::null_mut(),
    };

    to_json_string(&elementary::pow(
        base, exp,
    ))
}

/// Returns Pi as JSON.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_pi_json()
-> *mut c_char {

    to_json_string(&elementary::pi())
}

/// Returns Euler's number (e) as JSON.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_e_json()
-> *mut c_char {

    to_json_string(&elementary::e())
}

/// Expands a symbolic expression from JSON.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_expand_json(
    json_expr: *const c_char
) -> *mut c_char {

    let expr : Expr = match from_json_string(json_expr) {
        | Some(e) => e,
        | None => return std::ptr::null_mut(),
    };

    to_json_string(&elementary::expand(
        expr,
    ))
}

/// Computes binomial coefficient C(n, k) and returns as JSON string.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_binomial_coefficient_json(
    n: usize,
    k: usize,
) -> *mut c_char {

    let result = elementary::binomial_coefficient(n, k);

    to_c_string(format!(
        "{{\"result\":\"{result}\"}}"
    ))
}

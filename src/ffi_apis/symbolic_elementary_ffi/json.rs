//! JSON-based FFI API for symbolic elementary functions.

use crate::ffi_apis::common::*;
use crate::symbolic::core::Expr;
use crate::symbolic::elementary;
use std::os::raw::c_char;

/// Creates a sine expression from JSON: sin(expr).
/// 
/// # Arguments
/// * `json_expr` - JSON-serialized Expr
/// 
/// # Returns
/// JSON-serialized Expr or null on error
#[no_mangle]
pub unsafe extern "C" fn rssn_sin_json(json_expr: *const c_char) -> *mut c_char {
    let expr: Expr = match from_json_string(json_expr) {
        Some(e) => e,
        None => return std::ptr::null_mut(),
    };
    to_json_string(&elementary::sin(expr))
}

/// Creates a cosine expression from JSON: cos(expr).
#[no_mangle]
pub unsafe extern "C" fn rssn_cos_json(json_expr: *const c_char) -> *mut c_char {
    let expr: Expr = match from_json_string(json_expr) {
        Some(e) => e,
        None => return std::ptr::null_mut(),
    };
    to_json_string(&elementary::cos(expr))
}

/// Creates a tangent expression from JSON: tan(expr).
#[no_mangle]
pub unsafe extern "C" fn rssn_tan_json(json_expr: *const c_char) -> *mut c_char {
    let expr: Expr = match from_json_string(json_expr) {
        Some(e) => e,
        None => return std::ptr::null_mut(),
    };
    to_json_string(&elementary::tan(expr))
}

/// Creates an exponential expression from JSON: e^(expr).
#[no_mangle]
pub unsafe extern "C" fn rssn_exp_json(json_expr: *const c_char) -> *mut c_char {
    let expr: Expr = match from_json_string(json_expr) {
        Some(e) => e,
        None => return std::ptr::null_mut(),
    };
    to_json_string(&elementary::exp(expr))
}

/// Creates a natural logarithm expression from JSON: ln(expr).
#[no_mangle]
pub unsafe extern "C" fn rssn_ln_json(json_expr: *const c_char) -> *mut c_char {
    let expr: Expr = match from_json_string(json_expr) {
        Some(e) => e,
        None => return std::ptr::null_mut(),
    };
    to_json_string(&elementary::ln(expr))
}

/// Creates a square root expression from JSON: sqrt(expr).
#[no_mangle]
pub unsafe extern "C" fn rssn_sqrt_json(json_expr: *const c_char) -> *mut c_char {
    let expr: Expr = match from_json_string(json_expr) {
        Some(e) => e,
        None => return std::ptr::null_mut(),
    };
    to_json_string(&elementary::sqrt(expr))
}

/// Creates a power expression from JSON: base^exp.
/// 
/// # Arguments
/// * `json_base` - JSON-serialized base Expr
/// * `json_exp` - JSON-serialized exponent Expr
#[no_mangle]
pub unsafe extern "C" fn rssn_pow_json(
    json_base: *const c_char,
    json_exp: *const c_char,
) -> *mut c_char {
    let base: Expr = match from_json_string(json_base) {
        Some(e) => e,
        None => return std::ptr::null_mut(),
    };
    let exp: Expr = match from_json_string(json_exp) {
        Some(e) => e,
        None => return std::ptr::null_mut(),
    };
    to_json_string(&elementary::pow(base, exp))
}

/// Returns Pi as JSON.
#[no_mangle]
pub extern "C" fn rssn_pi_json() -> *mut c_char {
    to_json_string(&elementary::pi())
}

/// Returns Euler's number (e) as JSON.
#[no_mangle]
pub extern "C" fn rssn_e_json() -> *mut c_char {
    to_json_string(&elementary::e())
}

/// Expands a symbolic expression from JSON.
#[no_mangle]
pub unsafe extern "C" fn rssn_expand_json(json_expr: *const c_char) -> *mut c_char {
    let expr: Expr = match from_json_string(json_expr) {
        Some(e) => e,
        None => return std::ptr::null_mut(),
    };
    to_json_string(&elementary::expand(expr))
}

/// Computes binomial coefficient C(n, k) and returns as JSON string.
#[no_mangle]
pub extern "C" fn rssn_binomial_coefficient_json(n: usize, k: usize) -> *mut c_char {
    let result = elementary::binomial_coefficient(n, k);
    to_c_string(format!("{{\"result\":\"{}\"}}", result))
}

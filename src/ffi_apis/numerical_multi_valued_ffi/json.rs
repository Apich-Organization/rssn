//! JSON-based FFI API for numerical multi-valued functions.

use crate::ffi_apis::common::{from_json_string, to_c_string};
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::multi_valued;
use crate::symbolic::core::Expr;
use num_complex::Complex;
use serde::{Deserialize, Serialize};
use std::os::raw::c_char;

#[derive(Deserialize)]
struct NewtonInput {
    f: Expr,
    f_prime: Expr,
    start_re: f64,
    start_im: f64,
    tolerance: f64,
    max_iter: usize,
}

#[derive(Serialize)]
struct ComplexResult {
    re: f64,
    im: f64,
}

#[no_mangle]
pub unsafe extern "C" fn rssn_num_mv_newton_method_complex_json(input_json: *const c_char) -> *mut c_char {
    let input: NewtonInput = match from_json_string(input_json) {
        Some(i) => i,
        None => return to_c_string(serde_json::to_string(&FfiResult::<ComplexResult, String> { ok: None, err: Some("Invalid JSON input".to_string()) }).unwrap()),
    };

    let start_point = Complex::new(input.start_re, input.start_im);
    match multi_valued::newton_method_complex(&input.f, &input.f_prime, start_point, input.tolerance, input.max_iter) {
        Some(root) => {
            let res = ComplexResult { re: root.re, im: root.im };
            to_c_string(serde_json::to_string(&FfiResult { ok: Some(res), err: None::<String> }).unwrap())
        }
        None => to_c_string(serde_json::to_string(&FfiResult::<ComplexResult, String> { ok: None, err: Some("Newton's method failed to converge".to_string()) }).unwrap()),
    }
}

#[derive(Deserialize)]
struct LogSqrtInput {
    re: f64,
    im: f64,
    k: i32,
}

#[no_mangle]
pub unsafe extern "C" fn rssn_num_mv_complex_log_k_json(json: *const c_char) -> *mut c_char {
    let input: LogSqrtInput = match from_json_string(json) {
        Some(i) => i,
        None => return to_c_string(serde_json::to_string(&FfiResult::<ComplexResult, String> { ok: None, err: Some("Invalid JSON input".to_string()) }).unwrap()),
    };
    let z = Complex::new(input.re, input.im);
    let res = multi_valued::complex_log_k(z, input.k);
    let out = ComplexResult { re: res.re, im: res.im };
    to_c_string(serde_json::to_string(&FfiResult { ok: Some(out), err: None::<String> }).unwrap())
}

#[no_mangle]
pub unsafe extern "C" fn rssn_num_mv_complex_sqrt_k_json(json: *const c_char) -> *mut c_char {
    let input: LogSqrtInput = match from_json_string(json) {
        Some(i) => i,
        None => return to_c_string(serde_json::to_string(&FfiResult::<ComplexResult, String> { ok: None, err: Some("Invalid JSON input".to_string()) }).unwrap()),
    };
    let z = Complex::new(input.re, input.im);
    let res = multi_valued::complex_sqrt_k(z, input.k);
    let out = ComplexResult { re: res.re, im: res.im };
    to_c_string(serde_json::to_string(&FfiResult { ok: Some(out), err: None::<String> }).unwrap())
}

#[derive(Deserialize)]
struct PowInput {
    z_re: f64,
    z_im: f64,
    w_re: f64,
    w_im: f64,
    k: i32,
}

#[no_mangle]
pub unsafe extern "C" fn rssn_num_mv_complex_pow_k_json(json: *const c_char) -> *mut c_char {
    let input: PowInput = match from_json_string(json) {
        Some(i) => i,
        None => return to_c_string(serde_json::to_string(&FfiResult::<ComplexResult, String> { ok: None, err: Some("Invalid JSON input".to_string()) }).unwrap()),
    };
    let z = Complex::new(input.z_re, input.z_im);
    let w = Complex::new(input.w_re, input.w_im);
    let res = multi_valued::complex_pow_k(z, w, input.k);
    let out = ComplexResult { re: res.re, im: res.im };
    to_c_string(serde_json::to_string(&FfiResult { ok: Some(out), err: None::<String> }).unwrap())
}

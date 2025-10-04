//! FFI API for the rssn library.
//!
//! This module provides a C-compatible foreign function interface (FFI) for interacting
//! with the core data structures and functions of the `rssn` library.
//!
//! The primary design pattern used here is the "handle-based" interface. Instead of
//! exposing complex Rust structs directly, which is unsafe and unstable, we expose

use std::ffi::{CStr, CString};
use std::os::raw::c_char;
use std::ptr;

use serde::{Deserialize, Serialize};

use crate::symbolic::core::Expr;

/// A macro to generate the boilerplate for a handle-based FFI API.
///
/// This macro creates three functions for a given type `$T`:
/// - `_from_json`: Deserializes a JSON string into a `$T` object and returns it as a raw pointer (handle).
/// - `_to_json`: Serializes a `$T` object (given by its handle) into a JSON string.
/// - `_free`: Frees the memory of a `$T` object given by its handle.
macro_rules! impl_handle_api {
    ($T:ty, $from_json:ident, $to_json:ident, $free:ident) => {
        #[no_mangle]
        pub extern "C" fn $from_json(json_ptr: *const c_char) -> *mut $T {
            if json_ptr.is_null() {
                return ptr::null_mut();
            }
            let json_str = unsafe { CStr::from_ptr(json_ptr).to_str().unwrap() };
            let result: Result<$T, _> = serde_json::from_str(json_str);
            match result {
                Ok(obj) => Box::into_raw(Box::new(obj)),
                Err(_) => ptr::null_mut(),
            }
        }

        #[no_mangle]
        pub extern "C" fn $to_json(handle: *mut $T) -> *mut c_char {
            if handle.is_null() {
                return ptr::null_mut();
            }
            let obj = unsafe { &*handle };
            let json_str = serde_json::to_string(obj).unwrap();
            CString::new(json_str).unwrap().into_raw()
        }

        #[no_mangle]
        pub extern "C" fn $free(handle: *mut $T) {
            if !handle.is_null() {
                unsafe {
                    let _ = Box::from_raw(handle);
                }
            }
        }
    };
}

// Implement the handle API for `Expr` with the prefix "expr".
impl_handle_api!(Expr, expr_from_json, expr_to_json, expr_free);

/// Returns the string representation of an `Expr` handle.
///
/// The caller is responsible for freeing the returned string using `free_string`.
#[no_mangle]
pub extern "C" fn expr_to_string(handle: *mut Expr) -> *mut c_char {
    if handle.is_null() {
        return ptr::null_mut();
    }
    let expr = unsafe { &*handle };
    let s = CString::new(expr.to_string()).unwrap();
    s.into_raw()
}

use crate::symbolic::simplify::simplify;
use crate::symbolic::unit_unification::unify_expression;

#[derive(Serialize)]
struct FfiResult<T, E> {
    ok: Option<T>,
    err: Option<E>,
}

/// Simplifies an `Expr` and returns a handle to the new, simplified expression.
///
/// The caller is responsible for freeing the returned handle using `expr_free`.
#[no_mangle]
pub extern "C" fn expr_simplify(handle: *mut Expr) -> *mut Expr {
    if handle.is_null() {
        return ptr::null_mut();
    }
    let expr = unsafe { &*handle };
    let simplified_expr = simplify(expr.clone());
    Box::into_raw(Box::new(simplified_expr))
}

/// Attempts to unify the units within an expression.
///
/// Returns a JSON string representing a `FfiResult` which contains either the
/// new `Expr` object in the `ok` field or an error message in the `err` field.
/// The caller can then use `expr_from_json` to get a handle to the new expression.
/// The caller is responsible for freeing the returned string using `free_string`.
#[no_mangle]
pub extern "C" fn expr_unify_expression(handle: *mut Expr) -> *mut c_char {
    if handle.is_null() {
        let result = FfiResult {
            ok: None::<Expr>,
            err: Some("Null pointer passed to expr_unify_expression".to_string()),
        };
        let json_str = serde_json::to_string(&result).unwrap();
        return CString::new(json_str).unwrap().into_raw();
    }
    let expr = unsafe { &*handle };
    let unification_result = unify_expression(expr);

    let ffi_result = match unification_result {
        Ok(unified_expr) => FfiResult {
            ok: Some(unified_expr),
            err: None,
        },
        Err(e) => FfiResult {
            ok: None,
            err: Some(e),
        },
    };

    let json_str = serde_json::to_string(&ffi_result).unwrap();
    CString::new(json_str).unwrap().into_raw()
}

/// Allocates and returns a test string ("pong") to the caller.
///
/// This function serves as a more advanced health check for the FFI interface.
/// It allows the client to verify two things:
/// 1. That the FFI function can be called successfully.
/// 2. That memory allocated in Rust can be safely passed to and then freed by the client
///    by calling `free_string` on the returned pointer.
///
/// Returns a pointer to a null-terminated C string. The caller is responsible for freeing this string.
#[no_mangle]
pub extern "C" fn rssn_test_string_passing() -> *mut c_char {
    let s = CString::new("pong").unwrap();
    s.into_raw()
}

/// Frees a C string that was allocated by this library.
#[no_mangle]
pub extern "C" fn free_string(s: *mut c_char) {
    if !s.is_null() {
        unsafe {
            let _ = CString::from_raw(s);
        }
    }
}

// ===== Output FFI Functions =====

use crate::output::latex::to_latex;
use crate::output::pretty_print::pretty_print;

/// Converts an expression to a LaTeX string.
///
/// The caller is responsible for freeing the returned string using `free_string`.
#[no_mangle]
pub extern "C" fn expr_to_latex(handle: *mut Expr) -> *mut c_char {
    if handle.is_null() {
        return ptr::null_mut();
    }
    let expr = unsafe { &*handle };
    let latex_str = to_latex(expr);
    CString::new(latex_str).unwrap().into_raw()
}

/// Converts an expression to a formatted, pretty-printed string.
///
/// The caller is responsible for freeing the returned string using `free_string`.
#[no_mangle]
pub extern "C" fn expr_to_pretty_string(handle: *mut Expr) -> *mut c_char {
    if handle.is_null() {
        return ptr::null_mut();
    }
    let expr = unsafe { &*handle };
    let pretty_str = pretty_print(expr);
    CString::new(pretty_str).unwrap().into_raw()
}

// ===== Interpolation FFI Functions =====

// We need to see the definition of Polynomial to serialize it.
// Assuming it has a public `coeffs` field.
use crate::numerical::polynomial::Polynomial;

#[derive(Deserialize)]
struct LagrangeInput {
    points: Vec<(f64, f64)>,
}

#[derive(Serialize)]
struct FfiPolynomial {
    coeffs: Vec<f64>,
}

#[derive(Deserialize)]
struct BezierInput {
    control_points: Vec<Vec<f64>>,
    t: f64,
}

/// Computes a Lagrange interpolating polynomial and returns its coefficients as a JSON string.
#[no_mangle]
pub extern "C" fn interpolate_lagrange(json_ptr: *const c_char) -> *mut c_char {
    if json_ptr.is_null() {
        let result = FfiResult { ok: None::<FfiPolynomial>, err: Some("Null pointer passed to interpolate_lagrange".to_string()) };
        return CString::new(serde_json::to_string(&result).unwrap()).unwrap().into_raw();
    }
    let json_str = unsafe { CStr::from_ptr(json_ptr).to_str().unwrap() };
    let input: Result<LagrangeInput, _> = serde_json::from_str(json_str);
    let ffi_result = match input {
        Ok(input_data) => {
            match interp_module::lagrange_interpolation(&input_data.points) {
                Ok(poly) => {
                    let ffi_poly = FfiPolynomial { coeffs: poly.coeffs };
                    FfiResult { ok: Some(ffi_poly), err: None::<String> }
                }
                Err(e) => FfiResult { ok: None, err: Some(e) },
            }
        }
        Err(e) => FfiResult { ok: None, err: Some(format!("JSON deserialization error: {}", e)) },
    };
    CString::new(serde_json::to_string(&ffi_result).unwrap()).unwrap().into_raw()
}

/// Evaluates a point on a Bézier curve and returns the coordinates as a JSON string.
#[no_mangle]
pub extern "C" fn interpolate_bezier_curve(json_ptr: *const c_char) -> *mut c_char {
    if json_ptr.is_null() {
        let result = FfiResult { ok: None::<Vec<f64>>, err: Some("Null pointer passed to interpolate_bezier_curve".to_string()) };
        return CString::new(serde_json::to_string(&result).unwrap()).unwrap().into_raw();
    }
    let json_str = unsafe { CStr::from_ptr(json_ptr).to_str().unwrap() };
    let input: Result<BezierInput, _> = serde_json::from_str(json_str);
    let ffi_result = match input {
        Ok(input_data) => {
            let point = interp_module::bezier_curve(&input_data.control_points, input_data.t);
            FfiResult { ok: Some(point), err: None::<String> }
        }
        Err(e) => FfiResult { ok: None, err: Some(format!("JSON deserialization error: {}", e)) },
    };
    CString::new(serde_json::to_string(&ffi_result).unwrap()).unwrap().into_raw()
}

// ===== Vector & Combinatorics FFI Functions =====

use crate::numerical::{vector, combinatorics};

// --- Input Structs ---
#[derive(Deserialize)] struct VecInput { v: Vec<f64> }
#[derive(Deserialize)] struct TwoVecInput { v1: Vec<f64>, v2: Vec<f64> }
#[derive(Deserialize)] struct VecScalarInput { v: Vec<f64>, s: f64 }
#[derive(Deserialize)] struct U64Input { n: u64 }
#[derive(Deserialize)] struct TwoU64Input { n: u64, k: u64 }

// --- Macros for FFI Generation ---

macro_rules! impl_ffi_1_vec_in_f64_out {
    ($fn_name:ident, $wrapped_fn:ident) => {
        #[no_mangle]
        pub extern "C" fn $fn_name(json_ptr: *const c_char) -> *mut c_char {
            if json_ptr.is_null() { return CString::new("{\"err\":\"Null pointer passed to function\"}").unwrap().into_raw(); }
            let json_str = unsafe { CStr::from_ptr(json_ptr).to_str().unwrap() };
            let input: Result<VecInput, _> = serde_json::from_str(json_str);
            let ffi_result = match input {
                Ok(input_data) => FfiResult { ok: Some(vector::$wrapped_fn(&input_data.v)), err: None::<String> },
                Err(e) => FfiResult { ok: None, err: Some(format!("JSON error: {}", e)) },
            };
            CString::new(serde_json::to_string(&ffi_result).unwrap()).unwrap().into_raw()
        }
    };
}

macro_rules! impl_ffi_2_vec_in_f64_out {
    ($fn_name:ident, $wrapped_fn:ident) => {
        #[no_mangle]
        pub extern "C" fn $fn_name(json_ptr: *const c_char) -> *mut c_char {
            if json_ptr.is_null() { return CString::new("{\"err\":\"Null pointer passed to function\"}").unwrap().into_raw(); }
            let json_str = unsafe { CStr::from_ptr(json_ptr).to_str().unwrap() };
            let input: Result<TwoVecInput, _> = serde_json::from_str(json_str);
            let ffi_result = match input {
                Ok(input_data) => {
                    match vector::$wrapped_fn(&input_data.v1, &input_data.v2) {
                        Ok(val) => FfiResult { ok: Some(val), err: None },
                        Err(e) => FfiResult { ok: None, err: Some(e) },
                    }
                },
                Err(e) => FfiResult { ok: None, err: Some(format!("JSON error: {}", e)) },
            };
            CString::new(serde_json::to_string(&ffi_result).unwrap()).unwrap().into_raw()
        }
    };
}

macro_rules! impl_ffi_2_vec_in_vec_out {
    ($fn_name:ident, $wrapped_fn:ident) => {
        #[no_mangle]
        pub extern "C" fn $fn_name(json_ptr: *const c_char) -> *mut c_char {
            if json_ptr.is_null() { return CString::new("{\"err\":\"Null pointer passed to function\"}").unwrap().into_raw(); }
            let json_str = unsafe { CStr::from_ptr(json_ptr).to_str().unwrap() };
            let input: Result<TwoVecInput, _> = serde_json::from_str(json_str);
            let ffi_result = match input {
                Ok(input_data) => {
                    match vector::$wrapped_fn(&input_data.v1, &input_data.v2) {
                        Ok(val) => FfiResult { ok: Some(val), err: None },
                        Err(e) => FfiResult { ok: None, err: Some(e) },
                    }
                },
                Err(e) => FfiResult { ok: None, err: Some(format!("JSON error: {}", e)) },
            };
            CString::new(serde_json::to_string(&ffi_result).unwrap()).unwrap().into_raw()
        }
    };
}

macro_rules! impl_ffi_1_u64_in_f64_out {
    ($fn_name:ident, $wrapped_fn:ident) => {
        #[no_mangle]
        pub extern "C" fn $fn_name(json_ptr: *const c_char) -> *mut c_char {
            if json_ptr.is_null() { return CString::new("{\"err\":\"Null pointer passed to function\"}").unwrap().into_raw(); }
            let json_str = unsafe { CStr::from_ptr(json_ptr).to_str().unwrap() };
            let input: Result<U64Input, _> = serde_json::from_str(json_str);
            let ffi_result = match input {
                Ok(input_data) => FfiResult { ok: Some(combinatorics::$wrapped_fn(input_data.n)), err: None::<String> },
                Err(e) => FfiResult { ok: None, err: Some(format!("JSON error: {}", e)) },
            };
            CString::new(serde_json::to_string(&ffi_result).unwrap()).unwrap().into_raw()
        }
    };
}

macro_rules! impl_ffi_2_u64_in_f64_out {
    ($fn_name:ident, $wrapped_fn:ident) => {
        #[no_mangle]
        pub extern "C" fn $fn_name(json_ptr: *const c_char) -> *mut c_char {
            if json_ptr.is_null() { return CString::new("{\"err\":\"Null pointer passed to function\"}").unwrap().into_raw(); }
            let json_str = unsafe { CStr::from_ptr(json_ptr).to_str().unwrap() };
            let input: Result<TwoU64Input, _> = serde_json::from_str(json_str);
            let ffi_result = match input {
                Ok(input_data) => FfiResult { ok: Some(combinatorics::$wrapped_fn(input_data.n, input_data.k)), err: None::<String> },
                Err(e) => FfiResult { ok: None, err: Some(format!("JSON error: {}", e)) },
            };
            CString::new(serde_json::to_string(&ffi_result).unwrap()).unwrap().into_raw()
        }
    };
}

// --- Vector Functions ---
impl_ffi_1_vec_in_f64_out!(vector_norm, norm);
impl_ffi_2_vec_in_f64_out!(vector_dot_product, dot_product);
impl_ffi_2_vec_in_f64_out!(vector_distance, distance);
impl_ffi_2_vec_in_f64_out!(vector_angle, angle);
impl_ffi_2_vec_in_vec_out!(vector_add, vec_add);
impl_ffi_2_vec_in_vec_out!(vector_sub, vec_sub);
impl_ffi_2_vec_in_vec_out!(vector_cross_product, cross_product);

#[no_mangle]
pub extern "C" fn vector_scalar_mul(json_ptr: *const c_char) -> *mut c_char {
    if json_ptr.is_null() { return CString::new("{\"err\":\"Null pointer passed to function\"}").unwrap().into_raw(); }
    let json_str = unsafe { CStr::from_ptr(json_ptr).to_str().unwrap() };
    let input: Result<VecScalarInput, _> = serde_json::from_str(json_str);
    let ffi_result = match input {
        Ok(input_data) => FfiResult { ok: Some(vector::scalar_mul(&input_data.v, input_data.s)), err: None::<String> },
        Err(e) => FfiResult { ok: None, err: Some(format!("JSON error: {}", e)) },
    };
    CString::new(serde_json::to_string(&ffi_result).unwrap()).unwrap().into_raw()
}

// --- Combinatorics Functions ---
impl_ffi_1_u64_in_f64_out!(combinatorics_factorial, factorial);
impl_ffi_2_u64_in_f64_out!(combinatorics_permutations, permutations);
impl_ffi_2_u64_in_f64_out!(combinatorics_combinations, combinations);

// ===== Number Theory FFI Functions =====

use crate::numerical::number_theory as nt;

// --- Input Structs ---
#[derive(Deserialize)] struct TwoU64NtInput { a: u64, b: u64 }
#[derive(Deserialize)] struct ModPowInput { base: u64, exp: u64, modulus: u64 }
#[derive(Deserialize)] struct TwoI64NtInput { a: i64, b: i64 }
#[derive(Deserialize)] struct U64NtInput { n: u64 }

// --- Macros for FFI Generation ---

macro_rules! impl_ffi_2_u64_in_u64_out {
    ($fn_name:ident, $wrapped_fn:ident) => {
        #[no_mangle]
        pub extern "C" fn $fn_name(json_ptr: *const c_char) -> *mut c_char {
            if json_ptr.is_null() { return CString::new("{\"err\":\"Null pointer passed to function\"}").unwrap().into_raw(); }
            let json_str = unsafe { CStr::from_ptr(json_ptr).to_str().unwrap() };
            let input: Result<TwoU64NtInput, _> = serde_json::from_str(json_str);
            let ffi_result = match input {
                Ok(input_data) => FfiResult { ok: Some(nt::$wrapped_fn(input_data.a, input_data.b)), err: None::<String> },
                Err(e) => FfiResult { ok: None, err: Some(format!("JSON error: {}", e)) },
            };
            CString::new(serde_json::to_string(&ffi_result).unwrap()).unwrap().into_raw()
        }
    };
}

macro_rules! impl_ffi_1_u64_in_bool_out {
    ($fn_name:ident, $wrapped_fn:ident) => {
        #[no_mangle]
        pub extern "C" fn $fn_name(json_ptr: *const c_char) -> *mut c_char {
            if json_ptr.is_null() { return CString::new("{\"err\":\"Null pointer passed to function\"}").unwrap().into_raw(); }
            let json_str = unsafe { CStr::from_ptr(json_ptr).to_str().unwrap() };
            let input: Result<U64NtInput, _> = serde_json::from_str(json_str);
            let ffi_result = match input {
                Ok(input_data) => FfiResult { ok: Some(nt::$wrapped_fn(input_data.n)), err: None::<String> },
                Err(e) => FfiResult { ok: None, err: Some(format!("JSON error: {}", e)) },
            };
            CString::new(serde_json::to_string(&ffi_result).unwrap()).unwrap().into_raw()
        }
    };
}

// --- Function Implementations ---

impl_ffi_2_u64_in_u64_out!(nt_gcd, gcd);
impl_ffi_1_u64_in_bool_out!(nt_is_prime_miller_rabin, is_prime_miller_rabin);

#[no_mangle]
pub extern "C" fn nt_mod_pow(json_ptr: *const c_char) -> *mut c_char {
    if json_ptr.is_null() { return CString::new("{\"err\":\"Null pointer passed to function\"}").unwrap().into_raw(); }
    let json_str = unsafe { CStr::from_ptr(json_ptr).to_str().unwrap() };
    let input: Result<ModPowInput, _> = serde_json::from_str(json_str);
    let ffi_result = match input {
        Ok(d) => FfiResult { ok: Some(nt::mod_pow(d.base as u128, d.exp, d.modulus)), err: None::<String> },
        Err(e) => FfiResult { ok: None, err: Some(format!("JSON error: {}", e)) },
    };
    CString::new(serde_json::to_string(&ffi_result).unwrap()).unwrap().into_raw()
}

#[no_mangle]
pub extern "C" fn nt_mod_inverse(json_ptr: *const c_char) -> *mut c_char {
    if json_ptr.is_null() { return CString::new("{\"err\":\"Null pointer passed to function\"}").unwrap().into_raw(); }
    let json_str = unsafe { CStr::from_ptr(json_ptr).to_str().unwrap() };
    let input: Result<TwoI64NtInput, _> = serde_json::from_str(json_str);
    let ffi_result = match input {
        Ok(d) => FfiResult { ok: nt::mod_inverse(d.a, d.b), err: None::<String> },
        Err(e) => FfiResult { ok: None, err: Some(format!("JSON error: {}", e)) },
    };
    CString::new(serde_json::to_string(&ffi_result).unwrap()).unwrap().into_raw()
}

// ===== Special Functions FFI =====

use crate::numerical::special::{self as special_module};

#[derive(Deserialize)]
struct SpecialFunc1Input {
    x: f64,
}

#[derive(Deserialize)]
struct SpecialFunc2Input {
    a: f64,
    b: f64,
}

macro_rules! impl_special_fn_one_arg {
    ($fn_name:ident, $wrapped_fn:ident) => {
        #[no_mangle]
        pub extern "C" fn $fn_name(json_ptr: *const c_char) -> *mut c_char {
            if json_ptr.is_null() {
                let result = FfiResult { ok: None::<f64>, err: Some("Null pointer passed to function".to_string()) };
                return CString::new(serde_json::to_string(&result).unwrap()).unwrap().into_raw();
            }
            let json_str = unsafe { CStr::from_ptr(json_ptr).to_str().unwrap() };
            let input: Result<SpecialFunc1Input, _> = serde_json::from_str(json_str);
            let ffi_result = match input {
                Ok(input_data) => {
                    let result = special_module::$wrapped_fn(input_data.x);
                    FfiResult { ok: Some(result), err: None::<String> }
                }
                Err(e) => FfiResult { ok: None, err: Some(format!("JSON deserialization error: {}", e)) },
            };
            CString::new(serde_json::to_string(&ffi_result).unwrap()).unwrap().into_raw()
        }
    };
}

macro_rules! impl_special_fn_two_args {
    ($fn_name:ident, $wrapped_fn:ident) => {
        #[no_mangle]
        pub extern "C" fn $fn_name(json_ptr: *const c_char) -> *mut c_char {
            if json_ptr.is_null() {
                let result = FfiResult { ok: None::<f64>, err: Some("Null pointer passed to function".to_string()) };
                return CString::new(serde_json::to_string(&result).unwrap()).unwrap().into_raw();
            }
            let json_str = unsafe { CStr::from_ptr(json_ptr).to_str().unwrap() };
            let input: Result<SpecialFunc2Input, _> = serde_json::from_str(json_str);
            let ffi_result = match input {
                Ok(input_data) => {
                    let result = special_module::$wrapped_fn(input_data.a, input_data.b);
                    FfiResult { ok: Some(result), err: None::<String> }
                }
                Err(e) => FfiResult { ok: None, err: Some(format!("JSON deserialization error: {}", e)) },
            };
            CString::new(serde_json::to_string(&ffi_result).unwrap()).unwrap().into_raw()
        }
    };
}

impl_special_fn_one_arg!(special_gamma, gamma_numerical);
impl_special_fn_one_arg!(special_ln_gamma, ln_gamma_numerical);
impl_special_fn_one_arg!(special_erf, erf_numerical);
impl_special_fn_one_arg!(special_erfc, erfc_numerical);

impl_special_fn_two_args!(special_beta, beta_numerical);
impl_special_fn_two_args!(special_ln_beta, ln_beta_numerical);

// ===== Numerical Transforms FFI Functions =====

use crate::numerical::transforms::{fft, ifft};
use num_complex::Complex;

#[derive(Serialize, Deserialize)]
struct TransformsInput {
    data: Vec<Complex<f64>>,
}

/// Computes the Fast Fourier Transform (FFT) of a sequence of complex numbers.
/// The input is a JSON string representing a vector of complex numbers (e.g., `{"data":[{"re":1.0,"im":0.0},...]}
/// Returns the transformed data as a new JSON string.
#[no_mangle]
pub extern "C" fn transforms_fft(json_ptr: *const c_char) -> *mut c_char {
    if json_ptr.is_null() {
        let result = FfiResult { ok: None::<Vec<Complex<f64>>>, err: Some("Null pointer passed to transforms_fft".to_string()) };
        return CString::new(serde_json::to_string(&result).unwrap()).unwrap().into_raw();
    }
    let json_str = unsafe { CStr::from_ptr(json_ptr).to_str().unwrap() };
    let input: Result<TransformsInput, _> = serde_json::from_str(json_str);
    let ffi_result = match input {
        Ok(mut input_data) => {
            fft(&mut input_data.data);
            FfiResult { ok: Some(input_data.data), err: None::<String> }
        }
        Err(e) => FfiResult { ok: None, err: Some(format!("JSON deserialization error: {}", e)) },
    };
    CString::new(serde_json::to_string(&ffi_result).unwrap()).unwrap().into_raw()
}

/// Computes the Inverse Fast Fourier Transform (IFFT) of a sequence of complex numbers.
/// The input is a JSON string representing a vector of complex numbers (e.g., `{"data":[{"re":1.0,"im":0.0},...]}
/// Returns the transformed data as a new JSON string.
#[no_mangle]
pub extern "C" fn transforms_ifft(json_ptr: *const c_char) -> *mut c_char {
    if json_ptr.is_null() {
        let result = FfiResult { ok: None::<Vec<Complex<f64>>>, err: Some("Null pointer passed to transforms_ifft".to_string()) };
        return CString::new(serde_json::to_string(&result).unwrap()).unwrap().into_raw();
    }
    let json_str = unsafe { CStr::from_ptr(json_ptr).to_str().unwrap() };
    let input: Result<TransformsInput, _> = serde_json::from_str(json_str);
    let ffi_result = match input {
        Ok(mut input_data) => {
            ifft(&mut input_data.data);
            FfiResult { ok: Some(input_data.data), err: None::<String> }
        }
        Err(e) => FfiResult { ok: None, err: Some(format!("JSON deserialization error: {}", e)) },
    };
    CString::new(serde_json::to_string(&ffi_result).unwrap()).unwrap().into_raw()
}

// ===== Symbolic Polynomial FFI Functions =====

#[derive(Deserialize)]
struct PolyInput {
    expr: Expr,
    var: String,
}

#[derive(Deserialize)]
struct PolyDivInput {
    n: Expr,
    d: Expr,
    var: String,
}

#[derive(Deserialize)]
struct PolyFromCoeffsInput {
    coeffs: Vec<Expr>,
    var: String,
}

#[no_mangle]
pub extern "C" fn poly_is_polynomial(json_ptr: *const c_char) -> bool {
    if json_ptr.is_null() {
        return false;
    }
    let json_str = unsafe { CStr::from_ptr(json_ptr).to_str().unwrap() };
    let input: Result<PolyInput, _> = serde_json::from_str(json_str);
    match input {
        Ok(poly_input) => poly_module::is_polynomial(&poly_input.expr, &poly_input.var),
        Err(_) => false,
    }
}

#[no_mangle]
pub extern "C" fn poly_degree(json_ptr: *const c_char) -> i64 {
    if json_ptr.is_null() {
        return -1;
    }
    let json_str = unsafe { CStr::from_ptr(json_ptr).to_str().unwrap() };
    let input: Result<PolyInput, _> = serde_json::from_str(json_str);
    match input {
        Ok(poly_input) => poly_module::polynomial_degree(&poly_input.expr, &poly_input.var),
        Err(_) => -1,
    }
}

#[no_mangle]
pub extern "C" fn poly_leading_coefficient(handle: *mut Expr, var_ptr: *const c_char) -> *mut Expr {
    if handle.is_null() || var_ptr.is_null() {
        return ptr::null_mut();
    }
    let expr = unsafe { &*handle };
    let var = unsafe { CStr::from_ptr(var_ptr).to_str().unwrap() };
    let result_expr = poly_module::leading_coefficient(expr, var);
    Box::into_raw(Box::new(result_expr))
}

#[no_mangle]
pub extern "C" fn poly_long_division(json_ptr: *const c_char) -> *mut c_char {
    if json_ptr.is_null() {
        let result = FfiResult { ok: None::<MatrixPair>, err: Some("Null pointer passed to poly_long_division".to_string()) };
        return CString::new(serde_json::to_string(&result).unwrap()).unwrap().into_raw();
    }
    let json_str = unsafe { CStr::from_ptr(json_ptr).to_str().unwrap() };
    let input: Result<PolyDivInput, _> = serde_json::from_str(json_str);
    let ffi_result = match input {
        Ok(div_input) => {
            let (q, r) = poly_module::polynomial_long_division(&div_input.n, &div_input.d, &div_input.var);
            FfiResult { ok: Some(MatrixPair { p1: q, p2: r }), err: None::<String> }
        }
        Err(e) => FfiResult { ok: None, err: Some(format!("JSON deserialization error: {}", e)) },
    };
    CString::new(serde_json::to_string(&ffi_result).unwrap()).unwrap().into_raw()
}

#[no_mangle]
pub extern "C" fn poly_to_coeffs_vec(json_ptr: *const c_char) -> *mut c_char {
    if json_ptr.is_null() {
        let result = FfiResult { ok: None::<Vec<Expr>>, err: Some("Null pointer passed to poly_to_coeffs_vec".to_string()) };
        return CString::new(serde_json::to_string(&result).unwrap()).unwrap().into_raw();
    }
    let json_str = unsafe { CStr::from_ptr(json_ptr).to_str().unwrap() };
    let input: Result<PolyInput, _> = serde_json::from_str(json_str);
    let ffi_result = match input {
        Ok(poly_input) => {
            let coeffs = poly_module::to_polynomial_coeffs_vec(&poly_input.expr, &poly_input.var);
            FfiResult { ok: Some(coeffs), err: None::<String> }
        }
        Err(e) => FfiResult { ok: None, err: Some(format!("JSON deserialization error: {}", e)) },
    };
    CString::new(serde_json::to_string(&ffi_result).unwrap()).unwrap().into_raw()
}

#[no_mangle]
pub extern "C" fn poly_from_coeffs_vec(json_ptr: *const c_char) -> *mut Expr {
    if json_ptr.is_null() {
        return ptr::null_mut();
    }
    let json_str = unsafe { CStr::from_ptr(json_ptr).to_str().unwrap() };
    let input: Result<PolyFromCoeffsInput, _> = serde_json::from_str(json_str);
    match input {
        Ok(input_data) => {
            let result_expr = from_coeffs_to_expr(&input_data.coeffs, &input_data.var);
            Box::into_raw(Box::new(result_expr))
        }
        Err(_) => ptr::null_mut(),
    }
}

// ===== Statistics FFI Functions =====

#[derive(Deserialize)]
struct StatsDataInput {
    data: Vec<f64>,
}

#[derive(Deserialize)]
struct StatsTwoDataInput {
    data1: Vec<f64>,
    data2: Vec<f64>,
}

#[derive(Deserialize)]
struct PercentileInput {
    data: Vec<f64>,
    p: f64,
}

#[derive(Deserialize)]
struct RegressionInput {
    data: Vec<(f64, f64)>,
}

#[derive(Serialize)]
struct RegressionResult {
    slope: f64,
    intercept: f64,
}

macro_rules! impl_stats_fn_single_data {
    ($fn_name:ident, $wrapped_fn:ident) => {
        #[no_mangle]
        pub extern "C" fn $fn_name(json_ptr: *const c_char) -> *mut c_char {
            if json_ptr.is_null() {
                let result = FfiResult { ok: None::<f64>, err: Some("Null pointer passed to function".to_string()) };
                return CString::new(serde_json::to_string(&result).unwrap()).unwrap().into_raw();
            }
            let json_str = unsafe { CStr::from_ptr(json_ptr).to_str().unwrap() };
            let input: Result<StatsDataInput, _> = serde_json::from_str(json_str);
            let ffi_result = match input {
                Ok(mut input_data) => {
                    let result = stats_module::$wrapped_fn(&mut input_data.data);
                    FfiResult { ok: Some(result), err: None::<String> }
                }
                Err(e) => FfiResult { ok: None, err: Some(format!("JSON deserialization error: {}", e)) },
            };
            CString::new(serde_json::to_string(&ffi_result).unwrap()).unwrap().into_raw()
        }
    };
}

impl_stats_fn_single_data!(stats_mean, mean);
impl_stats_fn_single_data!(stats_variance, variance);
impl_stats_fn_single_data!(stats_std_dev, std_dev);
impl_stats_fn_single_data!(stats_median, median);
impl_stats_fn_single_data!(stats_min, min);
impl_stats_fn_single_data!(stats_max, max);
impl_stats_fn_single_data!(stats_skewness, skewness);
impl_stats_fn_single_data!(stats_kurtosis, kurtosis);
impl_stats_fn_single_data!(stats_shannon_entropy, shannon_entropy);

#[no_mangle]
pub extern "C" fn stats_percentile(json_ptr: *const c_char) -> *mut c_char {
    if json_ptr.is_null() {
        let result = FfiResult { ok: None::<f64>, err: Some("Null pointer passed to stats_percentile".to_string()) };
        return CString::new(serde_json::to_string(&result).unwrap()).unwrap().into_raw();
    }
    let json_str = unsafe { CStr::from_ptr(json_ptr).to_str().unwrap() };
    let input: Result<PercentileInput, _> = serde_json::from_str(json_str);
    let ffi_result = match input {
        Ok(mut input_data) => {
            let result = stats_module::percentile(&mut input_data.data, input_data.p);
            FfiResult { ok: Some(result), err: None::<String> }
        }
        Err(e) => FfiResult { ok: None, err: Some(format!("JSON deserialization error: {}", e)) },
    };
    CString::new(serde_json::to_string(&ffi_result).unwrap()).unwrap().into_raw()
}

macro_rules! impl_stats_fn_two_data {
    ($fn_name:ident, $wrapped_fn:ident) => {
        #[no_mangle]
        pub extern "C" fn $fn_name(json_ptr: *const c_char) -> *mut c_char {
            if json_ptr.is_null() {
                let result = FfiResult { ok: None::<f64>, err: Some("Null pointer passed to function".to_string()) };
                return CString::new(serde_json::to_string(&result).unwrap()).unwrap().into_raw();
            }
            let json_str = unsafe { CStr::from_ptr(json_ptr).to_str().unwrap() };
            let input: Result<StatsTwoDataInput, _> = serde_json::from_str(json_str);
            let ffi_result = match input {
                Ok(input_data) => {
                    let result = stats_module::$wrapped_fn(&input_data.data1, &input_data.data2);
                    FfiResult { ok: Some(result), err: None::<String> }
                }
                Err(e) => FfiResult { ok: None, err: Some(format!("JSON deserialization error: {}", e)) },
            };
            CString::new(serde_json::to_string(&ffi_result).unwrap()).unwrap().into_raw()
        }
    };
}

impl_stats_fn_two_data!(stats_covariance, covariance);
impl_stats_fn_two_data!(stats_correlation, correlation);

#[no_mangle]
pub extern "C" fn stats_simple_linear_regression(json_ptr: *const c_char) -> *mut c_char {
    if json_ptr.is_null() {
        let result = FfiResult { ok: None::<RegressionResult>, err: Some("Null pointer passed to stats_simple_linear_regression".to_string()) };
        return CString::new(serde_json::to_string(&result).unwrap()).unwrap().into_raw();
    }
    let json_str = unsafe { CStr::from_ptr(json_ptr).to_str().unwrap() };
    let input: Result<RegressionInput, _> = serde_json::from_str(json_str);
    let ffi_result = match input {
        Ok(input_data) => {
            let (slope, intercept) = simple_linear_regression(&input_data.data);
            let result = RegressionResult { slope, intercept };
            FfiResult { ok: Some(result), err: None::<String> }
        }
        Err(e) => FfiResult { ok: None, err: Some(format!("JSON deserialization error: {}", e)) },
    };
    CString::new(serde_json::to_string(&ffi_result).unwrap()).unwrap().into_raw()
}

// ===================================
// ===== ADDED FFI FUNCTIONS BELOW =====
// ===================================

// --- Imports ---
use crate::symbolic::calculus::{differentiate, substitute, integrate, definite_integrate, limit};
use crate::symbolic::solve::solve;
use crate::symbolic::matrix::{add_matrices, sub_matrices, mul_matrices, transpose_matrix, determinant, inverse_matrix, identity_matrix, scalar_mul_matrix, trace, characteristic_polynomial, rref, null_space, lu_decomposition, eigen_decomposition};
use crate::numerical::calculus::gradient;
use crate::numerical::integrate::{quadrature, QuadratureMethod};
use crate::physics::physics_sm::solve_advection_diffusion_1d;
use crate::symbolic::polynomial::{self as poly_module, from_coeffs_to_expr};
use crate::numerical::stats::{self as stats_module, simple_linear_regression};
use crate::numerical::interpolate::{self as interp_module};


// ===== Symbolic Calculus FFI Functions =====

/// Differentiates an `Expr` and returns a handle to the new, derivative expression.
#[no_mangle]
pub extern "C" fn expr_differentiate(handle: *mut Expr, var_ptr: *const c_char) -> *mut Expr {
    if handle.is_null() || var_ptr.is_null() {
        return ptr::null_mut();
    }
    let expr = unsafe { &*handle };
    let var = unsafe { CStr::from_ptr(var_ptr).to_str().unwrap() };
    let derivative_expr = differentiate(expr, var);
    Box::into_raw(Box::new(derivative_expr))
}

/// Substitutes a variable in an `Expr` with another `Expr` and returns a handle to the new expression.
#[no_mangle]
pub extern "C" fn expr_substitute(handle: *mut Expr, var_ptr: *const c_char, replacement_handle: *mut Expr) -> *mut Expr {
    if handle.is_null() || var_ptr.is_null() || replacement_handle.is_null() {
        return ptr::null_mut();
    }
    let expr = unsafe { &*handle };
    let var = unsafe { CStr::from_ptr(var_ptr).to_str().unwrap() };
    let replacement = unsafe { &*replacement_handle };
    let substituted_expr = substitute(expr, var, replacement);
    Box::into_raw(Box::new(substituted_expr))
}

/// Computes the indefinite integral of an `Expr` and returns a handle to the new expression.
#[no_mangle]
pub extern "C" fn expr_integrate(handle: *mut Expr, var_ptr: *const c_char) -> *mut Expr {
    if handle.is_null() || var_ptr.is_null() {
        return ptr::null_mut();
    }
    let expr = unsafe { &*handle };
    let var = unsafe { CStr::from_ptr(var_ptr).to_str().unwrap() };
    let integral_expr = integrate(expr, var, None, None);
    Box::into_raw(Box::new(integral_expr))
}

/// Computes the definite integral of an `Expr` and returns a handle to the new expression.
#[no_mangle]
pub extern "C" fn expr_definite_integrate(handle: *mut Expr, var_ptr: *const c_char, lower_handle: *mut Expr, upper_handle: *mut Expr) -> *mut Expr {
    if handle.is_null() || var_ptr.is_null() || lower_handle.is_null() || upper_handle.is_null() {
        return ptr::null_mut();
    }
    let expr = unsafe { &*handle };
    let var = unsafe { CStr::from_ptr(var_ptr).to_str().unwrap() };
    let lower = unsafe { &*lower_handle };
    let upper = unsafe { &*upper_handle };
    let integral_expr = definite_integrate(expr, var, lower, upper);
    Box::into_raw(Box::new(integral_expr))
}

/// Computes the limit of an `Expr` and returns a handle to the new expression.
#[no_mangle]
pub extern "C" fn expr_limit(handle: *mut Expr, var_ptr: *const c_char, to_handle: *mut Expr) -> *mut Expr {
    if handle.is_null() || var_ptr.is_null() || to_handle.is_null() {
        return ptr::null_mut();
    }
    let expr = unsafe { &*handle };
    let var = unsafe { CStr::from_ptr(var_ptr).to_str().unwrap() };
    let to = unsafe { &*to_handle };
    let limit_expr = limit(expr, var, to);
    Box::into_raw(Box::new(limit_expr))
}


// ===== Symbolic Solve & Matrix FFI Functions =====

/// Solves an equation for a given variable and returns the solutions as a JSON string.
#[no_mangle]
pub extern "C" fn expr_solve(handle: *mut Expr, var_ptr: *const c_char) -> *mut c_char {
    if handle.is_null() || var_ptr.is_null() {
        let result = FfiResult {
            ok: None::<Vec<Expr>>,
            err: Some("Null pointer passed to expr_solve".to_string()),
        };
        let json_str = serde_json::to_string(&result).unwrap();
        return CString::new(json_str).unwrap().into_raw();
    }
    let expr = unsafe { &*handle };
    let var = unsafe { CStr::from_ptr(var_ptr).to_str().unwrap() };
    let solutions = solve(expr, var);
    let ffi_result = FfiResult {
        ok: Some(solutions),
        err: None::<String>,
    };
    let json_str = serde_json::to_string(&ffi_result).unwrap();
    CString::new(json_str).unwrap().into_raw()
}

/// Adds two matrices and returns a handle to the new matrix expression.
#[no_mangle]
pub extern "C" fn matrix_add(h1: *mut Expr, h2: *mut Expr) -> *mut Expr {
    if h1.is_null() || h2.is_null() {
        return ptr::null_mut();
    }
    let m1 = unsafe { &*h1 };
    let m2 = unsafe { &*h2 };
    Box::into_raw(Box::new(add_matrices(m1, m2)))
}

/// Subtracts the second matrix from the first and returns a handle to the new matrix expression.
#[no_mangle]
pub extern "C" fn matrix_sub(h1: *mut Expr, h2: *mut Expr) -> *mut Expr {
    if h1.is_null() || h2.is_null() {
        return ptr::null_mut();
    }
    let m1 = unsafe { &*h1 };
    let m2 = unsafe { &*h2 };
    Box::into_raw(Box::new(sub_matrices(m1, m2)))
}

/// Multiplies two matrices and returns a handle to the new matrix expression.
#[no_mangle]
pub extern "C" fn matrix_mul(h1: *mut Expr, h2: *mut Expr) -> *mut Expr {
    if h1.is_null() || h2.is_null() {
        return ptr::null_mut();
    }
    let m1 = unsafe { &*h1 };
    let m2 = unsafe { &*h2 };
    Box::into_raw(Box::new(mul_matrices(m1, m2)))
}

/// Transposes a matrix and returns a handle to the new matrix expression.
#[no_mangle]
pub extern "C" fn matrix_transpose(handle: *mut Expr) -> *mut Expr {
    if handle.is_null() {
        return ptr::null_mut();
    }
    let m = unsafe { &*handle };
    Box::into_raw(Box::new(transpose_matrix(m)))
}

/// Computes the determinant of a matrix and returns a handle to the resulting expression.
#[no_mangle]
pub extern "C" fn matrix_determinant(handle: *mut Expr) -> *mut Expr {
    if handle.is_null() {
        return ptr::null_mut();
    }
    let m = unsafe { &*handle };
    Box::into_raw(Box::new(determinant(m)))
}

/// Inverts a matrix and returns a handle to the new matrix expression.
#[no_mangle]
pub extern "C" fn matrix_inverse(handle: *mut Expr) -> *mut Expr {
    if handle.is_null() {
        return ptr::null_mut();
    }
    let m = unsafe { &*handle };
    Box::into_raw(Box::new(inverse_matrix(m)))
}

/// Creates an identity matrix of a given size and returns a handle to it.
#[no_mangle]
pub extern "C" fn matrix_identity(size: usize) -> *mut Expr {
    Box::into_raw(Box::new(identity_matrix(size)))
}

/// Multiplies a matrix by a scalar and returns a handle to the new matrix expression.
#[no_mangle]
pub extern "C" fn matrix_scalar_mul(scalar_handle: *mut Expr, matrix_handle: *mut Expr) -> *mut Expr {
    if scalar_handle.is_null() || matrix_handle.is_null() {
        return ptr::null_mut();
    }
    let scalar = unsafe { &*scalar_handle };
    let matrix = unsafe { &*matrix_handle };
    Box::into_raw(Box::new(scalar_mul_matrix(scalar, matrix)))
}

/// Computes the trace of a matrix and returns the result as a JSON string.
#[no_mangle]
pub extern "C" fn matrix_trace(handle: *mut Expr) -> *mut c_char {
    if handle.is_null() {
        let result = FfiResult { ok: None::<Expr>, err: Some("Null pointer passed to matrix_trace".to_string()) };
        return CString::new(serde_json::to_string(&result).unwrap()).unwrap().into_raw();
    }
    let m = unsafe { &*handle };
    let result = trace(m);
    let ffi_result = match result {
        Ok(value) => FfiResult { ok: Some(value), err: None },
        Err(e) => FfiResult { ok: None, err: Some(e) },
    };
    CString::new(serde_json::to_string(&ffi_result).unwrap()).unwrap().into_raw()
}

/// Computes the characteristic polynomial of a matrix and returns the result as a JSON string.
#[no_mangle]
pub extern "C" fn matrix_characteristic_polynomial(handle: *mut Expr, var_ptr: *const c_char) -> *mut c_char {
    if handle.is_null() || var_ptr.is_null() {
        let result = FfiResult { ok: None::<Expr>, err: Some("Null pointer passed to matrix_characteristic_polynomial".to_string()) };
        return CString::new(serde_json::to_string(&result).unwrap()).unwrap().into_raw();
    }
    let m = unsafe { &*handle };
    let var = unsafe { CStr::from_ptr(var_ptr).to_str().unwrap() };
    let result = characteristic_polynomial(m, var);
    let ffi_result = match result {
        Ok(value) => FfiResult { ok: Some(value), err: None },
        Err(e) => FfiResult { ok: None, err: Some(e) },
    };
    CString::new(serde_json::to_string(&ffi_result).unwrap()).unwrap().into_raw()
}

/// Computes the Reduced Row Echelon Form (RREF) of a matrix and returns the result as a JSON string.
#[no_mangle]
pub extern "C" fn matrix_rref(handle: *mut Expr) -> *mut c_char {
    if handle.is_null() {
        let result = FfiResult { ok: None::<Expr>, err: Some("Null pointer passed to matrix_rref".to_string()) };
        return CString::new(serde_json::to_string(&result).unwrap()).unwrap().into_raw();
    }
    let m = unsafe { &*handle };
    let result = rref(m);
    let ffi_result = match result {
        Ok(value) => FfiResult { ok: Some(value), err: None },
        Err(e) => FfiResult { ok: None, err: Some(e) },
    };
    CString::new(serde_json::to_string(&ffi_result).unwrap()).unwrap().into_raw()
}

/// Computes the null space of a matrix and returns the result as a JSON string.
#[no_mangle]
pub extern "C" fn matrix_null_space(handle: *mut Expr) -> *mut c_char {
    if handle.is_null() {
        let result = FfiResult { ok: None::<Expr>, err: Some("Null pointer passed to matrix_null_space".to_string()) };
        return CString::new(serde_json::to_string(&result).unwrap()).unwrap().into_raw();
    }
    let m = unsafe { &*handle };
    let result = null_space(m);
    let ffi_result = match result {
        Ok(value) => FfiResult { ok: Some(value), err: None },
        Err(e) => FfiResult { ok: None, err: Some(e) },
    };
    CString::new(serde_json::to_string(&ffi_result).unwrap()).unwrap().into_raw()
}

// --- Decomposition Functions ---

#[derive(Serialize)]
struct MatrixPair {
    p1: Expr,
    p2: Expr,
}

/// Computes the LU decomposition of a matrix and returns the L and U matrices as a JSON string.
#[no_mangle]
pub extern "C" fn matrix_lu_decomposition(handle: *mut Expr) -> *mut c_char {
    if handle.is_null() {
        let result = FfiResult { ok: None::<MatrixPair>, err: Some("Null pointer passed to matrix_lu_decomposition".to_string()) };
        return CString::new(serde_json::to_string(&result).unwrap()).unwrap().into_raw();
    }
    let m = unsafe { &*handle };
    let result = lu_decomposition(m);
    let ffi_result = match result {
        Ok((l, u)) => FfiResult { ok: Some(MatrixPair { p1: l, p2: u }), err: None },
        Err(e) => FfiResult { ok: None, err: Some(e) },
    };
    CString::new(serde_json::to_string(&ffi_result).unwrap()).unwrap().into_raw()
}

/// Computes the eigenvalue decomposition of a matrix and returns the eigenvalues and eigenvectors as a JSON string.
#[no_mangle]
pub extern "C" fn matrix_eigen_decomposition(handle: *mut Expr) -> *mut c_char {
    if handle.is_null() {
        let result = FfiResult { ok: None::<MatrixPair>, err: Some("Null pointer passed to matrix_eigen_decomposition".to_string()) };
        return CString::new(serde_json::to_string(&result).unwrap()).unwrap().into_raw();
    }
    let m = unsafe { &*handle };
    let result = eigen_decomposition(m);
    let ffi_result = match result {
        Ok((eigenvalues, eigenvectors)) => FfiResult { ok: Some(MatrixPair { p1: eigenvalues, p2: eigenvectors }), err: None },
        Err(e) => FfiResult { ok: None, err: Some(e) },
    };
    CString::new(serde_json::to_string(&ffi_result).unwrap()).unwrap().into_raw()
}


// ===== Numerical FFI Functions =====

#[derive(Serialize, Deserialize)]
struct GradientInput {
    expr: Expr,
    vars: Vec<String>,
    point: Vec<f64>,
}

#[no_mangle]
pub extern "C" fn numerical_gradient(json_ptr: *const c_char) -> *mut c_char {
    if json_ptr.is_null() {
        let result = FfiResult {
            ok: None::<Vec<f64>>,
            err: Some("Null pointer passed to numerical_gradient".to_string()),
        };
        let json_str = serde_json::to_string(&result).unwrap();
        return CString::new(json_str).unwrap().into_raw();
    }
    let json_str = unsafe { CStr::from_ptr(json_ptr).to_str().unwrap() };
    let input: Result<GradientInput, _> = serde_json::from_str(json_str);
    let ffi_result = match input {
        Ok(grad_input) => {
            let vars_as_str: Vec<&str> = grad_input.vars.iter().map(|s| s.as_str()).collect();
            let grad_result = gradient(&grad_input.expr, &vars_as_str, &grad_input.point);
            match grad_result {
                Ok(grad_vec) => FfiResult {
                    ok: Some(grad_vec),
                    err: None,
                },
                Err(e) => FfiResult {
                    ok: None,
                    err: Some(e),
                },
            }
        }
        Err(e) => FfiResult {
            ok: None,
            err: Some(format!("JSON deserialization error: {}", e)),
        },
    };
    let json_str = serde_json::to_string(&ffi_result).unwrap();
    CString::new(json_str).unwrap().into_raw()
}

#[derive(Deserialize)]
enum FfiQuadratureMethod {
    Trapezoidal,
    Simpson,
}

#[derive(Deserialize)]
struct IntegrationInput {
    expr: Expr,
    var: String,
    start: f64,
    end: f64,
    n_steps: usize,
    method: FfiQuadratureMethod,
}

#[no_mangle]
pub extern "C" fn numerical_integrate(json_ptr: *const c_char) -> *mut c_char {
    if json_ptr.is_null() {
        let result = FfiResult {
            ok: None::<f64>,
            err: Some("Null pointer passed to numerical_integrate".to_string()),
        };
        let json_str = serde_json::to_string(&result).unwrap();
        return CString::new(json_str).unwrap().into_raw();
    }
    let json_str = unsafe { CStr::from_ptr(json_ptr).to_str().unwrap() };
    let input: Result<IntegrationInput, _> = serde_json::from_str(json_str);
    let ffi_result = match input {
        Ok(int_input) => {
            let method = match int_input.method {
                FfiQuadratureMethod::Trapezoidal => QuadratureMethod::Trapezoidal,
                FfiQuadratureMethod::Simpson => QuadratureMethod::Simpson,
            };
            let int_result = quadrature(
                &int_input.expr, 
                &int_input.var, 
                (int_input.start, int_input.end), 
                int_input.n_steps, 
                method
            );
            match int_result {
                Ok(val) => FfiResult {
                    ok: Some(val),
                    err: None,
                },
                Err(e) => FfiResult {
                    ok: None,
                    err: Some(e),
                },
            }
        }
        Err(e) => FfiResult {
            ok: None,
            err: Some(format!("JSON deserialization error: {}", e)),
        },
    };
    let json_str = serde_json::to_string(&ffi_result).unwrap();
    CString::new(json_str).unwrap().into_raw()
}


// ===== Physics FFI Functions =====

#[derive(Deserialize)]
struct AdvectionDiffusion1DInput {
    initial_condition: Vec<f64>,
    dx: f64,
    c: f64, // Advection speed
    d: f64, // Diffusion coefficient
    dt: f64,
    steps: usize,
}

#[no_mangle]
pub extern "C" fn physics_solve_advection_diffusion_1d(json_ptr: *const c_char) -> *mut c_char {
    if json_ptr.is_null() {
        let result = FfiResult {
            ok: None::<Vec<f64>>,
            err: Some("Null pointer passed to physics_solve_advection_diffusion_1d".to_string()),
        };
        let json_str = serde_json::to_string(&result).unwrap();
        return CString::new(json_str).unwrap().into_raw();
    }
    let json_str = unsafe { CStr::from_ptr(json_ptr).to_str().unwrap() };
    let input: Result<AdvectionDiffusion1DInput, _> = serde_json::from_str(json_str);
    let ffi_result = match input {
        Ok(sim_input) => {
            let result_vec = solve_advection_diffusion_1d(
                &sim_input.initial_condition,
                sim_input.dx,
                sim_input.c,
                sim_input.d,
                sim_input.dt,
                sim_input.steps,
            );
            FfiResult {
                ok: Some(result_vec),
                err: None,
            }
        }
        Err(e) => FfiResult {
            ok: None,
            err: Some(format!("JSON deserialization error: {}", e)),
        },
    };
    let json_str = serde_json::to_string(&ffi_result).unwrap();
    CString::new(json_str).unwrap().into_raw()
}

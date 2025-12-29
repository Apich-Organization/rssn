//! JSON-based FFI API for numerical series.

use std::os::raw::c_char;

use serde::Deserialize;

use crate::ffi_apis::common::from_json_string;
use crate::ffi_apis::common::to_c_string;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::series;
use crate::symbolic::core::Expr;

#[derive(Deserialize)]

struct TaylorInput {
    expr: Expr,
    var: String,
    at_point: f64,
    order: usize,
}

#[derive(Deserialize)]

struct SumInput {
    expr: Expr,
    var: String,
    start: i64,
    end: i64,
}

/// Computes Taylor series coefficients for a symbolic expression using JSON serialization.
///
/// Evaluates the derivatives of the expression at a point to obtain Taylor expansion coefficients:
/// f(x) ≈ Σ [fⁿⁿⁿ(a)/n!](x-a)ⁿ for n = 0 to order.
///
/// # Arguments
///
/// * `input_json` - A JSON string pointer containing:
///   - `expr`: Symbolic expression to expand
///   - `var`: Variable name for expansion
///   - `at_point`: Point a around which to expand
///   - `order`: Maximum order of Taylor expansion
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<Vec<f64>, String>` with
/// an array of Taylor coefficients [c₀, c₁, ..., cₙ].
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

/// # Panics
///
/// This function may panic if the FFI input is malformed, null where not expected,
/// or if internal state synchronization fails (e.g., poisoned locks).

pub unsafe extern "C" fn rssn_numerical_taylor_coefficients_json(
    input_json: *const c_char
) -> *mut c_char {

    let input : TaylorInput = match from_json_string(input_json) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<Vec<f64>, String> {
                        ok : None,
                        err : Some("Invalid JSON input".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let res =
        series::taylor_coefficients(
            &input.expr,
            &input.var,
            input.at_point,
            input.order,
        );

    let ffi_res = match res {
        | Ok(v) => {
            FfiResult {
                ok: Some(v),
                err: None,
            }
        },
        | Err(e) => {
            FfiResult {
                ok: None,
                err: Some(e),
            }
        },
    };

    to_c_string(
        serde_json::to_string(&ffi_res)
            .unwrap(),
    )
}

/// Computes the numerical sum of a symbolic series using JSON serialization.
///
/// Evaluates Σ f(var) for var from start to end, where f is a symbolic expression.
///
/// # Arguments
///
/// * `input_json` - A JSON string pointer containing:
///   - `expr`: Symbolic expression to sum
///   - `var`: Summation index variable name
///   - `start`: Lower limit of summation (inclusive)
///   - `end`: Upper limit of summation (inclusive)
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<f64, String>` with
/// the computed sum value.
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

/// # Panics
///
/// This function may panic if the FFI input is malformed, null where not expected,
/// or if internal state synchronization fails (e.g., poisoned locks).

pub unsafe extern "C" fn rssn_numerical_sum_series_json(
    input_json: *const c_char
) -> *mut c_char {

    let input : SumInput = match from_json_string(input_json) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok : None,
                        err : Some("Invalid JSON input".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let res = series::sum_series(
        &input.expr,
        &input.var,
        input.start,
        input.end,
    );

    let ffi_res = match res {
        | Ok(v) => {
            FfiResult {
                ok: Some(v),
                err: None,
            }
        },
        | Err(e) => {
            FfiResult {
                ok: None,
                err: Some(e),
            }
        },
    };

    to_c_string(
        serde_json::to_string(&ffi_res)
            .unwrap(),
    )
}

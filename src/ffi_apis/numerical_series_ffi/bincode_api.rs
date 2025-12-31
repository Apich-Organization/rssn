//! Bincode-based FFI API for numerical series.

use serde::Deserialize;

use crate::ffi_apis::common::BincodeBuffer;
use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
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

/// Computes Taylor series coefficients for a symbolic expression using bincode serialization.
///
/// Evaluates the derivatives of the expression at a point to obtain Taylor expansion coefficients:
/// f(x) ≈ Σ [fⁿⁿⁿ(a)/n!](x-a)ⁿ for n = 0 to order.
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing `TaylorInput` with:
///   - `expr`: Symbolic expression to expand
///   - `var`: Variable name for expansion
///   - `at_point`: Point a around which to expand
///   - `order`: Maximum order of Taylor expansion
///
/// # Returns
///
/// A bincode-encoded buffer containing `FfiResult<Vec<f64>, String>` with either:
/// - `ok`: Array of Taylor coefficients [c₀, c₁, ..., cₙ]
/// - `err`: Error message if computation failed
///
/// # Safety
///
/// This function is unsafe because it receives raw pointers through FFI.
/// The caller must ensure the input buffer contains valid bincode data.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_numerical_taylor_coefficients_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input : TaylorInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(
                &FfiResult::<Vec<f64>, String> {
                    ok : None,
                    err : Some("Invalid Bincode input".to_string()),
                },
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

    to_bincode_buffer(&ffi_res)
}

/// Computes the numerical sum of a symbolic series using bincode serialization.
///
/// Evaluates Σ f(var) for var from start to end, where f is a symbolic expression.
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing `SumInput` with:
///   - `expr`: Symbolic expression to sum
///   - `var`: Summation index variable name
///   - `start`: Lower limit of summation (inclusive)
///   - `end`: Upper limit of summation (inclusive)
///
/// # Returns
///
/// A bincode-encoded buffer containing `FfiResult<f64, String>` with either:
/// - `ok`: The computed sum value
/// - `err`: Error message if computation failed
///
/// # Safety
///
/// This function is unsafe because it receives raw pointers through FFI.
/// The caller must ensure the input buffer contains valid bincode data.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_numerical_sum_series_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input : SumInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(
                &FfiResult::<f64, String> {
                    ok : None,
                    err : Some("Invalid Bincode input".to_string()),
                },
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

    to_bincode_buffer(&ffi_res)
}

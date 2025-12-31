//! Bincode-based FFI API for numerical complex analysis.

use std::collections::HashMap;

use num_complex::Complex;
use serde::Deserialize;

use crate::ffi_apis::common::BincodeBuffer;
use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::complex_analysis;
use crate::symbolic::core::Expr;

#[derive(Deserialize)]

struct EvalInput {
    expr: Expr,
    vars: HashMap<String, Complex<f64>>,
}

/// Evaluates a complex expression using bincode for serialization.

#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_complex_eval_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input : EvalInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(
                &FfiResult::<Complex<f64>, String> {
                    ok : None,
                    err : Some("Invalid Bincode input".to_string()),
                },
            )
        },
    };

    match complex_analysis::eval_complex_expr(
        &input.expr,
        &input.vars,
    ) {
        | Ok(res) => {
            to_bincode_buffer(&FfiResult {
                ok : Some(res),
                err : None::<String>,
            })
        },
        | Err(e) => {
            to_bincode_buffer(
                &FfiResult::<Complex<f64>, String> {
                    ok : None,
                    err : Some(e),
                },
            )
        },
    }
}

#[derive(Deserialize)]

struct ContourInput {
    expr: Expr,
    var: String,
    path: Vec<Complex<f64>>,
}

/// Computes the contour integral of a complex expression using bincode for serialization.

#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_complex_contour_integral_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input : ContourInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(
                &FfiResult::<Complex<f64>, String> {
                    ok : None,
                    err : Some("Invalid Bincode input".to_string()),
                },
            )
        },
    };

    match complex_analysis::contour_integral_expr(
        &input.expr,
        &input.var,
        &input.path,
    ) {
        | Ok(res) => {
            to_bincode_buffer(&FfiResult {
                ok : Some(res),
                err : None::<String>,
            })
        },
        | Err(e) => {
            to_bincode_buffer(
                &FfiResult::<Complex<f64>, String> {
                    ok : None,
                    err : Some(e),
                },
            )
        },
    }
}

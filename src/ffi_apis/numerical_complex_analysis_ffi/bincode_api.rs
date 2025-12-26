//! Bincode-based FFI API for numerical complex analysis.

use crate::ffi_apis::common::{from_bincode_buffer, to_bincode_buffer, BincodeBuffer};
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::complex_analysis;
use crate::symbolic::core::Expr;
use num_complex::Complex;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

#[derive(Deserialize)]
struct EvalInput {
    expr: Expr,
    vars: HashMap<String, Complex<f64>>,
}

#[no_mangle]
pub unsafe extern "C" fn rssn_num_complex_eval_bincode(buffer: BincodeBuffer) -> BincodeBuffer {
    let input: EvalInput = match from_bincode_buffer(&buffer) {
        Some(i) => i,
        None => {
            return to_bincode_buffer(&FfiResult::<Complex<f64>, String> {
                ok: None,
                err: Some("Invalid Bincode input".to_string()),
            })
        }
    };

    match complex_analysis::eval_complex_expr(&input.expr, &input.vars) {
        Ok(res) => to_bincode_buffer(&FfiResult {
            ok: Some(res),
            err: None::<String>,
        }),
        Err(e) => to_bincode_buffer(&FfiResult::<Complex<f64>, String> {
            ok: None,
            err: Some(e),
        }),
    }
}

#[derive(Deserialize)]
struct ContourInput {
    expr: Expr,
    var: String,
    path: Vec<Complex<f64>>,
}

#[no_mangle]
pub unsafe extern "C" fn rssn_num_complex_contour_integral_bincode(
    buffer: BincodeBuffer,
) -> BincodeBuffer {
    let input: ContourInput = match from_bincode_buffer(&buffer) {
        Some(i) => i,
        None => {
            return to_bincode_buffer(&FfiResult::<Complex<f64>, String> {
                ok: None,
                err: Some("Invalid Bincode input".to_string()),
            })
        }
    };

    match complex_analysis::contour_integral_expr(&input.expr, &input.var, &input.path) {
        Ok(res) => to_bincode_buffer(&FfiResult {
            ok: Some(res),
            err: None::<String>,
        }),
        Err(e) => to_bincode_buffer(&FfiResult::<Complex<f64>, String> {
            ok: None,
            err: Some(e),
        }),
    }
}

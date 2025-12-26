//! Bincode-based FFI API for numerical series.

use crate::ffi_apis::common::{from_bincode_buffer, to_bincode_buffer, BincodeBuffer};
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::series;
use crate::symbolic::core::Expr;
use serde::{Deserialize, Serialize};

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

#[no_mangle]

pub unsafe extern "C" fn rssn_numerical_taylor_coefficients_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input: TaylorInput = match from_bincode_buffer(&buffer) {
        Some(i) => i,
        None => {
            return to_bincode_buffer(
                &FfiResult::<Vec<f64>, String> {
                    ok: None,
                    err: Some("Invalid Bincode input".to_string()),
                },
            )
        }
    };

    let res = series::taylor_coefficients(
        &input.expr,
        &input.var,
        input.at_point,
        input.order,
    );

    let ffi_res = match res {
        Ok(v) => FfiResult {
            ok: Some(v),
            err: None,
        },
        Err(e) => FfiResult {
            ok: None,
            err: Some(e),
        },
    };

    to_bincode_buffer(&ffi_res)
}

#[no_mangle]

pub unsafe extern "C" fn rssn_numerical_sum_series_bincode(buffer: BincodeBuffer) -> BincodeBuffer {

    let input: SumInput = match from_bincode_buffer(&buffer) {
        Some(i) => i,
        None => {
            return to_bincode_buffer(
                &FfiResult::<f64, String> {
                    ok: None,
                    err: Some("Invalid Bincode input".to_string()),
                },
            )
        }
    };

    let res = series::sum_series(
        &input.expr,
        &input.var,
        input.start,
        input.end,
    );

    let ffi_res = match res {
        Ok(v) => FfiResult {
            ok: Some(v),
            err: None,
        },
        Err(e) => FfiResult {
            ok: None,
            err: Some(e),
        },
    };

    to_bincode_buffer(&ffi_res)
}

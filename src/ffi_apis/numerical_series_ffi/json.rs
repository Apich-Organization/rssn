//! JSON-based FFI API for numerical series.

use std::os::raw::c_char;

use serde::Deserialize;
use serde::Serialize;

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

#[no_mangle]

pub unsafe extern "C" fn rssn_numerical_taylor_coefficients_json(
    input_json: *const c_char
) -> *mut c_char {

    let input: TaylorInput =
        match from_json_string(input_json) {
            | Some(i) => i,
            | None => {
                return to_c_string(
                    serde_json::to_string(&FfiResult::<
                        Vec<f64>,
                        String,
                    > {
                        ok: None,
                        err: Some(
                            "Invalid JSON input"
                                .to_string(),
                        ),
                    })
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

#[no_mangle]

pub unsafe extern "C" fn rssn_numerical_sum_series_json(
    input_json: *const c_char
) -> *mut c_char {

    let input: SumInput = match from_json_string(input_json)
    {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(&FfiResult::<
                    f64,
                    String,
                > {
                    ok: None,
                    err: Some(
                        "Invalid JSON input".to_string(),
                    ),
                })
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

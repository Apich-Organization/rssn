//! JSON-based FFI API for numerical complex analysis.

use std::collections::HashMap;
use std::os::raw::c_char;

use num_complex::Complex;
use serde::Deserialize;
use serde::Serialize;

use crate::ffi_apis::common::from_json_string;
use crate::ffi_apis::common::to_c_string;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::complex_analysis;
use crate::symbolic::core::Expr;

#[derive(Deserialize)]

struct EvalInput {
    expr : Expr,
    vars :
        HashMap<String, Complex<f64>>,
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_complex_eval_json(
    input_json : *const c_char
) -> *mut c_char {

    let input: EvalInput =
        match from_json_string(input_json) {
            | Some(i) => i,
            | None => {
                return to_c_string(
                    serde_json::to_string(&FfiResult::<
                        Complex<f64>,
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

    match complex_analysis::eval_complex_expr(
        &input.expr,
        &input.vars,
    ) {
        | Ok(res) => {
            to_c_string(
                serde_json::to_string(&FfiResult {
                    ok: Some(res),
                    err: None::<String>,
                })
                .unwrap(),
            )
        },
        | Err(e) => {
            to_c_string(
                serde_json::to_string(&FfiResult::<
                    Complex<f64>,
                    String,
                > {
                    ok: None,
                    err: Some(e),
                })
                .unwrap(),
            )
        },
    }
}

#[derive(Deserialize)]

struct ContourInput {
    expr : Expr,
    var : String,
    path : Vec<Complex<f64>>,
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_complex_contour_integral_json(
    input_json : *const c_char
) -> *mut c_char {

    let input: ContourInput =
        match from_json_string(input_json) {
            | Some(i) => i,
            | None => {
                return to_c_string(
                    serde_json::to_string(&FfiResult::<
                        Complex<f64>,
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

    match complex_analysis::contour_integral_expr(
        &input.expr,
        &input.var,
        &input.path,
    ) {
        | Ok(res) => {
            to_c_string(
                serde_json::to_string(&FfiResult {
                    ok: Some(res),
                    err: None::<String>,
                })
                .unwrap(),
            )
        },
        | Err(e) => {
            to_c_string(
                serde_json::to_string(&FfiResult::<
                    Complex<f64>,
                    String,
                > {
                    ok: None,
                    err: Some(e),
                })
                .unwrap(),
            )
        },
    }
}

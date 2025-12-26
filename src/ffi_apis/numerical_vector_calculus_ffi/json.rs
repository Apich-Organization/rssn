//! JSON-based FFI API for numerical vector calculus.

use std::os::raw::c_char;

use serde::Deserialize;
use serde::Serialize;

use crate::ffi_apis::common::from_json_string;
use crate::ffi_apis::common::to_c_string;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::vector_calculus;
use crate::symbolic::core::Expr;

#[derive(Deserialize)]

struct DivergenceInput {
    funcs : Vec<Expr>,
    vars : Vec<String>,
    point : Vec<f64>,
}

#[derive(Deserialize)]

struct CurlInput {
    funcs : Vec<Expr>,
    vars : Vec<String>,
    point : Vec<f64>,
}

#[derive(Deserialize)]

struct LaplacianInput {
    f : Expr,
    vars : Vec<String>,
    point : Vec<f64>,
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_vector_calculus_divergence_json(
    input_json : *const c_char
) -> *mut c_char {

    let input : DivergenceInput = match from_json_string(input_json) {
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

    let vars_refs : Vec<&str> = input
        .vars
        .iter()
        .map(|s| s.as_str())
        .collect();

    let res = vector_calculus::divergence_expr(
        &input.funcs,
        &vars_refs,
        &input.point,
    );

    let ffi_res = match res {
        | Ok(v) => {
            FfiResult {
                ok : Some(v),
                err : None::<String>,
            }
        },
        | Err(e) => {
            FfiResult {
                ok : None,
                err : Some(e),
            }
        },
    };

    to_c_string(
        serde_json::to_string(&ffi_res)
            .unwrap(),
    )
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_vector_calculus_curl_json(
    input_json : *const c_char
) -> *mut c_char {

    let input : CurlInput = match from_json_string(input_json) {
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

    let vars_refs : Vec<&str> = input
        .vars
        .iter()
        .map(|s| s.as_str())
        .collect();

    let res =
        vector_calculus::curl_expr(
            &input.funcs,
            &vars_refs,
            &input.point,
        );

    let ffi_res = match res {
        | Ok(v) => {
            FfiResult {
                ok : Some(v),
                err : None::<String>,
            }
        },
        | Err(e) => {
            FfiResult {
                ok : None,
                err : Some(e),
            }
        },
    };

    to_c_string(
        serde_json::to_string(&ffi_res)
            .unwrap(),
    )
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_vector_calculus_laplacian_json(
    input_json : *const c_char
) -> *mut c_char {

    let input : LaplacianInput = match from_json_string(input_json) {
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

    let vars_refs : Vec<&str> = input
        .vars
        .iter()
        .map(|s| s.as_str())
        .collect();

    let res =
        vector_calculus::laplacian(
            &input.f,
            &vars_refs,
            &input.point,
        );

    let ffi_res = match res {
        | Ok(v) => {
            FfiResult {
                ok : Some(v),
                err : None::<String>,
            }
        },
        | Err(e) => {
            FfiResult {
                ok : None,
                err : Some(e),
            }
        },
    };

    to_c_string(
        serde_json::to_string(&ffi_res)
            .unwrap(),
    )
}

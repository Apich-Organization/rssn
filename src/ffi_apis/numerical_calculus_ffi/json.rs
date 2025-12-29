//! JSON-based FFI API for numerical calculus.

use std::os::raw::c_char;

use serde::Deserialize;

use crate::ffi_apis::common::from_json_string;
use crate::ffi_apis::common::to_c_string;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::calculus;
use crate::symbolic::core::Expr;

#[derive(Deserialize)]

struct GradientInput {
    expr: Expr,
    vars: Vec<String>,
    point: Vec<f64>,
}

#[derive(Deserialize)]

struct JacobianInput {
    funcs: Vec<Expr>,
    vars: Vec<String>,
    point: Vec<f64>,
}

#[derive(Deserialize)]

struct HessianInput {
    expr: Expr,
    vars: Vec<String>,
    point: Vec<f64>,
}

/// Computes the gradient of an expression at a given point using JSON for serialization.

#[no_mangle]

pub unsafe extern "C" fn rssn_numerical_gradient_json(
    input_json: *const c_char
) -> *mut c_char {

    let input : GradientInput = match from_json_string(input_json) {
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

    let vars_refs: Vec<&str> = input
        .vars
        .iter()
        .map(std::string::String::as_str)
        .collect();

    let res = calculus::gradient(
        &input.expr,
        &vars_refs,
        &input.point,
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

/// Computes the Jacobian matrix of a set of expressions at a given point using JSON for serialization.

#[no_mangle]

pub unsafe extern "C" fn rssn_numerical_jacobian_json(
    input_json: *const c_char
) -> *mut c_char {

    let input : JacobianInput = match from_json_string(input_json) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<Vec<Vec<f64>>, String> {
                        ok : None,
                        err : Some("Invalid JSON input".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let vars_refs: Vec<&str> = input
        .vars
        .iter()
        .map(std::string::String::as_str)
        .collect();

    let res = calculus::jacobian(
        &input.funcs,
        &vars_refs,
        &input.point,
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

/// Computes the Hessian matrix of an expression at a given point using JSON for serialization.

#[no_mangle]

pub unsafe extern "C" fn rssn_numerical_hessian_json(
    input_json: *const c_char
) -> *mut c_char {

    let input : HessianInput = match from_json_string(input_json) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<Vec<Vec<f64>>, String> {
                        ok : None,
                        err : Some("Invalid JSON input".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let vars_refs: Vec<&str> = input
        .vars
        .iter()
        .map(std::string::String::as_str)
        .collect();

    let res = calculus::hessian(
        &input.expr,
        &vars_refs,
        &input.point,
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

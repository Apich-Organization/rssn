//! Bincode-based FFI API for numerical calculus.

use serde::Deserialize;
use serde::Serialize;

use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::ffi_apis::common::BincodeBuffer;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::calculus;
use crate::symbolic::core::Expr;

#[derive(Deserialize)]

struct GradientInput {
    expr : Expr,
    vars : Vec<String>,
    point : Vec<f64>,
}

#[derive(Deserialize)]

struct JacobianInput {
    funcs : Vec<Expr>,
    vars : Vec<String>,
    point : Vec<f64>,
}

#[derive(Deserialize)]

struct HessianInput {
    expr : Expr,
    vars : Vec<String>,
    point : Vec<f64>,
}

#[no_mangle]

pub unsafe extern "C" fn rssn_numerical_gradient_bincode(
    buffer : BincodeBuffer
) -> BincodeBuffer {

    let input: GradientInput =
        match from_bincode_buffer(&buffer) {
            | Some(i) => i,
            | None => {
                return to_bincode_buffer(&FfiResult::<
                    Vec<f64>,
                    String,
                > {
                    ok: None,
                    err: Some(
                        "Invalid Bincode input".to_string(),
                    ),
                })
            },
        };

    let vars_refs : Vec<&str> = input
        .vars
        .iter()
        .map(|s| s.as_str())
        .collect();

    let res = calculus::gradient(
        &input.expr,
        &vars_refs,
        &input.point,
    );

    let ffi_res = match res {
        | Ok(v) => {
            FfiResult {
                ok : Some(v),
                err : None,
            }
        },
        | Err(e) => {
            FfiResult {
                ok : None,
                err : Some(e),
            }
        },
    };

    to_bincode_buffer(&ffi_res)
}

#[no_mangle]

pub unsafe extern "C" fn rssn_numerical_jacobian_bincode(
    buffer : BincodeBuffer
) -> BincodeBuffer {

    let input: JacobianInput =
        match from_bincode_buffer(&buffer) {
            | Some(i) => i,
            | None => {
                return to_bincode_buffer(&FfiResult::<
                    Vec<Vec<f64>>,
                    String,
                > {
                    ok: None,
                    err: Some(
                        "Invalid Bincode input".to_string(),
                    ),
                })
            },
        };

    let vars_refs : Vec<&str> = input
        .vars
        .iter()
        .map(|s| s.as_str())
        .collect();

    let res = calculus::jacobian(
        &input.funcs,
        &vars_refs,
        &input.point,
    );

    let ffi_res = match res {
        | Ok(v) => {
            FfiResult {
                ok : Some(v),
                err : None,
            }
        },
        | Err(e) => {
            FfiResult {
                ok : None,
                err : Some(e),
            }
        },
    };

    to_bincode_buffer(&ffi_res)
}

#[no_mangle]

pub unsafe extern "C" fn rssn_numerical_hessian_bincode(
    buffer : BincodeBuffer
) -> BincodeBuffer {

    let input: HessianInput =
        match from_bincode_buffer(&buffer) {
            | Some(i) => i,
            | None => {
                return to_bincode_buffer(&FfiResult::<
                    Vec<Vec<f64>>,
                    String,
                > {
                    ok: None,
                    err: Some(
                        "Invalid Bincode input".to_string(),
                    ),
                })
            },
        };

    let vars_refs : Vec<&str> = input
        .vars
        .iter()
        .map(|s| s.as_str())
        .collect();

    let res = calculus::hessian(
        &input.expr,
        &vars_refs,
        &input.point,
    );

    let ffi_res = match res {
        | Ok(v) => {
            FfiResult {
                ok : Some(v),
                err : None,
            }
        },
        | Err(e) => {
            FfiResult {
                ok : None,
                err : Some(e),
            }
        },
    };

    to_bincode_buffer(&ffi_res)
}

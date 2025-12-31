//! Bincode-based FFI API for numerical calculus.

use serde::Deserialize;

use crate::ffi_apis::common::BincodeBuffer;
use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
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

/// Computes the gradient of an expression at a given point using bincode for serialization.

#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_numerical_gradient_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input : GradientInput = match from_bincode_buffer(&buffer) {
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

    let vars_refs: Vec<&str> = input
        .vars
        .iter()
        .map(
            std::string::String::as_str,
        )
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

    to_bincode_buffer(&ffi_res)
}

/// Computes the Jacobian matrix of a set of expressions at a given point using bincode for serialization.

#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_numerical_jacobian_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input : JacobianInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(
                &FfiResult::<Vec<Vec<f64>>, String> {
                    ok : None,
                    err : Some("Invalid Bincode input".to_string()),
                },
            )
        },
    };

    let vars_refs: Vec<&str> = input
        .vars
        .iter()
        .map(
            std::string::String::as_str,
        )
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

    to_bincode_buffer(&ffi_res)
}

/// Computes the Hessian matrix of an expression at a given point using bincode for serialization.

#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_numerical_hessian_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input : HessianInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(
                &FfiResult::<Vec<Vec<f64>>, String> {
                    ok : None,
                    err : Some("Invalid Bincode input".to_string()),
                },
            )
        },
    };

    let vars_refs: Vec<&str> = input
        .vars
        .iter()
        .map(
            std::string::String::as_str,
        )
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

    to_bincode_buffer(&ffi_res)
}

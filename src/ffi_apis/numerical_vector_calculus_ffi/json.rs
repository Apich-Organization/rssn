//! JSON-based FFI API for numerical vector calculus.

use std::os::raw::c_char;

use serde::Deserialize;

use crate::ffi_apis::common::from_json_string;
use crate::ffi_apis::common::to_c_string;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::vector_calculus;
use crate::symbolic::core::Expr;

#[derive(Deserialize)]

struct DivergenceInput {
    funcs: Vec<Expr>,
    vars: Vec<String>,
    point: Vec<f64>,
}

#[derive(Deserialize)]

struct CurlInput {
    funcs: Vec<Expr>,
    vars: Vec<String>,
    point: Vec<f64>,
}

#[derive(Deserialize)]

struct LaplacianInput {
    f: Expr,
    vars: Vec<String>,
    point: Vec<f64>,
}

/// Computes the divergence of a vector field at a point via JSON serialization.
///
/// The divergence measures the net outward flux of a vector field:
/// div(F) = ∂F₁/∂x₁ + ∂F₂/∂x₂ + ... + ∂Fₙ/∂xₙ.
///
/// # Arguments
///
/// * `input_json` - A JSON string pointer containing:
///   - `funcs`: Vector field components as symbolic expressions
///   - `vars`: Variable names corresponding to coordinates
///   - `point`: Point at which to evaluate divergence
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<f64, String>` with
/// the divergence value (scalar).
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

/// # Panics
///
/// This function may panic if the FFI input is malformed, null where not expected,
/// or if internal state synchronization fails (e.g., poisoned locks).

pub unsafe extern "C" fn rssn_num_vector_calculus_divergence_json(
    input_json: *const c_char
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

    let vars_refs: Vec<&str> = input
        .vars
        .iter()
        .map(
            std::string::String::as_str,
        )
        .collect();

    let res = vector_calculus::divergence_expr(
        &input.funcs,
        &vars_refs,
        &input.point,
    );

    let ffi_res = match res {
        | Ok(v) => {
            FfiResult {
                ok: Some(v),
                err: None::<String>,
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

/// Computes the curl of a vector field at a point via JSON serialization.
///
/// The curl measures the rotational tendency of a vector field. In 3D:
/// curl(F) = (∂F₃/∂x₂ - ∂F₂/∂x₃, ∂F₁/∂x₃ - ∂F₃/∂x₁, ∂F₂/∂x₁ - ∂F₁/∂x₂).
///
/// # Arguments
///
/// * `input_json` - A JSON string pointer containing:
///   - `funcs`: Vector field components as symbolic expressions
///   - `vars`: Variable names corresponding to coordinates
///   - `point`: Point at which to evaluate curl
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<Vec<f64>, String>` with
/// the curl vector.
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

/// # Panics
///
/// This function may panic if the FFI input is malformed, null where not expected,
/// or if internal state synchronization fails (e.g., poisoned locks).

pub unsafe extern "C" fn rssn_num_vector_calculus_curl_json(
    input_json: *const c_char
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

    let vars_refs: Vec<&str> = input
        .vars
        .iter()
        .map(
            std::string::String::as_str,
        )
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
                ok: Some(v),
                err: None::<String>,
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

/// Computes the Laplacian of a scalar field at a point via JSON serialization.
///
/// The Laplacian is the divergence of the gradient:
/// ∇²f = ∂²f/∂x₁² + ∂²f/∂x₂² + ... + ∂²f/∂xₙ².
///
/// # Arguments
///
/// * `input_json` - A JSON string pointer containing:
///   - `f`: Scalar field as a symbolic expression
///   - `vars`: Variable names corresponding to coordinates
///   - `point`: Point at which to evaluate Laplacian
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<f64, String>` with
/// the Laplacian value (scalar).
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

/// # Panics
///
/// This function may panic if the FFI input is malformed, null where not expected,
/// or if internal state synchronization fails (e.g., poisoned locks).

pub unsafe extern "C" fn rssn_num_vector_calculus_laplacian_json(
    input_json: *const c_char
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

    let vars_refs: Vec<&str> = input
        .vars
        .iter()
        .map(
            std::string::String::as_str,
        )
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
                ok: Some(v),
                err: None::<String>,
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

//! Bincode-based FFI API for numerical vector calculus.

use serde::Deserialize;
use serde::Serialize;

use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::ffi_apis::common::BincodeBuffer;
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

/// Computes the divergence of a vector field at a point using bincode serialization.
///
/// The divergence measures the net outward flux of a vector field:
/// div(F) = ∂F₁/∂x₁ + ∂F₂/∂x₂ + ... + ∂Fₙ/∂xₙ.
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing `DivergenceInput` with:
///   - `funcs`: Vector field components as symbolic expressions
///   - `vars`: Variable names corresponding to coordinates
///   - `point`: Point at which to evaluate divergence
///
/// # Returns
///
/// A bincode-encoded buffer containing `FfiResult<f64, String>` with either:
/// - `ok`: The divergence value (scalar)
/// - `err`: Error message if computation failed
///
/// # Safety
///
/// This function is unsafe because it receives a raw bincode buffer that must be
/// valid and properly encoded.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_vector_calculus_divergence_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input : DivergenceInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(
                &FfiResult::<f64, String> {
                    ok : None,
                    err : Some("Invalid Bincode input".to_string()),
                },
            )
        },
    };

    let vars_refs: Vec<&str> = input
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

    to_bincode_buffer(&ffi_res)
}

/// Computes the curl of a vector field at a point using bincode serialization.
///
/// The curl measures the rotational tendency of a vector field. In 3D:
/// curl(F) = (∂F₃/∂x₂ - ∂F₂/∂x₃, ∂F₁/∂x₃ - ∂F₃/∂x₁, ∂F₂/∂x₁ - ∂F₁/∂x₂).
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing `CurlInput` with:
///   - `funcs`: Vector field components as symbolic expressions
///   - `vars`: Variable names corresponding to coordinates
///   - `point`: Point at which to evaluate curl
///
/// # Returns
///
/// A bincode-encoded buffer containing `FfiResult<Vec<f64>, String>` with either:
/// - `ok`: The curl vector
/// - `err`: Error message if computation failed
///
/// # Safety
///
/// This function is unsafe because it receives a raw bincode buffer that must be
/// valid and properly encoded.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_vector_calculus_curl_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input : CurlInput = match from_bincode_buffer(&buffer) {
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

    to_bincode_buffer(&ffi_res)
}

/// Computes the Laplacian of a scalar field at a point using bincode serialization.
///
/// The Laplacian is the divergence of the gradient:
/// ∇²f = ∂²f/∂x₁² + ∂²f/∂x₂² + ... + ∂²f/∂xₙ².
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing `LaplacianInput` with:
///   - `f`: Scalar field as a symbolic expression
///   - `vars`: Variable names corresponding to coordinates
///   - `point`: Point at which to evaluate Laplacian
///
/// # Returns
///
/// A bincode-encoded buffer containing `FfiResult<f64, String>` with either:
/// - `ok`: The Laplacian value (scalar)
/// - `err`: Error message if computation failed
///
/// # Safety
///
/// This function is unsafe because it receives a raw bincode buffer that must be
/// valid and properly encoded.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_vector_calculus_laplacian_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input : LaplacianInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(
                &FfiResult::<f64, String> {
                    ok : None,
                    err : Some("Invalid Bincode input".to_string()),
                },
            )
        },
    };

    let vars_refs: Vec<&str> = input
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

    to_bincode_buffer(&ffi_res)
}

//! Bincode-based FFI API for numerical vector calculus.

use crate::ffi_apis::common::{from_bincode_buffer, to_bincode_buffer, BincodeBuffer};
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::vector_calculus;
use crate::symbolic::core::Expr;
use serde::{Deserialize, Serialize};

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

#[no_mangle]
pub unsafe extern "C" fn rssn_num_vector_calculus_divergence_bincode(
    buffer: BincodeBuffer,
) -> BincodeBuffer {
    let input: DivergenceInput = match from_bincode_buffer(&buffer) {
        Some(i) => i,
        None => {
            return to_bincode_buffer(&FfiResult::<f64, String> {
                ok: None,
                err: Some("Invalid Bincode input".to_string()),
            })
        }
    };
    let vars_refs: Vec<&str> = input.vars.iter().map(|s| s.as_str()).collect();
    let res = vector_calculus::divergence_expr(&input.funcs, &vars_refs, &input.point);
    let ffi_res = match res {
        Ok(v) => FfiResult {
            ok: Some(v),
            err: None::<String>,
        },
        Err(e) => FfiResult {
            ok: None,
            err: Some(e),
        },
    };
    to_bincode_buffer(&ffi_res)
}

#[no_mangle]
pub unsafe extern "C" fn rssn_num_vector_calculus_curl_bincode(
    buffer: BincodeBuffer,
) -> BincodeBuffer {
    let input: CurlInput = match from_bincode_buffer(&buffer) {
        Some(i) => i,
        None => {
            return to_bincode_buffer(&FfiResult::<Vec<f64>, String> {
                ok: None,
                err: Some("Invalid Bincode input".to_string()),
            })
        }
    };
    let vars_refs: Vec<&str> = input.vars.iter().map(|s| s.as_str()).collect();
    let res = vector_calculus::curl_expr(&input.funcs, &vars_refs, &input.point);
    let ffi_res = match res {
        Ok(v) => FfiResult {
            ok: Some(v),
            err: None::<String>,
        },
        Err(e) => FfiResult {
            ok: None,
            err: Some(e),
        },
    };
    to_bincode_buffer(&ffi_res)
}

#[no_mangle]
pub unsafe extern "C" fn rssn_num_vector_calculus_laplacian_bincode(
    buffer: BincodeBuffer,
) -> BincodeBuffer {
    let input: LaplacianInput = match from_bincode_buffer(&buffer) {
        Some(i) => i,
        None => {
            return to_bincode_buffer(&FfiResult::<f64, String> {
                ok: None,
                err: Some("Invalid Bincode input".to_string()),
            })
        }
    };
    let vars_refs: Vec<&str> = input.vars.iter().map(|s| s.as_str()).collect();
    let res = vector_calculus::laplacian(&input.f, &vars_refs, &input.point);
    let ffi_res = match res {
        Ok(v) => FfiResult {
            ok: Some(v),
            err: None::<String>,
        },
        Err(e) => FfiResult {
            ok: None,
            err: Some(e),
        },
    };
    to_bincode_buffer(&ffi_res)
}

//! JSON-based FFI API for numerical sparse matrix operations.

use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::sparse::{
    self,
    SparseMatrixData,
};
use serde::{
    Deserialize,
    Serialize,
};
use std::ffi::{
    CStr,
    CString,
};
use std::os::raw::c_char;

#[derive(Deserialize)]

struct SpMvRequest {
    matrix: SparseMatrixData,
    vector: Vec<f64>,
}

/// Sparse matrix-vector multiplication from JSON.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_sparse_spmv_json(json_ptr: *const c_char) -> *mut c_char {

    if json_ptr.is_null() {

        return std::ptr::null_mut();
    }

    let json_str = match unsafe {

        CStr::from_ptr(json_ptr).to_str()
    } {
        | Ok(s) => s,
        | Err(_) => return std::ptr::null_mut(),
    };

    let req: SpMvRequest = match serde_json::from_str(json_str) {
        | Ok(r) => r,
        | Err(e) => {

            let res: FfiResult<Vec<f64>, String> = FfiResult {
                ok: None,
                err: Some(e.to_string()),
            };

            return CString::new(serde_json::to_string(&res).unwrap())
                .unwrap()
                .into_raw();
        },
    };

    let mat = req
        .matrix
        .to_csmat();

    match sparse::sp_mat_vec_mul(&mat, &req.vector) {
        | Ok(res) => {

            let ffi_res: FfiResult<Vec<f64>, String> = FfiResult {
                ok: Some(res),
                err: None,
            };

            CString::new(serde_json::to_string(&ffi_res).unwrap())
                .unwrap()
                .into_raw()
        },
        | Err(e) => {

            let ffi_res: FfiResult<Vec<f64>, String> = FfiResult {
                ok: None,
                err: Some(e),
            };

            CString::new(serde_json::to_string(&ffi_res).unwrap())
                .unwrap()
                .into_raw()
        },
    }
}

/// Solves Ax=b using Conjugate Gradient from JSON.
#[derive(Deserialize)]

struct CgRequest {
    a: SparseMatrixData,
    b: Vec<f64>,
    x0: Option<Vec<f64>>,
    max_iter: usize,
    tolerance: f64,
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_sparse_solve_cg_json(json_ptr: *const c_char) -> *mut c_char {

    if json_ptr.is_null() {

        return std::ptr::null_mut();
    }

    let json_str = match unsafe {

        CStr::from_ptr(json_ptr).to_str()
    } {
        | Ok(s) => s,
        | Err(_) => return std::ptr::null_mut(),
    };

    let req: CgRequest = match serde_json::from_str(json_str) {
        | Ok(r) => r,
        | Err(e) => {

            let res: FfiResult<Vec<f64>, String> = FfiResult {
                ok: None,
                err: Some(e.to_string()),
            };

            return CString::new(serde_json::to_string(&res).unwrap())
                .unwrap()
                .into_raw();
        },
    };

    let a = req.a.to_csmat();

    let b = ndarray::Array1::from_vec(req.b);

    let x0 = req
        .x0
        .map(ndarray::Array1::from_vec);

    match sparse::solve_conjugate_gradient(
        &a,
        &b,
        x0.as_ref(),
        req.max_iter,
        req.tolerance,
    ) {
        | Ok(res) => {

            let ffi_res: FfiResult<Vec<f64>, String> = FfiResult {
                ok: Some(res.to_vec()),
                err: None,
            };

            CString::new(serde_json::to_string(&ffi_res).unwrap())
                .unwrap()
                .into_raw()
        },
        | Err(e) => {

            let ffi_res: FfiResult<Vec<f64>, String> = FfiResult {
                ok: None,
                err: Some(e),
            };

            CString::new(serde_json::to_string(&ffi_res).unwrap())
                .unwrap()
                .into_raw()
        },
    }
}

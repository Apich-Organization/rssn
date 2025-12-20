//! JSON-based FFI API for numerical matrix operations.

use crate::numerical::matrix::Matrix;
use crate::ffi_apis::ffi_api::FfiResult;
use serde::{Deserialize, Serialize};
use std::os::raw::c_char;
use std::ffi::{CStr, CString};

#[derive(Deserialize)]
struct MatrixOpRequest {
    m1: Matrix<f64>,
    m2: Option<Matrix<f64>>,
}

/// Evaluates a matrix addition from JSON.
#[no_mangle]
pub unsafe extern "C" fn rssn_num_matrix_add_json(json_ptr: *const c_char) -> *mut c_char {
    if json_ptr.is_null() { return std::ptr::null_mut(); }
    let json_str = match unsafe { CStr::from_ptr(json_ptr).to_str() } {
        Ok(s) => s,
        Err(_) => return std::ptr::null_mut(),
    };

    let req: MatrixOpRequest = match serde_json::from_str(json_str) {
        Ok(r) => r,
        Err(e) => {
            let res: FfiResult<Matrix<f64>, String> = FfiResult { ok: None, err: Some(e.to_string()) };
            return CString::new(serde_json::to_string(&res).unwrap()).unwrap().into_raw();
        }
    };

    let m2 = match req.m2 {
        Some(m) => m,
        None => {
            let res: FfiResult<Matrix<f64>, String> = FfiResult { ok: None, err: Some("Second matrix m2 is required".to_string()) };
            return CString::new(serde_json::to_string(&res).unwrap()).unwrap().into_raw();
        }
    };

    if req.m1.rows() != m2.rows() || req.m1.cols() != m2.cols() {
        let res: FfiResult<Matrix<f64>, String> = FfiResult { ok: None, err: Some("Dimension mismatch".to_string()) };
        return CString::new(serde_json::to_string(&res).unwrap()).unwrap().into_raw();
    }

    let result = req.m1 + m2;
    let ffi_res: FfiResult<Matrix<f64>, String> = FfiResult { ok: Some(result), err: None };
    CString::new(serde_json::to_string(&ffi_res).unwrap()).unwrap().into_raw()
}

/// Evaluates a matrix multiplication from JSON.
#[no_mangle]
pub unsafe extern "C" fn rssn_num_matrix_mul_json(json_ptr: *const c_char) -> *mut c_char {
    if json_ptr.is_null() { return std::ptr::null_mut(); }
    let json_str = match unsafe { CStr::from_ptr(json_ptr).to_str() } {
        Ok(s) => s,
        Err(_) => return std::ptr::null_mut(),
    };

    let req: MatrixOpRequest = match serde_json::from_str(json_str) {
        Ok(r) => r,
        Err(e) => {
            let res: FfiResult<Matrix<f64>, String> = FfiResult { ok: None, err: Some(e.to_string()) };
            return CString::new(serde_json::to_string(&res).unwrap()).unwrap().into_raw();
        }
    };

    let m2 = match req.m2 {
        Some(m) => m,
        None => {
            let res: FfiResult<Matrix<f64>, String> = FfiResult { ok: None, err: Some("Second matrix m2 is required".to_string()) };
            return CString::new(serde_json::to_string(&res).unwrap()).unwrap().into_raw();
        }
    };

    if req.m1.cols() != m2.rows() {
        let res: FfiResult<Matrix<f64>, String> = FfiResult { ok: None, err: Some("Dimension mismatch".to_string()) };
        return CString::new(serde_json::to_string(&res).unwrap()).unwrap().into_raw();
    }

    let result = req.m1 * m2;
    let ffi_res: FfiResult<Matrix<f64>, String> = FfiResult { ok: Some(result), err: None };
    CString::new(serde_json::to_string(&ffi_res).unwrap()).unwrap().into_raw()
}

/// Computes determinant from JSON.
#[no_mangle]
pub unsafe extern "C" fn rssn_num_matrix_det_json(json_ptr: *const c_char) -> *mut c_char {
    if json_ptr.is_null() { return std::ptr::null_mut(); }
    let json_str = match unsafe { CStr::from_ptr(json_ptr).to_str() } {
        Ok(s) => s,
        Err(_) => return std::ptr::null_mut(),
    };

    let matrix: Matrix<f64> = match serde_json::from_str(json_str) {
        Ok(m) => m,
        Err(e) => {
            let res: FfiResult<f64, String> = FfiResult { ok: None, err: Some(e.to_string()) };
            return CString::new(serde_json::to_string(&res).unwrap()).unwrap().into_raw();
        }
    };

    match matrix.determinant() {
        Ok(d) => {
            let ffi_res: FfiResult<f64, String> = FfiResult { ok: Some(d), err: None };
            CString::new(serde_json::to_string(&ffi_res).unwrap()).unwrap().into_raw()
        }
        Err(e) => {
            let ffi_res: FfiResult<f64, String> = FfiResult { ok: None, err: Some(e) };
            CString::new(serde_json::to_string(&ffi_res).unwrap()).unwrap().into_raw()
        }
    }
}

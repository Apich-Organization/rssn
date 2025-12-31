//! JSON-based FFI API for numerical matrix operations.

use std::ffi::CStr;
use std::ffi::CString;
use std::os::raw::c_char;

use serde::Deserialize;

use crate::ffi_apis::ffi_api::FfiResult;
use crate::nightly::matrix::Backend;
use crate::nightly::matrix::FaerDecompositionResult;
use crate::nightly::matrix::Matrix;

#[derive(Deserialize)]

struct MatrixOpRequest {
    m1: Matrix<f64>,
    m2: Option<Matrix<f64>>,
}

#[derive(Deserialize)]

struct MatrixBackendRequest {
    matrix: Matrix<f64>,
    backend_id: i32, /* 0: Native, 1: Faer */
}

#[derive(Deserialize)]

struct MatrixDecompositionRequest {
    matrix: Matrix<f64>,
    kind: crate::nightly::matrix::FaerDecompositionType,
}

/// Evaluates a matrix addition from JSON.
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

pub unsafe extern "C" fn rssn_num_matrix_add_json_nightly(
    json_ptr: *const c_char
) -> *mut c_char {

    if json_ptr.is_null() {

        return std::ptr::null_mut();
    }

    let json_str = match unsafe {

        CStr::from_ptr(json_ptr)
            .to_str()
    } {
        | Ok(s) => s,
        | Err(_) => {
            return std::ptr::null_mut();
        },
    };

    let req: MatrixOpRequest =
        match serde_json::from_str(
            json_str,
        ) {
            | Ok(r) => r,
            | Err(e) => {

                let res: FfiResult<
                    Matrix<f64>,
                    String,
                > = FfiResult {
                    ok: None,
                    err: Some(
                        e.to_string(),
                    ),
                };

                return CString::new(serde_json::to_string(&res).unwrap())
                .unwrap()
                .into_raw();
            },
        };

    let m2 = match req.m2 {
        | Some(m) => m,
        | None => {

            let res: FfiResult<
                Matrix<f64>,
                String,
            > = FfiResult {
                ok: None,
                err: Some(
                    "Second matrix m2 \
                     is required"
                        .to_string(),
                ),
            };

            return CString::new(
                serde_json::to_string(
                    &res,
                )
                .unwrap(),
            )
            .unwrap()
            .into_raw();
        },
    };

    if req.m1.rows() != m2.rows()
        || req.m1.cols() != m2.cols()
    {

        let res: FfiResult<
            Matrix<f64>,
            String,
        > = FfiResult {
            ok: None,
            err: Some(
                "Dimension mismatch"
                    .to_string(),
            ),
        };

        return CString::new(
            serde_json::to_string(&res)
                .unwrap(),
        )
        .unwrap()
        .into_raw();
    }

    let result = req.m1 + m2;

    let ffi_res: FfiResult<
        Matrix<f64>,
        String,
    > = FfiResult {
        ok: Some(result),
        err: None,
    };

    CString::new(
        serde_json::to_string(&ffi_res)
            .unwrap(),
    )
    .unwrap()
    .into_raw()
}

/// Evaluates a matrix multiplication from JSON.
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

pub unsafe extern "C" fn rssn_num_matrix_mul_json_nightly(
    json_ptr: *const c_char
) -> *mut c_char {

    if json_ptr.is_null() {

        return std::ptr::null_mut();
    }

    let json_str = match unsafe {

        CStr::from_ptr(json_ptr)
            .to_str()
    } {
        | Ok(s) => s,
        | Err(_) => {
            return std::ptr::null_mut();
        },
    };

    let req: MatrixOpRequest =
        match serde_json::from_str(
            json_str,
        ) {
            | Ok(r) => r,
            | Err(e) => {

                let res: FfiResult<
                    Matrix<f64>,
                    String,
                > = FfiResult {
                    ok: None,
                    err: Some(
                        e.to_string(),
                    ),
                };

                return CString::new(serde_json::to_string(&res).unwrap())
                .unwrap()
                .into_raw();
            },
        };

    let m2 = match req.m2 {
        | Some(m) => m,
        | None => {

            let res: FfiResult<
                Matrix<f64>,
                String,
            > = FfiResult {
                ok: None,
                err: Some(
                    "Second matrix m2 \
                     is required"
                        .to_string(),
                ),
            };

            return CString::new(
                serde_json::to_string(
                    &res,
                )
                .unwrap(),
            )
            .unwrap()
            .into_raw();
        },
    };

    if req.m1.cols() != m2.rows() {

        let res: FfiResult<
            Matrix<f64>,
            String,
        > = FfiResult {
            ok: None,
            err: Some(
                "Dimension mismatch"
                    .to_string(),
            ),
        };

        return CString::new(
            serde_json::to_string(&res)
                .unwrap(),
        )
        .unwrap()
        .into_raw();
    }

    let result = req.m1 * m2;

    let ffi_res: FfiResult<
        Matrix<f64>,
        String,
    > = FfiResult {
        ok: Some(result),
        err: None,
    };

    CString::new(
        serde_json::to_string(&ffi_res)
            .unwrap(),
    )
    .unwrap()
    .into_raw()
}

/// Computes determinant from JSON.
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

pub unsafe extern "C" fn rssn_num_matrix_det_json_nightly(
    json_ptr: *const c_char
) -> *mut c_char {

    if json_ptr.is_null() {

        return std::ptr::null_mut();
    }

    let json_str = match unsafe {

        CStr::from_ptr(json_ptr)
            .to_str()
    } {
        | Ok(s) => s,
        | Err(_) => {
            return std::ptr::null_mut();
        },
    };

    let matrix: Matrix<f64> =
        match serde_json::from_str(
            json_str,
        ) {
            | Ok(m) => m,
            | Err(e) => {

                let res: FfiResult<
                    f64,
                    String,
                > = FfiResult {
                    ok: None,
                    err: Some(
                        e.to_string(),
                    ),
                };

                return CString::new(serde_json::to_string(&res).unwrap())
                .unwrap()
                .into_raw();
            },
        };

    match matrix.determinant() {
        | Ok(d) => {

            let ffi_res: FfiResult<
                f64,
                String,
            > = FfiResult {
                ok: Some(d),
                err: None,
            };

            CString::new(
                serde_json::to_string(
                    &ffi_res,
                )
                .unwrap(),
            )
            .unwrap()
            .into_raw()
        },
        | Err(e) => {

            let ffi_res: FfiResult<
                f64,
                String,
            > = FfiResult {
                ok: None,
                err: Some(e),
            };

            CString::new(
                serde_json::to_string(
                    &ffi_res,
                )
                .unwrap(),
            )
            .unwrap()
            .into_raw()
        },
    }
}

/// Sets backend for a matrix (returns new matrix with backend set) from JSON.
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

pub unsafe extern "C" fn rssn_num_matrix_set_backend_json_nightly(
    json_ptr: *const c_char
) -> *mut c_char {

    if json_ptr.is_null() {

        return std::ptr::null_mut();
    }

    let json_str = match unsafe {

        CStr::from_ptr(json_ptr)
            .to_str()
    } {
        | Ok(s) => s,
        | Err(_) => {
            return std::ptr::null_mut();
        },
    };

    let req: MatrixBackendRequest =
        match serde_json::from_str(
            json_str,
        ) {
            | Ok(r) => r,
            | Err(e) => {

                let res: FfiResult<
                    Matrix<f64>,
                    String,
                > = FfiResult {
                    ok: None,
                    err: Some(
                        e.to_string(),
                    ),
                };

                return CString::new(serde_json::to_string(&res).unwrap()).unwrap().into_raw();
            },
        };

    let backend = match req.backend_id {
        | 1 => Backend::Faer,
        | _ => Backend::Native,
    };

    let m = req
        .matrix
        .with_backend(backend);

    let ffi_res: FfiResult<
        Matrix<f64>,
        String,
    > = FfiResult {
        ok: Some(m),
        err: None,
    };

    CString::new(
        serde_json::to_string(&ffi_res)
            .unwrap(),
    )
    .unwrap()
    .into_raw()
}

/// Decomposes a matrix from JSON.
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

pub unsafe extern "C" fn rssn_num_matrix_decompose_json_nightly(
    json_ptr: *const c_char
) -> *mut c_char {

    if json_ptr.is_null() {

        return std::ptr::null_mut();
    }

    let json_str = match unsafe {

        CStr::from_ptr(json_ptr)
            .to_str()
    } {
        | Ok(s) => s,
        | Err(_) => {
            return std::ptr::null_mut();
        },
    };

    let req: MatrixDecompositionRequest = match serde_json::from_str(json_str) {
        Ok(r) => r,
        Err(e) => {
             let res: FfiResult<FaerDecompositionResult<f64>, String> = FfiResult { ok: None, err: Some(e.to_string()) };
             return CString::new(serde_json::to_string(&res).unwrap()).unwrap().into_raw();
        }
    };

    match req
        .matrix
        .decompose(req.kind)
    {
        | Some(result) => {

            let ffi_res: FfiResult<
                FaerDecompositionResult<
                    f64,
                >,
                String,
            > = FfiResult {
                ok: Some(result),
                err: None,
            };

            CString::new(
                serde_json::to_string(
                    &ffi_res,
                )
                .unwrap(),
            )
            .unwrap()
            .into_raw()
        },
        | None => {

            let ffi_res: FfiResult<
                FaerDecompositionResult<
                    f64,
                >,
                String,
            > = FfiResult {
                ok: None,
                err: Some(
                    "Decomposition \
                     failed or backend \
                     not supported"
                        .to_string(),
                ),
            };

            CString::new(
                serde_json::to_string(
                    &ffi_res,
                )
                .unwrap(),
            )
            .unwrap()
            .into_raw()
        },
    }
}

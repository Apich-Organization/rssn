//! JSON-based FFI API for numerical tensor operations.

use std::ffi::CStr;
use std::ffi::CString;
use std::os::raw::c_char;

use serde::Deserialize;

use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::tensor::TensorData;
use crate::numerical::tensor::{
    self,
};

#[derive(Deserialize)]

struct TensordotRequest {
    a: TensorData,
    b: TensorData,
    axes_a: Vec<usize>,
    axes_b: Vec<usize>,
}

/// Tensor contraction from JSON.
#[no_mangle]

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

pub unsafe extern "C" fn rssn_num_tensor_tensordot_json(
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
            return std::ptr::null_mut()
        },
    };

    let req: TensordotRequest =
        match serde_json::from_str(
            json_str,
        ) {
            | Ok(r) => r,
            | Err(e) => {

                let res: FfiResult<
                    TensorData,
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

    let a = match req.a.to_arrayd() {
        | Ok(arr) => arr,
        | Err(e) => {

            let res: FfiResult<
                TensorData,
                String,
            > = FfiResult {
                ok: None,
                err: Some(e),
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

    let b = match req.b.to_arrayd() {
        | Ok(arr) => arr,
        | Err(e) => {

            let res: FfiResult<
                TensorData,
                String,
            > = FfiResult {
                ok: None,
                err: Some(e),
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

    match tensor::tensordot(
        &a,
        &b,
        &req.axes_a,
        &req.axes_b,
    ) {
        | Ok(res) => {

            let ffi_res: FfiResult<
                TensorData,
                String,
            > = FfiResult {
                ok: Some(
                    TensorData::from(
                        &res,
                    ),
                ),
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
                TensorData,
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

#[derive(Deserialize)]

struct OuterProductRequest {
    a: TensorData,
    b: TensorData,
}

/// Outer product from JSON.
#[no_mangle]

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

pub unsafe extern "C" fn rssn_num_tensor_outer_product_json(
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
            return std::ptr::null_mut()
        },
    };

    let req: OuterProductRequest =
        match serde_json::from_str(
            json_str,
        ) {
            | Ok(r) => r,
            | Err(e) => {

                let res: FfiResult<
                    TensorData,
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

    let a = match req.a.to_arrayd() {
        | Ok(arr) => arr,
        | Err(e) => {

            let res: FfiResult<
                TensorData,
                String,
            > = FfiResult {
                ok: None,
                err: Some(e),
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

    let b = match req.b.to_arrayd() {
        | Ok(arr) => arr,
        | Err(e) => {

            let res: FfiResult<
                TensorData,
                String,
            > = FfiResult {
                ok: None,
                err: Some(e),
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

    match tensor::outer_product(&a, &b)
    {
        | Ok(res) => {

            let ffi_res: FfiResult<
                TensorData,
                String,
            > = FfiResult {
                ok: Some(
                    TensorData::from(
                        &res,
                    ),
                ),
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
                TensorData,
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

//! JSON-based FFI API for numerical coordinate transformations.

use std::ffi::CStr;
use std::ffi::CString;
use std::os::raw::c_char;

use serde::Deserialize;

use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::coordinates as nc;
use crate::symbolic::coordinates::CoordinateSystem;

#[derive(Deserialize)]

struct CoordinateTransformRequest {
    point: Vec<f64>,
    from: CoordinateSystem,
    to: CoordinateSystem,
}

/// Transforms a point using JSON.
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

pub unsafe extern "C" fn rssn_num_coord_transform_json(
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

    let req : CoordinateTransformRequest = match serde_json::from_str(json_str) {
        | Ok(r) => r,
        | Err(e) => {

            let res : FfiResult<Vec<f64>, String> = FfiResult {
                ok : None,
                err : Some(e.to_string()),
            };

            return CString::new(serde_json::to_string(&res).unwrap())
                .unwrap()
                .into_raw();
        },
    };

    match nc::transform_point(
        &req.point,
        req.from,
        req.to,
    ) {
        | Ok(res) => {

            let ffi_res: FfiResult<
                Vec<f64>,
                String,
            > = FfiResult {
                ok: Some(res),
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
                Vec<f64>,
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

/// Transforms a point (pure numerical) using JSON.
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

pub unsafe extern "C" fn rssn_num_coord_transform_pure_json(
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

    let req : CoordinateTransformRequest = match serde_json::from_str(json_str) {
        | Ok(r) => r,
        | Err(e) => {

            let res : FfiResult<Vec<f64>, String> = FfiResult {
                ok : None,
                err : Some(e.to_string()),
            };

            return CString::new(serde_json::to_string(&res).unwrap())
                .unwrap()
                .into_raw();
        },
    };

    match nc::transform_point_pure(
        &req.point,
        req.from,
        req.to,
    ) {
        | Ok(res) => {

            let ffi_res: FfiResult<
                Vec<f64>,
                String,
            > = FfiResult {
                ok: Some(res),
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
                Vec<f64>,
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

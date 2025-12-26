//! JSON-based FFI API for numerical geometric algebra operations.

use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::geometric_algebra::Multivector3D;
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

struct GaInput {
    mv: Multivector3D,
}

#[derive(Deserialize)]

struct TwoGaInput {
    mv1: Multivector3D,
    mv2: Multivector3D,
}

/// JSON FFI for ga_add.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_ga_add_json(json_ptr: *const c_char) -> *mut c_char {

    let json_str = match CStr::from_ptr(json_ptr).to_str() {
        Ok(s) => s,
        Err(_) => return std::ptr::null_mut(),
    };

    let input: TwoGaInput = match serde_json::from_str(json_str) {
        Ok(v) => v,
        Err(e) => {
            return CString::new(format!(
                "{{\"err\": \"{}\"}}",
                e
            ))
            .unwrap()
            .into_raw()
        }
    };

    let res = FfiResult {
        ok: Some(input.mv1 + input.mv2),
        err: None::<String>,
    };

    CString::new(serde_json::to_string(&res).unwrap())
        .unwrap()
        .into_raw()
}

/// JSON FFI for ga_sub.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_ga_sub_json(json_ptr: *const c_char) -> *mut c_char {

    let json_str = match CStr::from_ptr(json_ptr).to_str() {
        Ok(s) => s,
        Err(_) => return std::ptr::null_mut(),
    };

    let input: TwoGaInput = match serde_json::from_str(json_str) {
        Ok(v) => v,
        Err(e) => {
            return CString::new(format!(
                "{{\"err\": \"{}\"}}",
                e
            ))
            .unwrap()
            .into_raw()
        }
    };

    let res = FfiResult {
        ok: Some(input.mv1 - input.mv2),
        err: None::<String>,
    };

    CString::new(serde_json::to_string(&res).unwrap())
        .unwrap()
        .into_raw()
}

/// JSON FFI for ga_mul.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_ga_mul_json(json_ptr: *const c_char) -> *mut c_char {

    let json_str = match CStr::from_ptr(json_ptr).to_str() {
        Ok(s) => s,
        Err(_) => return std::ptr::null_mut(),
    };

    let input: TwoGaInput = match serde_json::from_str(json_str) {
        Ok(v) => v,
        Err(e) => {
            return CString::new(format!(
                "{{\"err\": \"{}\"}}",
                e
            ))
            .unwrap()
            .into_raw()
        }
    };

    let res = FfiResult {
        ok: Some(input.mv1 * input.mv2),
        err: None::<String>,
    };

    CString::new(serde_json::to_string(&res).unwrap())
        .unwrap()
        .into_raw()
}

/// JSON FFI for ga_wedge.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_ga_wedge_json(json_ptr: *const c_char) -> *mut c_char {

    let json_str = match CStr::from_ptr(json_ptr).to_str() {
        Ok(s) => s,
        Err(_) => return std::ptr::null_mut(),
    };

    let input: TwoGaInput = match serde_json::from_str(json_str) {
        Ok(v) => v,
        Err(e) => {
            return CString::new(format!(
                "{{\"err\": \"{}\"}}",
                e
            ))
            .unwrap()
            .into_raw()
        }
    };

    let res = FfiResult {
        ok: Some(
            input
                .mv1
                .wedge(input.mv2),
        ),
        err: None::<String>,
    };

    CString::new(serde_json::to_string(&res).unwrap())
        .unwrap()
        .into_raw()
}

/// JSON FFI for ga_dot.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_ga_dot_json(json_ptr: *const c_char) -> *mut c_char {

    let json_str = match CStr::from_ptr(json_ptr).to_str() {
        Ok(s) => s,
        Err(_) => return std::ptr::null_mut(),
    };

    let input: TwoGaInput = match serde_json::from_str(json_str) {
        Ok(v) => v,
        Err(e) => {
            return CString::new(format!(
                "{{\"err\": \"{}\"}}",
                e
            ))
            .unwrap()
            .into_raw()
        }
    };

    let res = FfiResult {
        ok: Some(
            input
                .mv1
                .dot(input.mv2),
        ),
        err: None::<String>,
    };

    CString::new(serde_json::to_string(&res).unwrap())
        .unwrap()
        .into_raw()
}

/// JSON FFI for ga_reverse.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_ga_reverse_json(json_ptr: *const c_char) -> *mut c_char {

    let json_str = match CStr::from_ptr(json_ptr).to_str() {
        Ok(s) => s,
        Err(_) => return std::ptr::null_mut(),
    };

    let input: GaInput = match serde_json::from_str(json_str) {
        Ok(v) => v,
        Err(e) => {
            return CString::new(format!(
                "{{\"err\": \"{}\"}}",
                e
            ))
            .unwrap()
            .into_raw()
        }
    };

    let res = FfiResult {
        ok: Some(input.mv.reverse()),
        err: None::<String>,
    };

    CString::new(serde_json::to_string(&res).unwrap())
        .unwrap()
        .into_raw()
}

/// JSON FFI for ga_norm.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_ga_norm_json(json_ptr: *const c_char) -> *mut c_char {

    let json_str = match CStr::from_ptr(json_ptr).to_str() {
        Ok(s) => s,
        Err(_) => return std::ptr::null_mut(),
    };

    let input: GaInput = match serde_json::from_str(json_str) {
        Ok(v) => v,
        Err(e) => {
            return CString::new(format!(
                "{{\"err\": \"{}\"}}",
                e
            ))
            .unwrap()
            .into_raw()
        }
    };

    let res = FfiResult {
        ok: Some(input.mv.norm()),
        err: None::<String>,
    };

    CString::new(serde_json::to_string(&res).unwrap())
        .unwrap()
        .into_raw()
}

/// JSON FFI for ga_inv.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_ga_inv_json(json_ptr: *const c_char) -> *mut c_char {

    let json_str = match CStr::from_ptr(json_ptr).to_str() {
        Ok(s) => s,
        Err(_) => return std::ptr::null_mut(),
    };

    let input: GaInput = match serde_json::from_str(json_str) {
        Ok(v) => v,
        Err(e) => {
            return CString::new(format!(
                "{{\"err\": \"{}\"}}",
                e
            ))
            .unwrap()
            .into_raw()
        }
    };

    let res = match input.mv.inv() {
        Some(v) => {
            FfiResult {
                ok: Some(v),
                err: None::<String>,
            }
        }
        None => {
            FfiResult {
                ok: None,
                err: Some("Multivector is not invertible".to_string()),
            }
        }
    };

    CString::new(serde_json::to_string(&res).unwrap())
        .unwrap()
        .into_raw()
}

//! JSON-based FFI API for numerical number theory operations.

use std::ffi::CStr;
use std::ffi::CString;
use std::os::raw::c_char;

use serde::Deserialize;
use serde::Serialize;

use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::number_theory as nt;

#[derive(Deserialize)]

struct FactorizeRequest {
    n : u64,
}

/// Factorizes a number from JSON.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_nt_factorize_json(json_ptr : *const c_char) -> *mut c_char {

    if json_ptr.is_null() {

        return std::ptr::null_mut();
    }

    let json_str = match unsafe {

        CStr::from_ptr(json_ptr).to_str()
    } {
        | Ok(s) => s,
        | Err(_) => return std::ptr::null_mut(),
    };

    let req : FactorizeRequest = match serde_json::from_str(json_str) {
        | Ok(r) => r,
        | Err(e) => {

            let res : FfiResult<Vec<u64>, String> = FfiResult {
                ok : None,
                err : Some(e.to_string()),
            };

            return CString::new(serde_json::to_string(&res).unwrap())
                .unwrap()
                .into_raw();
        },
    };

    let factors = nt::factorize(req.n);

    let ffi_res : FfiResult<Vec<u64>, String> = FfiResult {
        ok : Some(factors),
        err : None,
    };

    CString::new(serde_json::to_string(&ffi_res).unwrap())
        .unwrap()
        .into_raw()
}

#[derive(Deserialize)]

struct ModInverseRequest {
    a : i64,
    m : i64,
}

/// Modular inverse from JSON.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_nt_mod_inverse_json(json_ptr : *const c_char) -> *mut c_char {

    if json_ptr.is_null() {

        return std::ptr::null_mut();
    }

    let json_str = match unsafe {

        CStr::from_ptr(json_ptr).to_str()
    } {
        | Ok(s) => s,
        | Err(_) => return std::ptr::null_mut(),
    };

    let req : ModInverseRequest = match serde_json::from_str(json_str) {
        | Ok(r) => r,
        | Err(e) => {

            let res : FfiResult<i64, String> = FfiResult {
                ok : None,
                err : Some(e.to_string()),
            };

            return CString::new(serde_json::to_string(&res).unwrap())
                .unwrap()
                .into_raw();
        },
    };

    match nt::mod_inverse(req.a, req.m) {
        | Some(inv) => {

            let ffi_res : FfiResult<i64, String> = FfiResult {
                ok : Some(inv),
                err : None,
            };

            CString::new(serde_json::to_string(&ffi_res).unwrap())
                .unwrap()
                .into_raw()
        },
        | None => {

            let ffi_res : FfiResult<i64, String> = FfiResult {
                ok : None,
                err : Some("No modular inverse exists".to_string()),
            };

            CString::new(serde_json::to_string(&ffi_res).unwrap())
                .unwrap()
                .into_raw()
        },
    }
}

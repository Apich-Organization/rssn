//! JSON-based FFI API for numerical finite field arithmetic.

use std::ffi::CStr;
use std::ffi::CString;
use std::os::raw::c_char;

use serde::Deserialize;
use serde::Serialize;

use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::finite_field::PrimeFieldElement;
use crate::numerical::finite_field::{
    self,
};

#[derive(Deserialize)]

struct PfeBinaryOpRequest {
    a : PrimeFieldElement,
    b : PrimeFieldElement,
}

/// GF(p) addition from JSON.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_ff_pfe_add_json(
    json_ptr : *const c_char
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

    let req : PfeBinaryOpRequest =
        match serde_json::from_str(
            json_str,
        ) {
            | Ok(r) => r,
            | Err(e) => {

                let res : FfiResult<
                    PrimeFieldElement,
                    String,
                > = FfiResult {
                    ok : None,
                    err : Some(
                        e.to_string(),
                    ),
                };

                return CString::new(
                    serde_json::to_string(&res).unwrap(),
                )
                .unwrap()
                .into_raw();
            },
        };

    let res_pfe = req.a + req.b;

    let ffi_res : FfiResult<
        PrimeFieldElement,
        String,
    > = FfiResult {
        ok : Some(res_pfe),
        err : None,
    };

    CString::new(
        serde_json::to_string(&ffi_res)
            .unwrap(),
    )
    .unwrap()
    .into_raw()
}

/// GF(p) multiplication from JSON.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_ff_pfe_mul_json(
    json_ptr : *const c_char
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

    let req : PfeBinaryOpRequest =
        match serde_json::from_str(
            json_str,
        ) {
            | Ok(r) => r,
            | Err(e) => {

                let res : FfiResult<
                    PrimeFieldElement,
                    String,
                > = FfiResult {
                    ok : None,
                    err : Some(
                        e.to_string(),
                    ),
                };

                return CString::new(
                    serde_json::to_string(&res).unwrap(),
                )
                .unwrap()
                .into_raw();
            },
        };

    let res_pfe = req.a * req.b;

    let ffi_res : FfiResult<
        PrimeFieldElement,
        String,
    > = FfiResult {
        ok : Some(res_pfe),
        err : None,
    };

    CString::new(
        serde_json::to_string(&ffi_res)
            .unwrap(),
    )
    .unwrap()
    .into_raw()
}

#[derive(Deserialize)]

struct Gf256OpRequest {
    a : u8,
    b : u8,
}

/// GF(2^8) multiplication from JSON.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_ff_gf256_mul_json(
    json_ptr : *const c_char
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

    let req : Gf256OpRequest =
        match serde_json::from_str(
            json_str,
        ) {
            | Ok(r) => r,
            | Err(e) => {

                let res : FfiResult<
                    u8,
                    String,
                > = FfiResult {
                    ok : None,
                    err : Some(
                        e.to_string(),
                    ),
                };

                return CString::new(
                    serde_json::to_string(&res).unwrap(),
                )
                .unwrap()
                .into_raw();
            },
        };

    let res = finite_field::gf256_mul(
        req.a, req.b,
    );

    let ffi_res : FfiResult<
        u8,
        String,
    > = FfiResult {
        ok : Some(res),
        err : None,
    };

    CString::new(
        serde_json::to_string(&ffi_res)
            .unwrap(),
    )
    .unwrap()
    .into_raw()
}

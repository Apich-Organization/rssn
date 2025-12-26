//! JSON-based FFI API for numerical real root finding.

use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::polynomial::Polynomial;
use crate::numerical::real_roots;
use serde::Deserialize;
use std::ffi::{CStr, CString};
use std::os::raw::c_char;

#[derive(Deserialize)]

struct FindRootsInput {
    coeffs: Vec<f64>,
    tolerance: f64,
}

#[no_mangle]

pub unsafe extern "C" fn rssn_real_roots_find_roots_json(json_ptr: *const c_char) -> *mut c_char {

    let json_str = match CStr::from_ptr(json_ptr).to_str() {
        Ok(s) => s,
        Err(_) => return std::ptr::null_mut(),
    };

    let input: FindRootsInput = match serde_json::from_str(json_str) {
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

    let poly = Polynomial::new(input.coeffs);

    let result = real_roots::find_roots(
        &poly,
        input.tolerance,
    );

    let res = match result {
        Ok(roots) => {
            FfiResult {
                ok: Some(roots),
                err: None::<String>,
            }
        }
        Err(e) => {
            FfiResult {
                ok: None,
                err: Some(e),
            }
        }
    };

    CString::new(serde_json::to_string(&res).unwrap())
        .unwrap()
        .into_raw()
}

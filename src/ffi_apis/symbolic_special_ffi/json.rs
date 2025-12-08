use crate::symbolic::special;
use crate::ffi_apis::common::*;
use std::os::raw::c_char;

#[no_mangle]
pub unsafe extern "C" fn rssn_json_gamma_numerical(x_json: *const c_char) -> *mut c_char {
    let x: Option<f64> = from_json_string(x_json);
    if let Some(val) = x {
        to_json_string(&special::gamma_numerical(val))
    } else {
        std::ptr::null_mut()
    }
}

#[no_mangle]
pub unsafe extern "C" fn rssn_json_ln_gamma_numerical(x_json: *const c_char) -> *mut c_char {
    let x: Option<f64> = from_json_string(x_json);
    if let Some(val) = x {
        to_json_string(&special::ln_gamma_numerical(val))
    } else {
        std::ptr::null_mut()
    }
}

#[no_mangle]
pub unsafe extern "C" fn rssn_json_beta_numerical(a_json: *const c_char, b_json: *const c_char) -> *mut c_char {
    let a: Option<f64> = from_json_string(a_json);
    let b: Option<f64> = from_json_string(b_json);
    if let (Some(val_a), Some(val_b)) = (a, b) {
        to_json_string(&special::beta_numerical(val_a, val_b))
    } else {
        std::ptr::null_mut()
    }
}

#[no_mangle]
pub unsafe extern "C" fn rssn_json_ln_beta_numerical(a_json: *const c_char, b_json: *const c_char) -> *mut c_char {
    let a: Option<f64> = from_json_string(a_json);
    let b: Option<f64> = from_json_string(b_json);
    if let (Some(val_a), Some(val_b)) = (a, b) {
        to_json_string(&special::ln_beta_numerical(val_a, val_b))
    } else {
        std::ptr::null_mut()
    }
}

#[no_mangle]
pub unsafe extern "C" fn rssn_json_erf_numerical(x_json: *const c_char) -> *mut c_char {
    let x: Option<f64> = from_json_string(x_json);
    if let Some(val) = x {
        to_json_string(&special::erf_numerical(val))
    } else {
        std::ptr::null_mut()
    }
}

#[no_mangle]
pub unsafe extern "C" fn rssn_json_erfc_numerical(x_json: *const c_char) -> *mut c_char {
    let x: Option<f64> = from_json_string(x_json);
    if let Some(val) = x {
        to_json_string(&special::erfc_numerical(val))
    } else {
        std::ptr::null_mut()
    }
}

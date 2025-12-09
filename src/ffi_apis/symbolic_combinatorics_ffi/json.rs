use crate::ffi_apis::common::*;
use crate::symbolic::combinatorics::*;
use crate::symbolic::core::Expr;
use std::os::raw::c_char;

#[no_mangle]
pub unsafe extern "C" fn rssn_json_permutations(
    n_json: *const c_char,
    k_json: *const c_char,
) -> *mut c_char {
    let n: Expr = match from_json_string(n_json) {
        Some(e) => e,
        None => return std::ptr::null_mut(),
    };
    let k: Expr = match from_json_string(k_json) {
        Some(e) => e,
        None => return std::ptr::null_mut(),
    };
    let result = permutations(n, k);
    to_json_string(&result)
}

#[no_mangle]
pub unsafe extern "C" fn rssn_json_combinations(
    n_json: *const c_char,
    k_json: *const c_char,
) -> *mut c_char {
    let n: Expr = match from_json_string(n_json) {
        Some(e) => e,
        None => return std::ptr::null_mut(),
    };
    let k: Expr = match from_json_string(k_json) {
        Some(e) => e,
        None => return std::ptr::null_mut(),
    };
    let result = combinations(&n, k);
    to_json_string(&result)
}

#[no_mangle]
pub unsafe extern "C" fn rssn_json_catalan_number(n: usize) -> *mut c_char {
    let result = catalan_number(n);
    to_json_string(&result)
}

#[no_mangle]
pub unsafe extern "C" fn rssn_json_stirling_number_second_kind(n: usize, k: usize) -> *mut c_char {
    let result = stirling_number_second_kind(n, k);
    to_json_string(&result)
}

#[no_mangle]
pub unsafe extern "C" fn rssn_json_bell_number(n: usize) -> *mut c_char {
    let result = bell_number(n);
    to_json_string(&result)
}

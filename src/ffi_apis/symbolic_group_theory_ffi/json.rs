use crate::ffi_apis::common::*;
use crate::symbolic::core::Expr;
use crate::symbolic::group_theory::*;
use std::collections::HashMap;
use std::os::raw::c_char;

#[no_mangle]
pub unsafe extern "C" fn rssn_json_group_create(
    json_str: *const c_char,
) -> *mut c_char {
    let group: Group = match from_json_string(json_str) {
        Some(g) => g,
        None => return std::ptr::null_mut(),
    };
    to_json_string(&group)
}

#[no_mangle]
pub unsafe extern "C" fn rssn_json_group_multiply(
    group_json: *const c_char,
    a_json: *const c_char,
    b_json: *const c_char,
) -> *mut c_char {
    let group: Group = match from_json_string(group_json) {
        Some(g) => g,
        None => return std::ptr::null_mut(),
    };
    let a: GroupElement = match from_json_string(a_json) {
        Some(e) => e,
        None => return std::ptr::null_mut(),
    };
    let b: GroupElement = match from_json_string(b_json) {
        Some(e) => e,
        None => return std::ptr::null_mut(),
    };
    
    let result = group.multiply(&a, &b);
    to_json_string(&result)
}

#[no_mangle]
pub unsafe extern "C" fn rssn_json_group_inverse(
    group_json: *const c_char,
    a_json: *const c_char,
) -> *mut c_char {
    let group: Group = match from_json_string(group_json) {
        Some(g) => g,
        None => return std::ptr::null_mut(),
    };
    let a: GroupElement = match from_json_string(a_json) {
        Some(e) => e,
        None => return std::ptr::null_mut(),
    };
    
    let result = group.inverse(&a);
    to_json_string(&result)
}

#[no_mangle]
pub unsafe extern "C" fn rssn_json_group_is_abelian(
    group_json: *const c_char,
) -> bool {
    let group: Group = match from_json_string(group_json) {
        Some(g) => g,
        None => return false,
    };
    group.is_abelian()
}

#[no_mangle]
pub unsafe extern "C" fn rssn_json_group_element_order(
    group_json: *const c_char,
    a_json: *const c_char,
) -> usize {
    let group: Group = match from_json_string(group_json) {
        Some(g) => g,
        None => return 0,
    };
    let a: GroupElement = match from_json_string(a_json) {
        Some(e) => e,
        None => return 0,
    };
    group.element_order(&a).unwrap_or(0)
}

#[no_mangle]
pub unsafe extern "C" fn rssn_json_group_conjugacy_classes(
    group_json: *const c_char,
) -> *mut c_char {
    let group: Group = match from_json_string(group_json) {
        Some(g) => g,
        None => return std::ptr::null_mut(),
    };
    let classes = group.conjugacy_classes();
    to_json_string(&classes)
}

#[no_mangle]
pub unsafe extern "C" fn rssn_json_group_center(
    group_json: *const c_char,
) -> *mut c_char {
    let group: Group = match from_json_string(group_json) {
        Some(g) => g,
        None => return std::ptr::null_mut(),
    };
    let center = group.center();
    to_json_string(&center)
}

#[no_mangle]
pub unsafe extern "C" fn rssn_json_representation_create(
    json_str: *const c_char,
) -> *mut c_char {
    let rep: Representation = match from_json_string(json_str) {
        Some(r) => r,
        None => return std::ptr::null_mut(),
    };
    to_json_string(&rep)
}

#[no_mangle]
pub unsafe extern "C" fn rssn_json_representation_is_valid(
    rep_json: *const c_char,
    group_json: *const c_char,
) -> bool {
    let rep: Representation = match from_json_string(rep_json) {
        Some(r) => r,
        None => return false,
    };
    let group: Group = match from_json_string(group_json) {
        Some(g) => g,
        None => return false,
    };
    rep.is_valid(&group)
}

#[no_mangle]
pub unsafe extern "C" fn rssn_json_character(
    rep_json: *const c_char,
) -> *mut c_char {
    let rep: Representation = match from_json_string(rep_json) {
        Some(r) => r,
        None => return std::ptr::null_mut(),
    };
    let chars = character(&rep);
    to_json_string(&chars)
}

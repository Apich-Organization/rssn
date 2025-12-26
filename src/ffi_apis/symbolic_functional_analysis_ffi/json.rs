use crate::ffi_apis::common::*;
use crate::symbolic::core::Expr;
use crate::symbolic::functional_analysis::*;
use std::os::raw::c_char;

#[no_mangle]

pub unsafe extern "C" fn rssn_json_hilbert_space_create(
    json_str: *const c_char
) -> *mut c_char {

    let space: HilbertSpace =
        match from_json_string(json_str) {
            | Some(s) => s,
            | None => return std::ptr::null_mut(),
        };

    to_json_string(&space)
}

#[no_mangle]

pub unsafe extern "C" fn rssn_json_inner_product(
    space_json: *const c_char,
    f_json: *const c_char,
    g_json: *const c_char,
) -> *mut c_char {

    let space: HilbertSpace =
        match from_json_string(space_json) {
            | Some(s) => s,
            | None => return std::ptr::null_mut(),
        };

    let f: Expr = match from_json_string(
        f_json,
    ) {
        | Some(e) => e,
        | None => {
            return std::ptr::null_mut()
        },
    };

    let g: Expr = match from_json_string(
        g_json,
    ) {
        | Some(e) => e,
        | None => {
            return std::ptr::null_mut()
        },
    };

    let result =
        inner_product(&space, &f, &g);

    to_json_string(&result)
}

#[no_mangle]

pub unsafe extern "C" fn rssn_json_norm(
    space_json: *const c_char,
    f_json: *const c_char,
) -> *mut c_char {

    let space: HilbertSpace =
        match from_json_string(space_json) {
            | Some(s) => s,
            | None => return std::ptr::null_mut(),
        };

    let f: Expr = match from_json_string(
        f_json,
    ) {
        | Some(e) => e,
        | None => {
            return std::ptr::null_mut()
        },
    };

    let result = norm(&space, &f);

    to_json_string(&result)
}

#[no_mangle]

pub unsafe extern "C" fn rssn_json_gram_schmidt(
    space_json: *const c_char,
    basis_json: *const c_char,
) -> *mut c_char {

    let space: HilbertSpace =
        match from_json_string(space_json) {
            | Some(s) => s,
            | None => return std::ptr::null_mut(),
        };

    let basis: Vec<Expr> =
        match from_json_string(basis_json) {
            | Some(b) => b,
            | None => return std::ptr::null_mut(),
        };

    let result =
        gram_schmidt(&space, &basis);

    to_json_string(&result)
}

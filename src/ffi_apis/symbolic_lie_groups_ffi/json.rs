use std::os::raw::c_char;

use crate::ffi_apis::common::*;
use crate::symbolic::core::Expr;
use crate::symbolic::lie_groups_and_algebras::*;

// --- LieAlgebra Creation ---

#[no_mangle]

pub unsafe extern "C" fn rssn_json_lie_algebra_so3() -> *mut c_char {

    let algebra = so3();

    to_json_string(&algebra)
}

#[no_mangle]

pub unsafe extern "C" fn rssn_json_lie_algebra_su2() -> *mut c_char {

    let algebra = su2();

    to_json_string(&algebra)
}

// --- Lie Bracket ---

#[no_mangle]

pub unsafe extern "C" fn rssn_json_lie_bracket(
    x_json : *const c_char,
    y_json : *const c_char,
) -> *mut c_char {

    let x : Expr = match from_json_string(x_json) {
        | Some(e) => e,
        | None => return std::ptr::null_mut(),
    };

    let y : Expr = match from_json_string(y_json) {
        | Some(e) => e,
        | None => return std::ptr::null_mut(),
    };

    match lie_bracket(&x, &y) {
        | Ok(result) => to_json_string(&result),
        | Err(_) => std::ptr::null_mut(),
    }
}

// --- Exponential Map ---

#[no_mangle]

pub unsafe extern "C" fn rssn_json_exponential_map(
    x_json : *const c_char,
    order : usize,
) -> *mut c_char {

    let x : Expr = match from_json_string(x_json) {
        | Some(e) => e,
        | None => return std::ptr::null_mut(),
    };

    match exponential_map(&x, order) {
        | Ok(result) => to_json_string(&result),
        | Err(_) => std::ptr::null_mut(),
    }
}

// --- Adjoint Representations ---

#[no_mangle]

pub unsafe extern "C" fn rssn_json_adjoint_representation_group(
    g_json : *const c_char,
    x_json : *const c_char,
) -> *mut c_char {

    let g : Expr = match from_json_string(g_json) {
        | Some(e) => e,
        | None => return std::ptr::null_mut(),
    };

    let x : Expr = match from_json_string(x_json) {
        | Some(e) => e,
        | None => return std::ptr::null_mut(),
    };

    match adjoint_representation_group(&g, &x) {
        | Ok(result) => to_json_string(&result),
        | Err(_) => std::ptr::null_mut(),
    }
}

#[no_mangle]

pub unsafe extern "C" fn rssn_json_adjoint_representation_algebra(
    x_json : *const c_char,
    y_json : *const c_char,
) -> *mut c_char {

    let x : Expr = match from_json_string(x_json) {
        | Some(e) => e,
        | None => return std::ptr::null_mut(),
    };

    let y : Expr = match from_json_string(y_json) {
        | Some(e) => e,
        | None => return std::ptr::null_mut(),
    };

    match adjoint_representation_algebra(&x, &y) {
        | Ok(result) => to_json_string(&result),
        | Err(_) => std::ptr::null_mut(),
    }
}

// --- Commutator Table ---

#[no_mangle]

pub unsafe extern "C" fn rssn_json_commutator_table(algebra_json : *const c_char) -> *mut c_char {

    let algebra : LieAlgebra = match from_json_string(algebra_json) {
        | Some(a) => a,
        | None => return std::ptr::null_mut(),
    };

    match commutator_table(&algebra) {
        | Ok(table) => to_json_string(&table),
        | Err(_) => std::ptr::null_mut(),
    }
}

// --- Jacobi Identity Check ---

#[no_mangle]

pub unsafe extern "C" fn rssn_json_check_jacobi_identity(algebra_json : *const c_char) -> bool {

    let algebra : LieAlgebra = match from_json_string(algebra_json) {
        | Some(a) => a,
        | None => return false,
    };

    match check_jacobi_identity(&algebra) {
        | Ok(result) => result,
        | Err(_) => false,
    }
}

// --- Generators ---

#[no_mangle]

pub unsafe extern "C" fn rssn_json_so3_generators() -> *mut c_char {

    let generators = so3_generators();

    let exprs : Vec<Expr> = generators
        .into_iter()
        .map(|g| g.0)
        .collect();

    to_json_string(&exprs)
}

#[no_mangle]

pub unsafe extern "C" fn rssn_json_su2_generators() -> *mut c_char {

    let generators = su2_generators();

    let exprs : Vec<Expr> = generators
        .into_iter()
        .map(|g| g.0)
        .collect();

    to_json_string(&exprs)
}

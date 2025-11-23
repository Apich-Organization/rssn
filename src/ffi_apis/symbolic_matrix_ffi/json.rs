use crate::ffi_apis::common::{from_json_string, rssn_free_string, to_json_string};
use crate::symbolic::core::Expr;
use crate::symbolic::matrix::*;
use std::ffi::{c_char, CStr, CString};

#[no_mangle]
pub extern "C" fn rssn_matrix_add_json(m1_json: *const c_char, m2_json: *const c_char) -> *mut c_char {
    let m1: Expr = from_json_string(m1_json);
    let m2: Expr = from_json_string(m2_json);
    let result = add_matrices(&m1, &m2);
    to_json_string(&result)
}

#[no_mangle]
pub extern "C" fn rssn_matrix_mul_json(m1_json: *const c_char, m2_json: *const c_char) -> *mut c_char {
    let m1: Expr = from_json_string(m1_json);
    let m2: Expr = from_json_string(m2_json);
    let result = mul_matrices(&m1, &m2);
    to_json_string(&result)
}

#[no_mangle]
pub extern "C" fn rssn_matrix_transpose_json(matrix_json: *const c_char) -> *mut c_char {
    let matrix: Expr = from_json_string(matrix_json);
    let result = transpose_matrix(&matrix);
    to_json_string(&result)
}

#[no_mangle]
pub extern "C" fn rssn_matrix_determinant_json(matrix_json: *const c_char) -> *mut c_char {
    let matrix: Expr = from_json_string(matrix_json);
    let result = determinant(&matrix);
    to_json_string(&result)
}

#[no_mangle]
pub extern "C" fn rssn_matrix_inverse_json(matrix_json: *const c_char) -> *mut c_char {
    let matrix: Expr = from_json_string(matrix_json);
    let result = inverse_matrix(&matrix);
    to_json_string(&result)
}

#[no_mangle]
pub extern "C" fn rssn_matrix_solve_linear_system_json(a_json: *const c_char, b_json: *const c_char) -> *mut c_char {
    let a: Expr = from_json_string(a_json);
    let b: Expr = from_json_string(b_json);
    match solve_linear_system(&a, &b) {
        Ok(result) => to_json_string(&result),
        Err(e) => to_json_string(&Expr::Variable(format!("Error: {}", e))),
    }
}

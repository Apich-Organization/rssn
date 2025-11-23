use crate::ffi_apis::common::*;
use crate::symbolic::core::Expr;
use crate::symbolic::matrix::*;
use std::ffi::c_char;

#[no_mangle]
pub extern "C" fn rssn_json_matrix_add(m1_json: *const c_char, m2_json: *const c_char) -> *mut c_char {
    let m1: Option<Expr> = from_json_string(m1_json);
    let m2: Option<Expr> = from_json_string(m2_json);
    
    if let (Some(matrix1), Some(matrix2)) = (m1, m2) {
        let result = add_matrices(&matrix1, &matrix2);
        to_json_string(&result)
    } else {
        std::ptr::null_mut()
    }
}

#[no_mangle]
pub extern "C" fn rssn_json_matrix_mul(m1_json: *const c_char, m2_json: *const c_char) -> *mut c_char {
    let m1: Option<Expr> = from_json_string(m1_json);
    let m2: Option<Expr> = from_json_string(m2_json);
    
    if let (Some(matrix1), Some(matrix2)) = (m1, m2) {
        let result = mul_matrices(&matrix1, &matrix2);
        to_json_string(&result)
    } else {
        std::ptr::null_mut()
    }
}

#[no_mangle]
pub extern "C" fn rssn_json_matrix_transpose(matrix_json: *const c_char) -> *mut c_char {
    let matrix: Option<Expr> = from_json_string(matrix_json);
    
    if let Some(m) = matrix {
        let result = transpose_matrix(&m);
        to_json_string(&result)
    } else {
        std::ptr::null_mut()
    }
}

#[no_mangle]
pub extern "C" fn rssn_json_matrix_determinant(matrix_json: *const c_char) -> *mut c_char {
    let matrix: Option<Expr> = from_json_string(matrix_json);
    
    if let Some(m) = matrix {
        let result = determinant(&m);
        to_json_string(&result)
    } else {
        std::ptr::null_mut()
    }
}

#[no_mangle]
pub extern "C" fn rssn_json_matrix_inverse(matrix_json: *const c_char) -> *mut c_char {
    let matrix: Option<Expr> = from_json_string(matrix_json);
    
    if let Some(m) = matrix {
        let result = inverse_matrix(&m);
        to_json_string(&result)
    } else {
        std::ptr::null_mut()
    }
}

#[no_mangle]
pub extern "C" fn rssn_json_matrix_solve_linear_system(a_json: *const c_char, b_json: *const c_char) -> *mut c_char {
    let a: Option<Expr> = from_json_string(a_json);
    let b: Option<Expr> = from_json_string(b_json);
    
    if let (Some(matrix_a), Some(vector_b)) = (a, b) {
        match solve_linear_system(&matrix_a, &vector_b) {
            Ok(result) => to_json_string(&result),
            Err(_) => std::ptr::null_mut(),
        }
    } else {
        std::ptr::null_mut()
    }
}

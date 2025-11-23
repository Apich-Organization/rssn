use crate::ffi_apis::common::{from_bincode_buffer, rssn_free_bincode_buffer, to_bincode_buffer, BincodeBuffer};
use crate::symbolic::core::Expr;
use crate::symbolic::matrix::*;

#[no_mangle]
pub extern "C" fn rssn_matrix_add_bincode(m1_buf: BincodeBuffer, m2_buf: BincodeBuffer) -> BincodeBuffer {
    let m1: Expr = from_bincode_buffer(m1_buf);
    let m2: Expr = from_bincode_buffer(m2_buf);
    let result = add_matrices(&m1, &m2);
    to_bincode_buffer(&result)
}

#[no_mangle]
pub extern "C" fn rssn_matrix_mul_bincode(m1_buf: BincodeBuffer, m2_buf: BincodeBuffer) -> BincodeBuffer {
    let m1: Expr = from_bincode_buffer(m1_buf);
    let m2: Expr = from_bincode_buffer(m2_buf);
    let result = mul_matrices(&m1, &m2);
    to_bincode_buffer(&result)
}

#[no_mangle]
pub extern "C" fn rssn_matrix_transpose_bincode(matrix_buf: BincodeBuffer) -> BincodeBuffer {
    let matrix: Expr = from_bincode_buffer(matrix_buf);
    let result = transpose_matrix(&matrix);
    to_bincode_buffer(&result)
}

#[no_mangle]
pub extern "C" fn rssn_matrix_determinant_bincode(matrix_buf: BincodeBuffer) -> BincodeBuffer {
    let matrix: Expr = from_bincode_buffer(matrix_buf);
    let result = determinant(&matrix);
    to_bincode_buffer(&result)
}

#[no_mangle]
pub extern "C" fn rssn_matrix_inverse_bincode(matrix_buf: BincodeBuffer) -> BincodeBuffer {
    let matrix: Expr = from_bincode_buffer(matrix_buf);
    let result = inverse_matrix(&matrix);
    to_bincode_buffer(&result)
}

#[no_mangle]
pub extern "C" fn rssn_matrix_solve_linear_system_bincode(a_buf: BincodeBuffer, b_buf: BincodeBuffer) -> BincodeBuffer {
    let a: Expr = from_bincode_buffer(a_buf);
    let b: Expr = from_bincode_buffer(b_buf);
    match solve_linear_system(&a, &b) {
        Ok(result) => to_bincode_buffer(&result),
        Err(e) => to_bincode_buffer(&Expr::Variable(format!("Error: {}", e))),
    }
}

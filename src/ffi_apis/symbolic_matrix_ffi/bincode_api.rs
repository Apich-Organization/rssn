use crate::ffi_apis::common::*;
use crate::symbolic::core::Expr;
use crate::symbolic::matrix::*;

#[no_mangle]
pub extern "C" fn rssn_bincode_matrix_add(
    m1_buf: BincodeBuffer,
    m2_buf: BincodeBuffer,
) -> BincodeBuffer {
    let m1: Option<Expr> = from_bincode_buffer(&m1_buf);
    let m2: Option<Expr> = from_bincode_buffer(&m2_buf);

    if let (Some(matrix1), Some(matrix2)) = (m1, m2) {
        let result = add_matrices(&matrix1, &matrix2);
        to_bincode_buffer(&result)
    } else {
        BincodeBuffer::empty()
    }
}

#[no_mangle]
pub extern "C" fn rssn_bincode_matrix_mul(
    m1_buf: BincodeBuffer,
    m2_buf: BincodeBuffer,
) -> BincodeBuffer {
    let m1: Option<Expr> = from_bincode_buffer(&m1_buf);
    let m2: Option<Expr> = from_bincode_buffer(&m2_buf);

    if let (Some(matrix1), Some(matrix2)) = (m1, m2) {
        let result = mul_matrices(&matrix1, &matrix2);
        to_bincode_buffer(&result)
    } else {
        BincodeBuffer::empty()
    }
}

#[no_mangle]
pub extern "C" fn rssn_bincode_matrix_transpose(matrix_buf: BincodeBuffer) -> BincodeBuffer {
    let matrix: Option<Expr> = from_bincode_buffer(&matrix_buf);

    if let Some(m) = matrix {
        let result = transpose_matrix(&m);
        to_bincode_buffer(&result)
    } else {
        BincodeBuffer::empty()
    }
}

#[no_mangle]
pub extern "C" fn rssn_bincode_matrix_determinant(matrix_buf: BincodeBuffer) -> BincodeBuffer {
    let matrix: Option<Expr> = from_bincode_buffer(&matrix_buf);

    if let Some(m) = matrix {
        let result = determinant(&m);
        to_bincode_buffer(&result)
    } else {
        BincodeBuffer::empty()
    }
}

#[no_mangle]
pub extern "C" fn rssn_bincode_matrix_inverse(matrix_buf: BincodeBuffer) -> BincodeBuffer {
    let matrix: Option<Expr> = from_bincode_buffer(&matrix_buf);

    if let Some(m) = matrix {
        let result = inverse_matrix(&m);
        to_bincode_buffer(&result)
    } else {
        BincodeBuffer::empty()
    }
}

#[no_mangle]
pub extern "C" fn rssn_bincode_matrix_solve_linear_system(
    a_buf: BincodeBuffer,
    b_buf: BincodeBuffer,
) -> BincodeBuffer {
    let a: Option<Expr> = from_bincode_buffer(&a_buf);
    let b: Option<Expr> = from_bincode_buffer(&b_buf);

    if let (Some(matrix_a), Some(vector_b)) = (a, b) {
        match solve_linear_system(&matrix_a, &vector_b) {
            Ok(result) => to_bincode_buffer(&result),
            Err(_) => BincodeBuffer::empty(),
        }
    } else {
        BincodeBuffer::empty()
    }
}

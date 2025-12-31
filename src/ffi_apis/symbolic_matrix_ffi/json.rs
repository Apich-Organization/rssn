use std::ffi::c_char;

use crate::ffi_apis::common::from_json_string;
use crate::ffi_apis::common::to_json_string;
use crate::symbolic::core::Expr;
use crate::symbolic::matrix::add_matrices;
use crate::symbolic::matrix::determinant;
use crate::symbolic::matrix::inverse_matrix;
use crate::symbolic::matrix::mul_matrices;
use crate::symbolic::matrix::solve_linear_system;
use crate::symbolic::matrix::transpose_matrix;

/// Performs matrix addition.

///

/// Takes two JSON strings representing `Expr` (matrices) as input,

/// and returns a JSON string representing their sum.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_json_matrix_add(
    m1_json: *const c_char,
    m2_json: *const c_char,
) -> *mut c_char {

    let m1: Option<Expr> =
        from_json_string(m1_json);

    let m2: Option<Expr> =
        from_json_string(m2_json);

    match (m1, m2) {
        | (
            Some(matrix1),
            Some(matrix2),
        ) => {

            let result = add_matrices(
                &matrix1,
                &matrix2,
            );

            to_json_string(&result)
        },
        | _ => std::ptr::null_mut(),
    }
}

/// Performs matrix multiplication.

///

/// Takes two JSON strings representing `Expr` (matrices) as input,

/// and returns a JSON string representing their product.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_json_matrix_mul(
    m1_json: *const c_char,
    m2_json: *const c_char,
) -> *mut c_char {

    let m1: Option<Expr> =
        from_json_string(m1_json);

    let m2: Option<Expr> =
        from_json_string(m2_json);

    match (m1, m2) {
        | (
            Some(matrix1),
            Some(matrix2),
        ) => {

            let result = mul_matrices(
                &matrix1,
                &matrix2,
            );

            to_json_string(&result)
        },
        | _ => std::ptr::null_mut(),
    }
}

/// Performs matrix transposition.

///

/// Takes a JSON string representing an `Expr` (matrix) as input,

/// and returns a JSON string representing its transpose.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_json_matrix_transpose(
    matrix_json: *const c_char
) -> *mut c_char {

    let matrix: Option<Expr> =
        from_json_string(matrix_json);

    if let Some(m) = matrix {

        let result =
            transpose_matrix(&m);

        to_json_string(&result)
    } else {

        std::ptr::null_mut()
    }
}

/// Computes the determinant of a matrix.

///

/// Takes a JSON string representing an `Expr` (matrix) as input,

/// and returns a JSON string representing its determinant.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_json_matrix_determinant(
    matrix_json: *const c_char
) -> *mut c_char {

    let matrix: Option<Expr> =
        from_json_string(matrix_json);

    if let Some(m) = matrix {

        let result = determinant(&m);

        to_json_string(&result)
    } else {

        std::ptr::null_mut()
    }
}

/// Computes the inverse of a matrix.

///

/// Takes a JSON string representing an `Expr` (matrix) as input,

/// and returns a JSON string representing its inverse.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_json_matrix_inverse(
    matrix_json: *const c_char
) -> *mut c_char {

    let matrix: Option<Expr> =
        from_json_string(matrix_json);

    if let Some(m) = matrix {

        let result = inverse_matrix(&m);

        to_json_string(&result)
    } else {

        std::ptr::null_mut()
    }
}

/// Solves a linear system of equations AX = B.

///

/// Takes two JSON strings representing `Expr` (matrix A and vector B) as input,

/// and returns a JSON string representing the solution vector X.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_json_matrix_solve_linear_system(
    a_json: *const c_char,
    b_json: *const c_char,
) -> *mut c_char {

    let a: Option<Expr> =
        from_json_string(a_json);

    let b: Option<Expr> =
        from_json_string(b_json);

    match (a, b) {
        | (
            Some(matrix_a),
            Some(vector_b),
        ) => {

            match solve_linear_system(
                &matrix_a,
                &vector_b,
            ) {
                | Ok(result) => {
                    to_json_string(
                        &result,
                    )
                },
                | Err(_) => {
                    std::ptr::null_mut()
                },
            }
        },
        | _ => std::ptr::null_mut(),
    }
}

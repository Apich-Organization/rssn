use crate::symbolic::core::Expr;
use crate::symbolic::matrix::add_matrices;
use crate::symbolic::matrix::determinant;
use crate::symbolic::matrix::inverse_matrix;
use crate::symbolic::matrix::mul_matrices;
use crate::symbolic::matrix::solve_linear_system;
use crate::symbolic::matrix::transpose_matrix;

/// Performs matrix addition using raw pointers to `Expr` objects.

///

/// Takes two raw pointers to `Expr` (representing matrices) as input,

/// and returns a raw pointer to a new `Expr` representing their sum.

#[no_mangle]

pub extern "C" fn rssn_matrix_add_handle(
    m1: *const Expr,
    m2: *const Expr,
) -> *mut Expr {

    let m1 = unsafe {

        &*m1
    };

    let m2 = unsafe {

        &*m2
    };

    let result = add_matrices(m1, m2);

    Box::into_raw(Box::new(result))
}

/// Performs matrix multiplication using raw pointers to `Expr` objects.

///

/// Takes two raw pointers to `Expr` (representing matrices) as input,

/// and returns a raw pointer to a new `Expr` representing their product.

#[no_mangle]

pub extern "C" fn rssn_matrix_mul_handle(
    m1: *const Expr,
    m2: *const Expr,
) -> *mut Expr {

    let m1 = unsafe {

        &*m1
    };

    let m2 = unsafe {

        &*m2
    };

    let result = mul_matrices(m1, m2);

    Box::into_raw(Box::new(result))
}

/// Performs matrix transposition using a raw pointer to an `Expr` object.

///

/// Takes a raw pointer to an `Expr` (representing a matrix) as input,

/// and returns a raw pointer to a new `Expr` representing its transpose.

#[no_mangle]

pub extern "C" fn rssn_matrix_transpose_handle(
    matrix: *const Expr
) -> *mut Expr {

    let matrix = unsafe {

        &*matrix
    };

    let result =
        transpose_matrix(matrix);

    Box::into_raw(Box::new(result))
}

/// Computes the determinant of a matrix using a raw pointer to an `Expr` object.

///

/// Takes a raw pointer to an `Expr` (representing a matrix) as input,

/// and returns a raw pointer to a new `Expr` representing its determinant.

#[no_mangle]

pub extern "C" fn rssn_matrix_determinant_handle(
    matrix: *const Expr
) -> *mut Expr {

    let matrix = unsafe {

        &*matrix
    };

    let result = determinant(matrix);

    Box::into_raw(Box::new(result))
}

/// Computes the inverse of a matrix using a raw pointer to an `Expr` object.

///

/// Takes a raw pointer to an `Expr` (representing a matrix) as input,

/// and returns a raw pointer to a new `Expr` representing its inverse.

#[no_mangle]

pub extern "C" fn rssn_matrix_inverse_handle(
    matrix: *const Expr
) -> *mut Expr {

    let matrix = unsafe {

        &*matrix
    };

    let result = inverse_matrix(matrix);

    Box::into_raw(Box::new(result))
}

/// Solves a linear system of equations AX = B using raw pointers to `Expr` objects.

///

/// Takes two raw pointers to `Expr` (representing matrix A and vector B) as input,

/// and returns a raw pointer to a new `Expr` representing the solution vector X.

#[no_mangle]

pub extern "C" fn rssn_matrix_solve_linear_system_handle(
    a: *const Expr,
    b: *const Expr,
) -> *mut Expr {

    let a = unsafe {

        &*a
    };

    let b = unsafe {

        &*b
    };

    match solve_linear_system(a, b) {
        | Ok(result) => {
            Box::into_raw(Box::new(
                result,
            ))
        },
        | Err(e) => {
            Box::into_raw(Box::new(
                Expr::Variable(
                    format!(
                        "Error: {e}"
                    ),
                ),
            ))
        },
    }
}

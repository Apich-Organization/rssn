use crate::symbolic::core::Expr;
use crate::symbolic::matrix::*;

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

#[no_mangle]

pub extern "C" fn rssn_matrix_transpose_handle(matrix: *const Expr) -> *mut Expr {

    let matrix = unsafe {

        &*matrix
    };

    let result = transpose_matrix(matrix);

    Box::into_raw(Box::new(result))
}

#[no_mangle]

pub extern "C" fn rssn_matrix_determinant_handle(matrix: *const Expr) -> *mut Expr {

    let matrix = unsafe {

        &*matrix
    };

    let result = determinant(matrix);

    Box::into_raw(Box::new(result))
}

#[no_mangle]

pub extern "C" fn rssn_matrix_inverse_handle(matrix: *const Expr) -> *mut Expr {

    let matrix = unsafe {

        &*matrix
    };

    let result = inverse_matrix(matrix);

    Box::into_raw(Box::new(result))
}

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
        | Ok(result) => Box::into_raw(Box::new(result)),
        | Err(e) => {
            Box::into_raw(Box::new(
                Expr::Variable(format!(
                    "Error: {}",
                    e
                )),
            ))
        },
    }
}

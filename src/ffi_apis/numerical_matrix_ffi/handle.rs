//! Handle-based FFI API for numerical matrix operations.

use crate::ffi_apis::ffi_api::update_last_error;
use crate::numerical::matrix::Matrix;
use std::ptr;

/// Creates a new f64 matrix from dimensions and a raw data array.
///
/// # Arguments
/// * `rows` - Number of rows.
/// * `cols` - Number of columns.
/// * `data` - Pointer to an array of doubles in row-major order.
///
/// # Returns
/// A raw pointer to the Matrix object, or null on error.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_matrix_create(
    rows: usize,
    cols: usize,
    data: *const f64,
) -> *mut Matrix<f64> {

    if data.is_null() {

        update_last_error(
            "Null pointer passed to \
             rssn_num_matrix_create"
                .to_string(),
        );

        return ptr::null_mut();
    }

    let len = rows * cols;

    let slice = unsafe {

        std::slice::from_raw_parts(
            data, len,
        )
    };

    let matrix = Matrix::new(
        rows,
        cols,
        slice.to_vec(),
    );

    Box::into_raw(Box::new(matrix))
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_matrix_free(
    matrix: *mut Matrix<f64>
) {

    if !matrix.is_null() {

        unsafe {

            let _ =
                Box::from_raw(matrix);
        }
    }
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_matrix_get_rows(
    matrix: *const Matrix<f64>
) -> usize {

    if matrix.is_null() {

        return 0;
    }

    unsafe {

        (*matrix).rows()
    }
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_matrix_get_cols(
    matrix: *const Matrix<f64>
) -> usize {

    if matrix.is_null() {

        return 0;
    }

    unsafe {

        (*matrix).cols()
    }
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_matrix_get_data(
    matrix: *const Matrix<f64>,
    buffer: *mut f64,
) -> i32 {

    if matrix.is_null()
        || buffer.is_null()
    {

        update_last_error(
            "Null pointer passed to \
             rssn_num_matrix_get_data"
                .to_string(),
        );

        return -1;
    }

    let m = unsafe {

        &*matrix
    };

    let data = m.data();

    unsafe {

        ptr::copy_nonoverlapping(
            data.as_ptr(),
            buffer,
            data.len(),
        );
    }

    0
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_matrix_add(
    m1: *const Matrix<f64>,
    m2: *const Matrix<f64>,
) -> *mut Matrix<f64> {

    if m1.is_null() || m2.is_null() {

        return ptr::null_mut();
    }

    let v1 = unsafe {

        &*m1
    };

    let v2 = unsafe {

        &*m2
    };

    if v1.rows() != v2.rows()
        || v1.cols() != v2.cols()
    {

        update_last_error(
            "Dimension mismatch in \
             matrix addition"
                .to_string(),
        );

        return ptr::null_mut();
    }

    let res = v1.clone() + v2.clone();

    Box::into_raw(Box::new(res))
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_matrix_mul(
    m1: *const Matrix<f64>,
    m2: *const Matrix<f64>,
) -> *mut Matrix<f64> {

    if m1.is_null() || m2.is_null() {

        return ptr::null_mut();
    }

    let v1 = unsafe {

        &*m1
    };

    let v2 = unsafe {

        &*m2
    };

    if v1.cols() != v2.rows() {

        update_last_error(
            "Dimension mismatch in \
             matrix multiplication"
                .to_string(),
        );

        return ptr::null_mut();
    }

    let res = v1.clone() * v2.clone();

    Box::into_raw(Box::new(res))
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_matrix_transpose(
    matrix: *const Matrix<f64>
) -> *mut Matrix<f64> {

    if matrix.is_null() {

        return ptr::null_mut();
    }

    let m = unsafe {

        &*matrix
    };

    let res = m.transpose();

    Box::into_raw(Box::new(res))
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_matrix_determinant(
    matrix: *const Matrix<f64>,
    result: *mut f64,
) -> i32 {

    if matrix.is_null()
        || result.is_null()
    {

        return -1;
    }

    let m = unsafe {

        &*matrix
    };

    match m.determinant() {
        | Ok(d) => {

            unsafe {

                *result = d
            };

            0
        },
        | Err(e) => {

            update_last_error(e);

            -1
        },
    }
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_matrix_inverse(
    matrix: *const Matrix<f64>
) -> *mut Matrix<f64> {

    if matrix.is_null() {

        return ptr::null_mut();
    }

    let m = unsafe {

        &*matrix
    };

    match m.inverse() {
        | Some(inv) => {
            Box::into_raw(Box::new(inv))
        },
        | None => {

            update_last_error(
                "Matrix is singular \
                 or not square"
                    .to_string(),
            );

            ptr::null_mut()
        },
    }
}

/// Creates an identity matrix.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_matrix_identity(
    size: usize
) -> *mut Matrix<f64> {

    let m = Matrix::identity(size);

    Box::into_raw(Box::new(m))
}

/// Checks if it's identity.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_matrix_is_identity(
    matrix: *const Matrix<f64>,
    epsilon: f64,
) -> i32 {

    if matrix.is_null() {

        return 0;
    }

    let m = unsafe {

        &*matrix
    };

    if m.is_identity(epsilon) {

        1
    } else {

        0
    }
}

/// Checks if it's orthogonal.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_matrix_is_orthogonal(
    matrix: *const Matrix<f64>,
    epsilon: f64,
) -> i32 {

    if matrix.is_null() {

        return 0;
    }

    let m = unsafe {

        &*matrix
    };

    if m.is_orthogonal(epsilon) {

        1
    } else {

        0
    }
}

/// Returns the rank.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_matrix_rank(
    matrix: *const Matrix<f64>,
    out_rank: *mut usize,
) -> i32 {

    if matrix.is_null()
        || out_rank.is_null()
    {

        return -1;
    }

    let m = unsafe {

        &*matrix
    };

    match m.rank() {
        | Ok(r) => {

            unsafe {

                *out_rank = r
            };

            0
        },
        | Err(e) => {

            update_last_error(e);

            -1
        },
    }
}

/// Returns the trace.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_matrix_trace(
    matrix: *const Matrix<f64>,
    out_trace: *mut f64,
) -> i32 {

    if matrix.is_null()
        || out_trace.is_null()
    {

        return -1;
    }

    let m = unsafe {

        &*matrix
    };

    match m.trace() {
        | Ok(t) => {

            unsafe {

                *out_trace = t
            };

            0
        },
        | Err(e) => {

            update_last_error(e);

            -1
        },
    }
}

/// Returns the Frobenius norm.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_matrix_frobenius_norm(
    matrix: *const Matrix<f64>
) -> f64 {

    if matrix.is_null() {

        return 0.0;
    }

    let m = unsafe {

        &*matrix
    };

    m.frobenius_norm()
}

//! Handle-based FFI API for numerical matrix operations.

use std::ptr;

use crate::ffi_apis::ffi_api::update_last_error;
use crate::nightly::matrix::Backend;
use crate::nightly::matrix::FaerDecompositionResult;
use crate::nightly::matrix::FaerDecompositionType;
use crate::nightly::matrix::Matrix;

/// Opaque handle for a Matrix<f64>.

pub struct RssnMatrixHandle;

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

pub unsafe extern "C" fn rssn_num_matrix_create_nightly(
    rows: usize,
    cols: usize,
    data: *const f64,
) -> *mut RssnMatrixHandle {

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

    Box::into_raw(Box::new(matrix)).cast::<RssnMatrixHandle>()
}

/// Frees a previously allocated Matrix.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_matrix_free_nightly(
    matrix: *mut RssnMatrixHandle
) {

    if !matrix.is_null() {

        unsafe {

            let _ = Box::from_raw(
                matrix.cast::<Matrix<f64>>(),
            );
        }
    }
}

/// Returns the number of rows of a given matrix.

#[no_mangle]

pub const unsafe extern "C" fn rssn_num_matrix_get_rows_nightly(
    matrix: *const RssnMatrixHandle
) -> usize {

    if matrix.is_null() {

        return 0;
    }

    unsafe {

        (*matrix.cast::<Matrix<f64>>())
            .rows()
    }
}

/// Returns the number of columns of a given matrix.

#[no_mangle]

pub const unsafe extern "C" fn rssn_num_matrix_get_cols_nightly(
    matrix: *const RssnMatrixHandle
) -> usize {

    if matrix.is_null() {

        return 0;
    }

    unsafe {

        (*matrix.cast::<Matrix<f64>>())
            .cols()
    }
}

/// Retrieves the raw data of a given matrix.

#[no_mangle]

pub unsafe extern "C" fn rssn_num_matrix_get_data_nightly(
    matrix: *const RssnMatrixHandle,

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

        &*matrix.cast::<Matrix<f64>>()
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

/// Adds two matrices.

#[no_mangle]

pub unsafe extern "C" fn rssn_num_matrix_add_nightly(
    m1: *const RssnMatrixHandle,

    m2: *const RssnMatrixHandle,
) -> *mut RssnMatrixHandle {

    if m1.is_null() || m2.is_null() {

        return ptr::null_mut();
    }

    let v1 = unsafe {

        &*m1.cast::<Matrix<f64>>()
    };

    let v2 = unsafe {

        &*m2.cast::<Matrix<f64>>()
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

    Box::into_raw(Box::new(res)).cast::<RssnMatrixHandle>()
}

/// Multiplies two matrices.

#[no_mangle]

pub unsafe extern "C" fn rssn_num_matrix_mul_nightly(
    m1: *const RssnMatrixHandle,

    m2: *const RssnMatrixHandle,
) -> *mut RssnMatrixHandle {

    if m1.is_null() || m2.is_null() {

        return ptr::null_mut();
    }

    let v1 = unsafe {

        &*m1.cast::<Matrix<f64>>()
    };

    let v2 = unsafe {

        &*m2.cast::<Matrix<f64>>()
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

    Box::into_raw(Box::new(res)).cast::<RssnMatrixHandle>()
}

/// Transposes a matrix.

#[no_mangle]

pub unsafe extern "C" fn rssn_num_matrix_transpose_nightly(
    matrix: *const RssnMatrixHandle
) -> *mut RssnMatrixHandle {

    if matrix.is_null() {

        return ptr::null_mut();
    }

    let m = unsafe {

        &*matrix.cast::<Matrix<f64>>()
    };

    let res = m.transpose();

    Box::into_raw(Box::new(res)).cast::<RssnMatrixHandle>()
}

/// Computes the determinant of a matrix.

#[no_mangle]

pub unsafe extern "C" fn rssn_num_matrix_determinant_nightly(
    matrix: *const RssnMatrixHandle,

    result: *mut f64,
) -> i32 {

    if matrix.is_null()
        || result.is_null()
    {

        return -1;
    }

    let m = unsafe {

        &*matrix.cast::<Matrix<f64>>()
    };

    match m.determinant() {
        | Ok(d) => {

            unsafe {

                *result = d;
            };

            0
        },
        | Err(e) => {

            update_last_error(e);

            -1
        },
    }
}

/// Computes the inverse of a matrix.

#[no_mangle]

pub unsafe extern "C" fn rssn_num_matrix_inverse_nightly(
    matrix: *const RssnMatrixHandle
) -> *mut RssnMatrixHandle {

    if matrix.is_null() {

        return ptr::null_mut();
    }

    let m = unsafe {

        &*matrix.cast::<Matrix<f64>>()
    };

    match m.inverse() {
        | Some(inv) => {
            Box::into_raw(Box::new(inv)).cast::<RssnMatrixHandle>()
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

pub unsafe extern "C" fn rssn_num_matrix_identity_nightly(
    size: usize
) -> *mut RssnMatrixHandle {

    let m =
        Matrix::<f64>::identity(size);

    Box::into_raw(Box::new(m)).cast::<RssnMatrixHandle>()
}

/// Checks if it's identity.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_matrix_is_identity_nightly(
    matrix: *const RssnMatrixHandle,
    epsilon: f64,
) -> i32 {

    if matrix.is_null() {

        return 0;
    }

    let m = unsafe {

        &*matrix.cast::<Matrix<f64>>()
    };

    i32::from(m.is_identity(epsilon))
}

/// Checks if it's orthogonal.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_matrix_is_orthogonal_nightly(
    matrix: *const RssnMatrixHandle,
    epsilon: f64,
) -> i32 {

    if matrix.is_null() {

        return 0;
    }

    let m = unsafe {

        &*matrix.cast::<Matrix<f64>>()
    };

    i32::from(m.is_orthogonal(epsilon))
}

/// Returns the rank.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_matrix_rank_nightly(
    matrix: *const RssnMatrixHandle,
    out_rank: *mut usize,
) -> i32 {

    if matrix.is_null()
        || out_rank.is_null()
    {

        return -1;
    }

    let m = unsafe {

        &*matrix.cast::<Matrix<f64>>()
    };

    match m.rank() {
        | Ok(r) => {

            unsafe {

                *out_rank = r;
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

pub unsafe extern "C" fn rssn_num_matrix_trace_nightly(
    matrix: *const RssnMatrixHandle,
    out_trace: *mut f64,
) -> i32 {

    if matrix.is_null()
        || out_trace.is_null()
    {

        return -1;
    }

    let m = unsafe {

        &*matrix.cast::<Matrix<f64>>()
    };

    match m.trace() {
        | Ok(t) => {

            unsafe {

                *out_trace = t;
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

pub unsafe extern "C" fn rssn_num_matrix_frobenius_norm_nightly(
    matrix: *const RssnMatrixHandle
) -> f64 {

    if matrix.is_null() {

        return 0.0;
    }

    let m = unsafe {

        &*matrix.cast::<Matrix<f64>>()
    };

    m.frobenius_norm()
}

/// Sets the backend for the matrix.
/// 0: Native, 1: Faer
#[no_mangle]

pub unsafe extern "C" fn rssn_num_matrix_set_backend_nightly(
    matrix: *mut RssnMatrixHandle,
    backend_id: i32,
) -> i32 {

    if matrix.is_null() {

        return -1;
    }

    let m = &mut *matrix.cast::<Matrix<f64>>();

    match backend_id {
        | 0 => {
            m.set_backend(
                Backend::Native,
            );
        },
        | 1 => {
            m.set_backend(Backend::Faer);
        },
        | _ => return -1,
    }

    0
}

/// Computes SVD decomposition: A = U * S * V^T.
/// Returns 0 on success, -1 on error.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_matrix_decompose_svd_nightly(
    matrix: *const RssnMatrixHandle,
    out_u: *mut *mut RssnMatrixHandle,
    out_s: *mut *mut RssnMatrixHandle,
    out_v: *mut *mut RssnMatrixHandle,
) -> i32 {

    if matrix.is_null()
        || out_u.is_null()
        || out_s.is_null()
        || out_v.is_null()
    {

        return -1;
    }

    let m = &*matrix.cast::<Matrix<f64>>();

    if let Some(res) = m.decompose(
        FaerDecompositionType::Svd,
    ) {

        if let FaerDecompositionResult::Svd { u, s, v } = res {
            let s_mat = Matrix::new(s.len(), 1, s).with_backend(m.backend);

            *out_u = Box::into_raw(Box::new(u)).cast::<RssnMatrixHandle>();
            *out_s = Box::into_raw(Box::new(s_mat)).cast::<RssnMatrixHandle>();
            *out_v = Box::into_raw(Box::new(v)).cast::<RssnMatrixHandle>();
            return 0;
        }
    }

    -1
}

/// Computes Cholesky decomposition: A = L * L^T.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_matrix_decompose_cholesky_nightly(
    matrix: *const RssnMatrixHandle,
    out_l: *mut *mut RssnMatrixHandle,
) -> i32 {

    if matrix.is_null()
        || out_l.is_null()
    {

        return -1;
    }

    let m = &*matrix.cast::<Matrix<f64>>();

    if let Some(res) = m.decompose(
        FaerDecompositionType::Cholesky,
    ) {

        if let FaerDecompositionResult::Cholesky { l } = res {
             *out_l = Box::into_raw(Box::new(l)).cast::<RssnMatrixHandle>();
             return 0;
        }
    }

    -1
}

/// Computes Symmetric Eigendecomposition: A = V * D * V^T.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_matrix_decompose_eigen_symmetric_nightly(
    matrix: *const RssnMatrixHandle,
    out_values: *mut *mut RssnMatrixHandle,
    out_vectors: *mut *mut RssnMatrixHandle,
) -> i32 {

    if matrix.is_null()
        || out_values.is_null()
        || out_vectors.is_null()
    {

        return -1;
    }

    let m = &*matrix.cast::<Matrix<f64>>();

    if let Some(res) = m.decompose(FaerDecompositionType::EigenSymmetric) {
        if let FaerDecompositionResult::EigenSymmetric { values, vectors } = res {
             let val_mat = Matrix::new(values.len(), 1, values).with_backend(m.backend);
             *out_values = Box::into_raw(Box::new(val_mat)).cast::<RssnMatrixHandle>();
             *out_vectors = Box::into_raw(Box::new(vectors)).cast::<RssnMatrixHandle>();
             return 0;
        }
    }

    -1
}

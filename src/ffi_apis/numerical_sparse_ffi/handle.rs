//! Handle-based FFI API for numerical sparse matrix operations.

use std::ptr;

use sprs_rssn::CsMat;

use crate::ffi_apis::ffi_api::update_last_error;
use crate::numerical::sparse::{
    self,
};

/// Creates a new sparse CSR matrix from dimensions and triplet data.
///
/// # Arguments
/// * `rows` - Number of rows.
/// * `cols` - Number of columns.
/// * `row_indices` - Array of row indices.
/// * `col_indices` - Array of column indices.
/// * `values` - Array of values.
/// * `nnz` - Number of non-zero elements (length of the input arrays).
///
/// # Returns
/// A raw pointer to the `CsMat` object, or null on error.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_sparse_create(
    rows: usize,
    cols: usize,
    row_indices: *const usize,
    col_indices: *const usize,
    values: *const f64,
    nnz: usize,
) -> *mut CsMat<f64> { unsafe {

    if row_indices.is_null()
        || col_indices.is_null()
        || values.is_null()
    {

        update_last_error(
            "Null pointer passed to \
             rssn_num_sparse_create"
                .to_string(),
        );

        return ptr::null_mut();
    }

    let r_idx =  {

        std::slice::from_raw_parts(
            row_indices,
            nnz,
        )
    };

    let c_idx =  {

        std::slice::from_raw_parts(
            col_indices,
            nnz,
        )
    };

    let vals =  {

        std::slice::from_raw_parts(
            values,
            nnz,
        )
    };

    let mut triplets =
        Vec::with_capacity(nnz);

    for i in 0 .. nnz {

        triplets.push((
            r_idx[i],
            c_idx[i],
            vals[i],
        ));
    }

    let mat = sparse::csr_from_triplets(
        rows,
        cols,
        &triplets,
    );

    Box::into_raw(Box::new(mat))
}}

/// Frees a sparse matrix object.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_sparse_free(
    matrix: *mut CsMat<f64>
) {

    if !matrix.is_null() {

        unsafe {

            let _ =
                Box::from_raw(matrix);
        }
    }
}

/// Returns the number of rows.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_sparse_get_rows(
    matrix: *const CsMat<f64>
) -> usize {

    if matrix.is_null() {

        return 0;
    }

    unsafe {

        (*matrix).rows()
    }
}

/// Returns the number of columns.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_sparse_get_cols(
    matrix: *const CsMat<f64>
) -> usize {

    if matrix.is_null() {

        return 0;
    }

    unsafe {

        (*matrix).cols()
    }
}

/// Returns the number of non-zero elements.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_sparse_get_nnz(
    matrix: *const CsMat<f64>
) -> usize {

    if matrix.is_null() {

        return 0;
    }

    unsafe {

        (*matrix).nnz()
    }
}

/// Sparse matrix-vector multiplication.
///
/// result = matrix * vector
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_sparse_spmv(
    matrix: *const CsMat<f64>,
    vector: *const f64,
    vec_len: usize,
    result: *mut f64,
) -> i32 { unsafe {

    if matrix.is_null()
        || vector.is_null()
        || result.is_null()
    {

        update_last_error(
            "Null pointer passed to \
             rssn_num_sparse_spmv"
                .to_string(),
        );

        return -1;
    }

    let m =  {

        &*matrix
    };

    let v =  {

        std::slice::from_raw_parts(
            vector,
            vec_len,
        )
    };

    match sparse::sp_mat_vec_mul(m, v) {
        | Ok(res) => {

            if res.len() != m.rows() {

                update_last_error(
                    "Internal error: \
                     result length \
                     mismatch"
                        .to_string(),
                );

                return -1;
            }

            ptr::copy_nonoverlapping(
                res.as_ptr(),
                result,
                res.len(),
            );

            0
        },
        | Err(e) => {

            update_last_error(e);

            -1
        },
    }
}}

/// Computes the Frobenius norm.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_sparse_frobenius_norm(
    matrix: *const CsMat<f64>
) -> f64 {

    if matrix.is_null() {

        return 0.0;
    }

    let m = unsafe {

        &*matrix
    };

    sparse::frobenius_norm(m)
}

/// Computes the trace.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_sparse_trace(
    matrix: *const CsMat<f64>,
    out_trace: *mut f64,
) -> i32 { unsafe {

    if matrix.is_null()
        || out_trace.is_null()
    {

        return -1;
    }

    let m =  {

        &*matrix
    };

    match sparse::trace(m) {
        | Ok(t) => {

             {

                *out_trace = t;
            };

            0
        },
        | Err(e) => {

            update_last_error(e);

            -1
        },
    }
}}

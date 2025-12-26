//! Handle-based FFI API for numerical equation solvers.

use crate::numerical::matrix::Matrix;
use crate::numerical::solve::{self, LinearSolution};
use std::ptr;
use std::slice;

/// Solves a linear system Ax = b.
///
/// # Arguments
/// * `matrix_ptr` - Pointer to the coefficient matrix A.
/// * `vector_data` - Pointer to the constant vector b.
/// * `vector_len` - Length of the vector b.
///
/// # Returns
/// A pointer to a `LinearSolution` object, or null on error.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_solve_linear_system_handle(
    matrix_ptr: *const Matrix<f64>,
    vector_data: *const f64,
    vector_len: usize,
) -> *mut LinearSolution {

    if matrix_ptr.is_null() || vector_data.is_null() {

        return ptr::null_mut();
    }

    let matrix = &*matrix_ptr;

    let vector = slice::from_raw_parts(vector_data, vector_len);

    match solve::solve_linear_system(matrix, vector) {
        Ok(solution) => Box::into_raw(Box::new(solution)),
        Err(_) => ptr::null_mut(),
    }
}

/// Frees a `LinearSolution` object.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_solve_free_solution(ptr: *mut LinearSolution) {

    if !ptr.is_null() {

        let _ = Box::from_raw(ptr);
    }
}

/// Checks if the solution is unique.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_solve_is_unique(ptr: *const LinearSolution) -> bool {

    if ptr.is_null() {

        return false;
    }

    matches!(*ptr, LinearSolution::Unique(_))
}

/// Helper to copy vector data.

unsafe fn copy_vec_to_buffer(
    vec: &[f64],
    buffer: *mut f64,
) {

    ptr::copy_nonoverlapping(vec.as_ptr(), buffer, vec.len());
}

/// Gets the data of a unique solution.
///
/// # Arguments
/// * `ptr` - Pointer to the `LinearSolution`.
/// * `buffer` - Buffer to store the solution vector.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_solve_get_unique_solution(
    ptr: *const LinearSolution,
    buffer: *mut f64,
) {

    if ptr.is_null() || buffer.is_null() {

        return;
    }

    if let LinearSolution::Unique(ref sol) = *ptr {

        copy_vec_to_buffer(sol, buffer);
    }
}

/// Gets the length of the unique solution vector.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_solve_get_unique_solution_len(
    ptr: *const LinearSolution
) -> usize {

    if ptr.is_null() {

        return 0;
    }

    if let LinearSolution::Unique(ref sol) = *ptr {

        sol.len()
    } else {

        0
    }
}

// TODO: Add accessors for Parametric solutions if needed.
// For now, Unique and checking for NoSolution (via !is_unique and check types) is a good start.

/// Checks if there is no solution.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_solve_is_no_solution(ptr: *const LinearSolution) -> bool {

    if ptr.is_null() {

        return false;
    }

    matches!(*ptr, LinearSolution::NoSolution)
}

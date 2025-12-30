//! Handle-based FFI API for numerical equation solvers.

use std::ptr;
use std::slice;

use crate::numerical::matrix::Matrix;
use crate::numerical::solve::LinearSolution;
use crate::numerical::solve::{
    self,
};

/// Solves a linear system Ax = b.
///
/// # Arguments
/// * `matrix_ptr` - Pointer to the coefficient matrix A.
/// * `vector_data` - Pointer to the constant vector b.
/// * `vector_len` - Length of the vector b.
///
/// # Returns
/// A pointer to a `LinearSolution` object, or null on error.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_solve_linear_system_handle(
    matrix_ptr: *const Matrix<f64>,
    vector_data: *const f64,
    vector_len: usize,
) -> *mut LinearSolution { unsafe {

    if matrix_ptr.is_null()
        || vector_data.is_null()
    {

        return ptr::null_mut();
    }

    let matrix = &*matrix_ptr;

    let vector = slice::from_raw_parts(
        vector_data,
        vector_len,
    );

    match solve::solve_linear_system(
        matrix,
        vector,
    ) {
        | Ok(solution) => {
            Box::into_raw(Box::new(
                solution,
            ))
        },
        | Err(_) => ptr::null_mut(),
    }
}}

/// Frees a `LinearSolution` object.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_solve_free_solution(
    ptr: *mut LinearSolution
) { unsafe {

    if !ptr.is_null() {

        let _ = Box::from_raw(ptr);
    }
}}

/// Checks if the solution is unique.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub const unsafe extern "C" fn rssn_num_solve_is_unique(
    ptr: *const LinearSolution
) -> bool { unsafe {

    if ptr.is_null() {

        return false;
    }

    matches!(
        *ptr,
        LinearSolution::Unique(_)
    )
}}

/// Helper to copy vector data.

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

const unsafe fn copy_vec_to_buffer(
    vec: &[f64],
    buffer: *mut f64,
) { unsafe {

    ptr::copy_nonoverlapping(
        vec.as_ptr(),
        buffer,
        vec.len(),
    );
}}

/// Gets the data of a unique solution.
///
/// # Arguments
/// * `ptr` - Pointer to the `LinearSolution`.
/// * `buffer` - Buffer to store the solution vector.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_solve_get_unique_solution(
    ptr: *const LinearSolution,
    buffer: *mut f64,
) { unsafe {

    if ptr.is_null() || buffer.is_null()
    {

        return;
    }

    if let LinearSolution::Unique(
        ref sol,
    ) = *ptr
    {

        copy_vec_to_buffer(sol, buffer);
    }
}}

/// Gets the length of the unique solution vector.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub const unsafe extern "C" fn rssn_num_solve_get_unique_solution_len(
    ptr: *const LinearSolution
) -> usize { unsafe {

    if ptr.is_null() {

        return 0;
    }

    if let LinearSolution::Unique(
        ref sol,
    ) = *ptr
    {

        sol.len()
    } else {

        0
    }
}}

// TODO: Add accessors for Parametric solutions if needed.
// For now, Unique and checking for NoSolution (via !is_unique and check types) is a good start.

/// Checks if there is no solution.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub const unsafe extern "C" fn rssn_num_solve_is_no_solution(
    ptr: *const LinearSolution
) -> bool { unsafe {

    if ptr.is_null() {

        return false;
    }

    matches!(
        *ptr,
        LinearSolution::NoSolution
    )
}}

//! Handle-based FFI API for physics CNM functions.

use std::ptr;

use crate::physics::physics_cnm;

/// Solves 1D heat equation using CN and returns a flat array of doubles.
/// The caller is responsible for freeing the memory using `rssn_free_f64_array`.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_physics_cnm_solve_heat_1d(
    initial_condition: *const f64,
    n: usize,
    dx: f64,
    dt: f64,
    d_coeff: f64,
    steps: usize,
    out_size: *mut usize,
) -> *mut f64 { unsafe {

    if initial_condition.is_null() {

        return ptr::null_mut();
    }

    let initial_slice =
        std::slice::from_raw_parts(
            initial_condition,
            n,
        );

    let res = physics_cnm::solve_heat_equation_1d_cn(
        initial_slice,
        dx,
        dt,
        d_coeff,
        steps,
    );

    *out_size = res.len();

    let mut res_boxed =
        res.into_boxed_slice();

    let ptr = res_boxed.as_mut_ptr();

    std::mem::forget(res_boxed);

    ptr
}}

/// Frees a float64 array allocated by the CNM FFI.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_free_f64_cnm_array(
    ptr: *mut f64,
    size: usize,
) { unsafe {

    if !ptr.is_null() {

        let _ = Box::from_raw(std::ptr::slice_from_raw_parts_mut(ptr, size));
    }
}}

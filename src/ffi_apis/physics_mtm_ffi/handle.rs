//! Handle-based FFI API for physics MTM (Multigrid) functions.

use std::ptr;

use crate::physics::physics_mtm;

/// Solves 1D Poisson using Multigrid and returns a flat array of doubles.
///
/// The `out_size` will be set to `n + 2` (including boundaries).
/// The caller is responsible for freeing the memory using `rssn_free_f64_mtm_array`.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_physics_mtm_solve_poisson_1d(
    n_interior: usize,
    f: *const f64,
    num_cycles: usize,
    out_size: *mut usize,
) -> *mut f64 { unsafe {

    if f.is_null() {

        return ptr::null_mut();
    }

    let f_slice =
        std::slice::from_raw_parts(
            f,
            n_interior,
        );

    match physics_mtm::solve_poisson_1d_multigrid(
        n_interior,
        f_slice,
        num_cycles,
    ) {
        | Ok(res) => {

            *out_size = res.len();

            let mut res_boxed = res.into_boxed_slice();

            let ptr = res_boxed.as_mut_ptr();

            std::mem::forget(res_boxed);

            ptr
        },
        | Err(_) => {

            *out_size = 0;

            ptr::null_mut()
        },
    }
}}

/// Solves 2D Poisson using Multigrid and returns a flat array of doubles.
/// The `out_size` will be set to `n * n`.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_physics_mtm_solve_poisson_2d(
    n: usize,
    f: *const f64,
    num_cycles: usize,
    out_size: *mut usize,
) -> *mut f64 { unsafe {

    if f.is_null() {

        return ptr::null_mut();
    }

    let f_slice =
        std::slice::from_raw_parts(
            f,
            n * n,
        );

    match physics_mtm::solve_poisson_2d_multigrid(
        n,
        f_slice,
        num_cycles,
    ) {
        | Ok(res) => {

            *out_size = res.len();

            let mut res_boxed = res.into_boxed_slice();

            let ptr = res_boxed.as_mut_ptr();

            std::mem::forget(res_boxed);

            ptr
        },
        | Err(_) => {

            *out_size = 0;

            ptr::null_mut()
        },
    }
}}

/// Frees a float64 array allocated by the MTM FFI.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_free_f64_mtm_array(
    ptr: *mut f64,
    size: usize,
) { unsafe {

    if !ptr.is_null() {

        let _ = Box::from_raw(std::ptr::slice_from_raw_parts_mut(ptr, size));
    }
}}

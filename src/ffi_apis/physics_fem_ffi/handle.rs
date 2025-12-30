//! Handle-based FFI API for physics FEM functions.

use std::ptr;

use crate::physics::physics_fem;

/// Solves 1D Poisson using FEM and returns a flat array of doubles.
/// The caller is responsible for freeing the memory using `rssn_free_f64_array`.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_physics_fem_solve_poisson_1d(
    n_elements: usize,
    domain_length: f64,
    out_size: *mut usize,
) -> *mut f64 {

    match physics_fem::solve_poisson_1d(
        n_elements,
        domain_length,
        |_| 2.0,
    ) {
        | Ok(res) => {

            unsafe {

                *out_size = res.len();
            }

            let mut res =
                res.into_boxed_slice();

            let ptr = res.as_mut_ptr();

            std::mem::forget(res);

            ptr
        },
        | Err(_) => {

            unsafe {

                *out_size = 0;
            }

            ptr::null_mut()
        },
    }
}

/// Frees a float64 array allocated by the FEM FFI.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_free_f64_array(
    ptr: *mut f64,
    size: usize,
) { unsafe {

    if !ptr.is_null() {

        let _ = Box::from_raw(std::ptr::slice_from_raw_parts_mut(ptr, size));
    }
}}

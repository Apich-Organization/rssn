//! Handle-based FFI API for physics FEM functions.

use crate::physics::physics_fem;
use std::ptr;

/// Solves 1D Poisson using FEM and returns a flat array of doubles.
/// The caller is responsible for freeing the memory using rssn_free_f64_array.
#[no_mangle]
pub extern "C" fn rssn_physics_fem_solve_poisson_1d(
    n_elements: usize,
    domain_length: f64,
    out_size: *mut usize,
) -> *mut f64 {
    match physics_fem::solve_poisson_1d(n_elements, domain_length, |_| 2.0) {
        Ok(res) => {
            unsafe { *out_size = res.len(); }
            let mut res = res.into_boxed_slice();
            let ptr = res.as_mut_ptr();
            std::mem::forget(res);
            ptr
        }
        Err(_) => {
            unsafe { *out_size = 0; }
            ptr::null_mut()
        }
    }
}

/// Frees a float64 array allocated by the FEM FFI.
#[no_mangle]
pub unsafe extern "C" fn rssn_free_f64_array(ptr: *mut f64, size: usize) {
    if !ptr.is_null() {
        let _ = Box::from_raw(std::slice::from_raw_parts_mut(ptr, size));
    }
}

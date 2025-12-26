//! Handle-based FFI API for physics CNM functions.

use std::ptr;

use crate::physics::physics_cnm;

/// Solves 1D heat equation using CN and returns a flat array of doubles.
/// The caller is responsible for freeing the memory using rssn_free_f64_array.
#[no_mangle]

pub unsafe extern "C" fn rssn_physics_cnm_solve_heat_1d(
    initial_condition : *const f64,
    n : usize,
    dx : f64,
    dt : f64,
    d_coeff : f64,
    steps : usize,
    out_size : *mut usize,
) -> *mut f64 {

    if initial_condition.is_null() {

        return ptr::null_mut();
    }

    let initial_slice = std::slice::from_raw_parts(initial_condition, n);

    let res = physics_cnm::solve_heat_equation_1d_cn(
        initial_slice,
        dx,
        dt,
        d_coeff,
        steps,
    );

    *out_size = res.len();

    let mut res_boxed = res.into_boxed_slice();

    let ptr = res_boxed.as_mut_ptr();

    std::mem::forget(res_boxed);

    ptr
}

/// Frees a float64 array allocated by the CNM FFI.
#[no_mangle]

pub unsafe extern "C" fn rssn_free_f64_cnm_array(
    ptr : *mut f64,
    size : usize,
) {

    if !ptr.is_null() {

        let _ = Box::from_raw(std::slice::from_raw_parts_mut(ptr, size));
    }
}

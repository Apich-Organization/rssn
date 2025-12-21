//! Handle-based FFI API for numerical coordinate transformations.

use crate::ffi_apis::ffi_api::update_last_error;
use crate::numerical::coordinates as nc;
use crate::symbolic::coordinates::CoordinateSystem;
use std::os::raw::c_double;
use std::ptr;

/// Transforms a point from one coordinate system to another.
///
/// # Arguments
/// * `point_ptr` - Pointer to an array of doubles.
/// * `point_len` - Number of elements in the point.
/// * `from` - Source coordinate system.
/// * `to` - Target coordinate system.
/// * `out_len` - Pointer to store the number of elements in the resulting vector.
///
/// # Returns
/// A pointer to the transformed point (array of doubles), or null on error.
#[no_mangle]
pub unsafe extern "C" fn rssn_num_coord_transform_point(
    point_ptr: *const c_double,
    point_len: usize,
    from: CoordinateSystem,
    to: CoordinateSystem,
    out_len: *mut usize,
) -> *mut c_double {
    if point_ptr.is_null() || out_len.is_null() {
        return ptr::null_mut();
    }
    let point = unsafe { std::slice::from_raw_parts(point_ptr, point_len) };

    match nc::transform_point(point, from, to) {
        Ok(res) => {
            unsafe {
                *out_len = res.len();
            }
            let mut res_vec = res;
            let ptr = res_vec.as_mut_ptr();
            std::mem::forget(res_vec);
            ptr
        }
        Err(e) => {
            unsafe {
                update_last_error(e);
            }
            ptr::null_mut()
        }
    }
}

/// Transforms a point from one coordinate system to another (pure numerical).
#[no_mangle]
pub unsafe extern "C" fn rssn_num_coord_transform_point_pure(
    point_ptr: *const c_double,
    point_len: usize,
    from: CoordinateSystem,
    to: CoordinateSystem,
    out_len: *mut usize,
) -> *mut c_double {
    if point_ptr.is_null() || out_len.is_null() {
        return ptr::null_mut();
    }
    let point = unsafe { std::slice::from_raw_parts(point_ptr, point_len) };

    match nc::transform_point_pure(point, from, to) {
        Ok(res) => {
            unsafe {
                *out_len = res.len();
            }
            let mut res_vec = res;
            let ptr = res_vec.as_mut_ptr();
            std::mem::forget(res_vec);
            ptr
        }
        Err(e) => {
            unsafe {
                update_last_error(e);
            }
            ptr::null_mut()
        }
    }
}

/// Computes the numerical Jacobian matrix.
/// Returns a pointer to a flat array of doubles (row-major).
#[no_mangle]
pub unsafe extern "C" fn rssn_num_coord_jacobian(
    from: CoordinateSystem,
    to: CoordinateSystem,
    at_point_ptr: *const c_double,
    point_len: usize,
    out_rows: *mut usize,
    out_cols: *mut usize,
) -> *mut c_double {
    if at_point_ptr.is_null() || out_rows.is_null() || out_cols.is_null() {
        return ptr::null_mut();
    }
    let point = unsafe { std::slice::from_raw_parts(at_point_ptr, point_len) };

    match nc::numerical_jacobian(from, to, point) {
        Ok(matrix) => {
            let rows = matrix.rows();
            let cols = matrix.cols();
            unsafe {
                *out_rows = rows;
                *out_cols = cols;
            }
            let mut data = matrix.into_data();
            let ptr = data.as_mut_ptr();
            std::mem::forget(data);
            ptr
        }
        Err(e) => {
            unsafe {
                update_last_error(e);
            }
            ptr::null_mut()
        }
    }
}

/// Frees a pointer allocated by the coordinate transformation functions.
#[no_mangle]
pub unsafe extern "C" fn rssn_num_coord_free(ptr: *mut c_double, len: usize) {
    if !ptr.is_null() {
        unsafe {
            let _ = Vec::from_raw_parts(ptr, len, len);
        }
    }
}

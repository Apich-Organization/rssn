//! Raw pointer-based FFI API for numerical vector operations.

use crate::ffi_apis::ffi_api::update_last_error;
use crate::numerical::vector;
use std::os::raw::c_char;

/// Creates a new numerical vector from a raw array of doubles.
/// The caller is responsible for freeing the returned pointer using rssn_num_vec_free.
#[no_mangle]
pub unsafe extern "C" fn rssn_num_vec_create(data: *const f64, len: usize) -> *mut Vec<f64> {
    if data.is_null() {
        update_last_error("Null pointer passed to rssn_num_vec_create".to_string());
        return std::ptr::null_mut();
    }
    let v = unsafe { std::slice::from_raw_parts(data, len) }.to_vec();
    Box::into_raw(Box::new(v))
}

/// Frees a numerical vector allocated by the library.
#[no_mangle]
pub unsafe extern "C" fn rssn_num_vec_free(v: *mut Vec<f64>) {
    if !v.is_null() {
        unsafe {
            let _ = Box::from_raw(v);
        }
    }
}

/// Returns the length of a numerical vector.
#[no_mangle]
pub unsafe extern "C" fn rssn_num_vec_len(v: *const Vec<f64>) -> usize {
    if v.is_null() {
        return 0;
    }
    unsafe { (*v).len() }
}

/// Returns a pointer to the underlying data of a numerical vector.
#[no_mangle]
pub unsafe extern "C" fn rssn_num_vec_data(v: *const Vec<f64>) -> *const f64 {
    if v.is_null() {
        return std::ptr::null();
    }
    unsafe { (*v).as_ptr() }
}

/// Computes the sum of two vectors.
#[no_mangle]
pub unsafe extern "C" fn rssn_num_vec_add(
    v1: *const Vec<f64>,
    v2: *const Vec<f64>,
) -> *mut Vec<f64> {
    if v1.is_null() || v2.is_null() {
        return std::ptr::null_mut();
    }
    let res = vector::vec_add(unsafe { &*v1 }, unsafe { &*v2 });
    match res {
        Ok(v) => Box::into_raw(Box::new(v)),
        Err(e) => {
            update_last_error(e);
            std::ptr::null_mut()
        }
    }
}

/// Computes the difference of two vectors.
#[no_mangle]
pub unsafe extern "C" fn rssn_num_vec_sub(
    v1: *const Vec<f64>,
    v2: *const Vec<f64>,
) -> *mut Vec<f64> {
    if v1.is_null() || v2.is_null() {
        return std::ptr::null_mut();
    }
    let res = vector::vec_sub(unsafe { &*v1 }, unsafe { &*v2 });
    match res {
        Ok(v) => Box::into_raw(Box::new(v)),
        Err(e) => {
            update_last_error(e);
            std::ptr::null_mut()
        }
    }
}

/// Multiplies a vector by a scalar.
#[no_mangle]
pub unsafe extern "C" fn rssn_num_vec_scalar_mul(v: *const Vec<f64>, s: f64) -> *mut Vec<f64> {
    if v.is_null() {
        return std::ptr::null_mut();
    }
    let res = vector::scalar_mul(unsafe { &*v }, s);
    Box::into_raw(Box::new(res))
}

/// Computes the dot product of two vectors.
#[no_mangle]
pub unsafe extern "C" fn rssn_num_vec_dot_product(
    v1: *const Vec<f64>,
    v2: *const Vec<f64>,
    result: *mut f64,
) -> i32 {
    if v1.is_null() || v2.is_null() || result.is_null() {
        return -1;
    }
    match vector::dot_product(unsafe { &*v1 }, unsafe { &*v2 }) {
        Ok(val) => {
            unsafe {
                *result = val;
            }
            0
        }
        Err(e) => {
            update_last_error(e);
            -1
        }
    }
}

/// Computes the L2 norm of a vector.
#[no_mangle]
pub unsafe extern "C" fn rssn_num_vec_norm(v: *const Vec<f64>, result: *mut f64) -> i32 {
    if v.is_null() || result.is_null() {
        return -1;
    }
    unsafe {
        *result = vector::norm(&*v);
    }
    0
}

/// Computes the Lp norm of a vector.
#[no_mangle]
pub unsafe extern "C" fn rssn_num_vec_lp_norm(v: *const Vec<f64>, p: f64, result: *mut f64) -> i32 {
    if v.is_null() || result.is_null() {
        return -1;
    }
    unsafe {
        *result = vector::lp_norm(&*v, p);
    }
    0
}

/// Normalizes a vector.
#[no_mangle]
pub unsafe extern "C" fn rssn_num_vec_normalize(v: *const Vec<f64>) -> *mut Vec<f64> {
    if v.is_null() {
        return std::ptr::null_mut();
    }
    match vector::normalize(unsafe { &*v }) {
        Ok(res) => Box::into_raw(Box::new(res)),
        Err(e) => {
            update_last_error(e);
            std::ptr::null_mut()
        }
    }
}

/// Computes the cross product of two 3D vectors.
#[no_mangle]
pub unsafe extern "C" fn rssn_num_vec_cross_product(
    v1: *const Vec<f64>,
    v2: *const Vec<f64>,
) -> *mut Vec<f64> {
    if v1.is_null() || v2.is_null() {
        return std::ptr::null_mut();
    }
    match vector::cross_product(unsafe { &*v1 }, unsafe { &*v2 }) {
        Ok(res) => Box::into_raw(Box::new(res)),
        Err(e) => {
            update_last_error(e);
            std::ptr::null_mut()
        }
    }
}

/// Computes the angle between two vectors.
#[no_mangle]
pub unsafe extern "C" fn rssn_num_vec_angle(
    v1: *const Vec<f64>,
    v2: *const Vec<f64>,
    result: *mut f64,
) -> i32 {
    if v1.is_null() || v2.is_null() || result.is_null() {
        return -1;
    }
    match vector::angle(unsafe { &*v1 }, unsafe { &*v2 }) {
        Ok(val) => {
            unsafe {
                *result = val;
            }
            0
        }
        Err(e) => {
            update_last_error(e);
            -1
        }
    }
}

/// Projects v1 onto v2.
#[no_mangle]
pub unsafe extern "C" fn rssn_num_vec_project(
    v1: *const Vec<f64>,
    v2: *const Vec<f64>,
) -> *mut Vec<f64> {
    if v1.is_null() || v2.is_null() {
        return std::ptr::null_mut();
    }
    match vector::project(unsafe { &*v1 }, unsafe { &*v2 }) {
        Ok(res) => Box::into_raw(Box::new(res)),
        Err(e) => {
            update_last_error(e);
            std::ptr::null_mut()
        }
    }
}

/// Reflects v about n.
#[no_mangle]
pub unsafe extern "C" fn rssn_num_vec_reflect(
    v: *const Vec<f64>,
    n: *const Vec<f64>,
) -> *mut Vec<f64> {
    if v.is_null() || n.is_null() {
        return std::ptr::null_mut();
    }
    match vector::reflect(unsafe { &*v }, unsafe { &*n }) {
        Ok(res) => Box::into_raw(Box::new(res)),
        Err(e) => {
            update_last_error(e);
            std::ptr::null_mut()
        }
    }
}
// Correction for reflect: vector::reflect(&*v, &*n)

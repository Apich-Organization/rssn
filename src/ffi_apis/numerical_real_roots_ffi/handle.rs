//! Handle-based FFI API for numerical real root finding.

use crate::numerical::polynomial::Polynomial;
use crate::numerical::real_roots;
use std::ptr;
use std::slice;

/// Findings roots of a polynomial from coefficients.
///
/// # Arguments
/// * `coeffs_ptr` - Pointer to the coefficients array (f64).
/// * `len` - Number of coefficients.
/// * `tolerance` - The tolerance for root finding.
///
/// # Returns
/// A pointer to a `Vec<f64>` containing the sorted real roots, or null on error.
#[no_mangle]

pub unsafe extern "C" fn rssn_real_roots_find_roots(
    coeffs_ptr: *const f64,
    len: usize,
    tolerance: f64,
) -> *mut Vec<f64> {

    if coeffs_ptr.is_null() || len == 0 {

        return ptr::null_mut();
    }

    let coeffs = slice::from_raw_parts(coeffs_ptr, len).to_vec();

    let poly = Polynomial::new(coeffs);

    match real_roots::find_roots(&poly, tolerance) {
        Ok(roots) => Box::into_raw(Box::new(roots)),
        Err(_) => ptr::null_mut(),
    }
}

/// Frees a roots vector.
#[no_mangle]

pub unsafe extern "C" fn rssn_real_roots_free_vec(ptr: *mut Vec<f64>) {

    if !ptr.is_null() {

        let _ = Box::from_raw(ptr);
    }
}

/// Gets the length of the roots vector.
#[no_mangle]

pub unsafe extern "C" fn rssn_real_roots_get_vec_len(ptr: *const Vec<f64>) -> usize {

    if ptr.is_null() {

        return 0;
    }

    (*ptr).len()
}

/// Gets the data of the roots vector.
#[no_mangle]

pub unsafe extern "C" fn rssn_real_roots_get_vec_data(
    ptr: *const Vec<f64>,
    buffer: *mut f64,
) {

    if ptr.is_null() || buffer.is_null() {

        return;
    }

    let vec = &*ptr;

    ptr::copy_nonoverlapping(
        vec.as_ptr(),
        buffer,
        vec.len(),
    );
}

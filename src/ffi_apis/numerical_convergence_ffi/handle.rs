//! Handle-based FFI API for numerical convergence operations.

use std::ptr;
use std::slice;

use crate::numerical::convergence;

/// Applies Aitken's acceleration to the input sequence.
///
/// # Arguments
/// * `data` - Pointer to the input sequence array.
/// * `len` - Length of the input sequence.
///
/// # Returns
/// A pointer to a new `Vec<f64>` containing the accelerated sequence, or null on error.
#[no_mangle]

pub unsafe extern "C" fn rssn_convergence_aitken(
    data: *const f64,
    len: usize,
) -> *mut Vec<f64> {

    if data.is_null() {

        return ptr::null_mut();
    }

    let slice = slice::from_raw_parts(
        data, len,
    );

    let res = convergence::aitken_acceleration(slice);

    Box::into_raw(Box::new(res))
}

/// Applies Richardson extrapolation to the input sequence.
///
/// # Arguments
/// * `data` - Pointer to the input sequence array (assumed approximations with halving steps).
/// * `len` - Length of the input sequence.
///
/// # Returns
/// A pointer to a new `Vec<f64>` containing the extrapolated sequence, or null on error.
#[no_mangle]

pub unsafe extern "C" fn rssn_convergence_richardson(
    data: *const f64,
    len: usize,
) -> *mut Vec<f64> {

    if data.is_null() {

        return ptr::null_mut();
    }

    let slice = slice::from_raw_parts(
        data, len,
    );

    let res = convergence::richardson_extrapolation(slice);

    Box::into_raw(Box::new(res))
}

/// Applies Wynn's epsilon algorithm to the input sequence.
///
/// # Arguments
/// * `data` - Pointer to the input sequence array.
/// * `len` - Length of the input sequence.
///
/// # Returns
/// A pointer to a new `Vec<f64>` containing the accelerated sequence, or null on error.
#[no_mangle]

pub unsafe extern "C" fn rssn_convergence_wynn(
    data: *const f64,
    len: usize,
) -> *mut Vec<f64> {

    if data.is_null() {

        return ptr::null_mut();
    }

    let slice = slice::from_raw_parts(
        data, len,
    );

    let res = convergence::wynn_epsilon(
        slice,
    );

    Box::into_raw(Box::new(res))
}

/// Frees a generic `Vec<f64>` pointer created by convergence functions.
#[no_mangle]

pub unsafe extern "C" fn rssn_convergence_free_vec(
    vec: *mut Vec<f64>
) {

    if !vec.is_null() {

        unsafe {

            let _ = Box::from_raw(vec);
        }
    }
}

/// Returns the length of the vector.
#[no_mangle]

pub unsafe extern "C" fn rssn_convergence_get_vec_len(
    vec: *const Vec<f64>
) -> usize {

    if vec.is_null() {

        return 0;
    }

    unsafe {

        (*vec).len()
    }
}

/// Copies the vector data into a provided buffer.
/// buffer must have size at least `len * sizeof(f64)`.
#[no_mangle]

pub unsafe extern "C" fn rssn_convergence_get_vec_data(
    vec: *const Vec<f64>,
    buffer: *mut f64,
) {

    if vec.is_null() || buffer.is_null()
    {

        return;
    }

    let v = unsafe {

        &*vec
    };

    unsafe {

        ptr::copy_nonoverlapping(
            v.as_ptr(),
            buffer,
            v.len(),
        );
    }
}

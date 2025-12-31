//! Handle-based FFI API for numerical real root finding.

use std::ptr;
use std::slice;

use crate::numerical::polynomial::Polynomial;
use crate::numerical::real_roots;

/// Findings roots of a polynomial from coefficients.
///
/// # Arguments
/// * `coeffs_ptr` - Pointer to the coefficients array (f64).
/// * `len` - Number of coefficients.
/// * `tolerance` - The tolerance for root finding.
///
/// # Returns
/// A pointer to a `Vec<f64>` containing the sorted real roots, or null on error.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_real_roots_find_roots(
    coeffs_ptr: *const f64,
    len: usize,
    tolerance: f64,
) -> *mut Vec<f64> {

    unsafe {

        if coeffs_ptr.is_null()
            || len == 0
        {

            return ptr::null_mut();
        }

        let coeffs =
            slice::from_raw_parts(
                coeffs_ptr,
                len,
            )
            .to_vec();

        let poly =
            Polynomial::new(coeffs);

        match real_roots::find_roots(
            &poly,
            tolerance,
        ) {
            | Ok(roots) => {
                Box::into_raw(Box::new(
                    roots,
                ))
            },
            | Err(_) => ptr::null_mut(),
        }
    }
}

/// Frees a roots vector.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_real_roots_free_vec(
    ptr: *mut Vec<f64>
) {

    unsafe {

        if !ptr.is_null() {

            let _ = Box::from_raw(ptr);
        }
    }
}

/// Gets the length of the roots vector.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub const unsafe extern "C" fn rssn_real_roots_get_vec_len(
    ptr: *const Vec<f64>
) -> usize {

    unsafe {

        if ptr.is_null() {

            return 0;
        }

        (*ptr).len()
    }
}

/// Gets the data of the roots vector.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub const unsafe extern "C" fn rssn_real_roots_get_vec_data(
    ptr: *const Vec<f64>,
    buffer: *mut f64,
) {

    unsafe {

        if ptr.is_null()
            || buffer.is_null()
        {

            return;
        }

        let vec = &*ptr;

        ptr::copy_nonoverlapping(
            vec.as_ptr(),
            buffer,
            vec.len(),
        );
    }
}

//! Handle-based FFI API for numerical polynomial operations.

use std::ptr;

use crate::ffi_apis::ffi_api::update_last_error;
use crate::numerical::polynomial::Polynomial;

/// Creates a new polynomial from coefficients.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_poly_create(
    coeffs: *const f64,
    len: usize,
) -> *mut Polynomial { unsafe {

    if coeffs.is_null() {

        update_last_error(
            "Null pointer passed to \
             rssn_num_poly_create"
                .to_string(),
        );

        return ptr::null_mut();
    }

    let c = unsafe {

        std::slice::from_raw_parts(
            coeffs,
            len,
        )
    };

    Box::into_raw(Box::new(
        Polynomial {
            coeffs: c.to_vec(),
        },
    ))
}}

/// Frees a polynomial object.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_poly_free(
    poly: *mut Polynomial
) {

    if !poly.is_null() {

        unsafe {

            let _ = Box::from_raw(poly);
        }
    }
}

/// Evaluates a polynomial at x.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_poly_eval(
    poly: *const Polynomial,
    x: f64,
) -> f64 {

    if poly.is_null() {

        return 0.0;
    }

    unsafe {

        (*poly).eval(x)
    }
}

/// Returns the degree of a polynomial.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub const unsafe extern "C" fn rssn_num_poly_degree(
    poly: *const Polynomial
) -> usize {

    if poly.is_null() {

        return 0;
    }

    unsafe {

        (*poly).degree()
    }
}

/// Computes the derivative.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_poly_derivative(
    poly: *const Polynomial
) -> *mut Polynomial {

    if poly.is_null() {

        return ptr::null_mut();
    }

    let res = unsafe {

        (*poly).derivative()
    };

    Box::into_raw(Box::new(res))
}

/// Computes the integral.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_poly_integral(
    poly: *const Polynomial
) -> *mut Polynomial {

    if poly.is_null() {

        return ptr::null_mut();
    }

    let res = unsafe {

        (*poly).integral()
    };

    Box::into_raw(Box::new(res))
}

/// Adds two polynomials.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_poly_add(
    a: *const Polynomial,
    b: *const Polynomial,
) -> *mut Polynomial {

    if a.is_null() || b.is_null() {

        return ptr::null_mut();
    }

    let res = unsafe {

        (*a).clone() + (*b).clone()
    };

    Box::into_raw(Box::new(res))
}

/// Subtracts two polynomials.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_poly_sub(
    a: *const Polynomial,
    b: *const Polynomial,
) -> *mut Polynomial {

    if a.is_null() || b.is_null() {

        return ptr::null_mut();
    }

    let res = unsafe {

        (*a).clone() - (*b).clone()
    };

    Box::into_raw(Box::new(res))
}

/// Multiplies two polynomials.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_poly_mul(
    a: *const Polynomial,
    b: *const Polynomial,
) -> *mut Polynomial {

    if a.is_null() || b.is_null() {

        return ptr::null_mut();
    }

    let res = unsafe {

        (*a).clone() * (*b).clone()
    };

    Box::into_raw(Box::new(res))
}

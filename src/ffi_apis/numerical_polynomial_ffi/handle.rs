//! Handle-based FFI API for numerical polynomial operations.

use std::ptr;

use crate::ffi_apis::ffi_api::update_last_error;
use crate::numerical::polynomial::Polynomial;

/// Creates a new polynomial from coefficients.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_poly_create(
    coeffs: *const f64,
    len: usize,
) -> *mut Polynomial {

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
}

/// Frees a polynomial object.
#[no_mangle]

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
#[no_mangle]

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
#[no_mangle]

pub unsafe extern "C" fn rssn_num_poly_degree(
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
#[no_mangle]

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
#[no_mangle]

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
#[no_mangle]

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
#[no_mangle]

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
#[no_mangle]

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

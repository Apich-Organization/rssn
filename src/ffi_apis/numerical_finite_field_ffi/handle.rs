//! Handle-based FFI API for numerical finite field arithmetic.

use crate::ffi_apis::ffi_api::update_last_error;
use crate::numerical::finite_field::{self, PrimeFieldElement};
use std::ptr;

/// Creates a new PrimeFieldElement.
#[no_mangle]

pub extern "C" fn rssn_num_ff_pfe_new(value: u64, modulus: u64) -> *mut PrimeFieldElement {

    Box::into_raw(Box::new(PrimeFieldElement::new(value, modulus)))
}

/// Frees a PrimeFieldElement.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_ff_pfe_free(pfe: *mut PrimeFieldElement) {

    if !pfe.is_null() {

        unsafe {

            let _ = Box::from_raw(pfe);
        }
    }
}

/// Computes the inverse of a PrimeFieldElement.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_ff_pfe_inverse(
    pfe: *const PrimeFieldElement,
) -> *mut PrimeFieldElement {

    if pfe.is_null() {

        return ptr::null_mut();
    }

    match unsafe {

        (*pfe).inverse()
    } {
        Some(inv) => Box::into_raw(Box::new(inv)),
        None => {

            update_last_error("Element is not invertible".to_string());

            ptr::null_mut()
        }
    }
}

/// Computes (pfe^exp) mod modulus.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_ff_pfe_pow(
    pfe: *const PrimeFieldElement,
    exp: u64,
) -> *mut PrimeFieldElement {

    if pfe.is_null() {

        return ptr::null_mut();
    }

    let res = unsafe {

        (*pfe).pow(exp)
    };

    Box::into_raw(Box::new(res))
}

/// Performs addition in GF(p).
#[no_mangle]

pub unsafe extern "C" fn rssn_num_ff_pfe_add(
    a: *const PrimeFieldElement,
    b: *const PrimeFieldElement,
) -> *mut PrimeFieldElement {

    if a.is_null() || b.is_null() {

        return ptr::null_mut();
    }

    let res = unsafe {

        *a + *b
    };

    Box::into_raw(Box::new(res))
}

/// Performs multiplication in GF(p).
#[no_mangle]

pub unsafe extern "C" fn rssn_num_ff_pfe_mul(
    a: *const PrimeFieldElement,
    b: *const PrimeFieldElement,
) -> *mut PrimeFieldElement {

    if a.is_null() || b.is_null() {

        return ptr::null_mut();
    }

    let res = unsafe {

        *a * *b
    };

    Box::into_raw(Box::new(res))
}

/// GF(2^8) addition.
#[no_mangle]

pub extern "C" fn rssn_num_ff_gf256_add(a: u8, b: u8) -> u8 {

    finite_field::gf256_add(a, b)
}

/// GF(2^8) multiplication.
#[no_mangle]

pub extern "C" fn rssn_num_ff_gf256_mul(a: u8, b: u8) -> u8 {

    finite_field::gf256_mul(a, b)
}

/// GF(2^8) division.
/// Returns 0 and sets error if divisor is 0.
#[no_mangle]

pub extern "C" fn rssn_num_ff_gf256_div(a: u8, b: u8) -> u8 {

    match finite_field::gf256_div(a, b) {
        Ok(res) => res,
        Err(e) => {

            unsafe {

                update_last_error(e);
            }

            0
        }
    }
}

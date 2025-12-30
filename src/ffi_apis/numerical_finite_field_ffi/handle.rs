//! Handle-based FFI API for numerical finite field arithmetic.

use std::ptr;

use crate::ffi_apis::ffi_api::update_last_error;
use crate::numerical::finite_field::PrimeFieldElement;
use crate::numerical::finite_field::{
    self,
};

/// Creates a new `PrimeFieldElement`.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_num_ff_pfe_new(
    value: u64,
    modulus: u64,
) -> *mut PrimeFieldElement {

    Box::into_raw(Box::new(
        PrimeFieldElement::new(
            value,
            modulus,
        ),
    ))
}

/// Frees a `PrimeFieldElement`.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_ff_pfe_free(
    pfe: *mut PrimeFieldElement
) {

    if !pfe.is_null() {

        unsafe {

            let _ = Box::from_raw(pfe);
        }
    }
}

/// Computes the inverse of a `PrimeFieldElement`.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_ff_pfe_inverse(
    pfe: *const PrimeFieldElement
) -> *mut PrimeFieldElement { unsafe {

    if pfe.is_null() {

        return ptr::null_mut();
    }

    match (*pfe).inverse() {
        | Some(inv) => {
            Box::into_raw(Box::new(inv))
        },
        | None => {

            update_last_error(
                "Element is not \
                 invertible"
                    .to_string(),
            );

            ptr::null_mut()
        },
    }
}}

/// Computes (pfe^exp) mod modulus.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

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
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

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
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

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
#[unsafe(no_mangle)]

pub const extern "C" fn rssn_num_ff_gf256_add(
    a: u8,
    b: u8,
) -> u8 {

    finite_field::gf256_add(a, b)
}

/// GF(2^8) multiplication.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_num_ff_gf256_mul(
    a: u8,
    b: u8,
) -> u8 {

    finite_field::gf256_mul(a, b)
}

/// GF(2^8) division.
/// Returns 0 and sets error if divisor is 0.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_num_ff_gf256_div(
    a: u8,
    b: u8,
) -> u8 {

    match finite_field::gf256_div(a, b)
    {
        | Ok(res) => res,
        | Err(e) => {

            unsafe {

                update_last_error(e);
            }

            0
        },
    }
}

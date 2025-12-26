use crate::symbolic::finite_field::*;
use num_bigint::BigInt;

/// Creates a new prime field element (Handle)
/// Returns a boxed pointer to the element
#[no_mangle]

pub extern "C" fn rssn_prime_field_element_new_handle(
    value: *const BigInt,
    modulus: *const BigInt,
) -> *mut PrimeFieldElement {

    let value_ref = unsafe {

        if value.is_null() {

            return std::ptr::null_mut();
        }

        &*value
    };

    let modulus_ref = unsafe {

        if modulus.is_null() {

            return std::ptr::null_mut();
        }

        &*modulus
    };

    let field = PrimeField::new(modulus_ref.clone());

    let element = PrimeFieldElement::new(
        value_ref.clone(),
        field,
    );

    Box::into_raw(Box::new(element))
}

/// Adds two prime field elements (Handle)
#[no_mangle]

pub extern "C" fn rssn_prime_field_element_add_handle(
    a: *const PrimeFieldElement,
    b: *const PrimeFieldElement,
) -> *mut PrimeFieldElement {

    let a_ref = unsafe {

        if a.is_null() {

            return std::ptr::null_mut();
        }

        &*a
    };

    let b_ref = unsafe {

        if b.is_null() {

            return std::ptr::null_mut();
        }

        &*b
    };

    let result = a_ref.clone() + b_ref.clone();

    Box::into_raw(Box::new(result))
}

/// Multiplies two prime field elements (Handle)
#[no_mangle]

pub extern "C" fn rssn_prime_field_element_mul_handle(
    a: *const PrimeFieldElement,
    b: *const PrimeFieldElement,
) -> *mut PrimeFieldElement {

    let a_ref = unsafe {

        if a.is_null() {

            return std::ptr::null_mut();
        }

        &*a
    };

    let b_ref = unsafe {

        if b.is_null() {

            return std::ptr::null_mut();
        }

        &*b
    };

    let result = a_ref.clone() * b_ref.clone();

    Box::into_raw(Box::new(result))
}

/// Computes the inverse of a prime field element (Handle)
#[no_mangle]

pub extern "C" fn rssn_prime_field_element_inverse_handle(
    elem: *const PrimeFieldElement
) -> *mut PrimeFieldElement {

    let elem_ref = unsafe {

        if elem.is_null() {

            return std::ptr::null_mut();
        }

        &*elem
    };

    match elem_ref.inverse() {
        | Some(inv) => Box::into_raw(Box::new(inv)),
        | None => std::ptr::null_mut(),
    }
}

/// Frees a prime field element (Handle)
#[no_mangle]

pub extern "C" fn rssn_prime_field_element_free_handle(elem: *mut PrimeFieldElement) {

    if !elem.is_null() {

        unsafe {

            let _ = Box::from_raw(elem);
        }
    }
}

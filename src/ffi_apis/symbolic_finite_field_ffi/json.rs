use crate::ffi_apis::common::*;
use crate::symbolic::finite_field::*;
use num_bigint::BigInt;
use std::ffi::c_char;

/// Creates a new prime field element (JSON)
#[no_mangle]

pub extern "C" fn rssn_json_prime_field_element_new(
    value_json: *const c_char,
    modulus_json: *const c_char,
) -> *mut c_char {

    let value: Option<BigInt> =
        from_json_string(value_json);

    let modulus: Option<BigInt> =
        from_json_string(modulus_json);

    if let (Some(v), Some(m)) =
        (value, modulus)
    {

        let field = PrimeField::new(m);

        let element =
            PrimeFieldElement::new(
                v, field,
            );

        to_json_string(&element)
    } else {

        std::ptr::null_mut()
    }
}

/// Adds two prime field elements (JSON)
#[no_mangle]

pub extern "C" fn rssn_json_prime_field_element_add(
    a_json: *const c_char,
    b_json: *const c_char,
) -> *mut c_char {

    let a: Option<PrimeFieldElement> =
        from_json_string(a_json);

    let b: Option<PrimeFieldElement> =
        from_json_string(b_json);

    if let (
        Some(a_elem),
        Some(b_elem),
    ) = (a, b)
    {

        let result = a_elem + b_elem;

        to_json_string(&result)
    } else {

        std::ptr::null_mut()
    }
}

/// Subtracts two prime field elements (JSON)
#[no_mangle]

pub extern "C" fn rssn_json_prime_field_element_sub(
    a_json: *const c_char,
    b_json: *const c_char,
) -> *mut c_char {

    let a: Option<PrimeFieldElement> =
        from_json_string(a_json);

    let b: Option<PrimeFieldElement> =
        from_json_string(b_json);

    if let (
        Some(a_elem),
        Some(b_elem),
    ) = (a, b)
    {

        let result = a_elem - b_elem;

        to_json_string(&result)
    } else {

        std::ptr::null_mut()
    }
}

/// Multiplies two prime field elements (JSON)
#[no_mangle]

pub extern "C" fn rssn_json_prime_field_element_mul(
    a_json: *const c_char,
    b_json: *const c_char,
) -> *mut c_char {

    let a: Option<PrimeFieldElement> =
        from_json_string(a_json);

    let b: Option<PrimeFieldElement> =
        from_json_string(b_json);

    if let (
        Some(a_elem),
        Some(b_elem),
    ) = (a, b)
    {

        let result = a_elem * b_elem;

        to_json_string(&result)
    } else {

        std::ptr::null_mut()
    }
}

/// Divides two prime field elements (JSON)
#[no_mangle]

pub extern "C" fn rssn_json_prime_field_element_div(
    a_json: *const c_char,
    b_json: *const c_char,
) -> *mut c_char {

    let a: Option<PrimeFieldElement> =
        from_json_string(a_json);

    let b: Option<PrimeFieldElement> =
        from_json_string(b_json);

    if let (
        Some(a_elem),
        Some(b_elem),
    ) = (a, b)
    {

        let result = a_elem / b_elem;

        to_json_string(&result)
    } else {

        std::ptr::null_mut()
    }
}

/// Computes the inverse of a prime field element (JSON)
#[no_mangle]

pub extern "C" fn rssn_json_prime_field_element_inverse(
    elem_json: *const c_char
) -> *mut c_char {

    let elem: Option<
        PrimeFieldElement,
    > = from_json_string(elem_json);

    if let Some(e) = elem {

        match e.inverse() {
            | Some(inv) => {
                to_json_string(&inv)
            },
            | None => {
                std::ptr::null_mut()
            },
        }
    } else {

        std::ptr::null_mut()
    }
}

/// Creates a new finite field polynomial (JSON)
#[no_mangle]

pub extern "C" fn rssn_json_finite_field_polynomial_new(
    coeffs_json: *const c_char,
    modulus_json: *const c_char,
) -> *mut c_char {

    let coeffs: Option<
        Vec<PrimeFieldElement>,
    > = from_json_string(coeffs_json);

    let modulus: Option<BigInt> =
        from_json_string(modulus_json);

    if let (Some(c), Some(m)) =
        (coeffs, modulus)
    {

        let field = PrimeField::new(m);

        let poly =
            FiniteFieldPolynomial::new(
                c, field,
            );

        to_json_string(&poly)
    } else {

        std::ptr::null_mut()
    }
}

/// Gets the degree of a finite field polynomial (JSON)
#[no_mangle]

pub extern "C" fn rssn_json_finite_field_polynomial_degree(
    poly_json: *const c_char
) -> i64 {

    let poly: Option<
        FiniteFieldPolynomial,
    > = from_json_string(poly_json);

    if let Some(p) = poly {

        p.degree() as i64
    } else {

        -1
    }
}

/// Performs polynomial long division (JSON)
#[no_mangle]

pub extern "C" fn rssn_json_finite_field_polynomial_long_division(
    dividend_json: *const c_char,
    divisor_json: *const c_char,
) -> *mut c_char {

    let dividend: Option<
        FiniteFieldPolynomial,
    > = from_json_string(dividend_json);

    let divisor: Option<
        FiniteFieldPolynomial,
    > = from_json_string(divisor_json);

    if let (
        Some(div),
        Some(divisor_poly),
    ) = (dividend, divisor)
    {

        match div.long_division(
            &divisor_poly,
        ) {
            | Ok((
                quotient,
                remainder,
            )) => {

                let result = serde_json::json!({
                    "quotient": quotient,
                    "remainder": remainder
                });

                to_json_string(&result)
            },
            | Err(_) => {
                std::ptr::null_mut()
            },
        }
    } else {

        std::ptr::null_mut()
    }
}

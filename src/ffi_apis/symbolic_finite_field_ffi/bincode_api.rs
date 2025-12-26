use num_bigint::BigInt;

use crate::ffi_apis::common::*;
use crate::symbolic::finite_field::*;

/// Creates a new prime field element (Bincode)
#[no_mangle]

pub extern "C" fn rssn_bincode_prime_field_element_new(
    value_buf : BincodeBuffer,
    modulus_buf : BincodeBuffer,
) -> BincodeBuffer {

    let value : Option<BigInt> =
        from_bincode_buffer(&value_buf);

    let modulus : Option<BigInt> =
        from_bincode_buffer(
            &modulus_buf,
        );

    if let (Some(v), Some(m)) =
        (value, modulus)
    {

        let field = PrimeField::new(m);

        let element =
            PrimeFieldElement::new(
                v, field,
            );

        to_bincode_buffer(&element)
    } else {

        BincodeBuffer::empty()
    }
}

/// Adds two prime field elements (Bincode)
#[no_mangle]

pub extern "C" fn rssn_bincode_prime_field_element_add(
    a_buf : BincodeBuffer,
    b_buf : BincodeBuffer,
) -> BincodeBuffer {

    let a : Option<PrimeFieldElement> =
        from_bincode_buffer(&a_buf);

    let b : Option<PrimeFieldElement> =
        from_bincode_buffer(&b_buf);

    if let (
        Some(a_elem),
        Some(b_elem),
    ) = (a, b)
    {

        let result = a_elem + b_elem;

        to_bincode_buffer(&result)
    } else {

        BincodeBuffer::empty()
    }
}

/// Subtracts two prime field elements (Bincode)
#[no_mangle]

pub extern "C" fn rssn_bincode_prime_field_element_sub(
    a_buf : BincodeBuffer,
    b_buf : BincodeBuffer,
) -> BincodeBuffer {

    let a : Option<PrimeFieldElement> =
        from_bincode_buffer(&a_buf);

    let b : Option<PrimeFieldElement> =
        from_bincode_buffer(&b_buf);

    if let (
        Some(a_elem),
        Some(b_elem),
    ) = (a, b)
    {

        let result = a_elem - b_elem;

        to_bincode_buffer(&result)
    } else {

        BincodeBuffer::empty()
    }
}

/// Multiplies two prime field elements (Bincode)
#[no_mangle]

pub extern "C" fn rssn_bincode_prime_field_element_mul(
    a_buf : BincodeBuffer,
    b_buf : BincodeBuffer,
) -> BincodeBuffer {

    let a : Option<PrimeFieldElement> =
        from_bincode_buffer(&a_buf);

    let b : Option<PrimeFieldElement> =
        from_bincode_buffer(&b_buf);

    if let (
        Some(a_elem),
        Some(b_elem),
    ) = (a, b)
    {

        let result = a_elem * b_elem;

        to_bincode_buffer(&result)
    } else {

        BincodeBuffer::empty()
    }
}

/// Divides two prime field elements (Bincode)
#[no_mangle]

pub extern "C" fn rssn_bincode_prime_field_element_div(
    a_buf : BincodeBuffer,
    b_buf : BincodeBuffer,
) -> BincodeBuffer {

    let a : Option<PrimeFieldElement> =
        from_bincode_buffer(&a_buf);

    let b : Option<PrimeFieldElement> =
        from_bincode_buffer(&b_buf);

    if let (
        Some(a_elem),
        Some(b_elem),
    ) = (a, b)
    {

        let result = a_elem / b_elem;

        to_bincode_buffer(&result)
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes the inverse of a prime field element (Bincode)
#[no_mangle]

pub extern "C" fn rssn_bincode_prime_field_element_inverse(
    elem_buf : BincodeBuffer
) -> BincodeBuffer {

    let elem : Option<
        PrimeFieldElement,
    > = from_bincode_buffer(&elem_buf);

    if let Some(e) = elem {

        match e.inverse() {
            | Some(inv) => {
                to_bincode_buffer(&inv)
            },
            | None => {
                BincodeBuffer::empty()
            },
        }
    } else {

        BincodeBuffer::empty()
    }
}

/// Creates a new finite field polynomial (Bincode)
#[no_mangle]

pub extern "C" fn rssn_bincode_finite_field_polynomial_new(
    coeffs_buf : BincodeBuffer,
    modulus_buf : BincodeBuffer,
) -> BincodeBuffer {

    let coeffs : Option<
        Vec<PrimeFieldElement>,
    > = from_bincode_buffer(
        &coeffs_buf,
    );

    let modulus : Option<BigInt> =
        from_bincode_buffer(
            &modulus_buf,
        );

    if let (Some(c), Some(m)) =
        (coeffs, modulus)
    {

        let field = PrimeField::new(m);

        let poly =
            FiniteFieldPolynomial::new(
                c, field,
            );

        to_bincode_buffer(&poly)
    } else {

        BincodeBuffer::empty()
    }
}

/// Gets the degree of a finite field polynomial (Bincode)
#[no_mangle]

pub extern "C" fn rssn_bincode_finite_field_polynomial_degree(
    poly_buf : BincodeBuffer
) -> i64 {

    let poly : Option<
        FiniteFieldPolynomial,
    > = from_bincode_buffer(&poly_buf);

    if let Some(p) = poly {

        p.degree() as i64
    } else {

        -1
    }
}

/// Performs polynomial long division (Bincode)
#[no_mangle]

pub extern "C" fn rssn_bincode_finite_field_polynomial_long_division(
    dividend_buf : BincodeBuffer,
    divisor_buf : BincodeBuffer,
) -> BincodeBuffer {

    let dividend : Option<
        FiniteFieldPolynomial,
    > = from_bincode_buffer(
        &dividend_buf,
    );

    let divisor : Option<
        FiniteFieldPolynomial,
    > = from_bincode_buffer(
        &divisor_buf,
    );

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

                #[derive(
                    serde::Serialize,
                )]

                struct DivisionResult {
                    quotient: FiniteFieldPolynomial,
                    remainder: FiniteFieldPolynomial,
                }

                let result =
                    DivisionResult {
                        quotient,
                        remainder,
                    };

                to_bincode_buffer(
                    &result,
                )
            },
            | Err(_) => {
                BincodeBuffer::empty()
            },
        }
    } else {

        BincodeBuffer::empty()
    }
}

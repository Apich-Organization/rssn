use crate::ffi_apis::common::*;
use crate::symbolic::finite_field::FiniteFieldPolynomial;
use crate::symbolic::poly_factorization::*;

/// Factors a polynomial over a finite field (Bincode)
#[no_mangle]

pub extern "C" fn rssn_bincode_factor_gf(
    poly_buf: BincodeBuffer
) -> BincodeBuffer {

    let poly: Option<
        FiniteFieldPolynomial,
    > = from_bincode_buffer(&poly_buf);

    if let Some(p) = poly {

        match factor_gf(&p) {
            | Ok(factors) => {
                to_bincode_buffer(
                    &factors,
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

/// Computes square-free factorization (Bincode)
#[no_mangle]

pub extern "C" fn rssn_bincode_square_free_factorization_gf(
    poly_buf: BincodeBuffer
) -> BincodeBuffer {

    let poly: Option<
        FiniteFieldPolynomial,
    > = from_bincode_buffer(&poly_buf);

    if let Some(p) = poly {

        match square_free_factorization_gf(p) {
            | Ok(factors) => to_bincode_buffer(&factors),
            | Err(_) => BincodeBuffer::empty(),
        }
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes polynomial GCD over finite field (Bincode)
#[no_mangle]

pub extern "C" fn rssn_bincode_poly_gcd_gf(
    a_buf: BincodeBuffer,
    b_buf: BincodeBuffer,
) -> BincodeBuffer {

    let a: Option<
        FiniteFieldPolynomial,
    > = from_bincode_buffer(&a_buf);

    let b: Option<
        FiniteFieldPolynomial,
    > = from_bincode_buffer(&b_buf);

    if let (
        Some(poly_a),
        Some(poly_b),
    ) = (a, b)
    {

        match poly_gcd_gf(
            poly_a, poly_b,
        ) {
            | Ok(gcd) => {
                to_bincode_buffer(&gcd)
            },
            | Err(_) => {
                BincodeBuffer::empty()
            },
        }
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes polynomial derivative over finite field (Bincode)
#[no_mangle]

pub extern "C" fn rssn_bincode_poly_derivative_gf(
    poly_buf: BincodeBuffer
) -> BincodeBuffer {

    let poly: Option<
        FiniteFieldPolynomial,
    > = from_bincode_buffer(&poly_buf);

    if let Some(p) = poly {

        let derivative =
            poly_derivative_gf(&p);

        to_bincode_buffer(&derivative)
    } else {

        BincodeBuffer::empty()
    }
}

use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::ffi_apis::common::BincodeBuffer;
use crate::symbolic::core::SparsePolynomial;
use crate::symbolic::grobner::buchberger;
use crate::symbolic::grobner::poly_division_multivariate;
use crate::symbolic::grobner::MonomialOrder;

#[no_mangle]

/// Computes a Gröbner basis using Buchberger's algorithm and returns it via bincode serialization.
///
/// Given a basis of multivariate polynomials and a monomial order, this runs
/// Buchberger's algorithm to produce a Gröbner basis for the ideal they generate.
///
/// # Arguments
///
/// * `basis_buf` - `BincodeBuffer` encoding `Vec<SparsePolynomial>` for the initial basis.
/// * `order_buf` - `BincodeBuffer` encoding the [`MonomialOrder`] to use.
///
/// # Returns
///
/// A `BincodeBuffer` encoding `Vec<SparsePolynomial>` forming a Gröbner basis, or an
/// empty buffer if deserialization fails or the computation encounters an error.
///
/// # Safety
///
/// This function is an FFI entry point; callers must treat the returned buffer as
/// opaque and only pass it to compatible APIs.

pub extern "C" fn rssn_bincode_buchberger(
    basis_buf: BincodeBuffer,
    order_buf: BincodeBuffer,
) -> BincodeBuffer {

    let basis: Option<
        Vec<SparsePolynomial>,
    > = from_bincode_buffer(&basis_buf);

    let order: Option<MonomialOrder> =
        from_bincode_buffer(&order_buf);

    if let (Some(b), Some(o)) =
        (basis, order)
    {

        match buchberger(&b, o) {
            | Ok(result) => {
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

#[no_mangle]

/// Divides a multivariate polynomial by a list of divisors under a given monomial order,
/// returning the quotients and remainder via bincode serialization.
///
/// # Arguments
///
/// * `dividend_buf` - `BincodeBuffer` encoding the dividend `SparsePolynomial`.
/// * `divisors_buf` - `BincodeBuffer` encoding `Vec<SparsePolynomial>` of divisors.
/// * `order_buf` - `BincodeBuffer` encoding the [`MonomialOrder`] to use.
///
/// # Returns
///
/// A `BincodeBuffer` encoding `(Vec<SparsePolynomial>, SparsePolynomial)` containing
/// the quotient polynomials and the remainder, or an empty buffer if deserialization
/// fails or the division fails.
///
/// # Safety
///
/// This function is an FFI entry point; callers must treat the returned buffer as
/// opaque and only pass it to compatible APIs.

pub extern "C" fn rssn_bincode_poly_division_multivariate(
    dividend_buf: BincodeBuffer,
    divisors_buf: BincodeBuffer,
    order_buf: BincodeBuffer,
) -> BincodeBuffer {

    let dividend: Option<
        SparsePolynomial,
    > = from_bincode_buffer(
        &dividend_buf,
    );

    let divisors: Option<
        Vec<SparsePolynomial>,
    > = from_bincode_buffer(
        &divisors_buf,
    );

    let order: Option<MonomialOrder> =
        from_bincode_buffer(&order_buf);

    if let (
        Some(d),
        Some(divs),
        Some(o),
    ) = (
        dividend,
        divisors,
        order,
    ) {

        match poly_division_multivariate(
            &d, &divs, o,
        ) {
            | Ok(result) => {
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

use std::ffi::c_char;

use crate::ffi_apis::common::*;
use crate::symbolic::core::SparsePolynomial;
use crate::symbolic::grobner::buchberger;
use crate::symbolic::grobner::poly_division_multivariate;
use crate::symbolic::grobner::MonomialOrder;

#[no_mangle]

/// Computes a Gröbner basis using Buchberger's algorithm and returns it as JSON-encoded polynomials.
///
/// # Arguments
///
/// * `basis_json` - C string pointer with JSON-encoded `Vec<SparsePolynomial>` for the initial basis.
/// * `order_json` - C string pointer with JSON-encoded [`MonomialOrder`].
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `Vec<SparsePolynomial>` representing a
/// Gröbner basis, or null if deserialization fails or the computation encounters an error.
///
/// # Safety
///
/// This function is unsafe because it dereferences raw C string pointers and returns
/// ownership of a heap-allocated C string that must be freed by the caller.

pub extern "C" fn rssn_json_buchberger(
    basis_json: *const c_char,
    order_json: *const c_char,
) -> *mut c_char {

    let basis: Option<
        Vec<SparsePolynomial>,
    > = from_json_string(basis_json);

    let order: Option<MonomialOrder> =
        from_json_string(order_json);

    if let (Some(b), Some(o)) =
        (basis, order)
    {

        match buchberger(&b, o) {
            | Ok(result) => {
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

#[no_mangle]

/// Divides a multivariate polynomial by a list of divisors under a given monomial order
/// and returns quotients and remainder as JSON-encoded polynomials.
///
/// # Arguments
///
/// * `dividend_json` - C string pointer with JSON-encoded dividend `SparsePolynomial`.
/// * `divisors_json` - C string pointer with JSON-encoded `Vec<SparsePolynomial>` of divisors.
/// * `order_json` - C string pointer with JSON-encoded [`MonomialOrder`].
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `(Vec<SparsePolynomial>, SparsePolynomial)`
/// with quotient polynomials and remainder, or null if deserialization fails or the
/// computation encounters an error.
///
/// # Safety
///
/// This function is unsafe because it dereferences raw C string pointers and returns
/// ownership of a heap-allocated C string that must be freed by the caller.

pub extern "C" fn rssn_json_poly_division_multivariate(
    dividend_json: *const c_char,
    divisors_json: *const c_char,
    order_json: *const c_char,
) -> *mut c_char {

    let dividend: Option<
        SparsePolynomial,
    > = from_json_string(dividend_json);

    let divisors: Option<
        Vec<SparsePolynomial>,
    > = from_json_string(divisors_json);

    let order: Option<MonomialOrder> =
        from_json_string(order_json);

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

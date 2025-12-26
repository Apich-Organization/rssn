use std::ffi::c_char;

use crate::ffi_apis::common::*;
use crate::symbolic::core::SparsePolynomial;
use crate::symbolic::grobner::buchberger;
use crate::symbolic::grobner::poly_division_multivariate;
use crate::symbolic::grobner::MonomialOrder;

#[no_mangle]

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

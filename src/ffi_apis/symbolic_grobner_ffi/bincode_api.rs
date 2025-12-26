use crate::ffi_apis::common::*;
use crate::symbolic::core::SparsePolynomial;
use crate::symbolic::grobner::buchberger;
use crate::symbolic::grobner::poly_division_multivariate;
use crate::symbolic::grobner::MonomialOrder;

#[no_mangle]

pub extern "C" fn rssn_bincode_buchberger(
    basis_buf : BincodeBuffer,
    order_buf : BincodeBuffer,
) -> BincodeBuffer {

    let basis : Option<Vec<SparsePolynomial>> = from_bincode_buffer(&basis_buf);

    let order : Option<MonomialOrder> = from_bincode_buffer(&order_buf);

    if let (Some(b), Some(o)) = (basis, order) {

        match buchberger(&b, o) {
            | Ok(result) => to_bincode_buffer(&result),
            | Err(_) => BincodeBuffer::empty(),
        }
    } else {

        BincodeBuffer::empty()
    }
}

#[no_mangle]

pub extern "C" fn rssn_bincode_poly_division_multivariate(
    dividend_buf : BincodeBuffer,
    divisors_buf : BincodeBuffer,
    order_buf : BincodeBuffer,
) -> BincodeBuffer {

    let dividend : Option<SparsePolynomial> = from_bincode_buffer(&dividend_buf);

    let divisors : Option<Vec<SparsePolynomial>> = from_bincode_buffer(&divisors_buf);

    let order : Option<MonomialOrder> = from_bincode_buffer(&order_buf);

    if let (Some(d), Some(divs), Some(o)) = (
        dividend,
        divisors,
        order,
    ) {

        match poly_division_multivariate(&d, &divs, o) {
            | Ok(result) => to_bincode_buffer(&result),
            | Err(_) => BincodeBuffer::empty(),
        }
    } else {

        BincodeBuffer::empty()
    }
}

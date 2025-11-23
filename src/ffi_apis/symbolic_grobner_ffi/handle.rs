use crate::symbolic::core::SparsePolynomial;
use crate::symbolic::grobner::{buchberger, poly_division_multivariate, MonomialOrder};

#[no_mangle]
pub extern "C" fn rssn_buchberger_handle(
    basis: *const Vec<SparsePolynomial>,
    order: MonomialOrder
) -> *mut Vec<SparsePolynomial> {
    let basis_ref = unsafe { &*basis };
    
    match buchberger(basis_ref, order) {
        Ok(result) => Box::into_raw(Box::new(result)),
        Err(_) => std::ptr::null_mut(),
    }
}

#[no_mangle]
pub extern "C" fn rssn_poly_division_multivariate_handle(
    dividend: *const SparsePolynomial,
    divisors: *const Vec<SparsePolynomial>,
    order: MonomialOrder
) -> *mut (Vec<SparsePolynomial>, SparsePolynomial) {
    let dividend_ref = unsafe { &*dividend };
    let divisors_ref = unsafe { &*divisors };
    
    match poly_division_multivariate(dividend_ref, divisors_ref, order) {
        Ok(result) => Box::into_raw(Box::new(result)),
        Err(_) => std::ptr::null_mut(),
    }
}

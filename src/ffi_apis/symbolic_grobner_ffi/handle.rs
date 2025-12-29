use crate::symbolic::core::SparsePolynomial;
use crate::symbolic::grobner::buchberger;
use crate::symbolic::grobner::poly_division_multivariate;
use crate::symbolic::grobner::MonomialOrder;

#[no_mangle]

/// Computes a Gröbner basis using Buchberger's algorithm and returns it as a heap-allocated vector handle.
///
/// # Arguments
///
/// * `basis` - Pointer to a `Vec<SparsePolynomial>` representing the initial basis.
/// * `order` - [`MonomialOrder`] specifying the term ordering.
///
/// # Returns
///
/// A pointer to a heap-allocated `Vec<SparsePolynomial>` representing a Gröbner
/// basis, or null if the computation fails.
///
/// # Safety
///
/// This function is unsafe because it dereferences a raw pointer and returns
/// ownership of a heap-allocated vector that must be freed by the caller.

pub extern "C" fn rssn_buchberger_handle(
    basis: *const Vec<SparsePolynomial>,
    order: MonomialOrder,
) -> *mut Vec<SparsePolynomial> {

    let basis_ref = unsafe {

        &*basis
    };

    match buchberger(basis_ref, order) {
        | Ok(result) => {
            Box::into_raw(Box::new(
                result,
            ))
        },
        | Err(_) => {
            std::ptr::null_mut()
        },
    }
}

#[no_mangle]

/// Divides a multivariate polynomial by a list of divisors under a given monomial order
/// and returns quotient polynomials and remainder as a heap-allocated tuple.
///
/// # Arguments
///
/// * `dividend` - Pointer to the dividend `SparsePolynomial`.
/// * `divisors` - Pointer to a `Vec<SparsePolynomial>` of divisor polynomials.
/// * `order` - [`MonomialOrder`] specifying the term ordering.
///
/// # Returns
///
/// A pointer to a heap-allocated `(Vec<SparsePolynomial>, SparsePolynomial)` containing
/// the quotient polynomials and the remainder, or null if the computation fails.
///
/// # Safety
///
/// This function is unsafe because it dereferences raw pointers and returns ownership
/// of a heap-allocated tuple that must be freed by the caller.

pub extern "C" fn rssn_poly_division_multivariate_handle(
    dividend: *const SparsePolynomial,
    divisors: *const Vec<
        SparsePolynomial,
    >,
    order: MonomialOrder,
) -> *mut (
    Vec<SparsePolynomial>,
    SparsePolynomial,
) {

    let dividend_ref = unsafe {

        &*dividend
    };

    let divisors_ref = unsafe {

        &*divisors
    };

    match poly_division_multivariate(
        dividend_ref,
        divisors_ref,
        order,
    ) {
        | Ok(result) => {
            Box::into_raw(Box::new(
                result,
            ))
        },
        | Err(_) => {
            std::ptr::null_mut()
        },
    }
}

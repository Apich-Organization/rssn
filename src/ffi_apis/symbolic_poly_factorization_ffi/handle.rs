use crate::symbolic::finite_field::FiniteFieldPolynomial;
use crate::symbolic::poly_factorization::*;

/// Factors a polynomial over a finite field (Handle)
#[no_mangle]

pub extern "C" fn rssn_factor_gf_handle(
    poly: *const FiniteFieldPolynomial
) -> *mut Vec<FiniteFieldPolynomial> {

    let poly_ref = unsafe {

        if poly.is_null() {

            return std::ptr::null_mut(
            );
        }

        &*poly
    };

    match factor_gf(poly_ref) {
        | Ok(factors) => {
            Box::into_raw(Box::new(
                factors,
            ))
        },
        | Err(_) => {
            std::ptr::null_mut()
        },
    }
}

/// Computes square-free factorization (Handle)
#[no_mangle]

pub extern "C" fn rssn_square_free_factorization_gf_handle(
    poly: *const FiniteFieldPolynomial
) -> *mut Vec<(
    FiniteFieldPolynomial,
    usize,
)> {

    let poly_ref = unsafe {

        if poly.is_null() {

            return std::ptr::null_mut(
            );
        }

        &*poly
    };

    match square_free_factorization_gf(
        poly_ref.clone(),
    ) {
        | Ok(factors) => {
            Box::into_raw(Box::new(
                factors,
            ))
        },
        | Err(_) => {
            std::ptr::null_mut()
        },
    }
}

/// Computes polynomial GCD over finite field (Handle)
#[no_mangle]

pub extern "C" fn rssn_poly_gcd_gf_handle(
    a: *const FiniteFieldPolynomial,
    b: *const FiniteFieldPolynomial,
) -> *mut FiniteFieldPolynomial {

    let a_ref = unsafe {

        if a.is_null() {

            return std::ptr::null_mut(
            );
        }

        &*a
    };

    let b_ref = unsafe {

        if b.is_null() {

            return std::ptr::null_mut(
            );
        }

        &*b
    };

    match poly_gcd_gf(
        a_ref.clone(),
        b_ref.clone(),
    ) {
        | Ok(gcd) => {
            Box::into_raw(Box::new(gcd))
        },
        | Err(_) => {
            std::ptr::null_mut()
        },
    }
}

/// Computes polynomial derivative over finite field (Handle)
#[no_mangle]

pub extern "C" fn rssn_poly_derivative_gf_handle(
    poly: *const FiniteFieldPolynomial
) -> *mut FiniteFieldPolynomial {

    let poly_ref = unsafe {

        if poly.is_null() {

            return std::ptr::null_mut(
            );
        }

        &*poly
    };

    let derivative =
        poly_derivative_gf(poly_ref);

    Box::into_raw(Box::new(derivative))
}

/// Frees a vector of polynomials (Handle)
#[no_mangle]

pub extern "C" fn rssn_free_poly_vec_handle(
    ptr: *mut Vec<
        FiniteFieldPolynomial,
    >
) {

    if !ptr.is_null() {

        unsafe {

            let _ = Box::from_raw(ptr);
        }
    }
}

/// Frees a vector of polynomial-multiplicity pairs (Handle)
#[no_mangle]

pub extern "C" fn rssn_free_poly_mult_vec_handle(
    ptr: *mut Vec<(
        FiniteFieldPolynomial,
        usize,
    )>
) {

    if !ptr.is_null() {

        unsafe {

            let _ = Box::from_raw(ptr);
        }
    }
}

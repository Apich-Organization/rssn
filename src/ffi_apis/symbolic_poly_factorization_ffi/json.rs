use crate::ffi_apis::common::*;
use crate::symbolic::finite_field::FiniteFieldPolynomial;
use crate::symbolic::poly_factorization::*;
use std::ffi::c_char;

/// Factors a polynomial over a finite field (JSON)
#[no_mangle]
pub extern "C" fn rssn_json_factor_gf(poly_json: *const c_char) -> *mut c_char {
    let poly: Option<FiniteFieldPolynomial> = from_json_string(poly_json);

    if let Some(p) = poly {
        match factor_gf(&p) {
            Ok(factors) => to_json_string(&factors),
            Err(_) => std::ptr::null_mut(),
        }
    } else {
        std::ptr::null_mut()
    }
}

/// Computes square-free factorization (JSON)
#[no_mangle]
pub extern "C" fn rssn_json_square_free_factorization_gf(
    poly_json: *const c_char,
) -> *mut c_char {
    let poly: Option<FiniteFieldPolynomial> = from_json_string(poly_json);

    if let Some(p) = poly {
        match square_free_factorization_gf(p) {
            Ok(factors) => to_json_string(&factors),
            Err(_) => std::ptr::null_mut(),
        }
    } else {
        std::ptr::null_mut()
    }
}

/// Computes polynomial GCD over finite field (JSON)
#[no_mangle]
pub extern "C" fn rssn_json_poly_gcd_gf(
    a_json: *const c_char,
    b_json: *const c_char,
) -> *mut c_char {
    let a: Option<FiniteFieldPolynomial> = from_json_string(a_json);
    let b: Option<FiniteFieldPolynomial> = from_json_string(b_json);

    if let (Some(poly_a), Some(poly_b)) = (a, b) {
        match poly_gcd_gf(poly_a, poly_b) {
            Ok(gcd) => to_json_string(&gcd),
            Err(_) => std::ptr::null_mut(),
        }
    } else {
        std::ptr::null_mut()
    }
}

/// Computes polynomial derivative over finite field (JSON)
#[no_mangle]
pub extern "C" fn rssn_json_poly_derivative_gf(poly_json: *const c_char) -> *mut c_char {
    let poly: Option<FiniteFieldPolynomial> = from_json_string(poly_json);

    if let Some(p) = poly {
        let derivative = poly_derivative_gf(&p);
        to_json_string(&derivative)
    } else {
        std::ptr::null_mut()
    }
}

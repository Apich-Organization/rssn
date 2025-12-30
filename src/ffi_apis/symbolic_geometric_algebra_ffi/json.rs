use std::ffi::c_char;

use crate::ffi_apis::common::{from_json_string, to_json_string};
use crate::symbolic::core::Expr;
use crate::symbolic::geometric_algebra::Multivector;

/// Creates a new scalar multivector (JSON)
#[unsafe(no_mangle)]

pub extern "C" fn rssn_json_multivector_scalar(
    p: u32,
    q: u32,
    r: u32,
    value_json: *const c_char,
) -> *mut c_char {

    let value: Option<Expr> =
        from_json_string(value_json);

    if let Some(val) = value {

        let mv = Multivector::scalar(
            (p, q, r),
            val,
        );

        to_json_string(&mv)
    } else {

        std::ptr::null_mut()
    }
}

/// Computes geometric product (JSON)
#[unsafe(no_mangle)]

pub extern "C" fn rssn_json_multivector_geometric_product(
    a_json: *const c_char,
    b_json: *const c_char,
) -> *mut c_char {

    let a: Option<Multivector> =
        from_json_string(a_json);

    let b: Option<Multivector> =
        from_json_string(b_json);

    match (a, b)
    { (Some(mv_a), Some(mv_b)) => {

        let result = mv_a
            .geometric_product(&mv_b);

        to_json_string(&result)
    } _ => {

        std::ptr::null_mut()
    }}
}

/// Computes outer product (JSON)
#[unsafe(no_mangle)]

pub extern "C" fn rssn_json_multivector_outer_product(
    a_json: *const c_char,
    b_json: *const c_char,
) -> *mut c_char {

    let a: Option<Multivector> =
        from_json_string(a_json);

    let b: Option<Multivector> =
        from_json_string(b_json);

    match (a, b)
    { (Some(mv_a), Some(mv_b)) => {

        let result =
            mv_a.outer_product(&mv_b);

        to_json_string(&result)
    } _ => {

        std::ptr::null_mut()
    }}
}

/// Computes inner product (JSON)
#[unsafe(no_mangle)]

pub extern "C" fn rssn_json_multivector_inner_product(
    a_json: *const c_char,
    b_json: *const c_char,
) -> *mut c_char {

    let a: Option<Multivector> =
        from_json_string(a_json);

    let b: Option<Multivector> =
        from_json_string(b_json);

    match (a, b)
    { (Some(mv_a), Some(mv_b)) => {

        let result =
            mv_a.inner_product(&mv_b);

        to_json_string(&result)
    } _ => {

        std::ptr::null_mut()
    }}
}

/// Computes reverse (JSON)
#[unsafe(no_mangle)]

pub extern "C" fn rssn_json_multivector_reverse(
    mv_json: *const c_char
) -> *mut c_char {

    let mv: Option<Multivector> =
        from_json_string(mv_json);

    if let Some(multivector) = mv {

        let result =
            multivector.reverse();

        to_json_string(&result)
    } else {

        std::ptr::null_mut()
    }
}

/// Computes grade projection (JSON)
#[unsafe(no_mangle)]

pub extern "C" fn rssn_json_multivector_grade_projection(
    mv_json: *const c_char,
    grade: u32,
) -> *mut c_char {

    let mv: Option<Multivector> =
        from_json_string(mv_json);

    if let Some(multivector) = mv {

        let result = multivector
            .grade_projection(grade);

        to_json_string(&result)
    } else {

        std::ptr::null_mut()
    }
}

/// Computes magnitude (JSON)
#[unsafe(no_mangle)]

pub extern "C" fn rssn_json_multivector_magnitude(
    mv_json: *const c_char
) -> *mut c_char {

    let mv: Option<Multivector> =
        from_json_string(mv_json);

    if let Some(multivector) = mv {

        let result =
            multivector.magnitude();

        to_json_string(&result)
    } else {

        std::ptr::null_mut()
    }
}

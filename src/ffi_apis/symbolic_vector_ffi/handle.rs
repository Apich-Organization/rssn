use crate::symbolic::vector::*;

#[no_mangle]

pub extern "C" fn rssn_vector_magnitude_handle(
    v: *const Vector,
) -> *mut crate::symbolic::core::Expr {

    let v_ref = unsafe {

        &*v
    };

    let result = v_ref.magnitude();

    Box::into_raw(Box::new(result))
}

#[no_mangle]

pub extern "C" fn rssn_vector_dot_handle(
    v1: *const Vector,
    v2: *const Vector,
) -> *mut crate::symbolic::core::Expr {

    let v1_ref = unsafe {

        &*v1
    };

    let v2_ref = unsafe {

        &*v2
    };

    let result = v1_ref.dot(v2_ref);

    Box::into_raw(Box::new(result))
}

#[no_mangle]

pub extern "C" fn rssn_vector_cross_handle(v1: *const Vector, v2: *const Vector) -> *mut Vector {

    let v1_ref = unsafe {

        &*v1
    };

    let v2_ref = unsafe {

        &*v2
    };

    let result = v1_ref.cross(v2_ref);

    Box::into_raw(Box::new(result))
}

#[no_mangle]

pub extern "C" fn rssn_vector_normalize_handle(v: *const Vector) -> *mut Vector {

    let v_ref = unsafe {

        &*v
    };

    let result = v_ref.normalize();

    Box::into_raw(Box::new(result))
}

use crate::symbolic::core::Expr;
use crate::symbolic::geometric_algebra::Multivector;

/// Creates a new scalar multivector (Handle)
#[no_mangle]

pub extern "C" fn rssn_multivector_scalar_handle(
    p : u32,
    q : u32,
    r : u32,
    value : *const Expr,
) -> *mut Multivector {

    let value_ref = unsafe {

        if value.is_null() {

            return std::ptr::null_mut(
            );
        }

        &*value
    };

    let mv = Multivector::scalar(
        (p, q, r),
        value_ref.clone(),
    );

    Box::into_raw(Box::new(mv))
}

/// Computes geometric product (Handle)
#[no_mangle]

pub extern "C" fn rssn_multivector_geometric_product_handle(
    a : *const Multivector,
    b : *const Multivector,
) -> *mut Multivector {

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

    let result =
        a_ref.geometric_product(b_ref);

    Box::into_raw(Box::new(result))
}

/// Computes outer product (Handle)
#[no_mangle]

pub extern "C" fn rssn_multivector_outer_product_handle(
    a : *const Multivector,
    b : *const Multivector,
) -> *mut Multivector {

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

    let result =
        a_ref.outer_product(b_ref);

    Box::into_raw(Box::new(result))
}

/// Computes inner product (Handle)
#[no_mangle]

pub extern "C" fn rssn_multivector_inner_product_handle(
    a : *const Multivector,
    b : *const Multivector,
) -> *mut Multivector {

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

    let result =
        a_ref.inner_product(b_ref);

    Box::into_raw(Box::new(result))
}

/// Computes reverse (Handle)
#[no_mangle]

pub extern "C" fn rssn_multivector_reverse_handle(
    mv : *const Multivector
) -> *mut Multivector {

    let mv_ref = unsafe {

        if mv.is_null() {

            return std::ptr::null_mut(
            );
        }

        &*mv
    };

    let result = mv_ref.reverse();

    Box::into_raw(Box::new(result))
}

/// Computes grade projection (Handle)
#[no_mangle]

pub extern "C" fn rssn_multivector_grade_projection_handle(
    mv : *const Multivector,
    grade : u32,
) -> *mut Multivector {

    let mv_ref = unsafe {

        if mv.is_null() {

            return std::ptr::null_mut(
            );
        }

        &*mv
    };

    let result =
        mv_ref.grade_projection(grade);

    Box::into_raw(Box::new(result))
}

/// Computes magnitude (Handle)
#[no_mangle]

pub extern "C" fn rssn_multivector_magnitude_handle(
    mv : *const Multivector
) -> *mut Expr {

    let mv_ref = unsafe {

        if mv.is_null() {

            return std::ptr::null_mut(
            );
        }

        &*mv
    };

    let result = mv_ref.magnitude();

    Box::into_raw(Box::new(result))
}

/// Frees a multivector (Handle)
#[no_mangle]

pub extern "C" fn rssn_free_multivector_handle(
    ptr : *mut Multivector
) {

    if !ptr.is_null() {

        unsafe {

            let _ = Box::from_raw(ptr);
        }
    }
}

//! Handle-based FFI API for relativity functions.

use crate::symbolic::core::Expr;
use crate::symbolic::relativity;

/// Structure to hold two expressions (e.g., transformed x and t).
#[repr(C)]

pub struct ExprPair {
    pub first: *mut Expr,
    pub second: *mut Expr,
}

/// Calculates the Lorentz factor.
#[no_mangle]

pub unsafe extern "C" fn rssn_lorentz_factor(
    velocity: *const Expr
) -> *mut Expr {

    if velocity.is_null() {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(
        relativity::lorentz_factor(
            &*velocity,
        ),
    ))
}

/// Performs a Lorentz transformation in the x-direction.
#[no_mangle]

pub unsafe extern "C" fn rssn_lorentz_transformation_x(
    x: *const Expr,
    t: *const Expr,
    v: *const Expr,
) -> ExprPair {

    if x.is_null()
        || t.is_null()
        || v.is_null()
    {

        return ExprPair {
            first: std::ptr::null_mut(),
            second: std::ptr::null_mut(
            ),
        };
    }

    let (xp, tp) =
        relativity::lorentz_transformation_x(&*x, &*t, &*v);

    ExprPair {
        first: Box::into_raw(Box::new(
            xp,
        )),
        second: Box::into_raw(
            Box::new(tp),
        ),
    }
}

/// Calculates mass-energy equivalence.
#[no_mangle]

pub unsafe extern "C" fn rssn_mass_energy_equivalence(
    mass: *const Expr
) -> *mut Expr {

    if mass.is_null() {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(
        relativity::mass_energy_equivalence(&*mass),
    ))
}

/// Calculates Schwarzschild radius.
#[no_mangle]

pub unsafe extern "C" fn rssn_schwarzschild_radius(
    mass: *const Expr
) -> *mut Expr {

    if mass.is_null() {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(
        relativity::schwarzschild_radius(&*mass),
    ))
}

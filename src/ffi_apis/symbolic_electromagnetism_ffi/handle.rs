//! Handle-based FFI API for electromagnetism functions.

use std::ffi::CStr;
use std::os::raw::c_char;

use crate::symbolic::core::Expr;
use crate::symbolic::electromagnetism;
use crate::symbolic::vector::Vector;

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

unsafe fn c_str_to_str<'a>(
    s: *const c_char
) -> Option<&'a str> {

    if s.is_null() {

        None
    } else {

        CStr::from_ptr(s)
            .to_str()
            .ok()
    }
}

/// Calculates the Lorentz force.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_lorentz_force(
    charge: *const Expr,
    e_field: *const Vector,
    velocity: *const Vector,
    b_field: *const Vector,
) -> *mut Vector {

    if charge.is_null()
        || e_field.is_null()
        || velocity.is_null()
        || b_field.is_null()
    {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(
        electromagnetism::lorentz_force(
            &*charge,
            &*e_field,
            &*velocity,
            &*b_field,
        ),
    ))
}

/// Calculates the Poynting vector.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_poynting_vector(
    e_field: *const Vector,
    b_field: *const Vector,
) -> *mut Vector {

    if e_field.is_null()
        || b_field.is_null()
    {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(
        electromagnetism::poynting_vector(&*e_field, &*b_field),
    ))
}

/// Calculates energy density.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_electromagnetic_energy_density(
    e_field: *const Vector,
    b_field: *const Vector,
) -> *mut Expr {

    if e_field.is_null()
        || b_field.is_null()
    {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(
        electromagnetism::energy_density(&*e_field, &*b_field),
    ))
}

/// Computes magnetic field from vector potential.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_magnetic_field_from_vector_potential(
    a: *const Vector,
    x: *const c_char,
    y: *const c_char,
    z: *const c_char,
) -> *mut Vector {

    if a.is_null()
        || x.is_null()
        || y.is_null()
        || z.is_null()
    {

        return std::ptr::null_mut();
    }

    let xs = match c_str_to_str(x) {
        | Some(s) => s,
        | None => {
            return std::ptr::null_mut()
        },
    };

    let ys = match c_str_to_str(y) {
        | Some(s) => s,
        | None => {
            return std::ptr::null_mut()
        },
    };

    let zs = match c_str_to_str(z) {
        | Some(s) => s,
        | None => {
            return std::ptr::null_mut()
        },
    };

    Box::into_raw(Box::new(
        electromagnetism::magnetic_field_from_vector_potential(&*a, (xs, ys, zs)),
    ))
}

/// Computes electric field from scalar and vector potentials.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_electric_field_from_potentials(
    v: *const Expr,
    a: *const Vector,
    x: *const c_char,
    y: *const c_char,
    z: *const c_char,
    t: *const c_char,
) -> *mut Vector {

    if v.is_null()
        || a.is_null()
        || x.is_null()
        || y.is_null()
        || z.is_null()
        || t.is_null()
    {

        return std::ptr::null_mut();
    }

    let xs = match c_str_to_str(x) {
        | Some(s) => s,
        | None => {
            return std::ptr::null_mut()
        },
    };

    let ys = match c_str_to_str(y) {
        | Some(s) => s,
        | None => {
            return std::ptr::null_mut()
        },
    };

    let zs = match c_str_to_str(z) {
        | Some(s) => s,
        | None => {
            return std::ptr::null_mut()
        },
    };

    let ts = match c_str_to_str(t) {
        | Some(s) => s,
        | None => {
            return std::ptr::null_mut()
        },
    };

    Box::into_raw(Box::new(
        electromagnetism::electric_field_from_potentials(
            &*v,
            &*a,
            (xs, ys, zs),
            ts,
        ),
    ))
}

/// Calculates Coulomb's Law field.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_coulombs_law(
    charge: *const Expr,
    r: *const Vector,
) -> *mut Vector {

    if charge.is_null() || r.is_null() {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(
        electromagnetism::coulombs_law(
            &*charge,
            &*r,
        ),
    ))
}

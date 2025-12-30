//! Handle-based FFI API for classical mechanics functions.

use std::ffi::CStr;
use std::os::raw::c_char;

use crate::symbolic::classical_mechanics;
use crate::symbolic::core::Expr;
use crate::symbolic::vector::Vector;
use crate::symbolic::vector_calculus::ParametricCurve;

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

unsafe fn c_str_to_str<'a>(
    s: *const c_char
) -> Option<&'a str> { unsafe {

    if s.is_null() {

        None
    } else {

        CStr::from_ptr(s)
            .to_str()
            .ok()
    }
}}

/// Calculates kinetic energy: 1/2 * m * v^2.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_kinetic_energy(
    mass: *const Expr,
    velocity: *const Expr,
) -> *mut Expr { unsafe {

    if mass.is_null()
        || velocity.is_null()
    {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(
        classical_mechanics::kinetic_energy(&*mass, &*velocity),
    ))
}}

/// Calculates Lagrangian: T - V.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_lagrangian(
    t: *const Expr,
    v: *const Expr,
) -> *mut Expr { unsafe {

    if t.is_null() || v.is_null() {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(
        classical_mechanics::lagrangian(
            &*t, &*v,
        ),
    ))
}}

/// Calculates Hamiltonian: T + V.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_hamiltonian(
    t: *const Expr,
    v: *const Expr,
) -> *mut Expr { unsafe {

    if t.is_null() || v.is_null() {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(
        classical_mechanics::hamiltonian(&*t, &*v),
    ))
}}

/// Computes Euler-Lagrange equation.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_euler_lagrange_equation(
    lagrangian: *const Expr,
    q: *const c_char,
    q_dot: *const c_char,
    t_var: *const c_char,
) -> *mut Expr { unsafe {

    if lagrangian.is_null()
        || q.is_null()
        || q_dot.is_null()
        || t_var.is_null()
    {

        return std::ptr::null_mut();
    }

    let q_str = match c_str_to_str(q) {
        | Some(s) => s,
        | None => {
            return std::ptr::null_mut()
        },
    };

    let q_dot_str = match c_str_to_str(
        q_dot,
    ) {
        | Some(s) => s,
        | None => {
            return std::ptr::null_mut()
        },
    };

    let t_str = match c_str_to_str(
        t_var,
    ) {
        | Some(s) => s,
        | None => {
            return std::ptr::null_mut()
        },
    };

    Box::into_raw(Box::new(
        classical_mechanics::euler_lagrange_equation(
            &*lagrangian,
            q_str,
            q_dot_str,
            t_str,
        ),
    ))
}}

/// Calculates torque: r x F.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_torque(
    r: *const Vector,
    force: *const Vector,
) -> *mut Vector { unsafe {

    if r.is_null() || force.is_null() {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(
        classical_mechanics::torque(
            &*r,
            &*force,
        ),
    ))
}}

/// Calculates power: F . v.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_power(
    force: *const Vector,
    velocity: *const Vector,
) -> *mut Expr { unsafe {

    if force.is_null()
        || velocity.is_null()
    {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(
        classical_mechanics::power(
            &*force,
            &*velocity,
        ),
    ))
}}

/// Calculates work done by a variable force field along a path.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_work_line_integral(
    force_field: *const Vector,
    path: *const ParametricCurve,
) -> *mut Expr { unsafe {

    if force_field.is_null()
        || path.is_null()
    {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(
        classical_mechanics::work_line_integral(
            &*force_field,
            &*path,
        ),
    ))
}}

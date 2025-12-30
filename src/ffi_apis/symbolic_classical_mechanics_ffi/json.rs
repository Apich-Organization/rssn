//! JSON-based FFI API for classical mechanics functions.

use std::os::raw::c_char;

use crate::ffi_apis::common::from_json_string;
use crate::ffi_apis::common::to_json_string;
use crate::symbolic::classical_mechanics;
use crate::symbolic::core::Expr;

/// Calculates kinetic energy using JSON.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_json_kinetic_energy(
    mass_json: *const c_char,
    velocity_json: *const c_char,
) -> *mut c_char {

    let mass: Option<Expr> =
        from_json_string(mass_json);

    let velocity: Option<Expr> =
        from_json_string(velocity_json);

    match (mass, velocity)
    { (Some(m), Some(v)) => {

        to_json_string(&classical_mechanics::kinetic_energy(&m, &v))
    } _ => {

        std::ptr::null_mut()
    }}
}

/// Computes Euler-Lagrange equation using JSON.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_json_euler_lagrange_equation(
    lagrangian_json: *const c_char,
    q: *const c_char,
    q_dot: *const c_char,
    t_var: *const c_char,
) -> *mut c_char {

    let lagrangian: Option<Expr> =
        from_json_string(
            lagrangian_json,
        );

    let q_str = unsafe {

        if q.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(q)
                .to_str()
                .ok()
        }
    };

    let q_dot_str = unsafe {

        if q_dot.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(
                q_dot,
            )
            .to_str()
            .ok()
        }
    };

    let t_str = unsafe {

        if t_var.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(
                t_var,
            )
            .to_str()
            .ok()
        }
    };

    match (
        lagrangian,
        q_str,
        q_dot_str,
        t_str,
    ) { (
        Some(l),
        Some(qs),
        Some(qds),
        Some(ts),
    ) => {

        to_json_string(&classical_mechanics::euler_lagrange_equation(&l, qs, qds, ts))
    } _ => {

        std::ptr::null_mut()
    }}
}

//! Bincode-based FFI API for classical mechanics functions.

use crate::ffi_apis::common::*;
use crate::symbolic::classical_mechanics;
use crate::symbolic::core::Expr;
use std::os::raw::c_char;

/// Calculates kinetic energy using Bincode.
#[no_mangle]

pub extern "C" fn rssn_bincode_kinetic_energy(
    mass_buf: BincodeBuffer,
    velocity_buf: BincodeBuffer,
) -> BincodeBuffer {

    let mass: Option<Expr> = from_bincode_buffer(&mass_buf);

    let velocity: Option<Expr> = from_bincode_buffer(&velocity_buf);

    if let (Some(m), Some(v)) = (mass, velocity) {

        to_bincode_buffer(&classical_mechanics::kinetic_energy(&m, &v))
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes Euler-Lagrange equation using Bincode.
#[no_mangle]

pub extern "C" fn rssn_bincode_euler_lagrange_equation(
    lagrangian_buf: BincodeBuffer,
    q: *const c_char,
    q_dot: *const c_char,
    t_var: *const c_char,
) -> BincodeBuffer {

    let lagrangian: Option<Expr> = from_bincode_buffer(&lagrangian_buf);

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

            std::ffi::CStr::from_ptr(q_dot)
                .to_str()
                .ok()
        }
    };

    let t_str = unsafe {

        if t_var.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(t_var)
                .to_str()
                .ok()
        }
    };

    if let (Some(l), Some(qs), Some(qds), Some(ts)) = (lagrangian, q_str, q_dot_str, t_str) {

        to_bincode_buffer(&classical_mechanics::euler_lagrange_equation(
            &l, qs, qds, ts,
        ))
    } else {

        BincodeBuffer::empty()
    }
}

//! Bincode-based FFI API for calculus of variations functions.

use crate::ffi_apis::common::*;
use crate::symbolic::calculus_of_variations;
use crate::symbolic::core::Expr;
use std::os::raw::c_char;

/// Computes the Euler-Lagrange equation using Bincode.
#[no_mangle]

pub extern "C" fn rssn_bincode_euler_lagrange(
    lagrangian_buf: BincodeBuffer,
    func: *const c_char,
    var: *const c_char,
) -> BincodeBuffer {

    let lagrangian: Option<Expr> = from_bincode_buffer(&lagrangian_buf);

    let func_str = unsafe {

        if func.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(func).to_str().ok()
        }
    };

    let var_str = unsafe {

        if var.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(var).to_str().ok()
        }
    };

    if let (Some(l), Some(f), Some(v)) = (lagrangian, func_str, var_str) {

        let result = calculus_of_variations::euler_lagrange(&l, f, v);

        to_bincode_buffer(&result)
    } else {

        BincodeBuffer::empty()
    }
}

/// Generates and attempts to solve the Euler-Lagrange equation using Bincode.
#[no_mangle]

pub extern "C" fn rssn_bincode_solve_euler_lagrange(
    lagrangian_buf: BincodeBuffer,
    func: *const c_char,
    var: *const c_char,
) -> BincodeBuffer {

    let lagrangian: Option<Expr> = from_bincode_buffer(&lagrangian_buf);

    let func_str = unsafe {

        if func.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(func).to_str().ok()
        }
    };

    let var_str = unsafe {

        if var.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(var).to_str().ok()
        }
    };

    if let (Some(l), Some(f), Some(v)) = (lagrangian, func_str, var_str) {

        let result = calculus_of_variations::solve_euler_lagrange(&l, f, v);

        to_bincode_buffer(&result)
    } else {

        BincodeBuffer::empty()
    }
}

/// Applies Hamilton's Principle using Bincode.
#[no_mangle]

pub extern "C" fn rssn_bincode_hamiltons_principle(
    lagrangian_buf: BincodeBuffer,
    func: *const c_char,
    var: *const c_char,
) -> BincodeBuffer {

    let lagrangian: Option<Expr> = from_bincode_buffer(&lagrangian_buf);

    let func_str = unsafe {

        if func.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(func).to_str().ok()
        }
    };

    let var_str = unsafe {

        if var.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(var).to_str().ok()
        }
    };

    if let (Some(l), Some(f), Some(v)) = (lagrangian, func_str, var_str) {

        let result = calculus_of_variations::hamiltons_principle(&l, f, v);

        to_bincode_buffer(&result)
    } else {

        BincodeBuffer::empty()
    }
}

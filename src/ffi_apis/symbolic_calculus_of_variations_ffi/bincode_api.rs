//! Bincode-based FFI API for calculus of variations functions.

use std::os::raw::c_char;

use crate::ffi_apis::common::BincodeBuffer;
use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::symbolic::calculus_of_variations;
use crate::symbolic::core::Expr;

/// Computes the Euler-Lagrange equation using Bincode.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_euler_lagrange(
    lagrangian_buf: BincodeBuffer,
    func: *const c_char,
    var: *const c_char,
) -> BincodeBuffer {

    let lagrangian: Option<Expr> =
        from_bincode_buffer(
            &lagrangian_buf,
        );

    let func_str = unsafe {

        if func.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(
                func,
            )
            .to_str()
            .ok()
        }
    };

    let var_str = unsafe {

        if var.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(
                var,
            )
            .to_str()
            .ok()
        }
    };

    match (
        lagrangian,
        func_str,
        var_str,
    ) {
        | (
            Some(l),
            Some(f),
            Some(v),
        ) => {

            let result = calculus_of_variations::euler_lagrange(&l, f, v);

            to_bincode_buffer(&result)
        },
        | _ => BincodeBuffer::empty(),
    }
}

/// Generates and attempts to solve the Euler-Lagrange equation using Bincode.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_solve_euler_lagrange(
    lagrangian_buf: BincodeBuffer,
    func: *const c_char,
    var: *const c_char,
) -> BincodeBuffer {

    let lagrangian: Option<Expr> =
        from_bincode_buffer(
            &lagrangian_buf,
        );

    let func_str = unsafe {

        if func.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(
                func,
            )
            .to_str()
            .ok()
        }
    };

    let var_str = unsafe {

        if var.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(
                var,
            )
            .to_str()
            .ok()
        }
    };

    match (
        lagrangian,
        func_str,
        var_str,
    ) {
        | (
            Some(l),
            Some(f),
            Some(v),
        ) => {

            let result = calculus_of_variations::solve_euler_lagrange(&l, f, v);

            to_bincode_buffer(&result)
        },
        | _ => BincodeBuffer::empty(),
    }
}

/// Applies Hamilton's Principle using Bincode.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_hamiltons_principle(
    lagrangian_buf: BincodeBuffer,
    func: *const c_char,
    var: *const c_char,
) -> BincodeBuffer {

    let lagrangian: Option<Expr> =
        from_bincode_buffer(
            &lagrangian_buf,
        );

    let func_str = unsafe {

        if func.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(
                func,
            )
            .to_str()
            .ok()
        }
    };

    let var_str = unsafe {

        if var.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(
                var,
            )
            .to_str()
            .ok()
        }
    };

    match (
        lagrangian,
        func_str,
        var_str,
    ) {
        | (
            Some(l),
            Some(f),
            Some(v),
        ) => {

            let result = calculus_of_variations::hamiltons_principle(&l, f, v);

            to_bincode_buffer(&result)
        },
        | _ => BincodeBuffer::empty(),
    }
}

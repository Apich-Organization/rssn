//! Bincode-based FFI API for symbolic ODE functions.

use crate::ffi_apis::common::*;
use crate::symbolic::core::Expr;
use crate::symbolic::ode;
use std::os::raw::c_char;

/// Solves an ODE using Bincode.
#[no_mangle]
pub extern "C" fn rssn_bincode_solve_ode(
    ode_buf: BincodeBuffer,
    func: *const c_char,
    var: *const c_char,
) -> BincodeBuffer {
    let ode_expr: Option<Expr> = from_bincode_buffer(&ode_buf);
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

    if let (Some(ode), Some(f), Some(v)) = (ode_expr, func_str, var_str) {
        let result = ode::solve_ode(&ode, f, v, None);
        to_bincode_buffer(&result)
    } else {
        BincodeBuffer::empty()
    }
}

/// Solves a separable ODE using Bincode.
#[no_mangle]
pub extern "C" fn rssn_bincode_solve_separable_ode(
    equation_buf: BincodeBuffer,
    func: *const c_char,
    var: *const c_char,
) -> BincodeBuffer {
    let equation: Option<Expr> = from_bincode_buffer(&equation_buf);
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

    if let (Some(eq), Some(f), Some(v)) = (equation, func_str, var_str) {
        match ode::solve_separable_ode(&eq, f, v) {
            Some(result) => to_bincode_buffer(&result),
            None => BincodeBuffer::empty(),
        }
    } else {
        BincodeBuffer::empty()
    }
}

/// Solves a first-order linear ODE using Bincode.
#[no_mangle]
pub extern "C" fn rssn_bincode_solve_first_order_linear_ode(
    equation_buf: BincodeBuffer,
    func: *const c_char,
    var: *const c_char,
) -> BincodeBuffer {
    let equation: Option<Expr> = from_bincode_buffer(&equation_buf);
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

    if let (Some(eq), Some(f), Some(v)) = (equation, func_str, var_str) {
        match ode::solve_first_order_linear_ode(&eq, f, v) {
            Some(result) => to_bincode_buffer(&result),
            None => BincodeBuffer::empty(),
        }
    } else {
        BincodeBuffer::empty()
    }
}

/// Solves a Bernoulli ODE using Bincode.
#[no_mangle]
pub extern "C" fn rssn_bincode_solve_bernoulli_ode(
    equation_buf: BincodeBuffer,
    func: *const c_char,
    var: *const c_char,
) -> BincodeBuffer {
    let equation: Option<Expr> = from_bincode_buffer(&equation_buf);
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

    if let (Some(eq), Some(f), Some(v)) = (equation, func_str, var_str) {
        match ode::solve_bernoulli_ode(&eq, f, v) {
            Some(result) => to_bincode_buffer(&result),
            None => BincodeBuffer::empty(),
        }
    } else {
        BincodeBuffer::empty()
    }
}

/// Solves a Riccati ODE using Bincode.
#[no_mangle]
pub extern "C" fn rssn_bincode_solve_riccati_ode(
    equation_buf: BincodeBuffer,
    func: *const c_char,
    var: *const c_char,
    y1_buf: BincodeBuffer,
) -> BincodeBuffer {
    let equation: Option<Expr> = from_bincode_buffer(&equation_buf);
    let y1: Option<Expr> = from_bincode_buffer(&y1_buf);
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

    if let (Some(eq), Some(f), Some(v), Some(y)) = (equation, func_str, var_str, y1) {
        match ode::solve_riccati_ode(&eq, f, v, &y) {
            Some(result) => to_bincode_buffer(&result),
            None => BincodeBuffer::empty(),
        }
    } else {
        BincodeBuffer::empty()
    }
}

/// Solves a Cauchy-Euler ODE using Bincode.
#[no_mangle]
pub extern "C" fn rssn_bincode_solve_cauchy_euler_ode(
    equation_buf: BincodeBuffer,
    func: *const c_char,
    var: *const c_char,
) -> BincodeBuffer {
    let equation: Option<Expr> = from_bincode_buffer(&equation_buf);
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

    if let (Some(eq), Some(f), Some(v)) = (equation, func_str, var_str) {
        match ode::solve_cauchy_euler_ode(&eq, f, v) {
            Some(result) => to_bincode_buffer(&result),
            None => BincodeBuffer::empty(),
        }
    } else {
        BincodeBuffer::empty()
    }
}

/// Solves an exact ODE using Bincode.
#[no_mangle]
pub extern "C" fn rssn_bincode_solve_exact_ode(
    equation_buf: BincodeBuffer,
    func: *const c_char,
    var: *const c_char,
) -> BincodeBuffer {
    let equation: Option<Expr> = from_bincode_buffer(&equation_buf);
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

    if let (Some(eq), Some(f), Some(v)) = (equation, func_str, var_str) {
        match ode::solve_exact_ode(&eq, f, v) {
            Some(result) => to_bincode_buffer(&result),
            None => BincodeBuffer::empty(),
        }
    } else {
        BincodeBuffer::empty()
    }
}

/// Solves by reduction of order using Bincode.
#[no_mangle]
pub extern "C" fn rssn_bincode_solve_by_reduction_of_order(
    equation_buf: BincodeBuffer,
    func: *const c_char,
    var: *const c_char,
    y1_buf: BincodeBuffer,
) -> BincodeBuffer {
    let equation: Option<Expr> = from_bincode_buffer(&equation_buf);
    let y1: Option<Expr> = from_bincode_buffer(&y1_buf);
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

    if let (Some(eq), Some(f), Some(v), Some(y)) = (equation, func_str, var_str, y1) {
        match ode::solve_by_reduction_of_order(&eq, f, v, &y) {
            Some(result) => to_bincode_buffer(&result),
            None => BincodeBuffer::empty(),
        }
    } else {
        BincodeBuffer::empty()
    }
}

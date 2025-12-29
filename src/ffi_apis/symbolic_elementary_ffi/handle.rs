//! Handle-based FFI API for symbolic elementary functions.

use std::os::raw::c_char;

use crate::ffi_apis::common::to_c_string;
use crate::symbolic::core::Expr;
use crate::symbolic::elementary;

/// Creates a sine expression: sin(expr).
///
/// # Safety
/// The caller must ensure `expr` is a valid Expr pointer.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_sin(
    expr: *const Expr
) -> *mut Expr {

    if expr.is_null() {

        return std::ptr::null_mut();
    }

    let expr_ref = &*expr;

    Box::into_raw(Box::new(
        elementary::sin(
            expr_ref.clone(),
        ),
    ))
}

/// Creates a cosine expression: cos(expr).
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_cos(
    expr: *const Expr
) -> *mut Expr {

    if expr.is_null() {

        return std::ptr::null_mut();
    }

    let expr_ref = &*expr;

    Box::into_raw(Box::new(
        elementary::cos(
            expr_ref.clone(),
        ),
    ))
}

/// Creates a tangent expression: tan(expr).
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_tan(
    expr: *const Expr
) -> *mut Expr {

    if expr.is_null() {

        return std::ptr::null_mut();
    }

    let expr_ref = &*expr;

    Box::into_raw(Box::new(
        elementary::tan(
            expr_ref.clone(),
        ),
    ))
}

/// Creates an exponential expression: e^(expr).
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_exp(
    expr: *const Expr
) -> *mut Expr {

    if expr.is_null() {

        return std::ptr::null_mut();
    }

    let expr_ref = &*expr;

    Box::into_raw(Box::new(
        elementary::exp(
            expr_ref.clone(),
        ),
    ))
}

/// Creates a natural logarithm expression: ln(expr).
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_ln(
    expr: *const Expr
) -> *mut Expr {

    if expr.is_null() {

        return std::ptr::null_mut();
    }

    let expr_ref = &*expr;

    Box::into_raw(Box::new(
        elementary::ln(
            expr_ref.clone(),
        ),
    ))
}

/// Creates a square root expression: sqrt(expr).
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_sqrt(
    expr: *const Expr
) -> *mut Expr {

    if expr.is_null() {

        return std::ptr::null_mut();
    }

    let expr_ref = &*expr;

    Box::into_raw(Box::new(
        elementary::sqrt(
            expr_ref.clone(),
        ),
    ))
}

/// Creates a power expression: base^exp.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_pow(
    base: *const Expr,
    exp: *const Expr,
) -> *mut Expr {

    if base.is_null() || exp.is_null() {

        return std::ptr::null_mut();
    }

    let base_ref = &*base;

    let exp_ref = &*exp;

    Box::into_raw(Box::new(
        elementary::pow(
            base_ref.clone(),
            exp_ref.clone(),
        ),
    ))
}

/// Returns the symbolic representation of Pi.
#[no_mangle]

pub extern "C" fn rssn_pi() -> *mut Expr
{

    Box::into_raw(Box::new(
        elementary::pi(),
    ))
}

/// Returns the symbolic representation of Euler's number (e).
#[no_mangle]

pub extern "C" fn rssn_e() -> *mut Expr
{

    Box::into_raw(Box::new(
        elementary::e(),
    ))
}

/// Expands a symbolic expression.
///
/// # Safety
/// The caller must ensure `expr` is a valid Expr pointer.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_expand(
    expr: *const Expr
) -> *mut Expr {

    if expr.is_null() {

        return std::ptr::null_mut();
    }

    let expr_ref = &*expr;

    Box::into_raw(Box::new(
        elementary::expand(
            expr_ref.clone(),
        ),
    ))
}

/// Computes binomial coefficient C(n, k).
#[no_mangle]

pub extern "C" fn rssn_binomial_coefficient(
    n: usize,
    k: usize,
) -> *mut c_char {

    let result = elementary::binomial_coefficient(n, k);

    to_c_string(result.to_string())
}

/// Frees an Expr pointer created by this module.
///
/// # Safety
/// The caller must ensure `expr` was created by this module and hasn't been freed yet.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_free_expr(
    expr: *mut Expr
) {

    if !expr.is_null() {

        let _ = Box::from_raw(expr);
    }
}

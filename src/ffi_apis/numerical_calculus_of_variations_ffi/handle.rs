//! Handle-based FFI API for numerical calculus of variations.

use std::ffi::CStr;
use std::os::raw::c_char;
use std::ptr;

use crate::ffi_apis::ffi_api::update_last_error;
use crate::numerical::calculus_of_variations;
use crate::symbolic::core::Expr;

/// Evaluates the action for a given path.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_cov_evaluate_action(
    lagrangian: *const Expr,
    path: *const Expr,
    t_var: *const c_char,
    path_var: *const c_char,
    path_dot_var: *const c_char,
    t_start: f64,
    t_end: f64,
    result: *mut f64,
) -> i32 {

    unsafe {

        if lagrangian.is_null()
            || path.is_null()
            || t_var.is_null()
            || path_var.is_null()
            || path_dot_var.is_null()
            || result.is_null()
        {

            return -1;
        }

        let t_var_str =
            match CStr::from_ptr(t_var)
                .to_str()
            {
                | Ok(s) => s,
                | Err(_) => {

                    update_last_error(
                    "Invalid UTF-8 \
                     for t_var"
                        .to_string(),
                );

                    return -1;
                },
            };

        let path_var_str =
            match CStr::from_ptr(
                path_var,
            )
            .to_str()
            {
                | Ok(s) => s,
                | Err(_) => {

                    update_last_error(
                    "Invalid UTF-8 \
                     for path_var"
                        .to_string(),
                );

                    return -1;
                },
            };

        let path_dot_var_str =
            match CStr::from_ptr(
                path_dot_var,
            )
            .to_str()
            {
                | Ok(s) => s,
                | Err(_) => {

                    update_last_error(
                    "Invalid UTF-8 \
                     for path_dot_var"
                        .to_string(),
                );

                    return -1;
                },
            };

        match calculus_of_variations::evaluate_action(
        &*lagrangian,
        &*path,
        t_var_str,
        path_var_str,
        path_dot_var_str,
        (t_start, t_end),
    ) {
        | Ok(val) => {

            *result = val;

            0
        },
        | Err(e) => {

            update_last_error(e);

            -1
        },
    }
    }
}

/// Computes the Euler-Lagrange expression.
/// Returns a pointer to a new Expr.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_cov_euler_lagrange(
    lagrangian: *const Expr,
    t_var: *const c_char,
    path_var: *const c_char,
    path_dot_var: *const c_char,
) -> *mut Expr {

    unsafe {

        if lagrangian.is_null()
            || t_var.is_null()
            || path_var.is_null()
            || path_dot_var.is_null()
        {

            return ptr::null_mut();
        }

        let t_var_str =
            match CStr::from_ptr(t_var)
                .to_str()
            {
                | Ok(s) => s,
                | Err(_) => {

                    update_last_error(
                    "Invalid UTF-8 \
                     for t_var"
                        .to_string(),
                );

                    return ptr::null_mut();
                },
            };

        let path_var_str =
            match CStr::from_ptr(
                path_var,
            )
            .to_str()
            {
                | Ok(s) => s,
                | Err(_) => {

                    update_last_error(
                    "Invalid UTF-8 \
                     for path_var"
                        .to_string(),
                );

                    return ptr::null_mut();
                },
            };

        let path_dot_var_str =
            match CStr::from_ptr(
                path_dot_var,
            )
            .to_str()
            {
                | Ok(s) => s,
                | Err(_) => {

                    update_last_error(
                    "Invalid UTF-8 \
                     for path_dot_var"
                        .to_string(),
                );

                    return ptr::null_mut();
                },
            };

        let res = calculus_of_variations::euler_lagrange(
        &*lagrangian,
        t_var_str,
        path_var_str,
        path_dot_var_str,
    );

        Box::into_raw(Box::new(res))
    }
}

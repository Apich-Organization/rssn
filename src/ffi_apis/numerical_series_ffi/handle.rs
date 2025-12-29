//! Handle-based FFI API for numerical series operations.

use std::ffi::CStr;
use std::os::raw::c_char;
use std::ptr;

use crate::ffi_apis::ffi_api::update_last_error;
use crate::numerical::series;
use crate::symbolic::core::Expr;

/// Computes the numerical Taylor series coefficients.
/// Returns a pointer to a Vec<f64> containing the coefficients.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_numerical_taylor_coefficients(
    f: *const Expr,
    var: *const c_char,
    at_point: f64,
    order: usize,
) -> *mut Vec<f64> {

    if f.is_null() || var.is_null() {

        return ptr::null_mut();
    }

    let f_expr = &*f;

    let var_str =
        match CStr::from_ptr(var)
            .to_str()
        {
            | Ok(s) => s,
            | Err(_) => {

                update_last_error(
                    "Invalid UTF-8 \
                     string for \
                     variable name"
                        .to_string(),
                );

                return ptr::null_mut();
            },
        };

    match series::taylor_coefficients(
        f_expr,
        var_str,
        at_point,
        order,
    ) {
        | Ok(coeffs) => {
            Box::into_raw(Box::new(
                coeffs,
            ))
        },
        | Err(e) => {

            update_last_error(e);

            ptr::null_mut()
        },
    }
}

/// Evaluates a power series at a point.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_numerical_evaluate_power_series(
    coeffs: *const Vec<f64>,
    at_point: f64,
    x: f64,
) -> f64 {

    if coeffs.is_null() {

        return 0.0;
    }

    series::evaluate_power_series(
        &*coeffs,
        at_point,
        x,
    )
}

/// Computes the sum of a series.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_numerical_sum_series(
    f: *const Expr,
    var: *const c_char,
    start: i64,
    end: i64,
    result: *mut f64,
) -> i32 {

    if f.is_null()
        || var.is_null()
        || result.is_null()
    {

        return -1;
    }

    let f_expr = &*f;

    let var_str =
        match CStr::from_ptr(var)
            .to_str()
        {
            | Ok(s) => s,
            | Err(_) => {

                update_last_error(
                    "Invalid UTF-8 \
                     string for \
                     variable name"
                        .to_string(),
                );

                return -1;
            },
        };

    match series::sum_series(
        f_expr,
        var_str,
        start,
        end,
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

//! Handle-based FFI API for numerical integration.

use crate::ffi_apis::ffi_api::update_last_error;
use crate::numerical::integrate::{self, QuadratureMethod};
use crate::symbolic::core::Expr;
use std::ffi::CStr;
use std::os::raw::c_char;

/// Performs numerical integration (quadrature) of a function.
///
/// # Arguments
/// * `expr_ptr` - Pointer to the `Expr` to integrate.
/// * `var_ptr` - Pointer to the C string representing the variable of integration.
/// * `a` - Lower bound of integration.
/// * `b` - Upper bound of integration.
/// * `n_steps` - Number of steps for non-adaptive methods.
/// * `method` - Integration method:
///     0 - Trapezoidal
///     1 - Simpson
///     2 - Adaptive
///     3 - Romberg
///     4 - Gauss-Legendre
/// * `result` - Pointer to store the result.
///
/// # Returns
/// 0 on success, -1 on error.
#[no_mangle]

pub unsafe extern "C" fn rssn_numerical_quadrature(
    expr_ptr: *const Expr,
    var_ptr: *const c_char,
    a: f64,
    b: f64,
    n_steps: usize,
    method: i32,
    result: *mut f64,
) -> i32 {

    if expr_ptr.is_null() || var_ptr.is_null() || result.is_null() {

        update_last_error("Null pointer passed to rssn_numerical_quadrature".to_string());

        return -1;
    }

    let expr = &*expr_ptr;

    let var_str = match CStr::from_ptr(var_ptr).to_str() {
        Ok(s) => s,
        Err(e) => {

            update_last_error(format!(
                "Invalid UTF-8 in variable name: {}",
                e
            ));

            return -1;
        }
    };

    let q_method = match method {
        0 => QuadratureMethod::Trapezoidal,
        1 => QuadratureMethod::Simpson,
        2 => QuadratureMethod::Adaptive,
        3 => QuadratureMethod::Romberg,
        4 => QuadratureMethod::GaussLegendre,
        _ => {

            update_last_error(format!(
                "Invalid quadrature method: {}",
                method
            ));

            return -1;
        }
    };

    match integrate::quadrature(
        expr,
        var_str,
        (a, b),
        n_steps,
        &q_method,
    ) {
        Ok(val) => {

            *result = val;

            0
        }
        Err(e) => {

            update_last_error(e);

            -1
        }
    }
}

//! Handle-based FFI API for numerical elementary operations.

use crate::ffi_apis::ffi_api::update_last_error;
use crate::numerical::elementary;
use crate::symbolic::core::Expr;
use std::collections::HashMap;
use std::ffi::CStr;
use std::os::raw::c_char;

/// Evaluates an expression handle given variable values.
///
/// # Arguments
/// * `expr_ptr` - Pointer to the Expr object.
/// * `vars` - Array of C-strings for variable names.
/// * `vals` - Array of doubles for variable values.
/// * `num_vars` - Number of variables.
/// * `result` - Pointer to store the result.
///
/// # Returns
/// 0 on success, -1 on failure.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_eval_expr(
    expr_ptr: *const Expr,
    vars: *const *const c_char,
    vals: *const f64,
    num_vars: usize,
    result: *mut f64,
) -> i32 {

    if expr_ptr.is_null()
        || result.is_null()
        || (num_vars > 0 && (vars.is_null() || vals.is_null()))
    {

        update_last_error("Null pointer passed to rssn_num_eval_expr".to_string());

        return -1;
    }

    let mut vars_map = HashMap::new();

    for i in 0..num_vars {

        let name_ptr = unsafe {

            *vars.add(i)
        };

        if name_ptr.is_null() {

            update_last_error(format!("Variable name at index {} is null", i));

            return -1;
        }

        let name = match unsafe {

            CStr::from_ptr(name_ptr).to_str()
        } {
            Ok(s) => s.to_string(),
            Err(e) => {

                update_last_error(format!("Invalid UTF-8 in variable name {}: {}", i, e));

                return -1;
            }
        };

        let val = unsafe {

            *vals.add(i)
        };

        vars_map.insert(name, val);
    }

    match elementary::eval_expr(
        unsafe {

            &*expr_ptr
        },
        &vars_map,
    ) {
        Ok(v) => {

            unsafe {

                *result = v
            };

            0
        }
        Err(e) => {

            update_last_error(e);

            -1
        }
    }
}

/// Pure numerical functions exposed via FFI.
#[no_mangle]

pub extern "C" fn rssn_num_pure_sin(x: f64) -> f64 {

    elementary::pure::sin(x)
}

#[no_mangle]

pub extern "C" fn rssn_num_pure_cos(x: f64) -> f64 {

    elementary::pure::cos(x)
}

#[no_mangle]

pub extern "C" fn rssn_num_pure_tan(x: f64) -> f64 {

    elementary::pure::tan(x)
}

#[no_mangle]

pub extern "C" fn rssn_num_pure_asin(x: f64) -> f64 {

    elementary::pure::asin(x)
}

#[no_mangle]

pub extern "C" fn rssn_num_pure_acos(x: f64) -> f64 {

    elementary::pure::acos(x)
}

#[no_mangle]

pub extern "C" fn rssn_num_pure_atan(x: f64) -> f64 {

    elementary::pure::atan(x)
}

#[no_mangle]

pub extern "C" fn rssn_num_pure_atan2(y: f64, x: f64) -> f64 {

    elementary::pure::atan2(y, x)
}

#[no_mangle]

pub extern "C" fn rssn_num_pure_sinh(x: f64) -> f64 {

    elementary::pure::sinh(x)
}

#[no_mangle]

pub extern "C" fn rssn_num_pure_cosh(x: f64) -> f64 {

    elementary::pure::cosh(x)
}

#[no_mangle]

pub extern "C" fn rssn_num_pure_tanh(x: f64) -> f64 {

    elementary::pure::tanh(x)
}

#[no_mangle]

pub extern "C" fn rssn_num_pure_abs(x: f64) -> f64 {

    elementary::pure::abs(x)
}

#[no_mangle]

pub extern "C" fn rssn_num_pure_sqrt(x: f64) -> f64 {

    elementary::pure::sqrt(x)
}

#[no_mangle]

pub extern "C" fn rssn_num_pure_ln(x: f64) -> f64 {

    elementary::pure::ln(x)
}

#[no_mangle]

pub extern "C" fn rssn_num_pure_exp(x: f64) -> f64 {

    elementary::pure::exp(x)
}

#[no_mangle]

pub extern "C" fn rssn_num_pure_pow(base: f64, exp: f64) -> f64 {

    elementary::pure::pow(base, exp)
}

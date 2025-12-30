//! Handle-based FFI API for numerical elementary operations.

use std::collections::HashMap;
use std::ffi::CStr;
use std::os::raw::c_char;

use crate::ffi_apis::ffi_api::update_last_error;
use crate::numerical::elementary;
use crate::symbolic::core::Expr;

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
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_eval_expr(
    expr_ptr: *const Expr,
    vars: *const *const c_char,
    vals: *const f64,
    num_vars: usize,
    result: *mut f64,
) -> i32 { unsafe {

    if expr_ptr.is_null()
        || result.is_null()
        || (num_vars > 0
            && (vars.is_null()
                || vals.is_null()))
    {

        update_last_error(
            "Null pointer passed to \
             rssn_num_eval_expr"
                .to_string(),
        );

        return -1;
    }

    let mut vars_map = HashMap::new();

    for i in 0 .. num_vars {

        let name_ptr =  {

            *vars.add(i)
        };

        if name_ptr.is_null() {

            update_last_error(format!(
                "Variable name at \
                 index {i} is null"
            ));

            return -1;
        }

        let name = match CStr::from_ptr(name_ptr).to_str()
         {
            | Ok(s) => s.to_string(),
            | Err(e) => {

                update_last_error(format!(
                    "Invalid UTF-8 in variable name {i}: {e}"
                ));

                return -1;
            },
        };

        let val =  {

            *vals.add(i)
        };

        vars_map.insert(name, val);
    }

    match elementary::eval_expr(
         {

            &*expr_ptr
        },
        &vars_map,
    ) {
        | Ok(v) => {

             {

                *result = v;
            };

            0
        },
        | Err(e) => {

            update_last_error(e);

            -1
        },
    }
}}

/// Pure numerical functions exposed via FFI.
/// Computes the sine of a f64 value.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_num_pure_sin(
    x: f64
) -> f64 {

    elementary::pure::sin(x)
}

/// Computes the cosine of a f64 value.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_num_pure_cos(
    x: f64
) -> f64 {

    elementary::pure::cos(x)
}

/// Computes the tangent of a f64 value.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_num_pure_tan(
    x: f64
) -> f64 {

    elementary::pure::tan(x)
}

/// Computes the arcsine of a f64 value.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_num_pure_asin(
    x: f64
) -> f64 {

    elementary::pure::asin(x)
}

/// Computes the arccosine of a f64 value.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_num_pure_acos(
    x: f64
) -> f64 {

    elementary::pure::acos(x)
}

/// Computes the arctangent of a f64 value.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_num_pure_atan(
    x: f64
) -> f64 {

    elementary::pure::atan(x)
}

/// Computes the arctangent of y/x using the signs of the arguments to determine the correct quadrant.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_num_pure_atan2(
    y: f64,

    x: f64,
) -> f64 {

    elementary::pure::atan2(y, x)
}

/// Computes the hyperbolic sine of a f64 value.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_num_pure_sinh(
    x: f64
) -> f64 {

    elementary::pure::sinh(x)
}

/// Computes the hyperbolic cosine of a f64 value.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_num_pure_cosh(
    x: f64
) -> f64 {

    elementary::pure::cosh(x)
}

/// Computes the hyperbolic tangent of a f64 value.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_num_pure_tanh(
    x: f64
) -> f64 {

    elementary::pure::tanh(x)
}

/// Computes the absolute value of a f64 value.

#[unsafe(no_mangle)]

pub const extern "C" fn rssn_num_pure_abs(
    x: f64
) -> f64 {

    elementary::pure::abs(x)
}

/// Computes the square root of a f64 value.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_num_pure_sqrt(
    x: f64
) -> f64 {

    elementary::pure::sqrt(x)
}

/// Computes the natural logarithm of a f64 value.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_num_pure_ln(
    x: f64
) -> f64 {

    elementary::pure::ln(x)
}

/// Computes e raised to the power of a f64 value.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_num_pure_exp(
    x: f64
) -> f64 {

    elementary::pure::exp(x)
}

/// Computes `base` raised to the power of `exp`.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_num_pure_pow(
    base: f64,

    exp: f64,
) -> f64 {

    elementary::pure::pow(base, exp)
}

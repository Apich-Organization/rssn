use std::os::raw::c_char;

use crate::ffi_apis::common::from_json_string;
use crate::ffi_apis::common::to_json_string;
use crate::symbolic::core::Expr;
use crate::symbolic::stats_regression;

// JSON helpers usually take data in custom struct format or just generic List/Vector.
// Since data is [(Expr, Expr)], let's assume JSON input is List of List of 2 or object?
// Usually FFI JSON API takes serialization of input types.
// `&[(Expr, Expr)]` deserializes from `[[x1,y1], [x2,y2], ...]`.

/// Performs a simple linear regression.

///

/// Takes a JSON string representing `Vec<(Expr, Expr)>` (data points).

/// Returns a JSON string representing a `Vec<Expr>` containing the intercept and slope coefficients.

#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_simple_linear_regression(
    data_json: *const c_char
) -> *mut c_char {

    let data: Option<
        Vec<(Expr, Expr)>,
    > = from_json_string(data_json);

    if let Some(data) = data {

        let (b0, b1) = stats_regression::simple_linear_regression_symbolic(&data);

        to_json_string(&vec![b0, b1]) // serialized as [b0, b1]
    } else {

        std::ptr::null_mut()
    }
}

/// Performs a polynomial regression.

///

/// Takes a JSON string representing `Vec<(Expr, Expr)>` (data points) and a `usize` (degree).

/// Returns a JSON string representing a `Vec<Expr>` containing the coefficients of the polynomial.

#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_polynomial_regression(
    data_json: *const c_char,
    degree: usize,
) -> *mut c_char {

    let data: Option<
        Vec<(Expr, Expr)>,
    > = from_json_string(data_json);

    if let Some(data) = data {

        match stats_regression::polynomial_regression_symbolic(&data, degree) {
            | Ok(coeffs) => to_json_string(&coeffs),
            | Err(_) => std::ptr::null_mut(), // or separate error handling
        }
    } else {

        std::ptr::null_mut()
    }
}

/// Performs a nonlinear regression.

///

/// Takes JSON strings representing `Vec<(Expr, Expr)>` (data points), `Expr` (model),

/// `Vec<String>` (variables), and `Vec<String>` (parameters).

/// Returns a JSON string representing `Vec<(Expr, Expr)>` (optimized parameter values).

#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_nonlinear_regression(
    data_json: *const c_char,
    model_json: *const c_char,
    vars_json: *const c_char,
    params_json: *const c_char,
) -> *mut c_char {

    let data: Option<
        Vec<(Expr, Expr)>,
    > = from_json_string(data_json);

    let model: Option<Expr> =
        from_json_string(model_json);

    let vars: Option<Vec<String>> =
        from_json_string(vars_json);

    let params: Option<Vec<String>> =
        from_json_string(params_json);

    match (
        data,
        model,
        vars,
        params,
    ) { (
        Some(data),
        Some(model),
        Some(vars),
        Some(params),
    ) => {

        let vars_refs: Vec<&str> = vars
            .iter()
            .map(std::string::String::as_str)
            .collect();

        let params_refs: Vec<&str> =
            params
                .iter()
                .map(std::string::String::as_str)
                .collect();

        match stats_regression::nonlinear_regression_symbolic(
            &data,
            &model,
            &vars_refs,
            &params_refs,
        ) {
            | Some(solutions) => to_json_string(&solutions), // Vec<(Expr, Expr)>
            | None => std::ptr::null_mut(),
        }
    } _ => {

        std::ptr::null_mut()
    }}
}

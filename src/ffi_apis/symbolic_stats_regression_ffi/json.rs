use std::os::raw::c_char;
use std::os::raw::c_int;
use std::sync::Arc;

use crate::ffi_apis::common::*;
use crate::symbolic::core::Expr;
use crate::symbolic::stats_regression;

// JSON helpers usually take data in custom struct format or just generic List/Vector.
// Since data is [(Expr, Expr)], let's assume JSON input is List of List of 2 or object?
// Usually FFI JSON API takes serialization of input types.
// `&[(Expr, Expr)]` deserializes from `[[x1,y1], [x2,y2], ...]`.

#[no_mangle]

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

#[no_mangle]

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

#[no_mangle]

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

    if let (
        Some(data),
        Some(model),
        Some(vars),
        Some(params),
    ) = (
        data,
        model,
        vars,
        params,
    ) {

        let vars_refs: Vec<&str> = vars
            .iter()
            .map(|s| s.as_str())
            .collect();

        let params_refs: Vec<&str> =
            params
                .iter()
                .map(|s| s.as_str())
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
    } else {

        std::ptr::null_mut()
    }
}

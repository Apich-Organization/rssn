use std::os::raw::c_char;

use crate::ffi_apis::common::*;
use crate::symbolic::core::Expr;
use crate::symbolic::optimize::*;

/// Finds extrema of a function (JSON)
#[no_mangle]

pub extern "C" fn rssn_json_find_extrema(
    expr_json: *const c_char,
    vars_json: *const c_char,
) -> *mut c_char {

    let expr: Option<Expr> =
        from_json_string(expr_json);

    let vars: Option<Vec<String>> =
        from_json_string(vars_json);

    if let (Some(e), Some(v)) =
        (expr, vars)
    {

        let vars_refs: Vec<&str> = v
            .iter()
            .map(|s| s.as_str())
            .collect();

        match find_extrema(
            &e,
            &vars_refs,
        ) {
            | Ok(points) => {
                to_json_string(&points)
            },
            | Err(_) => {
                std::ptr::null_mut()
            },
        }
    } else {

        std::ptr::null_mut()
    }
}

/// Computes Hessian matrix (JSON)
#[no_mangle]

pub extern "C" fn rssn_json_hessian_matrix(
    expr_json: *const c_char,
    vars_json: *const c_char,
) -> *mut c_char {

    let expr: Option<Expr> =
        from_json_string(expr_json);

    let vars: Option<Vec<String>> =
        from_json_string(vars_json);

    if let (Some(e), Some(v)) =
        (expr, vars)
    {

        let vars_refs: Vec<&str> = v
            .iter()
            .map(|s| s.as_str())
            .collect();

        let hessian = hessian_matrix(
            &e,
            &vars_refs,
        );

        to_json_string(&hessian)
    } else {

        std::ptr::null_mut()
    }
}

/// Finds constrained extrema (JSON)
#[no_mangle]

pub extern "C" fn rssn_json_find_constrained_extrema(
    expr_json: *const c_char,
    constraints_json: *const c_char,
    vars_json: *const c_char,
) -> *mut c_char {

    let expr: Option<Expr> =
        from_json_string(expr_json);

    let constraints: Option<Vec<Expr>> =
        from_json_string(
            constraints_json,
        );

    let vars: Option<Vec<String>> =
        from_json_string(vars_json);

    if let (Some(e), Some(c), Some(v)) = (
        expr,
        constraints,
        vars,
    ) {

        let vars_refs: Vec<&str> = v
            .iter()
            .map(|s| s.as_str())
            .collect();

        match find_constrained_extrema(
            &e,
            &c,
            &vars_refs,
        ) {
            | Ok(solutions) => {
                to_json_string(
                    &solutions,
                )
            },
            | Err(_) => {
                std::ptr::null_mut()
            },
        }
    } else {

        std::ptr::null_mut()
    }
}

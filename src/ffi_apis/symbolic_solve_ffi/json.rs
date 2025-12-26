use crate::ffi_apis::common::*;
use crate::symbolic::core::Expr;
use crate::symbolic::solve::{solve, solve_linear_system, solve_system};
use std::ffi::c_char;

#[no_mangle]

pub extern "C" fn rssn_json_solve(
    expr_json: *const c_char,
    var_json: *const c_char,
) -> *mut c_char {

    let expr: Option<Expr> = from_json_string(expr_json);

    let var: Option<String> = from_json_string(var_json);

    if let (Some(e), Some(v)) = (expr, var) {

        let result = solve(&e, &v);

        to_json_string(&result)
    } else {

        std::ptr::null_mut()
    }
}

#[no_mangle]

pub extern "C" fn rssn_json_solve_system(
    equations_json: *const c_char,
    vars_json: *const c_char,
) -> *mut c_char {

    let equations: Option<Vec<Expr>> = from_json_string(equations_json);

    let vars: Option<Vec<String>> = from_json_string(vars_json);

    if let (Some(eqs), Some(vs)) = (equations, vars) {

        let vars_str: Vec<&str> = vs
            .iter()
            .map(|s| s.as_str())
            .collect();

        match solve_system(&eqs, &vars_str) {
            Some(result) => to_json_string(&result),
            None => std::ptr::null_mut(),
        }
    } else {

        std::ptr::null_mut()
    }
}

#[no_mangle]

pub extern "C" fn rssn_json_solve_linear_system(
    system_json: *const c_char,
    vars_json: *const c_char,
) -> *mut c_char {

    let system: Option<Expr> = from_json_string(system_json);

    let vars: Option<Vec<String>> = from_json_string(vars_json);

    if let (Some(sys), Some(vs)) = (system, vars) {

        match solve_linear_system(&sys, &vs) {
            Ok(result) => to_json_string(&result),
            Err(_) => std::ptr::null_mut(),
        }
    } else {

        std::ptr::null_mut()
    }
}

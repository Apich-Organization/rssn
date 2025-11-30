use crate::ffi_apis::common::*;
use crate::symbolic::core::Expr;
use crate::symbolic::number_theory::{
    chinese_remainder, extended_gcd, is_prime, solve_diophantine,
};
use std::ffi::c_char;

#[derive(serde::Deserialize)]
struct Congruence {
    remainder: Expr,
    modulus: Expr,
}

#[no_mangle]
pub extern "C" fn rssn_json_solve_diophantine(
    equation_json: *const c_char,
    vars_json: *const c_char,
) -> *mut c_char {
    let equation: Option<Expr> = from_json_string(equation_json);
    let vars: Option<Vec<String>> = from_json_string(vars_json);

    if let (Some(eq), Some(v)) = (equation, vars) {
        let v_str: Vec<&str> = v.iter().map(|s| s.as_str()).collect();
        match solve_diophantine(&eq, &v_str) {
            Ok(solutions) => to_json_string(&solutions),
            Err(_) => std::ptr::null_mut(),
        }
    } else {
        std::ptr::null_mut()
    }
}

#[no_mangle]
pub extern "C" fn rssn_json_extended_gcd(
    a_json: *const c_char,
    b_json: *const c_char,
) -> *mut c_char {
    let a: Option<Expr> = from_json_string(a_json);
    let b: Option<Expr> = from_json_string(b_json);

    if let (Some(a_expr), Some(b_expr)) = (a, b) {
        let (g, x, y) = extended_gcd(&a_expr, &b_expr);
        to_json_string(&(g, x, y))
    } else {
        std::ptr::null_mut()
    }
}

#[no_mangle]
pub extern "C" fn rssn_json_is_prime(n_json: *const c_char) -> *mut c_char {
    let n: Option<Expr> = from_json_string(n_json);

    if let Some(n_expr) = n {
        let result = is_prime(&n_expr);
        to_json_string(&result)
    } else {
        std::ptr::null_mut()
    }
}

#[no_mangle]
pub extern "C" fn rssn_json_chinese_remainder(congruences_json: *const c_char) -> *mut c_char {
    let congruences_input: Option<Vec<Congruence>> = from_json_string(congruences_json);

    if let Some(input) = congruences_input {
        let congruences: Vec<(Expr, Expr)> = input
            .into_iter()
            .map(|c| (c.remainder, c.modulus))
            .collect();
        match chinese_remainder(&congruences) {
            Some(result) => to_json_string(&result),
            None => std::ptr::null_mut(),
        }
    } else {
        std::ptr::null_mut()
    }
}

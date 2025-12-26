use crate::ffi_apis::common::*;
use crate::symbolic::core::Expr;
use crate::symbolic::logic::{
    is_satisfiable,
    simplify_logic,
    to_cnf,
    to_dnf,
};
use std::ffi::c_char;

/// Simplifies a logical expression using JSON-based FFI.
#[no_mangle]

pub extern "C" fn rssn_json_simplify_logic(expr_json: *const c_char) -> *mut c_char {

    let expr: Option<Expr> = from_json_string(expr_json);

    if let Some(e) = expr {

        let result = simplify_logic(&e);

        to_json_string(&result)
    } else {

        std::ptr::null_mut()
    }
}

/// Converts a logical expression to Conjunctive Normal Form (CNF) using JSON-based FFI.
#[no_mangle]

pub extern "C" fn rssn_json_to_cnf(expr_json: *const c_char) -> *mut c_char {

    let expr: Option<Expr> = from_json_string(expr_json);

    if let Some(e) = expr {

        let result = to_cnf(&e);

        to_json_string(&result)
    } else {

        std::ptr::null_mut()
    }
}

/// Converts a logical expression to Disjunctive Normal Form (DNF) using JSON-based FFI.
#[no_mangle]

pub extern "C" fn rssn_json_to_dnf(expr_json: *const c_char) -> *mut c_char {

    let expr: Option<Expr> = from_json_string(expr_json);

    if let Some(e) = expr {

        let result = to_dnf(&e);

        to_json_string(&result)
    } else {

        std::ptr::null_mut()
    }
}

/// Checks if a logical expression is satisfiable using JSON-based FFI.
///
/// Returns a JSON string containing:
/// - `{"result": "satisfiable"}` if satisfiable
/// - `{"result": "unsatisfiable"}` if unsatisfiable
/// - `{"result": "undecidable"}` if the expression contains quantifiers
#[no_mangle]

pub extern "C" fn rssn_json_is_satisfiable(expr_json: *const c_char) -> *mut c_char {

    let expr: Option<Expr> = from_json_string(expr_json);

    if let Some(e) = expr {

        let result = match is_satisfiable(&e) {
            Some(true) => "satisfiable",
            Some(false) => "unsatisfiable",
            None => "undecidable",
        };

        to_json_string(&result)
    } else {

        std::ptr::null_mut()
    }
}

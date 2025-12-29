//! JSON-based FFI API for symbolic CAS foundations.

use std::os::raw::c_char;

use serde::Deserialize;
use serde::Serialize;

use crate::ffi_apis::common::{from_json_string, to_json_string};
use crate::symbolic::cas_foundations;
use crate::symbolic::core::Expr;
use crate::symbolic::grobner::MonomialOrder;

#[derive(Serialize, Deserialize)]

struct SimplifyWithRelationsInput {
    expr: Expr,
    relations: Vec<Expr>,
    vars: Vec<String>,
    order: MonomialOrder,
}

/// Expands an expression using algebraic rules (JSON).
#[no_mangle]

pub extern "C" fn rssn_cas_expand_json(
    json_str: *const c_char
) -> *mut c_char {

    let expr: Option<Expr> =
        from_json_string(json_str);

    let expr = match expr {
        | Some(e) => e,
        | None => {
            return std::ptr::null_mut()
        },
    };

    let result =
        cas_foundations::expand(expr);

    to_json_string(&result)
}

/// Factorizes an expression (JSON).
#[no_mangle]

pub extern "C" fn rssn_cas_factorize_json(
    json_str: *const c_char
) -> *mut c_char {

    let expr: Option<Expr> =
        from_json_string(json_str);

    let expr = match expr {
        | Some(e) => e,
        | None => {
            return std::ptr::null_mut()
        },
    };

    let result =
        cas_foundations::factorize(
            expr,
        );

    to_json_string(&result)
}

/// Normalizes an expression to a canonical form (JSON).
#[no_mangle]

pub extern "C" fn rssn_cas_normalize_json(
    json_str: *const c_char
) -> *mut c_char {

    let expr: Option<Expr> =
        from_json_string(json_str);

    let expr = match expr {
        | Some(e) => e,
        | None => {
            return std::ptr::null_mut()
        },
    };

    let result =
        cas_foundations::normalize(
            expr,
        );

    to_json_string(&result)
}

/// Simplifies an expression using a set of polynomial side-relations (JSON).
#[no_mangle]

pub extern "C" fn rssn_cas_simplify_with_relations_json(
    json_str: *const c_char
) -> *mut c_char {

    let input: Option<
        SimplifyWithRelationsInput,
    > = from_json_string(json_str);

    let input = match input {
        | Some(i) => i,
        | None => {
            return std::ptr::null_mut()
        },
    };

    let vars_refs: Vec<&str> = input
        .vars
        .iter()
        .map(std::string::String::as_str)
        .collect();

    let result = cas_foundations::simplify_with_relations(
        &input.expr,
        &input.relations,
        &vars_refs,
        input.order,
    );

    to_json_string(&result)
}

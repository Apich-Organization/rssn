//! Bincode-based FFI API for symbolic CAS foundations.

use crate::ffi_apis::common::*;
use crate::symbolic::cas_foundations;
use crate::symbolic::core::Expr;
use crate::symbolic::grobner::MonomialOrder;
use serde::{
    Deserialize,
    Serialize,
};

#[derive(Serialize, Deserialize)]

struct SimplifyWithRelationsInput {
    expr: Expr,
    relations: Vec<Expr>,
    vars: Vec<String>,
    order: MonomialOrder,
}

/// Expands an expression using algebraic rules (Bincode).
#[no_mangle]

pub extern "C" fn rssn_cas_expand_bincode(input: BincodeBuffer) -> BincodeBuffer {

    let expr: Option<Expr> = from_bincode_buffer(&input);

    let expr = match expr {
        | Some(e) => e,
        | None => return BincodeBuffer::empty(),
    };

    let result = cas_foundations::expand(expr);

    to_bincode_buffer(&result)
}

/// Factorizes an expression (Bincode).
#[no_mangle]

pub extern "C" fn rssn_cas_factorize_bincode(input: BincodeBuffer) -> BincodeBuffer {

    let expr: Option<Expr> = from_bincode_buffer(&input);

    let expr = match expr {
        | Some(e) => e,
        | None => return BincodeBuffer::empty(),
    };

    let result = cas_foundations::factorize(expr);

    to_bincode_buffer(&result)
}

/// Normalizes an expression to a canonical form (Bincode).
#[no_mangle]

pub extern "C" fn rssn_cas_normalize_bincode(input: BincodeBuffer) -> BincodeBuffer {

    let expr: Option<Expr> = from_bincode_buffer(&input);

    let expr = match expr {
        | Some(e) => e,
        | None => return BincodeBuffer::empty(),
    };

    let result = cas_foundations::normalize(expr);

    to_bincode_buffer(&result)
}

/// Simplifies an expression using a set of polynomial side-relations (Bincode).
#[no_mangle]

pub extern "C" fn rssn_cas_simplify_with_relations_bincode(input: BincodeBuffer) -> BincodeBuffer {

    let input_data: Option<SimplifyWithRelationsInput> = from_bincode_buffer(&input);

    let input_data = match input_data {
        | Some(i) => i,
        | None => return BincodeBuffer::empty(),
    };

    let vars_refs: Vec<&str> = input_data
        .vars
        .iter()
        .map(|s| s.as_str())
        .collect();

    let result = cas_foundations::simplify_with_relations(
        &input_data.expr,
        &input_data.relations,
        &vars_refs,
        input_data.order,
    );

    to_bincode_buffer(&result)
}

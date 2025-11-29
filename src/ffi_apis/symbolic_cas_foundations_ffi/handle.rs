//! Handle-based FFI API for symbolic CAS foundations.

use crate::ffi_apis::common::*;
use crate::symbolic::cas_foundations;
use crate::symbolic::core::Expr;
use crate::symbolic::grobner::MonomialOrder;
use std::os::raw::c_char;

/// Expands an expression using algebraic rules.
///
/// # Safety
/// The caller must ensure `expr` is a valid Expr pointer.
#[no_mangle]
pub unsafe extern "C" fn rssn_cas_expand(expr: *const Expr) -> *mut Expr {
    if expr.is_null() {
        return std::ptr::null_mut();
    }
    let expr_ref = &*expr;
    Box::into_raw(Box::new(cas_foundations::expand(expr_ref.clone())))
}

/// Factorizes an expression.
///
/// # Safety
/// The caller must ensure `expr` is a valid Expr pointer.
#[no_mangle]
pub unsafe extern "C" fn rssn_cas_factorize(expr: *const Expr) -> *mut Expr {
    if expr.is_null() {
        return std::ptr::null_mut();
    }
    let expr_ref = &*expr;
    Box::into_raw(Box::new(cas_foundations::factorize(expr_ref.clone())))
}

/// Normalizes an expression to a canonical form.
///
/// # Safety
/// The caller must ensure `expr` is a valid Expr pointer.
#[no_mangle]
pub unsafe extern "C" fn rssn_cas_normalize(expr: *const Expr) -> *mut Expr {
    if expr.is_null() {
        return std::ptr::null_mut();
    }
    let expr_ref = &*expr;
    Box::into_raw(Box::new(cas_foundations::normalize(expr_ref.clone())))
}

/// Simplifies an expression using a set of polynomial side-relations.
///
/// # Arguments
/// * `expr` - The expression to simplify.
/// * `relations` - Array of pointers to relation expressions (e.g., `x^2 + y^2 - 1`).
/// * `relations_len` - Number of relations.
/// * `vars` - Array of C strings representing variable ordering.
/// * `vars_len` - Number of variables.
/// * `order_int` - Monomial ordering: 0=Lex, 1=GradedLex, 2=GradedReverseLex.
///
/// # Safety
/// The caller must ensure all pointers are valid.
#[no_mangle]
pub unsafe extern "C" fn rssn_cas_simplify_with_relations(
    expr: *const Expr,
    relations: *const *const Expr,
    relations_len: usize,
    vars: *const *const c_char,
    vars_len: usize,
    order_int: i32,
) -> *mut Expr {
    if expr.is_null()
        || (relations_len > 0 && relations.is_null())
        || (vars_len > 0 && vars.is_null())
    {
        return std::ptr::null_mut();
    }

    let expr_ref = &*expr;

    // Convert relations array
    let mut relations_vec = Vec::with_capacity(relations_len);
    if relations_len > 0 {
        let relations_slice = std::slice::from_raw_parts(relations, relations_len);
        for &rel_ptr in relations_slice {
            if rel_ptr.is_null() {
                return std::ptr::null_mut();
            }
            relations_vec.push((*rel_ptr).clone());
        }
    }

    // Convert vars array
    let mut vars_vec = Vec::with_capacity(vars_len);
    if vars_len > 0 {
        let vars_slice = std::slice::from_raw_parts(vars, vars_len);
        for &var_ptr in vars_slice {
            match c_str_to_str(var_ptr) {
                Some(s) => vars_vec.push(s),
                None => return std::ptr::null_mut(),
            }
        }
    }
    let vars_refs: Vec<&str> = vars_vec.iter().map(|s| *s).collect();

    // Convert order
    let order = match order_int {
        0 => MonomialOrder::Lexicographical,
        1 => MonomialOrder::GradedLexicographical,
        2 => MonomialOrder::GradedReverseLexicographical,
        _ => MonomialOrder::Lexicographical, // Default
    };

    Box::into_raw(Box::new(cas_foundations::simplify_with_relations(
        expr_ref,
        &relations_vec,
        &vars_refs,
        order,
    )))
}

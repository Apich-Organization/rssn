//! Handle-based FFI API for symbolic simplify functions.

use crate::symbolic::core::Expr;
use crate::symbolic::simplify;

/// Simplifies an expression using the heuristic simplifier.
///
/// # Safety
/// The caller must ensure `expr` is a valid Expr pointer.
#[no_mangle]

pub unsafe extern "C" fn rssn_heuristic_simplify(
    expr : *const Expr
) -> *mut Expr {

    if expr.is_null() {

        return std::ptr::null_mut();
    }

    let expr_ref = &*expr;

    #[allow(deprecated)]
    Box::into_raw(Box::new(
        simplify::heuristic_simplify(
            expr_ref.clone(),
        ),
    ))
}

/// Simplifies an expression using the legacy simplifier.
///
/// # Safety
/// The caller must ensure `expr` is a valid Expr pointer.
#[no_mangle]

pub unsafe extern "C" fn rssn_simplify(
    expr : *const Expr
) -> *mut Expr {

    if expr.is_null() {

        return std::ptr::null_mut();
    }

    let expr_ref = &*expr;

    #[allow(deprecated)]
    Box::into_raw(Box::new(
        simplify::simplify(
            expr_ref.clone(),
        ),
    ))
}

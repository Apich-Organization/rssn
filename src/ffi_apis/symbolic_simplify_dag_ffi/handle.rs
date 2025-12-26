//! Handle-based FFI API for symbolic simplify_dag functions.

use crate::symbolic::core::Expr;
use crate::symbolic::simplify_dag;

/// Simplifies an expression using the DAG-based simplifier.
///
/// # Safety
/// The caller must ensure `expr` is a valid Expr pointer.
#[no_mangle]

pub unsafe extern "C" fn rssn_simplify_dag(
    expr: *const Expr
) -> *mut Expr {

    if expr.is_null() {

        return std::ptr::null_mut();
    }

    let expr_ref = &*expr;

    Box::into_raw(Box::new(
        simplify_dag::simplify(
            expr_ref,
        ),
    ))
}

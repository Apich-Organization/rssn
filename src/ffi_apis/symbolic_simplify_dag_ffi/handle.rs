//! Handle-based FFI API for symbolic simplify_dag functions.

use crate::symbolic::core::Expr;
use crate::symbolic::simplify_dag;

/// Simplifies an expression using the DAG-based simplifier.
///
/// # Safety
/// The caller must ensure `expr` is a valid Expr pointer.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_simplify_dag(
    expr: *const Expr
) -> *mut Expr { unsafe {

    if expr.is_null() {

        return std::ptr::null_mut();
    }

    let expr_ref = &*expr;

    Box::into_raw(Box::new(
        simplify_dag::simplify(
            expr_ref,
        ),
    ))
}}

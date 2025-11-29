use crate::symbolic::core::Expr;
use crate::symbolic::logic::{is_satisfiable, simplify_logic, to_cnf, to_dnf};

/// Simplifies a logical expression using handle-based FFI.
///
/// # Safety
/// The caller must ensure that `expr` is a valid pointer to an `Expr`.
#[no_mangle]
pub extern "C" fn rssn_simplify_logic_handle(expr: *const Expr) -> *mut Expr {
    let expr_ref = unsafe {
        if expr.is_null() {
            return std::ptr::null_mut();
        }
        &*expr
    };

    let result = simplify_logic(expr_ref);
    Box::into_raw(Box::new(result))
}

/// Converts a logical expression to Conjunctive Normal Form (CNF) using handle-based FFI.
///
/// # Safety
/// The caller must ensure that `expr` is a valid pointer to an `Expr`.
#[no_mangle]
pub extern "C" fn rssn_to_cnf_handle(expr: *const Expr) -> *mut Expr {
    let expr_ref = unsafe {
        if expr.is_null() {
            return std::ptr::null_mut();
        }
        &*expr
    };

    let result = to_cnf(expr_ref);
    Box::into_raw(Box::new(result))
}

/// Converts a logical expression to Disjunctive Normal Form (DNF) using handle-based FFI.
///
/// # Safety
/// The caller must ensure that `expr` is a valid pointer to an `Expr`.
#[no_mangle]
pub extern "C" fn rssn_to_dnf_handle(expr: *const Expr) -> *mut Expr {
    let expr_ref = unsafe {
        if expr.is_null() {
            return std::ptr::null_mut();
        }
        &*expr
    };

    let result = to_dnf(expr_ref);
    Box::into_raw(Box::new(result))
}

/// Checks if a logical expression is satisfiable using handle-based FFI.
///
/// Returns:
/// - 1 if satisfiable
/// - 0 if unsatisfiable
/// - -1 if the expression contains quantifiers (undecidable)
///
/// # Safety
/// The caller must ensure that `expr` is a valid pointer to an `Expr`.
#[no_mangle]
pub extern "C" fn rssn_is_satisfiable_handle(expr: *const Expr) -> i32 {
    let expr_ref = unsafe {
        if expr.is_null() {
            return -1;
        }
        &*expr
    };

    match is_satisfiable(expr_ref) {
        Some(true) => 1,
        Some(false) => 0,
        None => -1,
    }
}

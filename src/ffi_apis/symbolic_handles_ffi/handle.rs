//! Handle-based FFI API for the HandleManager.
//!
//! This module provides C-compatible functions for managing expression handles.

use std::os::raw::c_char;

use crate::ffi_apis::common::{to_json_string, to_c_string};
use crate::symbolic::core::Expr;
use crate::symbolic::handles::HANDLE_MANAGER;

/// Inserts an expression into the handle manager and returns a unique handle.
///
/// # Safety
/// The caller must ensure `expr` is a valid Expr pointer.
#[no_mangle]

pub unsafe extern "C" fn rssn_handle_insert(
    expr: *const Expr
) -> usize {

    if expr.is_null() {

        return 0; // 0 represents null/invalid handle
    }

    let expr_ref = &*expr;

    HANDLE_MANAGER
        .insert(expr_ref.clone())
}

/// Retrieves an expression from the handle manager.
///
/// Returns a new owned Expr pointer that must be freed by the caller.
///
/// # Safety
/// The caller must ensure the returned pointer is freed using `rssn_free_expr`.
#[no_mangle]

pub unsafe extern "C" fn rssn_handle_get(
    handle: usize
) -> *mut Expr {

    match HANDLE_MANAGER.get(handle) {
        | Some(arc_expr) => {
            Box::into_raw(Box::new(
                (*arc_expr).clone(),
            ))
        },
        | None => std::ptr::null_mut(),
    }
}

/// Checks if a handle exists in the manager.
///
/// # Safety
/// This function is safe to call with any handle value.
#[no_mangle]

pub extern "C" fn rssn_handle_exists(
    handle: usize
) -> bool {

    HANDLE_MANAGER.exists(handle)
}

/// Frees a handle from the manager.
///
/// Returns true if the handle was found and freed, false otherwise.
///
/// # Safety
/// This function is safe to call with any handle value.
#[no_mangle]

pub extern "C" fn rssn_handle_free(
    handle: usize
) -> bool {

    HANDLE_MANAGER
        .free(handle)
        .is_some()
}

/// Returns the number of expressions currently managed.
///
/// # Safety
/// This function is always safe to call.
#[no_mangle]

pub extern "C" fn rssn_handle_count(
) -> usize {

    HANDLE_MANAGER.count()
}

/// Clears all expressions from the handle manager.
///
/// **Warning**: This invalidates all existing handles.
///
/// # Safety
/// This function is always safe to call, but will invalidate all handles.
#[no_mangle]

pub extern "C" fn rssn_handle_clear() {

    HANDLE_MANAGER.clear();
}

/// Returns a list of all active handles as a JSON array string.
///
/// The returned string must be freed using `rssn_free_string`.
///
/// # Safety
/// The caller must free the returned string.
#[no_mangle]

pub extern "C" fn rssn_handle_get_all(
) -> *mut c_char {

    let handles = HANDLE_MANAGER
        .get_all_handles();

    to_json_string(&handles)
}

/// Clones an expression handle, creating a new handle pointing to the same expression.
///
/// Returns 0 if the source handle doesn't exist.
///
/// # Safety
/// This function is safe to call with any handle value.
#[no_mangle]

pub extern "C" fn rssn_handle_clone(
    handle: usize
) -> usize {

    match HANDLE_MANAGER
        .clone_expr(handle)
    {
        | Some(expr) => {
            HANDLE_MANAGER.insert(expr)
        },
        | None => 0,
    }
}

/// Converts an expression handle to a string representation.
///
/// The returned string must be freed using `rssn_free_string`.
///
/// # Safety
/// The caller must free the returned string.
#[no_mangle]

pub extern "C" fn rssn_handle_to_string(
    handle: usize
) -> *mut c_char {

    match HANDLE_MANAGER.get(handle) {
        | Some(arc_expr) => {
            to_c_string(format!(
                "{arc_expr}"
            ))
        },
        | None => std::ptr::null_mut(),
    }
}

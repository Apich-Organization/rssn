//! JSON-based FFI API for compute cache module.

use crate::compute::cache::{ComputationResultCache, ParsingCache};
use crate::ffi_apis::common::{from_json_string, to_c_string};
use crate::symbolic::core::Expr;
use std::ffi::CStr;
use std::os::raw::c_char;
use std::sync::Arc;

// --- ParsingCache ---

/// Retrieves an expression from the ParsingCache as a JSON string.
/// Returns null if not found or error.
#[no_mangle]
pub extern "C" fn rssn_parsing_cache_get_json(
    cache: *mut ParsingCache,
    input: *const c_char,
) -> *mut c_char {
    if cache.is_null() || input.is_null() {
        return std::ptr::null_mut();
    }
    unsafe {
        let input_str = match CStr::from_ptr(input).to_str() {
            Ok(s) => s,
            Err(_) => return std::ptr::null_mut(),
        };

        match (*cache).get(input_str) {
            Some(expr) => match serde_json::to_string(&*expr) {
                Ok(json) => to_c_string(json),
                Err(_) => std::ptr::null_mut(),
            },
            None => std::ptr::null_mut(),
        }
    }
}

/// Stores an expression in the ParsingCache from a JSON string.
#[no_mangle]
pub extern "C" fn rssn_parsing_cache_set_json(
    cache: *mut ParsingCache,
    input: *const c_char,
    json_expr: *const c_char,
) {
    if cache.is_null() || input.is_null() || json_expr.is_null() {
        return;
    }
    unsafe {
        let input_str = match CStr::from_ptr(input).to_str() {
            Ok(s) => s.to_string(),
            Err(_) => return,
        };

        let expr: Option<Expr> = from_json_string(json_expr);
        if let Some(e) = expr {
            (*cache).set(input_str, Arc::new(e));
        }
    }
}

// --- ComputationResultCache ---

/// Retrieves a value from the ComputationResultCache using a JSON expression key.
/// Returns the value as a JSON string (e.g. "\"result\"").
#[no_mangle]
pub extern "C" fn rssn_computation_result_cache_get_json(
    cache: *mut ComputationResultCache,
    json_expr: *const c_char,
) -> *mut c_char {
    if cache.is_null() || json_expr.is_null() {
        return std::ptr::null_mut();
    }
    unsafe {
        let expr: Option<Expr> = from_json_string(json_expr);
        if let Some(e) = expr {
            match (*cache).get(&Arc::new(e)) {
                Some(value) => match serde_json::to_string(&value) {
                    Ok(json) => to_c_string(json),
                    Err(_) => std::ptr::null_mut(),
                },
                None => std::ptr::null_mut(),
            }
        } else {
            std::ptr::null_mut()
        }
    }
}

/// Stores a value in the ComputationResultCache using JSON strings.
#[no_mangle]
pub extern "C" fn rssn_computation_result_cache_set_json(
    cache: *mut ComputationResultCache,
    json_expr: *const c_char,
    json_value: *const c_char,
) {
    if cache.is_null() || json_expr.is_null() || json_value.is_null() {
        return;
    }
    unsafe {
        let expr: Option<Expr> = from_json_string(json_expr);
        let value: Option<String> = from_json_string(json_value);

        if let (Some(e), Some(v)) = (expr, value) {
            (*cache).set(Arc::new(e), v);
        }
    }
}

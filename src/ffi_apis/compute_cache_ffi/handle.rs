//! Handle-based FFI API for compute cache module.

use crate::compute::cache::{ComputationResultCache, ParsingCache};
use crate::compute::computation::Value;
use crate::symbolic::core::Expr;
use std::ffi::{CStr, CString};
use std::os::raw::c_char;
use std::sync::Arc;

// --- ParsingCache ---

/// Creates a new ParsingCache.
/// The caller is responsible for freeing the memory using rssn_parsing_cache_free.
#[no_mangle]
pub extern "C" fn rssn_parsing_cache_new() -> *mut ParsingCache {
    Box::into_raw(Box::new(ParsingCache::new()))
}

/// Frees a ParsingCache.
#[no_mangle]
pub extern "C" fn rssn_parsing_cache_free(cache: *mut ParsingCache) {
    if cache.is_null() {
        return;
    }
    unsafe {
        let _ = Box::from_raw(cache);
    }
}

/// Clears a ParsingCache.
#[no_mangle]
pub extern "C" fn rssn_parsing_cache_clear(cache: *mut ParsingCache) {
    if cache.is_null() {
        return;
    }
    unsafe {
        (*cache).clear();
    }
}

/// Retrieves an expression from the ParsingCache.
/// Returns a pointer to the Expr (Arc<Expr> with incremented refcount), or null if not found.
/// The caller is responsible for freeing the returned Expr (using the appropriate Expr free function).
#[no_mangle]
pub extern "C" fn rssn_parsing_cache_get(
    cache: *mut ParsingCache,
    input: *const c_char,
) -> *mut Expr {
    if cache.is_null() || input.is_null() {
        return std::ptr::null_mut();
    }
    unsafe {
        let input_str = match CStr::from_ptr(input).to_str() {
            Ok(s) => s,
            Err(_) => return std::ptr::null_mut(),
        };
        
        match (*cache).get(input_str) {
            Some(expr) => {
                // Return a raw pointer to the inner Expr, but we need to keep the Arc alive?
                // Actually, usually we pass Arc<Expr> across FFI as *const Expr if it's borrowed,
                // or we need a way to pass ownership.
                // For simplicity here, let's clone the inner Expr to a Box if we want to return a standalone pointer,
                // OR we return a pointer to an Arc<Expr>?
                //
                // If we return *mut Expr, we lose the Arc context unless we re-wrap it.
                // Given Expr is usually used inside Arc, maybe we should return *mut Expr and let the caller wrap it?
                // But Expr::Dag contains Arc<DagNode>.
                //
                // Let's assume for now we return a Box<Expr> that is a clone of the cached Expr.
                // This means we are giving a new copy to the caller.
                Box::into_raw(Box::new((*expr).clone()))
            }
            None => std::ptr::null_mut(),
        }
    }
}

/// Stores an expression in the ParsingCache.
/// The expr pointer is cloned (deep copy of the structure, but DAG nodes are shared).
#[no_mangle]
pub extern "C" fn rssn_parsing_cache_set(
    cache: *mut ParsingCache,
    input: *const c_char,
    expr: *const Expr,
) {
    if cache.is_null() || input.is_null() || expr.is_null() {
        return;
    }
    unsafe {
        let input_str = match CStr::from_ptr(input).to_str() {
            Ok(s) => s.to_string(),
            Err(_) => return,
        };
        
        let expr_arc = Arc::new((*expr).clone());
        (*cache).set(input_str, expr_arc);
    }
}

// --- ComputationResultCache ---

/// Creates a new ComputationResultCache.
/// The caller is responsible for freeing the memory using rssn_computation_result_cache_free.
#[no_mangle]
pub extern "C" fn rssn_computation_result_cache_new() -> *mut ComputationResultCache {
    Box::into_raw(Box::new(ComputationResultCache::new()))
}

/// Frees a ComputationResultCache.
#[no_mangle]
pub extern "C" fn rssn_computation_result_cache_free(cache: *mut ComputationResultCache) {
    if cache.is_null() {
        return;
    }
    unsafe {
        let _ = Box::from_raw(cache);
    }
}

/// Clears a ComputationResultCache.
#[no_mangle]
pub extern "C" fn rssn_computation_result_cache_clear(cache: *mut ComputationResultCache) {
    if cache.is_null() {
        return;
    }
    unsafe {
        (*cache).clear();
    }
}

/// Retrieves a value from the ComputationResultCache.
/// Returns a C string (char*) which must be freed by the caller using rssn_free_string.
/// Returns null if not found.
#[no_mangle]
pub extern "C" fn rssn_computation_result_cache_get(
    cache: *mut ComputationResultCache,
    expr: *const Expr,
) -> *mut c_char {
    if cache.is_null() || expr.is_null() {
        return std::ptr::null_mut();
    }
    unsafe {
        // We need to construct an Arc<Expr> to query the cache, but we only have a raw pointer.
        // The cache key is Arc<Expr>.
        // We can create a temporary Arc for the lookup if we clone the Expr.
        let expr_arc = Arc::new((*expr).clone());
        
        match (*cache).get(&expr_arc) {
            Some(value) => {
                match CString::new(value) {
                    Ok(c_str) => c_str.into_raw(),
                    Err(_) => std::ptr::null_mut(),
                }
            }
            None => std::ptr::null_mut(),
        }
    }
}

/// Stores a value in the ComputationResultCache.
#[no_mangle]
pub extern "C" fn rssn_computation_result_cache_set(
    cache: *mut ComputationResultCache,
    expr: *const Expr,
    value: *const c_char,
) {
    if cache.is_null() || expr.is_null() || value.is_null() {
        return;
    }
    unsafe {
        let value_str = match CStr::from_ptr(value).to_str() {
            Ok(s) => s.to_string(),
            Err(_) => return,
        };
        
        let expr_arc = Arc::new((*expr).clone());
        (*cache).set(expr_arc, value_str);
    }
}

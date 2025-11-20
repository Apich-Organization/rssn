use crate::compute::cache::{ParsingCache, ComputationResultCache};
use crate::symbolic::core::Expr;
use crate::compute::computation::Value;
use std::sync::Arc;
use std::ffi::{CString, CStr};
use std::os::raw::{c_char, c_void};

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

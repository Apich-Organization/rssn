//! Bincode-based FFI API for compute cache module.

use crate::compute::cache::{ComputationResultCache, ParsingCache};
use crate::symbolic::core::Expr;
use crate::ffi_apis::common::{BincodeBuffer, to_bincode_buffer, from_bincode_buffer};
use std::ffi::CStr;
use std::os::raw::c_char;
use std::sync::Arc;

// --- ParsingCache ---

/// Retrieves an expression from the ParsingCache as a bincode buffer.
#[no_mangle]
pub extern "C" fn rssn_parsing_cache_get_bincode(
    cache: *mut ParsingCache,
    input: *const c_char,
) -> BincodeBuffer {
    if cache.is_null() || input.is_null() {
        return BincodeBuffer::empty();
    }
    unsafe {
        let input_str = match CStr::from_ptr(input).to_str() {
            Ok(s) => s,
            Err(_) => return BincodeBuffer::empty(),
        };
        
        match (*cache).get(input_str) {
            Some(expr) => to_bincode_buffer(&*expr),
            None => BincodeBuffer::empty(),
        }
    }
}

/// Stores an expression in the ParsingCache from a bincode buffer.
#[no_mangle]
pub extern "C" fn rssn_parsing_cache_set_bincode(
    cache: *mut ParsingCache,
    input: *const c_char,
    buffer: BincodeBuffer,
) {
    if cache.is_null() || input.is_null() {
        return;
    }
    unsafe {
        let input_str = match CStr::from_ptr(input).to_str() {
            Ok(s) => s.to_string(),
            Err(_) => return,
        };
        
        let expr: Option<Expr> = from_bincode_buffer(&buffer);
        if let Some(e) = expr {
            (*cache).set(input_str, Arc::new(e));
        }
    }
}

// --- ComputationResultCache ---

/// Retrieves a value from the ComputationResultCache using a bincode expression key.
#[no_mangle]
pub extern "C" fn rssn_computation_result_cache_get_bincode(
    cache: *mut ComputationResultCache,
    expr_buffer: BincodeBuffer,
) -> BincodeBuffer {
    if cache.is_null() {
        return BincodeBuffer::empty();
    }
    unsafe {
        let expr: Option<Expr> = from_bincode_buffer(&expr_buffer);
        if let Some(e) = expr {
            match (*cache).get(&Arc::new(e)) {
                Some(value) => to_bincode_buffer(&value),
                None => BincodeBuffer::empty(),
            }
        } else {
            BincodeBuffer::empty()
        }
    }
}

/// Stores a value in the ComputationResultCache using bincode buffers.
#[no_mangle]
pub extern "C" fn rssn_computation_result_cache_set_bincode(
    cache: *mut ComputationResultCache,
    expr_buffer: BincodeBuffer,
    value_buffer: BincodeBuffer,
) {
    if cache.is_null() {
        return;
    }
    unsafe {
        let expr: Option<Expr> = from_bincode_buffer(&expr_buffer);
        let value: Option<String> = from_bincode_buffer(&value_buffer);
        
        if let (Some(e), Some(v)) = (expr, value) {
            (*cache).set(Arc::new(e), v);
        }
    }
}

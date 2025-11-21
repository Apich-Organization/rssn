//! Bincode-based FFI API for symbolic simplify functions.

use crate::ffi_apis::common::*;
use crate::symbolic::core::Expr;
use crate::symbolic::simplify;

/// Simplifies an expression using the heuristic simplifier (Bincode input/output).
#[no_mangle]
pub extern "C" fn rssn_bincode_heuristic_simplify(expr_buf: BincodeBuffer) -> BincodeBuffer {
    let expr: Option<Expr> = from_bincode_buffer(&expr_buf);
    
    if let Some(e) = expr {
        #[allow(deprecated)]
        let result = simplify::heuristic_simplify(e);
        to_bincode_buffer(&result)
    } else {
        BincodeBuffer::empty()
    }
}

/// Simplifies an expression using the legacy simplifier (Bincode input/output).
#[no_mangle]
pub extern "C" fn rssn_bincode_simplify(expr_buf: BincodeBuffer) -> BincodeBuffer {
    let expr: Option<Expr> = from_bincode_buffer(&expr_buf);
    
    if let Some(e) = expr {
        #[allow(deprecated)]
        let result = simplify::simplify(e);
        to_bincode_buffer(&result)
    } else {
        BincodeBuffer::empty()
    }
}

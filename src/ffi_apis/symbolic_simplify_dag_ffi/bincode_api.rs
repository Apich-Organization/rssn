//! Bincode-based FFI API for symbolic simplify_dag functions.

use crate::ffi_apis::common::*;
use crate::symbolic::core::Expr;
use crate::symbolic::simplify_dag;

/// Simplifies an expression using the DAG-based simplifier (Bincode input/output).
#[no_mangle]
pub extern "C" fn rssn_bincode_simplify_dag(expr_buf: BincodeBuffer) -> BincodeBuffer {
    let expr: Option<Expr> = from_bincode_buffer(&expr_buf);
    
    if let Some(e) = expr {
        let result = simplify_dag::simplify(&e);
        to_bincode_buffer(&result)
    } else {
        BincodeBuffer::empty()
    }
}

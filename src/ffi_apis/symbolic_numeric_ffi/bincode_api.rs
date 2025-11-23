use crate::ffi_apis::common::*;
use crate::symbolic::core::Expr;
use crate::symbolic::numeric::evaluate_numerical;

#[no_mangle]
pub extern "C" fn rssn_bincode_evaluate_numerical(expr_buf: BincodeBuffer) -> BincodeBuffer {
    let expr: Option<Expr> = from_bincode_buffer(&expr_buf);
    
    if let Some(e) = expr {
        if let Some(result) = evaluate_numerical(&e) {
            to_bincode_buffer(&result)
        } else {
            BincodeBuffer::empty()
        }
    } else {
        BincodeBuffer::empty()
    }
}

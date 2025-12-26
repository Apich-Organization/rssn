use crate::ffi_apis::common::*;
use crate::symbolic::core::Expr;
use crate::symbolic::unit_unification::unify_expression;

#[no_mangle]

pub extern "C" fn rssn_bincode_unify_expression(expr_buf : BincodeBuffer) -> BincodeBuffer {

    let expr : Option<Expr> = from_bincode_buffer(&expr_buf);

    if let Some(e) = expr {

        match unify_expression(&e) {
            | Ok(result) => to_bincode_buffer(&result),
            | Err(_) => BincodeBuffer::empty(),
        }
    } else {

        BincodeBuffer::empty()
    }
}

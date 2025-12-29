use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::ffi_apis::common::BincodeBuffer;
use crate::symbolic::core::Expr;
use crate::symbolic::numeric::evaluate_numerical;

/// Numerically evaluates a symbolic expression.

///

/// Takes a bincode-serialized `Expr` as input,

/// and returns a bincode-serialized numerical evaluation of that expression.

#[no_mangle]

pub extern "C" fn rssn_bincode_evaluate_numerical(
    expr_buf: BincodeBuffer
) -> BincodeBuffer {

    let expr: Option<Expr> =
        from_bincode_buffer(&expr_buf);

    if let Some(e) = expr {

        if let Some(result) =
            evaluate_numerical(&e)
        {

            to_bincode_buffer(&result)
        } else {

            BincodeBuffer::empty()
        }
    } else {

        BincodeBuffer::empty()
    }
}

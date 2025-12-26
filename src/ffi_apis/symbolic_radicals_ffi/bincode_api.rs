use crate::ffi_apis::common::*;
use crate::symbolic::core::Expr;
use crate::symbolic::radicals::*;

/// Simplifies radical expressions (Bincode)
#[no_mangle]

pub extern "C" fn rssn_bincode_simplify_radicals(
    expr_buf: BincodeBuffer
) -> BincodeBuffer {

    let expr: Option<Expr> =
        from_bincode_buffer(&expr_buf);

    if let Some(e) = expr {

        let result =
            simplify_radicals(&e);

        to_bincode_buffer(&result)
    } else {

        BincodeBuffer::empty()
    }
}

/// Denests a nested square root (Bincode)
#[no_mangle]

pub extern "C" fn rssn_bincode_denest_sqrt(
    expr_buf: BincodeBuffer
) -> BincodeBuffer {

    let expr: Option<Expr> =
        from_bincode_buffer(&expr_buf);

    if let Some(e) = expr {

        let result = denest_sqrt(&e);

        to_bincode_buffer(&result)
    } else {

        BincodeBuffer::empty()
    }
}

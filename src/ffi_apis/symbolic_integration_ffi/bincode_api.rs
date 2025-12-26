use crate::ffi_apis::common::*;
use crate::symbolic::core::Expr;
use crate::symbolic::integration::*;

/// Integrates an expression using the Risch-Norman algorithm (Bincode)
#[no_mangle]

pub extern "C" fn rssn_bincode_risch_norman_integrate(
    expr_buf: BincodeBuffer,
    x_buf: BincodeBuffer,
) -> BincodeBuffer {

    let expr: Option<Expr> = from_bincode_buffer(&expr_buf);

    let x: Option<String> = from_bincode_buffer(&x_buf);

    if let (Some(e), Some(var)) = (expr, x) {

        let result = risch_norman_integrate(&e, &var);

        to_bincode_buffer(&result)
    } else {

        BincodeBuffer::empty()
    }
}

/// Integrates a rational function (Bincode)
#[no_mangle]

pub extern "C" fn rssn_bincode_integrate_rational_function(
    expr_buf: BincodeBuffer,
    x_buf: BincodeBuffer,
) -> BincodeBuffer {

    let expr: Option<Expr> = from_bincode_buffer(&expr_buf);

    let x: Option<String> = from_bincode_buffer(&x_buf);

    if let (Some(e), Some(var)) = (expr, x) {

        match integrate_rational_function_expr(&e, &var) {
            Ok(result) => to_bincode_buffer(&result),
            Err(_) => BincodeBuffer::empty(),
        }
    } else {

        BincodeBuffer::empty()
    }
}

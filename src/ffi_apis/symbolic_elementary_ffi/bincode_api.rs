//! Bincode-based FFI API for symbolic elementary functions.

use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::ffi_apis::common::BincodeBuffer;
use crate::symbolic::core::Expr;
use crate::symbolic::elementary;

/// Creates a sine expression from bincode: sin(expr).
///
/// # Arguments
/// * `expr_buffer` - Bincode-serialized Expr
///
/// # Returns
/// Bincode-serialized Expr
#[no_mangle]

pub extern "C" fn rssn_sin_bincode(
    expr_buffer: BincodeBuffer
) -> BincodeBuffer {

    let expr : Expr = match from_bincode_buffer(&expr_buffer) {
        | Some(e) => e,
        | None => return BincodeBuffer::empty(),
    };

    to_bincode_buffer(&elementary::sin(
        expr,
    ))
}

/// Creates a cosine expression from bincode: cos(expr).
#[no_mangle]

pub extern "C" fn rssn_cos_bincode(
    expr_buffer: BincodeBuffer
) -> BincodeBuffer {

    let expr : Expr = match from_bincode_buffer(&expr_buffer) {
        | Some(e) => e,
        | None => return BincodeBuffer::empty(),
    };

    to_bincode_buffer(&elementary::cos(
        expr,
    ))
}

/// Creates a tangent expression from bincode: tan(expr).
#[no_mangle]

pub extern "C" fn rssn_tan_bincode(
    expr_buffer: BincodeBuffer
) -> BincodeBuffer {

    let expr : Expr = match from_bincode_buffer(&expr_buffer) {
        | Some(e) => e,
        | None => return BincodeBuffer::empty(),
    };

    to_bincode_buffer(&elementary::tan(
        expr,
    ))
}

/// Creates an exponential expression from bincode: e^(expr).
#[no_mangle]

pub extern "C" fn rssn_exp_bincode(
    expr_buffer: BincodeBuffer
) -> BincodeBuffer {

    let expr : Expr = match from_bincode_buffer(&expr_buffer) {
        | Some(e) => e,
        | None => return BincodeBuffer::empty(),
    };

    to_bincode_buffer(&elementary::exp(
        expr,
    ))
}

/// Creates a natural logarithm expression from bincode: ln(expr).
#[no_mangle]

pub extern "C" fn rssn_ln_bincode(
    expr_buffer: BincodeBuffer
) -> BincodeBuffer {

    let expr : Expr = match from_bincode_buffer(&expr_buffer) {
        | Some(e) => e,
        | None => return BincodeBuffer::empty(),
    };

    to_bincode_buffer(&elementary::ln(
        expr,
    ))
}

/// Creates a square root expression from bincode: sqrt(expr).
#[no_mangle]

pub extern "C" fn rssn_sqrt_bincode(
    expr_buffer: BincodeBuffer
) -> BincodeBuffer {

    let expr : Expr = match from_bincode_buffer(&expr_buffer) {
        | Some(e) => e,
        | None => return BincodeBuffer::empty(),
    };

    to_bincode_buffer(
        &elementary::sqrt(expr),
    )
}

/// Creates a power expression from bincode: base^exp.
///
/// # Arguments
/// * `base_buffer` - Bincode-serialized base Expr
/// * `exp_buffer` - Bincode-serialized exponent Expr
#[no_mangle]

pub extern "C" fn rssn_pow_bincode(
    base_buffer: BincodeBuffer,
    exp_buffer: BincodeBuffer,
) -> BincodeBuffer {

    let base : Expr = match from_bincode_buffer(&base_buffer) {
        | Some(e) => e,
        | None => return BincodeBuffer::empty(),
    };

    let exp : Expr = match from_bincode_buffer(&exp_buffer) {
        | Some(e) => e,
        | None => return BincodeBuffer::empty(),
    };

    to_bincode_buffer(&elementary::pow(
        base, exp,
    ))
}

/// Returns Pi as bincode.
#[no_mangle]

pub extern "C" fn rssn_pi_bincode(
) -> BincodeBuffer {

    to_bincode_buffer(&elementary::pi())
}

/// Returns Euler's number (e) as bincode.
#[no_mangle]

pub extern "C" fn rssn_e_bincode(
) -> BincodeBuffer {

    to_bincode_buffer(&elementary::e())
}

/// Expands a symbolic expression from bincode.
#[no_mangle]

pub extern "C" fn rssn_expand_bincode(
    expr_buffer: BincodeBuffer
) -> BincodeBuffer {

    let expr : Expr = match from_bincode_buffer(&expr_buffer) {
        | Some(e) => e,
        | None => return BincodeBuffer::empty(),
    };

    to_bincode_buffer(
        &elementary::expand(expr),
    )
}

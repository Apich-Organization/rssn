//! Bincode-based FFI API for numerical elementary operations.

use crate::ffi_apis::common::BincodeBuffer;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::elementary;
use crate::symbolic::core::Expr;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

#[derive(Deserialize)]

struct EvalRequest {
    expr: Expr,
    vars: HashMap<String, f64>,
}

fn decode<T: for<'de> Deserialize<'de>>(
    data: *const u8,
    len: usize,
) -> Option<T> {

    if data.is_null() {

        return None;
    }

    let slice = unsafe {

        std::slice::from_raw_parts(data, len)
    };

    bincode_next::serde::decode_from_slice(
        slice,
        bincode_next::config::standard(),
    )
    .ok()
    .map(|(v, _)| v)
}

fn encode<T: Serialize>(val: &T) -> BincodeBuffer {

    match bincode_next::serde::encode_to_vec(
        val,
        bincode_next::config::standard(),
    ) {
        Ok(bytes) => BincodeBuffer::from_vec(bytes),
        Err(_) => BincodeBuffer::empty(),
    }
}

/// Evaluates an expression from a Bincode buffer.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_eval_bincode(
    data: *const u8,
    len: usize,
) -> BincodeBuffer {

    let req: EvalRequest = match decode(data, len) {
        Some(r) => r,
        None => {

            let res: FfiResult<f64, String> = FfiResult {
                ok: None,
                err: Some("Bincode decode error".to_string()),
            };

            return encode(&res);
        }
    };

    let result = elementary::eval_expr(&req.expr, &req.vars);

    let ffi_res = match result {
        Ok(v) => FfiResult {
            ok: Some(v),
            err: None,
        },
        Err(e) => FfiResult {
            ok: None,
            err: Some(e),
        },
    };

    encode(&ffi_res)
}

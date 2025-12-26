//! Bincode-based FFI API for numerical integration.

use serde::Deserialize;

use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::ffi_apis::common::BincodeBuffer;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::integrate::QuadratureMethod;
use crate::numerical::integrate::{
    self,
};
use crate::symbolic::core::Expr;

#[derive(Deserialize)]

struct QuadratureInput {
    expr: Expr,
    var: String,
    a: f64,
    b: f64,
    n_steps: usize,
    method: QuadratureMethod,
}

/// Performs numerical integration via Bincode.
#[no_mangle]

pub unsafe extern "C" fn rssn_numerical_quadrature_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input: QuadratureInput =
        match from_bincode_buffer(
            &buffer,
        ) {
            | Some(v) => v,
            | None => {

                let res : FfiResult<f64, String> = FfiResult {
                ok : None,
                err : Some("Bincode decoding error".to_string()),
            };

                return to_bincode_buffer(&res);
            },
        };

    let result = integrate::quadrature(
        &input.expr,
        &input.var,
        (input.a, input.b),
        input.n_steps,
        &input.method,
    );

    let res = match result {
        | Ok(val) => {
            FfiResult {
                ok: Some(val),
                err: None::<String>,
            }
        },
        | Err(e) => {
            FfiResult {
                ok: None,
                err: Some(e),
            }
        },
    };

    to_bincode_buffer(&res)
}

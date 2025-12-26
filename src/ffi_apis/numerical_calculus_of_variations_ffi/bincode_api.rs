//! Bincode-based FFI API for numerical calculus of variations.

use serde::Deserialize;
use serde::Serialize;

use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::ffi_apis::common::BincodeBuffer;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::calculus_of_variations;
use crate::symbolic::core::Expr;

#[derive(Deserialize)]

struct ActionInput {
    lagrangian : Expr,
    path : Expr,
    t_var : String,
    path_var : String,
    path_dot_var : String,
    t_range : (f64, f64),
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_cov_evaluate_action_bincode(
    buffer : BincodeBuffer
) -> BincodeBuffer {

    let input : ActionInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(
                &FfiResult::<f64, String> {
                    ok : None,
                    err : Some("Invalid Bincode input".to_string()),
                },
            )
        },
    };

    match calculus_of_variations::evaluate_action(
        &input.lagrangian,
        &input.path,
        &input.t_var,
        &input.path_var,
        &input.path_dot_var,
        input.t_range,
    ) {
        | Ok(val) => {

            let ffi_res = FfiResult {
                ok : Some(val),
                err : None::<String>,
            };

            to_bincode_buffer(&ffi_res)
        },
        | Err(e) => {

            let ffi_res = FfiResult {
                ok : None::<f64>,
                err : Some(e),
            };

            to_bincode_buffer(&ffi_res)
        },
    }
}

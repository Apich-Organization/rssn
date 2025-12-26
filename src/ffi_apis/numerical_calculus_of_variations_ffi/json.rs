//! JSON-based FFI API for numerical calculus of variations.

use crate::ffi_apis::common::{
    from_json_string,
    to_c_string,
};
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::calculus_of_variations;
use crate::symbolic::core::Expr;
use serde::{
    Deserialize,
    Serialize,
};
use std::os::raw::c_char;

#[derive(Deserialize)]

struct ActionInput {
    lagrangian: Expr,
    path: Expr,
    t_var: String,
    path_var: String,
    path_dot_var: String,
    t_range: (f64, f64),
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_cov_evaluate_action_json(
    input_json: *const c_char
) -> *mut c_char {

    let input: ActionInput = match from_json_string(input_json) {
        Some(i) => i,
        None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok: None,
                        err: Some("Invalid JSON input".to_string()),
                    },
                )
                .unwrap(),
            )
        }
    };

    match calculus_of_variations::evaluate_action(
        &input.lagrangian,
        &input.path,
        &input.t_var,
        &input.path_var,
        &input.path_dot_var,
        input.t_range,
    ) {
        Ok(val) => {

            let ffi_res = FfiResult {
                ok: Some(val),
                err: None::<String>,
            };

            to_c_string(serde_json::to_string(&ffi_res).unwrap())
        }
        Err(e) => {

            let ffi_res = FfiResult {
                ok: None::<f64>,
                err: Some(e),
            };

            to_c_string(serde_json::to_string(&ffi_res).unwrap())
        }
    }
}

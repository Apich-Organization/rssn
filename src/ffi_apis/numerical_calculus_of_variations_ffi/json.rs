//! JSON-based FFI API for numerical calculus of variations.

use std::os::raw::c_char;

use serde::Deserialize;

use crate::ffi_apis::common::from_json_string;
use crate::ffi_apis::common::to_c_string;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::calculus_of_variations;
use crate::symbolic::core::Expr;

#[derive(Deserialize)]

struct ActionInput {
    lagrangian: Expr,
    path: Expr,
    t_var: String,
    path_var: String,
    path_dot_var: String,
    t_range: (f64, f64),
}

/// Evaluates the action for a given Lagrangian and path using JSON for serialization.

#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_cov_evaluate_action_json(
    input_json: *const c_char
) -> *mut c_char {

    let input : ActionInput = match from_json_string(input_json) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok : None,
                        err : Some("Invalid JSON input".to_string()),
                    },
                )
                .unwrap(),
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

            to_c_string(serde_json::to_string(&ffi_res).unwrap())
        },
        | Err(e) => {

            let ffi_res = FfiResult {
                ok : None::<f64>,
                err : Some(e),
            };

            to_c_string(serde_json::to_string(&ffi_res).unwrap())
        },
    }
}

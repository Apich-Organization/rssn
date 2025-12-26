//! JSON-based FFI API for numerical ODE solvers.

use crate::ffi_apis::common::{from_json_string, to_c_string};
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::ode::{self, OdeSolverMethod};
use crate::symbolic::core::Expr;
use serde::{Deserialize, Serialize};
use std::os::raw::c_char;

#[derive(Deserialize)]

struct OdeInput {
    funcs: Vec<Expr>,
    y0: Vec<f64>,
    x_range: (f64, f64),
    num_steps: usize,
    method: OdeSolverMethod,
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_ode_solve_json(input_json: *const c_char) -> *mut c_char {

    let input: OdeInput = match from_json_string(input_json) {
        Some(i) => i,
        None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<Vec<Vec<f64>>, String> {
                        ok: None,
                        err: Some("Invalid JSON input".to_string()),
                    },
                )
                .unwrap(),
            )
        }
    };

    let res = ode::solve_ode_system(
        &input.funcs,
        &input.y0,
        input.x_range,
        input.num_steps,
        input.method,
    );

    let ffi_res = match res {
        Ok(v) => FfiResult {
            ok: Some(v),
            err: None,
        },
        Err(e) => FfiResult {
            ok: None,
            err: Some(e),
        },
    };

    to_c_string(serde_json::to_string(&ffi_res).unwrap())
}

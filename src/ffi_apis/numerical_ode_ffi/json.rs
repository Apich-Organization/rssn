//! JSON-based FFI API for numerical ODE solvers.

use std::os::raw::c_char;

use serde::Deserialize;

use crate::ffi_apis::common::from_json_string;
use crate::ffi_apis::common::to_c_string;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::ode::OdeSolverMethod;
use crate::numerical::ode::{
    self,
};
use crate::symbolic::core::Expr;

#[derive(Deserialize)]

struct OdeInput {
    funcs: Vec<Expr>,
    y0: Vec<f64>,
    x_range: (f64, f64),
    num_steps: usize,
    method: OdeSolverMethod,
}

/// Solves a system of ordinary differential equations (ODEs) using JSON serialization.
///
/// # Arguments
///
/// * `input_json` - A JSON string pointer containing:
///   - `funcs`: Array of differential equations as `Expr`
///   - `y0`: Initial values for each function
///   - `x_range`: Integration range as [start, end]
///   - `num_steps`: Number of steps for the solver
///   - `method`: ODE solver method to use
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<Vec<Vec<f64>>, String>` with either:
/// - `ok`: Solution matrix where each row represents a step
/// - `err`: Error message if solving failed
///
/// # Safety
///
/// This function is unsafe because it:
/// - Receives a raw C string pointer that must be valid, null-terminated UTF-8
/// - Returns a raw pointer that the caller must free using `rssn_free_string`
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_ode_solve_json(
    input_json: *const c_char
) -> *mut c_char {

    let input : OdeInput = match from_json_string(input_json) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<Vec<Vec<f64>>, String> {
                        ok : None,
                        err : Some("Invalid JSON input".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let res = ode::solve_ode_system(
        &input.funcs,
        &input.y0,
        input.x_range,
        input.num_steps,
        input.method,
    );

    let ffi_res = match res {
        | Ok(v) => {
            FfiResult {
                ok: Some(v),
                err: None,
            }
        },
        | Err(e) => {
            FfiResult {
                ok: None,
                err: Some(e),
            }
        },
    };

    to_c_string(
        serde_json::to_string(&ffi_res)
            .unwrap(),
    )
}

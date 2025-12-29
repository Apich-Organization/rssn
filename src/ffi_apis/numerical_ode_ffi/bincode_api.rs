//! Bincode-based FFI API for numerical ODE solvers.

use serde::Deserialize;

use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::ffi_apis::common::BincodeBuffer;
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

/// Solves a system of ordinary differential equations (ODEs) using bincode serialization.
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing `OdeInput` with:
///   - `funcs`: Vector of differential equations as `Expr`
///   - `y0`: Initial values for each function
///   - `x_range`: Integration range as (start, end)
///   - `num_steps`: Number of steps for the solver
///   - `method`: ODE solver method to use
///
/// # Returns
///
/// A bincode-encoded buffer containing `FfiResult<Vec<Vec<f64>>, String>` with either:
/// - `ok`: Solution matrix where each row represents a step
/// - `err`: Error message if solving failed
///
/// # Safety
///
/// This function is unsafe because it receives raw pointers through FFI.
/// The caller must ensure the input buffer is valid bincode data.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_ode_solve_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input : OdeInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(
                &FfiResult::<Vec<Vec<f64>>, String> {
                    ok : None,
                    err : Some("Invalid Bincode input".to_string()),
                },
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

    to_bincode_buffer(&ffi_res)
}

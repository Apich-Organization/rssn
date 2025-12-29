//! JSON-based FFI API for physics CNM functions.

use std::os::raw::c_char;

use serde::Deserialize;

use crate::ffi_apis::common::from_json_string;
use crate::ffi_apis::common::to_c_string;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::physics::physics_cnm::HeatEquationSolverConfig;
use crate::physics::physics_cnm::{
    self,
};

#[derive(Deserialize)]

struct Heat2DInput {
    initial_condition: Vec<f64>,
    config: HeatEquationSolverConfig,
}

/// Solves the 2D heat equation using Crank-Nicolson ADI method via JSON serialization.
///
/// The heat equation ∂u/∂t = α∇²u is solved using the Crank-Nicolson Alternating
/// Direction Implicit (ADI) method, which is unconditionally stable and second-order
/// accurate in both space and time.
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `initial_condition`: Initial temperature distribution (flattened 2D grid)
///   - `config`: Solver configuration (grid size, time step, thermal diffusivity, etc.)
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<Vec<f64>, String>` with
/// the final temperature distribution after time evolution.
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[no_mangle]

pub unsafe extern "C" fn rssn_physics_cnm_solve_heat_2d_json(
    input: *const c_char
) -> *mut c_char {

    let input : Heat2DInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(&FfiResult::<
                    Vec<f64>,
                    String,
                >::err(
                    "Invalid JSON".to_string(),
                ))
                .unwrap(),
            )
        },
    };

    let res = physics_cnm::solve_heat_equation_2d_cn_adi(
        &input.initial_condition,
        &input.config,
    );

    to_c_string(
        serde_json::to_string(
            &FfiResult::<
                Vec<f64>,
                String,
            >::ok(res),
        )
        .unwrap(),
    )
}

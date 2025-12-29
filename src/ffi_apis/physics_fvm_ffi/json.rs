//! JSON-based FFI API for physics FVM functions.

use std::os::raw::c_char;

use serde::Deserialize;
use serde::Serialize;

use crate::ffi_apis::common::from_json_string;
use crate::ffi_apis::common::to_c_string;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::physics::physics_fvm::Mesh;
use crate::physics::physics_fvm::SweState;
use crate::physics::physics_fvm::{
    self,
};

#[derive(Deserialize)]

struct AdvectionInput {
    num_cells: usize,
    domain_size: f64,
    velocity: f64,
    dt: f64,
    steps: usize,
    initial_values: Vec<f64>,
}

#[derive(Deserialize)]

struct SweInput {
    h: Vec<f64>,
    hu: Vec<f64>,
    dx: f64,
    dt: f64,
    steps: usize,
    g: f64,
}

/// Solves the 1D advection equation using Finite Volume Method (FVM) via JSON serialization.
///
/// The advection equation ∂u/∂t + v∂u/∂x = 0 models conservative transport
/// of a scalar quantity u with constant velocity v.
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `num_cells`: Number of cells in the mesh
///   - `domain_size`: Total length of the 1D domain
///   - `velocity`: Advection velocity v
///   - `dt`: Time step size
///   - `steps`: Number of time steps to simulate
///   - `initial_values`: Initial field distribution
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<Vec<f64>, String>` with
/// the final field distribution.
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[no_mangle]

pub unsafe extern "C" fn rssn_physics_fvm_advection_json(
    input: *const c_char
) -> *mut c_char {

    let input : AdvectionInput = match from_json_string(input) {
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

    let mut mesh = Mesh::new(
        input.num_cells,
        input.domain_size,
        |_| 0.0,
    );

    for (i, &val) in input
        .initial_values
        .iter()
        .enumerate()
    {

        if i < mesh.cells.len() {

            mesh.cells[i].value = val;
        }
    }

    let result =
        physics_fvm::solve_advection_1d(
            &mut mesh,
            input.velocity,
            input.dt,
            input.steps,
            || (0.0, 0.0),
        );

    to_c_string(
        serde_json::to_string(
            &FfiResult::<
                Vec<f64>,
                String,
            >::ok(result),
        )
        .unwrap(),
    )
}

/// Solves the 1D shallow water equations using Finite Volume Method (FVM) via JSON serialization.
///
/// The shallow water equations model conservation of mass and momentum in free-surface flows:
/// ∂h/∂t + ∂(hu)/∂x = 0 and ∂(hu)/∂t + ∂(hu² + gh²/2)/∂x = 0.
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `h`: Initial water depth distribution
///   - `hu`: Initial momentum (h×velocity) distribution
///   - `dx`: Spatial step size
///   - `dt`: Time step size
///   - `steps`: Number of time steps to simulate
///   - `g`: Gravitational acceleration
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<Vec<SweState>, String>` with
/// the time series of shallow water states (h, hu).
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[no_mangle]

pub unsafe extern "C" fn rssn_physics_fvm_swe_json(
    input: *const c_char
) -> *mut c_char {

    let input : SweInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(&FfiResult::<
                    Vec<SweState>,
                    String,
                >::err(
                    "Invalid JSON".to_string(),
                ))
                .unwrap(),
            )
        },
    };

    let result = physics_fvm::solve_shallow_water_1d(
        input.h,
        input.hu,
        input.dx,
        input.dt,
        input.steps,
        input.g,
    );

    to_c_string(
        serde_json::to_string(
            &FfiResult::<
                Vec<SweState>,
                String,
            >::ok(result),
        )
        .unwrap(),
    )
}

//! JSON-based FFI API for physics FVM functions.

use crate::ffi_apis::common::{from_json_string, to_c_string};
use crate::ffi_apis::ffi_api::FfiResult;
use crate::physics::physics_fvm::{self, Mesh, SweState};
use serde::{Deserialize, Serialize};
use std::os::raw::c_char;

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

#[no_mangle]

pub unsafe extern "C" fn rssn_physics_fvm_advection_json(input: *const c_char) -> *mut c_char {

    let input: AdvectionInput = match from_json_string(input) {
        Some(i) => i,
        None => {
            return to_c_string(
                serde_json::to_string(&FfiResult::<Vec<f64>, String>::err(
                    "Invalid JSON".to_string(),
                ))
                .unwrap(),
            )
        }
    };

    let mut mesh = Mesh::new(input.num_cells, input.domain_size, |_| 0.0);

    for (i, &val) in input.initial_values.iter().enumerate() {

        if i < mesh.cells.len() {

            mesh.cells[i].value = val;
        }
    }

    let result =
        physics_fvm::solve_advection_1d(&mut mesh, input.velocity, input.dt, input.steps, || {

            (0.0, 0.0)
        });

    to_c_string(serde_json::to_string(&FfiResult::<Vec<f64>, String>::ok(result)).unwrap())
}

#[no_mangle]

pub unsafe extern "C" fn rssn_physics_fvm_swe_json(input: *const c_char) -> *mut c_char {

    let input: SweInput = match from_json_string(input) {
        Some(i) => i,
        None => {
            return to_c_string(
                serde_json::to_string(&FfiResult::<Vec<SweState>, String>::err(
                    "Invalid JSON".to_string(),
                ))
                .unwrap(),
            )
        }
    };

    let result = physics_fvm::solve_shallow_water_1d(
        input.h,
        input.hu,
        input.dx,
        input.dt,
        input.steps,
        input.g,
    );

    to_c_string(serde_json::to_string(&FfiResult::<Vec<SweState>, String>::ok(result)).unwrap())
}

//! Bincode-based FFI API for numerical CFD functions.

use crate::ffi_apis::common::{from_bincode_buffer, to_bincode_buffer, BincodeBuffer};
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::physics_cfd;
use serde::Deserialize;

#[derive(Deserialize)]
struct ReynoldsInput {
    velocity: f64,
    length: f64,
    kinematic_viscosity: f64,
}

#[derive(Deserialize)]
struct CflInput {
    velocity: f64,
    dt: f64,
    dx: f64,
}

#[derive(Deserialize)]
struct Advection1DInput {
    u0: Vec<f64>,
    c: f64,
    dx: f64,
    dt: f64,
    num_steps: usize,
}

#[no_mangle]
pub unsafe extern "C" fn rssn_num_cfd_reynolds_number_bincode(
    buffer: BincodeBuffer,
) -> BincodeBuffer {
    let input: ReynoldsInput = match from_bincode_buffer(&buffer) {
        Some(i) => i,
        None => {
            return to_bincode_buffer(&FfiResult::<f64, String> {
                ok: None,
                err: Some("Invalid Bincode".to_string()),
            })
        }
    };
    let re = physics_cfd::reynolds_number(input.velocity, input.length, input.kinematic_viscosity);
    to_bincode_buffer(&FfiResult {
        ok: Some(re),
        err: None::<String>,
    })
}

#[no_mangle]
pub unsafe extern "C" fn rssn_num_cfd_cfl_number_bincode(
    buffer: BincodeBuffer,
) -> BincodeBuffer {
    let input: CflInput = match from_bincode_buffer(&buffer) {
        Some(i) => i,
        None => {
            return to_bincode_buffer(&FfiResult::<f64, String> {
                ok: None,
                err: Some("Invalid Bincode".to_string()),
            })
        }
    };
    let cfl = physics_cfd::cfl_number(input.velocity, input.dt, input.dx);
    to_bincode_buffer(&FfiResult {
        ok: Some(cfl),
        err: None::<String>,
    })
}

#[no_mangle]
pub unsafe extern "C" fn rssn_num_cfd_solve_advection_1d_bincode(
    buffer: BincodeBuffer,
) -> BincodeBuffer {
    let input: Advection1DInput = match from_bincode_buffer(&buffer) {
        Some(i) => i,
        None => {
            return to_bincode_buffer(&FfiResult::<Vec<Vec<f64>>, String> {
                ok: None,
                err: Some("Invalid Bincode".to_string()),
            })
        }
    };
    let results = physics_cfd::solve_advection_1d(&input.u0, input.c, input.dx, input.dt, input.num_steps);
    to_bincode_buffer(&FfiResult {
        ok: Some(results),
        err: None::<String>,
    })
}

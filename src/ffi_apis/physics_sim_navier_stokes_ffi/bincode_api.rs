//! Bincode-based FFI API for physics sim Navier-Stokes functions.

use crate::ffi_apis::common::{from_bincode_buffer, to_bincode_buffer, BincodeBuffer};
use crate::ffi_apis::ffi_api::FfiResult;
use crate::physics::physics_sim::navier_stokes_fluid::{self, NavierStokesParameters};
use ndarray::Array2;
use serde::{Deserialize, Serialize};

#[derive(Serialize)]

struct NavierStokesOutputData {
    pub u: Array2<f64>,
    pub v: Array2<f64>,
    pub p: Array2<f64>,
}

#[no_mangle]

pub unsafe extern "C" fn rssn_physics_sim_navier_stokes_run_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let params: NavierStokesParameters = match from_bincode_buffer(&buffer) {
        Some(p) => p,
        None => {
            return to_bincode_buffer(&FfiResult::<
                NavierStokesOutputData,
                String,
            >::err(
                "Invalid Bincode".to_string(),
            ))
        }
    };

    match navier_stokes_fluid::run_lid_driven_cavity(&params) {
        Ok((u, v, p)) => {

            let out = NavierStokesOutputData { u, v, p };

            to_bincode_buffer(&FfiResult::<
                NavierStokesOutputData,
                String,
            >::ok(
                out
            ))
        }
        Err(e) => to_bincode_buffer(&FfiResult::<
            NavierStokesOutputData,
            String,
        >::err(e)),
    }
}

//! JSON-based FFI API for physics sim Navier-Stokes functions.

use std::os::raw::c_char;

use ndarray::Array2;
use serde::Deserialize;
use serde::Serialize;

use crate::ffi_apis::common::from_json_string;
use crate::ffi_apis::common::to_c_string;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::physics::physics_sim::navier_stokes_fluid::NavierStokesParameters;
use crate::physics::physics_sim::navier_stokes_fluid::{
    self,
};

#[derive(Serialize)]

struct NavierStokesOutputData {
    pub u : Array2<f64>,
    pub v : Array2<f64>,
    pub p : Array2<f64>,
}

#[no_mangle]

pub unsafe extern "C" fn rssn_physics_sim_navier_stokes_run_json(
    input : *const c_char
) -> *mut c_char {

    let params : NavierStokesParameters = match from_json_string(input) {
        | Some(p) => p,
        | None => {
            return to_c_string(
                serde_json::to_string(&FfiResult::<
                    NavierStokesOutputData,
                    String,
                >::err(
                    "Invalid JSON".to_string(),
                ))
                .unwrap(),
            )
        },
    };

    match navier_stokes_fluid::run_lid_driven_cavity(&params) {
        | Ok((u, v, p)) => {

            let out = NavierStokesOutputData {
                u,
                v,
                p,
            };

            to_c_string(
                serde_json::to_string(&FfiResult::<
                    NavierStokesOutputData,
                    String,
                >::ok(
                    out
                ))
                .unwrap(),
            )
        },
        | Err(e) => {
            to_c_string(
                serde_json::to_string(&FfiResult::<
                    NavierStokesOutputData,
                    String,
                >::err(
                    e
                ))
                .unwrap(),
            )
        },
    }
}

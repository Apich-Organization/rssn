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
    pub u: Array2<f64>,
    pub v: Array2<f64>,
    pub p: Array2<f64>,
}

/// Solves the incompressible Navier-Stokes equations for fluid flow in a lid-driven cavity via JSON serialization.
///
/// The Navier-Stokes equations ∂u/∂t + (u·∇)u = -∇p/ρ + ν∇²u with incompressibility
/// constraint ∇·u = 0 govern viscous fluid dynamics. This solver uses a projection method
/// to enforce divergence-free velocity fields.
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `nx`, `ny`: Grid dimensions
///   - `re`: Reynolds number Re = UL/ν (ratio of inertial to viscous forces)
///   - `dt`: Time step size
///   - `n_iter`: Number of time iterations
///   - `lid_velocity`: Velocity of the moving lid boundary
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<NavierStokesOutputData, String>` with:
/// - `u`: Horizontal velocity field u(x,y)
/// - `v`: Vertical velocity field v(x,y)
/// - `p`: Pressure field p(x,y)
///
/// # Safety
///
/// This function is unsafe because it dereferences a raw C string pointer.
#[no_mangle]

pub unsafe extern "C" fn rssn_physics_sim_navier_stokes_run_json(
    input: *const c_char
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

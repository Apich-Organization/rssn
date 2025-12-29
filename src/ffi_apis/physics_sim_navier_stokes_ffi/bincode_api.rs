//! Bincode-based FFI API for physics sim Navier-Stokes functions.

use ndarray::Array2;
use serde::Deserialize;
use serde::Serialize;

use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::ffi_apis::common::BincodeBuffer;
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

/// Solves the incompressible Navier-Stokes equations for fluid flow in a lid-driven cavity via bincode serialization.
///
/// The Navier-Stokes equations ∂u/∂t + (u·∇)u = -∇p/ρ + ν∇²u with incompressibility
/// constraint ∇·u = 0 govern viscous fluid dynamics. This solver uses a projection method
/// to enforce divergence-free velocity fields.
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing `NavierStokesParameters` with:
///   - `nx`, `ny`: Grid dimensions
///   - `re`: Reynolds number Re = UL/ν (ratio of inertial to viscous forces)
///   - `dt`: Time step size
///   - `n_iter`: Number of time iterations
///   - `lid_velocity`: Velocity of the moving lid boundary
///
/// # Returns
///
/// A bincode-encoded buffer containing `FfiResult<NavierStokesOutputData, String>` with either:
/// - `ok`: Object containing:
///   - `u`: Horizontal velocity field u(x,y)
///   - `v`: Vertical velocity field v(x,y)
///   - `p`: Pressure field p(x,y)
/// - `err`: Error message if computation failed
///
/// # Safety
///
/// This function is unsafe because it receives a raw bincode buffer that must be
/// valid and properly encoded.
#[no_mangle]

pub unsafe extern "C" fn rssn_physics_sim_navier_stokes_run_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let params : NavierStokesParameters = match from_bincode_buffer(&buffer) {
        | Some(p) => p,
        | None => {
            return to_bincode_buffer(&FfiResult::<
                NavierStokesOutputData,
                String,
            >::err(
                "Invalid Bincode".to_string(),
            ))
        },
    };

    match navier_stokes_fluid::run_lid_driven_cavity(&params) {
        | Ok((u, v, p)) => {

            let out = NavierStokesOutputData {
                u,
                v,
                p,
            };

            to_bincode_buffer(&FfiResult::<
                NavierStokesOutputData,
                String,
            >::ok(
                out
            ))
        },
        | Err(e) => {
            to_bincode_buffer(&FfiResult::<
                NavierStokesOutputData,
                String,
            >::err(
                e
            ))
        },
    }
}

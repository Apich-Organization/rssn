//! Bincode-based FFI API for numerical CFD functions.

use serde::Deserialize;

use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::ffi_apis::common::BincodeBuffer;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::physics_cfd;

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

/// Computes the Reynolds number for fluid flow using bincode serialization.
///
/// The Reynolds number is a dimensionless quantity characterizing the flow regime:
/// Re = (velocity × length) / `kinematic_viscosity`
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing `ReynoldsInput` with:
///   - `velocity`: Flow velocity (m/s)
///   - `length`: Characteristic length scale (m)
///   - `kinematic_viscosity`: Kinematic viscosity (ν = μ/ρ, m²/s)
///
/// # Returns
///
/// A bincode-encoded buffer containing `FfiResult<f64, String>` with either:
/// - `ok`: The computed Reynolds number (dimensionless)
/// - `err`: Error message if deserialization failed
///
/// # Safety
///
/// This function is unsafe because it receives raw pointers through FFI.
/// The caller must ensure the input buffer contains valid bincode data.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_cfd_reynolds_number_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input : ReynoldsInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(
                &FfiResult::<f64, String> {
                    ok : None,
                    err : Some("Invalid Bincode".to_string()),
                },
            )
        },
    };

    let re =
        physics_cfd::reynolds_number(
            input.velocity,
            input.length,
            input.kinematic_viscosity,
        );

    to_bincode_buffer(&FfiResult {
        ok: Some(re),
        err: None::<String>,
    })
}

/// Computes the Courant-Friedrichs-Lewy (CFL) number for numerical stability using bincode serialization.
///
/// The CFL number is a stability criterion for explicit time-stepping schemes:
/// CFL = (velocity × dt) / dx
///
/// For numerical stability in explicit schemes, CFL ≤ 1 is typically required.
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing `CflInput` with:
///   - `velocity`: Flow velocity or wave speed (m/s)
///   - `dt`: Time step size (s)
///   - `dx`: Spatial grid spacing (m)
///
/// # Returns
///
/// A bincode-encoded buffer containing `FfiResult<f64, String>` with either:
/// - `ok`: The computed CFL number (dimensionless)
/// - `err`: Error message if deserialization failed
///
/// # Safety
///
/// This function is unsafe because it receives raw pointers through FFI.
/// The caller must ensure the input buffer contains valid bincode data.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_cfd_cfl_number_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input : CflInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(
                &FfiResult::<f64, String> {
                    ok : None,
                    err : Some("Invalid Bincode".to_string()),
                },
            )
        },
    };

    let cfl = physics_cfd::cfl_number(
        input.velocity,
        input.dt,
        input.dx,
    );

    to_bincode_buffer(&FfiResult {
        ok: Some(cfl),
        err: None::<String>,
    })
}

/// Solves the 1D advection equation using a finite difference scheme and bincode serialization.
///
/// The advection equation describes the transport of a scalar quantity by a velocity field:
/// ∂u/∂t + c ∂u/∂x = 0
///
/// This function uses an explicit finite difference method to time-step the solution.
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing `Advection1DInput` with:
///   - `u0`: Initial condition as a spatial array of values
///   - `c`: Advection velocity (constant)
///   - `dx`: Spatial grid spacing
///   - `dt`: Time step size
///   - `num_steps`: Number of time steps to compute
///
/// # Returns
///
/// A bincode-encoded buffer containing `FfiResult<Vec<Vec<f64>>, String>` with either:
/// - `ok`: Solution history where each inner vector is the solution at one time step
/// - `err`: Error message if deserialization failed
///
/// # Safety
///
/// This function is unsafe because it receives raw pointers through FFI.
/// The caller must ensure the input buffer contains valid bincode data.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_cfd_solve_advection_1d_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input : Advection1DInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(
                &FfiResult::<Vec<Vec<f64>>, String> {
                    ok : None,
                    err : Some("Invalid Bincode".to_string()),
                },
            )
        },
    };

    let results =
        physics_cfd::solve_advection_1d(
            &input.u0,
            input.c,
            input.dx,
            input.dt,
            input.num_steps,
        );

    to_bincode_buffer(&FfiResult {
        ok: Some(results),
        err: None::<String>,
    })
}

//! JSON-based FFI API for physics SM functions.

use std::os::raw::c_char;

use serde::Deserialize;

use crate::ffi_apis::common::from_json_string;
use crate::ffi_apis::common::to_c_string;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::physics::physics_sm::AdvectionDiffusionConfig;
use crate::physics::physics_sm::{
    self,
};

#[derive(Deserialize)]

struct AdvectionDiffusion1DInput {
    initial_condition: Vec<f64>,
    dx: f64,
    c: f64,
    d: f64,
    dt: f64,
    steps: usize,
}

#[derive(Deserialize)]

struct AdvectionDiffusion2DInput {
    initial_condition: Vec<f64>,
    config: AdvectionDiffusionConfig,
}

/// Solves the 1D advection-diffusion equation using spectral methods via JSON serialization.
///
/// The 1D advection-diffusion equation ∂u/∂t + c∂u/∂x = D∂²u/∂x² models scalar transport
/// with velocity c and diffusivity D. Spectral methods provide exponential convergence
/// for smooth solutions.
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `initial_condition`: Initial scalar field u(x,0)
///   - `dx`: Spatial step size
///   - `c`: Advection velocity
///   - `d`: Diffusion coefficient D
///   - `dt`: Time step size
///   - `steps`: Number of time steps
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<Vec<f64>, String>` with
/// the final scalar field u(x,t).
///
/// # Safety
///
/// This function is unsafe because it dereferences a raw C string pointer.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_physics_sm_solve_advection_1d_json(
    input: *const c_char
) -> *mut c_char {

    let input : AdvectionDiffusion1DInput = match from_json_string(input) {
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

    let res = physics_sm::solve_advection_diffusion_1d(
        &input.initial_condition,
        input.dx,
        input.c,
        input.d,
        input.dt,
        input.steps,
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

/// Solves the 2D advection-diffusion equation using spectral methods via JSON serialization.
///
/// The 2D advection-diffusion equation ∂u/∂t + c·∇u = D∇²u models transport phenomena
/// combining convective transport (advection) and diffusive spreading. Spectral methods
/// use Fourier basis functions for high-order accuracy.
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `initial_condition`: Initial scalar field u(x,y,0) as flattened vector
///   - `config`: Configuration including:
///     - `nx`, `ny`: Grid dimensions
///     - `dx`, `dy`: Spatial steps
///     - `cx`, `cy`: Advection velocities in x and y directions
///     - `d`: Diffusion coefficient D
///     - `dt`: Time step size
///     - `steps`: Number of time steps
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<Vec<f64>, String>` with
/// the final scalar field u(x,y,t) as a flattened vector.
///
/// # Safety
///
/// This function is unsafe because it dereferences a raw C string pointer.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_physics_sm_solve_advection_2d_json(
    input: *const c_char
) -> *mut c_char {

    let input : AdvectionDiffusion2DInput = match from_json_string(input) {
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

    let res = physics_sm::solve_advection_diffusion_2d(
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

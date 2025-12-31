//! Bincode-based FFI API for physics SM functions.

use serde::Deserialize;

use crate::ffi_apis::common::BincodeBuffer;
use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::physics::physics_sm::{
    self,
}; // Alias for clarity if needed

#[derive(Deserialize)]

struct AdvectionDiffusion2DInput {
    initial_condition : Vec<f64>,
    config : physics_sm::AdvectionDiffusionConfig,
}

/// Solves the 2D advection-diffusion equation using spectral methods via bincode serialization.
///
/// The advection-diffusion equation ∂u/∂t + c·∇u = D∇²u models transport phenomena
/// combining convective transport (advection) and diffusive spreading. Spectral methods
/// use Fourier basis functions for high-order accuracy.
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing:
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
/// A bincode-encoded buffer containing `FfiResult<Vec<f64>, String>` with
/// the final scalar field u(x,y,t) as a flattened vector.
///
/// # Safety
///
/// This function is unsafe because it receives a raw bincode buffer that must be
/// valid and properly encoded.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_physics_sm_solve_advection_2d_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input : AdvectionDiffusion2DInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(&FfiResult::<
                Vec<f64>,
                String,
            >::err(
                "Invalid Bincode".to_string(),
            ))
        },
    };

    let res = physics_sm::solve_advection_diffusion_2d(
        &input.initial_condition,
        &input.config,
    );

    to_bincode_buffer(&FfiResult::<
        Vec<f64>,
        String,
    >::ok(res))
}

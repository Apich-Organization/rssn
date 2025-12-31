//! Bincode-based FFI API for physics FVM functions.

use serde::Deserialize;

use crate::ffi_apis::common::BincodeBuffer;
use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::physics::physics_fvm::SweState;
use crate::physics::physics_fvm::{
    self,
};

#[derive(Deserialize)]

struct SweInput {
    h: Vec<f64>,
    hu: Vec<f64>,
    dx: f64,
    dt: f64,
    steps: usize,
    g: f64,
}

/// Solves the 1D shallow water equations using Finite Volume Method (FVM) via bincode serialization.
///
/// The shallow water equations model conservation of mass and momentum in free-surface flows:
/// ∂h/∂t + ∂(hu)/∂x = 0 and ∂(hu)/∂t + ∂(hu² + gh²/2)/∂x = 0.
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing `SweInput` with:
///   - `h`: Initial water depth distribution
///   - `hu`: Initial momentum (h×velocity) distribution
///   - `dx`: Spatial step size
///   - `dt`: Time step size
///   - `steps`: Number of time steps to simulate
///   - `g`: Gravitational acceleration
///
/// # Returns
///
/// A bincode-encoded buffer containing `FfiResult<Vec<SweState>, String>` with either:
/// - `ok`: Time series of shallow water states (h, hu)
/// - `err`: Error message if computation failed
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

pub unsafe extern "C" fn rssn_physics_fvm_swe_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input : SweInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(&FfiResult::<
                Vec<SweState>,
                String,
            >::err(
                "Invalid Bincode".to_string(),
            ))
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

    to_bincode_buffer(&FfiResult::<
        Vec<SweState>,
        String,
    >::ok(result))
}

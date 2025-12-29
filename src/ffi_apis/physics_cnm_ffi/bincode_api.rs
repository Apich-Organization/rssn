//! Bincode-based FFI API for physics CNM functions.

use serde::Deserialize;

use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::ffi_apis::common::BincodeBuffer;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::physics::physics_cnm::HeatEquationSolverConfig;
use crate::physics::physics_cnm::{
    self,
};

#[derive(Deserialize)]

struct Heat2DInput {
    initial_condition: Vec<f64>,
    config: HeatEquationSolverConfig,
}

/// Solves the 2D heat equation using Crank-Nicolson ADI method via bincode serialization.
///
/// The heat equation ∂u/∂t = α∇²u is solved using the Crank-Nicolson Alternating
/// Direction Implicit (ADI) method, which is unconditionally stable and second-order
/// accurate in both space and time.
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing `Heat2DInput` with:
///   - `initial_condition`: Initial temperature distribution (flattened 2D grid)
///   - `config`: Solver configuration (grid size, time step, thermal diffusivity, etc.)
///
/// # Returns
///
/// A bincode-encoded buffer containing `FfiResult<Vec<f64>, String>` with either:
/// - `ok`: Final temperature distribution after time evolution
/// - `err`: Error message if computation failed
///
/// # Safety
///
/// This function is unsafe because it receives a raw bincode buffer that must be
/// valid and properly encoded.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_physics_cnm_solve_heat_2d_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input : Heat2DInput = match from_bincode_buffer(&buffer) {
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

    let res = physics_cnm::solve_heat_equation_2d_cn_adi(
        &input.initial_condition,
        &input.config,
    );

    to_bincode_buffer(&FfiResult::<
        Vec<f64>,
        String,
    >::ok(res))
}

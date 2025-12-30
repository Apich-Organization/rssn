//! Bincode-based FFI API for physics sim FDTD electrodynamics functions.

use ndarray::Array2;

use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::ffi_apis::common::BincodeBuffer;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::physics::physics_sim::fdtd_electrodynamics::FdtdParameters;
use crate::physics::physics_sim::fdtd_electrodynamics::{
    self,
};

/// Runs a Finite-Difference Time-Domain (FDTD) electromagnetic simulation via bincode serialization.
///
/// FDTD solves Maxwell's equations ∇×E = -∂B/∂t and ∇×H = ∂D/∂t + J using a staggered
/// Yee lattice grid, advancing the electric field Ez and magnetic field components in time.
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing `FdtdParameters` with:
///   - `width`, `height`: Grid dimensions
///   - `dx`, `dy`: Spatial discretization steps
///   - `dt`: Time step size (must satisfy Courant-Friedrichs-Lewy stability condition)
///   - `steps`: Number of time steps to simulate
///   - `source_x`, `source_y`: Position of electromagnetic source
///   - `source_frequency`: Angular frequency ω of the source
///
/// # Returns
///
/// A bincode-encoded buffer containing `FfiResult<Array2<f64>, String>` with either:
/// - `ok`: Final Ez field as a 2D array
/// - `err`: Error message if computation failed or no snapshots were produced
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

pub unsafe extern "C" fn rssn_physics_sim_fdtd_run_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let params : FdtdParameters = match from_bincode_buffer(&buffer) {
        | Some(p) => p,
        | None => {
            return to_bincode_buffer(&FfiResult::<
                Vec<f64>,
                String,
            >::err(
                "Invalid Bincode".to_string(),
            ))
        },
    };

    let snapshots = fdtd_electrodynamics::run_fdtd_simulation(&params);

    if let Some(final_ez) =
        snapshots.last()
    {

        to_bincode_buffer(&FfiResult::<
            Array2<f64>,
            String,
        >::ok(
            final_ez.clone(),
        ))
    } else {

        to_bincode_buffer(&FfiResult::<
            Array2<f64>,
            String,
        >::err(
            "No snapshots".to_string(),
        ))
    }
}

//! Bincode-based FFI API for physics sim GPE superfluidity functions.

use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::ffi_apis::common::BincodeBuffer;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::physics::physics_sim::gpe_superfluidity::GpeParameters;
use crate::physics::physics_sim::gpe_superfluidity::{
    self,
};

/// Solves the Gross-Pitaevskii equation (GPE) for Bose-Einstein condensate ground state via bincode serialization.
///
/// The GPE iℏ∂ψ/∂t = [-ℏ²∇²/(2m) + V(r) + g|ψ|²]ψ describes the macroscopic wavefunction
/// of a superfluid quantum gas. This solver finds the ground state using imaginary time
/// evolution or variational methods.
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing `GpeParameters` with:
///   - `n_points`: Number of spatial grid points
///   - `dx`: Spatial discretization step
///   - `g`: Nonlinear interaction strength (proportional to scattering length)
///   - `v_trap`: External trapping potential coefficients
///   - `tolerance`: Convergence tolerance for ground state search
///   - `max_iterations`: Maximum iterations for solver
///
/// # Returns
///
/// A bincode-encoded buffer containing `FfiResult<Vec<f64>, String>` with either:
/// - `ok`: Ground state wavefunction ψ(x) as probability density |ψ|²
/// - `err`: Error message if computation failed or did not converge
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

pub unsafe extern "C" fn rssn_physics_sim_gpe_run_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let params : GpeParameters = match from_bincode_buffer(&buffer) {
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

    match gpe_superfluidity::run_gpe_ground_state_finder(&params) {
        | Ok(res) => {
            to_bincode_buffer(&FfiResult::<
                Vec<f64>,
                String,
            >::ok(
                res.into_raw_vec_and_offset().0,
            ))
        },
        | Err(e) => {
            to_bincode_buffer(&FfiResult::<
                Vec<f64>,
                String,
            >::err(
                e
            ))
        },
    }
}

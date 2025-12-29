//! JSON-based FFI API for physics sim GPE superfluidity functions.

use std::os::raw::c_char;

use crate::ffi_apis::common::from_json_string;
use crate::ffi_apis::common::to_c_string;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::physics::physics_sim::gpe_superfluidity::GpeParameters;
use crate::physics::physics_sim::gpe_superfluidity::{
    self,
};

/// Solves the Gross-Pitaevskii equation (GPE) for Bose-Einstein condensate ground state via JSON serialization.
///
/// The GPE iℏ∂ψ/∂t = [-ℏ²∇²/(2m) + V(r) + g|ψ|²]ψ describes the macroscopic wavefunction
/// of a superfluid quantum gas. This solver finds the ground state using imaginary time
/// evolution or variational methods.
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `n_points`: Number of spatial grid points
///   - `dx`: Spatial discretization step
///   - `g`: Nonlinear interaction strength (proportional to scattering length)
///   - `v_trap`: External trapping potential coefficients
///   - `tolerance`: Convergence tolerance for ground state search
///   - `max_iterations`: Maximum iterations for solver
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<Vec<f64>, String>` with
/// the ground state wavefunction ψ(x) as probability density |ψ|².
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

/// # Panics
///
/// This function may panic if the FFI input is malformed, null where not expected,
/// or if internal state synchronization fails (e.g., poisoned locks).

pub unsafe extern "C" fn rssn_physics_sim_gpe_run_json(
    input: *const c_char
) -> *mut c_char {

    let params : GpeParameters = match from_json_string(input) {
        | Some(p) => p,
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

    match gpe_superfluidity::run_gpe_ground_state_finder(&params) {
        | Ok(res) => {
            to_c_string(
                serde_json::to_string(&FfiResult::<
                    Vec<f64>,
                    String,
                >::ok(
                    res.into_raw_vec_and_offset().0,
                ))
                .unwrap(),
            )
        },
        | Err(e) => {
            to_c_string(
                serde_json::to_string(&FfiResult::<
                    Vec<f64>,
                    String,
                >::err(
                    e
                ))
                .unwrap(),
            )
        },
    }
}

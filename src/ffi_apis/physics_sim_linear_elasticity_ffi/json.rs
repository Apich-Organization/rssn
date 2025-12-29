//! JSON-based FFI API for physics sim linear elasticity functions.

use std::os::raw::c_char;

use crate::ffi_apis::common::from_json_string;
use crate::ffi_apis::common::to_c_string;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::physics::physics_sim::linear_elasticity::ElasticityParameters;
use crate::physics::physics_sim::linear_elasticity::{
    self,
};

/// Solves the linear elasticity equations for stress and displacement in a deformable solid via JSON serialization.
///
/// The elasticity equations σᵢⱼ = Cᵢⱼₖₗεₖₗ with strain εᵢⱼ = ½(∂uᵢ/∂xⱼ + ∂uⱼ/∂xᵢ) describe
/// small deformations in elastic materials, where σ is stress tensor, ε is strain tensor,
/// u is displacement field, and C is the stiffness tensor.
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `nx`, `ny`: Grid dimensions for finite element discretization
///   - `young_modulus`: Young's modulus E (stiffness)
///   - `poisson_ratio`: Poisson's ratio ν (lateral contraction)
///   - `applied_force`: External force distribution
///   - `boundary_conditions`: Fixed displacement constraints
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<Vec<f64>, String>` with
/// the displacement field u(x,y) as a flattened vector.
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

pub unsafe extern "C" fn rssn_physics_sim_linear_elasticity_run_json(
    input: *const c_char
) -> *mut c_char {

    let params : ElasticityParameters = match from_json_string(input) {
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

    match linear_elasticity::run_elasticity_simulation(&params) {
        | Ok(res) => {
            to_c_string(
                serde_json::to_string(&FfiResult::<
                    Vec<f64>,
                    String,
                >::ok(
                    res
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

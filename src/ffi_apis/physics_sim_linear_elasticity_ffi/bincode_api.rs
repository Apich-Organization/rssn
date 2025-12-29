//! Bincode-based FFI API for physics sim linear elasticity functions.

use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::ffi_apis::common::BincodeBuffer;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::physics::physics_sim::linear_elasticity::ElasticityParameters;
use crate::physics::physics_sim::linear_elasticity::{
    self,
};

/// Solves the linear elasticity equations for stress and displacement in a deformable solid via bincode serialization.
///
/// The elasticity equations σᵢⱼ = Cᵢⱼₖₗεₖₗ with strain εᵢⱼ = ½(∂uᵢ/∂xⱼ + ∂uⱼ/∂xᵢ) describe
/// small deformations in elastic materials, where σ is stress tensor, ε is strain tensor,
/// u is displacement field, and C is the stiffness tensor.
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing `ElasticityParameters` with:
///   - `nx`, `ny`: Grid dimensions for finite element discretization
///   - `young_modulus`: Young's modulus E (stiffness)
///   - `poisson_ratio`: Poisson's ratio ν (lateral contraction)
///   - `applied_force`: External force distribution
///   - `boundary_conditions`: Fixed displacement constraints
///
/// # Returns
///
/// A bincode-encoded buffer containing `FfiResult<Vec<f64>, String>` with either:
/// - `ok`: Displacement field u(x,y) as flattened vector
/// - `err`: Error message if computation failed
///
/// # Safety
///
/// This function is unsafe because it receives a raw bincode buffer that must be
/// valid and properly encoded.
#[no_mangle]

pub unsafe extern "C" fn rssn_physics_sim_linear_elasticity_run_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let params : ElasticityParameters = match from_bincode_buffer(&buffer) {
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

    match linear_elasticity::run_elasticity_simulation(&params) {
        | Ok(res) => {
            to_bincode_buffer(&FfiResult::<
                Vec<f64>,
                String,
            >::ok(
                res
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

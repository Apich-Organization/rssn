//! Bincode-based FFI API for physics sim geodesic relativity functions.

use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::ffi_apis::common::BincodeBuffer;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::physics::physics_sim::geodesic_relativity::GeodesicParameters;
use crate::physics::physics_sim::geodesic_relativity::{
    self,
};

/// Computes a geodesic trajectory in curved spacetime using general relativity via bincode serialization.
///
/// Integrates the geodesic equation d²xᵘ/dτ² + Γᵘᵥᵨ(dxᵥ/dτ)(dxᵨ/dτ) = 0 where Γᵘᵥᵨ are
/// Christoffel symbols of the metric tensor, modeling particle motion in curved spacetime
/// (e.g., near a black hole using Schwarzschild metric).
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing `GeodesicParameters` with:
///   - `mass`: Central mass M (e.g., black hole mass)
///   - `r0`, `phi0`: Initial radial and angular coordinates
///   - `v_r`, `v_phi`: Initial radial and angular velocities
///   - `dt`: Time step for integration
///   - `steps`: Number of integration steps
///
/// # Returns
///
/// A bincode-encoded buffer containing `FfiResult<Vec<(f64, f64)>, String>` with
/// the geodesic path as (r, φ) coordinate pairs.
///
/// # Safety
///
/// This function is unsafe because it receives a raw bincode buffer that must be
/// valid and properly encoded.
#[no_mangle]

pub unsafe extern "C" fn rssn_physics_sim_geodesic_run_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let params : GeodesicParameters = match from_bincode_buffer(&buffer) {
        | Some(p) => p,
        | None => {
            return to_bincode_buffer(&FfiResult::<
                Vec<(f64, f64)>,
                String,
            >::err(
                "Invalid Bincode".to_string(),
            ))
        },
    };

    let path = geodesic_relativity::run_geodesic_simulation(&params);

    to_bincode_buffer(&FfiResult::<
        Vec<(f64, f64)>,
        String,
    >::ok(path))
}

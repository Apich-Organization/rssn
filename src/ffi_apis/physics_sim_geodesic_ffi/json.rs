//! JSON-based FFI API for physics sim geodesic relativity functions.

use std::os::raw::c_char;

use crate::ffi_apis::common::from_json_string;
use crate::ffi_apis::common::to_c_string;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::physics::physics_sim::geodesic_relativity::GeodesicParameters;
use crate::physics::physics_sim::geodesic_relativity::{
    self,
};

/// Computes a geodesic trajectory in curved spacetime using general relativity via JSON serialization.
///
/// Integrates the geodesic equation d²xᵘ/dτ² + Γᵘᵥᵨ(dxᵥ/dτ)(dxᵨ/dτ) = 0 where Γᵘᵥᵨ are
/// Christoffel symbols of the metric tensor, modeling particle motion in curved spacetime
/// (e.g., near a black hole using Schwarzschild metric).
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `mass`: Central mass M (e.g., black hole mass)
///   - `r0`, `phi0`: Initial radial and angular coordinates
///   - `v_r`, `v_phi`: Initial radial and angular velocities
///   - `dt`: Time step for integration
///   - `steps`: Number of integration steps
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<Vec<(f64, f64)>, String>` with
/// the geodesic path as (r, φ) coordinate pairs.
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

pub unsafe extern "C" fn rssn_physics_sim_geodesic_run_json(
    input: *const c_char
) -> *mut c_char {

    let params : GeodesicParameters = match from_json_string(input) {
        | Some(p) => p,
        | None => {
            return to_c_string(
                serde_json::to_string(&FfiResult::<
                    Vec<(f64, f64)>,
                    String,
                >::err(
                    "Invalid JSON".to_string(),
                ))
                .unwrap(),
            )
        },
    };

    let path = geodesic_relativity::run_geodesic_simulation(&params);

    to_c_string(
        serde_json::to_string(
            &FfiResult::<
                Vec<(f64, f64)>,
                String,
            >::ok(path),
        )
        .unwrap(),
    )
}

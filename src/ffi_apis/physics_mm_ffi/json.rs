//! JSON-based FFI API for physics MM functions.

use std::os::raw::c_char;

use serde::Deserialize;
use serde::Serialize;

use crate::ffi_apis::common::from_json_string;
use crate::ffi_apis::common::to_c_string;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::physics::physics_mm::Particle;
use crate::physics::physics_mm::SPHSystem;
use crate::physics::physics_mm::{
    self,
};

#[derive(Deserialize)]

struct SphInput {
    system: SPHSystem,
    dt: f64,
}

/// Updates a Smoothed Particle Hydrodynamics (SPH) system by one time step via JSON serialization.
///
/// SPH is a meshfree Lagrangian method for simulating fluid dynamics by representing
/// the continuum as a set of particles with smoothed properties.
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `system`: SPH system state (particles with positions, velocities, densities, etc.)
///   - `dt`: Time step size
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<SPHSystem, String>` with
/// the updated SPH system state after time step.
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[no_mangle]

pub unsafe extern "C" fn rssn_physics_mm_sph_update_json(
    input: *const c_char
) -> *mut c_char {

    let mut input : SphInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(&FfiResult::<
                    SPHSystem,
                    String,
                >::err(
                    "Invalid JSON".to_string(),
                ))
                .unwrap(),
            )
        },
    };

    input
        .system
        .update(input.dt);

    to_c_string(
        serde_json::to_string(
            &FfiResult::<
                SPHSystem,
                String,
            >::ok(
                input.system
            ),
        )
        .unwrap(),
    )
}

/// Simulates a 2D dam break scenario using SPH method via JSON serialization.
///
/// Models the collapse of a water column and its subsequent flow, a classical
/// validation case for SPH fluid simulation.
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<Vec<(f64, f64)>, String>` with
/// particle positions (x, y) after simulation.
///
/// # Safety
///
/// This function is unsafe because it returns a raw C string pointer that the
/// caller must free.
#[no_mangle]

pub unsafe extern "C" fn rssn_physics_mm_simulate_dam_break_json(
) -> *mut c_char {

    let res = physics_mm::simulate_dam_break_2d_scenario();

    to_c_string(
        serde_json::to_string(
            &FfiResult::<
                Vec<(f64, f64)>,
                String,
            >::ok(res),
        )
        .unwrap(),
    )
}

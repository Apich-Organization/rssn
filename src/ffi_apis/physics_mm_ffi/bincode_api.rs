//! Bincode-based FFI API for physics MM functions.

use serde::Deserialize;

use crate::ffi_apis::common::BincodeBuffer;
use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::physics::physics_mm::SPHSystem;

#[derive(Deserialize)]

struct SphInput {
    system: SPHSystem,
    dt: f64,
}

/// Updates a Smoothed Particle Hydrodynamics (SPH) system by one time step via bincode serialization.
///
/// SPH is a meshfree Lagrangian method for simulating fluid dynamics by representing
/// the continuum as a set of particles with smoothed properties.
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing `SphInput` with:
///   - `system`: SPH system state (particles with positions, velocities, densities, etc.)
///   - `dt`: Time step size
///
/// # Returns
///
/// A bincode-encoded buffer containing `FfiResult<SPHSystem, String>` with either:
/// - `ok`: Updated SPH system state after time step
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

pub unsafe extern "C" fn rssn_physics_mm_sph_update_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let mut input : SphInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(&FfiResult::<
                SPHSystem,
                String,
            >::err(
                "Invalid Bincode".to_string(),
            ))
        },
    };

    input
        .system
        .update(input.dt);

    to_bincode_buffer(&FfiResult::<
        SPHSystem,
        String,
    >::ok(
        input.system
    ))
}

//! Bincode-based FFI API for numerical MD functions.

use serde::Deserialize;
use serde::Serialize;

use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::ffi_apis::common::BincodeBuffer;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::physics_md;

#[derive(Deserialize)]

struct LennardJonesInput {
    p1_position: Vec<f64>,
    p2_position: Vec<f64>,
    epsilon: f64,
    sigma: f64,
}

#[derive(Serialize)]

struct InteractionOutput {
    potential: f64,
    force: Vec<f64>,
}

#[derive(Deserialize)]

struct PbcInput {
    position: Vec<f64>,
    box_size: Vec<f64>,
}

/// Computes the Lennard-Jones interaction potential and force between two particles using bincode serialization.
///
/// The Lennard-Jones potential models van der Waals interactions:
/// V(r) = 4ε[(σ/r)¹² - (σ/r)⁶], where ε is the depth of the potential well and σ is the finite distance at which the potential is zero.
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing `LennardJonesInput` with:
///   - `p1_position`: Position vector of first particle (x, y, z)
///   - `p2_position`: Position vector of second particle (x, y, z)
///   - `epsilon`: Well depth ε (energy units)
///   - `sigma`: Finite distance σ at which V=0 (length units)
///
/// # Returns
///
/// A bincode-encoded buffer containing `FfiResult<InteractionOutput, String>` with either:
/// - `ok`: Object containing:
///   - `potential`: Interaction potential energy V(r)
///   - `force`: Force vector acting on particle 1
/// - `err`: Error message if computation failed
///
/// # Safety
///
/// This function is unsafe because it receives raw pointers through FFI.
/// The caller must ensure the input buffer contains valid bincode data.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_md_lennard_jones_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input : LennardJonesInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(
                &FfiResult::<InteractionOutput, String> {
                    ok : None,
                    err : Some("Invalid Bincode".to_string()),
                },
            )
        },
    };

    let p1 = physics_md::Particle::new(
        0,
        1.0,
        input.p1_position,
        vec![0.0, 0.0, 0.0],
    );

    let p2 = physics_md::Particle::new(
        1,
        1.0,
        input.p2_position,
        vec![0.0, 0.0, 0.0],
    );

    match physics_md::lennard_jones_interaction(
        &p1,
        &p2,
        input.epsilon,
        input.sigma,
    ) {
        | Ok((potential, force)) => {

            let output = InteractionOutput {
                potential,
                force,
            };

            to_bincode_buffer(&FfiResult {
                ok : Some(output),
                err : None::<String>,
            })
        },
        | Err(e) => {
            to_bincode_buffer(
                &FfiResult::<InteractionOutput, String> {
                    ok : None,
                    err : Some(e),
                },
            )
        },
    }
}

/// Applies periodic boundary conditions (PBC) to wrap particle coordinates into the simulation box using bincode serialization.
///
/// Periodic boundary conditions create an infinite tiling of the simulation box,
/// ensuring particles that exit one side re-enter from the opposite side.
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing `PbcInput` with:
///   - `position`: Position vector to wrap (x, y, z)
///   - `box_size`: Simulation box dimensions (Lx, Ly, Lz)
///
/// # Returns
///
/// A bincode-encoded buffer containing `FfiResult<Vec<f64>, String>` with either:
/// - `ok`: Wrapped position vector within [0, `box_size`)
/// - `err`: Error message if input invalid
///
/// # Safety
///
/// This function is unsafe because it receives raw pointers through FFI.
/// The caller must ensure the input buffer contains valid bincode data.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_md_apply_pbc_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input : PbcInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(
                &FfiResult::<Vec<f64>, String> {
                    ok : None,
                    err : Some("Invalid Bincode".to_string()),
                },
            )
        },
    };

    let wrapped = physics_md::apply_pbc(
        &input.position,
        &input.box_size,
    );

    to_bincode_buffer(&FfiResult {
        ok: Some(wrapped),
        err: None::<String>,
    })
}

//! JSON-based FFI API for numerical MD functions.

use std::os::raw::c_char;

use serde::Deserialize;
use serde::Serialize;

use crate::ffi_apis::common::from_json_string;
use crate::ffi_apis::common::to_c_string;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::physics_md;

// ============================================================================
// Input/Output structs
// ============================================================================

#[derive(Deserialize)]

struct ParticleInput {
    id: usize,
    mass: f64,
    position: Vec<f64>,
    velocity: Vec<f64>,
}

#[derive(Serialize)]

struct ParticleOutput {
    id: usize,
    mass: f64,
    position: Vec<f64>,
    velocity: Vec<f64>,
    kinetic_energy: f64,
    speed: f64,
}

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

struct MorseInput {
    p1_position: Vec<f64>,
    p2_position: Vec<f64>,
    de: f64,
    a: f64,
    re: f64,
}

#[derive(Deserialize)]

struct HarmonicInput {
    p1_position: Vec<f64>,
    p2_position: Vec<f64>,
    k: f64,
    r0: f64,
}

#[derive(Deserialize)]

struct ParticleListInput {
    particles: Vec<ParticleInput>,
}

#[derive(Serialize)]

struct SystemPropertiesOutput {
    kinetic_energy: f64,
    temperature: f64,
    center_of_mass: Vec<f64>,
    total_momentum: Vec<f64>,
}

#[derive(Deserialize)]

struct PbcInput {
    position: Vec<f64>,
    box_size: Vec<f64>,
}

#[derive(Deserialize)]

struct LatticeInput {
    n_per_side: usize,
    lattice_constant: f64,
    mass: f64,
}

// ============================================================================
// Potential Functions
// ============================================================================

/// Computes the Lennard-Jones interaction potential and force between two particles using JSON serialization.
///
/// The Lennard-Jones potential models van der Waals interactions:
/// V(r) = 4ε[(σ/r)¹² - (σ/r)⁶].
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `p1_position`: Position vector of first particle [x, y, z]
///   - `p2_position`: Position vector of second particle [x, y, z]
///   - `epsilon`: Well depth ε (energy units)
///   - `sigma`: Finite distance σ at which V=0 (length units)
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult` with potential and force.
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_md_lennard_jones_json(
    input: *const c_char
) -> *mut c_char {

    let input : LennardJonesInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<InteractionOutput, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
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

            to_c_string(
                serde_json::to_string(&FfiResult {
                    ok : Some(output),
                    err : None::<String>,
                })
                .unwrap(),
            )
        },
        | Err(e) => {
            to_c_string(
                serde_json::to_string(
                    &FfiResult::<InteractionOutput, String> {
                        ok : None,
                        err : Some(e),
                    },
                )
                .unwrap(),
            )
        },
    }
}

/// Computes the Morse interaction potential and force between two particles using JSON serialization.
///
/// The Morse potential models chemical bonds:
/// V(r) = `D_e`[1 - e⁻ᵃʳʳ⁻ʳᵉ⁾]², where `D_e` is the dissociation energy.
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `p1_position`: Position vector of first particle [x, y, z]
///   - `p2_position`: Position vector of second particle [x, y, z]
///   - `de`: Dissociation energy `D_e`
///   - `a`: Width parameter a (controls potential curvature)
///   - `re`: Equilibrium bond distance `r_e`
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult` with potential and force.
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_md_morse_json(
    input: *const c_char
) -> *mut c_char {

    let input : MorseInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<InteractionOutput, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
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

    match physics_md::morse_interaction(
        &p1,
        &p2,
        input.de,
        input.a,
        input.re,
    ) {
        | Ok((potential, force)) => {

            let output =
                InteractionOutput {
                    potential,
                    force,
                };

            to_c_string(
                serde_json::to_string(
                    &FfiResult {
                        ok: Some(
                            output,
                        ),
                        err: None::<
                            String,
                        >,
                    },
                )
                .unwrap(),
            )
        },
        | Err(e) => to_c_string(
            serde_json::to_string(
                &FfiResult::<
                    InteractionOutput,
                    String,
                > {
                    ok: None,
                    err: Some(e),
                },
            )
            .unwrap(),
        ),
    }
}

/// Computes the harmonic interaction potential and force between two particles using JSON serialization.
///
/// The harmonic potential models elastic bonds:
/// V(r) = ½k(r - r₀)², where k is the spring constant.
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `p1_position`: Position vector of first particle [x, y, z]
///   - `p2_position`: Position vector of second particle [x, y, z]
///   - `k`: Spring constant k (force/length units)
///   - `r0`: Equilibrium distance r₀
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult` with potential and force.
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_md_harmonic_json(
    input: *const c_char
) -> *mut c_char {

    let input : HarmonicInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<InteractionOutput, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
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

    match physics_md::harmonic_interaction(
        &p1,
        &p2,
        input.k,
        input.r0,
    ) {
        | Ok((potential, force)) => {

            let output = InteractionOutput {
                potential,
                force,
            };

            to_c_string(
                serde_json::to_string(&FfiResult {
                    ok : Some(output),
                    err : None::<String>,
                })
                .unwrap(),
            )
        },
        | Err(e) => {
            to_c_string(
                serde_json::to_string(
                    &FfiResult::<InteractionOutput, String> {
                        ok : None,
                        err : Some(e),
                    },
                )
                .unwrap(),
            )
        },
    }
}

// ============================================================================
// System Properties
// ============================================================================

/// Computes system-level properties for a collection of particles using JSON serialization.
///
/// Calculates kinetic energy, temperature, center of mass, and total momentum.
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `particles`: Array of particle objects with id, mass, position, velocity
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<SystemPropertiesOutput, String>` with:
/// - `kinetic_energy`: Total kinetic energy of the system
/// - `temperature`: System temperature T = `2K/(3Nk_B)`
/// - `center_of_mass`: Center of mass position vector
/// - `total_momentum`: Total linear momentum vector
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_md_system_properties_json(
    input: *const c_char
) -> *mut c_char {

    let input : ParticleListInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<SystemPropertiesOutput, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let particles: Vec<
        physics_md::Particle,
    > = input
        .particles
        .into_iter()
        .map(|p| {

            physics_md::Particle::new(
                p.id,
                p.mass,
                p.position,
                p.velocity,
            )
        })
        .collect();

    let ke = physics_md::total_kinetic_energy(&particles);

    let temp = physics_md::temperature(
        &particles,
    );

    let com =
        physics_md::center_of_mass(
            &particles,
        )
        .unwrap_or_default();

    let momentum =
        physics_md::total_momentum(
            &particles,
        )
        .unwrap_or_default();

    let output =
        SystemPropertiesOutput {
            kinetic_energy: ke,
            temperature: temp,
            center_of_mass: com,
            total_momentum: momentum,
        };

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok: Some(output),
                err: None::<String>,
            },
        )
        .unwrap(),
    )
}

// ============================================================================
// Lattice Creation
// ============================================================================

/// Creates a simple cubic lattice of particles using JSON serialization.
///
/// Generates a 3D cubic lattice structure commonly used for initial configurations
/// in molecular dynamics simulations.
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `n_per_side`: Number of particles per lattice dimension
///   - `lattice_constant`: Spacing between adjacent lattice sites
///   - `mass`: Mass of each particle
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<Vec<ParticleOutput>, String>` with
/// an array of particles positioned on a cubic lattice.
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_md_create_cubic_lattice_json(
    input: *const c_char
) -> *mut c_char {

    let input : LatticeInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<Vec<ParticleOutput>, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let particles = physics_md::create_cubic_lattice(
        input.n_per_side,
        input.lattice_constant,
        input.mass,
    );

    let output: Vec<ParticleOutput> =
        particles
            .iter()
            .map(|p| {

                ParticleOutput {
                    id: p.id,
                    mass: p.mass,
                    position: p
                        .position
                        .clone(),
                    velocity: p
                        .velocity
                        .clone(),
                    kinetic_energy: p
                        .kinetic_energy(
                        ),
                    speed: p.speed(),
                }
            })
            .collect();

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok: Some(output),
                err: None::<String>,
            },
        )
        .unwrap(),
    )
}

// ============================================================================
// Periodic Boundary Conditions
// ============================================================================

/// Applies periodic boundary conditions (PBC) to wrap particle coordinates using JSON serialization.
///
/// Periodic boundary conditions create an infinite tiling of the simulation box.
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `position`: Position vector to wrap [x, y, z]
///   - `box_size`: Simulation box dimensions [Lx, Ly, Lz]
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<Vec<f64>, String>` with
/// the wrapped position vector within [0, `box_size`).
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_md_apply_pbc_json(
    input: *const c_char
) -> *mut c_char {

    let input : PbcInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<Vec<f64>, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let wrapped = physics_md::apply_pbc(
        &input.position,
        &input.box_size,
    );

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok: Some(wrapped),
                err: None::<String>,
            },
        )
        .unwrap(),
    )
}

/// Computes the minimum image distance vector under periodic boundary conditions using JSON serialization.
///
/// The minimum image convention finds the shortest distance between particles
/// considering all periodic images of the simulation box.
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `position`: Separation vector [x, y, z]
///   - `box_size`: Simulation box dimensions [Lx, Ly, Lz]
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<Vec<f64>, String>` with
/// the minimum image distance vector.
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_md_minimum_image_json(
    input: *const c_char
) -> *mut c_char {

    let input : PbcInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<Vec<f64>, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let result = physics_md::minimum_image_distance(
        &input.position,
        &input.box_size,
    );

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok: Some(result),
                err: None::<String>,
            },
        )
        .unwrap(),
    )
}

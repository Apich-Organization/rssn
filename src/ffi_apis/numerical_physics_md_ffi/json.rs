//! JSON-based FFI API for numerical MD functions.

use crate::ffi_apis::common::{from_json_string, to_c_string};
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::physics_md;
use serde::{Deserialize, Serialize};
use std::os::raw::c_char;

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

#[no_mangle]
pub unsafe extern "C" fn rssn_num_md_lennard_jones_json(input: *const c_char) -> *mut c_char {
    let input: LennardJonesInput = match from_json_string(input) {
        Some(i) => i,
        None => {
            return to_c_string(
                serde_json::to_string(&FfiResult::<InteractionOutput, String> {
                    ok: None,
                    err: Some("Invalid JSON".to_string()),
                })
                .unwrap(),
            )
        }
    };

    let p1 = physics_md::Particle::new(0, 1.0, input.p1_position, vec![0.0, 0.0, 0.0]);
    let p2 = physics_md::Particle::new(1, 1.0, input.p2_position, vec![0.0, 0.0, 0.0]);

    match physics_md::lennard_jones_interaction(&p1, &p2, input.epsilon, input.sigma) {
        Ok((potential, force)) => {
            let output = InteractionOutput { potential, force };
            to_c_string(
                serde_json::to_string(&FfiResult {
                    ok: Some(output),
                    err: None::<String>,
                })
                .unwrap(),
            )
        }
        Err(e) => to_c_string(
            serde_json::to_string(&FfiResult::<InteractionOutput, String> {
                ok: None,
                err: Some(e),
            })
            .unwrap(),
        ),
    }
}

#[no_mangle]
pub unsafe extern "C" fn rssn_num_md_morse_json(input: *const c_char) -> *mut c_char {
    let input: MorseInput = match from_json_string(input) {
        Some(i) => i,
        None => {
            return to_c_string(
                serde_json::to_string(&FfiResult::<InteractionOutput, String> {
                    ok: None,
                    err: Some("Invalid JSON".to_string()),
                })
                .unwrap(),
            )
        }
    };

    let p1 = physics_md::Particle::new(0, 1.0, input.p1_position, vec![0.0, 0.0, 0.0]);
    let p2 = physics_md::Particle::new(1, 1.0, input.p2_position, vec![0.0, 0.0, 0.0]);

    match physics_md::morse_interaction(&p1, &p2, input.de, input.a, input.re) {
        Ok((potential, force)) => {
            let output = InteractionOutput { potential, force };
            to_c_string(
                serde_json::to_string(&FfiResult {
                    ok: Some(output),
                    err: None::<String>,
                })
                .unwrap(),
            )
        }
        Err(e) => to_c_string(
            serde_json::to_string(&FfiResult::<InteractionOutput, String> {
                ok: None,
                err: Some(e),
            })
            .unwrap(),
        ),
    }
}

#[no_mangle]
pub unsafe extern "C" fn rssn_num_md_harmonic_json(input: *const c_char) -> *mut c_char {
    let input: HarmonicInput = match from_json_string(input) {
        Some(i) => i,
        None => {
            return to_c_string(
                serde_json::to_string(&FfiResult::<InteractionOutput, String> {
                    ok: None,
                    err: Some("Invalid JSON".to_string()),
                })
                .unwrap(),
            )
        }
    };

    let p1 = physics_md::Particle::new(0, 1.0, input.p1_position, vec![0.0, 0.0, 0.0]);
    let p2 = physics_md::Particle::new(1, 1.0, input.p2_position, vec![0.0, 0.0, 0.0]);

    match physics_md::harmonic_interaction(&p1, &p2, input.k, input.r0) {
        Ok((potential, force)) => {
            let output = InteractionOutput { potential, force };
            to_c_string(
                serde_json::to_string(&FfiResult {
                    ok: Some(output),
                    err: None::<String>,
                })
                .unwrap(),
            )
        }
        Err(e) => to_c_string(
            serde_json::to_string(&FfiResult::<InteractionOutput, String> {
                ok: None,
                err: Some(e),
            })
            .unwrap(),
        ),
    }
}

// ============================================================================
// System Properties
// ============================================================================

#[no_mangle]
pub unsafe extern "C" fn rssn_num_md_system_properties_json(input: *const c_char) -> *mut c_char {
    let input: ParticleListInput = match from_json_string(input) {
        Some(i) => i,
        None => {
            return to_c_string(
                serde_json::to_string(&FfiResult::<SystemPropertiesOutput, String> {
                    ok: None,
                    err: Some("Invalid JSON".to_string()),
                })
                .unwrap(),
            )
        }
    };

    let particles: Vec<physics_md::Particle> = input
        .particles
        .into_iter()
        .map(|p| physics_md::Particle::new(p.id, p.mass, p.position, p.velocity))
        .collect();

    let ke = physics_md::total_kinetic_energy(&particles);
    let temp = physics_md::temperature(&particles);
    let com = physics_md::center_of_mass(&particles).unwrap_or_default();
    let momentum = physics_md::total_momentum(&particles).unwrap_or_default();

    let output = SystemPropertiesOutput {
        kinetic_energy: ke,
        temperature: temp,
        center_of_mass: com,
        total_momentum: momentum,
    };

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok: Some(output),
            err: None::<String>,
        })
        .unwrap(),
    )
}

// ============================================================================
// Lattice Creation
// ============================================================================

#[no_mangle]
pub unsafe extern "C" fn rssn_num_md_create_cubic_lattice_json(
    input: *const c_char,
) -> *mut c_char {
    let input: LatticeInput = match from_json_string(input) {
        Some(i) => i,
        None => {
            return to_c_string(
                serde_json::to_string(&FfiResult::<Vec<ParticleOutput>, String> {
                    ok: None,
                    err: Some("Invalid JSON".to_string()),
                })
                .unwrap(),
            )
        }
    };

    let particles =
        physics_md::create_cubic_lattice(input.n_per_side, input.lattice_constant, input.mass);

    let output: Vec<ParticleOutput> = particles
        .iter()
        .map(|p| ParticleOutput {
            id: p.id,
            mass: p.mass,
            position: p.position.clone(),
            velocity: p.velocity.clone(),
            kinetic_energy: p.kinetic_energy(),
            speed: p.speed(),
        })
        .collect();

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok: Some(output),
            err: None::<String>,
        })
        .unwrap(),
    )
}

// ============================================================================
// Periodic Boundary Conditions
// ============================================================================

#[no_mangle]
pub unsafe extern "C" fn rssn_num_md_apply_pbc_json(input: *const c_char) -> *mut c_char {
    let input: PbcInput = match from_json_string(input) {
        Some(i) => i,
        None => {
            return to_c_string(
                serde_json::to_string(&FfiResult::<Vec<f64>, String> {
                    ok: None,
                    err: Some("Invalid JSON".to_string()),
                })
                .unwrap(),
            )
        }
    };

    let wrapped = physics_md::apply_pbc(&input.position, &input.box_size);

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok: Some(wrapped),
            err: None::<String>,
        })
        .unwrap(),
    )
}

#[no_mangle]
pub unsafe extern "C" fn rssn_num_md_minimum_image_json(input: *const c_char) -> *mut c_char {
    let input: PbcInput = match from_json_string(input) {
        Some(i) => i,
        None => {
            return to_c_string(
                serde_json::to_string(&FfiResult::<Vec<f64>, String> {
                    ok: None,
                    err: Some("Invalid JSON".to_string()),
                })
                .unwrap(),
            )
        }
    };

    let result = physics_md::minimum_image_distance(&input.position, &input.box_size);

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok: Some(result),
            err: None::<String>,
        })
        .unwrap(),
    )
}

//! Bincode-based FFI API for numerical physics functions.

use crate::ffi_apis::common::{from_bincode_buffer, to_bincode_buffer, BincodeBuffer};
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::physics;
use serde::Deserialize;

#[derive(Deserialize)]
struct HarmonicOscillatorInput {
    amplitude: f64,
    omega: f64,
    phase: f64,
    time: f64,
}

#[derive(Deserialize)]
struct TwoChargesInput {
    q1: f64,
    q2: f64,
    r: f64,
}

#[derive(Deserialize)]
struct IdealGasInput {
    n: f64,
    t: f64,
    v: f64,
}

#[derive(Deserialize)]
struct VelocityInput {
    velocity: f64,
}

#[derive(Deserialize)]
struct MassInput {
    mass: f64,
}

#[derive(Deserialize)]
struct QuantumHarmonicInput {
    n: u64,
    omega: f64,
}

#[derive(Deserialize)]
struct QuantumNumberInput {
    n: u64,
}

#[no_mangle]
pub unsafe extern "C" fn rssn_num_physics_simple_harmonic_oscillator_bincode(
    buffer: BincodeBuffer,
) -> BincodeBuffer {
    let input: HarmonicOscillatorInput = match from_bincode_buffer(&buffer) {
        Some(i) => i,
        None => {
            return to_bincode_buffer(&FfiResult::<f64, String> {
                ok: None,
                err: Some("Invalid Bincode".to_string()),
            })
        }
    };
    let result = physics::simple_harmonic_oscillator(input.amplitude, input.omega, input.phase, input.time);
    to_bincode_buffer(&FfiResult {
        ok: Some(result),
        err: None::<String>,
    })
}

#[no_mangle]
pub unsafe extern "C" fn rssn_num_physics_coulomb_force_bincode(
    buffer: BincodeBuffer,
) -> BincodeBuffer {
    let input: TwoChargesInput = match from_bincode_buffer(&buffer) {
        Some(i) => i,
        None => {
            return to_bincode_buffer(&FfiResult::<f64, String> {
                ok: None,
                err: Some("Invalid Bincode".to_string()),
            })
        }
    };
    let result = physics::coulomb_force(input.q1, input.q2, input.r);
    to_bincode_buffer(&FfiResult {
        ok: Some(result),
        err: None::<String>,
    })
}

#[no_mangle]
pub unsafe extern "C" fn rssn_num_physics_ideal_gas_pressure_bincode(
    buffer: BincodeBuffer,
) -> BincodeBuffer {
    let input: IdealGasInput = match from_bincode_buffer(&buffer) {
        Some(i) => i,
        None => {
            return to_bincode_buffer(&FfiResult::<f64, String> {
                ok: None,
                err: Some("Invalid Bincode".to_string()),
            })
        }
    };
    let result = physics::ideal_gas_pressure(input.n, input.t, input.v);
    to_bincode_buffer(&FfiResult {
        ok: Some(result),
        err: None::<String>,
    })
}

#[no_mangle]
pub unsafe extern "C" fn rssn_num_physics_lorentz_factor_bincode(
    buffer: BincodeBuffer,
) -> BincodeBuffer {
    let input: VelocityInput = match from_bincode_buffer(&buffer) {
        Some(i) => i,
        None => {
            return to_bincode_buffer(&FfiResult::<f64, String> {
                ok: None,
                err: Some("Invalid Bincode".to_string()),
            })
        }
    };
    let result = physics::lorentz_factor(input.velocity);
    to_bincode_buffer(&FfiResult {
        ok: Some(result),
        err: None::<String>,
    })
}

#[no_mangle]
pub unsafe extern "C" fn rssn_num_physics_mass_energy_bincode(
    buffer: BincodeBuffer,
) -> BincodeBuffer {
    let input: MassInput = match from_bincode_buffer(&buffer) {
        Some(i) => i,
        None => {
            return to_bincode_buffer(&FfiResult::<f64, String> {
                ok: None,
                err: Some("Invalid Bincode".to_string()),
            })
        }
    };
    let result = physics::mass_energy(input.mass);
    to_bincode_buffer(&FfiResult {
        ok: Some(result),
        err: None::<String>,
    })
}

#[no_mangle]
pub unsafe extern "C" fn rssn_num_physics_quantum_harmonic_oscillator_energy_bincode(
    buffer: BincodeBuffer,
) -> BincodeBuffer {
    let input: QuantumHarmonicInput = match from_bincode_buffer(&buffer) {
        Some(i) => i,
        None => {
            return to_bincode_buffer(&FfiResult::<f64, String> {
                ok: None,
                err: Some("Invalid Bincode".to_string()),
            })
        }
    };
    let result = physics::quantum_harmonic_oscillator_energy(input.n, input.omega);
    to_bincode_buffer(&FfiResult {
        ok: Some(result),
        err: None::<String>,
    })
}

#[no_mangle]
pub unsafe extern "C" fn rssn_num_physics_hydrogen_energy_level_bincode(
    buffer: BincodeBuffer,
) -> BincodeBuffer {
    let input: QuantumNumberInput = match from_bincode_buffer(&buffer) {
        Some(i) => i,
        None => {
            return to_bincode_buffer(&FfiResult::<f64, String> {
                ok: None,
                err: Some("Invalid Bincode".to_string()),
            })
        }
    };
    let result = physics::hydrogen_energy_level(input.n);
    to_bincode_buffer(&FfiResult {
        ok: Some(result),
        err: None::<String>,
    })
}

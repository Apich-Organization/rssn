//! JSON-based FFI API for numerical physics functions.

use crate::ffi_apis::common::{
    from_json_string,
    to_c_string,
};
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::physics;
use serde::Deserialize;
use std::os::raw::c_char;

// ============================================================================
// Input structs
// ============================================================================

#[derive(Deserialize)]

struct HarmonicOscillatorInput {
    amplitude: f64,
    omega: f64,
    phase: f64,
    time: f64,
}

#[derive(Deserialize)]

struct DampedOscillatorInput {
    amplitude: f64,
    omega0: f64,
    gamma: f64,
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

struct PointChargeInput {
    q: f64,
    r: f64,
}

#[derive(Deserialize)]

struct IdealGasInput {
    n: f64,
    t: f64,
    v: f64,
}

#[derive(Deserialize)]

struct MassTempInput {
    mass: f64,
    temperature: f64,
}

#[derive(Deserialize)]

struct BlackbodyInput {
    area: f64,
    temperature: f64,
}

#[derive(Deserialize)]

struct VelocityInput {
    velocity: f64,
}

#[derive(Deserialize)]

struct TimeDilationInput {
    proper_time: f64,
    velocity: f64,
}

#[derive(Deserialize)]

struct VelocityAdditionInput {
    v: f64,
    w: f64,
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

#[derive(Deserialize)]

struct MomentumInput {
    momentum: f64,
}

#[derive(Deserialize)]

struct WavelengthInput {
    wavelength: f64,
}

#[derive(Deserialize)]

struct EnergyInput {
    energy: f64,
}

#[derive(Deserialize)]

struct MassInput {
    mass: f64,
}

#[derive(Deserialize)]

struct TemperatureInput {
    temperature: f64,
}

// ============================================================================
// Classical Mechanics
// ============================================================================

#[no_mangle]

pub unsafe extern "C" fn rssn_num_physics_simple_harmonic_oscillator_json(
    input: *const c_char
) -> *mut c_char {

    let input: HarmonicOscillatorInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok: None,
                        err: Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let result = physics::simple_harmonic_oscillator(
        input.amplitude,
        input.omega,
        input.phase,
        input.time,
    );

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok: Some(result),
            err: None::<String>,
        })
        .unwrap(),
    )
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_physics_damped_harmonic_oscillator_json(
    input: *const c_char
) -> *mut c_char {

    let input: DampedOscillatorInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok: None,
                        err: Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let result = physics::damped_harmonic_oscillator(
        input.amplitude,
        input.omega0,
        input.gamma,
        input.phase,
        input.time,
    );

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok: Some(result),
            err: None::<String>,
        })
        .unwrap(),
    )
}

// ============================================================================
// Electromagnetism
// ============================================================================

#[no_mangle]

pub unsafe extern "C" fn rssn_num_physics_coulomb_force_json(input: *const c_char) -> *mut c_char {

    let input: TwoChargesInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok: None,
                        err: Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let result = physics::coulomb_force(
        input.q1, input.q2, input.r,
    );

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok: Some(result),
            err: None::<String>,
        })
        .unwrap(),
    )
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_physics_electric_field_point_charge_json(
    input: *const c_char
) -> *mut c_char {

    let input: PointChargeInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok: None,
                        err: Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let result = physics::electric_field_point_charge(input.q, input.r);

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok: Some(result),
            err: None::<String>,
        })
        .unwrap(),
    )
}

// ============================================================================
// Thermodynamics
// ============================================================================

#[no_mangle]

pub unsafe extern "C" fn rssn_num_physics_ideal_gas_pressure_json(
    input: *const c_char
) -> *mut c_char {

    let input: IdealGasInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok: None,
                        err: Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let result = physics::ideal_gas_pressure(
        input.n, input.t, input.v,
    );

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok: Some(result),
            err: None::<String>,
        })
        .unwrap(),
    )
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_physics_maxwell_boltzmann_mean_speed_json(
    input: *const c_char
) -> *mut c_char {

    let input: MassTempInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok: None,
                        err: Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let result = physics::maxwell_boltzmann_mean_speed(
        input.mass,
        input.temperature,
    );

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok: Some(result),
            err: None::<String>,
        })
        .unwrap(),
    )
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_physics_blackbody_power_json(
    input: *const c_char
) -> *mut c_char {

    let input: BlackbodyInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok: None,
                        err: Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let result = physics::blackbody_power(
        input.area,
        input.temperature,
    );

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok: Some(result),
            err: None::<String>,
        })
        .unwrap(),
    )
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_physics_wien_displacement_wavelength_json(
    input: *const c_char
) -> *mut c_char {

    let input: TemperatureInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok: None,
                        err: Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let result = physics::wien_displacement_wavelength(input.temperature);

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok: Some(result),
            err: None::<String>,
        })
        .unwrap(),
    )
}

// ============================================================================
// Special Relativity
// ============================================================================

#[no_mangle]

pub unsafe extern "C" fn rssn_num_physics_lorentz_factor_json(input: *const c_char) -> *mut c_char {

    let input: VelocityInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok: None,
                        err: Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let result = physics::lorentz_factor(input.velocity);

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok: Some(result),
            err: None::<String>,
        })
        .unwrap(),
    )
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_physics_time_dilation_json(input: *const c_char) -> *mut c_char {

    let input: TimeDilationInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok: None,
                        err: Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let result = physics::time_dilation(
        input.proper_time,
        input.velocity,
    );

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok: Some(result),
            err: None::<String>,
        })
        .unwrap(),
    )
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_physics_mass_energy_json(input: *const c_char) -> *mut c_char {

    let input: MassInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok: None,
                        err: Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let result = physics::mass_energy(input.mass);

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok: Some(result),
            err: None::<String>,
        })
        .unwrap(),
    )
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_physics_relativistic_velocity_addition_json(
    input: *const c_char
) -> *mut c_char {

    let input: VelocityAdditionInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok: None,
                        err: Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let result = physics::relativistic_velocity_addition(input.v, input.w);

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok: Some(result),
            err: None::<String>,
        })
        .unwrap(),
    )
}

// ============================================================================
// Quantum Mechanics
// ============================================================================

#[no_mangle]

pub unsafe extern "C" fn rssn_num_physics_quantum_harmonic_oscillator_energy_json(
    input: *const c_char
) -> *mut c_char {

    let input: QuantumHarmonicInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok: None,
                        err: Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let result = physics::quantum_harmonic_oscillator_energy(input.n, input.omega);

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok: Some(result),
            err: None::<String>,
        })
        .unwrap(),
    )
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_physics_hydrogen_energy_level_json(
    input: *const c_char
) -> *mut c_char {

    let input: QuantumNumberInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok: None,
                        err: Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let result = physics::hydrogen_energy_level(input.n);

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok: Some(result),
            err: None::<String>,
        })
        .unwrap(),
    )
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_physics_de_broglie_wavelength_json(
    input: *const c_char
) -> *mut c_char {

    let input: MomentumInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok: None,
                        err: Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let result = physics::de_broglie_wavelength(input.momentum);

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok: Some(result),
            err: None::<String>,
        })
        .unwrap(),
    )
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_physics_photon_energy_json(input: *const c_char) -> *mut c_char {

    let input: WavelengthInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok: None,
                        err: Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let result = physics::photon_energy(input.wavelength);

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok: Some(result),
            err: None::<String>,
        })
        .unwrap(),
    )
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_physics_photon_wavelength_json(
    input: *const c_char
) -> *mut c_char {

    let input: EnergyInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok: None,
                        err: Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let result = physics::photon_wavelength(input.energy);

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok: Some(result),
            err: None::<String>,
        })
        .unwrap(),
    )
}

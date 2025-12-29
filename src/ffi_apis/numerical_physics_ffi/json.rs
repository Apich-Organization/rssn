//! JSON-based FFI API for numerical physics functions.

use std::os::raw::c_char;

use serde::Deserialize;

use crate::ffi_apis::common::from_json_string;
use crate::ffi_apis::common::to_c_string;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::physics;

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

/// Computes the displacement of a simple harmonic oscillator using JSON serialization.
///
/// The displacement is x(t) = A cos(ωt + φ).
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `amplitude`: Oscillation amplitude A (m)
///   - `omega`: Angular frequency ω (rad/s)
///   - `phase`: Phase angle φ (radians)
///   - `time`: Time t at which to evaluate displacement (s)
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<f64, String>` with
/// the displacement x(t) (m).
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

pub unsafe extern "C" fn rssn_num_physics_simple_harmonic_oscillator_json(
    input: *const c_char
) -> *mut c_char {

    let input : HarmonicOscillatorInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
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
        serde_json::to_string(
            &FfiResult {
                ok: Some(result),
                err: None::<String>,
            },
        )
        .unwrap(),
    )
}

/// Computes the displacement of a damped harmonic oscillator using JSON serialization.
///
/// The displacement is x(t) = A e⁻γᵗ cos(ω't + φ), where ω' = √(ω₀² - γ²).
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `amplitude`: Initial amplitude A (m)
///   - `omega0`: Natural angular frequency ω₀ (rad/s)
///   - `gamma`: Damping coefficient γ (1/s)
///   - `phase`: Phase angle φ (radians)
///   - `time`: Time t at which to evaluate displacement (s)
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<f64, String>` with
/// the displacement x(t) for underdamped oscillation (m).
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

pub unsafe extern "C" fn rssn_num_physics_damped_harmonic_oscillator_json(
    input: *const c_char
) -> *mut c_char {

    let input : DampedOscillatorInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
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
        serde_json::to_string(
            &FfiResult {
                ok: Some(result),
                err: None::<String>,
            },
        )
        .unwrap(),
    )
}

// ============================================================================
// Electromagnetism
// ============================================================================

/// Computes the Coulomb electrostatic force between two point charges using JSON serialization.
///
/// The force is F = `k_e` q₁q₂ / r², where `k_e` = 8.99 × 10⁹ N·m²/C².
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `q1`: First charge q₁ (Coulombs)
///   - `q2`: Second charge q₂ (Coulombs)
///   - `r`: Separation distance r (m)
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<f64, String>` with
/// the Coulomb force magnitude F (Newtons, positive for repulsion).
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

pub unsafe extern "C" fn rssn_num_physics_coulomb_force_json(
    input: *const c_char
) -> *mut c_char {

    let input : TwoChargesInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let result = physics::coulomb_force(
        input.q1,
        input.q2,
        input.r,
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

/// Computes the electric field magnitude of a point charge using JSON serialization.
///
/// The electric field is E = `k_e` q / r².
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `q`: Point charge q (Coulombs)
///   - `r`: Distance from charge r (m)
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<f64, String>` with
/// the electric field magnitude E (N/C or V/m).
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

pub unsafe extern "C" fn rssn_num_physics_electric_field_point_charge_json(
    input: *const c_char
) -> *mut c_char {

    let input : PointChargeInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let result = physics::electric_field_point_charge(input.q, input.r);

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

// ============================================================================
// Thermodynamics
// ============================================================================

/// Computes the pressure of an ideal gas using JSON serialization.
///
/// Uses the ideal gas law PV = nRT, where R = 8.314 J/(mol·K).
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `n`: Amount of substance (moles)
///   - `t`: Absolute temperature T (Kelvin)
///   - `v`: Volume V (m³)
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<f64, String>` with
/// the gas pressure P (Pascals).
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

pub unsafe extern "C" fn rssn_num_physics_ideal_gas_pressure_json(
    input: *const c_char
) -> *mut c_char {

    let input : IdealGasInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let result =
        physics::ideal_gas_pressure(
            input.n,
            input.t,
            input.v,
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

/// Computes the mean speed of particles in a Maxwell-Boltzmann distribution using JSON serialization.
///
/// The mean speed is ⟨v⟩ = √(`8k_BT/(πm)`).
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `mass`: Particle mass m (kg)
///   - `temperature`: Absolute temperature T (Kelvin)
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<f64, String>` with
/// the mean speed ⟨v⟩ (m/s).
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

pub unsafe extern "C" fn rssn_num_physics_maxwell_boltzmann_mean_speed_json(
    input: *const c_char
) -> *mut c_char {

    let input : MassTempInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
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
        serde_json::to_string(
            &FfiResult {
                ok: Some(result),
                err: None::<String>,
            },
        )
        .unwrap(),
    )
}

/// Computes the total power radiated by a blackbody using JSON serialization.
///
/// Uses the Stefan-Boltzmann law P = `σAT⁴`, where σ = 5.67 × 10⁻⁸ W/(m²·K⁴).
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `area`: Surface area A (m²)
///   - `temperature`: Absolute temperature T (Kelvin)
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<f64, String>` with
/// the radiated power P (Watts).
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

pub unsafe extern "C" fn rssn_num_physics_blackbody_power_json(
    input: *const c_char
) -> *mut c_char {

    let input : BlackbodyInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let result =
        physics::blackbody_power(
            input.area,
            input.temperature,
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

/// Computes the peak wavelength of blackbody radiation using Wien's displacement law and JSON serialization.
///
/// Wien's law states `λ_max` = b/T, where b = 2.898 × 10⁻³ m·K.
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `temperature`: Absolute temperature T (Kelvin)
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<f64, String>` with
/// the peak wavelength `λ_max` (meters).
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

pub unsafe extern "C" fn rssn_num_physics_wien_displacement_wavelength_json(
    input: *const c_char
) -> *mut c_char {

    let input : TemperatureInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let result = physics::wien_displacement_wavelength(input.temperature);

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

// ============================================================================
// Special Relativity
// ============================================================================

/// Computes the Lorentz factor for relativistic transformations using JSON serialization.
///
/// The Lorentz factor is γ = 1 / √(1 - v²/c²).
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `velocity`: Velocity v (m/s)
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<f64, String>` with
/// the Lorentz factor γ (dimensionless, ≥ 1).
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

pub unsafe extern "C" fn rssn_num_physics_lorentz_factor_json(
    input: *const c_char
) -> *mut c_char {

    let input : VelocityInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let result =
        physics::lorentz_factor(
            input.velocity,
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

/// Computes relativistic time dilation using JSON serialization.
///
/// The dilated time is t = γt₀, where γ is the Lorentz factor.
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `proper_time`: Proper time t₀ (s)
///   - `velocity`: Relative velocity v (m/s)
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<f64, String>` with
/// the dilated time t (s).
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

pub unsafe extern "C" fn rssn_num_physics_time_dilation_json(
    input: *const c_char
) -> *mut c_char {

    let input : TimeDilationInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
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
        serde_json::to_string(
            &FfiResult {
                ok: Some(result),
                err: None::<String>,
            },
        )
        .unwrap(),
    )
}

/// Computes rest mass energy using Einstein's mass-energy equivalence and JSON serialization.
///
/// The rest energy is E = mc².
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `mass`: Rest mass m (kg)
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<f64, String>` with
/// the rest energy E (Joules).
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

pub unsafe extern "C" fn rssn_num_physics_mass_energy_json(
    input: *const c_char
) -> *mut c_char {

    let input : MassInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let result = physics::mass_energy(
        input.mass,
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

/// Computes relativistic velocity addition using JSON serialization.
///
/// The combined velocity is u = (v + w) / (1 + vw/c²).
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `v`: First velocity (m/s)
///   - `w`: Second velocity (m/s)
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<f64, String>` with
/// the combined velocity u (m/s).
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

pub unsafe extern "C" fn rssn_num_physics_relativistic_velocity_addition_json(
    input: *const c_char
) -> *mut c_char {

    let input : VelocityAdditionInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let result = physics::relativistic_velocity_addition(input.v, input.w);

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

// ============================================================================
// Quantum Mechanics
// ============================================================================

/// Computes the energy eigenvalue of a quantum harmonic oscillator using JSON serialization.
///
/// The energy is `E_n` = ℏω(n + 1/2), where ℏ is the reduced Planck constant.
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `n`: Quantum number n (non-negative integer, ground state = 0)
///   - `omega`: Angular frequency ω (rad/s)
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<f64, String>` with
/// the energy eigenvalue `E_n` (Joules).
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

pub unsafe extern "C" fn rssn_num_physics_quantum_harmonic_oscillator_energy_json(
    input: *const c_char
) -> *mut c_char {

    let input : QuantumHarmonicInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let result = physics::quantum_harmonic_oscillator_energy(input.n, input.omega);

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

/// Computes the energy level of the hydrogen atom using the Bohr model and JSON serialization.
///
/// The energy is `E_n` = -13.6 eV / n², where n is the principal quantum number.
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `n`: Principal quantum number n (positive integer, ground state = 1)
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<f64, String>` with
/// the energy level `E_n` (Joules, negative for bound states).
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

pub unsafe extern "C" fn rssn_num_physics_hydrogen_energy_level_json(
    input: *const c_char
) -> *mut c_char {

    let input : QuantumNumberInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let result =
        physics::hydrogen_energy_level(
            input.n,
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

/// Computes the de Broglie wavelength of a particle using JSON serialization.
///
/// The wavelength is λ = h/p, where h is Planck's constant and p is momentum.
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `momentum`: Particle momentum p (kg·m/s)
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<f64, String>` with
/// the de Broglie wavelength λ (meters).
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

pub unsafe extern "C" fn rssn_num_physics_de_broglie_wavelength_json(
    input: *const c_char
) -> *mut c_char {

    let input : MomentumInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let result =
        physics::de_broglie_wavelength(
            input.momentum,
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

/// Computes the energy of a photon from its wavelength using JSON serialization.
///
/// The energy is E = hc/λ.
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `wavelength`: Photon wavelength λ (meters)
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<f64, String>` with
/// the photon energy E (Joules).
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

pub unsafe extern "C" fn rssn_num_physics_photon_energy_json(
    input: *const c_char
) -> *mut c_char {

    let input : WavelengthInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let result = physics::photon_energy(
        input.wavelength,
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

/// Computes the wavelength of a photon from its energy using JSON serialization.
///
/// The wavelength is λ = hc/E.
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `energy`: Photon energy E (Joules)
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<f64, String>` with
/// the photon wavelength λ (meters).
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

pub unsafe extern "C" fn rssn_num_physics_photon_wavelength_json(
    input: *const c_char
) -> *mut c_char {

    let input : EnergyInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let result =
        physics::photon_wavelength(
            input.energy,
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

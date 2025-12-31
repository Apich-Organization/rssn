//! Bincode-based FFI API for numerical physics functions.

use serde::Deserialize;

use crate::ffi_apis::common::BincodeBuffer;
use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::physics;

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

/// Computes the displacement of a simple harmonic oscillator using bincode serialization.
///
/// The simple harmonic oscillator describes oscillatory motion with displacement:
/// x(t) = A cos(ωt + φ), where A is amplitude, ω is angular frequency, and φ is phase.
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing `HarmonicOscillatorInput` with:
///   - `amplitude`: Oscillation amplitude A (m)
///   - `omega`: Angular frequency ω (rad/s)
///   - `phase`: Phase angle φ (radians)
///   - `time`: Time t at which to evaluate displacement (s)
///
/// # Returns
///
/// A bincode-encoded buffer containing `FfiResult<f64, String>` with either:
/// - `ok`: The displacement x(t) (m)
/// - `err`: Error message if deserialization failed
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

pub unsafe extern "C" fn rssn_num_physics_simple_harmonic_oscillator_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input : HarmonicOscillatorInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(
                &FfiResult::<f64, String> {
                    ok : None,
                    err : Some("Invalid Bincode".to_string()),
                },
            )
        },
    };

    let result = physics::simple_harmonic_oscillator(
        input.amplitude,
        input.omega,
        input.phase,
        input.time,
    );

    to_bincode_buffer(&FfiResult {
        ok: Some(result),
        err: None::<String>,
    })
}

/// Computes the Coulomb electrostatic force between two point charges using bincode serialization.
///
/// The Coulomb force is the electrostatic interaction between charged particles:
/// F = `k_e` q₁q₂ / r², where `k_e` is Coulomb's constant (8.99×10⁹ N·m²/C²).
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing `TwoChargesInput` with:
///   - `q1`: First charge q₁ (Coulombs)
///   - `q2`: Second charge q₂ (Coulombs)
///   - `r`: Separation distance r (m)
///
/// # Returns
///
/// A bincode-encoded buffer containing `FfiResult<f64, String>` with either:
/// - `ok`: The Coulomb force magnitude F (Newtons, positive for repulsion)
/// - `err`: Error message if deserialization failed
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

pub unsafe extern "C" fn rssn_num_physics_coulomb_force_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input : TwoChargesInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(
                &FfiResult::<f64, String> {
                    ok : None,
                    err : Some("Invalid Bincode".to_string()),
                },
            )
        },
    };

    let result = physics::coulomb_force(
        input.q1,
        input.q2,
        input.r,
    );

    to_bincode_buffer(&FfiResult {
        ok: Some(result),
        err: None::<String>,
    })
}

/// Computes the pressure of an ideal gas using the ideal gas law and bincode serialization.
///
/// The ideal gas law relates pressure, volume, temperature, and amount of substance:
/// PV = nRT, where R is the universal gas constant (8.314 J/(mol·K)).
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing `IdealGasInput` with:
///   - `n`: Amount of substance (moles)
///   - `t`: Absolute temperature T (Kelvin)
///   - `v`: Volume V (m³)
///
/// # Returns
///
/// A bincode-encoded buffer containing `FfiResult<f64, String>` with either:
/// - `ok`: The gas pressure P (Pascals)
/// - `err`: Error message if deserialization failed
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

pub unsafe extern "C" fn rssn_num_physics_ideal_gas_pressure_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input : IdealGasInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(
                &FfiResult::<f64, String> {
                    ok : None,
                    err : Some("Invalid Bincode".to_string()),
                },
            )
        },
    };

    let result =
        physics::ideal_gas_pressure(
            input.n,
            input.t,
            input.v,
        );

    to_bincode_buffer(&FfiResult {
        ok: Some(result),
        err: None::<String>,
    })
}

/// Computes the Lorentz factor for relativistic time dilation and length contraction using bincode serialization.
///
/// The Lorentz factor γ appears in special relativity transformations:
/// γ = 1 / √(1 - v²/c²), where v is velocity and c is the speed of light.
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing `VelocityInput` with:
///   - `velocity`: Velocity v (m/s, typically as a fraction of c)
///
/// # Returns
///
/// A bincode-encoded buffer containing `FfiResult<f64, String>` with either:
/// - `ok`: The Lorentz factor γ (dimensionless, ≥ 1)
/// - `err`: Error message if deserialization failed
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

pub unsafe extern "C" fn rssn_num_physics_lorentz_factor_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input : VelocityInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(
                &FfiResult::<f64, String> {
                    ok : None,
                    err : Some("Invalid Bincode".to_string()),
                },
            )
        },
    };

    let result =
        physics::lorentz_factor(
            input.velocity,
        );

    to_bincode_buffer(&FfiResult {
        ok: Some(result),
        err: None::<String>,
    })
}

/// Computes the rest mass energy using Einstein's mass-energy equivalence and bincode serialization.
///
/// Einstein's mass-energy relation is one of the most famous equations in physics:
/// E = mc², where m is rest mass and c is the speed of light (2.998×10⁸ m/s).
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing `MassInput` with:
///   - `mass`: Rest mass m (kg)
///
/// # Returns
///
/// A bincode-encoded buffer containing `FfiResult<f64, String>` with either:
/// - `ok`: The rest energy E (Joules)
/// - `err`: Error message if deserialization failed
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

pub unsafe extern "C" fn rssn_num_physics_mass_energy_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input : MassInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(
                &FfiResult::<f64, String> {
                    ok : None,
                    err : Some("Invalid Bincode".to_string()),
                },
            )
        },
    };

    let result = physics::mass_energy(
        input.mass,
    );

    to_bincode_buffer(&FfiResult {
        ok: Some(result),
        err: None::<String>,
    })
}

/// Computes the energy eigenvalue of a quantum harmonic oscillator using bincode serialization.
///
/// The quantum harmonic oscillator has discrete energy levels given by:
/// `E_n` = ℏω(n + 1/2), where n is the quantum number, ℏ is the reduced Planck constant,
/// and ω is the angular frequency.
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing `QuantumHarmonicInput` with:
///   - `n`: Quantum number n (non-negative integer, ground state = 0)
///   - `omega`: Angular frequency ω (rad/s)
///
/// # Returns
///
/// A bincode-encoded buffer containing `FfiResult<f64, String>` with either:
/// - `ok`: The energy eigenvalue `E_n` (Joules)
/// - `err`: Error message if deserialization failed
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

pub unsafe extern "C" fn rssn_num_physics_quantum_harmonic_oscillator_energy_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input : QuantumHarmonicInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(
                &FfiResult::<f64, String> {
                    ok : None,
                    err : Some("Invalid Bincode".to_string()),
                },
            )
        },
    };

    let result = physics::quantum_harmonic_oscillator_energy(input.n, input.omega);

    to_bincode_buffer(&FfiResult {
        ok: Some(result),
        err: None::<String>,
    })
}

/// Computes the energy level of the hydrogen atom using the Bohr model and bincode serialization.
///
/// The Bohr model gives hydrogen atom energy levels as:
/// `E_n` = -13.6 eV / n², where n is the principal quantum number (n ≥ 1).
/// Energy is negative, indicating a bound state.
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing `QuantumNumberInput` with:
///   - `n`: Principal quantum number n (positive integer, ground state = 1)
///
/// # Returns
///
/// A bincode-encoded buffer containing `FfiResult<f64, String>` with either:
/// - `ok`: The energy level `E_n` (Joules, negative for bound states)
/// - `err`: Error message if deserialization failed
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

pub unsafe extern "C" fn rssn_num_physics_hydrogen_energy_level_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input : QuantumNumberInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(
                &FfiResult::<f64, String> {
                    ok : None,
                    err : Some("Invalid Bincode".to_string()),
                },
            )
        },
    };

    let result =
        physics::hydrogen_energy_level(
            input.n,
        );

    to_bincode_buffer(&FfiResult {
        ok: Some(result),
        err: None::<String>,
    })
}

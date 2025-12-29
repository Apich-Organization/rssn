//! Handle-based FFI API for numerical physics functions.

use crate::numerical::physics;

// ============================================================================
// Physical Constants
// ============================================================================

/// Returns the speed of light constant c via handle-based FFI.
///
/// # Returns
///
/// The speed of light c = 2.998 × 10⁸ m/s.
#[no_mangle]

pub extern "C" fn rssn_num_physics_speed_of_light(
) -> f64 {

    physics::SPEED_OF_LIGHT
}

/// Returns Planck's constant h via handle-based FFI.
///
/// # Returns
///
/// Planck's constant h = 6.626 × 10⁻³⁴ J·s.
#[no_mangle]

pub extern "C" fn rssn_num_physics_planck_constant(
) -> f64 {

    physics::PLANCK_CONSTANT
}

/// Returns the gravitational constant G via handle-based FFI.
///
/// # Returns
///
/// The gravitational constant G = 6.674 × 10⁻¹¹ N·m²/kg².
#[no_mangle]

pub extern "C" fn rssn_num_physics_gravitational_constant(
) -> f64 {

    physics::GRAVITATIONAL_CONSTANT
}

/// Returns the Boltzmann constant k_B via handle-based FFI.
///
/// # Returns
///
/// The Boltzmann constant k_B = 1.381 × 10⁻²³ J/K.
#[no_mangle]

pub extern "C" fn rssn_num_physics_boltzmann_constant(
) -> f64 {

    physics::BOLTZMANN_CONSTANT
}

/// Returns the elementary charge e via handle-based FFI.
///
/// # Returns
///
/// The elementary charge e = 1.602 × 10⁻¹⁹ Coulombs.
#[no_mangle]

pub extern "C" fn rssn_num_physics_elementary_charge(
) -> f64 {

    physics::ELEMENTARY_CHARGE
}

/// Returns the electron rest mass m_e via handle-based FFI.
///
/// # Returns
///
/// The electron rest mass m_e = 9.109 × 10⁻³¹ kg.
#[no_mangle]

pub extern "C" fn rssn_num_physics_electron_mass(
) -> f64 {

    physics::ELECTRON_MASS
}

// ============================================================================
// Classical Mechanics
// ============================================================================

/// Computes the displacement of a simple harmonic oscillator via handle-based FFI.
///
/// The displacement is given by x(t) = A cos(ωt + φ).
///
/// # Arguments
///
/// * `amplitude` - Oscillation amplitude A (m)
/// * `omega` - Angular frequency ω (rad/s)
/// * `phase` - Phase angle φ (radians)
/// * `time` - Time t at which to evaluate displacement (s)
///
/// # Returns
///
/// The displacement x(t) (m).
#[no_mangle]

pub extern "C" fn rssn_num_physics_simple_harmonic_oscillator(
    amplitude: f64,
    omega: f64,
    phase: f64,
    time: f64,
) -> f64 {

    physics::simple_harmonic_oscillator(
        amplitude,
        omega,
        phase,
        time,
    )
}

/// Computes the displacement of a damped harmonic oscillator via handle-based FFI.
///
/// The displacement is x(t) = A e⁻γᵗ cos(ω't + φ), where ω' = √(ω₀² - γ²).
///
/// # Arguments
///
/// * `amplitude` - Initial amplitude A (m)
/// * `omega0` - Natural angular frequency ω₀ (rad/s)
/// * `gamma` - Damping coefficient γ (1/s)
/// * `phase` - Phase angle φ (radians)
/// * `time` - Time t at which to evaluate displacement (s)
///
/// # Returns
///
/// The displacement x(t) for underdamped oscillation (m).
#[no_mangle]

pub extern "C" fn rssn_num_physics_damped_harmonic_oscillator(
    amplitude: f64,
    omega0: f64,
    gamma: f64,
    phase: f64,
    time: f64,
) -> f64 {

    physics::damped_harmonic_oscillator(
        amplitude,
        omega0,
        gamma,
        phase,
        time,
    )
}

// ============================================================================
// Electromagnetism
// ============================================================================

/// Computes the Coulomb electrostatic force between two point charges via handle-based FFI.
///
/// The force is F = k_e q₁q₂ / r², where k_e = 8.99 × 10⁹ N·m²/C².
///
/// # Arguments
///
/// * `q1` - First charge q₁ (Coulombs)
/// * `q2` - Second charge q₂ (Coulombs)
/// * `r` - Separation distance r (m)
///
/// # Returns
///
/// The Coulomb force magnitude F (Newtons, positive for repulsion).
#[no_mangle]

pub extern "C" fn rssn_num_physics_coulomb_force(
    q1: f64,
    q2: f64,
    r: f64,
) -> f64 {

    physics::coulomb_force(q1, q2, r)
}

/// Computes the electric field magnitude of a point charge via handle-based FFI.
///
/// The electric field is E = k_e q / r².
///
/// # Arguments
///
/// * `q` - Point charge q (Coulombs)
/// * `r` - Distance from charge r (m)
///
/// # Returns
///
/// The electric field magnitude E (N/C or V/m).
#[no_mangle]

pub extern "C" fn rssn_num_physics_electric_field_point_charge(
    q: f64,
    r: f64,
) -> f64 {

    physics::electric_field_point_charge(
        q, r,
    )
}

/// Computes the electric potential of a point charge via handle-based FFI.
///
/// The potential is V = k_e q / r.
///
/// # Arguments
///
/// * `q` - Point charge q (Coulombs)
/// * `r` - Distance from charge r (m)
///
/// # Returns
///
/// The electric potential V (Volts).
#[no_mangle]

pub extern "C" fn rssn_num_physics_electric_potential_point_charge(
    q: f64,
    r: f64,
) -> f64 {

    physics::electric_potential_point_charge(q, r)
}

/// Computes the magnetic field magnitude around an infinite straight current-carrying wire via handle-based FFI.
///
/// The magnetic field is B = (μ₀ I) / (2πr), where μ₀ is the permeability of free space.
///
/// # Arguments
///
/// * `current` - Electric current I (Amperes)
/// * `r` - Perpendicular distance from wire r (m)
///
/// # Returns
///
/// The magnetic field magnitude B (Tesla).
#[no_mangle]

pub extern "C" fn rssn_num_physics_magnetic_field_infinite_wire(
    current: f64,
    r: f64,
) -> f64 {

    physics::magnetic_field_infinite_wire(current, r)
}

/// Computes the Lorentz force on a charged particle in electromagnetic fields via handle-based FFI.
///
/// The force is F = q(E + v × B), simplified for parallel fields.
///
/// # Arguments
///
/// * `charge` - Particle charge q (Coulombs)
/// * `velocity` - Particle velocity v (m/s)
/// * `e_field` - Electric field magnitude E (V/m)
/// * `b_field` - Magnetic field magnitude B (Tesla)
///
/// # Returns
///
/// The Lorentz force magnitude F (Newtons).
#[no_mangle]

pub extern "C" fn rssn_num_physics_lorentz_force(
    charge: f64,
    velocity: f64,
    e_field: f64,
    b_field: f64,
) -> f64 {

    physics::lorentz_force(
        charge,
        velocity,
        e_field,
        b_field,
    )
}

/// Computes the cyclotron radius for a charged particle in a magnetic field via handle-based FFI.
///
/// The radius is r = (mv) / (qB), representing circular motion radius in a uniform B-field.
///
/// # Arguments
///
/// * `mass` - Particle mass m (kg)
/// * `velocity` - Particle velocity v (m/s)
/// * `charge` - Particle charge q (Coulombs)
/// * `b_field` - Magnetic field magnitude B (Tesla)
///
/// # Returns
///
/// The cyclotron radius r (m).
#[no_mangle]

pub extern "C" fn rssn_num_physics_cyclotron_radius(
    mass: f64,
    velocity: f64,
    charge: f64,
    b_field: f64,
) -> f64 {

    physics::cyclotron_radius(
        mass,
        velocity,
        charge,
        b_field,
    )
}

// ============================================================================
// Thermodynamics
// ============================================================================

/// Computes the pressure of an ideal gas via handle-based FFI.
///
/// Uses the ideal gas law PV = nRT, where R = 8.314 J/(mol·K).
///
/// # Arguments
///
/// * `n` - Amount of substance (moles)
/// * `t` - Absolute temperature T (Kelvin)
/// * `v` - Volume V (m³)
///
/// # Returns
///
/// The gas pressure P (Pascals).
#[no_mangle]

pub extern "C" fn rssn_num_physics_ideal_gas_pressure(
    n: f64,
    t: f64,
    v: f64,
) -> f64 {

    physics::ideal_gas_pressure(n, t, v)
}

/// Computes the volume of an ideal gas via handle-based FFI.
///
/// Uses the ideal gas law V = nRT/P.
///
/// # Arguments
///
/// * `n` - Amount of substance (moles)
/// * `t` - Absolute temperature T (Kelvin)
/// * `p` - Pressure P (Pascals)
///
/// # Returns
///
/// The gas volume V (m³).
#[no_mangle]

pub extern "C" fn rssn_num_physics_ideal_gas_volume(
    n: f64,
    t: f64,
    p: f64,
) -> f64 {

    physics::ideal_gas_volume(n, t, p)
}

/// Computes the temperature of an ideal gas via handle-based FFI.
///
/// Uses the ideal gas law T = PV/(nR).
///
/// # Arguments
///
/// * `p` - Pressure P (Pascals)
/// * `v` - Volume V (m³)
/// * `n` - Amount of substance (moles)
///
/// # Returns
///
/// The absolute temperature T (Kelvin).
#[no_mangle]

pub extern "C" fn rssn_num_physics_ideal_gas_temperature(
    p: f64,
    v: f64,
    n: f64,
) -> f64 {

    physics::ideal_gas_temperature(
        p, v, n,
    )
}

/// Computes the Maxwell-Boltzmann speed distribution probability density via handle-based FFI.
///
/// The distribution is f(v) = 4π(m/(2πk_BT))³¹² v² exp(-mv²/(2k_BT)).
///
/// # Arguments
///
/// * `v` - Particle speed v (m/s)
/// * `mass` - Particle mass m (kg)
/// * `temperature` - Absolute temperature T (Kelvin)
///
/// # Returns
///
/// The probability density f(v) (s/m).
#[no_mangle]

pub extern "C" fn rssn_num_physics_maxwell_boltzmann_speed_distribution(
    v: f64,
    mass: f64,
    temperature: f64,
) -> f64 {

    physics::maxwell_boltzmann_speed_distribution(v, mass, temperature)
}

/// Computes the mean speed of particles in a Maxwell-Boltzmann distribution via handle-based FFI.
///
/// The mean speed is ⟨v⟩ = √(8k_BT/(πm)).
///
/// # Arguments
///
/// * `mass` - Particle mass m (kg)
/// * `temperature` - Absolute temperature T (Kelvin)
///
/// # Returns
///
/// The mean speed ⟨v⟩ (m/s).
#[no_mangle]

pub extern "C" fn rssn_num_physics_maxwell_boltzmann_mean_speed(
    mass: f64,
    temperature: f64,
) -> f64 {

    physics::maxwell_boltzmann_mean_speed(mass, temperature)
}

/// Computes the root-mean-square speed of particles in a Maxwell-Boltzmann distribution via handle-based FFI.
///
/// The RMS speed is v_rms = √(3k_BT/m).
///
/// # Arguments
///
/// * `mass` - Particle mass m (kg)
/// * `temperature` - Absolute temperature T (Kelvin)
///
/// # Returns
///
/// The RMS speed v_rms (m/s).
#[no_mangle]

pub extern "C" fn rssn_num_physics_maxwell_boltzmann_rms_speed(
    mass: f64,
    temperature: f64,
) -> f64 {

    physics::maxwell_boltzmann_rms_speed(
        mass,
        temperature,
    )
}

/// Computes the total power radiated by a blackbody via handle-based FFI.
///
/// Uses the Stefan-Boltzmann law P = σAT⁴, where σ = 5.67 × 10⁻⁸ W/(m²·K⁴).
///
/// # Arguments
///
/// * `area` - Surface area A (m²)
/// * `temperature` - Absolute temperature T (Kelvin)
///
/// # Returns
///
/// The radiated power P (Watts).
#[no_mangle]

pub extern "C" fn rssn_num_physics_blackbody_power(
    area: f64,
    temperature: f64,
) -> f64 {

    physics::blackbody_power(
        area,
        temperature,
    )
}

/// Computes the peak wavelength of blackbody radiation using Wien's displacement law via handle-based FFI.
///
/// Wien's law states λ_max = b/T, where b = 2.898 × 10⁻³ m·K.
///
/// # Arguments
///
/// * `temperature` - Absolute temperature T (Kelvin)
///
/// # Returns
///
/// The peak wavelength λ_max (meters).
#[no_mangle]

pub extern "C" fn rssn_num_physics_wien_displacement_wavelength(
    temperature: f64
) -> f64 {

    physics::wien_displacement_wavelength(temperature)
}

// ============================================================================
// Special Relativity
// ============================================================================

/// Computes the Lorentz factor for relativistic transformations via handle-based FFI.
///
/// The Lorentz factor is γ = 1 / √(1 - v²/c²).
///
/// # Arguments
///
/// * `velocity` - Velocity v (m/s)
///
/// # Returns
///
/// The Lorentz factor γ (dimensionless, ≥ 1).
#[no_mangle]

pub extern "C" fn rssn_num_physics_lorentz_factor(
    velocity: f64
) -> f64 {

    physics::lorentz_factor(velocity)
}

/// Computes relativistic time dilation via handle-based FFI.
///
/// The dilated time is t = γt₀, where γ is the Lorentz factor.
///
/// # Arguments
///
/// * `proper_time` - Proper time t₀ (s)
/// * `velocity` - Relative velocity v (m/s)
///
/// # Returns
///
/// The dilated time t (s).
#[no_mangle]

pub extern "C" fn rssn_num_physics_time_dilation(
    proper_time: f64,
    velocity: f64,
) -> f64 {

    physics::time_dilation(
        proper_time,
        velocity,
    )
}

/// Computes relativistic length contraction via handle-based FFI.
///
/// The contracted length is L = L₀/γ, where γ is the Lorentz factor.
///
/// # Arguments
///
/// * `proper_length` - Proper length L₀ (m)
/// * `velocity` - Relative velocity v (m/s)
///
/// # Returns
///
/// The contracted length L (m).
#[no_mangle]

pub extern "C" fn rssn_num_physics_length_contraction(
    proper_length: f64,
    velocity: f64,
) -> f64 {

    physics::length_contraction(
        proper_length,
        velocity,
    )
}

/// Computes relativistic momentum via handle-based FFI.
///
/// The momentum is p = γmv, where γ is the Lorentz factor.
///
/// # Arguments
///
/// * `mass` - Rest mass m (kg)
/// * `velocity` - Velocity v (m/s)
///
/// # Returns
///
/// The relativistic momentum p (kg·m/s).
#[no_mangle]

pub extern "C" fn rssn_num_physics_relativistic_momentum(
    mass: f64,
    velocity: f64,
) -> f64 {

    physics::relativistic_momentum(
        mass,
        velocity,
    )
}

/// Computes relativistic kinetic energy via handle-based FFI.
///
/// The kinetic energy is K = (γ - 1)mc².
///
/// # Arguments
///
/// * `mass` - Rest mass m (kg)
/// * `velocity` - Velocity v (m/s)
///
/// # Returns
///
/// The kinetic energy K (Joules).
#[no_mangle]

pub extern "C" fn rssn_num_physics_relativistic_kinetic_energy(
    mass: f64,
    velocity: f64,
) -> f64 {

    physics::relativistic_kinetic_energy(
        mass,
        velocity,
    )
}

/// Computes rest mass energy using Einstein's mass-energy equivalence via handle-based FFI.
///
/// The rest energy is E = mc².
///
/// # Arguments
///
/// * `mass` - Rest mass m (kg)
///
/// # Returns
///
/// The rest energy E (Joules).
#[no_mangle]

pub extern "C" fn rssn_num_physics_mass_energy(
    mass: f64
) -> f64 {

    physics::mass_energy(mass)
}

/// Computes relativistic velocity addition via handle-based FFI.
///
/// The combined velocity is u = (v + w) / (1 + vw/c²).
///
/// # Arguments
///
/// * `v` - First velocity (m/s)
/// * `w` - Second velocity (m/s)
///
/// # Returns
///
/// The combined velocity u (m/s).
#[no_mangle]

pub extern "C" fn rssn_num_physics_relativistic_velocity_addition(
    v: f64,
    w: f64,
) -> f64 {

    physics::relativistic_velocity_addition(v, w)
}

// ============================================================================
// Quantum Mechanics
// ============================================================================

/// Computes the energy eigenvalue of a quantum harmonic oscillator via handle-based FFI.
///
/// The energy is E_n = ℏω(n + 1/2), where ℏ is the reduced Planck constant.
///
/// # Arguments
///
/// * `n` - Quantum number n (non-negative integer, ground state = 0)
/// * `omega` - Angular frequency ω (rad/s)
///
/// # Returns
///
/// The energy eigenvalue E_n (Joules).
#[no_mangle]

pub extern "C" fn rssn_num_physics_quantum_harmonic_oscillator_energy(
    n: u64,
    omega: f64,
) -> f64 {

    physics::quantum_harmonic_oscillator_energy(n, omega)
}

/// Computes the energy level of the hydrogen atom using the Bohr model via handle-based FFI.
///
/// The energy is E_n = -13.6 eV / n², where n is the principal quantum number.
///
/// # Arguments
///
/// * `n` - Principal quantum number n (positive integer, ground state = 1)
///
/// # Returns
///
/// The energy level E_n (Joules, negative for bound states).
#[no_mangle]

pub extern "C" fn rssn_num_physics_hydrogen_energy_level(
    n: u64
) -> f64 {

    physics::hydrogen_energy_level(n)
}

/// Computes the de Broglie wavelength of a particle via handle-based FFI.
///
/// The wavelength is λ = h/p, where h is Planck's constant and p is momentum.
///
/// # Arguments
///
/// * `momentum` - Particle momentum p (kg·m/s)
///
/// # Returns
///
/// The de Broglie wavelength λ (meters).
#[no_mangle]

pub extern "C" fn rssn_num_physics_de_broglie_wavelength(
    momentum: f64
) -> f64 {

    physics::de_broglie_wavelength(
        momentum,
    )
}

/// Computes the energy of a photon from its wavelength via handle-based FFI.
///
/// The energy is E = hc/λ.
///
/// # Arguments
///
/// * `wavelength` - Photon wavelength λ (meters)
///
/// # Returns
///
/// The photon energy E (Joules).
#[no_mangle]

pub extern "C" fn rssn_num_physics_photon_energy(
    wavelength: f64
) -> f64 {

    physics::photon_energy(wavelength)
}

/// Computes the wavelength of a photon from its energy via handle-based FFI.
///
/// The wavelength is λ = hc/E.
///
/// # Arguments
///
/// * `energy` - Photon energy E (Joules)
///
/// # Returns
///
/// The photon wavelength λ (meters).
#[no_mangle]

pub extern "C" fn rssn_num_physics_photon_wavelength(
    energy: f64
) -> f64 {

    physics::photon_wavelength(energy)
}

/// Computes the Compton wavelength of a particle via handle-based FFI.
///
/// The Compton wavelength is λ_C = h/(mc), characteristic of quantum scattering.
///
/// # Arguments
///
/// * `mass` - Particle rest mass m (kg)
///
/// # Returns
///
/// The Compton wavelength λ_C (meters).
#[no_mangle]

pub extern "C" fn rssn_num_physics_compton_wavelength(
    mass: f64
) -> f64 {

    physics::compton_wavelength(mass)
}

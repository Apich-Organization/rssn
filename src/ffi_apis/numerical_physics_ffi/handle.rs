//! Handle-based FFI API for numerical physics functions.

use crate::numerical::physics;

// ============================================================================
// Physical Constants
// ============================================================================

#[no_mangle]

pub extern "C" fn rssn_num_physics_speed_of_light() -> f64 {

    physics::SPEED_OF_LIGHT
}

#[no_mangle]

pub extern "C" fn rssn_num_physics_planck_constant() -> f64 {

    physics::PLANCK_CONSTANT
}

#[no_mangle]

pub extern "C" fn rssn_num_physics_gravitational_constant() -> f64 {

    physics::GRAVITATIONAL_CONSTANT
}

#[no_mangle]

pub extern "C" fn rssn_num_physics_boltzmann_constant() -> f64 {

    physics::BOLTZMANN_CONSTANT
}

#[no_mangle]

pub extern "C" fn rssn_num_physics_elementary_charge() -> f64 {

    physics::ELEMENTARY_CHARGE
}

#[no_mangle]

pub extern "C" fn rssn_num_physics_electron_mass() -> f64 {

    physics::ELECTRON_MASS
}

// ============================================================================
// Classical Mechanics
// ============================================================================

#[no_mangle]

pub extern "C" fn rssn_num_physics_simple_harmonic_oscillator(
    amplitude : f64,
    omega : f64,
    phase : f64,
    time : f64,
) -> f64 {

    physics::simple_harmonic_oscillator(
        amplitude,
        omega,
        phase,
        time,
    )
}

#[no_mangle]

pub extern "C" fn rssn_num_physics_damped_harmonic_oscillator(
    amplitude : f64,
    omega0 : f64,
    gamma : f64,
    phase : f64,
    time : f64,
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

#[no_mangle]

pub extern "C" fn rssn_num_physics_coulomb_force(
    q1 : f64,
    q2 : f64,
    r : f64,
) -> f64 {

    physics::coulomb_force(q1, q2, r)
}

#[no_mangle]

pub extern "C" fn rssn_num_physics_electric_field_point_charge(
    q : f64,
    r : f64,
) -> f64 {

    physics::electric_field_point_charge(q, r)
}

#[no_mangle]

pub extern "C" fn rssn_num_physics_electric_potential_point_charge(
    q : f64,
    r : f64,
) -> f64 {

    physics::electric_potential_point_charge(q, r)
}

#[no_mangle]

pub extern "C" fn rssn_num_physics_magnetic_field_infinite_wire(
    current : f64,
    r : f64,
) -> f64 {

    physics::magnetic_field_infinite_wire(current, r)
}

#[no_mangle]

pub extern "C" fn rssn_num_physics_lorentz_force(
    charge : f64,
    velocity : f64,
    e_field : f64,
    b_field : f64,
) -> f64 {

    physics::lorentz_force(
        charge,
        velocity,
        e_field,
        b_field,
    )
}

#[no_mangle]

pub extern "C" fn rssn_num_physics_cyclotron_radius(
    mass : f64,
    velocity : f64,
    charge : f64,
    b_field : f64,
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

#[no_mangle]

pub extern "C" fn rssn_num_physics_ideal_gas_pressure(
    n : f64,
    t : f64,
    v : f64,
) -> f64 {

    physics::ideal_gas_pressure(n, t, v)
}

#[no_mangle]

pub extern "C" fn rssn_num_physics_ideal_gas_volume(
    n : f64,
    t : f64,
    p : f64,
) -> f64 {

    physics::ideal_gas_volume(n, t, p)
}

#[no_mangle]

pub extern "C" fn rssn_num_physics_ideal_gas_temperature(
    p : f64,
    v : f64,
    n : f64,
) -> f64 {

    physics::ideal_gas_temperature(p, v, n)
}

#[no_mangle]

pub extern "C" fn rssn_num_physics_maxwell_boltzmann_speed_distribution(
    v : f64,
    mass : f64,
    temperature : f64,
) -> f64 {

    physics::maxwell_boltzmann_speed_distribution(v, mass, temperature)
}

#[no_mangle]

pub extern "C" fn rssn_num_physics_maxwell_boltzmann_mean_speed(
    mass : f64,
    temperature : f64,
) -> f64 {

    physics::maxwell_boltzmann_mean_speed(mass, temperature)
}

#[no_mangle]

pub extern "C" fn rssn_num_physics_maxwell_boltzmann_rms_speed(
    mass : f64,
    temperature : f64,
) -> f64 {

    physics::maxwell_boltzmann_rms_speed(mass, temperature)
}

#[no_mangle]

pub extern "C" fn rssn_num_physics_blackbody_power(
    area : f64,
    temperature : f64,
) -> f64 {

    physics::blackbody_power(area, temperature)
}

#[no_mangle]

pub extern "C" fn rssn_num_physics_wien_displacement_wavelength(temperature : f64) -> f64 {

    physics::wien_displacement_wavelength(temperature)
}

// ============================================================================
// Special Relativity
// ============================================================================

#[no_mangle]

pub extern "C" fn rssn_num_physics_lorentz_factor(velocity : f64) -> f64 {

    physics::lorentz_factor(velocity)
}

#[no_mangle]

pub extern "C" fn rssn_num_physics_time_dilation(
    proper_time : f64,
    velocity : f64,
) -> f64 {

    physics::time_dilation(
        proper_time,
        velocity,
    )
}

#[no_mangle]

pub extern "C" fn rssn_num_physics_length_contraction(
    proper_length : f64,
    velocity : f64,
) -> f64 {

    physics::length_contraction(
        proper_length,
        velocity,
    )
}

#[no_mangle]

pub extern "C" fn rssn_num_physics_relativistic_momentum(
    mass : f64,
    velocity : f64,
) -> f64 {

    physics::relativistic_momentum(mass, velocity)
}

#[no_mangle]

pub extern "C" fn rssn_num_physics_relativistic_kinetic_energy(
    mass : f64,
    velocity : f64,
) -> f64 {

    physics::relativistic_kinetic_energy(mass, velocity)
}

#[no_mangle]

pub extern "C" fn rssn_num_physics_mass_energy(mass : f64) -> f64 {

    physics::mass_energy(mass)
}

#[no_mangle]

pub extern "C" fn rssn_num_physics_relativistic_velocity_addition(
    v : f64,
    w : f64,
) -> f64 {

    physics::relativistic_velocity_addition(v, w)
}

// ============================================================================
// Quantum Mechanics
// ============================================================================

#[no_mangle]

pub extern "C" fn rssn_num_physics_quantum_harmonic_oscillator_energy(
    n : u64,
    omega : f64,
) -> f64 {

    physics::quantum_harmonic_oscillator_energy(n, omega)
}

#[no_mangle]

pub extern "C" fn rssn_num_physics_hydrogen_energy_level(n : u64) -> f64 {

    physics::hydrogen_energy_level(n)
}

#[no_mangle]

pub extern "C" fn rssn_num_physics_de_broglie_wavelength(momentum : f64) -> f64 {

    physics::de_broglie_wavelength(momentum)
}

#[no_mangle]

pub extern "C" fn rssn_num_physics_photon_energy(wavelength : f64) -> f64 {

    physics::photon_energy(wavelength)
}

#[no_mangle]

pub extern "C" fn rssn_num_physics_photon_wavelength(energy : f64) -> f64 {

    physics::photon_wavelength(energy)
}

#[no_mangle]

pub extern "C" fn rssn_num_physics_compton_wavelength(mass : f64) -> f64 {

    physics::compton_wavelength(mass)
}

//! Comprehensive tests for numerical physics module.
//!
//! Tests for physical constants, classical mechanics, electromagnetism,
//! thermodynamics, special relativity, and quantum mechanics.

use rssn::numerical::physics::*;
use std::f64::consts::PI;

// ============================================================================
// Physical Constants Tests
// ============================================================================

#[test]

fn test_speed_of_light() {

    assert!((SPEED_OF_LIGHT - 299_792_458.0).abs() < 1.0);
}

#[test]

fn test_planck_constant() {

    assert!(PLANCK_CONSTANT > 6.62e-34 && PLANCK_CONSTANT < 6.63e-34);
}

#[test]

fn test_gravitational_constant() {

    assert!(GRAVITATIONAL_CONSTANT > 6.67e-11 && GRAVITATIONAL_CONSTANT < 6.68e-11);
}

#[test]

fn test_boltzmann_constant() {

    assert!(BOLTZMANN_CONSTANT > 1.38e-23 && BOLTZMANN_CONSTANT < 1.39e-23);
}

#[test]

fn test_elementary_charge() {

    assert!(ELEMENTARY_CHARGE > 1.60e-19 && ELEMENTARY_CHARGE < 1.61e-19);
}

// ============================================================================
// Particle3D Tests
// ============================================================================

#[test]

fn test_particle3d_new() {

    let p = Particle3D::new(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0);

    assert_eq!(p.mass, 1.0);

    assert_eq!(p.vx, 1.0);
}

#[test]

fn test_particle3d_kinetic_energy() {

    let p = Particle3D::new(2.0, 0.0, 0.0, 0.0, 3.0, 4.0, 0.0);

    // KE = 0.5 * m * v^2 = 0.5 * 2 * (9+16) = 25
    assert!((p.kinetic_energy() - 25.0).abs() < 1e-10);
}

#[test]

fn test_particle3d_momentum() {

    let p = Particle3D::new(2.0, 0.0, 0.0, 0.0, 3.0, 4.0, 0.0);

    // p = m * v = 2 * 5 = 10
    assert!((p.momentum() - 10.0).abs() < 1e-10);
}

// ============================================================================
// Classical Mechanics Tests
// ============================================================================

#[test]

fn test_simple_harmonic_oscillator() {

    // At t=0, x = A * cos(φ)
    let x = simple_harmonic_oscillator(2.0, 1.0, 0.0, 0.0);

    assert!((x - 2.0).abs() < 1e-10);

    // At t=π/(2ω), x = A * cos(π/2) = 0
    let x = simple_harmonic_oscillator(2.0, 1.0, 0.0, PI / 2.0);

    assert!(x.abs() < 1e-10);
}

#[test]

fn test_damped_harmonic_oscillator() {

    // At t=0, should be at amplitude
    let x = damped_harmonic_oscillator(2.0, 10.0, 0.5, 0.0, 0.0);

    assert!((x - 2.0).abs() < 1e-10);

    // After some time, amplitude should decay
    let x1 = damped_harmonic_oscillator(2.0, 10.0, 0.5, 0.0, 0.0);

    let x2 = damped_harmonic_oscillator(2.0, 10.0, 0.5, 0.0, 1.0);

    assert!(x2.abs() < x1.abs());
}

#[test]

fn test_projectile_motion_with_drag() {

    let trajectory = projectile_motion_with_drag(
        10.0,     // v0
        PI / 4.0, // 45 degrees
        1.0,      // mass
        0.47,     // drag coefficient
        0.01,     // area
        1.225,    // air density
        0.001,    // dt
        10.0,     // max time
    );

    assert!(!trajectory.is_empty());

    // First point should be at origin
    let (_, x, y, _, _) = trajectory[0];

    assert!(x.abs() < 1e-10);

    assert!(y.abs() < 1e-10);
}

// ============================================================================
// N-Body Tests
// ============================================================================

#[test]

fn test_simulate_n_body() {

    let particles = vec![
        Particle3D::new(1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        Particle3D::new(1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0),
    ];

    let snapshots = simulate_n_body(particles, 0.01, 10, GRAVITATIONAL_CONSTANT);

    assert_eq!(snapshots.len(), 11); // Initial + 10 steps
}

#[test]

fn test_gravitational_potential_energy() {

    let particles = vec![
        Particle3D::new(1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        Particle3D::new(1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0),
    ];

    let pe = gravitational_potential_energy(&particles, GRAVITATIONAL_CONSTANT);

    // Should be negative (bound system)
    assert!(pe < 0.0);
}

#[test]

fn test_total_kinetic_energy() {

    let particles = vec![
        Particle3D::new(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0),
        Particle3D::new(1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0),
    ];

    let ke = total_kinetic_energy(&particles);

    // Each has KE = 0.5 * 1 * 1 = 0.5
    assert!((ke - 1.0).abs() < 1e-10);
}

// ============================================================================
// Electromagnetism Tests
// ============================================================================

#[test]

fn test_coulomb_force() {

    // Two electrons at 1m apart
    let f = coulomb_force(ELEMENTARY_CHARGE, ELEMENTARY_CHARGE, 1.0);

    assert!(f > 0.0); // Repulsive
    assert!(f > 2.3e-28 && f < 2.4e-28);
}

#[test]

fn test_coulomb_force_sign() {

    // Opposite charges should attract (negative force)
    let f = coulomb_force(ELEMENTARY_CHARGE, -ELEMENTARY_CHARGE, 1.0);

    assert!(f < 0.0);
}

#[test]

fn test_electric_field_point_charge() {

    let e = electric_field_point_charge(ELEMENTARY_CHARGE, 1.0);

    assert!(e > 0.0);
}

#[test]

fn test_electric_potential_point_charge() {

    let v = electric_potential_point_charge(ELEMENTARY_CHARGE, 1.0);

    assert!(v > 0.0);
}

#[test]

fn test_magnetic_field_infinite_wire() {

    let b = magnetic_field_infinite_wire(1.0, 0.1);

    // B = μ₀I/(2πr) = 4π×10⁻⁷ * 1 / (2π * 0.1) = 2×10⁻⁶ T
    assert!(b > 1e-6 && b < 3e-6);
}

#[test]

fn test_lorentz_force() {

    let f = lorentz_force(ELEMENTARY_CHARGE, 1e6, 0.0, 1.0);

    // F = qvB = 1.6e-19 * 1e6 * 1 = 1.6e-13 N
    assert!(f > 1e-13 && f < 2e-13);
}

#[test]

fn test_cyclotron_radius() {

    let r = cyclotron_radius(ELECTRON_MASS, 1e6, ELEMENTARY_CHARGE, 1.0);

    // r = mv/(qB) for electron at 1e6 m/s in 1T field
    assert!(r > 5e-6 && r < 6e-6);
}

// ============================================================================
// Thermodynamics Tests
// ============================================================================

#[test]

fn test_ideal_gas_pressure() {

    // 1 mol at 300K in 1m³
    let p = ideal_gas_pressure(1.0, 300.0, 1.0);

    // P = nRT/V = 1 * 8.314 * 300 / 1 ≈ 2494 Pa
    assert!(p > 2400.0 && p < 2600.0);
}

#[test]

fn test_ideal_gas_volume() {

    let v = ideal_gas_volume(1.0, 300.0, 101325.0);

    // V = nRT/P ≈ 0.0246 m³
    assert!(v > 0.024 && v < 0.025);
}

#[test]

fn test_ideal_gas_temperature() {

    let t = ideal_gas_temperature(101325.0, 0.0224, 1.0);

    // T = PV/(nR) ≈ 273K
    assert!(t > 270.0 && t < 280.0);
}

#[test]

fn test_maxwell_boltzmann_mean_speed() {

    // Mean speed of N2 at 300K
    let v = maxwell_boltzmann_mean_speed(NEUTRON_MASS * 28.0, 300.0);

    // Should be around 475 m/s
    assert!(v > 400.0 && v < 550.0);
}

#[test]

fn test_maxwell_boltzmann_rms_speed() {

    let v_rms = maxwell_boltzmann_rms_speed(NEUTRON_MASS * 28.0, 300.0);

    let v_mean = maxwell_boltzmann_mean_speed(NEUTRON_MASS * 28.0, 300.0);

    // RMS speed should be slightly higher than mean speed
    assert!(v_rms > v_mean);
}

#[test]

fn test_blackbody_power() {

    // Sun surface (5778K, r≈7e8m)
    let power_per_m2 = blackbody_power(1.0, 5778.0);

    // σT⁴ ≈ 6.3×10⁷ W/m²
    assert!(power_per_m2 > 6e7 && power_per_m2 < 7e7);
}

#[test]

fn test_wien_displacement_wavelength() {

    // Sun surface temperature
    let wavelength = wien_displacement_wavelength(5778.0);

    // Should be around 500nm (visible light)
    assert!(wavelength > 4e-7 && wavelength < 6e-7);
}

// ============================================================================
// Special Relativity Tests
// ============================================================================

#[test]

fn test_lorentz_factor_low_speed() {

    // At low speeds, γ ≈ 1
    let gamma = lorentz_factor(1000.0);

    assert!((gamma - 1.0).abs() < 1e-10);
}

#[test]

fn test_lorentz_factor_high_speed() {

    // At 0.9c, γ ≈ 2.29
    let gamma = lorentz_factor(0.9 * SPEED_OF_LIGHT);

    assert!(gamma > 2.2 && gamma < 2.4);
}

#[test]

fn test_lorentz_factor_very_high_speed() {

    // At 0.99c, γ ≈ 7.09
    let gamma = lorentz_factor(0.99 * SPEED_OF_LIGHT);

    assert!(gamma > 7.0 && gamma < 7.2);
}

#[test]

fn test_time_dilation() {

    // 1 second proper time at 0.9c
    let dilated = time_dilation(1.0, 0.9 * SPEED_OF_LIGHT);

    assert!(dilated > 2.2 && dilated < 2.4);
}

#[test]

fn test_length_contraction() {

    // 1 meter at 0.9c
    let contracted = length_contraction(1.0, 0.9 * SPEED_OF_LIGHT);

    assert!(contracted > 0.4 && contracted < 0.5);
}

#[test]

fn test_relativistic_momentum() {

    // At low speed, p ≈ mv
    let p = relativistic_momentum(1.0, 1.0);

    assert!((p - 1.0).abs() < 1e-10);
}

#[test]

fn test_mass_energy() {

    // 1 kg -> E = mc² ≈ 9×10¹⁶ J
    let e = mass_energy(1.0);

    assert!(e > 8e16 && e < 1e17);
}

#[test]

fn test_relativistic_velocity_addition() {

    // Two velocities of 0.5c should give < c
    let u = relativistic_velocity_addition(0.5 * SPEED_OF_LIGHT, 0.5 * SPEED_OF_LIGHT);

    assert!(u < SPEED_OF_LIGHT);

    // Should be 0.8c
    assert!((u / SPEED_OF_LIGHT - 0.8).abs() < 0.01);
}

// ============================================================================
// Quantum Mechanics Tests
// ============================================================================

#[test]

fn test_quantum_harmonic_oscillator_energy() {

    // Ground state (n=0): E = ħω/2
    let e0 = quantum_harmonic_oscillator_energy(0, 1.0);

    assert!((e0 - HBAR * 0.5).abs() < 1e-40);

    // First excited state (n=1): E = 3ħω/2
    let e1 = quantum_harmonic_oscillator_energy(1, 1.0);

    assert!((e1 - HBAR * 1.5).abs() < 1e-40);
}

#[test]

fn test_hydrogen_energy_level() {

    // Ground state (n=1): E = -13.6 eV
    let e1 = hydrogen_energy_level(1);

    assert!(e1 < 0.0);

    // E1 in Joules ≈ -2.18×10⁻¹⁸ J
    assert!(e1 > -2.2e-18 && e1 < -2.1e-18);

    // n=2: E = E1/4
    let e2 = hydrogen_energy_level(2);

    assert!((e2 - e1 / 4.0).abs() < 1e-20);
}

#[test]

fn test_de_broglie_wavelength() {

    // Electron at 1 eV kinetic energy
    // p = sqrt(2mE) ≈ 5.4×10⁻²⁵ kg⋅m/s
    // λ = h/p ≈ 1.2 nm
    let p = (2.0 * ELECTRON_MASS * 1.0 * ELEMENTARY_CHARGE).sqrt();

    let lambda = de_broglie_wavelength(p);

    assert!(lambda > 1e-9 && lambda < 2e-9);
}

#[test]

fn test_photon_energy() {

    // 500nm light (green)
    let e = photon_energy(500e-9);

    // E ≈ 4×10⁻¹⁹ J ≈ 2.5 eV
    assert!(e > 3e-19 && e < 5e-19);
}

#[test]

fn test_photon_wavelength() {

    // Energy of 2.5 eV
    let lambda = photon_wavelength(2.5 * ELEMENTARY_CHARGE);

    // Should be around 500nm
    assert!(lambda > 4e-7 && lambda < 6e-7);
}

#[test]

fn test_compton_wavelength() {

    // Electron Compton wavelength ≈ 2.43×10⁻¹² m
    let lambda_c = compton_wavelength(ELECTRON_MASS);

    assert!(lambda_c > 2.4e-12 && lambda_c < 2.5e-12);
}

// ============================================================================
// Property Tests
// ============================================================================

mod proptests {

    use super::*;
    use proptest::prelude::*;

    proptest! {
        #[test]
        fn prop_lorentz_factor_positive(v in 0.0..(0.99 * SPEED_OF_LIGHT)) {
            let gamma = lorentz_factor(v);
            prop_assert!(gamma >= 1.0);
        }

        #[test]
        fn prop_harmonic_oscillator_bounded(
            amplitude in 0.1..10.0f64,
            omega in 0.1..10.0f64,
            phase in 0.0..2.0*PI,
            time in 0.0..10.0f64
        ) {
            let x = simple_harmonic_oscillator(amplitude, omega, phase, time);
            prop_assert!(x.abs() <= amplitude + 1e-10);
        }

        #[test]
        fn prop_ideal_gas_consistency(
            n in 0.1..10.0f64,
            t in 100.0..500.0f64,
            v in 0.01..1.0f64
        ) {
            let p = ideal_gas_pressure(n, t, v);
            let v2 = ideal_gas_volume(n, t, p);
            prop_assert!((v - v2).abs() < 1e-6);
        }

        #[test]
        fn prop_photon_energy_wavelength_inverse(wavelength in 1e-9..1e-6f64) {
            let e = photon_energy(wavelength);
            let lambda2 = photon_wavelength(e);
            prop_assert!((wavelength - lambda2).abs() < 1e-15);
        }

        #[test]
        fn prop_coulomb_inverse_square(r in 0.01..10.0f64) {
            let f1 = coulomb_force(1.0, 1.0, r);
            let f2 = coulomb_force(1.0, 1.0, 2.0 * r);
            // Force at 2r should be 1/4 of force at r
            prop_assert!((f1 / 4.0 - f2).abs() < 1e-10 * f1.abs());
        }
    }
}

//! Classical physics helpers and small PDE/eigenvalue solvers (ported from
//! `numerical_physics_test.rs`; the old `rssn::numerical::physics` is now `rssn::sim::classical`).
//!
//! Tests for physical constants, classical mechanics, electromagnetism,
//! thermodynamics, special relativity, and quantum mechanics.

use std::f64::consts::PI;

use rssn::sim::classical::*;

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
    let params = ProjectileParams {
        v0: 10.0,
        angle: PI / 4.0,
        mass: 1.0,
        drag_coeff: 0.47,
        area: 0.01,
        air_density: 1.225,
        dt: 0.001,
        max_time: 10.0,
    };

    let trajectory = projectile_motion_with_drag(params);

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

    use proptest::prelude::*;

    use super::*;

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

// ============================================================================
// Added: closed-form checks (constants, relativity, thermodynamics, quantum)
// ============================================================================

mod strengthened {
    use std::f64::consts::PI;

    use proptest::prelude::*;
    use proptest::test_runner::RngSeed;
    use rssn::kernels::matrix::Matrix;
    use rssn::sim::classical::*;

    fn cfg() -> ProptestConfig {
        ProptestConfig {
            cases: 48,
            rng_seed: RngSeed::Fixed(0x5EED),
            failure_persistence: None,
            ..ProptestConfig::default()
        }
    }

    fn rel(
        a: f64,
        b: f64,
    ) -> f64 {
        (a - b).abs() / b.abs().max(f64::MIN_POSITIVE)
    }

    #[test]
    fn constants_are_mutually_consistent() {
        // c^2 = 1 / (eps0 mu0)
        assert!(
            rel(
                1.0 / (VACUUM_PERMITTIVITY * VACUUM_PERMEABILITY),
                SPEED_OF_LIGHT.powi(2)
            ) < 1e-8
        );
        // k_e = 1 / (4 pi eps0)
        assert!(rel(COULOMB_CONSTANT, 1.0 / (4.0 * PI * VACUUM_PERMITTIVITY)) < 1e-8);
        // hbar = h / 2 pi
        assert!(rel(HBAR, PLANCK_CONSTANT / (2.0 * PI)) < 1e-9);
        // R = N_A k_B
        assert!(rel(GAS_CONSTANT, AVOGADRO_NUMBER * BOLTZMANN_CONSTANT) < 1e-9);
        // Bohr radius = hbar / (m_e c alpha)
        assert!(
            rel(
                BOHR_RADIUS,
                HBAR / (ELECTRON_MASS * SPEED_OF_LIGHT * FINE_STRUCTURE_CONSTANT)
            ) < 1e-6
        );
        // Stefan-Boltzmann sigma = 2 pi^5 k^4 / (15 h^3 c^2)
        let sigma = 2.0 * PI.powi(5) * BOLTZMANN_CONSTANT.powi(4)
            / (15.0 * PLANCK_CONSTANT.powi(3) * SPEED_OF_LIGHT.powi(2));
        assert!(rel(STEFAN_BOLTZMANN, sigma) < 1e-8);
        assert!((STANDARD_GRAVITY - 9.80665).abs() < 1e-12);
    }

    #[test]
    fn coulomb_electric_field_and_potential_agree() {
        let q = ELEMENTARY_CHARGE;
        let r = 0.37;
        assert!(rel(coulomb_force(q, q, r), COULOMB_CONSTANT * q * q / (r * r)) < 1e-9);
        // F = q E for a test charge q and V = E * r for a point charge.
        assert!(
            rel(
                coulomb_force(q, q, r),
                q * electric_field_point_charge(q, r)
            ) < 1e-9
        );
        assert!(
            rel(
                electric_potential_point_charge(q, r),
                electric_field_point_charge(q, r) * r
            ) < 1e-9
        );
        assert_eq!(coulomb_force(1.0, 1.0, 0.0), f64::INFINITY);
    }

    #[test]
    fn magnetic_helpers_reference_values() {
        // B = mu0 I / (2 pi r)
        assert!(
            rel(
                magnetic_field_infinite_wire(1.0, 0.1),
                VACUUM_PERMEABILITY / (2.0 * PI * 0.1)
            ) < 1e-9
        );
        // Lorentz force is q v B at 90 degrees and vanishes for parallel motion.
        assert!(
            rel(
                lorentz_force(ELEMENTARY_CHARGE, 1e6, 0.0, 1.0),
                ELEMENTARY_CHARGE * 1e6
            ) < 1e-9
        );
        assert!(lorentz_force(ELEMENTARY_CHARGE, 1e6, PI / 2.0, 1.0).is_finite());
        // r = m v / (q B)
        assert!(
            rel(
                cyclotron_radius(ELECTRON_MASS, 1e6, ELEMENTARY_CHARGE, 1.0),
                ELECTRON_MASS * 1e6 / ELEMENTARY_CHARGE
            ) < 1e-9
        );
    }

    #[test]
    fn ideal_gas_law_molar_volume() {
        // 1 mol at 273.15 K and 101325 Pa occupies 22.414 L.
        let v = ideal_gas_volume(1.0, 273.15, 101_325.0);
        assert!((v - 0.022_414).abs() < 2e-5, "V = {v}");
        assert!(
            rel(
                ideal_gas_pressure(2.0, 300.0, 0.5),
                2.0 * GAS_CONSTANT * 300.0 / 0.5
            ) < 1e-12
        );
        assert!(
            rel(
                ideal_gas_temperature(ideal_gas_pressure(1.5, 350.0, 0.2), 0.2, 1.5),
                350.0
            ) < 1e-12
        );
    }

    #[test]
    fn kinetic_theory_speed_ratios() {
        let (m, t) = (28.0 * ATOMIC_MASS_UNIT, 300.0);
        let (mean, rms) = (
            maxwell_boltzmann_mean_speed(m, t),
            maxwell_boltzmann_rms_speed(m, t),
        );
        // v_rms / v_mean = sqrt(3 pi / 8)
        assert!(rel(rms / mean, (3.0 * PI / 8.0).sqrt()) < 1e-12);
        // N2 at 300 K: mean speed about 476 m/s
        assert!((mean - 476.0).abs() < 2.0, "mean speed {mean}");
        // 1/2 m v_rms^2 = 3/2 k T
        assert!(rel(0.5 * m * rms * rms, 1.5 * BOLTZMANN_CONSTANT * t) < 1e-12);
    }

    #[test]
    fn maxwell_boltzmann_distribution_is_normalised_and_peaks_at_most_probable_speed() {
        let (m, t) = (28.0 * ATOMIC_MASS_UNIT, 300.0);
        let v_mp = (2.0 * BOLTZMANN_CONSTANT * t / m).sqrt();
        // Trapezoid integral over [0, 8 v_mp]
        let n = 4000;
        let h = 8.0 * v_mp / n as f64;
        let integral: f64 = (0..n)
            .map(|i| {
                let (a, b) = (i as f64 * h, (i + 1) as f64 * h);
                0.5 * h
                    * (maxwell_boltzmann_speed_distribution(a, m, t)
                        + maxwell_boltzmann_speed_distribution(b, m, t))
            })
            .sum();
        assert!((integral - 1.0).abs() < 1e-6, "integral {integral}");
        let peak = maxwell_boltzmann_speed_distribution(v_mp, m, t);
        assert!(peak > maxwell_boltzmann_speed_distribution(0.9 * v_mp, m, t));
        assert!(peak > maxwell_boltzmann_speed_distribution(1.1 * v_mp, m, t));
        assert_eq!(maxwell_boltzmann_speed_distribution(100.0, m, 0.0), 0.0);
    }

    #[test]
    fn blackbody_and_wien() {
        // Sun: sigma T^4 = 6.32e7 W/m^2 ; Wien peak 501.5 nm
        let flux = blackbody_power(1.0, 5778.0);
        assert!((flux - 6.32e7).abs() < 0.02e7, "flux {flux}");
        assert!((wien_displacement_wavelength(5778.0) - 501.5e-9).abs() < 1e-9);
        assert_eq!(wien_displacement_wavelength(0.0), f64::INFINITY);
        assert!(
            rel(
                blackbody_power(4.0, 300.0),
                4.0 * blackbody_power(1.0, 300.0)
            ) < 1e-12
        );
    }

    #[test]
    fn relativity_reference_values() {
        let c = SPEED_OF_LIGHT;
        assert!((lorentz_factor(0.6 * c) - 1.25).abs() < 1e-12);
        assert!((lorentz_factor(0.8 * c) - 5.0 / 3.0).abs() < 1e-12);
        assert_eq!(lorentz_factor(c), f64::INFINITY);
        assert_eq!(lorentz_factor(2.0 * c), f64::INFINITY);
        assert!((time_dilation(1.0, 0.6 * c) - 1.25).abs() < 1e-12);
        assert!((length_contraction(1.0, 0.6 * c) - 0.8).abs() < 1e-12);
        // E^2 = (pc)^2 + (m c^2)^2
        let (m, v) = (2.5, 0.7 * c);
        let e = relativistic_total_energy(m, v);
        let p = relativistic_momentum(m, v);
        assert!(rel(e * e, (p * c).powi(2) + mass_energy(m).powi(2)) < 1e-12);
        assert!(rel(relativistic_kinetic_energy(m, v), e - mass_energy(m)) < 1e-12);
        // Kinetic energy tends to (1/2) m v^2 at low speed.
        // (gamma - 1 loses about 5 digits to cancellation at 1 km/s, hence the loose bound.)
        assert!(rel(relativistic_kinetic_energy(1.0, 1000.0), 0.5 * 1e6) < 1e-3);
        // Velocity addition: 0.5c (+) 0.5c = 0.8c exactly, and c (+) anything = c.
        assert!((relativistic_velocity_addition(0.5 * c, 0.5 * c) / c - 0.8).abs() < 1e-12);
        assert!((relativistic_velocity_addition(c, 0.9 * c) - c).abs() < 1e-6);
    }

    #[test]
    fn quantum_helpers_reference_values() {
        // Hydrogen ground state 13.6057 eV; Lyman-alpha = 10.204 eV.
        let e1 = hydrogen_energy_level(1) / ELEMENTARY_CHARGE;
        assert!((e1 + 13.605_69).abs() < 1e-3, "{e1} eV");
        let lyman_alpha = (hydrogen_energy_level(2) - hydrogen_energy_level(1)) / ELEMENTARY_CHARGE;
        assert!((lyman_alpha - 10.204).abs() < 1e-3);
        assert_eq!(hydrogen_energy_level(0), f64::NEG_INFINITY);
        // Photon: 1 eV <-> 1239.84 nm
        let lambda = photon_wavelength(ELEMENTARY_CHARGE);
        assert!((lambda - 1239.84e-9).abs() < 1e-11);
        // Compton wavelength of the electron: 2.42631e-12 m
        assert!((compton_wavelength(ELECTRON_MASS) - 2.426_31e-12).abs() < 1e-16);
        // Electron at 1 eV: lambda = 1.2264 nm
        let p = (2.0 * ELECTRON_MASS * ELEMENTARY_CHARGE).sqrt();
        assert!((de_broglie_wavelength(p) - 1.2264e-9).abs() < 1e-13);
        // Energy spacing of the oscillator is hbar omega.
        assert!(
            rel(
                quantum_harmonic_oscillator_energy(5, 2.0)
                    - quantum_harmonic_oscillator_energy(4, 2.0),
                2.0 * HBAR
            ) < 1e-9
        );
        assert!(rel(heisenberg_position_uncertainty(1e-24), HBAR / 2e-24) < 1e-12);
        assert_eq!(heisenberg_position_uncertainty(0.0), 0.0);
        assert_eq!(de_broglie_wavelength(0.0), f64::INFINITY);
        assert_eq!(photon_energy(0.0), f64::INFINITY);
    }

    #[test]
    fn oscillators_analytic() {
        assert!((simple_harmonic_oscillator(3.0, 2.0, PI / 2.0, 0.0)).abs() < 1e-12);
        assert!((simple_harmonic_oscillator(3.0, 2.0, 0.0, PI) - 3.0).abs() < 1e-12);
        // Underdamped: A e^{-g t} cos(sqrt(w0^2 - g^2) t)
        let x = damped_harmonic_oscillator(1.0, 5.0, 3.0, 0.0, 0.2);
        assert!((x - (-0.6f64).exp() * (0.8f64).cos()).abs() < 1e-12);
        // Critically damped: A (1 + g t) e^{-g t}
        let x = damped_harmonic_oscillator(2.0, 3.0, 3.0, 0.0, 0.5);
        assert!((x - 2.0 * (1.0 + 1.5) * (-1.5f64).exp()).abs() < 1e-12);
        // Overdamped decays monotonically towards zero.
        let a = damped_harmonic_oscillator(1.0, 1.0, 3.0, 0.0, 1.0);
        let b = damped_harmonic_oscillator(1.0, 1.0, 3.0, 0.0, 3.0);
        assert!(a.is_finite() && b.is_finite() && b.abs() < a.abs() * 10.0);
    }

    #[test]
    fn particle_helpers() {
        let p = Particle3D::new(2.0, 1.0, 2.0, 3.0, 3.0, 4.0, 0.0);
        assert_eq!(
            (p.x, p.y, p.z, p.vx, p.vy, p.vz),
            (1.0, 2.0, 3.0, 3.0, 4.0, 0.0)
        );
        let ps = [p, Particle3D::new(1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0)];
        assert!(rel(total_kinetic_energy(&ps), 25.0 + 2.0) < 1e-12);
    }

    #[test]
    fn gravitational_potential_energy_reference_value() {
        let ps = [
            Particle3D::new(2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            Particle3D::new(3.0, 3.0, 4.0, 0.0, 0.0, 0.0, 0.0),
        ];
        // U = -G m1 m2 / r = -G * 6 / 5
        assert!(rel(gravitational_potential_energy(&ps, 1.0), -1.2) < 1e-12);
        assert_eq!(gravitational_potential_energy(&ps[..1], 1.0), 0.0);
    }

    #[test]
    fn n_body_circular_orbit_conserves_energy_and_momentum() {
        // Two unit masses (G = 1) on a circular orbit about their centre of mass: v = sqrt(G m / (2 d)).
        let v = (1.0f64 / 2.0).sqrt();
        let start = vec![
            Particle3D::new(1.0, -0.5, 0.0, 0.0, 0.0, -v, 0.0),
            Particle3D::new(1.0, 0.5, 0.0, 0.0, 0.0, v, 0.0),
        ];
        let e =
            |ps: &[Particle3D]| total_kinetic_energy(ps) + gravitational_potential_energy(ps, 1.0);
        let e0 = e(&start);
        // KE = 2 * (1/2) * (1/2) = 1/2 and PE = -1, so E = -1/2.
        assert!((e0 + 0.5).abs() < 1e-12, "E0 = {e0}");
        let snaps = simulate_n_body(start, 1e-3, 3000, 1.0);
        assert_eq!(snaps.len(), 3001);
        let (mut max_dev, mut max_p) = (0.0f64, 0.0f64);
        for s in &snaps {
            max_dev = max_dev.max(((e(s) - e0) / e0).abs());
            let px: f64 = s.iter().map(|p| p.mass * p.vx).sum();
            let py: f64 = s.iter().map(|p| p.mass * p.vy).sum();
            max_p = max_p.max(px.abs().max(py.abs()));
            // separation stays 1 on a circular orbit
            let d = ((s[1].x - s[0].x).powi(2) + (s[1].y - s[0].y).powi(2)).sqrt();
            assert!((d - 1.0).abs() < 5e-3, "separation drifted to {d}");
        }
        assert!(max_dev < 5e-3, "relative energy drift {max_dev}");
        assert!(max_p < 1e-10, "momentum drift {max_p}");
    }

    #[test]
    fn projectile_without_drag_matches_the_range_formula() {
        let params = ProjectileParams {
            v0: 20.0,
            angle: PI / 4.0,
            mass: 1.0,
            drag_coeff: 0.0,
            area: 0.01,
            air_density: 1.225,
            dt: 1e-4,
            max_time: 10.0,
        };
        let traj = projectile_motion_with_drag(params);
        let range = 20.0f64.powi(2) * (2.0 * PI / 4.0).sin() / STANDARD_GRAVITY;
        let (t_end, x_end, y_end, ..) =
            *traj.last().unwrap_or(&(0.0, f64::NAN, f64::NAN, 0.0, 0.0));
        assert!((x_end - range).abs() < 0.05, "range {x_end} vs {range}");
        assert!(y_end >= 0.0 && y_end < 0.01);
        // Flight time 2 v0 sin(theta) / g
        assert!((t_end - 2.0 * 20.0 * (PI / 4.0).sin() / STANDARD_GRAVITY).abs() < 0.01);
        // Peak height v0^2 sin^2 / (2 g)
        let hmax = traj.iter().map(|p| p.2).fold(0.0, f64::max);
        assert!((hmax - 400.0 * 0.5 / (2.0 * STANDARD_GRAVITY)).abs() < 0.01);
    }

    #[test]
    fn drag_shortens_the_range_and_lowers_the_peak() {
        let base = ProjectileParams {
            v0: 30.0,
            angle: PI / 4.0,
            mass: 0.5,
            drag_coeff: 0.0,
            area: 0.02,
            air_density: 1.225,
            dt: 1e-3,
            max_time: 20.0,
        };
        let no_drag = projectile_motion_with_drag(base);
        let drag = projectile_motion_with_drag(ProjectileParams {
            drag_coeff: 0.5,
            ..base
        });
        let far = |t: &Vec<(f64, f64, f64, f64, f64)>| t.last().map_or(f64::NAN, |p| p.1);
        assert!(far(&drag) < 0.95 * far(&no_drag));
        // With drag, the descent is steeper than the ascent: horizontal velocity only decreases.
        assert!(drag.windows(2).all(|w| w[1].3 <= w[0].3 + 1e-12));
    }

    #[test]
    fn particle_motion_in_uniform_gravity_is_exact() {
        // F = (0, 0, -m g); RK4 is exact for quadratic trajectories.
        let (m, g) = (2.0, 9.81);
        let traj = simulate_particle_motion(
            |_, _| [0.0, 0.0, -m * g],
            m,
            (1.0, 2.0, 10.0),
            (3.0, 0.0, 5.0),
            (0.0, 2.0),
            20,
        )
        .unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(traj.len(), 21);
        let last = &traj[20];
        assert_eq!(last.len(), 6);
        assert!((last[0] - (1.0 + 3.0 * 2.0)).abs() < 1e-12);
        assert!((last[1] - 2.0).abs() < 1e-12);
        assert!((last[2] - (10.0 + 5.0 * 2.0 - 0.5 * g * 4.0)).abs() < 1e-10);
        assert!((last[5] - (5.0 - g * 2.0)).abs() < 1e-10);
    }

    #[test]
    fn particle_motion_in_a_harmonic_trap_matches_cosine_and_conserves_energy() {
        let (m, k) = (1.5_f64, 6.0_f64);
        let omega = (k / m).sqrt();
        let traj = simulate_particle_motion(
            |_, s| [-k * s[0], -k * s[1], -k * s[2]],
            m,
            (1.0, 0.0, 0.0),
            (0.0, 0.0, 0.0),
            (0.0, 5.0),
            2000,
        )
        .unwrap_or_else(|e| panic!("{e}"));
        let e0 = 0.5 * k;
        for (i, s) in traj.iter().enumerate() {
            let t = 5.0 * i as f64 / 2000.0;
            assert!((s[0] - (omega * t).cos()).abs() < 1e-8, "step {i}");
            let energy = 0.5 * m * (s[3] * s[3] + s[4] * s[4] + s[5] * s[5])
                + 0.5 * k * (s[0] * s[0] + s[1] * s[1] + s[2] * s[2]);
            assert!((energy - e0).abs() < 1e-8, "energy at step {i}: {energy}");
        }
    }

    #[test]
    fn particle_motion_rejects_zero_steps_and_blowup() {
        assert!(
            simulate_particle_motion(
                |_, _| [0.0; 3],
                1.0,
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 1.0),
                0
            )
            .is_err()
        );
        let r = simulate_particle_motion(
            |_, s| [s[0].powi(3) * 1e3, 0.0, 0.0],
            1.0,
            (10.0, 0.0, 0.0),
            (0.0, 0.0, 0.0),
            (0.0, 50.0),
            50,
        );
        assert!(r.is_err());
    }

    fn sorted(mut v: Vec<f64>) -> Vec<f64> {
        v.sort_by(f64::total_cmp);
        v
    }

    #[test]
    fn schrodinger_1d_particle_in_a_box_matches_the_discrete_spectrum() {
        // The n x n tridiagonal Hamiltonian (V = 0) has exact eigenvalues (1 - cos(k pi / (n+1))) / dx^2.
        let n = 40;
        let (x0, x1) = (0.0, 1.0);
        let dx = (x1 - x0) / (n as f64 - 1.0);
        let (vals, vecs) =
            solve_1d_schrodinger(|_| 0.0, (x0, x1), n).unwrap_or_else(|e| panic!("{e}"));
        let vals = sorted(vals);
        assert_eq!(vals.len(), n);
        for (k, e) in vals.iter().enumerate().take(6) {
            let exact = (1.0 - ((k + 1) as f64 * PI / (n as f64 + 1.0)).cos()) / (dx * dx);
            assert!(
                (e - exact).abs() < 1e-8 * exact.max(1.0),
                "E_{k} = {e} vs {exact}"
            );
        }
        assert!(vecs.is_orthogonal(1e-8), "eigenvectors must be orthonormal");
        // Continuum limit for the lowest level: pi^2 / (2 L^2) with L = (n+1) dx.
        let l = (n as f64 + 1.0) * dx;
        assert!((vals[0] - PI * PI / (2.0 * l * l)).abs() / vals[0] < 5e-3);
    }

    #[test]
    fn schrodinger_1d_harmonic_oscillator_levels() {
        // V = x^2 / 2 (hbar = m = omega = 1): E_n = n + 1/2.
        let (vals, _) = solve_1d_schrodinger(|x| 0.5 * x * x, (-6.0, 6.0), 60)
            .unwrap_or_else(|e| panic!("{e}"));
        let vals = sorted(vals);
        for (n, e) in vals.iter().enumerate().take(3) {
            assert!((e - (n as f64 + 0.5)).abs() < 0.03, "E_{n} = {e}");
        }
    }

    #[test]
    fn schrodinger_1d_ground_state_has_no_nodes() {
        let n = 30;
        let (vals, vecs) =
            solve_1d_schrodinger(|x| 0.5 * x * x, (-5.0, 5.0), n).unwrap_or_else(|e| panic!("{e}"));
        let (i0, _) = vals
            .iter()
            .enumerate()
            .min_by(|a, b| a.1.total_cmp(b.1))
            .unwrap_or((0, &0.0));
        let col: Vec<f64> = (0..n).map(|r| *vecs.get(r, i0)).collect();
        assert!(
            col.iter().all(|&v| v >= 0.0) || col.iter().all(|&v| v <= 0.0),
            "ground state must not change sign"
        );
        assert!((col.iter().map(|v| v * v).sum::<f64>() - 1.0).abs() < 1e-8);
    }

    #[test]
    fn schrodinger_solvers_validate_grid_size() {
        assert!(solve_1d_schrodinger(|_| 0.0, (0.0, 1.0), 2).is_err());
        assert!(solve_2d_schrodinger(|_, _| 0.0, (0.0, 1.0, 0.0, 1.0), (2, 5)).is_err());
        assert!(
            solve_3d_schrodinger(|_, _, _| 0.0, (0.0, 1.0, 0.0, 1.0, 0.0, 1.0), (5, 5, 2)).is_err()
        );
        // Dense 3D solver refuses grids with more than 25000 points.
        assert!(
            solve_3d_schrodinger(|_, _, _| 0.0, (0.0, 1.0, 0.0, 1.0, 0.0, 1.0), (30, 30, 30))
                .is_err()
        );
    }

    #[test]
    fn schrodinger_2d_box_levels_are_sums_of_1d_levels() {
        // Interior 4 x 4 (grid 6 x 6): E = e_a + e_b with e_k = (1 - cos(k pi / 5)) / dx^2.
        let (nx, ny) = (6usize, 6usize);
        let (vals, _) = solve_2d_schrodinger(|_, _| 0.0, (0.0, 1.0, 0.0, 1.0), (nx, ny))
            .unwrap_or_else(|e| panic!("{e}"));
        let vals = sorted(vals);
        assert_eq!(vals.len(), nx * ny);
        let d = 1.0 / (nx as f64 - 1.0);
        let e1 = |k: usize| (1.0 - (k as f64 * PI / (nx as f64 - 1.0)).cos()) / (d * d);
        assert!(
            (vals[0] - 2.0 * e1(1)).abs() < 1e-6,
            "ground {} vs {}",
            vals[0],
            2.0 * e1(1)
        );
        // (1,2) and (2,1) are degenerate.
        assert!((vals[1] - (e1(1) + e1(2))).abs() < 1e-6);
        assert!((vals[2] - (e1(1) + e1(2))).abs() < 1e-6);
        // Boundary points carry a huge potential, so exactly (nx-2)(ny-2) levels stay low.
        let low = vals.iter().filter(|&&v| v < 1e6).count();
        assert_eq!(low, (nx - 2) * (ny - 2));
    }

    #[test]
    fn schrodinger_3d_box_ground_state() {
        let (n, l) = (5usize, 1.0);
        let d = l / (n as f64 - 1.0);
        let (vals, _) = solve_3d_schrodinger(|_, _, _| 0.0, (0.0, l, 0.0, l, 0.0, l), (n, n, n))
            .unwrap_or_else(|e| panic!("{e}"));
        let vals = sorted(vals);
        let e1 = (1.0 - (PI / (n as f64 - 1.0)).cos()) / (d * d);
        assert!(
            (vals[0] - 3.0 * e1).abs() < 1e-5,
            "ground {} vs {}",
            vals[0],
            3.0 * e1
        );
        // Next level is triply degenerate: 2 e1(1) + e1(2).
        let e2 = (1.0 - (2.0 * PI / (n as f64 - 1.0)).cos()) / (d * d);
        for v in &vals[1..4] {
            assert!((v - (2.0 * e1 + e2)).abs() < 1e-5, "{v}");
        }
    }

    #[test]
    fn heat_equation_matches_the_decaying_sine_mode() {
        // u_t = alpha u_xx, u(x,0) = sin(pi x), u(0)=u(1)=0  =>  u = exp(-alpha pi^2 t) sin(pi x)
        let (alpha, nx, nt) = (0.5, 101, 100);
        let res = solve_heat_equation_1d_crank_nicolson(
            &|x| (PI * x).sin(),
            alpha,
            (0.0, 1.0),
            nx,
            (0.0, 0.2),
            nt,
        )
        .unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(res.len(), nt + 1);
        let decay = (-alpha * PI * PI * 0.2).exp();
        for (i, &u) in res[nt].iter().enumerate() {
            let x = i as f64 / (nx as f64 - 1.0);
            assert!((u - decay * (PI * x).sin()).abs() < 2e-4, "x = {x}: {u}");
        }
        // Boundary values stay zero and the maximum decreases monotonically.
        assert!(res.iter().all(|u| u[0] == 0.0 && u[nx - 1] == 0.0));
        let peaks: Vec<f64> = res
            .iter()
            .map(|u| u.iter().cloned().fold(f64::MIN, f64::max))
            .collect();
        assert!(peaks.windows(2).all(|w| w[1] <= w[0] + 1e-15));
    }

    #[test]
    fn heat_equation_validates_arguments() {
        assert!(
            solve_heat_equation_1d_crank_nicolson(&|_| 1.0, 1.0, (0.0, 1.0), 2, (0.0, 1.0), 10)
                .is_err()
        );
        assert!(
            solve_heat_equation_1d_crank_nicolson(&|_| 1.0, 1.0, (0.0, 1.0), 10, (0.0, 1.0), 0)
                .is_err()
        );
    }

    #[test]
    fn wave_equation_standing_wave() {
        // u(x, 0) = sin(pi x), u_t = 0, c = 1  =>  u(x, t) = sin(pi x) cos(pi t)
        let n = 101;
        let dx = 1.0 / (n as f64 - 1.0);
        let u0: Vec<f64> = (0..n).map(|i| (PI * i as f64 * dx).sin()).collect();
        let ut = vec![0.0; n];
        let dt = 0.005;
        let steps = 200;
        let snaps =
            solve_wave_equation_1d(&u0, &ut, 1.0, dx, dt, steps).unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(snaps.len(), steps + 1);
        for &(k, tol) in &[(100usize, 2e-3), (200usize, 3e-3)] {
            let t = k as f64 * dt;
            for (i, &u) in snaps[k].iter().enumerate() {
                let want = (PI * i as f64 * dx).sin() * (PI * t).cos();
                assert!((u - want).abs() < tol, "step {k}, i = {i}: {u} vs {want}");
            }
        }
        assert!(snaps[1..].iter().all(|s| s[0] == 0.0 && s[n - 1] == 0.0));
    }

    #[test]
    fn wave_equation_validates_arguments() {
        let u = vec![0.0; 10];
        assert!(solve_wave_equation_1d(&u, &vec![0.0; 9], 1.0, 0.1, 0.05, 5).is_err());
        assert!(solve_wave_equation_1d(&[0.0, 0.0], &[0.0, 0.0], 1.0, 0.1, 0.05, 5).is_err());
        // CFL number 2 is rejected, 1 is accepted.
        assert!(solve_wave_equation_1d(&u, &u, 1.0, 0.1, 0.2, 5).is_err());
        assert!(solve_wave_equation_1d(&u, &u, 1.0, 0.1, 0.1, 5).is_ok());
    }

    #[test]
    fn ising_model_returns_a_spin_lattice() {
        // steps == size keeps the (currently size-coupled) initialisation in range.
        let lattice = simulate_ising_model(8, 2.0, 8);
        assert_eq!(lattice.len(), 8);
        for row in &lattice {
            assert_eq!(row.len(), 8);
            assert!(
                row.iter().all(|&s| s == 1 || s == -1),
                "spins must be +-1: {row:?}"
            );
        }
    }

    #[test]
    fn ising_model_with_more_steps_than_rows() {
        let lattice = simulate_ising_model(8, 2.0, 1000);
        assert_eq!(lattice.len(), 8);
        assert!(lattice.iter().flatten().all(|&s| s == 1 || s == -1));
    }

    #[test]
    fn ising_model_with_fewer_steps_than_rows_has_only_unit_spins() {
        let lattice = simulate_ising_model(8, 2.0, 3);
        assert!(lattice.iter().flatten().all(|&s| s == 1 || s == -1));
    }

    #[test]
    fn ising_model_cools_towards_low_energy_and_handles_degenerate_sizes() {
        // At T = 0.05 almost only energy-lowering flips are accepted, so
        // the energy per site must fall far below the random-start value (~0).
        let n = 8;
        let lat = simulate_ising_model(n, 0.05, 100_000);
        let mut e = 0.0;
        for i in 0..n {
            for j in 0..n {
                let s = f64::from(lat[i][j]);
                e -= s * f64::from(lat[(i + 1) % n][j] + lat[i][(j + 1) % n]);
            }
        }
        let per_site = e / ((n * n) as f64);
        assert!(per_site < -0.5, "E/N = {per_site}");
        assert!(simulate_ising_model(0, 1.0, 10).is_empty());
        assert!(
            simulate_ising_model(4, 1.0, 0)
                .iter()
                .flatten()
                .all(|&s| s.abs() == 1)
        );
    }

    proptest! {
        #![proptest_config(cfg())]

        #[test]
        fn prop_relativistic_velocity_addition_never_exceeds_c(a in -0.999..0.999f64, b in -0.999..0.999f64) {
            let c = SPEED_OF_LIGHT;
            let u = relativistic_velocity_addition(a * c, b * c);
            prop_assert!(u.abs() < c);
            prop_assert!((u - relativistic_velocity_addition(b * c, a * c)).abs() < 1e-6);
        }

        #[test]
        fn prop_energy_momentum_relation(m in 0.1..100.0f64, beta in 0.0..0.99f64) {
            let c = SPEED_OF_LIGHT;
            let v = beta * c;
            let e = relativistic_total_energy(m, v);
            let p = relativistic_momentum(m, v);
            prop_assert!(((e * e - (p * c).powi(2)) / (m * c * c).powi(2) - 1.0).abs() < 1e-9);
        }

        #[test]
        fn prop_heat_solution_obeys_the_maximum_principle(
            a in 0.0..1.0f64, b in 0.0..1.0f64, c in 0.0..1.0f64,
        ) {
            let init = move |x: f64| a * (PI * x).sin() + b * (2.0 * PI * x).sin() + c * (3.0 * PI * x).sin();
            let res = solve_heat_equation_1d_crank_nicolson(&init, 0.1, (0.0, 1.0), 41, (0.0, 0.5), 25)
                .map_err(TestCaseError::fail)?;
            let peak0 = res[0].iter().fold(0.0f64, |m, v| m.max(v.abs()));
            for u in &res {
                prop_assert!(u.iter().all(|v| v.abs() <= peak0 + 1e-9));
            }
        }

        #[test]
        fn prop_uniform_gravity_trajectory_is_quadratic(z0 in -5.0..5.0f64, vz in -5.0..5.0f64, tend in 0.5..3.0f64) {
            let traj = simulate_particle_motion(|_, _| [0.0, 0.0, -9.81], 1.0, (0.0, 0.0, z0), (0.0, 0.0, vz), (0.0, tend), 10)
                .map_err(TestCaseError::fail)?;
            let z = traj[10][2];
            prop_assert!((z - (z0 + vz * tend - 0.5 * 9.81 * tend * tend)).abs() < 1e-9);
        }
    }

    // Keep the `Matrix` import used even when the eigenvector checks above change.
    #[allow(dead_code)]
    fn _assert_matrix_type(_: &Matrix<f64>) {}
}

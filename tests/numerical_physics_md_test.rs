//! Comprehensive tests for numerical MD module.
//!
//! Tests for particles, potentials, thermodynamics, and analysis functions.

use rssn::numerical::physics_md::*;

// ============================================================================
// Particle Tests
// ============================================================================

#[test]

fn test_particle_new() {

    let p = Particle::new(0, 1.0, vec![0.0, 0.0, 0.0], vec![1.0, 0.0, 0.0]);

    assert_eq!(p.id, 0);

    assert_eq!(p.mass, 1.0);

    assert_eq!(p.position, vec![0.0, 0.0, 0.0]);

    assert_eq!(p.velocity, vec![1.0, 0.0, 0.0]);

    assert_eq!(p.charge, 0.0);
}

#[test]

fn test_particle_with_charge() {

    let p = Particle::with_charge(0, 1.0, vec![0.0, 0.0, 0.0], vec![0.0, 0.0, 0.0], 1.0);

    assert_eq!(p.charge, 1.0);
}

#[test]

fn test_particle_kinetic_energy() {

    let p = Particle::new(0, 2.0, vec![0.0, 0.0, 0.0], vec![3.0, 4.0, 0.0]);

    // KE = 0.5 * 2 * (3² + 4² + 0²) = 1 * 25 = 25
    assert!((p.kinetic_energy() - 25.0).abs() < 1e-10);
}

#[test]

fn test_particle_momentum() {

    let p = Particle::new(0, 2.0, vec![0.0, 0.0, 0.0], vec![1.0, 2.0, 3.0]);

    let mom = p.momentum();

    assert_eq!(mom, vec![2.0, 4.0, 6.0]);
}

#[test]

fn test_particle_speed() {

    let p = Particle::new(0, 1.0, vec![0.0, 0.0, 0.0], vec![3.0, 4.0, 0.0]);

    assert!((p.speed() - 5.0).abs() < 1e-10);
}

#[test]

fn test_particle_distance_to() {

    let p1 = Particle::new(0, 1.0, vec![0.0, 0.0, 0.0], vec![0.0, 0.0, 0.0]);

    let p2 = Particle::new(1, 1.0, vec![3.0, 4.0, 0.0], vec![0.0, 0.0, 0.0]);

    let d = p1
        .distance_to(&p2)
        .unwrap();

    assert!((d - 5.0).abs() < 1e-10);
}

// ============================================================================
// Lennard-Jones Tests
// ============================================================================

#[test]

fn test_lennard_jones_equilibrium() {

    // At r = 2^(1/6) * σ, force should be zero (equilibrium)
    let r_eq = 2.0_f64.powf(1.0 / 6.0);

    let p1 = Particle::new(0, 1.0, vec![0.0, 0.0, 0.0], vec![0.0, 0.0, 0.0]);

    let p2 = Particle::new(1, 1.0, vec![r_eq, 0.0, 0.0], vec![0.0, 0.0, 0.0]);

    let (potential, force) = lennard_jones_interaction(&p1, &p2, 1.0, 1.0).unwrap();

    // Potential at equilibrium is -ε
    assert!((potential + 1.0).abs() < 1e-10);

    // Force should be approximately zero
    let force_mag: f64 = force
        .iter()
        .map(|f| f * f)
        .sum::<f64>()
        .sqrt();

    assert!(force_mag.abs() < 1e-10);
}

#[test]

fn test_lennard_jones_repulsive() {

    // At r < equilibrium, potential is positive and force is repulsive
    let p1 = Particle::new(0, 1.0, vec![0.0, 0.0, 0.0], vec![0.0, 0.0, 0.0]);

    let p2 = Particle::new(1, 1.0, vec![0.9, 0.0, 0.0], vec![0.0, 0.0, 0.0]);

    let (potential, _force) = lennard_jones_interaction(&p1, &p2, 1.0, 1.0).unwrap();

    assert!(potential > 0.0);
}

// ============================================================================
// Morse Potential Tests
// ============================================================================

#[test]

fn test_morse_at_equilibrium() {

    let p1 = Particle::new(0, 1.0, vec![0.0, 0.0, 0.0], vec![0.0, 0.0, 0.0]);

    let p2 = Particle::new(1, 1.0, vec![1.0, 0.0, 0.0], vec![0.0, 0.0, 0.0]);

    let (potential, _force) = morse_interaction(&p1, &p2, 1.0, 1.0, 1.0).unwrap();

    // At r = re, potential should be 0
    assert!(potential.abs() < 1e-10);
}

#[test]

fn test_morse_stretched() {

    let p1 = Particle::new(0, 1.0, vec![0.0, 0.0, 0.0], vec![0.0, 0.0, 0.0]);

    let p2 = Particle::new(1, 1.0, vec![2.0, 0.0, 0.0], vec![0.0, 0.0, 0.0]);

    let (potential, _force) = morse_interaction(&p1, &p2, 1.0, 1.0, 1.0).unwrap();

    // At r > re, potential should be positive
    assert!(potential > 0.0);
}

// ============================================================================
// Harmonic Potential Tests
// ============================================================================

#[test]

fn test_harmonic_at_equilibrium() {

    let p1 = Particle::new(0, 1.0, vec![0.0, 0.0, 0.0], vec![0.0, 0.0, 0.0]);

    let p2 = Particle::new(1, 1.0, vec![1.0, 0.0, 0.0], vec![0.0, 0.0, 0.0]);

    let (potential, _force) = harmonic_interaction(&p1, &p2, 100.0, 1.0).unwrap();

    // At r = r0, potential should be 0
    assert!(potential.abs() < 1e-10);
}

#[test]

fn test_harmonic_stretched() {

    let p1 = Particle::new(0, 1.0, vec![0.0, 0.0, 0.0], vec![0.0, 0.0, 0.0]);

    let p2 = Particle::new(1, 1.0, vec![1.5, 0.0, 0.0], vec![0.0, 0.0, 0.0]);

    let (potential, _force) = harmonic_interaction(&p1, &p2, 100.0, 1.0).unwrap();

    // V = 0.5 * 100 * (1.5 - 1.0)² = 50 * 0.25 = 12.5
    assert!((potential - 12.5).abs() < 1e-10);
}

// ============================================================================
// Coulomb Potential Tests
// ============================================================================

#[test]

fn test_coulomb_like_charges() {

    let p1 = Particle::with_charge(0, 1.0, vec![0.0, 0.0, 0.0], vec![0.0, 0.0, 0.0], 1.0);

    let p2 = Particle::with_charge(1, 1.0, vec![1.0, 0.0, 0.0], vec![0.0, 0.0, 0.0], 1.0);

    let (potential, _force) = coulomb_interaction(&p1, &p2, 1.0).unwrap();

    // V = k * q1 * q2 / r = 1 * 1 * 1 / 1 = 1
    assert!((potential - 1.0).abs() < 1e-10);
}

#[test]

fn test_coulomb_opposite_charges() {

    let p1 = Particle::with_charge(0, 1.0, vec![0.0, 0.0, 0.0], vec![0.0, 0.0, 0.0], 1.0);

    let p2 = Particle::with_charge(1, 1.0, vec![1.0, 0.0, 0.0], vec![0.0, 0.0, 0.0], -1.0);

    let (potential, _force) = coulomb_interaction(&p1, &p2, 1.0).unwrap();

    // V = k * q1 * q2 / r = 1 * 1 * (-1) / 1 = -1
    assert!((potential + 1.0).abs() < 1e-10);
}

// ============================================================================
// System Properties Tests
// ============================================================================

#[test]

fn test_total_kinetic_energy() {

    let particles = vec![
        Particle::new(0, 1.0, vec![0.0, 0.0, 0.0], vec![1.0, 0.0, 0.0]),
        Particle::new(1, 1.0, vec![1.0, 0.0, 0.0], vec![0.0, 1.0, 0.0]),
    ];

    let ke = total_kinetic_energy(&particles);

    // KE = 0.5 * 1 * 1 + 0.5 * 1 * 1 = 1.0
    assert!((ke - 1.0).abs() < 1e-10);
}

#[test]

fn test_total_momentum() {

    let particles = vec![
        Particle::new(0, 1.0, vec![0.0, 0.0, 0.0], vec![1.0, 0.0, 0.0]),
        Particle::new(1, 1.0, vec![1.0, 0.0, 0.0], vec![-1.0, 0.0, 0.0]),
    ];

    let mom = total_momentum(&particles).unwrap();

    // Opposite velocities should cancel
    assert!(mom[0].abs() < 1e-10);
}

#[test]

fn test_center_of_mass() {

    let particles = vec![
        Particle::new(0, 1.0, vec![0.0, 0.0, 0.0], vec![0.0, 0.0, 0.0]),
        Particle::new(1, 1.0, vec![2.0, 0.0, 0.0], vec![0.0, 0.0, 0.0]),
    ];

    let com = center_of_mass(&particles).unwrap();

    assert!((com[0] - 1.0).abs() < 1e-10);
}

#[test]

fn test_temperature() {

    // With 3D particles and KE = (3/2) * N * T, T = 2 * KE / (3 * N)
    let particles = vec![Particle::new(
        0,
        1.0,
        vec![0.0, 0.0, 0.0],
        vec![1.0, 1.0, 1.0],
    )];

    let t = temperature(&particles);

    // KE = 0.5 * 1 * 3 = 1.5, T = 2 * 1.5 / 3 = 1.0
    assert!((t - 1.0).abs() < 1e-10);
}

// ============================================================================
// Thermostat Tests
// ============================================================================

#[test]

fn test_velocity_rescale() {

    let mut particles = vec![Particle::new(
        0,
        1.0,
        vec![0.0, 0.0, 0.0],
        vec![1.0, 1.0, 1.0],
    )];

    let initial_temp = temperature(&particles);

    let target_temp = 2.0 * initial_temp;

    velocity_rescale(&mut particles, target_temp);

    let final_temp = temperature(&particles);

    assert!((final_temp - target_temp).abs() < 1e-10);
}

#[test]

fn test_berendsen_thermostat() {

    let mut particles = vec![Particle::new(
        0,
        1.0,
        vec![0.0, 0.0, 0.0],
        vec![1.0, 1.0, 1.0],
    )];

    let initial_temp = temperature(&particles);

    let target_temp = 2.0 * initial_temp;

    berendsen_thermostat(&mut particles, target_temp, 0.1, 0.01);

    let final_temp = temperature(&particles);

    // Temperature should move toward target
    assert!(final_temp > initial_temp);
}

// ============================================================================
// Periodic Boundary Conditions Tests
// ============================================================================

#[test]

fn test_apply_pbc() {

    let position = vec![11.0, -1.0, 5.0];

    let box_size = vec![10.0, 10.0, 10.0];

    let wrapped = apply_pbc(&position, &box_size);

    assert!((wrapped[0] - 1.0).abs() < 1e-10);

    assert!((wrapped[1] - 9.0).abs() < 1e-10);

    assert!((wrapped[2] - 5.0).abs() < 1e-10);
}

#[test]

fn test_minimum_image_distance() {

    let r = vec![8.0, -8.0, 0.0];

    let box_size = vec![10.0, 10.0, 10.0];

    let r_mic = minimum_image_distance(&r, &box_size);

    assert!((r_mic[0] + 2.0).abs() < 1e-10);

    assert!((r_mic[1] - 2.0).abs() < 1e-10);
}

// ============================================================================
// Analysis Tests
// ============================================================================

#[test]

fn test_mean_square_displacement() {

    let initial = vec![Particle::new(
        0,
        1.0,
        vec![0.0, 0.0, 0.0],
        vec![0.0, 0.0, 0.0],
    )];

    let current = vec![Particle::new(
        0,
        1.0,
        vec![3.0, 4.0, 0.0],
        vec![0.0, 0.0, 0.0],
    )];

    let msd = mean_square_displacement(&initial, &current);

    // MSD = 3² + 4² = 25
    assert!((msd - 25.0).abs() < 1e-10);
}

#[test]

fn test_radial_distribution_function() {

    let particles = vec![
        Particle::new(0, 1.0, vec![0.0, 0.0, 0.0], vec![0.0, 0.0, 0.0]),
        Particle::new(1, 1.0, vec![1.0, 0.0, 0.0], vec![0.0, 0.0, 0.0]),
    ];

    let box_size = vec![10.0, 10.0, 10.0];

    let (r_values, g_r) = radial_distribution_function(&particles, &box_size, 10, 5.0);

    assert_eq!(r_values.len(), 10);

    assert_eq!(g_r.len(), 10);
}

// ============================================================================
// Lattice Creation Tests
// ============================================================================

#[test]

fn test_create_cubic_lattice() {

    let particles = create_cubic_lattice(3, 1.0, 1.0);

    assert_eq!(particles.len(), 27); // 3³ = 27
}

#[test]

fn test_create_fcc_lattice() {

    let particles = create_fcc_lattice(2, 1.0, 1.0);

    assert_eq!(particles.len(), 32); // 4 × 2³ = 32
}

#[test]

fn test_remove_com_velocity() {

    let mut particles = vec![
        Particle::new(0, 1.0, vec![0.0, 0.0, 0.0], vec![1.0, 0.0, 0.0]),
        Particle::new(1, 1.0, vec![1.0, 0.0, 0.0], vec![3.0, 0.0, 0.0]),
    ];

    remove_com_velocity(&mut particles).unwrap();

    let mom = total_momentum(&particles).unwrap();

    assert!(mom[0].abs() < 1e-10);
}

// ============================================================================
// Property Tests
// ============================================================================

mod proptests {

    use super::*;
    use proptest::prelude::*;

    proptest! {
        #[test]
        fn prop_kinetic_energy_positive(vx in -10.0..10.0f64, vy in -10.0..10.0f64, vz in -10.0..10.0f64) {
            let p = Particle::new(0, 1.0, vec![0.0, 0.0, 0.0], vec![vx, vy, vz]);
            prop_assert!(p.kinetic_energy() >= 0.0);
        }

        #[test]
        fn prop_temperature_positive(v in 0.1..10.0f64) {
            let particles = vec![
                Particle::new(0, 1.0, vec![0.0, 0.0, 0.0], vec![v, 0.0, 0.0]),
            ];
            prop_assert!(temperature(&particles) >= 0.0);
        }

        #[test]
        fn prop_pbc_in_box(x in -100.0..100.0f64) {
            let position = vec![x];
            let box_size = vec![10.0];
            let wrapped = apply_pbc(&position, &box_size);
            prop_assert!(wrapped[0] >= 0.0 && wrapped[0] < 10.0);
        }

        #[test]
        fn prop_minimum_image_bounded(dx in -50.0..50.0f64) {
            let r = vec![dx];
            let box_size = vec![10.0];
            let r_mic = minimum_image_distance(&r, &box_size);
            prop_assert!(r_mic[0] >= -5.0 && r_mic[0] <= 5.0);
        }
    }
}

//! Molecular-dynamics helpers (ported from `numerical_physics_md_test.rs`).
//!
//! Tests for particles, potentials, thermodynamics, and analysis functions.

use rssn::sim::physics_md::*;

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

    let d = p1.distance_to(&p2).unwrap_or_else(|e| panic!("{e}"));

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

    let (potential, force) =
        lennard_jones_interaction(&p1, &p2, 1.0, 1.0).unwrap_or_else(|e| panic!("{e}"));

    // Potential at equilibrium is -ε
    assert!((potential + 1.0).abs() < 1e-10);

    // Force should be approximately zero
    let force_mag: f64 = force.iter().map(|f| f * f).sum::<f64>().sqrt();

    assert!(force_mag.abs() < 1e-10);
}

#[test]

fn test_lennard_jones_repulsive() {
    // At r < equilibrium, potential is positive and force is repulsive
    let p1 = Particle::new(0, 1.0, vec![0.0, 0.0, 0.0], vec![0.0, 0.0, 0.0]);

    let p2 = Particle::new(1, 1.0, vec![0.9, 0.0, 0.0], vec![0.0, 0.0, 0.0]);

    let (potential, _force) =
        lennard_jones_interaction(&p1, &p2, 1.0, 1.0).unwrap_or_else(|e| panic!("{e}"));

    assert!(potential > 0.0);
}

// ============================================================================
// Morse Potential Tests
// ============================================================================

#[test]

fn test_morse_at_equilibrium() {
    let p1 = Particle::new(0, 1.0, vec![0.0, 0.0, 0.0], vec![0.0, 0.0, 0.0]);

    let p2 = Particle::new(1, 1.0, vec![1.0, 0.0, 0.0], vec![0.0, 0.0, 0.0]);

    let (potential, _force) =
        morse_interaction(&p1, &p2, 1.0, 1.0, 1.0).unwrap_or_else(|e| panic!("{e}"));

    // At r = re, potential should be 0
    assert!(potential.abs() < 1e-10);
}

#[test]

fn test_morse_stretched() {
    let p1 = Particle::new(0, 1.0, vec![0.0, 0.0, 0.0], vec![0.0, 0.0, 0.0]);

    let p2 = Particle::new(1, 1.0, vec![2.0, 0.0, 0.0], vec![0.0, 0.0, 0.0]);

    let (potential, _force) =
        morse_interaction(&p1, &p2, 1.0, 1.0, 1.0).unwrap_or_else(|e| panic!("{e}"));

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

    let (potential, _force) =
        harmonic_interaction(&p1, &p2, 100.0, 1.0).unwrap_or_else(|e| panic!("{e}"));

    // At r = r0, potential should be 0
    assert!(potential.abs() < 1e-10);
}

#[test]

fn test_harmonic_stretched() {
    let p1 = Particle::new(0, 1.0, vec![0.0, 0.0, 0.0], vec![0.0, 0.0, 0.0]);

    let p2 = Particle::new(1, 1.0, vec![1.5, 0.0, 0.0], vec![0.0, 0.0, 0.0]);

    let (potential, _force) =
        harmonic_interaction(&p1, &p2, 100.0, 1.0).unwrap_or_else(|e| panic!("{e}"));

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

    let (potential, _force) = coulomb_interaction(&p1, &p2, 1.0).unwrap_or_else(|e| panic!("{e}"));

    // V = k * q1 * q2 / r = 1 * 1 * 1 / 1 = 1
    assert!((potential - 1.0).abs() < 1e-10);
}

#[test]

fn test_coulomb_opposite_charges() {
    let p1 = Particle::with_charge(0, 1.0, vec![0.0, 0.0, 0.0], vec![0.0, 0.0, 0.0], 1.0);

    let p2 = Particle::with_charge(1, 1.0, vec![1.0, 0.0, 0.0], vec![0.0, 0.0, 0.0], -1.0);

    let (potential, _force) = coulomb_interaction(&p1, &p2, 1.0).unwrap_or_else(|e| panic!("{e}"));

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

    let mom = total_momentum(&particles).unwrap_or_else(|e| panic!("{e}"));

    // Opposite velocities should cancel
    assert!(mom[0].abs() < 1e-10);
}

#[test]

fn test_center_of_mass() {
    let particles = vec![
        Particle::new(0, 1.0, vec![0.0, 0.0, 0.0], vec![0.0, 0.0, 0.0]),
        Particle::new(1, 1.0, vec![2.0, 0.0, 0.0], vec![0.0, 0.0, 0.0]),
    ];

    let com = center_of_mass(&particles).unwrap_or_else(|e| panic!("{e}"));

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

    remove_com_velocity(&mut particles).unwrap_or_else(|e| panic!("{e}"));

    let mom = total_momentum(&particles).unwrap_or_else(|e| panic!("{e}"));

    assert!(mom[0].abs() < 1e-10);
}

// ============================================================================
// Property Tests
// ============================================================================

mod proptests {

    use proptest::prelude::*;

    use super::*;

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

// ============================================================================
// Added: physics checks (forces are gradients of potentials, Verlet, RDF, ...)
// ============================================================================

mod strengthened {
    use proptest::prelude::*;
    use proptest::test_runner::RngSeed;
    use rssn::sim::physics_md::*;

    fn cfg() -> ProptestConfig {
        ProptestConfig {
            cases: 48,
            rng_seed: RngSeed::Fixed(0x5EED),
            failure_persistence: None,
            ..ProptestConfig::default()
        }
    }

    fn at(
        x: f64,
        y: f64,
        z: f64,
    ) -> Particle {
        Particle::new(0, 1.0, vec![x, y, z], vec![0.0; 3])
    }

    /// Finite-difference check that `force == -grad_{p1} V` for a pair potential.
    fn check_force_is_minus_gradient(
        name: &str,
        pair: impl Fn(&Particle, &Particle) -> Result<(f64, Vec<f64>), String>,
        p2: &Particle,
        pos1: [f64; 3],
    ) -> Result<(), String> {
        let h = 1e-6;
        let (_, force) = pair(&at(pos1[0], pos1[1], pos1[2]), p2)?;
        for axis in 0..3 {
            let (mut plus, mut minus) = (pos1, pos1);
            plus[axis] += h;
            minus[axis] -= h;
            let vp = pair(&at(plus[0], plus[1], plus[2]), p2)?.0;
            let vm = pair(&at(minus[0], minus[1], minus[2]), p2)?.0;
            let grad = (vp - vm) / (2.0 * h);
            let scale = force[axis].abs().max(grad.abs()).max(1.0);
            if (force[axis] + grad).abs() > 1e-5 * scale {
                return Err(format!(
                    "{name}: F[{axis}] = {} but -dV/dx = {}",
                    force[axis], -grad
                ));
            }
        }
        Ok(())
    }

    #[test]
    fn particle_accessors_and_dimension_checks() {
        let p = Particle::new(7, 3.0, vec![1.0, 2.0], vec![0.5, -0.5]);
        assert_eq!((p.id, p.mass), (7, 3.0));
        assert_eq!(p.force, vec![0.0, 0.0]);
        assert!((p.speed() - 0.5f64.hypot(0.5)).abs() < 1e-15);
        assert!(
            p.distance_to(&at(0.0, 0.0, 0.0)).is_err(),
            "2D and 3D positions cannot be compared"
        );
        assert!(lennard_jones_interaction(&p, &at(0.0, 0.0, 0.0), 1.0, 1.0).is_err());
    }

    #[test]
    fn lennard_jones_reference_values() {
        // r = sigma: V = 0, F = 24 epsilon / sigma (repulsive, along +r from p2 to p1)
        let (v, f) = lennard_jones_interaction(&at(1.0, 0.0, 0.0), &at(0.0, 0.0, 0.0), 1.0, 1.0)
            .unwrap_or_else(|e| panic!("{e}"));
        assert!(v.abs() < 1e-12);
        assert!(
            (f[0] - 24.0).abs() < 1e-9 && f[1] == 0.0 && f[2] == 0.0,
            "{f:?}"
        );
        // r = 2: V = 4 (1/4096 - 1/64)
        let (v, _) = lennard_jones_interaction(&at(2.0, 0.0, 0.0), &at(0.0, 0.0, 0.0), 1.0, 1.0)
            .unwrap_or_else(|e| panic!("{e}"));
        assert!((v - 4.0 * (1.0 / 4096.0 - 1.0 / 64.0)).abs() < 1e-12);
        // Newton's third law: swapping the particles reverses the force.
        let (_, f12) = lennard_jones_interaction(&at(1.3, 0.2, 0.0), &at(0.0, 0.0, 0.1), 0.7, 1.1)
            .unwrap_or_else(|e| panic!("{e}"));
        let (_, f21) = lennard_jones_interaction(&at(0.0, 0.0, 0.1), &at(1.3, 0.2, 0.0), 0.7, 1.1)
            .unwrap_or_else(|e| panic!("{e}"));
        for k in 0..3 {
            assert!((f12[k] + f21[k]).abs() < 1e-12);
        }
        // Coincident particles: infinite potential, zero force.
        let (v, f) = lennard_jones_interaction(&at(0.0, 0.0, 0.0), &at(0.0, 0.0, 0.0), 1.0, 1.0)
            .unwrap_or_else(|e| panic!("{e}"));
        assert_eq!((v, f), (f64::INFINITY, vec![0.0; 3]));
    }

    #[test]
    fn morse_reference_values() {
        // V(re + ln 2 / a) = De / 4 ; well depth is De below the dissociation limit.
        let (a, de, re) = (1.5, 2.0, 1.2);
        let r = re + 2f64.ln() / a;
        let (v, _) = morse_interaction(&at(r, 0.0, 0.0), &at(0.0, 0.0, 0.0), de, a, re)
            .unwrap_or_else(|e| panic!("{e}"));
        assert!((v - de / 4.0).abs() < 1e-12);
        let (v, _) = morse_interaction(&at(1e3, 0.0, 0.0), &at(0.0, 0.0, 0.0), de, a, re)
            .unwrap_or_else(|e| panic!("{e}"));
        assert!((v - de).abs() < 1e-12);
        // Equilibrium: zero force.
        let (_, f) = morse_interaction(&at(re, 0.0, 0.0), &at(0.0, 0.0, 0.0), de, a, re)
            .unwrap_or_else(|e| panic!("{e}"));
        assert!(f.iter().all(|c| c.abs() < 1e-12));
    }

    #[test]
    fn morse_force_is_minus_the_gradient_of_the_potential() {
        let p2 = at(0.0, 0.0, 0.0);
        let morse = |a: &Particle, b: &Particle| morse_interaction(a, b, 1.0, 1.0, 1.0);
        for pos in [[1.5, 0.0, 0.0], [0.8, 0.3, 0.1], [2.0, 1.0, -0.5]] {
            if let Err(e) = check_force_is_minus_gradient("morse", morse, &p2, pos) {
                panic!("{e}");
            }
        }
        // Stretched bond pulls p1 back toward p2 (negative x component).
        let (_, f) = morse_interaction(&at(1.5, 0.0, 0.0), &p2, 1.0, 1.0, 1.0)
            .unwrap_or_else(|e| panic!("{e}"));
        assert!(f[0] < 0.0, "expected attraction, got {f:?}");
    }

    #[test]
    fn morse_force_matches_analytic_value_and_is_antisymmetric() {
        // r = re + 0.5, De = a = 1: F = -2 (1 - e^-0.5) e^-0.5 = -0.47730...
        let (_, f) =
            morse_interaction(&at(1.5, 0.0, 0.0), &at(0.0, 0.0, 0.0), 1.0, 1.0, 1.0).unwrap();
        let e = (-0.5f64).exp();
        assert!((f[0] + 2.0 * (1.0 - e) * e).abs() < 1e-12);
        let (_, g) =
            morse_interaction(&at(0.0, 0.0, 0.0), &at(1.5, 0.0, 0.0), 1.0, 1.0, 1.0).unwrap();
        assert!((f[0] + g[0]).abs() < 1e-12, "Newton's third law");
        // Compressed bond repels: force on p1 points away from p2.
        let (_, h) =
            morse_interaction(&at(0.7, 0.0, 0.0), &at(0.0, 0.0, 0.0), 1.0, 1.0, 1.0).unwrap();
        assert!(h[0] > 0.0);
    }

    #[test]
    fn lennard_jones_harmonic_coulomb_forces_are_minus_the_gradient() {
        let p2 = Particle::with_charge(1, 1.0, vec![0.1, -0.2, 0.3], vec![0.0; 3], -2.0);
        let mk = |q: f64| {
            move |a: &Particle, b: &Particle| {
                let mut a = a.clone();
                a.charge = q;
                coulomb_interaction(&a, b, 1.7)
            }
        };
        for pos in [[1.4, 0.2, 0.1], [0.9, 0.9, 0.9], [-1.0, 0.5, 0.7]] {
            let lj = |a: &Particle, b: &Particle| lennard_jones_interaction(a, b, 0.8, 1.1);
            let ha = |a: &Particle, b: &Particle| harmonic_interaction(a, b, 40.0, 1.3);
            let ss = |a: &Particle, b: &Particle| soft_sphere_interaction(a, b, 0.5, 2.5, 6);
            for (name, res) in [
                ("lj", check_force_is_minus_gradient("lj", lj, &p2, pos)),
                (
                    "harmonic",
                    check_force_is_minus_gradient("harmonic", ha, &p2, pos),
                ),
                (
                    "soft sphere",
                    check_force_is_minus_gradient("soft sphere", ss, &p2, pos),
                ),
                (
                    "coulomb",
                    check_force_is_minus_gradient("coulomb", mk(3.0), &p2, pos),
                ),
            ] {
                if let Err(e) = res {
                    panic!("{name} at {pos:?}: {e}");
                }
            }
        }
    }

    #[test]
    fn coulomb_harmonic_and_soft_sphere_reference_values() {
        // Like charges repel: force on p1 points away from p2.
        let a = Particle::with_charge(0, 1.0, vec![2.0, 0.0, 0.0], vec![0.0; 3], 3.0);
        let b = Particle::with_charge(1, 1.0, vec![0.0, 0.0, 0.0], vec![0.0; 3], 2.0);
        let (v, f) = coulomb_interaction(&a, &b, 5.0).unwrap_or_else(|e| panic!("{e}"));
        assert!((v - 5.0 * 6.0 / 2.0).abs() < 1e-12);
        assert!((f[0] - 5.0 * 6.0 / 4.0).abs() < 1e-12);
        // Harmonic bond stretched by 0.5 pulls p1 back: F = -k dr.
        let (v, f) = harmonic_interaction(&at(1.5, 0.0, 0.0), &at(0.0, 0.0, 0.0), 100.0, 1.0)
            .unwrap_or_else(|e| panic!("{e}"));
        assert!((v - 12.5).abs() < 1e-12 && (f[0] + 50.0).abs() < 1e-12);
        // Soft sphere is zero beyond sigma, epsilon (sigma/r)^n inside.
        let (v, f) = soft_sphere_interaction(&at(3.0, 0.0, 0.0), &at(0.0, 0.0, 0.0), 1.0, 2.0, 12)
            .unwrap_or_else(|e| panic!("{e}"));
        assert_eq!((v, f), (0.0, vec![0.0; 3]));
        let (v, f) = soft_sphere_interaction(&at(1.0, 0.0, 0.0), &at(0.0, 0.0, 0.0), 1.0, 2.0, 12)
            .unwrap_or_else(|e| panic!("{e}"));
        assert!((v - 4096.0).abs() < 1e-9 && (f[0] - 12.0 * 4096.0).abs() < 1e-6);
    }

    #[test]
    fn velocity_verlet_reproduces_a_harmonic_dimer_and_conserves_energy() {
        // Two unit masses joined by a spring (k = 8, r0 = 1): reduced mass 1/2 => omega = sqrt(k / mu) = 4.
        let (k, r0) = (8.0, 1.0);
        let mut ps = vec![
            Particle::new(0, 1.0, vec![0.0, 0.0, 0.0], vec![0.0; 3]),
            Particle::new(1, 1.0, vec![1.25, 0.0, 0.0], vec![0.0; 3]),
        ];
        let calc = |p: &mut Vec<Particle>| -> Result<(), String> {
            let (_, f) = harmonic_interaction(&p[0], &p[1], k, r0)?;
            p[0].force = f.clone();
            p[1].force = f.iter().map(|c| -c).collect();
            Ok(())
        };
        let (dt, steps) = (1e-3, 2000);
        let traj =
            integrate_velocity_verlet(&mut ps, dt, steps, calc).unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(traj.len(), steps + 1);
        let omega = 4.0f64;
        let energy = |s: &[Particle]| {
            let r = s[0].distance_to(&s[1]).unwrap_or(f64::NAN);
            total_kinetic_energy(s) + 0.5 * k * (r - r0) * (r - r0)
        };
        let e0 = energy(&traj[0]);
        for (i, s) in traj.iter().enumerate() {
            let t = i as f64 * dt;
            let r = s[1].position[0] - s[0].position[0];
            assert!(
                (r - (r0 + 0.25 * (omega * t).cos())).abs() < 2e-5,
                "step {i}: r = {r}"
            );
            assert!(
                (energy(s) - e0).abs() < 1e-4 * e0,
                "energy drift at step {i}"
            );
        }
        // Total momentum is exactly conserved and the centre of mass stays put.
        for s in &traj {
            assert!(
                total_momentum(s)
                    .unwrap_or_default()
                    .iter()
                    .all(|p| p.abs() < 1e-12)
            );
            assert!((center_of_mass(s).unwrap_or_default()[0] - 0.625).abs() < 1e-12);
        }
    }

    #[test]
    fn velocity_verlet_propagates_force_calculator_errors() {
        let mut ps = vec![at(0.0, 0.0, 0.0)];
        let r = integrate_velocity_verlet(&mut ps, 0.1, 3, |_| Err("boom".to_string()));
        assert_eq!(r.err().as_deref(), Some("boom"));
    }

    #[test]
    fn lennard_jones_dimer_conserves_energy_near_the_minimum() {
        let mut ps = vec![
            Particle::new(0, 1.0, vec![0.0, 0.0, 0.0], vec![0.0; 3]),
            Particle::new(1, 1.0, vec![1.05, 0.0, 0.0], vec![0.0; 3]),
        ];
        let calc = |p: &mut Vec<Particle>| -> Result<(), String> {
            let (_, f) = lennard_jones_interaction(&p[0], &p[1], 1.0, 1.0)?;
            p[0].force = f.clone();
            p[1].force = f.iter().map(|c| -c).collect();
            Ok(())
        };
        let traj =
            integrate_velocity_verlet(&mut ps, 2e-3, 3000, calc).unwrap_or_else(|e| panic!("{e}"));
        let energy = |s: &[Particle]| {
            let (v, _) =
                lennard_jones_interaction(&s[0], &s[1], 1.0, 1.0).unwrap_or((f64::NAN, vec![]));
            total_kinetic_energy(s) + v
        };
        let e0 = energy(&traj[0]);
        assert!(e0 < 0.0, "bound state");
        assert!(
            traj.iter()
                .all(|s| (energy(s) - e0).abs() < 1e-4 * e0.abs())
        );
        // The bond oscillates around the minimum 2^(1/6).
        let rs: Vec<f64> = traj
            .iter()
            .map(|s| s[1].position[0] - s[0].position[0])
            .collect();
        let (lo, hi) = (
            rs.iter().cloned().fold(f64::MAX, f64::min),
            rs.iter().cloned().fold(f64::MIN, f64::max),
        );
        assert!(lo < 2f64.powf(1.0 / 6.0) && hi > 2f64.powf(1.0 / 6.0));
    }

    #[test]
    fn system_properties() {
        let ps = vec![
            Particle::new(0, 2.0, vec![0.0, 0.0, 0.0], vec![1.0, 0.0, 0.0]),
            Particle::new(1, 6.0, vec![4.0, 0.0, 0.0], vec![0.0, 1.0, 0.0]),
        ];
        assert_eq!(total_momentum(&ps), Ok(vec![2.0, 6.0, 0.0]));
        assert_eq!(center_of_mass(&ps), Ok(vec![3.0, 0.0, 0.0]));
        assert_eq!(total_kinetic_energy(&ps), 1.0 + 3.0);
        // T = 2 KE / (3 N) with k_B = 1
        assert!((temperature(&ps) - 2.0 * 4.0 / 6.0).abs() < 1e-15);
        assert_eq!(temperature(&[]), 0.0);
        assert!(total_momentum(&[]).is_err() && center_of_mass(&[]).is_err());
        // Ideal-gas pressure P V = N T (+ virial).
        assert!((pressure(&ps, 10.0, 0.0) - 2.0 * temperature(&ps) / 10.0).abs() < 1e-15);
        assert!((pressure(&ps, 10.0, 5.0) - (2.0 * temperature(&ps) + 5.0) / 10.0).abs() < 1e-15);
        assert_eq!(pressure(&ps, 0.0, 1.0), 0.0);
    }

    #[test]
    fn thermostats_hit_or_approach_the_target() {
        let mk = || {
            vec![
                Particle::new(0, 1.0, vec![0.0; 3], vec![1.0, 2.0, 2.0]),
                Particle::new(1, 2.0, vec![1.0; 3], vec![-1.0, 0.5, 0.0]),
            ]
        };
        let mut ps = mk();
        velocity_rescale(&mut ps, 3.5);
        assert!((temperature(&ps) - 3.5).abs() < 1e-12);
        // Berendsen with dt = tau reduces to a full rescale; with dt << tau it moves only slightly.
        let mut full = mk();
        let t0 = temperature(&full);
        berendsen_thermostat(&mut full, 3.5, 0.1, 0.1);
        assert!((temperature(&full) - 3.5).abs() < 1e-12);
        let mut weak = mk();
        berendsen_thermostat(&mut weak, 3.5, 1.0, 0.01);
        let t1 = temperature(&weak);
        assert!(t1 > t0 && t1 < 3.5, "{t0} -> {t1}");
        // A zero-temperature system stays frozen.
        let mut cold = vec![at(0.0, 0.0, 0.0)];
        velocity_rescale(&mut cold, 5.0);
        berendsen_thermostat(&mut cold, 5.0, 1.0, 0.1);
        assert_eq!(cold[0].velocity, vec![0.0; 3]);
    }

    #[test]
    fn maxwell_boltzmann_initialisation_is_seeded_and_exact() {
        let make = || create_cubic_lattice(3, 1.0, 2.0);
        let mut a = make();
        let mut b = make();
        let mut c = make();
        initialize_velocities_maxwell_boltzmann(&mut a, 1.5, 42);
        initialize_velocities_maxwell_boltzmann(&mut b, 1.5, 42);
        initialize_velocities_maxwell_boltzmann(&mut c, 1.5, 43);
        assert!(
            a.iter().zip(&b).all(|(p, q)| p.velocity == q.velocity),
            "same seed => same velocities"
        );
        assert!(
            a.iter().zip(&c).any(|(p, q)| p.velocity != q.velocity),
            "different seed => different velocities"
        );
        assert!((temperature(&a) - 1.5).abs() < 1e-12);
        assert!(
            total_momentum(&a)
                .unwrap_or_default()
                .iter()
                .all(|p| p.abs() < 1e-10)
        );
    }

    #[test]
    fn periodic_boundary_helpers() {
        assert_eq!(
            apply_pbc(&[10.0, 25.0, -25.0, 0.0], &[10.0; 4]),
            vec![0.0, 5.0, 5.0, 0.0]
        );
        assert_eq!(apply_pbc(&[1.0, 2.0], &[10.0, 1.5]), vec![1.0, 0.5]);
        assert_eq!(
            minimum_image_distance(&[9.0, -9.0, 5.0, 4.0], &[10.0; 4]),
            vec![-1.0, 1.0, 5.0, 4.0]
        );
        assert_eq!(minimum_image_distance(&[23.0], &[10.0]), vec![3.0]);
    }

    #[test]
    fn mean_square_displacement_and_com_removal() {
        let a = vec![at(0.0, 0.0, 0.0), at(1.0, 1.0, 1.0)];
        let b = vec![at(1.0, 0.0, 0.0), at(1.0, 3.0, 1.0)];
        assert!((mean_square_displacement(&a, &b) - (1.0 + 4.0) / 2.0).abs() < 1e-15);
        assert_eq!(mean_square_displacement(&a, &b[..1]), 0.0);
        assert_eq!(mean_square_displacement(&[], &[]), 0.0);
        let mut ps = vec![
            Particle::new(0, 1.0, vec![0.0; 3], vec![1.0, 0.0, 0.0]),
            Particle::new(1, 3.0, vec![1.0; 3], vec![5.0, 2.0, 0.0]),
        ];
        remove_com_velocity(&mut ps).unwrap_or_else(|e| panic!("{e}"));
        // v_com = (16/4, 6/4, 0): velocities shift by it and momentum vanishes.
        assert_eq!(ps[0].velocity, vec![-3.0, -1.5, 0.0]);
        assert_eq!(ps[1].velocity, vec![1.0, 0.5, 0.0]);
        assert!(
            total_momentum(&ps)
                .unwrap_or_default()
                .iter()
                .all(|p| p.abs() < 1e-12)
        );
        assert!(remove_com_velocity(&mut []).is_ok());
    }

    #[test]
    fn lattices_have_the_expected_geometry() {
        let cubic = create_cubic_lattice(3, 1.5, 2.0);
        assert_eq!(cubic.len(), 27);
        assert!(
            cubic
                .iter()
                .enumerate()
                .all(|(i, p)| p.id == i && p.mass == 2.0 && p.velocity == vec![0.0; 3])
        );
        assert_eq!(cubic[26].position, vec![3.0, 3.0, 3.0]);
        let fcc = create_fcc_lattice(2, 2.0, 1.0);
        assert_eq!(fcc.len(), 32);
        // Nearest-neighbour distance a / sqrt(2), with 12 neighbours in the periodic box.
        let box_size = vec![4.0; 3];
        let d_nn = 2.0 / 2f64.sqrt();
        let neighbours = fcc
            .iter()
            .skip(1)
            .filter(|p| {
                let r = minimum_image_distance(
                    &fcc[0]
                        .position
                        .iter()
                        .zip(&p.position)
                        .map(|(a, b)| a - b)
                        .collect::<Vec<_>>(),
                    &box_size,
                );
                (r.iter().map(|x| x * x).sum::<f64>().sqrt() - d_nn).abs() < 1e-9
            })
            .count();
        assert_eq!(neighbours, 12);
        assert!(create_cubic_lattice(0, 1.0, 1.0).is_empty());
    }

    #[test]
    fn radial_distribution_of_a_single_pair_and_of_an_ideal_gas() {
        // Two particles at distance 1 in a 10^3 box: one pair in bin 2 (dr = 0.5).
        let pair = vec![at(0.0, 0.0, 0.0), at(1.0, 0.0, 0.0)];
        let (r, g) = radial_distribution_function(&pair, &[10.0; 3], 10, 5.0);
        assert_eq!(r.len(), 10);
        assert!((r[0] - 0.25).abs() < 1e-12 && (r[9] - 4.75).abs() < 1e-12);
        let rho = 2.0 / 1000.0;
        let shell = 4.0 / 3.0 * std::f64::consts::PI * (1.5f64.powi(3) - 1.0);
        assert!((g[2] - 1.0 / (rho * shell)).abs() < 1e-9, "g[2] = {}", g[2]);
        assert!(g.iter().enumerate().all(|(i, &v)| i == 2 || v == 0.0));
        assert!(
            radial_distribution_function(&pair[..1], &[10.0; 3], 10, 5.0)
                .0
                .is_empty()
        );

        // Uniformly distributed points (seeded LCG): g(r) fluctuates around 1.
        let mut state = 12345u64;
        let mut next = || {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (state >> 11) as f64 / (1u64 << 53) as f64 * 10.0
        };
        let gas: Vec<Particle> = (0..600).map(|_| at(next(), next(), next())).collect();
        let (_, g) = radial_distribution_function(&gas, &[10.0; 3], 10, 5.0);
        for &v in &g[3..] {
            assert!((v - 1.0).abs() < 0.15, "g = {v}");
        }
    }

    proptest! {
        #![proptest_config(cfg())]

        #[test]
        fn prop_pair_forces_obey_newtons_third_law(
            x in 0.7..3.0f64, y in -1.0..1.0f64, z in -1.0..1.0f64,
        ) {
            let (a, b) = (at(x, y, z), at(0.0, 0.1, -0.1));
            type Pair = fn(&Particle, &Particle) -> Result<(f64, Vec<f64>), String>;
            let pairs: [Pair; 3] = [
                |p, q| lennard_jones_interaction(p, q, 1.0, 1.0),
                |p, q| harmonic_interaction(p, q, 5.0, 1.5),
                |p, q| soft_sphere_interaction(p, q, 1.0, 3.0, 8),
            ];
            for pair in pairs {
                let (v_ab, f_ab) = pair(&a, &b).map_err(TestCaseError::fail)?;
                let (v_ba, f_ba) = pair(&b, &a).map_err(TestCaseError::fail)?;
                prop_assert!((v_ab - v_ba).abs() < 1e-9 * (1.0 + v_ab.abs()), "potential must be symmetric");
                for k in 0..3 {
                    prop_assert!((f_ab[k] + f_ba[k]).abs() < 1e-9 * (1.0 + f_ab[k].abs()));
                }
            }
        }

        #[test]
        fn prop_minimum_image_is_shorter_than_half_the_box_and_periodic(dx in -100.0..100.0f64, l in 1.0..20.0f64) {
            let m = minimum_image_distance(&[dx], &[l])[0];
            prop_assert!(m.abs() <= l / 2.0 + 1e-9);
            // Differs from dx by an integer number of box lengths.
            let k = (dx - m) / l;
            prop_assert!((k - k.round()).abs() < 1e-9);
        }

        #[test]
        fn prop_temperature_scales_with_velocity_squared(v in 0.1..10.0f64, s in 0.1..5.0f64) {
            let mk = |scale: f64| vec![Particle::new(0, 1.5, vec![0.0; 3], vec![v * scale, 0.0, 0.0])];
            let ratio = temperature(&mk(s)) / temperature(&mk(1.0));
            prop_assert!((ratio - s * s).abs() < 1e-9 * (1.0 + s * s));
        }
    }
}

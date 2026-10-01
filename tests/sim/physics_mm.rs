//! Meshless SPH method (ported from `physics_mm_test.rs`).

use assert_approx_eq::assert_approx_eq;
use proptest::prelude::*;
use proptest::test_runner::RngSeed;
use rssn::sim::physics_mm::*;
use std::f64::consts::PI;

fn cfg() -> ProptestConfig {
    ProptestConfig {
        cases: 24,
        rng_seed: RngSeed::Fixed(0x5EED),
        failure_persistence: None,
        ..ProptestConfig::default()
    }
}

fn particle(
    x: f64,
    y: f64,
) -> Particle {
    Particle {
        pos: Vector2D::new(x, y),
        vel: Vector2D::default(),
        force: Vector2D::default(),
        density: 0.0,
        pressure: 0.0,
        mass: 1.0,
    }
}

fn system(
    h: f64,
    particles: Vec<Particle>,
) -> SPHSystem {
    SPHSystem {
        particles,
        poly6: Poly6Kernel::new(h),
        spiky: SpikyKernel::new(h),
        gravity: Vector2D::new(0.0, 0.0),
        viscosity: 0.1,
        gas_const: 1000.0,
        rest_density: 1000.0,
        bounds: Vector2D::new(10.0, 10.0),
    }
}

#[test]
fn vector2d_arithmetic() {
    let a = Vector2D::new(1.0, 2.0);
    let b = Vector2D::new(3.0, -4.0);
    let s = a + b;
    assert_eq!((s.x, s.y), (4.0, -2.0));
    let d = a - b;
    assert_eq!((d.x, d.y), (-2.0, 6.0));
    let m = a * 2.0;
    assert_eq!((m.x, m.y), (2.0, 4.0));
    let q = b / 2.0;
    assert_eq!((q.x, q.y), (1.5, -2.0));
    let z = Vector2D::default();
    assert_eq!((z.x, z.y), (0.0, 0.0));
}

#[test]
fn kernel_constructors_store_normalisations() {
    let h = 0.5;
    let p = Poly6Kernel::new(h);
    assert_approx_eq!(p.h_sq, 0.25);
    assert_approx_eq!(p.factor, 315.0 / (64.0 * PI * h.powi(9)), 1e-9 * p.factor);
    let s = SpikyKernel::new(h);
    assert_approx_eq!(s.h, h);
    assert_approx_eq!(s.factor, -45.0 / (PI * h.powi(6)), 1e-9 * s.factor.abs());
}

#[test]
fn isolated_particle_density_is_the_kernel_peak() {
    let h = 0.5;
    let mut sys = system(h, vec![particle(1.0, 1.0)]);
    sys.compute_density_pressure();
    let peak = 315.0 / (64.0 * PI * h.powi(3));
    assert_approx_eq!(sys.particles[0].density, peak, 1e-9 * peak);
    // Density below the rest density gives zero pressure.
    assert_eq!(sys.particles[0].pressure, 0.0);
}

#[test]
fn density_pressure_for_a_close_pair_is_symmetric_and_larger() {
    let h = 0.5;
    let mut single = system(h, vec![particle(0.0, 0.0)]);
    single.compute_density_pressure();
    let mut pair = system(h, vec![particle(0.0, 0.0), particle(0.05, 0.0)]);
    pair.compute_density_pressure();
    assert!(pair.particles[0].density > single.particles[0].density);
    assert_approx_eq!(pair.particles[0].density, pair.particles[1].density, 1e-9);
    let p = pair.particles[0].pressure;
    assert_approx_eq!(
        p,
        1000.0 * (pair.particles[0].density - 1000.0).max(0.0),
        1e-9 * (1.0 + p)
    );
}

#[test]
fn particles_outside_the_support_do_not_interact() {
    let h = 0.5;
    let mut sys = system(h, vec![particle(0.0, 0.0), particle(2.0, 0.0)]);
    sys.compute_density_pressure();
    sys.compute_forces();
    for p in &sys.particles {
        assert_eq!((p.force.x, p.force.y), (0.0, 0.0));
    }
}

#[test]
fn pair_forces_are_equal_and_opposite_without_gravity() {
    let h = 0.5;
    let mut sys = system(h, vec![particle(1.0, 1.0), particle(1.2, 1.0)]);
    sys.rest_density = 1.0; // ensure non-zero pressure
    sys.compute_density_pressure();
    sys.compute_forces();
    let (f0, f1) = (sys.particles[0].force, sys.particles[1].force);
    assert!(f0.x.abs() > 0.0);
    assert_approx_eq!(f0.x + f1.x, 0.0, 1e-9 * (1.0 + f0.x.abs()));
    assert_approx_eq!(f0.y + f1.y, 0.0, 1e-9);
    // Pressure pushes the pair apart: the left particle is pushed toward -x.
    assert!(f0.x < 0.0 && f1.x > 0.0);
}

#[test]
fn integrate_clamps_to_the_bounds_and_damps_velocity() {
    let mut sys = system(0.5, vec![particle(0.05, 5.0)]);
    sys.bounds = Vector2D::new(1.0, 1.0);
    sys.particles[0].density = 1.0;
    sys.particles[0].vel = Vector2D::new(-10.0, 0.0);
    sys.integrate(0.1);
    let p = &sys.particles[0];
    assert_eq!(p.pos.x, 0.0);
    assert_approx_eq!(p.vel.x, 5.0);
    assert_eq!(p.pos.y, 1.0);
}

#[test]
fn free_fall_accelerates_with_gravity() {
    let mut sys = system(0.1, vec![particle(5.0, 5.0)]);
    sys.gravity = Vector2D::new(0.0, -9.8);
    let dt = 0.01;
    for _ in 0..10 {
        sys.update(dt);
    }
    let p = &sys.particles[0];
    // Velocity after 10 explicit steps of constant acceleration.
    assert_approx_eq!(p.vel.y, -9.8 * 0.1, 1e-9);
    assert!(p.pos.y < 5.0);
    assert_approx_eq!(p.pos.x, 5.0);
}

#[test]
fn dam_break_scenario_stays_inside_the_box_and_spreads() {
    let res = simulate_dam_break_2d_scenario();
    assert_eq!(res.len(), 200);
    assert!(res.iter().all(|&(x, y)| x.is_finite()
        && y.is_finite()
        && (0.0..=4.0).contains(&x)
        && (0.0..=4.0).contains(&y)));
    let max_x = res.iter().map(|p| p.0).fold(f64::MIN, f64::max);
    assert!(max_x > 0.72, "fluid did not spread: max x = {max_x}");
}

proptest! {
    #![proptest_config(cfg())]

    #[test]
    fn prop_update_keeps_particles_finite_and_inside(dt in 0.001f64..0.01) {
        let mut sys = system(0.1, vec![particle(0.5, 0.5)]);
        sys.gravity = Vector2D::new(0.0, -9.8);
        sys.bounds = Vector2D::new(1.0, 1.0);
        sys.particles[0].vel = Vector2D::new(0.1, 0.1);
        sys.update(dt);
        for p in &sys.particles {
            prop_assert!(p.pos.x.is_finite() && p.pos.y.is_finite());
            prop_assert!((0.0..=1.0).contains(&p.pos.x));
            prop_assert!((0.0..=1.0).contains(&p.pos.y));
        }
    }

    #[test]
    fn prop_density_is_symmetric_for_two_equal_particles(d in 0.0f64..0.45) {
        let mut sys = system(0.5, vec![particle(1.0, 1.0), particle(1.0 + d, 1.0)]);
        sys.compute_density_pressure();
        let (a, b) = (sys.particles[0].density, sys.particles[1].density);
        prop_assert!((a - b).abs() <= 1e-9 * a.abs());
    }
}

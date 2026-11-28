//! Ordered neighbor reuse for the reduced-unit SPH model.

use rssn::sim::physics_mm::Particle;
use rssn::sim::physics_mm::Poly6Kernel;
use rssn::sim::physics_mm::SPHSystem;
use rssn::sim::physics_mm::SpikyKernel;
use rssn::sim::physics_mm::Vector2D;

fn system(n: usize) -> SPHSystem {
    let particles = (0..n)
        .map(|i| {
            Particle {
                pos: Vector2D::new((i % 16) as f64 * 0.0625, (i / 16) as f64 * 0.0625),
                vel: Vector2D::new((i % 7) as f64 * 0.007, (i % 11) as f64 * -0.01),
                force: Vector2D::default(),
                density: 1.0,
                pressure: 0.0,
                mass: 0.075 + (i % 11) as f64 * 0.002,
            }
        })
        .collect();
    SPHSystem {
        particles,
        poly6: Poly6Kernel::new(0.125),
        spiky: SpikyKernel::new(0.125),
        gravity: Vector2D::new(0.1, -0.5),
        viscosity: 0.003,
        gas_const: 0.05,
        rest_density: 30.0,
        bounds: Vector2D::new(8.0, 8.0),
    }
}

fn weight(
    square: f64,
    kernel: &Poly6Kernel,
) -> f64 {
    if square >= kernel.h_sq {
        0.0
    } else {
        let delta = kernel.h_sq - square;
        kernel.factor * delta * delta * delta
    }
}

fn square(delta: Vector2D) -> f64 {
    delta.x.mul_add(delta.x, delta.y * delta.y)
}

fn direct_density(sys: &SPHSystem) -> Vec<(f64, f64)> {
    sys.particles
        .iter()
        .map(|outer| {
            let density = sys.particles.iter().fold(0.0, |sum, inner| {
                sum + inner.mass * weight(square(outer.pos - inner.pos), &sys.poly6)
            });
            (
                density,
                sys.gas_const * (density - sys.rest_density).max(0.0),
            )
        })
        .collect()
}

fn direct_forces(
    sys: &SPHSystem,
    states: &[(f64, f64)],
) -> Vec<Vector2D> {
    sys.particles
        .iter()
        .enumerate()
        .map(|(i, outer)| {
            let mut force = Vector2D::default();
            for (j, inner) in sys.particles.iter().enumerate() {
                let delta = outer.pos - inner.pos;
                let r = square(delta).sqrt();
                if i != j && r < sys.spiky.h {
                    let gradient = if r == 0.0 {
                        Vector2D::default()
                    } else {
                        let q = sys.spiky.h - r;
                        delta * (sys.spiky.factor * q * q / r)
                    };
                    let pressure = f64::midpoint(states[i].1, states[j].1) / states[j].0;
                    force = force - gradient * pressure;
                    force = force
                        + (inner.vel - outer.vel)
                            * (sys.viscosity * weight(square(delta), &sys.poly6));
                }
            }
            force + sys.gravity * states[i].0
        })
        .collect()
}

fn bits(
    actual: f64,
    expected: f64,
) {
    if expected.is_nan() {
        assert!(actual.is_nan());
    } else {
        assert_eq!(actual.to_bits(), expected.to_bits());
    }
}

fn assert_snapshot(
    sys: &SPHSystem,
    states: &[(f64, f64)],
    forces: &[Vector2D],
) {
    for ((p, &(rho, pressure)), f) in sys.particles.iter().zip(states).zip(forces) {
        bits(p.density, rho);
        bits(p.pressure, pressure);
        bits(p.force.x, f.x);
        bits(p.force.y, f.y);
    }
}

#[test]
fn support_edges_and_both_cutoff_orders_agree_with_the_pair_sum() {
    for h in [0.0625_f64, 0.125, 0.25] {
        let mut sys = system(1440);
        sys.spiky = SpikyKernel::new(h);
        sys.particles[1].pos = Vector2D::new(h.next_down(), 0.0);
        sys.particles[2].pos = Vector2D::new(h, 0.0);
        sys.particles[3].pos = Vector2D::new(h.next_up(), 0.0);
        sys.particles[4].pos = sys.particles[0].pos;
        sys.particles[5].pos.x = -0.03;
        sys.particles.swap(7, 123);
        sys.particles.rotate_left(53);
        let states = direct_density(&sys);
        let forces = direct_forces(&sys, &states);

        let mut separate = sys.clone();
        separate.compute_density_pressure();
        separate.compute_forces();
        assert_snapshot(&separate, &states, &forces);
        sys.compute_density_and_forces();
        assert_snapshot(&sys, &states, &forces);
    }
}

#[test]
fn each_step_rebuilds_from_moved_particles_and_matches_explicit_integration() {
    let one = rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .build()
        .unwrap();
    let two = rayon::ThreadPoolBuilder::new()
        .num_threads(2)
        .build()
        .unwrap();
    let mut serial = system(1152);
    let mut parallel = serial.clone();
    let dt = 0.0005;
    for step in 0..3 {
        if step == 1 {
            serial.particles[5].pos.x += 0.45;
            parallel.particles[5].pos = serial.particles[5].pos;
        }
        let states = direct_density(&serial);
        let forces = direct_forces(&serial, &states);
        let expected: Vec<_> = serial
            .particles
            .iter()
            .zip(&states)
            .zip(&forces)
            .map(|((p, &(density, _)), &force)| {
                let mut velocity = p.vel + force * (dt / density);
                let mut position = p.pos + velocity * dt;
                if position.x < 0.0 || position.x > serial.bounds.x {
                    velocity.x *= -0.5;
                    position.x = position.x.clamp(0.0, serial.bounds.x);
                }
                if position.y < 0.0 || position.y > serial.bounds.y {
                    velocity.y *= -0.5;
                    position.y = position.y.clamp(0.0, serial.bounds.y);
                }
                (position, velocity)
            })
            .collect();
        one.install(|| serial.update(dt));
        two.install(|| parallel.update(dt));
        for sys in [&serial, &parallel] {
            assert_snapshot(sys, &states, &forces);
            for (p, &(pos, vel)) in sys.particles.iter().zip(&expected) {
                bits(p.pos.x, pos.x);
                bits(p.pos.y, pos.y);
                bits(p.vel.x, vel.x);
                bits(p.vel.y, vel.y);
            }
        }
    }
}

#[test]
fn preserves_nonfinite_masses_and_underflowed_force_distances() {
    for translated in [2.0_f64.powi(44), 2.0_f64.powi(48)] {
        let mut sys = system(1152);
        sys.particles[0].pos.x = translated;
        sys.particles[1].pos.x = translated + 0.0625;
        let states = direct_density(&sys);
        let forces = direct_forces(&sys, &states);
        assert!(states[0].0.is_finite());
        sys.compute_density_and_forces();
        assert_snapshot(&sys, &states, &forces);
    }

    let mut nonfinite = system(1152);
    nonfinite.particles[1151].mass = f64::INFINITY;
    let states = direct_density(&nonfinite);
    let forces = direct_forces(&nonfinite, &states);
    assert!(states[0].0.is_nan());
    nonfinite.compute_density_and_forces();
    assert_snapshot(&nonfinite, &states, &forces);

    let mut tiny = system(1152);
    tiny.spiky = SpikyKernel {
        h: 1e-175,
        factor: -1.0,
    };
    for (i, p) in tiny.particles.iter_mut().enumerate() {
        p.pos = Vector2D::new(i as f64 * 1e-166, 0.0);
        p.vel = Vector2D::new(i as f64, 0.0);
    }
    let states = vec![(1.0, 0.0); tiny.particles.len()];
    let forces = direct_forces(&tiny, &states);
    assert!(forces[0].x > 1.0);
    tiny.compute_forces();
    assert_snapshot(&tiny, &states, &forces);
}

#[test]
fn modern_density_reuses_its_points_and_falls_back_for_extreme_open_positions() {
    use rssn::sim::sph::KernelKind;
    use rssn::sim::sph::SphBoundary;
    use rssn::sim::sph::SphDensityPlan;

    let h = 0.125;
    let positions = [[0.0], [0.0625], [0.1875]];
    let mut moved = positions;
    moved[1][0] += 0.25;
    let expanded = [[0.0], [0.0625], [0.1875], [2.0_f64.powi(52)]];
    let build =
        |pos: &[[f64; 1]]| SphDensityPlan::new(pos, [0.0], [1.0], [SphBoundary::Open], h).unwrap();
    let initial = build(&positions);
    let moved_plan = build(&moved);
    let fallback = build(&expanded);
    for mass in [[1.0, 2.0, 3.0, 4.0], [0.125, 0.0, 7.0, 0.25]] {
        let kernel = KernelKind::WendlandC2;
        let original = initial.evaluate(&mass[..3], &[h; 3], kernel).unwrap();
        let rebuild = moved_plan.evaluate(&mass[..3], &[h; 3], kernel).unwrap();
        let fallback = fallback.evaluate(&mass, &[h; 4], kernel).unwrap();
        for (field, points) in [(&original, &positions), (&rebuild, &moved)] {
            for (actual, &[x]) in field.iter().zip(points) {
                let direct = points
                    .iter()
                    .zip(mass)
                    .fold(0.0, |s, (&[y], m)| s + m * kernel.w(1, (x - y).abs(), h));
                assert!((actual / direct - 1.0).abs() < 1e-14);
            }
        }
        for (&cached, &all_pairs) in original.iter().zip(&fallback) {
            assert!((cached / all_pairs - 1.0).abs() < 1e-14);
        }
        bits(fallback[3], mass[3] * kernel.w(1, 0.0, h));
    }
}

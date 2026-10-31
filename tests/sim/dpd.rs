//! `sim::dpd`: the Groot–Warren DPD fluid. Direct checks of the pair force
//! (conservative repulsion, minimum image, virial), momentum conservation,
//! the free-particle and ideal-gas limits, thermostat relaxation and
//! equipartition, determinism and input validation.

use rssn::sim::dpd::Dpd;
use rssn::sim::dpd::DpdParams;

fn params(
    dim: usize,
    box_len: f64,
    a: f64,
    gamma: f64,
    kt: f64,
) -> DpdParams {
    DpdParams {
        dim,
        box_len: [box_len, box_len, box_len],
        a,
        gamma,
        kt,
        rc: 1.0,
        mass: 1.0,
    }
}

/// Two particles at the given positions at rest, `kt = 0` (no random force)
/// and `gamma = 0`: the force is the conservative repulsion alone.
fn pair(
    a: f64,
    x1: [f64; 3],
    x2: [f64; 3],
    box_len: f64,
) -> Dpd {
    let mut s = Dpd::new(2, params(3, box_len, a, 0.0, 0.0), 1).unwrap();
    s.pos = vec![x1, x2];
    s.vel = vec![[0.0; 3]; 2];
    s
}

#[test]
fn constructor_validates_and_initialises() {
    let p = params(3, 4.0, 25.0, 4.5, 1.0);
    assert!(Dpd::new(1, p, 1).is_err(), "fewer than two particles");
    assert!(
        Dpd::new(10, DpdParams { dim: 1, ..p }, 1).is_err(),
        "dimension 1"
    );
    assert!(
        Dpd::new(10, DpdParams { dim: 4, ..p }, 1).is_err(),
        "dimension 4"
    );
    assert!(
        Dpd::new(10, params(3, 2.9, 25.0, 4.5, 1.0), 1).is_err(),
        "box below three cut-offs"
    );
    let s = Dpd::new(500, p, 7).unwrap();
    assert_eq!(s.len(), 500);
    assert!(!s.is_empty());
    assert!((s.density() - 500.0 / 64.0).abs() < 1e-14);
    // Zero total momentum to round-off, positions inside the box, and the
    // Maxwell temperature to sampling accuracy (relative sd = sqrt(2/(d N))).
    let m = s.momentum();
    assert!(m.iter().all(|x| x.abs() < 1e-11), "{m:?}");
    assert!(
        s.pos
            .iter()
            .all(|r| r.iter().all(|&x| (0.0..4.0).contains(&x)))
    );
    assert!(
        (s.temperature() - 1.0).abs() < 5.0 * (2.0 / 1500.0_f64).sqrt(),
        "{}",
        s.temperature()
    );
    // In two dimensions the third component is unused.
    let two = Dpd::new(200, params(2, 6.0, 25.0, 4.5, 1.0), 3).unwrap();
    assert!(two.pos.iter().chain(&two.vel).all(|v| v[2] == 0.0));
    assert!((two.density() - 200.0 / 36.0).abs() < 1e-14);
}

/// One step to refresh the force stored at construction (computed at the
/// random initial positions), then the velocities are reset to rest.
fn settle(
    s: &mut Dpd,
    dt: f64,
) {
    s.step(dt);
    s.vel.fill([0.0; 3]);
}

#[test]
fn pair_force_is_the_soft_conservative_repulsion() {
    // a (1 - r/rc) along the line of centres: the velocity gained in a short
    // step from rest is dt F / m.
    let (a, r, dt) = (25.0, 0.4, 1e-4);
    let mut s = pair(a, [1.0, 1.0, 1.0], [1.0 + r, 1.0, 1.0], 3.0);
    settle(&mut s, dt);
    s.step(dt);
    let want = a * (1.0 - r);
    let (g0, g1) = (s.vel[0][0] / dt, s.vel[1][0] / dt);
    assert!(
        (g0 + g1).abs() < 1e-9 * want,
        "equal and opposite: {g0} {g1}"
    );
    assert!(
        g0 < 0.0 && g1 > 0.0,
        "repulsion pushes the particles apart: {g0} {g1}"
    );
    assert!((g1 - want).abs() < 1e-3 * want, "{g1} vs {want}");
    // The force falls linearly to zero at the cut-off: r = 0.7 gives 0.3 a.
    let mut t = pair(a, [1.0, 1.0, 1.0], [1.7, 1.0, 1.0], 3.0);
    settle(&mut t, dt);
    t.step(dt);
    assert!(
        (t.vel[1][0] / dt - 0.3 * a).abs() < 1e-3 * a,
        "{}",
        t.vel[1][0] / dt
    );
    // Transverse components stay zero, also for an oblique pair.
    assert!(
        s.vel
            .iter()
            .all(|v| v[1].abs() < 1e-12 && v[2].abs() < 1e-12)
    );
    let mut o = pair(a, [1.0, 1.0, 1.0], [1.3, 1.3, 1.3], 3.0);
    settle(&mut o, dt);
    o.step(dt);
    let distance = 0.3 * 3.0_f64.sqrt();
    let force = a * (1.0 - distance);
    for axis in 0..3 {
        assert!(
            (o.vel[1][axis] / dt - force / 3.0_f64.sqrt()).abs() < 2e-3 * force,
            "axis {axis}"
        );
    }
    // Beyond the cut-off there is no force at all.
    let mut far = pair(a, [1.0, 1.0, 1.0], [2.2, 1.0, 1.0], 4.0);
    settle(&mut far, dt);
    far.run(dt, 3);
    assert!(far.vel.iter().all(|v| v.iter().all(|x| x.abs() < 1e-12)));
}

#[test]
fn virial_pressure_of_a_pair_matches_the_conservative_virial() {
    // p = rho T_kin + (a w r) / (d V) with the last force evaluation.
    let (a, r, box_len) = (25.0, 0.4, 3.0);
    let mut s = pair(a, [1.0, 1.0, 1.0], [1.0 + r, 1.0, 1.0], box_len);
    settle(&mut s, 1e-6);
    s.step(1e-6);
    let volume = box_len.powi(3);
    let want = a * (1.0 - r) * r / (3.0 * volume);
    assert!(
        (s.pressure() - want).abs() < 1e-4 * want,
        "{} vs {want}",
        s.pressure()
    );
}

#[test]
fn minimum_image_acts_across_the_periodic_boundary() {
    let dt = 1e-4;
    // Particles 0.2 apart through the boundary of a box of edge 3: they
    // repel through it, so the one at x = 0.1 is pushed to +x.
    let mut s = pair(25.0, [0.1, 1.0, 1.0], [2.9, 1.0, 1.0], 3.0);
    settle(&mut s, dt);
    s.step(dt);
    assert!(s.vel[0][0] > 0.0 && s.vel[1][0] < 0.0, "{:?}", s.vel);
    // The same pair without a boundary in between pushes the other way
    // round, with the same magnitude 25 (1 - 0.2).
    let mut t = pair(25.0, [1.0, 1.0, 1.0], [1.2, 1.0, 1.0], 3.0);
    settle(&mut t, dt);
    t.step(dt);
    assert!(t.vel[0][0] < 0.0 && t.vel[1][0] > 0.0);
    assert!(
        (t.vel[0][0] + s.vel[0][0]).abs() < 1e-9,
        "{} vs {}",
        t.vel[0][0],
        s.vel[0][0]
    );
    assert!(
        (s.vel[0][0] / dt - 20.0).abs() < 1e-2,
        "{}",
        s.vel[0][0] / dt
    );
    // Positions are wrapped back into the box.
    let mut w = pair(0.0, [0.01, 1.0, 1.0], [1.5, 1.0, 1.0], 3.0);
    w.vel[0][0] = -2.0;
    w.run(0.05, 5);
    assert!(
        w.pos
            .iter()
            .all(|r| r.iter().all(|&x| (0.0..3.0).contains(&x))),
        "{:?}",
        w.pos
    );
    assert!(
        w.pos[0][0] > 2.0,
        "wrapped to the far side: {}",
        w.pos[0][0]
    );
}

#[test]
fn free_particles_without_forces_move_ballistically() {
    // a = 0, gamma = 0, kt = 0: no force at all, x(t) = x0 + v t (mod L).
    let mut s = Dpd::new(60, params(3, 4.0, 0.0, 0.0, 0.0), 11).unwrap();
    let mut rng = rssn::kernels::random::Rng::new(5);
    for v in &mut s.vel {
        *v = [rng.normal(), rng.normal(), rng.normal()];
    }
    let start = s.pos.clone();
    let velocity = s.vel.clone();
    s.run(0.02, 50);
    assert!((s.t - 1.0).abs() < 1e-12);
    for i in 0..s.len() {
        for a in 0..3 {
            assert_eq!(s.vel[i][a], velocity[i][a], "velocity is conserved exactly");
            let want = (start[i][a] + velocity[i][a]).rem_euclid(4.0);
            let gap = (s.pos[i][a] - want + 2.0).rem_euclid(4.0) - 2.0;
            assert!(
                gap.abs() < 1e-10,
                "particle {i} axis {a}: {} vs {want}",
                s.pos[i][a]
            );
        }
    }
    // The pressure of this ideal gas is rho T_kin.
    let p = s.pressure();
    assert!((p - s.density() * s.temperature()).abs() < 1e-12);
}

#[test]
fn momentum_is_conserved_and_positions_stay_in_the_box() {
    for dim in [2_usize, 3] {
        let box_len = if dim == 3 { 4.0 } else { 6.0 };
        let n = if dim == 3 { 192 } else { 108 };
        let mut s = Dpd::new(n, params(dim, box_len, 25.0, 4.5, 1.0), 21 + dim as u64).unwrap();
        s.run(0.03, 120);
        let m = s.momentum();
        assert!(
            m.iter().all(|x| x.abs() < 1e-9),
            "dim {dim}: momentum {m:?}"
        );
        assert!(
            s.pos
                .iter()
                .all(|r| (0..dim).all(|a| (0.0..box_len).contains(&r[a])))
        );
        assert!(
            s.pos
                .iter()
                .chain(&s.vel)
                .all(|v| v.iter().all(|x| x.is_finite()))
        );
        if dim == 2 {
            assert!(s.pos.iter().chain(&s.vel).all(|v| v[2] == 0.0));
        }
    }
}

#[test]
fn the_thermostat_relaxes_a_hot_start_and_equipartition_holds() {
    // Velocities scaled to four times the target temperature relax to kT by
    // the pair friction and noise (fluctuation-dissipation, sigma^2 = 2 gamma kT).
    let kt = 1.0;
    let mut s = Dpd::new(750, params(3, 6.3, 25.0, 4.5, kt), 31).unwrap();
    for v in &mut s.vel {
        for x in v.iter_mut() {
            *x *= 2.0;
        }
    }
    let hot = s.temperature();
    assert!((hot - 4.0 * kt).abs() < 0.3, "{hot}");
    // A moderate step: the modified velocity-Verlet temperature error grows
    // with dt (a few per cent at dt = 0.03), and the time average of a
    // finite system fluctuates by about 1 %, so the check is at 5 %.
    s.run(0.02, 700);
    let mid = s.temperature();
    assert!(mid < 0.6 * hot, "temperature decreases: {hot} -> {mid}");
    s.run(0.02, 800);
    // Time-average over a stretch at equilibrium.
    let mut samples = Vec::new();
    for _ in 0..100 {
        s.run(0.02, 10);
        samples.push(s.temperature());
    }
    let mean = samples.iter().sum::<f64>() / samples.len() as f64;
    assert!((mean - kt).abs() < 0.05, "equilibrium temperature {mean}");
    // Each Cartesian component carries kT/2 and the velocity distribution is
    // Gaussian (kurtosis 3).
    let n = s.len() as f64;
    for a in 0..3 {
        let m2 = s.vel.iter().map(|v| v[a] * v[a]).sum::<f64>() / n;
        assert!((m2 - kt).abs() < 0.25, "axis {a}: <v^2> = {m2}");
    }
    let m2 = s.vel.iter().flatten().map(|x| x * x).sum::<f64>() / (3.0 * n);
    let m4 = s.vel.iter().flatten().map(|x| x.powi(4)).sum::<f64>() / (3.0 * n);
    assert!(
        (m4 / (m2 * m2) - 3.0).abs() < 0.5,
        "kurtosis {}",
        m4 / (m2 * m2)
    );
}

#[test]
fn the_particle_mass_enters_the_temperature_and_the_dynamics() {
    let heavy = DpdParams {
        mass: 2.5,
        ..params(3, 4.0, 25.0, 4.5, 1.0)
    };
    let mut s = Dpd::new(200, heavy, 41).unwrap();
    assert!((s.temperature() - 1.0).abs() < 0.2, "{}", s.temperature());
    s.run(0.03, 400);
    assert!((s.temperature() - 1.0).abs() < 0.2, "{}", s.temperature());
    // Momentum uses the mass: sum m v = 0.
    assert!(s.momentum().iter().all(|x| x.abs() < 1e-9));
}

#[test]
fn runs_are_deterministic_in_the_seed() {
    let p = params(3, 4.0, 25.0, 4.5, 1.0);
    let run = |seed: u64| {
        let mut s = Dpd::new(100, p, seed).unwrap();
        s.run(0.03, 30);
        (s.pos, s.vel)
    };
    assert_eq!(run(5), run(5));
    assert_ne!(run(5).0, run(6).0);
}

mod walls {
    use rssn::sim::dpd::DpdWall;
    use rssn::sim::dpd::WallKind;

    use super::*;

    fn fluid(
        wall: Option<DpdWall>,
        seed: u64,
    ) -> Dpd {
        let mut s = Dpd::new(375, params(3, 5.0, 25.0, 4.5, 1.0), seed).unwrap();
        if let Some(w) = wall {
            s.set_walls(w).unwrap();
        }
        s
    }

    #[test]
    fn particles_stay_between_the_walls_and_tangential_momentum_survives_specular_walls() {
        let wall = DpdWall {
            axis: 2,
            a_wall: 10.0,
            kind: WallKind::Specular,
        };
        let mut s = fluid(Some(wall), 7);
        for v in &mut s.vel {
            v[0] += 0.8;
        }
        let m0 = s.momentum();
        s.run(0.02, 1000);
        assert!(
            s.pos.iter().all(|p| (0.0..=5.0).contains(&p[2])),
            "left the slab"
        );
        let m1 = s.momentum();
        // the wall force and the reflection act along z only: x and y momentum are conserved
        assert!(
            (m1[0] - m0[0]).abs() < 1e-8 * 375.0,
            "{} vs {}",
            m1[0],
            m0[0]
        );
        assert!((m1[1] - m0[1]).abs() < 1e-8 * 375.0);
        // temperature stays near the thermostat value plus the kinetic energy of the drift 0.8
        let expected = 1.0 + 0.8 * 0.8 / 3.0;
        assert!(
            (s.temperature() - expected).abs() < 0.1,
            "T = {}",
            s.temperature()
        );
    }

    #[test]
    fn bounce_back_walls_remove_the_flow_and_specular_walls_keep_it() {
        let run = |kind: WallKind| {
            let mut s = fluid(
                Some(DpdWall {
                    axis: 2,
                    a_wall: 0.0,
                    kind,
                }),
                11,
            );
            for v in &mut s.vel {
                v[0] += 0.8;
            }
            let start = s.mean_velocity(0);
            s.run(0.02, 1500);
            (start, s.mean_velocity(0))
        };
        let (start, specular) = run(WallKind::Specular);
        let (_, bounce) = run(WallKind::BounceBack);
        println!("mean vx: start {start:.3}, specular {specular:.3}, bounce-back {bounce:.3}");
        assert!((specular - start).abs() < 1e-8 * 10.0);
        assert!(
            bounce.abs() < 0.6 * start,
            "bounce-back kept {bounce} of {start}"
        );
    }

    #[test]
    fn a_wall_leaves_the_bulk_density_and_depletes_next_to_it() {
        let wall = DpdWall {
            axis: 2,
            a_wall: 25.0,
            kind: WallKind::Specular,
        };
        let mut s = fluid(Some(wall), 3);
        s.run(0.02, 800);
        let (mut mid, mut edge, mut samples) = (0.0, 0.0, 0.0);
        for _ in 0..600 {
            s.step(0.02);
            mid += s.pos.iter().filter(|p| (1.5..3.5).contains(&p[2])).count() as f64;
            edge += s.pos.iter().filter(|p| p[2] < 0.25).count() as f64;
            samples += 1.0;
        }
        // densities per unit slab thickness
        let (rho_mid, rho_edge) = (mid / samples / 2.0, edge / samples / 0.25);
        println!("slab density: middle {rho_mid:.1}, wall layer {rho_edge:.1}");
        assert!(
            (rho_mid / (375.0 / 5.0) - 1.0).abs() < 0.15,
            "bulk {rho_mid}"
        );
        assert!(
            rho_edge < 0.8 * rho_mid,
            "no depletion at the wall: {rho_edge} vs {rho_mid}"
        );
    }

    #[test]
    fn invalid_walls_are_rejected_and_no_wall_is_unchanged() {
        let mut s = fluid(None, 2);
        assert!(
            s.set_walls(DpdWall {
                axis: 3,
                a_wall: 1.0,
                kind: WallKind::Specular
            })
            .is_err()
        );
        assert!(
            s.set_walls(DpdWall {
                axis: 2,
                a_wall: -1.0,
                kind: WallKind::Specular
            })
            .is_err()
        );
        let (mut a, mut b) = (fluid(None, 5), fluid(None, 5));
        a.run(0.02, 20);
        b.run(0.02, 20);
        assert_eq!(a.pos, b.pos);
    }
}

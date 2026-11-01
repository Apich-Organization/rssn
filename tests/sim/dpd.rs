//! `sim::dpd`: the Groot–Warren DPD fluid. Direct checks of the pair force
//! (conservative repulsion, minimum image, virial), momentum conservation,
//! the free-particle and ideal-gas limits, thermostat relaxation and
//! equipartition, determinism and input validation, and of the walled slab:
//! confinement, momentum and energy balance at the reflection, the wall
//! contribution to the virial pressure, and the friction of a no-slip wall.

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

    /// 375 particles at density 3 in a slab `0 <= z <= 5`, optionally walled along `z`.
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

    fn specular(a_wall: f64) -> DpdWall {
        DpdWall {
            axis: 2,
            a_wall,
            kind: WallKind::Specular,
        }
    }

    fn kinetic(s: &Dpd) -> f64 {
        let sum = s
            .vel
            .iter()
            .flat_map(|v| v.iter())
            .map(|x| x * x)
            .sum::<f64>();
        0.5 * s.params.mass * sum
    }

    /// Separations of every pair, with the minimum image on the periodic axes and the plain
    /// distance on the walled axis 2, exactly as the force loop uses them.
    fn separations(
        s: &Dpd,
        walled: bool,
    ) -> Vec<f64> {
        let p = s.params;
        let mut out = Vec::new();
        for i in 0..s.len() {
            for j in (i + 1)..s.len() {
                let mut r2 = 0.0;
                for a in 0..p.dim {
                    let mut x = s.pos[i][a] - s.pos[j][a];
                    if !(walled && a == 2) {
                        x -= p.box_len[a] * (x / p.box_len[a]).round();
                    }
                    r2 += x * x;
                }
                out.push(r2.sqrt());
            }
        }
        out
    }

    /// Conservative pair potential `Σ a r_c (1 - r/r_c)² / 2` and the virial `Σ a (1 - r/r_c) r`.
    fn pair_terms(
        s: &Dpd,
        walled: bool,
    ) -> (f64, f64) {
        let p = s.params;
        let (mut energy, mut virial) = (0.0, 0.0);
        for r in separations(s, walled).iter().copied().filter(|r| *r < p.rc) {
            let w = 1.0 - r / p.rc;
            energy += 0.5 * p.a * p.rc * w * w;
            virial += p.a * w * r;
        }
        (energy, virial)
    }

    /// Wall potential `a_wall (r_c - d)² / (2 r_c)` and virial `f_wall d` summed over both
    /// planes for every particle within `r_c` of one.
    fn wall_terms(
        s: &Dpd,
        a_wall: f64,
    ) -> (f64, f64) {
        let p = s.params;
        let len = p.box_len[2];
        let (mut energy, mut virial) = (0.0, 0.0);
        for r in &s.pos {
            for d in [r[2], len - r[2]] {
                if d < p.rc {
                    let f = a_wall * (1.0 - d / p.rc);
                    energy += 0.5 * f * (p.rc - d);
                    virial += f * d;
                }
            }
        }
        (energy, virial)
    }

    /// `m Σ (v - V)² / (d (N - 1))`: the temperature with the streaming velocity removed, which
    /// is what the DPD thermostat targets and what a drift leaves untouched.
    fn peculiar_temperature(s: &Dpd) -> f64 {
        let n = s.len() as f64;
        let mut mean = [0.0; 3];
        for v in &s.vel {
            for a in 0..3 {
                mean[a] += v[a] / n;
            }
        }
        let mut sum = 0.0;
        for v in &s.vel {
            for a in 0..3 {
                let d = v[a] - mean[a];
                sum += d * d;
            }
        }
        s.params.mass * sum / (s.params.dim as f64 * (n - 1.0))
    }

    fn total_energy(
        s: &Dpd,
        a_wall: f64,
    ) -> f64 {
        kinetic(s) + pair_terms(s, true).0 + wall_terms(s, a_wall).0
    }

    #[test]
    fn particles_stay_between_the_walls_and_tangential_momentum_survives_specular_walls() {
        let mut s = fluid(Some(specular(10.0)), 7);
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
        // The temperature is the thermostat value plus the kinetic energy of the imposed drift.
        // A single instantaneous value is a poor estimator: `m Σ v² / (d (N - 1))` has a relative
        // standard deviation of sqrt(2 / (d N)) = 4.2 % here (plus the 0.8² drift), so the old
        // `±0.1` on one sample was a two-sigma test. The measured autocorrelation of the
        // peculiar temperature is 0.10 after t = 2 and 0.00 after t = 4, so averages over blocks
        // of t = 2 are independent to within 10 %; ten of them give a standard error of 1.3 %,
        // and the bound is three of those.
        s.run(0.02, 1000);
        let mut blocks = Vec::new();
        for _ in 0..10 {
            s.run(0.02, 100);
            blocks.push(s.temperature());
        }
        let n = blocks.len() as f64;
        let mean = blocks.iter().sum::<f64>() / n;
        let var = blocks.iter().map(|x| (x - mean) * (x - mean)).sum::<f64>() / (n - 1.0);
        let se = (var / n).sqrt();
        let expected = 1.0 + 0.8 * 0.8 / 3.0;
        println!("mean T = {mean:.4} +- {se:.4} over {n} blocks, want {expected:.4}");
        // + 0.002 covers the N / (N - 1) of the estimator, which puts the exact expectation
        // 0.0006 above 1 + 0.8²/3.
        assert!(
            (mean - expected).abs() < 3.0 * se + 0.002,
            "mean T = {mean} +- {se} over {n} blocks, want {expected}"
        );
    }

    #[test]
    fn a_specular_wall_neither_injects_nor_removes_energy() {
        // With the thermostat switched off the wall force is conservative and the reflection is
        // elastic, so KE + U_pair + U_wall must be constant. Applying the reflection to the whole
        // post-step velocity, and so mirroring the post-collision force impulse as well, added
        // 0.075 per crossing (+4.6 % of the total energy over t = 20); mirroring only the
        // pre-collision half of the step leaves 4e-4 to 1.8e-3 over eight seeds, which is the
        // same order as the periodic integrator error (5e-4 at this step).
        let wall = specular(10.0);
        let mut s = fluid(Some(wall), 5);
        s.run(0.02, 2000);
        s.params.gamma = 0.0;
        s.params.kt = 0.0;
        s.set_walls(wall).unwrap();
        let e0 = total_energy(&s, 10.0);
        s.run(0.01, 2000);
        let e1 = total_energy(&s, 10.0);
        assert!(
            (e1 - e0).abs() < 5e-3 * e0,
            "energy {} -> {} ({:e} relative)",
            e0,
            e1,
            (e1 - e0) / e0
        );
    }

    #[test]
    fn the_wall_force_enters_the_virial_pressure() {
        // p = rho T + (pair virial + wall virial) / (d V), with both virials recomputed here from
        // the positions of the last force evaluation.
        let mut s = fluid(Some(specular(25.0)), 3);
        s.run(0.02, 600);
        let (_, pair) = pair_terms(&s, true);
        let (_, wall) = wall_terms(&s, 25.0);
        assert!(
            wall > 0.05 * pair,
            "the wall virial must matter here: {wall} vs {pair}"
        );
        let volume = 125.0;
        let want = s.density() * s.temperature() + (pair + wall) / (3.0 * volume);
        assert!(
            (s.pressure() - want).abs() < 1e-9 * want,
            "{} vs {want}",
            s.pressure()
        );
    }

    #[test]
    fn a_wall_out_of_range_leaves_the_pressure_alone() {
        // Every particle sits at least r_c away from both planes, so neither the wall force nor
        // its virial is non-zero and the pressure is the periodic one: adding the wall cannot
        // make it jump.
        let mut s = Dpd::new(64, params(3, 8.0, 25.0, 4.5, 1.0), 3).unwrap();
        for (i, p) in s.pos.iter_mut().enumerate() {
            for (a, x) in p.iter_mut().enumerate() {
                *x = 3.2 + 0.4 * ((i >> (2 * a)) & 3) as f64;
            }
        }
        s.set_walls(specular(0.0)).unwrap();
        let cold = s.pressure();
        s.set_walls(specular(25.0)).unwrap();
        assert_eq!(s.pressure(), cold);
        assert!(
            s.pos.iter().all(|p| (1.0..7.0).contains(&p[2])),
            "the cluster must stay clear of the walls"
        );
    }

    #[test]
    fn a_frictional_wall_decays_a_drift_and_a_specular_wall_keeps_it() {
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
            let drift = s.mean_velocity(0);
            // time average of the peculiar temperature over the last t = 20 in blocks of t = 4
            let mut t = 0.0;
            for _ in 0..5 {
                s.run(0.02, 200);
                t += peculiar_temperature(&s);
            }
            (start, drift, t / 5.0)
        };
        let (start, specular_drift, specular_t) = run(WallKind::Specular);
        let (_, frictional_drift, frictional_t) = run(WallKind::NoSlip { friction: 1.0 });
        assert!((specular_drift - start).abs() < 1e-8 * 10.0);
        assert!(
            frictional_drift.abs() < 0.1 * start,
            "frictional wall kept {frictional_drift} of {start}"
        );
        // The Langevin collision removes the drift without cooling the fluid; a purely
        // dissipative rule leaves it at about 0.85 kT.
        assert!(
            (frictional_t - 1.0).abs() < 0.1,
            "frictional wall temperature {frictional_t}"
        );
        assert!(
            (specular_t - 1.0).abs() < 0.15,
            "specular wall temperature {specular_t}"
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
        let mut s = fluid(Some(specular(25.0)), 3);
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
        for friction in [-0.1, 1.5, f64::NAN] {
            assert!(
                s.set_walls(DpdWall {
                    axis: 2,
                    a_wall: 1.0,
                    kind: WallKind::NoSlip { friction }
                })
                .is_err(),
                "friction {friction}"
            );
        }
        let (mut a, mut b) = (fluid(None, 5), fluid(None, 5));
        a.run(0.02, 20);
        b.run(0.02, 20);
        assert_eq!(a.pos, b.pos);
    }
}

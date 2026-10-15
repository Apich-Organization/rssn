//! Dimension-generic MD engine: neighbour lists, boundaries, thermostats, barostats, constraints, observables.

use rssn::kernels::random::Rng;
use rssn::sim::md_engine::*;

fn lj() -> Vec<PairTerm> {
    vec![PairTerm::LennardJones { eps: 1.0, sigma: 1.0 }]
}

fn lj_fluid_3d(
    cells: usize,
    rho: f64,
    t: f64,
    dt: f64,
    seed: u64,
) -> MdSimulation<3> {
    let n = 4 * cells * cells * cells;
    let l = (n as f64 / rho).cbrt();
    let a = l / cells as f64;
    let bx = SimBox::periodic([l; 3]);
    let mut sys = MdSystem::new(lattice_fcc(cells, a), bx);
    let mut rng = Rng::new(seed);
    sys.maxwell_boltzmann(t, &mut rng);
    let mut ff = ForceField::new(lj(), 2.5);
    ff.shift_force = true;
    MdSimulation::new(sys, ff, 0.3, dt, seed + 1).unwrap()
}

/// Brute-force reference of pair energy and forces.
fn brute<const D: usize>(
    sys: &MdSystem<D>,
    rc: f64,
    shift: bool,
) -> (f64, Vec<[f64; D]>) {
    let term = PairTerm::LennardJones { eps: 1.0, sigma: 1.0 };
    let us = if shift { term.eval(rc * rc, 0.0).0 } else { 0.0 };
    let mut e = 0.0;
    let mut f = vec![[0.0; D]; sys.n()];
    for i in 0..sys.n() {
        for j in i + 1..sys.n() {
            let mut d = [0.0; D];
            for a in 0..D {
                d[a] = sys.pos[i][a] - sys.pos[j][a];
            }
            sys.bx.min_image(&mut d);
            let r2: f64 = d.iter().map(|x| x * x).sum();
            if r2 < rc * rc {
                let (u, fr) = term.eval(r2, 0.0);
                e += u - us;
                for a in 0..D {
                    f[i][a] += fr * d[a];
                    f[j][a] -= fr * d[a];
                }
            }
        }
    }
    (e, f)
}

fn check_cell_list<const D: usize>(
    bx: SimBox<D>,
    n: usize,
    dmin: f64,
    seed: u64,
) {
    let mut rng = Rng::new(seed);
    let pos = random_positions(n, &bx, dmin, &mut rng).unwrap();
    let sys = MdSystem::new(pos, bx);
    let mut ff = ForceField::<D>::new(lj(), 2.5);
    ff.shift = true;
    let mut nl = NeighbourList::new(2.5, 0.3);
    let out = compute_forces(&sys, &ff, &mut nl).unwrap();
    let (e, f) = brute(&sys, 2.5, true);
    assert!((out.epot - e).abs() < 1e-9 * e.abs().max(1.0), "D={D}: {} vs {e}", out.epot);
    for i in 0..n {
        for a in 0..D {
            assert!((out.f[i][a] - f[i][a]).abs() < 1e-8 * f[i][a].abs().max(1.0), "D={D} i={i}");
        }
    }
}

#[test]
fn cell_list_matches_brute_force_in_one_two_and_three_dimensions() {
    check_cell_list::<1>(SimBox::periodic([30.0]), 25, 0.9, 1);
    check_cell_list::<2>(SimBox::periodic([14.0, 15.0]), 120, 0.85, 2);
    check_cell_list::<3>(SimBox::periodic([9.0, 9.5, 10.0]), 300, 0.85, 3);
    // mixed boundaries: periodic, reflecting, open
    check_cell_list::<3>(SimBox::new([9.0, 9.5, 10.0], [Boundary::Periodic, Boundary::Reflecting, Boundary::Open]), 250, 0.85, 4);
    // very small box (fewer than three cells per axis)
    check_cell_list::<2>(SimBox::periodic([5.0, 5.2]), 14, 0.9, 5);
}

#[test]
fn verlet_skin_avoids_rebuilds_and_stays_correct() {
    let mut sim = lj_fluid_3d(3, 0.8, 1.0, 0.004, 7);
    sim.run(200).unwrap();
    assert!(sim.nl.builds < 200 / 4, "rebuilds {}", sim.nl.builds);
    let mut fresh = NeighbourList::new(2.5, 0.0);
    let reference = compute_forces(&sim.sys, &sim.ff, &mut fresh).unwrap();
    assert!((sim.forces.epot - reference.epot).abs() < 1e-8 * reference.epot.abs());
    for i in 0..sim.sys.n() {
        for a in 0..3 {
            assert!((sim.forces.f[i][a] - reference.f[i][a]).abs() < 1e-7 * reference.f[i][a].abs().max(1.0));
        }
    }
}

#[test]
fn pair_terms_have_consistent_energy_and_force() {
    let terms = [
        PairTerm::LennardJones { eps: 0.7, sigma: 1.1 },
        PairTerm::Morse { d: 1.3, alpha: 1.7, r0: 1.2 },
        PairTerm::Coulomb { k: 1.4, alpha: 0.0 },
        PairTerm::Coulomb { k: 1.4, alpha: 0.9 },
    ];
    for t in terms {
        for r in [0.9, 1.3, 2.1] {
            let h = 1e-6;
            let up = t.eval((r + h) * (r + h), -1.0).0;
            let um = t.eval((r - h) * (r - h), -1.0).0;
            let (_, fr) = t.eval(r * r, -1.0);
            assert!((fr * r + (up - um) / (2.0 * h)).abs() < 1e-6, "{t:?} at {r}");
        }
    }
}

fn numeric_forces<const D: usize>(
    sys: &MdSystem<D>,
    ff: &ForceField<D>,
) -> Vec<[f64; D]> {
    let mut out = vec![[0.0; D]; sys.n()];
    let h = 1e-6;
    for i in 0..sys.n() {
        for a in 0..D {
            let mut p = sys.clone();
            p.pos[i][a] += h;
            let mut nl = NeighbourList::new(ff.rcut, 0.3);
            let ep = compute_forces(&p, ff, &mut nl).unwrap().epot;
            p.pos[i][a] -= 2.0 * h;
            let mut nl = NeighbourList::new(ff.rcut, 0.3);
            let em = compute_forces(&p, ff, &mut nl).unwrap().epot;
            out[i][a] = -(ep - em) / (2.0 * h);
        }
    }
    out
}

#[test]
fn bonded_forces_match_numerical_gradients() {
    // 3-D molecule with bonds, angles and a dihedral
    let pos = vec![[0.1, 0.0, 0.2], [1.2, 0.3, 0.0], [1.9, 1.4, 0.4], [3.0, 1.6, 1.3], [2.5, 2.9, 0.1]];
    let sys = MdSystem::new(pos, SimBox::new([20.0; 3], [Boundary::Open; 3]));
    let mut ff = ForceField::<3>::new(vec![], 1.0);
    ff.bonds = vec![Bond { i: 0, j: 1, k: 30.0, r0: 1.1 }, Bond { i: 1, j: 2, k: 25.0, r0: 1.3 }];
    ff.angles = vec![Angle { i: 0, j: 1, k: 2, ktheta: 12.0, theta0: 1.9 }, Angle { i: 2, j: 3, k: 4, ktheta: 9.0, theta0: 2.1 }];
    ff.dihedrals = vec![Dihedral { atoms: [0, 1, 2, 3], k: 2.0, n: 3.0, delta: 0.4 }, Dihedral { atoms: [1, 2, 3, 4], k: 1.5, n: 2.0, delta: 0.0 }];
    let mut nl = NeighbourList::new(1.0, 0.3);
    let out = compute_forces(&sys, &ff, &mut nl).unwrap();
    let num = numeric_forces(&sys, &ff);
    for i in 0..5 {
        for a in 0..3 {
            assert!((out.f[i][a] - num[i][a]).abs() < 1e-6, "atom {i} comp {a}: {} vs {}", out.f[i][a], num[i][a]);
        }
    }
    let net: f64 = (0..3).map(|a| out.f.iter().map(|f| f[a]).sum::<f64>().abs()).sum();
    assert!(net < 1e-10);
    // dihedrals need D = 3
    let mut ff2 = ForceField::<2>::new(vec![], 1.0);
    ff2.dihedrals = ff.dihedrals.clone();
    assert!(ff2.validate().is_err());
}

#[test]
fn angle_forces_work_in_two_dimensions() {
    let sys = MdSystem::new(vec![[0.0, 0.1], [1.0, 0.0], [1.4, 1.2]], SimBox::new([20.0; 2], [Boundary::Open; 2]));
    let mut ff = ForceField::<2>::new(vec![], 1.0);
    ff.angles = vec![Angle { i: 0, j: 1, k: 2, ktheta: 5.0, theta0: 1.8 }];
    let mut nl = NeighbourList::new(1.0, 0.3);
    let out = compute_forces(&sys, &ff, &mut nl).unwrap();
    let num = numeric_forces(&sys, &ff);
    for i in 0..3 {
        for a in 0..2 {
            assert!((out.f[i][a] - num[i][a]).abs() < 1e-6);
        }
    }
}

#[test]
fn virial_trace_equals_the_scaling_derivative_of_the_energy() {
    let pos = vec![[0.1, 0.0, 0.2], [1.2, 0.3, 0.0], [1.9, 1.4, 0.4], [3.0, 1.6, 1.3], [2.5, 2.9, 0.1], [1.0, 2.0, 2.5]];
    let sys = MdSystem::new(pos, SimBox::new([40.0; 3], [Boundary::Open; 3]));
    let mut ff = ForceField::<3>::new(lj(), 3.0);
    ff.bonds = vec![Bond { i: 0, j: 1, k: 30.0, r0: 1.1 }];
    ff.angles = vec![Angle { i: 0, j: 1, k: 2, ktheta: 12.0, theta0: 1.9 }];
    ff.dihedrals = vec![Dihedral { atoms: [1, 2, 3, 4], k: 1.5, n: 2.0, delta: 0.3 }];
    ff.exclude_bonded = true;
    let energy_at = |lam: f64| {
        let mut p = sys.clone();
        for q in &mut p.pos {
            for c in q.iter_mut() {
                *c *= lam;
            }
        }
        let mut f2 = ForceField::<3>::new(lj(), 3.0);
        f2.bonds = ff.bonds.clone();
        f2.angles = ff.angles.clone();
        f2.dihedrals = ff.dihedrals.clone();
        f2.rebuild_exclusions();
        let mut nl = NeighbourList::new(3.0, 0.3);
        compute_forces(&p, &f2, &mut nl).unwrap().epot
    };
    let h = 1e-6;
    let de = (energy_at(1.0 + h) - energy_at(1.0 - h)) / (2.0 * h);
    let mut f2 = ForceField::<3>::new(lj(), 3.0);
    f2.bonds = ff.bonds.clone();
    f2.angles = ff.angles.clone();
    f2.dihedrals = ff.dihedrals.clone();
    f2.rebuild_exclusions();
    let mut nl = NeighbourList::new(3.0, 0.3);
    let out = compute_forces(&sys, &f2, &mut nl).unwrap();
    let tr: f64 = (0..3).map(|a| out.virial[a][a]).sum();
    assert!((tr + de).abs() < 1e-5 * de.abs().max(1.0), "trace {tr} vs -dE/dlambda {}", -de);
}

#[test]
fn nve_energy_error_scales_with_the_square_of_the_time_step() {
    let err = |dt: f64| {
        let mut sim = lj_fluid_3d(3, 0.8, 1.0, dt, 11);
        sim.run((0.8 / dt).round() as usize).unwrap();
        let e0 = sim.total_energy();
        let steps = (0.4 / dt).round() as usize;
        let mut worst = 0.0_f64;
        for _ in 0..steps {
            sim.step().unwrap();
            worst = worst.max((sim.total_energy() - e0).abs());
        }
        worst
    };
    let e1 = err(0.008);
    let e2 = err(0.004);
    assert!(e1 < 0.05, "energy error {e1}");
    assert!(e1 / e2 > 2.5 && e1 / e2 < 6.0, "ratio {}", e1 / e2);
    // momentum stays zero
    let mut sim = lj_fluid_3d(3, 0.8, 1.0, 0.004, 12);
    sim.run(300).unwrap();
    for p in sim.sys.momentum() {
        assert!(p.abs() < 1e-10);
    }
}

#[test]
fn lj_liquid_rdf_has_its_first_peak_near_the_potential_minimum() {
    let mut sim = lj_fluid_3d(4, 0.8442, 1.0, 0.004, 21);
    sim.thermostat = Thermostat::Berendsen { t0: 1.0, tau: 0.2 };
    sim.run(1500).unwrap();
    sim.thermostat = Thermostat::None;
    let mut rdf = Rdf::<3>::new(80, sim.sys.bx.len[0] / 2.0);
    for k in 0..800 {
        sim.step().unwrap();
        if k % 20 == 0 {
            rdf.add_frame(&sim.sys);
        }
    }
    let g = rdf.result();
    let (rpk, gpk) = g.iter().copied().fold((0.0, 0.0), |b, x| if x.1 > b.1 { x } else { b });
    assert!(rpk > 1.0 && rpk < 1.18, "first peak at {rpk}");
    assert!(gpk > 2.2 && gpk < 3.6, "peak height {gpk}");
    assert!(g[0].1 < 1e-3, "core is excluded: {}", g[0].1);
    let tail: f64 = g.iter().filter(|x| x.0 > 2.8).map(|x| x.1).sum::<f64>() / g.iter().filter(|x| x.0 > 2.8).count() as f64;
    assert!((tail - 1.0).abs() < 0.08, "long-range limit {tail}");
}

#[test]
fn rdf_of_an_ideal_gas_is_one_in_every_dimension() {
    fn run<const D: usize>(bx: SimBox<D>, n: usize) {
        let mut rng = Rng::new(5 + D as u64);
        let mut rdf = Rdf::<D>::new(12, 0.5 * bx.len[0].min(bx.len[D - 1]));
        let mut sys = MdSystem::new(vec![[0.0; D]; n], bx);
        for _ in 0..60 {
            for p in &mut sys.pos {
                for a in 0..D {
                    p[a] = rng.uniform() * sys.bx.len[a];
                }
            }
            rdf.add_frame(&sys);
        }
        for (r, g) in rdf.result().iter().skip(2) {
            assert!((g - 1.0).abs() < 0.12, "D={D} r={r} g={g}");
        }
    }
    run::<1>(SimBox::periodic([100.0]), 400);
    run::<2>(SimBox::periodic([20.0, 20.0]), 400);
    run::<3>(SimBox::periodic([8.0, 8.0, 8.0]), 400);
}

#[test]
fn langevin_baoab_samples_the_canonical_ensemble() {
    // 64 independent 1-D harmonic oscillators, k = 1, T = 1.3
    let n = 64;
    let bx = SimBox::new([100.0], [Boundary::Open]);
    let mut sys = MdSystem::new((0..n).map(|i| [0.01 * i as f64]).collect(), bx);
    sys.vel = vec![[0.0]; n];
    let ff = {
        let mut ff = ForceField::<1>::new(vec![], 1.0);
        ff.extra.push(Box::new(|sys: &MdSystem<1>, out: &mut Forces<1>| {
            for i in 0..sys.n() {
                out.f[i][0] -= sys.pos[i][0];
                out.epot += 0.5 * sys.pos[i][0] * sys.pos[i][0];
            }
        }));
        ff
    };
    let mut sim = MdSimulation::new(sys, ff, 0.3, 0.1, 99).unwrap();
    sim.thermostat = Thermostat::Langevin { t0: 1.3, gamma: 1.0 };
    sim.run(500).unwrap();
    let (mut x2, mut v2, mut cnt) = (0.0, 0.0, 0.0);
    for _ in 0..4000 {
        sim.step().unwrap();
        for i in 0..n {
            x2 += sim.sys.pos[i][0].powi(2);
            v2 += sim.sys.vel[i][0].powi(2);
            cnt += 1.0;
        }
    }
    assert!((x2 / cnt - 1.3).abs() < 0.06, "<x^2> = {}", x2 / cnt);
    assert!((v2 / cnt - 1.3).abs() < 0.06, "<v^2> = {}", v2 / cnt);
}

#[test]
fn nose_hoover_chain_thermostats_a_fluid_and_conserves_its_extended_energy() {
    let mut sim = lj_fluid_3d(3, 0.8, 1.2, 0.004, 31);
    let nf = sim.dof();
    sim.thermostat = Thermostat::NoseHooverChain(NoseHooverChain::new(1.2, 0.5, 3, nf));
    sim.run(1500).unwrap();
    let h0 = sim.conserved_energy();
    let (mut t_sum, mut worst) = (0.0, 0.0_f64);
    let steps = 6000;
    for _ in 0..steps {
        sim.step().unwrap();
        t_sum += sim.temperature();
        worst = worst.max((sim.conserved_energy() - h0).abs());
    }
    assert!((t_sum / steps as f64 - 1.2).abs() < 0.04, "mean T {}", t_sum / steps as f64);
    assert!(worst < 3e-3 * sim.sys.kinetic_energy(), "extended energy drift {worst}");
}

#[test]
fn berendsen_thermostat_relaxes_to_the_target() {
    let mut sim = lj_fluid_3d(3, 0.8, 0.5, 0.004, 41);
    sim.thermostat = Thermostat::Berendsen { t0: 1.5, tau: 0.1 };
    sim.run(1500).unwrap();
    assert!((sim.temperature() - 1.5).abs() < 0.1, "T {}", sim.temperature());
}

#[test]
fn berendsen_barostat_drives_the_pressure_to_its_target() {
    let mut sim = lj_fluid_3d(3, 0.6, 1.5, 0.004, 51);
    sim.thermostat = Thermostat::Berendsen { t0: 1.5, tau: 0.1 };
    sim.barostat = Barostat::Berendsen { p0: 2.0, tau: 1.0, beta: 0.05 };
    let v0 = sim.sys.bx.volume();
    sim.run(3000).unwrap();
    let mut p = 0.0;
    let n = 2000;
    for _ in 0..n {
        sim.step().unwrap();
        p += sim.pressure();
    }
    assert!((p / n as f64 - 2.0).abs() < 0.25, "mean pressure {}", p / n as f64);
    assert!((sim.sys.bx.volume() - v0).abs() > 0.5, "the box must have changed");
}

#[test]
fn mtk_barostat_conserves_the_extended_enthalpy_with_second_order_error() {
    let err = |dt: f64| {
        let mut sim = lj_fluid_3d(3, 0.75, 1.2, dt, 61);
        sim.barostat = Barostat::Mtk { p0: 1.0, w: 50.0, p_eps: 0.0 };
        sim.run((0.4 / dt) as usize).unwrap();
        let h0 = sim.conserved_energy();
        let mut worst = 0.0_f64;
        for _ in 0..(0.8 / dt) as usize {
            sim.step().unwrap();
            worst = worst.max((sim.conserved_energy() - h0).abs());
        }
        (worst, (sim.sys.bx.volume() - 27.0 * 4.0 / 0.75).abs())
    };
    let (e1, dv) = err(0.008);
    let (e2, _) = err(0.004);
    assert!(e1 < 6e-3 * 160.0, "enthalpy error {e1}");
    assert!(e1 / e2 > 2.0, "ratio {}", e1 / e2);
    assert!(dv > 1e-3, "volume must fluctuate");
    // combining the MTK barostat with a thermostat is rejected
    let mut sim = lj_fluid_3d(3, 0.75, 1.2, 0.004, 62);
    sim.barostat = Barostat::Mtk { p0: 1.0, w: 50.0, p_eps: 0.0 };
    sim.thermostat = Thermostat::Berendsen { t0: 1.0, tau: 1.0 };
    assert!(sim.step().is_err());
}

#[test]
fn shake_rattle_keep_rigid_dimers_rigid_and_conserve_energy() {
    let nmol = 32;
    let l = 7.0;
    let bx = SimBox::periodic([l; 3]);
    let mut rng = Rng::new(71);
    let centres = random_positions(nmol, &bx, 1.6, &mut rng).unwrap();
    let d = 0.7;
    let mut pos = Vec::new();
    for c in &centres {
        let mut u = [rng.normal(), rng.normal(), rng.normal()];
        let nrm = (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]).sqrt();
        for x in &mut u {
            *x /= nrm;
        }
        pos.push([c[0] + 0.5 * d * u[0], c[1] + 0.5 * d * u[1], c[2] + 0.5 * d * u[2]]);
        pos.push([c[0] - 0.5 * d * u[0], c[1] - 0.5 * d * u[1], c[2] - 0.5 * d * u[2]]);
    }
    let mut sys = MdSystem::new(pos, bx);
    sys.maxwell_boltzmann(1.0, &mut rng);
    let mut ff = ForceField::<3>::new(lj(), 2.5);
    ff.shift = true;
    // zero-strength bonds exclude the intramolecular pair from the non-bonded sum
    ff.bonds = (0..nmol).map(|m| Bond { i: 2 * m, j: 2 * m + 1, k: 0.0, r0: d }).collect();
    let mut sim = MdSimulation::new(sys, ff, 0.3, 0.004, 72).unwrap();
    sim.constraints = (0..nmol).map(|m| Constraint { i: 2 * m, j: 2 * m + 1, d }).collect();
    // start from velocities consistent with the constraints
    for _ in 0..3 {
        sim.step().unwrap();
    }
    assert_eq!(sim.dof(), (3 * 2 * nmol - 3 - nmol) as f64);
    let e0 = sim.total_energy();
    let mut worst_len = 0.0_f64;
    let mut worst_e = 0.0_f64;
    for _ in 0..1500 {
        sim.step().unwrap();
        for c in &sim.constraints {
            let mut dv = [0.0; 3];
            for a in 0..3 {
                dv[a] = sim.sys.pos[c.i][a] - sim.sys.pos[c.j][a];
            }
            sim.sys.bx.min_image(&mut dv);
            worst_len = worst_len.max((dv.iter().map(|x| x * x).sum::<f64>().sqrt() - d).abs());
            let mut vr = 0.0;
            for a in 0..3 {
                vr += dv[a] * (sim.sys.vel[c.i][a] - sim.sys.vel[c.j][a]);
            }
            assert!(vr.abs() < 1e-9, "relative velocity along the bond {vr}");
        }
        worst_e = worst_e.max((sim.total_energy() - e0).abs());
    }
    assert!(worst_len < 1e-9, "bond length drift {worst_len}");
    assert!(worst_e < 0.05 * sim.sys.kinetic_energy(), "energy drift {worst_e}");
}

#[test]
fn reflecting_walls_conserve_energy_and_bounce_in_every_dimension() {
    // 1-D particle bouncing between walls: triangular-wave trajectory
    let bx = SimBox::new([1.0], [Boundary::Reflecting]);
    let mut sys = MdSystem::new(vec![[0.2]], bx);
    sys.vel[0] = [0.7];
    let ff = ForceField::<1>::new(vec![], 1.0);
    let mut sim = MdSimulation::new(sys, ff, 0.3, 0.01, 1).unwrap();
    sim.run(1000).unwrap();
    // unfolded position 0.2 + 0.7*10 = 7.2 -> reflected into [0,1]: 7.2 mod 2 = 1.2 -> 0.8
    assert!((sim.sys.pos[0][0] - 0.8).abs() < 1e-9, "{}", sim.sys.pos[0][0]);
    assert!((sim.sys.kinetic_energy() - 0.245).abs() < 1e-12);
    // 3-D gas with LJ interaction in a reflecting box keeps energy and stays inside
    let bx = SimBox::new([6.0; 3], [Boundary::Reflecting; 3]);
    let mut rng = Rng::new(3);
    let pos = random_positions(60, &bx, 1.0, &mut rng).unwrap();
    let mut sys = MdSystem::new(pos, bx);
    sys.maxwell_boltzmann(1.0, &mut rng);
    let mut ff = ForceField::<3>::new(lj(), 2.5);
    ff.shift = true;
    let mut sim = MdSimulation::new(sys, ff, 0.3, 0.002, 4).unwrap();
    let e0 = sim.total_energy();
    sim.run(2000).unwrap();
    assert!((sim.total_energy() - e0).abs() < 0.02 * e0.abs().max(10.0));
    for p in &sim.sys.pos {
        assert!(p.iter().all(|&c| (0.0..=6.0).contains(&c)));
    }
}

#[test]
fn thermal_walls_impose_their_temperature_on_a_gas() {
    let bx = SimBox::new([4.0, 4.0], [Boundary::Thermal(2.0), Boundary::Thermal(2.0)]);
    let mut rng = Rng::new(8);
    let pos = random_positions(300, &bx, 0.0, &mut rng).unwrap();
    let mut sys = MdSystem::new(pos, bx);
    sys.maxwell_boltzmann(0.2, &mut rng);
    let ff = ForceField::<2>::new(vec![], 1.0);
    let mut sim = MdSimulation::new(sys, ff, 0.3, 0.02, 9).unwrap();
    sim.run(3000).unwrap();
    let (mut ke, mut cnt) = (0.0, 0.0);
    for _ in 0..3000 {
        sim.step().unwrap();
        ke += sim.sys.kinetic_energy();
        cnt += 1.0;
    }
    // equipartition: K = N D T / 2
    let t = 2.0 * ke / cnt / (300.0 * 2.0);
    assert!((t - 2.0).abs() < 0.12, "gas temperature {t}");
}

#[test]
fn lees_edwards_minimum_image_and_shear_stress() {
    let mut bx = SimBox::<2>::periodic([1.0, 1.0]);
    bx.shear_rate = 0.5;
    bx.offset = 0.3;
    let mut d = [0.1, 0.9];
    bx.min_image(&mut d);
    assert!((d[0] + 0.2).abs() < 1e-12 && (d[1] + 0.1).abs() < 1e-12, "{d:?}");
    // a particle crossing the top boundary is shifted by -offset and loses gamma*L of velocity
    let mut sys = MdSystem::new(vec![[0.5, 0.99]], bx);
    sys.vel[0] = [0.0, 1.0];
    let ff = ForceField::<2>::new(vec![], 1.0);
    let mut sim = MdSimulation::new(sys, ff, 0.3, 0.02, 1).unwrap();
    sim.step().unwrap();
    assert!(sim.sys.pos[0][1] < 0.05);
    assert!((sim.sys.vel[0][0] + 0.5).abs() < 1e-12);
    // sheared LJ fluid: negative P_xy for positive shear rate
    let cells = 3;
    let n = 4 * cells * cells * cells;
    let l = (n as f64 / 0.8442).cbrt();
    let mut bx = SimBox::<3>::periodic([l; 3]);
    bx.shear_rate = 0.5;
    let mut sys = MdSystem::new(lattice_fcc(cells, l / cells as f64), bx);
    let mut rng = Rng::new(3);
    sys.maxwell_boltzmann(1.0, &mut rng);
    for i in 0..sys.n() {
        sys.vel[i][0] += 0.5 * (sys.pos[i][1] - 0.5 * l);
    }
    let mut ff = ForceField::<3>::new(lj(), 2.5);
    ff.shift = true;
    let mut sim = MdSimulation::new(sys, ff, 0.3, 0.003, 4).unwrap();
    sim.run(600).unwrap();
    let mut pxy = 0.0;
    for _ in 0..600 {
        sim.step().unwrap();
        pxy += sim.pressure_tensor()[0][1];
    }
    assert!(pxy / 600.0 < -0.1, "mean P_xy = {}", pxy / 600.0);
}

#[test]
fn ideal_gas_pressure_msd_and_initial_conditions() {
    fn msd_check<const D: usize>() {
        let bx = SimBox::<D>::periodic([5.0; D]);
        let mut rng = Rng::new(2);
        let pos = random_positions(80, &bx, 0.0, &mut rng).unwrap();
        let mut sys = MdSystem::new(pos, bx);
        sys.maxwell_boltzmann(1.3, &mut rng);
        let reference = sys.unwrapped.clone();
        let v2: f64 = sys.vel.iter().map(|v| v.iter().map(|x| x * x).sum::<f64>()).sum::<f64>() / 80.0;
        let ff = ForceField::<D>::new(vec![], 1.0);
        let mut sim = MdSimulation::new(sys, ff, 0.3, 0.05, 1).unwrap();
        sim.run(200).unwrap();
        // free flight: msd = <v^2> t^2, far beyond the box size thanks to unwrapping
        let t = 10.0;
        assert!((msd(&sim.sys, &reference) - v2 * t * t).abs() < 1e-8 * v2 * t * t, "D={D}");
        // ideal-gas pressure: P V = (N-1) T for exact temperature and zero momentum
        let p = sim.pressure();
        assert!((p * sim.sys.bx.volume() - 79.0 * 1.3).abs() < 1e-8, "D={D}: {p}");
        for q in sim.sys.momentum() {
            assert!(q.abs() < 1e-10);
        }
        assert!((sim.temperature() - 1.3).abs() < 1e-10);
    }
    msd_check::<1>();
    msd_check::<2>();
    msd_check::<3>();
    // lattices
    let fcc = lattice_fcc(3, 1.6);
    assert_eq!(fcc.len(), 108);
    let bx = SimBox::<3>::periodic([4.8; 3]);
    let mut dmin = f64::INFINITY;
    for i in 0..fcc.len() {
        for j in i + 1..fcc.len() {
            let mut d = [fcc[i][0] - fcc[j][0], fcc[i][1] - fcc[j][1], fcc[i][2] - fcc[j][2]];
            bx.min_image(&mut d);
            dmin = dmin.min(d.iter().map(|x| x * x).sum::<f64>().sqrt());
        }
    }
    assert!((dmin - 1.6 / 2.0_f64.sqrt()).abs() < 1e-12);
    let tri = lattice_triangular(6, 6, 1.0);
    let bx = SimBox::<2>::periodic([6.0, 6.0 * 3.0_f64.sqrt() / 2.0]);
    let mut dmin = f64::INFINITY;
    for i in 0..tri.len() {
        for j in i + 1..tri.len() {
            let mut d = [tri[i][0] - tri[j][0], tri[i][1] - tri[j][1]];
            bx.min_image(&mut d);
            dmin = dmin.min(d.iter().map(|x| x * x).sum::<f64>().sqrt());
        }
    }
    assert!((dmin - 1.0).abs() < 1e-12, "triangular nearest neighbour {dmin}");
    assert_eq!(lattice_cubic([3, 4], 1.0).len(), 12);
    let bx = SimBox::<3>::periodic([3.0; 3]);
    assert!(random_positions(500, &bx, 1.5, &mut Rng::new(1)).is_none());
}

#[test]
fn vacf_and_energy_drift_helpers() {
    // harmonic oscillator velocity: v(t) = cos(t)
    let frames: Vec<Vec<[f64; 1]>> = (0..4000).map(|k| vec![[(0.01 * k as f64 * 2.0).cos()]]).collect();
    let c = vacf(&frames, 400);
    assert!((c[0] - 0.5).abs() < 0.02);
    // <v(0) v(tau)> = cos(0.02 lag) / 2: lag 39 gives cos(0.78)/2, a quarter period (lag 78.5) gives ~0
    assert!((c[39] - 0.5 * 0.78_f64.cos()).abs() < 0.03, "{}", c[39]);
    assert!(c[78].abs() < 0.03, "{}", c[78]);
    let t: Vec<f64> = (0..50).map(f64::from).collect();
    let e: Vec<f64> = t.iter().map(|x| 3.0 + 0.25 * x).collect();
    assert!((energy_drift_slope(&t, &e) - 0.25).abs() < 1e-12);
}

/// LJ fluid in a skewed (triclinic) cell with the volume of the cubic fcc box.
fn lj_triclinic(
    cells: usize,
    rho: f64,
    t: f64,
    dt: f64,
    seed: u64,
) -> MdSimulation<3> {
    let n = 4 * cells * cells * cells;
    let l = (n as f64 / rho).cbrt();
    let a = l / cells as f64;
    let h = [[l, 0.0, 0.0], [0.15 * l, l, 0.0], [0.1 * l, 0.12 * l, l]];
    let bx = SimBox::triclinic(h).unwrap();
    let pos: Vec<[f64; 3]> = lattice_fcc(cells, a)
        .iter()
        .map(|p| {
            let s = [p[0] / l, p[1] / l, p[2] / l];
            bx.from_fractional(&s)
        })
        .collect();
    let mut sys = MdSystem::new(pos, bx);
    let mut rng = Rng::new(seed);
    sys.maxwell_boltzmann(t, &mut rng);
    let mut ff = ForceField::new(lj(), 2.5);
    ff.shift_force = true;
    MdSimulation::new(sys, ff, 0.3, dt, seed + 1).unwrap()
}

#[test]
fn triclinic_cell_energy_matches_an_explicit_image_sum() {
    let sim = lj_triclinic(3, 0.7, 1.0, 0.004, 5);
    let h = sim.sys.bx.matrix();
    assert!((sim.sys.bx.volume() - 108.0 / 0.7).abs() < 1e-9);
    let rc = 2.5;
    let term = PairTerm::LennardJones { eps: 1.0, sigma: 1.0 };
    // shifted-force reference: u(r) - u(rc) + F(rc) (r - rc)
    let (ec, fcr) = term.eval(rc * rc, 0.0);
    let fcut = fcr * rc; // -dU/dr at rc
    let n = sim.sys.n();
    let mut e = 0.0;
    for i in 0..n {
        for j in 0..n {
            for a in -2i32..=2 {
                for b in -2i32..=2 {
                    for c in -2i32..=2 {
                        if i == j && a == 0 && b == 0 && c == 0 {
                            continue;
                        }
                        let mut d = [0.0; 3];
                        for k in 0..3 {
                            d[k] = sim.sys.pos[i][k] - sim.sys.pos[j][k] + a as f64 * h[0][k] + b as f64 * h[1][k] + c as f64 * h[2][k];
                        }
                        let r = d.iter().map(|x| x * x).sum::<f64>().sqrt();
                        if r < rc {
                            e += 0.5 * (term.eval(r * r, 0.0).0 - ec + fcut * (r - rc));
                        }
                    }
                }
            }
        }
    }
    assert!((sim.forces.epot - e).abs() < 1e-8 * e.abs().max(1.0), "{} vs {e}", sim.forces.epot);
    let net: f64 = (0..3).map(|k| sim.forces.f.iter().map(|f| f[k]).sum::<f64>().abs()).sum();
    assert!(net < 1e-9, "net force {net}");
    for p in &sim.sys.pos {
        let s = sim.sys.bx.to_fractional(p);
        assert!(s.iter().all(|x| (-1e-12..1.0 + 1e-12).contains(x)));
    }
}

#[test]
fn flexible_cell_nph_conserves_its_extended_enthalpy_with_second_order_error() {
    let err = |dt: f64| {
        let mut sim = lj_triclinic(3, 0.75, 1.2, dt, 81);
        sim.barostat = Barostat::Flexible(FlexibleBarostat::hydrostatic(3, 1.0, 40.0, CellMode::Full));
        sim.run((0.4 / dt) as usize).unwrap();
        let h0 = sim.conserved_energy();
        let mut worst = 0.0_f64;
        for _ in 0..(0.8 / dt) as usize {
            sim.step().unwrap();
            worst = worst.max((sim.conserved_energy() - h0).abs());
        }
        let h = sim.sys.bx.matrix();
        (worst, h)
    };
    let (e1, h) = err(0.008);
    let (e2, _) = err(0.004);
    assert!(e1 < 0.5, "enthalpy error {e1}");
    assert!(e1 / e2 > 2.0, "ratio {}", e1 / e2);
    let l = (108.0_f64 / 0.75).cbrt();
    assert!((h[0][0] - l).abs() + (h[1][0] - 0.15 * l).abs() > 1e-4);
}

#[test]
fn flexible_cell_npt_with_nose_hoover_chains_reaches_the_target_stress_and_temperature() {
    let (t0, p0) = (1.2, 1.0);
    let mut sim = lj_fluid_3d(3, 0.7, t0, 0.005, 91);
    let dof = sim.dof();
    sim.thermostat = Thermostat::NoseHooverChain(NoseHooverChain::new(t0, 0.5, 3, dof));
    sim.barostat = Barostat::Flexible(FlexibleBarostat::hydrostatic(3, p0, 20.0, CellMode::Full).with_chain(t0, 1.0, 3));
    sim.run(3000).unwrap();
    let (mut tm, mut pt) = (0.0, [[0.0; 3]; 3]);
    let ns = 6000;
    let c0 = sim.conserved_energy();
    let mut cw = 0.0_f64;
    for _ in 0..ns {
        sim.step().unwrap();
        tm += sim.temperature();
        let p = sim.pressure_tensor();
        for a in 0..3 {
            for b in 0..3 {
                pt[a][b] += p[a][b];
            }
        }
        cw = cw.max((sim.conserved_energy() - c0).abs());
    }
    let nsf = ns as f64;
    assert!((tm / nsf - t0).abs() < 0.05, "T = {}", tm / nsf);
    for a in 0..3 {
        assert!((pt[a][a] / nsf - p0).abs() < 0.12, "P[{a}][{a}] = {}", pt[a][a] / nsf);
        for b in 0..3 {
            if a != b {
                assert!((pt[a][b] / nsf).abs() < 0.08, "P[{a}][{b}] = {}", pt[a][b] / nsf);
            }
        }
    }
    assert!(cw < 0.5, "extended energy wandered by {cw}");
}

#[test]
fn orthorhombic_mode_follows_an_anisotropic_target_pressure() {
    let t0 = 1.0;
    let mut sim = lj_fluid_3d(3, 0.6, t0, 0.005, 95);
    let dof = sim.dof();
    sim.thermostat = Thermostat::NoseHooverChain(NoseHooverChain::new(t0, 0.5, 3, dof));
    let pext = vec![0.8, 0.0, 0.0, 0.0, 0.4, 0.0, 0.0, 0.0, 0.2];
    sim.barostat = Barostat::Flexible(FlexibleBarostat::with_tensor(3, pext, 20.0, CellMode::Orthorhombic).with_chain(t0, 1.0, 3));
    sim.run(4000).unwrap();
    let mut p = [0.0; 3];
    let ns = 6000;
    for _ in 0..ns {
        sim.step().unwrap();
        let pt = sim.pressure_tensor();
        for a in 0..3 {
            p[a] += pt[a][a] / ns as f64;
        }
    }
    for (a, target) in [0.8, 0.4, 0.2].iter().enumerate() {
        assert!((p[a] - target).abs() < 0.1, "P[{a}] = {} vs {target}", p[a]);
    }
    assert!(sim.sys.bx.cell.is_none());
}

fn dimer_system(
    nmol: usize,
    l: f64,
    d: f64,
    seed: u64,
) -> MdSimulation<3> {
    let bx = SimBox::periodic([l; 3]);
    let mut rng = Rng::new(seed);
    let centres = random_positions(nmol, &bx, 1.6, &mut rng).unwrap();
    let mut pos = Vec::new();
    for c in &centres {
        let mut u = [rng.normal(), rng.normal(), rng.normal()];
        let nrm = (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]).sqrt();
        for x in &mut u {
            *x /= nrm;
        }
        pos.push([c[0] + 0.5 * d * u[0], c[1] + 0.5 * d * u[1], c[2] + 0.5 * d * u[2]]);
        pos.push([c[0] - 0.5 * d * u[0], c[1] - 0.5 * d * u[1], c[2] - 0.5 * d * u[2]]);
    }
    let mut sys = MdSystem::new(pos, bx);
    sys.maxwell_boltzmann(1.0, &mut rng);
    let mut ff = ForceField::<3>::new(lj(), 2.5);
    ff.shift = true;
    ff.bonds = (0..nmol).map(|m| Bond { i: 2 * m, j: 2 * m + 1, k: 0.0, r0: d }).collect();
    let mut sim = MdSimulation::new(sys, ff, 0.3, 0.004, seed + 1).unwrap();
    sim.constraints = (0..nmol).map(|m| Constraint { i: 2 * m, j: 2 * m + 1, d }).collect();
    sim
}

fn max_bond_error(sim: &MdSimulation<3>) -> f64 {
    let mut worst = 0.0_f64;
    for c in &sim.constraints {
        let mut dv = [0.0; 3];
        for a in 0..3 {
            dv[a] = sim.sys.pos[c.i][a] - sim.sys.pos[c.j][a];
        }
        sim.sys.bx.min_image(&mut dv);
        worst = worst.max((dv.iter().map(|x| x * x).sum::<f64>().sqrt() - c.d).abs());
    }
    worst
}

#[test]
fn flexible_cell_with_constraints_keeps_dimers_rigid_and_conserves_the_extended_energy() {
    let d = 0.7;
    let mut sim = dimer_system(32, 7.0, d, 101);
    for _ in 0..3 {
        sim.step().unwrap();
    }
    let dof = sim.dof();
    sim.thermostat = Thermostat::NoseHooverChain(NoseHooverChain::new(1.0, 0.5, 3, dof));
    sim.barostat = Barostat::Flexible(FlexibleBarostat::hydrostatic(3, 0.5, 30.0, CellMode::Full).with_chain(1.0, 1.0, 3));
    sim.run(200).unwrap();
    let c0 = sim.conserved_energy();
    let (mut worst_len, mut worst_e) = (0.0_f64, 0.0_f64);
    for _ in 0..1500 {
        sim.step().unwrap();
        worst_len = worst_len.max(max_bond_error(&sim));
        worst_e = worst_e.max((sim.conserved_energy() - c0).abs());
    }
    assert!(worst_len < 1e-9, "bond length drift {worst_len}");
    assert!(worst_e < 0.02 * sim.sys.kinetic_energy(), "extended energy drift {worst_e}");
    let tr: f64 = (0..3).map(|a| sim.constraint_virial[a][a]).sum();
    assert!(tr.is_finite() && tr != 0.0);
}

#[test]
fn constrained_langevin_keeps_bonds_rigid_and_samples_the_target_temperature() {
    let d = 0.7;
    let mut sim = dimer_system(32, 7.0, d, 111);
    sim.thermostat = Thermostat::Langevin { t0: 1.5, gamma: 2.0 };
    sim.run(1000).unwrap();
    let (mut tm, mut worst) = (0.0, 0.0_f64);
    let ns = 6000;
    for _ in 0..ns {
        sim.step().unwrap();
        tm += sim.temperature();
        worst = worst.max(max_bond_error(&sim));
    }
    assert!(worst < 1e-9, "bond drift {worst}");
    assert!((tm / ns as f64 - 1.5).abs() < 0.06, "T = {}", tm / ns as f64);
    sim.barostat = Barostat::Flexible(FlexibleBarostat::hydrostatic(3, 1.0, 30.0, CellMode::Isotropic));
    assert!(sim.step().is_err());
}

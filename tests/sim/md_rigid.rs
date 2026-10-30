//! Rigid-molecule MD: geometry, conservation laws, agreement with SHAKE/RATTLE, molecular virial.

use rssn::kernels::random::Rng;
use rssn::sim::md_engine::*;
use rssn::sim::md_rigid::*;

fn lj_ff() -> ForceField<3> {
    let mut ff = ForceField::new(vec![PairTerm::LennardJones { eps: 1.0, sigma: 0.6 }], 2.0);
    ff.shift_force = true;
    ff
}

/// Random rotation (rows of a rotation matrix) from a random unit quaternion.
fn random_rotation(rng: &mut Rng) -> [[f64; 3]; 3] {
    let mut q = [rng.normal(), rng.normal(), rng.normal(), rng.normal()];
    let n = q.iter().map(|x| x * x).sum::<f64>().sqrt();
    for x in &mut q {
        *x /= n;
    }
    let [w, x, y, z] = q;
    [
        [
            1.0 - 2.0 * (y * y + z * z),
            2.0 * (x * y - w * z),
            2.0 * (x * z + w * y),
        ],
        [
            2.0 * (x * y + w * z),
            1.0 - 2.0 * (x * x + z * z),
            2.0 * (y * z - w * x),
        ],
        [
            2.0 * (x * z - w * y),
            2.0 * (y * z + w * x),
            1.0 - 2.0 * (x * x + y * y),
        ],
    ]
}

/// Fluid of rigid molecules with the given site shape and masses.
fn molecular_fluid(
    shape: &[[f64; 3]],
    masses: &[f64],
    nmol: usize,
    l: f64,
    temp: f64,
    seed: u64,
) -> (MdSystem<3>, Vec<Vec<usize>>) {
    let bx = SimBox::periodic([l; 3]);
    let mut rng = Rng::new(seed);
    let centres = random_positions(
        nmol,
        &bx,
        if shape.len() == 1 { 1.0 } else { 2.0 },
        &mut rng,
    )
    .unwrap();
    let ns = shape.len();
    let mut pos = Vec::new();
    let mut mass = Vec::new();
    let mut mols = Vec::new();
    for (m, c) in centres.iter().enumerate() {
        let r = random_rotation(&mut rng);
        let mut ids = Vec::new();
        for (s, b) in shape.iter().enumerate() {
            let mut p = *c;
            for i in 0..3 {
                for j in 0..3 {
                    p[i] += r[i][j] * b[j];
                }
            }
            pos.push(p);
            mass.push(masses[s]);
            ids.push(m * ns + s);
        }
        mols.push(ids);
    }
    let mut sys = MdSystem::new(pos, bx);
    sys.mass = mass;
    for v in &mut sys.vel {
        for c in v.iter_mut() {
            *c = temp.sqrt() * rng.normal();
        }
    }
    (sys, mols)
}

fn water_like() -> ([[f64; 3]; 3], [f64; 3]) {
    (
        [[0.0, 0.0, 0.0], [0.8, 0.1, 0.0], [-0.3, 0.7, 0.2]],
        [1.0, 0.5, 0.7],
    )
}

#[test]
fn rigid_molecules_stay_rigid_and_conserve_energy_and_momentum() {
    let (shape, masses) = water_like();
    let run = |dt: f64| {
        let (sys, mols) = molecular_fluid(&shape, &masses, 27, 7.0, 1.0, 5);
        let mut sim = RigidMdSimulation::new(sys, lj_ff(), 0.3, dt, &mols).unwrap();
        let p0 = sim.momentum();
        sim.run((0.5 / dt) as usize).unwrap();
        let e0 = sim.total_energy();
        let mut worst = 0.0_f64;
        for _ in 0..(1.5 / dt) as usize {
            sim.step().unwrap();
            worst = worst.max((sim.total_energy() - e0).abs());
        }
        let p = sim.momentum();
        assert!(
            (0..3).all(|k| (p[k] - p0[k]).abs() < 1e-10),
            "momentum {p:?} vs {p0:?}"
        );
        assert!(sim.max_geometry_error() < 1e-12);
        assert_eq!(sim.dof(), (27 * 6 - 3) as f64);
        worst
    };
    let e1 = run(0.004);
    let e2 = run(0.002);
    assert!(e1 < 0.05, "energy error {e1}");
    assert!(e1 / e2 > 2.5, "ratio {}", e1 / e2);
}

#[test]
fn rigid_integrator_agrees_with_shake_rattle_for_dimers() {
    let d = 0.7;
    let shape = [[-0.5 * d, 0.0, 0.0], [0.5 * d, 0.0, 0.0]];
    let nmol = 20;
    let compare = |dt: f64| -> f64 {
        let (sys, mols) = molecular_fluid(&shape, &[1.0, 1.0], nmol, 6.5, 1.0, 9);
        let mut rigid = RigidMdSimulation::new(sys.clone(), lj_ff(), 0.3, dt, &mols).unwrap();
        // SHAKE/RATTLE system with the rigid-consistent velocities
        let mut ff = lj_ff();
        ff.bonds = mols
            .iter()
            .map(|m| Bond {
                i: m[0],
                j: m[1],
                k: 0.0,
                r0: d,
            })
            .collect();
        let mut s2 = sys;
        s2.vel = rigid.sys.vel.clone();
        let mut shake = MdSimulation::new(s2, ff, 0.3, dt, 3).unwrap();
        shake.constraints = mols
            .iter()
            .map(|m| Constraint { i: m[0], j: m[1], d })
            .collect();
        let steps = (0.3 / dt) as usize;
        rigid.run(steps).unwrap();
        shake.run(steps).unwrap();
        let mut worst = 0.0_f64;
        for i in 0..rigid.sys.n() {
            let mut dv = [0.0; 3];
            for a in 0..3 {
                dv[a] = rigid.sys.pos[i][a] - shake.sys.pos[i][a];
            }
            rigid.sys.bx.min_image(&mut dv);
            worst = worst.max(dv.iter().map(|x| x * x).sum::<f64>().sqrt());
        }
        worst
    };
    let e1 = compare(0.004);
    let e2 = compare(0.002);
    assert!(e1 < 5e-3, "trajectory difference {e1}");
    assert!(
        e2 < e1 / 2.5,
        "difference should shrink as dt^2: {e1} -> {e2}"
    );
}

#[test]
fn single_site_molecules_reduce_to_plain_velocity_verlet() {
    let (sys, mols) = molecular_fluid(&[[0.0; 3]], &[1.0], 40, 7.5, 1.0, 13);
    let mut rigid = RigidMdSimulation::new(sys.clone(), lj_ff(), 0.3, 0.003, &mols).unwrap();
    let mut plain = MdSimulation::new(sys, lj_ff(), 0.3, 0.003, 1).unwrap();
    // the rigid driver removes nothing from the velocities of point particles
    plain.sys.vel = rigid.sys.vel.clone();
    rigid.run(150).unwrap();
    plain.run(150).unwrap();
    for i in 0..40 {
        for a in 0..3 {
            let mut d = [0.0; 3];
            d[a] = rigid.sys.pos[i][a] - plain.sys.pos[i][a];
            rigid.sys.bx.min_image(&mut d);
            assert!(d[a].abs() < 1e-9, "site {i}: {}", d[a]);
        }
    }
    assert_eq!(rigid.dof(), 117.0);
}

#[test]
fn molecular_virial_equals_the_affine_scaling_derivative_of_the_energy() {
    let (shape, masses) = water_like();
    let (sys, mols) = molecular_fluid(&shape, &masses, 27, 7.0, 1.0, 21);
    let mut sim = RigidMdSimulation::new(sys, lj_ff(), 0.5, 0.004, &mols).unwrap();
    for m in &mut sim.mols {
        m.vcom = [0.0; 3];
    }
    let l0 = sim.sys.bx.len[0];
    let com0: Vec<[f64; 3]> = sim.mols.iter().map(|m| m.com).collect();
    let mut energy_at = |eps: f64| -> f64 {
        sim.sys.bx.len = [l0 * (1.0 + eps); 3];
        for (m, c) in sim.mols.iter_mut().zip(&com0) {
            m.com = [c[0] * (1.0 + eps), c[1] * (1.0 + eps), c[2] * (1.0 + eps)];
        }
        sim.sync_sites();
        sim.forces = compute_forces(&sim.sys, &sim.ff, &mut sim.nl).unwrap();
        sim.forces.epot
    };
    let h = 1e-5;
    let dudeps = (energy_at(h) - energy_at(-h)) / (2.0 * h);
    energy_at(0.0);
    let p = sim.pressure_tensor();
    let tr_w = (p[0][0] + p[1][1] + p[2][2]) * sim.sys.bx.volume();
    assert!(
        (tr_w + dudeps).abs() < 1e-5 * dudeps.abs().max(1.0),
        "{tr_w} vs {}",
        -dudeps
    );
}

#[test]
fn rigid_driver_validates_its_input_and_supports_both_rotor_splittings() {
    let (shape, masses) = water_like();
    let (sys, mols) = molecular_fluid(&shape, &masses, 8, 6.0, 1.0, 31);
    assert!(RigidMdSimulation::new(sys.clone(), lj_ff(), 0.3, 0.004, &mols[..4]).is_err());
    let mut dup = mols.clone();
    dup[1][0] = dup[0][0];
    assert!(RigidMdSimulation::new(sys.clone(), lj_ff(), 0.3, 0.004, &dup).is_err());
    let mut energies = Vec::new();
    for (top, order) in [(true, 2), (false, 2), (true, 4)] {
        let mut sim = RigidMdSimulation::new(sys.clone(), lj_ff(), 0.3, 0.004, &mols).unwrap();
        sim.top_split = top;
        sim.order = order;
        sim.run(200).unwrap();
        energies.push(sim.total_energy());
        assert!(sim.max_geometry_error() < 1e-12);
    }
    assert!((energies[0] - energies[1]).abs() < 0.05 && (energies[0] - energies[2]).abs() < 0.05);
}

#[test]
fn exact_rotor_flow_keeps_geometry_and_has_smaller_energy_drift_than_splitting() {
    let (shape, masses) = water_like();
    let run = |exact: bool, dt: f64| {
        let (sys, mols) = molecular_fluid(&shape, &masses, 27, 7.0, 1.0, 5);
        let mut sim = RigidMdSimulation::new(sys, lj_ff(), 0.3, dt, &mols).unwrap();
        sim.exact_rotor = exact;
        sim.run((0.5 / dt) as usize).unwrap();
        let e0 = sim.total_energy();
        let mut worst = 0.0_f64;
        for _ in 0..(1.0 / dt) as usize {
            sim.step().unwrap();
            worst = worst.max((sim.total_energy() - e0).abs());
        }
        assert!(sim.max_geometry_error() < 1e-10);
        worst
    };
    let split = run(false, 0.008);
    let exact = run(true, 0.008);
    assert!(exact < 0.3, "exact-rotor energy error {exact}");
    assert!(exact <= split * 1.2, "exact {exact} vs splitting {split}");
    // second order in dt
    let e2 = run(true, 0.004);
    assert!(exact / e2 > 2.0, "ratio {}", exact / e2);
}

#[test]
fn exact_constraint_virial_reproduces_the_molecular_pressure_without_a_time_step_lag() {
    // rigid dimers: the atomistic pressure with the exact instantaneous constraint forces has the
    // same trace as the molecular (centre-of-mass) pressure of the rigid-body driver
    let d = 0.7;
    let shape = [[-0.5 * d, 0.0, 0.0], [0.5 * d, 0.0, 0.0]];
    let (sys, mols) = molecular_fluid(&shape, &[1.0, 1.0], 20, 6.5, 1.0, 9);
    let mut rigid = RigidMdSimulation::new(sys, lj_ff(), 0.3, 0.004, &mols).unwrap();
    rigid.run(60).unwrap();
    let mut ff = lj_ff();
    ff.bonds = mols
        .iter()
        .map(|m| Bond {
            i: m[0],
            j: m[1],
            k: 0.0,
            r0: d,
        })
        .collect();
    let mut sim = MdSimulation::new(rigid.sys.clone(), ff, 0.3, 0.004, 3).unwrap();
    sim.constraints = mols
        .iter()
        .map(|m| Constraint { i: m[0], j: m[1], d })
        .collect();
    sim.refresh_constraint_virial();
    let p_atom = sim.pressure_tensor();
    let p_mol = rigid.pressure_tensor();
    let tr = |p: &[[f64; 3]; 3]| p[0][0] + p[1][1] + p[2][2];
    assert!(
        (tr(&p_atom) - tr(&p_mol)).abs() < 1e-9 * (1.0 + tr(&p_mol).abs()),
        "{} vs {}",
        tr(&p_atom),
        tr(&p_mol)
    );
    // the constraint force is a genuine contribution
    let tc =
        sim.constraint_virial[0][0] + sim.constraint_virial[1][1] + sim.constraint_virial[2][2];
    assert!(tc.abs() > 1e-3);
}

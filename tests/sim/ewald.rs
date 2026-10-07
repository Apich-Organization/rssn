//! Ewald summation and PME: Madelung constant of NaCl, consistency and mesh accuracy.

use rssn::sim::ewald::*;
use rssn::sim::md_engine::*;

const MADELUNG_NACL: f64 = 1.747_564_594_633_18;

/// Rock salt in an `n x n x n` arrangement of conventional cells (edge 2, nearest-neighbour distance 1).
fn rock_salt(n: usize) -> (Vec<[f64; 3]>, Vec<f64>, [f64; 3]) {
    let mut pos = Vec::new();
    let mut q = Vec::new();
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                let o = [2.0 * i as f64, 2.0 * j as f64, 2.0 * k as f64];
                for (b, sign) in [
                    ([0.0, 0.0, 0.0], 1.0),
                    ([1.0, 1.0, 0.0], 1.0),
                    ([1.0, 0.0, 1.0], 1.0),
                    ([0.0, 1.0, 1.0], 1.0),
                    ([1.0, 0.0, 0.0], -1.0),
                    ([0.0, 1.0, 0.0], -1.0),
                    ([0.0, 0.0, 1.0], -1.0),
                    ([1.0, 1.0, 1.0], -1.0),
                ] {
                    pos.push([o[0] + b[0], o[1] + b[1], o[2] + b[2]]);
                    q.push(sign);
                }
            }
        }
    }
    let l = 2.0 * n as f64;
    (pos, q, [l; 3])
}

fn perturb(
    pos: &mut [[f64; 3]],
    amp: f64,
) {
    let mut s = 12345u64;
    for p in pos.iter_mut() {
        for c in p.iter_mut() {
            s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            *c += amp * ((s >> 33) as f64 / (1u64 << 31) as f64 - 0.5);
        }
    }
}

#[test]
fn ewald_reproduces_the_madelung_constant_of_rock_salt() {
    let (pos, q, len) = rock_salt(2);
    let (e, f) = ewald_energy_forces(&pos, &q, len, 2.2, 14, 1.99, 1.0);
    // E per ion = -M q^2 / (2 r0) with r0 = 1
    let per_ion = e / pos.len() as f64;
    assert!((per_ion + MADELUNG_NACL / 2.0).abs() < 1e-6, "E/ion = {per_ion}");
    // the perfect crystal is in equilibrium
    for fi in &f {
        for c in fi {
            assert!(c.abs() < 1e-6);
        }
    }
}

#[test]
fn ewald_energy_is_independent_of_the_splitting_parameter() {
    let (mut pos, q, len) = rock_salt(2);
    perturb(&mut pos, 0.2);
    let (e1, f1) = ewald_energy_forces(&pos, &q, len, 2.0, 14, 1.99, 1.0);
    let (e2, f2) = ewald_energy_forces(&pos, &q, len, 2.5, 18, 1.99, 1.0);
    assert!((e1 - e2).abs() < 1e-5 * e1.abs(), "{e1} vs {e2}");
    for i in 0..pos.len() {
        for a in 0..3 {
            assert!((f1[i][a] - f2[i][a]).abs() < 1e-4, "force {i}/{a}: {} vs {}", f1[i][a], f2[i][a]);
        }
    }
}

#[test]
fn ewald_forces_are_the_gradient_of_the_energy() {
    let (mut pos, q, len) = rock_salt(1);
    perturb(&mut pos, 0.3);
    let alpha = 1.4;
    let (_, f) = ewald_energy_forces(&pos, &q, len, alpha, 12, 0.999, 1.0);
    // a real-space cutoff below L/2 = 1 is too short for alpha = 1.4; compare against the same truncated model
    let h = 1e-5;
    for i in [0, 3, 5] {
        for a in 0..3 {
            let mut p = pos.clone();
            p[i][a] += h;
            let (ep, _) = ewald_energy_forces(&p, &q, len, alpha, 12, 0.999, 1.0);
            p[i][a] -= 2.0 * h;
            let (em, _) = ewald_energy_forces(&p, &q, len, alpha, 12, 0.999, 1.0);
            let num = -(ep - em) / (2.0 * h);
            assert!((f[i][a] - num).abs() < 1e-5 * f[i][a].abs().max(1.0), "{i}/{a}: {} vs {num}", f[i][a]);
        }
    }
}

#[test]
fn pme_matches_the_direct_reciprocal_sum() {
    let (mut pos, q, len) = rock_salt(2);
    perturb(&mut pos, 0.3);
    let alpha = 1.1;
    let direct = ewald_reciprocal(&pos, &q, len, alpha, 10, 1.0);
    let mesh = pme_reciprocal(&pos, &q, len, alpha, [32, 32, 32], 6, 1.0);
    assert!((direct.energy - mesh.energy).abs() < 1e-4 * direct.energy.abs(), "{} vs {}", direct.energy, mesh.energy);
    let scale = direct.forces.iter().flat_map(|f| f.iter()).map(|c| c.abs()).fold(0.0, f64::max);
    for i in 0..pos.len() {
        for a in 0..3 {
            assert!((direct.forces[i][a] - mesh.forces[i][a]).abs() < 2e-3 * scale, "{i}/{a}: {} vs {}", direct.forces[i][a], mesh.forces[i][a]);
        }
    }
    for a in 0..3 {
        for b in 0..3 {
            assert!((direct.virial[a][b] - mesh.virial[a][b]).abs() < 1e-3 * direct.virial[0][0].abs().max(1.0));
        }
    }
    // refining the mesh reduces the error
    let coarse = pme_reciprocal(&pos, &q, len, alpha, [16, 16, 16], 4, 1.0);
    let fine = pme_reciprocal(&pos, &q, len, alpha, [32, 32, 32], 6, 1.0);
    assert!((coarse.energy - direct.energy).abs() > (fine.energy - direct.energy).abs());
    // non power-of-two meshes work through the Bluestein FFT
    let odd = pme_reciprocal(&pos, &q, len, alpha, [30, 24, 28], 6, 1.0);
    assert!((odd.energy - direct.energy).abs() < 2e-4 * direct.energy.abs());
}

#[test]
fn md_force_field_with_ewald_and_pme_matches_the_reference() {
    let (mut pos, q, len) = rock_salt(2);
    perturb(&mut pos, 0.25);
    let alpha = 2.2;
    let (e_ref, f_ref) = ewald_energy_forces(&pos, &q, len, alpha, 14, 1.99, 1.0);
    for pme in [false, true] {
        let mut sys = MdSystem::new(pos.clone(), SimBox::periodic(len));
        sys.charge = q.clone();
        let mut ff = ForceField::<3>::new(vec![PairTerm::Coulomb { k: 1.0, alpha }], 1.99);
        ff.extra.push(if pme { pme_extra(alpha, [48, 48, 48], 6, 1.0) } else { ewald_extra(alpha, 14, 1.0) });
        let mut nl = NeighbourList::new(1.99, 0.2);
        let out = compute_forces(&sys, &ff, &mut nl).unwrap();
        let tol = if pme { 1e-3 } else { 1e-8 };
        assert!((out.epot - e_ref).abs() < tol * e_ref.abs(), "pme={pme}: {} vs {e_ref}", out.epot);
        for i in 0..pos.len() {
            for a in 0..3 {
                assert!((out.f[i][a] - f_ref[i][a]).abs() < 50.0 * tol, "pme={pme} {i}/{a}");
            }
        }
    }
}

#[test]
fn charged_cell_gets_the_neutralising_background() {
    // a single unit charge in a cubic box: E = -alpha q^2/sqrt(pi) - pi q^2/(2 V alpha^2) + recip; independent of alpha
    let pos = vec![[1.0, 1.0, 1.0]];
    let q = vec![1.0];
    let len = [4.0; 3];
    let (e1, _) = ewald_energy_forces(&pos, &q, len, 1.5, 12, 1.99, 1.0);
    let (e2, _) = ewald_energy_forces(&pos, &q, len, 2.0, 16, 1.99, 1.0);
    assert!((e1 - e2).abs() < 1e-5, "{e1} vs {e2}");
    // known Wigner-lattice value: -2.837297 q^2 / (2 L)
    assert!((e1 + 2.837_297 / (2.0 * 4.0)).abs() < 1e-4, "{e1}");
}

#[test]
fn bspline_partition_of_unity() {
    for u in [0.0, 0.3, 0.77] {
        let s: f64 = (0..6).map(|j| bspline(6, u + j as f64)).sum();
        assert!((s - 1.0).abs() < 1e-13);
    }
}

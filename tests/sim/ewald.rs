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
            s = s
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
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
    assert!(
        (per_ion + MADELUNG_NACL / 2.0).abs() < 1e-6,
        "E/ion = {per_ion}"
    );
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
            assert!(
                (f1[i][a] - f2[i][a]).abs() < 1e-4,
                "force {i}/{a}: {} vs {}",
                f1[i][a],
                f2[i][a]
            );
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
            assert!(
                (f[i][a] - num).abs() < 1e-5 * f[i][a].abs().max(1.0),
                "{i}/{a}: {} vs {num}",
                f[i][a]
            );
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
    assert!(
        (direct.energy - mesh.energy).abs() < 1e-4 * direct.energy.abs(),
        "{} vs {}",
        direct.energy,
        mesh.energy
    );
    let scale = direct
        .forces
        .iter()
        .flat_map(|f| f.iter())
        .map(|c| c.abs())
        .fold(0.0, f64::max);
    for i in 0..pos.len() {
        for a in 0..3 {
            assert!(
                (direct.forces[i][a] - mesh.forces[i][a]).abs() < 2e-3 * scale,
                "{i}/{a}: {} vs {}",
                direct.forces[i][a],
                mesh.forces[i][a]
            );
        }
    }
    for a in 0..3 {
        for b in 0..3 {
            assert!(
                (direct.virial[a][b] - mesh.virial[a][b]).abs()
                    < 1e-3 * direct.virial[0][0].abs().max(1.0)
            );
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
        ff.extra.push(if pme {
            pme_extra(alpha, [48, 48, 48], 6, 1.0)
        } else {
            ewald_extra(alpha, 14, 1.0)
        });
        let mut nl = NeighbourList::new(1.99, 0.2);
        let out = compute_forces(&sys, &ff, &mut nl).unwrap();
        let tol = if pme { 1e-3 } else { 1e-8 };
        assert!(
            (out.epot - e_ref).abs() < tol * e_ref.abs(),
            "pme={pme}: {} vs {e_ref}",
            out.epot
        );
        for i in 0..pos.len() {
            for a in 0..3 {
                assert!(
                    (out.f[i][a] - f_ref[i][a]).abs() < 50.0 * tol,
                    "pme={pme} {i}/{a}"
                );
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

fn triclinic_cell() -> [[f64; 3]; 3] {
    [[5.0, 0.0, 0.0], [1.3, 4.6, 0.0], [-0.8, 1.1, 5.2]]
}

fn neutral_charges(h: &[[f64; 3]; 3]) -> (Vec<[f64; 3]>, Vec<f64>) {
    let mut rng = rssn::kernels::random::Rng::new(77);
    let n = 8;
    let mut pos = Vec::new();
    let mut q = Vec::new();
    for i in 0..n {
        let s = [rng.uniform(), rng.uniform(), rng.uniform()];
        let mut p = [0.0; 3];
        for c in 0..3 {
            p[c] = s[0] * h[0][c] + s[1] * h[1][c] + s[2] * h[2][c];
        }
        pos.push(p);
        q.push(if i % 2 == 0 { 1.0 } else { -1.0 });
    }
    (pos, q)
}

#[test]
fn triclinic_ewald_gives_the_madelung_constant_from_the_fcc_primitive_cell() {
    // rock salt: primitive fcc cell (rhombohedral, 60 degrees), nearest-neighbour distance 1
    let h = [[0.0, 1.0, 1.0], [1.0, 0.0, 1.0], [1.0, 1.0, 0.0]];
    let pos = vec![[0.0; 3], [1.0, 0.0, 0.0]];
    let q = vec![1.0, -1.0];
    let (e, f) = ewald_energy_forces_cell(&pos, &q, &h, 2.4, 14, 3.9, 1.0).unwrap();
    assert!(
        (e + MADELUNG_NACL).abs() < 1e-8,
        "{e} vs {}",
        -MADELUNG_NACL
    );
    for fi in &f {
        assert!(fi.iter().all(|c| c.abs() < 1e-8));
    }
    let r = pme_reciprocal_cell(&pos, &q, &h, 2.4, [32, 32, 32], 6, 1.0).unwrap();
    let d = ewald_reciprocal_cell(&pos, &q, &h, 2.4, 14, 1.0).unwrap();
    assert!(
        (r.energy - d.energy).abs() < 1e-6,
        "{} vs {}",
        r.energy,
        d.energy
    );
}

#[test]
fn triclinic_ewald_is_independent_of_alpha_and_forces_are_gradients() {
    let h = triclinic_cell();
    let (pos, q) = neutral_charges(&h);
    let (e1, f1) = ewald_energy_forces_cell(&pos, &q, &h, 1.0, 12, 8.5, 1.0).unwrap();
    let (e2, f2) = ewald_energy_forces_cell(&pos, &q, &h, 1.4, 16, 6.5, 1.0).unwrap();
    assert!((e1 - e2).abs() < 1e-6 * e1.abs().max(1.0), "{e1} vs {e2}");
    for i in 0..pos.len() {
        for a in 0..3 {
            assert!(
                (f1[i][a] - f2[i][a]).abs() < 1e-5,
                "force {i},{a}: {} vs {}",
                f1[i][a],
                f2[i][a]
            );
        }
    }
    let alpha = 1.2;
    let (_, f) = ewald_energy_forces_cell(&pos, &q, &h, alpha, 14, 7.5, 1.0).unwrap();
    for (i, a) in [(0usize, 0usize), (3, 1), (5, 2)] {
        let eps = 1e-5;
        let mut p = pos.clone();
        p[i][a] += eps;
        let (ep, _) = ewald_energy_forces_cell(&p, &q, &h, alpha, 14, 7.5, 1.0).unwrap();
        p[i][a] -= 2.0 * eps;
        let (em, _) = ewald_energy_forces_cell(&p, &q, &h, alpha, 14, 7.5, 1.0).unwrap();
        let fd = -(ep - em) / (2.0 * eps);
        assert!((fd - f[i][a]).abs() < 1e-6, "{fd} vs {}", f[i][a]);
    }
}

#[test]
fn triclinic_pme_matches_the_direct_sum_and_the_virial_is_the_strain_derivative() {
    let h = triclinic_cell();
    let (pos, q) = neutral_charges(&h);
    let alpha = 0.9;
    let direct = ewald_reciprocal_cell(&pos, &q, &h, alpha, 12, 1.0).unwrap();
    let mesh = pme_reciprocal_cell(&pos, &q, &h, alpha, [40, 40, 40], 6, 1.0).unwrap();
    assert!(
        (direct.energy - mesh.energy).abs() < 2e-5 * direct.energy.abs().max(1.0),
        "{} vs {}",
        direct.energy,
        mesh.energy
    );
    for i in 0..pos.len() {
        for a in 0..3 {
            assert!(
                (direct.forces[i][a] - mesh.forces[i][a]).abs() < 2e-4,
                "{} vs {}",
                direct.forces[i][a],
                mesh.forces[i][a]
            );
        }
    }
    // virial = -dE/d(strain) for symmetric strains (affine deformation of cell and positions)
    let strained = |eps: [[f64; 3]; 3]| {
        let tf = |v: &[f64; 3]| -> [f64; 3] {
            std::array::from_fn(|c| v[c] + (0..3).map(|d| v[d] * eps[d][c]).sum::<f64>())
        };
        let h2 = [tf(&h[0]), tf(&h[1]), tf(&h[2])];
        let p2: Vec<[f64; 3]> = pos.iter().map(tf).collect();
        ewald_reciprocal_cell(&p2, &q, &h2, alpha, 12, 1.0)
            .unwrap()
            .energy
    };
    for (x, y) in [(0usize, 0usize), (1, 1), (0, 1), (2, 0)] {
        let d = 1e-6;
        let mut e = [[0.0; 3]; 3];
        e[x][y] = d;
        e[y][x] = d;
        let ep = strained(e);
        let mut e = [[0.0; 3]; 3];
        e[x][y] = -d;
        e[y][x] = -d;
        let em = strained(e);
        let mut fd = -(ep - em) / (2.0 * d);
        if x != y {
            fd /= 2.0;
        }
        assert!(
            (fd - direct.virial[x][y]).abs() < 1e-5 * (1.0 + fd.abs()),
            "virial {x}{y}: {fd} vs {}",
            direct.virial[x][y]
        );
    }
    assert!(pme_reciprocal_cell(&pos, &q, &h, alpha, [4, 4, 4], 6, 1.0).is_err());
    assert!(
        ewald_reciprocal_cell(
            &pos,
            &q,
            &[[1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
            alpha,
            4,
            1.0
        )
        .is_err()
    );
}

#[test]
fn p3m_with_the_optimal_influence_function_beats_pme_on_a_coarse_mesh() {
    let h = triclinic_cell();
    let (pos, q) = neutral_charges(&h);
    let alpha = 0.9;
    let direct = ewald_reciprocal_cell(&pos, &q, &h, alpha, 12, 1.0).unwrap();
    let rms = |f: &[[f64; 3]]| -> f64 {
        let s: f64 = f
            .iter()
            .zip(&direct.forces)
            .map(|(a, b)| (0..3).map(|c| (a[c] - b[c]).powi(2)).sum::<f64>())
            .sum();
        (s / f.len() as f64).sqrt()
    };
    for (grid, order) in [
        ([12usize, 12, 12], 4usize),
        ([16, 16, 16], 4),
        ([16, 16, 16], 6),
    ] {
        let pme = pme_reciprocal_cell(&pos, &q, &h, alpha, grid, order, 1.0).unwrap();
        let p3m = p3m_reciprocal_cell(&pos, &q, &h, alpha, grid, order, 1.0).unwrap();
        let (e_pme, e_p3m) = (
            (pme.energy - direct.energy).abs(),
            (p3m.energy - direct.energy).abs(),
        );
        let (f_pme, f_p3m) = (rms(&pme.forces), rms(&p3m.forces));
        assert!(
            f_p3m < f_pme,
            "grid {grid:?} order {order}: P3M force error {f_p3m} vs PME {f_pme}"
        );
        assert!(
            e_p3m < 3.0 * e_pme + 1e-8,
            "energy errors {e_p3m} vs {e_pme}"
        );
        assert!(f_p3m < 0.05, "P3M force error {f_p3m}");
    }
    // fine mesh: accurate forces, energy and virial
    let p3m = p3m_reciprocal_cell(&pos, &q, &h, alpha, [32, 32, 32], 6, 1.0).unwrap();
    assert!((p3m.energy - direct.energy).abs() < 1e-6 * direct.energy.abs().max(1.0));
    for i in 0..pos.len() {
        for a in 0..3 {
            assert!(
                (p3m.forces[i][a] - direct.forces[i][a]).abs() < 1e-3,
                "{} vs {}",
                p3m.forces[i][a],
                direct.forces[i][a]
            );
        }
    }
    for a in 0..3 {
        for b in 0..3 {
            assert!(
                (p3m.virial[a][b] - direct.virial[a][b]).abs()
                    < 1e-3 * (1.0 + direct.virial[a][b].abs())
            );
        }
    }
    assert!(p3m_reciprocal_cell(&pos, &q, &h, alpha, [4, 4, 4], 6, 1.0).is_err());
    assert!(p3m_reciprocal_cell(&pos, &q, &h, alpha, [16, 16, 16], 5, 1.0).is_err());
}

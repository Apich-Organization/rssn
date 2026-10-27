//! SSE directed loops vs exact diagonalisation and Bethe-ansatz / literature energies.

use rssn::sim::sse::*;

/// Dense XXZ Hamiltonian (`Jxy = 1`, antiferromagnetic) on the given bonds.
fn xxz_matrix(n: usize, bonds: &[(usize, usize)], jz: f64) -> Vec<Vec<f64>> {
    let dim = 1usize << n;
    let mut h = vec![vec![0.0; dim]; dim];
    for s in 0..dim {
        for &(i, j) in bonds {
            let (a, b) = ((s >> i) & 1, (s >> j) & 1);
            let (za, zb) = (a as f64 - 0.5, b as f64 - 0.5);
            h[s][s] += jz * za * zb;
            if a != b {
                let t = s ^ (1 << i) ^ (1 << j);
                h[t][s] += 0.5;
            }
        }
    }
    h
}

/// Eigenvalues of a symmetric matrix by cyclic Jacobi rotations.
fn eigenvalues(mut a: Vec<Vec<f64>>) -> Vec<f64> {
    let n = a.len();
    for _ in 0..60 {
        let mut off = 0.0;
        for i in 0..n {
            for j in i + 1..n {
                off += a[i][j] * a[i][j];
            }
        }
        if off < 1e-22 {
            break;
        }
        for p in 0..n {
            for q in p + 1..n {
                if a[p][q].abs() < 1e-14 {
                    continue;
                }
                let theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
                let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
                let t = if theta == 0.0 { 1.0 } else { t };
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                for k in 0..n {
                    let (akp, akq) = (a[k][p], a[k][q]);
                    a[k][p] = c * akp - s * akq;
                    a[k][q] = s * akp + c * akq;
                }
                for k in 0..n {
                    let (apk, aqk) = (a[p][k], a[q][k]);
                    a[p][k] = c * apk - s * aqk;
                    a[q][k] = s * apk + c * aqk;
                }
            }
        }
    }
    (0..n).map(|i| a[i][i]).collect()
}

fn thermal(evs: &[f64], beta: f64) -> (f64, f64) {
    let e0 = evs.iter().copied().fold(f64::INFINITY, f64::min);
    let (mut z, mut e1, mut e2) = (0.0, 0.0, 0.0);
    for &e in evs {
        let w = (-beta * (e - e0)).exp();
        z += w;
        e1 += w * e;
        e2 += w * e * e;
    }
    let (e1, e2) = (e1 / z, e2 / z);
    (e1, beta * beta * (e2 - e1 * e1))
}

#[test]
fn sse_reproduces_exact_diagonalisation_of_small_xxz_systems() {
    type Case = (usize, Vec<(usize, usize)>, f64, f64);
    let cases: [Case; 4] = [
        (6, chain_bonds(6, true), 1.0, 2.0),
        (6, chain_bonds(6, true), 0.5, 1.5),
        (6, chain_bonds(6, true), 0.0, 2.0),
        (6, square_bonds(3, 2, false), 1.0, 2.0),
    ];
    for (n, bonds, jz, beta) in cases {
        let evs = eigenvalues(xxz_matrix(n, &bonds, jz));
        let (e_ex, cv_ex) = thermal(&evs, beta);
        let mut sim = SseXxz::new(n, bonds, jz, beta, 17).unwrap();
        let r = sim.run(3000, 40_000);
        assert!((r.energy - e_ex).abs() < 5.0 * r.energy_error + 0.01, "n={n} jz={jz}: E = {} +- {} vs exact {e_ex}", r.energy, r.energy_error);
        assert!((r.specific_heat - cv_ex).abs() < 0.12 * cv_ex.abs() + 0.08, "Cv {} vs {cv_ex}", r.specific_heat);
    }
}

#[test]
fn heisenberg_chain_ground_state_energy_matches_the_bethe_ansatz() {
    // E0/N -> 1/4 - ln 2 = -0.443147 ; finite size: - pi^2/(12 L^2)
    let l = 24;
    let mut sim = SseXxz::new(l, chain_bonds(l, true), 1.0, 2.0 * l as f64, 5).unwrap();
    let r = sim.run(2000, 20_000);
    let e_per_site = r.energy / l as f64;
    let expected = 0.25 - 2.0_f64.ln() - std::f64::consts::PI.powi(2) / (12.0 * (l * l) as f64);
    assert!((e_per_site - expected).abs() < 6e-4, "E/N = {e_per_site} vs {expected}");
    // staggered structure factor is large for the critical AF chain
    assert!(r.staggered_sq / l as f64 > 0.2);
}

#[test]
fn square_lattice_heisenberg_matches_the_known_4x4_ground_state_energy() {
    // exact 4x4 periodic cluster: E0 = -11.228483 (-0.701780 per site)
    let mut sim = SseXxz::new(16, square_bonds(4, 4, true), 1.0, 16.0, 9).unwrap();
    let r = sim.run(2000, 10_000);
    assert!((r.energy / 16.0 + 0.701_780).abs() < 2.5e-3, "E/N = {}", r.energy / 16.0);
    assert!(r.magnetisation_sq < 1e-2 * 16.0);
}

#[test]
fn sse_rejects_frustrated_or_invalid_lattices() {
    assert!(SseXxz::new(3, vec![(0, 1), (1, 2), (2, 0)], 1.0, 1.0, 1).is_err());
    assert!(SseXxz::new(2, vec![], 1.0, 1.0, 1).is_err());
    assert!(SseXxz::new(2, vec![(0, 1)], 1.0, -1.0, 1).is_err());
}

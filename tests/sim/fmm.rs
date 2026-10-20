//! Softened FMM far field and periodic FMM vs Ewald summation.

use rssn::kernels::random::Rng;
use rssn::sim::ewald::ewald_energy_forces_cell;
use rssn::sim::fmm::*;

fn random_system(n: usize, h: &[[f64; 3]; 3], seed: u64, neutral: bool) -> (Vec<[f64; 3]>, Vec<f64>) {
    let mut rng = Rng::new(seed);
    let mut pos = Vec::new();
    let mut q = Vec::new();
    for i in 0..n {
        let s = [rng.uniform(), rng.uniform(), rng.uniform()];
        pos.push(std::array::from_fn(|c| s[0] * h[0][c] + s[1] * h[1][c] + s[2] * h[2][c]));
        q.push(if neutral {
            if i % 2 == 0 {
                1.0
            } else {
                -1.0
            }
        } else {
            0.5 + rng.uniform()
        });
    }
    (pos, q)
}

#[test]
fn softened_far_field_matches_the_direct_softened_sum() {
    let mut rng = Rng::new(5);
    let n = 400;
    let pos: Vec<[f64; 3]> = (0..n).map(|_| [rng.uniform(), rng.uniform(), rng.uniform()]).collect();
    let src: Vec<f64> = (0..n).map(|_| rng.uniform() - 0.2).collect();
    let eps = 0.15;
    // direct softened sums
    let mut pot = vec![0.0; n];
    let mut fld = vec![[0.0; 3]; n];
    for i in 0..n {
        for j in 0..n {
            if i != j {
                let d = [pos[j][0] - pos[i][0], pos[j][1] - pos[i][1], pos[j][2] - pos[i][2]];
                let r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2] + eps * eps;
                pot[i] += src[j] / r2.sqrt();
                for k in 0..3 {
                    fld[i][k] += src[j] * d[k] / (r2 * r2.sqrt());
                }
            }
        }
    }
    let mut prev = f64::INFINITY;
    for order in [4, 6, 8] {
        let r = fmm3d(&pos, &src, &FmmOptions { order, leaf_size: 20, theta: 0.5, eps }).unwrap();
        let fscale = fld.iter().map(|f| f[0].abs() + f[1].abs() + f[2].abs()).fold(0.0, f64::max);
        let err = (0..n).map(|i| (0..3).map(|k| (r.field[i][k] - fld[i][k]).abs()).fold(0.0, f64::max)).fold(0.0, f64::max) / fscale;
        let perr = (0..n).map(|i| (r.potential[i] - pot[i]).abs()).fold(0.0, f64::max) / pot.iter().map(|x| x.abs()).fold(0.0, f64::max);
        assert!(err < prev, "error must decrease with the order: {err} vs {prev}");
        prev = err;
        if order == 8 {
            assert!(err < 2e-4 && perr < 2e-4, "softened FMM error {err} {perr}");
        }
    }
}

fn check_periodic(h: [[f64; 3]; 3], neutral: bool, tol_f: f64, tol_e: f64) {
    let (pos, q) = random_system(24, &h, 3, neutral);
    let opts = PeriodicFmmOptions { fmm: FmmOptions { order: 8, leaf_size: 40, theta: 0.5, eps: 0.0 }, near_shells: 1, far_half_order: 8 };
    let r = periodic_fmm3d(&pos, &q, &h, &opts).unwrap();
    let (e_ref, f_ref) = ewald_energy_forces_cell(&pos, &q, &h, 0.9, 14, 9.0, 1.0).unwrap();
    let fmax = f_ref.iter().map(|f| f[0].abs().max(f[1].abs()).max(f[2].abs())).fold(0.0, f64::max);
    for i in 0..pos.len() {
        for a in 0..3 {
            let f = -q[i] * r.field[i][a];
            assert!((f - f_ref[i][a]).abs() < tol_f * fmax, "force {i},{a}: {f} vs {}", f_ref[i][a]);
        }
    }
    let e: f64 = 0.5 * pos.iter().enumerate().map(|(i, _)| q[i] * r.potential[i]).sum::<f64>();
    assert!((e - e_ref).abs() < tol_e * e_ref.abs().max(1.0), "energy {e} vs {e_ref}");
}

#[test]
fn periodic_fmm_matches_ewald_in_cubic_and_triclinic_cells() {
    check_periodic([[5.0, 0.0, 0.0], [0.0, 5.0, 0.0], [0.0, 0.0, 5.0]], true, 2e-4, 2e-4);
    check_periodic([[5.0, 0.0, 0.0], [1.3, 4.6, 0.0], [-0.8, 1.1, 5.2]], true, 5e-4, 5e-4);
    // charged cell: neutralising background in the potential and energy
    check_periodic([[5.0, 0.0, 0.0], [0.0, 5.0, 0.0], [0.0, 0.0, 5.0]], false, 5e-4, 5e-4);
}

#[test]
fn periodic_fmm_is_invariant_under_lattice_translations_and_converges_with_the_far_order() {
    let h = [[5.0, 0.0, 0.0], [1.0, 4.5, 0.0], [0.5, -0.7, 5.1]];
    let (pos, q) = random_system(16, &h, 8, true);
    let shifted: Vec<[f64; 3]> = pos
        .iter()
        .enumerate()
        .map(|(i, x)| if i % 3 == 0 { [x[0] + h[0][0] + 2.0 * h[1][0], x[1] + 2.0 * h[1][1], x[2] + h[2][2] * 0.0] } else { *x })
        .collect();
    let base = PeriodicFmmOptions { fmm: FmmOptions { order: 8, leaf_size: 40, theta: 0.5, eps: 0.0 }, near_shells: 1, far_half_order: 8 };
    let a = periodic_fmm3d(&pos, &q, &h, &base).unwrap();
    let b = periodic_fmm3d(&shifted, &q, &h, &base).unwrap();
    for i in 0..pos.len() {
        for k in 0..3 {
            assert!((a.field[i][k] - b.field[i][k]).abs() < 1e-8 * (1.0 + a.field[i][k].abs()));
        }
    }
    // reference with a very high far order
    let reference = periodic_fmm3d(&pos, &q, &h, &PeriodicFmmOptions { far_half_order: 10, ..base }).unwrap();
    let err = |ho: usize| {
        let r = periodic_fmm3d(&pos, &q, &h, &PeriodicFmmOptions { far_half_order: ho, ..base }).unwrap();
        (0..pos.len()).map(|i| (0..3).map(|k| (r.field[i][k] - reference.field[i][k]).abs()).fold(0.0, f64::max)).fold(0.0, f64::max)
    };
    let (e2, e4) = (err(1), err(6));
    assert!(e4 < e2, "far-order convergence: {e2} -> {e4}");
    assert!(periodic_fmm3d(&[], &[], &h, &base).is_err());
    assert!(periodic_fmm3d(&pos, &q, &[[1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 0.0, 1.0]], &base).is_err());
}

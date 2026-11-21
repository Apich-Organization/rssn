//! Softened FMM far field and periodic FMM vs Ewald summation.

use rssn::kernels::random::Rng;
use rssn::sim::ewald::ewald_energy_forces_cell;
use rssn::sim::fmm::*;

fn random_system(
    n: usize,
    h: &[[f64; 3]; 3],
    seed: u64,
    neutral: bool,
) -> (Vec<[f64; 3]>, Vec<f64>) {
    let mut rng = Rng::new(seed);
    let mut pos = Vec::new();
    let mut q = Vec::new();
    for i in 0..n {
        let s = [rng.uniform(), rng.uniform(), rng.uniform()];
        pos.push(std::array::from_fn(|c| {
            s[0] * h[0][c] + s[1] * h[1][c] + s[2] * h[2][c]
        }));
        q.push(if neutral {
            if i % 2 == 0 { 1.0 } else { -1.0 }
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
    let pos: Vec<[f64; 3]> = (0..n)
        .map(|_| [rng.uniform(), rng.uniform(), rng.uniform()])
        .collect();
    let src: Vec<f64> = (0..n).map(|_| rng.uniform() - 0.2).collect();
    let eps = 0.15;
    // direct softened sums
    let mut pot = vec![0.0; n];
    let mut fld = vec![[0.0; 3]; n];
    for i in 0..n {
        for j in 0..n {
            if i != j {
                let d = [
                    pos[j][0] - pos[i][0],
                    pos[j][1] - pos[i][1],
                    pos[j][2] - pos[i][2],
                ];
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
        let r = fmm3d(
            &pos,
            &src,
            &FmmOptions {
                order,
                leaf_size: 20,
                theta: 0.5,
                eps,
            },
        )
        .unwrap();
        let fscale = fld
            .iter()
            .map(|f| f[0].abs() + f[1].abs() + f[2].abs())
            .fold(0.0, f64::max);
        let err = (0..n)
            .map(|i| {
                (0..3)
                    .map(|k| (r.field[i][k] - fld[i][k]).abs())
                    .fold(0.0, f64::max)
            })
            .fold(0.0, f64::max)
            / fscale;
        let perr = (0..n)
            .map(|i| (r.potential[i] - pot[i]).abs())
            .fold(0.0, f64::max)
            / pot.iter().map(|x| x.abs()).fold(0.0, f64::max);
        assert!(
            err < prev,
            "error must decrease with the order: {err} vs {prev}"
        );
        prev = err;
        if order == 8 {
            assert!(err < 2e-4 && perr < 2e-4, "softened FMM error {err} {perr}");
        }
    }
}

fn check_periodic(
    h: [[f64; 3]; 3],
    neutral: bool,
    tol_f: f64,
    tol_e: f64,
) {
    let (pos, q) = random_system(24, &h, 3, neutral);
    let opts = PeriodicFmmOptions {
        fmm: FmmOptions {
            order: 8,
            leaf_size: 40,
            theta: 0.5,
            eps: 0.0,
        },
        near_shells: 1,
        far_half_order: 8,
        soften_far: true,
    };
    let r = periodic_fmm3d(&pos, &q, &h, &opts).unwrap();
    let (e_ref, f_ref) = ewald_energy_forces_cell(&pos, &q, &h, 0.9, 14, 9.0, 1.0).unwrap();
    let fmax = f_ref
        .iter()
        .map(|f| f[0].abs().max(f[1].abs()).max(f[2].abs()))
        .fold(0.0, f64::max);
    for i in 0..pos.len() {
        for a in 0..3 {
            let f = -q[i] * r.field[i][a];
            assert!(
                (f - f_ref[i][a]).abs() < tol_f * fmax,
                "force {i},{a}: {f} vs {}",
                f_ref[i][a]
            );
        }
    }
    let e: f64 = 0.5
        * pos
            .iter()
            .enumerate()
            .map(|(i, _)| q[i] * r.potential[i])
            .sum::<f64>();
    assert!(
        (e - e_ref).abs() < tol_e * e_ref.abs().max(1.0),
        "energy {e} vs {e_ref}"
    );
}

#[test]
fn periodic_fmm_matches_ewald_in_cubic_and_triclinic_cells() {
    check_periodic(
        [[5.0, 0.0, 0.0], [0.0, 5.0, 0.0], [0.0, 0.0, 5.0]],
        true,
        2e-4,
        2e-4,
    );
    check_periodic(
        [[5.0, 0.0, 0.0], [1.3, 4.6, 0.0], [-0.8, 1.1, 5.2]],
        true,
        5e-4,
        5e-4,
    );
    // charged cell: neutralising background in the potential and energy
    check_periodic(
        [[5.0, 0.0, 0.0], [0.0, 5.0, 0.0], [0.0, 0.0, 5.0]],
        false,
        5e-4,
        5e-4,
    );
}

#[test]
fn periodic_fmm_is_invariant_under_lattice_translations_and_converges_with_the_far_order() {
    let h = [[5.0, 0.0, 0.0], [1.0, 4.5, 0.0], [0.5, -0.7, 5.1]];
    let (pos, q) = random_system(16, &h, 8, true);
    let shifted: Vec<[f64; 3]> = pos
        .iter()
        .enumerate()
        .map(|(i, x)| {
            if i % 3 == 0 {
                [
                    x[0] + h[0][0] + 2.0 * h[1][0],
                    x[1] + 2.0 * h[1][1],
                    x[2] + h[2][2] * 0.0,
                ]
            } else {
                *x
            }
        })
        .collect();
    let base = PeriodicFmmOptions {
        fmm: FmmOptions {
            order: 8,
            leaf_size: 40,
            theta: 0.5,
            eps: 0.0,
        },
        near_shells: 1,
        far_half_order: 8,
        soften_far: true,
    };
    let a = periodic_fmm3d(&pos, &q, &h, &base).unwrap();
    let b = periodic_fmm3d(&shifted, &q, &h, &base).unwrap();
    for i in 0..pos.len() {
        for k in 0..3 {
            assert!((a.field[i][k] - b.field[i][k]).abs() < 1e-8 * (1.0 + a.field[i][k].abs()));
        }
    }
    // reference with a very high far order
    let reference = periodic_fmm3d(
        &pos,
        &q,
        &h,
        &PeriodicFmmOptions {
            far_half_order: 10,
            ..base
        },
    )
    .unwrap();
    let err = |ho: usize| {
        let r = periodic_fmm3d(
            &pos,
            &q,
            &h,
            &PeriodicFmmOptions {
                far_half_order: ho,
                ..base
            },
        )
        .unwrap();
        (0..pos.len())
            .map(|i| {
                (0..3)
                    .map(|k| (r.field[i][k] - reference.field[i][k]).abs())
                    .fold(0.0, f64::max)
            })
            .fold(0.0, f64::max)
    };
    let (e2, e4) = (err(1), err(6));
    assert!(e4 < e2, "far-order convergence: {e2} -> {e4}");
    assert!(periodic_fmm3d(&[], &[], &h, &base).is_err());
    assert!(
        periodic_fmm3d(
            &pos,
            &q,
            &[[1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
            &base
        )
        .is_err()
    );
}

#[test]
fn periodic_fmm_far_lattice_is_consistently_softened() {
    // With a large softening the unsoftened far lattice is visibly inconsistent with an
    // explicit sum over more near images; the Plummer form factor in the far lattice removes it.
    let h = [[5.0, 0.0, 0.0], [0.0, 5.0, 0.0], [0.0, 0.0, 5.0]];
    let (pos, q) = random_system(8, &h, 21, true);
    let fmm = FmmOptions {
        order: 8,
        leaf_size: 40,
        theta: 0.5,
        eps: 1.2,
    };
    let reference = periodic_fmm3d(
        &pos,
        &q,
        &h,
        &PeriodicFmmOptions {
            fmm,
            near_shells: 3,
            far_half_order: 8,
            soften_far: true,
        },
    )
    .unwrap();
    let err = |soften: bool| {
        let r = periodic_fmm3d(
            &pos,
            &q,
            &h,
            &PeriodicFmmOptions {
                fmm,
                near_shells: 1,
                far_half_order: 8,
                soften_far: soften,
            },
        )
        .unwrap();
        (0..pos.len())
            .map(|i| {
                (0..3)
                    .map(|k| (r.field[i][k] - reference.field[i][k]).abs())
                    .fold(0.0, f64::max)
            })
            .fold(0.0, f64::max)
    };
    let (e_soft, e_plain) = (err(true), err(false));
    assert!(
        e_soft < 0.2 * e_plain,
        "softened far lattice {e_soft} vs unsoftened {e_plain}"
    );
    assert!(e_soft < 3e-4, "residual {e_soft}");
}

// ---------------------------------------------------------------------------
// Expansion order: arbitrary p, geometric truncation error, linear cost
// ---------------------------------------------------------------------------

/// Direct all-pairs potential `sum_j s_j / r_ij` and field
/// `sum_j s_j (x_j - x_i) / r_ij^3`, self excluded.
fn direct_sum(
    pos: &[[f64; 3]],
    src: &[f64],
) -> (Vec<f64>, Vec<[f64; 3]>) {
    let n = pos.len();
    let mut pot = vec![0.0; n];
    let mut fld = vec![[0.0; 3]; n];
    for i in 0..n {
        for j in 0..n {
            if i == j {
                continue;
            }
            let d = [
                pos[j][0] - pos[i][0],
                pos[j][1] - pos[i][1],
                pos[j][2] - pos[i][2],
            ];
            let r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
            let ir = 1.0 / r2.sqrt();
            pot[i] += src[j] * ir;
            for k in 0..3 {
                fld[i][k] += src[j] * d[k] * ir / r2;
            }
        }
    }
    (pot, fld)
}

/// Relative `L2` error of `a` against the reference `b`.
fn rel_err(
    a: &[f64],
    b: &[f64],
) -> f64 {
    let num: f64 = a
        .iter()
        .zip(b)
        .map(|(x, y)| (x - y) * (x - y))
        .sum::<f64>()
        .sqrt();
    let den: f64 = b.iter().map(|x| x * x).sum::<f64>().sqrt();
    num / den
}

fn rel_err_vec(
    a: &[[f64; 3]],
    b: &[[f64; 3]],
) -> f64 {
    let af: Vec<f64> = a.iter().flat_map(|x| x.iter().copied()).collect();
    let bf: Vec<f64> = b.iter().flat_map(|x| x.iter().copied()).collect();
    rel_err(&af, &bf)
}

/// Least-squares slope of `ln err` against `p`: the measured decay rate per order.
fn fitted_slope(
    ps: &[usize],
    errs: &[f64],
) -> f64 {
    let n = ps.len() as f64;
    let mx = ps.iter().map(|&p| p as f64).sum::<f64>() / n;
    let my = errs.iter().map(|&e| e.ln()).sum::<f64>() / n;
    let num: f64 = ps
        .iter()
        .zip(errs)
        .map(|(&p, &e)| (p as f64 - mx) * (e.ln() - my))
        .sum();
    let den: f64 = ps.iter().map(|&p| (p as f64 - mx).powi(2)).sum();
    num / den
}

/// Two compact clusters sitting at the centres of two octant boxes of the unit
/// cube; only the cluster in the `(0, 0, 0)` octant carries charge, so the field
/// at the other cluster is exactly the far field of one well-separated source box.
/// The two corner points (zero charge) pin the root box to `[0, 1]^3` so that the
/// octant centres coincide with the child box centres.
fn far_pair() -> (Vec<[f64; 3]>, Vec<f64>, Vec<usize>) {
    let mut rng = Rng::new(32);
    let mut pos = Vec::new();
    let mut src = Vec::new();
    let mut targets = Vec::new();
    for o in 0..8usize {
        let c = [
            if o & 1 == 0 { 0.25 } else { 0.75 },
            if o & 2 == 0 { 0.25 } else { 0.75 },
            if o & 4 == 0 { 0.25 } else { 0.75 },
        ];
        for k in 0..12 {
            let a = 0.05 * f64::from(k as u32);
            pos.push([
                c[0] + 0.05 * (a * 3.1).cos() * (a * 0.7).sin(),
                c[1] + 0.05 * (a * 2.3).sin() * (a * 1.3).cos(),
                c[2] + 0.05 * (a * 1.7).cos(),
            ]);
            // charge only the (0,0,0) octant
            src.push(if o == 0 {
                0.5 + rng.uniform()
            } else {
                0.0
            });
            if o == 7 {
                targets.push(pos.len() - 1);
            }
        }
    }
    for corner in [[0.0, 0.0, 1.0], [1.0, 1.0, 0.0]] {
        pos.push(corner);
        src.push(0.0);
    }
    (pos, src, targets)
}

#[test]
fn expansion_error_decays_geometrically_in_p_at_the_mac_angle() {
    // fixed source distribution and target set (uniform cloud, unit cube)
    let mut rng = Rng::new(11);
    let n = 800;
    let pos: Vec<[f64; 3]> = (0..n)
        .map(|_| [rng.uniform(), rng.uniform(), rng.uniform()])
        .collect();
    let src: Vec<f64> = (0..n).map(|_| rng.normal()).collect();
    let (_, exact) = direct_sum(&pos, &src);
    let ps: Vec<usize> = (2..=12).collect();
    for theta in [0.5, 0.35] {
        let mut errs = Vec::new();
        let mut prev = f64::INFINITY;
        for &p in &ps {
            let r = fmm3d(
                &pos,
                &src,
                &FmmOptions {
                    order: p,
                    leaf_size: 16,
                    theta,
                    eps: 0.0,
                },
            )
            .unwrap();
            let e = rel_err_vec(&r.field, &exact);
            errs.push(e);
            if p <= 10 {
                assert!(
                    e < prev,
                    "theta={theta} p={p}: error {e:.3e} did not improve on {prev:.3e}"
                );
                prev = e;
            }
        }
        // measured decay per unit order, and the per-order error ratio
        let slope = fitted_slope(&ps[..9], &errs[..9]);
        println!("theta={theta}: p  error  ratio");
        for (i, (&p, &e)) in ps.iter().zip(&errs).enumerate() {
            let ratio = if i == 0 {
                f64::NAN
            } else {
                e / errs[i - 1]
            };
            println!("  {p:2}  {e:.3e}  {ratio:.3}");
        }
        println!(
            "  fitted slope d(ln err)/dp = {slope:.3}  (ln theta = {:.3}, ln(theta/2) = {:.3})",
            theta.ln(),
            (theta / 2.0).ln()
        );
        // The MAC guarantees (r_A + r_B) <= theta d, so the truncation decays at
        // least as fast as theta^p and no faster than the (theta/2)^p of two boxes
        // of equal radius sitting on the acceptance boundary.
        assert!(
            slope <= theta.ln() + 0.02,
            "theta={theta}: slope {slope:.3} decays slower than theta^p"
        );
        assert!(
            slope >= (theta / 2.0).ln() - 0.05,
            "theta={theta}: slope {slope:.3} decays faster than any accepted pair allows"
        );
        assert!(
            errs[8] < 1e-5,
            "theta={theta}: p = 10 error {:.3e}",
            errs[8]
        );
    }
}

#[test]
fn m2l_translation_matches_direct_summation_for_one_well_separated_pair() {
    let (pos, src, targets) = far_pair();
    let (pot, fld) = direct_sum(&pos, &src);
    let q: f64 = src.iter().sum();
    let c_a: [f64; 3] = [0.25, 0.25, 0.25];
    let c_b: [f64; 3] = [0.75, 0.75, 0.75];
    let d: [f64; 3] = [c_b[0] - c_a[0], c_b[1] - c_a[1], c_b[2] - c_a[2]];
    let r: f64 = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
    // Geometry of the accepted pair: both boxes have radius 0.05, so each expansion is
    // truncated at the ratio r/d of its own radius to the centre distance. The pair is
    // accepted because r_A + r_B = 0.1 <= theta d = 0.43.
    let ratio = 0.05 / r;
    let mut errs = Vec::new();
    let ps: Vec<usize> = vec![2, 3, 4, 5, 6, 8, 10];
    for &p in &ps {
        let res = fmm3d(
            &pos,
            &src,
            &FmmOptions {
                order: p,
                leaf_size: 12,
                theta: 0.5,
                eps: 0.0,
            },
        )
        .unwrap();
        let sel_f: Vec<[f64; 3]> = targets.iter().map(|&i| res.field[i]).collect();
        let ref_f: Vec<[f64; 3]> = targets.iter().map(|&i| fld[i]).collect();
        let sel_p: Vec<f64> = targets.iter().map(|&i| res.potential[i]).collect();
        let ref_p: Vec<f64> = targets.iter().map(|&i| pot[i]).collect();
        let e = rel_err_vec(&sel_f, &ref_f);
        println!(
            "p={p:2}  field {e:.3e}  potential {:.3e}  (ratio {ratio:.4}^{p} = {:.3e})",
            rel_err(&sel_p, &ref_p),
            ratio.powi(p as i32)
        );
        errs.push(e);
    }
    let slope = fitted_slope(&ps[..6], &errs[..6]);
    println!("fitted slope {slope:.3} (ln r/d = {:.3})", ratio.ln());
    assert!(
        (slope - ratio.ln()).abs() < 0.4,
        "slope {slope:.3} vs ln(r/d) = {:.3}",
        ratio.ln()
    );
    assert!(errs[6] < 1e-8, "p = 10 error {:.3e}", errs[6]);

    // The lower orders are the truncations of the same expansion: at order 0 the local
    // expansion is constant (the monopole, potential q/r, zero field) and at order 1 it is
    // affine, pot(y) - field(y) . (y - c_B) = const with a constant field (the field is the
    // gradient of the local expansion), which no higher order can satisfy.
    let w_of = |i: usize| [pos[i][0] - c_b[0], pos[i][1] - c_b[1], pos[i][2] - c_b[2]];
    let at = |order: usize| {
        fmm3d(
            &pos,
            &src,
            &FmmOptions {
                order,
                leaf_size: 12,
                theta: 0.5,
                eps: 0.0,
            },
        )
        .unwrap()
    };
    let mono = at(0);
    for &i in &targets {
        assert!(
            (mono.potential[i] - q / r).abs() < 2e-5 * (q / r).abs(),
            "order 0: potential {} vs monopole {}",
            mono.potential[i],
            q / r
        );
        for k in 0..3 {
            assert!(
                mono.field[i][k].abs() < 1e-12,
                "order 0: field {} is not zero",
                mono.field[i][k]
            );
        }
    }
    let lin = at(1);
    let affine: Vec<f64> = targets
        .iter()
        .map(|&i| {
            let w = w_of(i);
            lin.potential[i]
                - lin.field[i][0] * w[0]
                - lin.field[i][1] * w[1]
                - lin.field[i][2] * w[2]
        })
        .collect();
    let lo = affine.iter().copied().fold(f64::INFINITY, f64::min);
    let hi = affine.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    assert!(
        hi - lo < 1e-9 * (q / r).abs(),
        "order 1 is not affine: spread {:.3e}",
        hi - lo
    );
    for &i in &targets {
        for k in 0..3 {
            assert!(
                (lin.field[i][k] - lin.field[targets[0]][k]).abs() < 1e-12 * (q / (r * r)).abs(),
                "order 1: field is not constant"
            );
        }
    }
}

#[test]
fn cost_per_particle_is_bounded_from_one_thousand_to_one_hundred_thousand() {
    // Small leaves keep the near region (a fixed number of leaf widths around a
    // leaf) smaller than the domain at every N here, so the tree does the same work
    // per particle at each size. Order 2 keeps the run short; the cost structure is
    // the same at any order (M2L is O(p^6) per box, independent of N).
    let sizes = [1_000usize, 10_000, 100_000];
    let mut per_particle = Vec::new();
    let mut times = Vec::new();
    for &nn in &sizes {
        let mut rng = Rng::new(13);
        let pos: Vec<[f64; 3]> = (0..nn)
            .map(|_| [rng.uniform(), rng.uniform(), rng.uniform()])
            .collect();
        let src: Vec<f64> = (0..nn).map(|_| rng.normal()).collect();
        let opts = FmmOptions {
            order: 2,
            leaf_size: 4,
            theta: 0.5,
            eps: 0.0,
        };
        // min of two runs: the host is shared, so the minimum is the least noisy
        // estimate of the algorithmic cost
        let mut best = f64::INFINITY;
        let mut checksum = 0.0;
        for _ in 0..2 {
            let t0 = std::time::Instant::now();
            let r = fmm3d(&pos, &src, &opts).unwrap();
            best = best.min(t0.elapsed().as_secs_f64());
            checksum = r.field[nn - 1][0];
        }
        assert!(checksum.is_finite());
        times.push(best);
        per_particle.push(best / nn as f64);
        println!(
            "N = {nn:6}  t = {best:8.3} s  t/N = {:7.2} us  (N^2 would be {:8.1} s)",
            per_particle.last().copied().unwrap_or(f64::NAN) * 1e6,
            times[0] * (nn as f64 / sizes[0] as f64).powi(2)
        );
    }
    let total = times[2] / times[0];
    let per_head = per_particle[2] / per_particle[0];
    println!(
        "t(1e5)/t(1e3) = {total:.1} (linear 100, quadratic 10000); cost per particle x{per_head:.2}"
    );
    // A hundredfold larger problem must not cost a hundredfold more per particle:
    // the near field is O(1) per particle and M2L is O(1) per box, so the residual
    // growth is the interaction list and the cache, not the particle count.
    assert!(
        total < 400.0,
        "t(1e5)/t(1e3) = {total:.1}, close to quadratic"
    );
    assert!(
        per_head < 4.0,
        "cost per particle grew {per_head:.2}x from 1e3 to 1e5"
    );
}

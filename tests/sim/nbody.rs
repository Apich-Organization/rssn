//! Barnes-Hut and 2-D FMM against direct summation; gravity with symplectic integration.

use rssn::kernels::numeric_common::NumericError;
use rssn::kernels::random::Rng;
use rssn::sim::nbody::*;
use rssn::sim::symplectic::*;

fn random_points<const D: usize>(
    n: usize,
    seed: u64,
) -> Vec<[f64; D]> {
    let mut rng = Rng::new(seed);
    (0..n)
        .map(|_| {
            let mut p = [0.0; D];
            for c in &mut p {
                *c = rng.normal();
            }
            p
        })
        .collect()
}

fn rel_error<const D: usize>(
    a: &[[f64; D]],
    b: &[[f64; D]],
) -> f64 {
    let mut num = 0.0;
    let mut den = 0.0;
    for (x, y) in a.iter().zip(b) {
        for k in 0..D {
            num += (x[k] - y[k]).powi(2);
            den += y[k].powi(2);
        }
    }
    (num / den).sqrt()
}

#[test]
fn direct_field_matches_the_softened_two_body_formula() {
    let pos = [[0.0, 0.0, 0.0], [3.0, -4.0, 0.0]];
    let src = [2.0, -5.0];
    let field = direct_field(&pos, &src, 12.0);
    let denominator = 13.0_f64.powi(3);
    let expected = [
        [-15.0 / denominator, 20.0 / denominator, 0.0],
        [-6.0 / denominator, 8.0 / denominator, 0.0],
    ];
    assert!(rel_error(&field, &expected) < 1e-14);
    assert_eq!(field, direct_field(&pos, &src, -12.0));

    assert_eq!(direct_field(&[[0.0], [2.0]], &src, 0.0), [[-5.0], [-2.0]]);
    assert_eq!(
        direct_field(&[[0.0; 4], [2.0, 0.0, 0.0, 0.0]], &src, 0.0),
        [[-0.625, 0.0, 0.0, 0.0], [-0.25, 0.0, 0.0, 0.0]]
    );
    let coincident = [[1.0; 3]; 2];
    assert_eq!(direct_field(&coincident, &src, 0.5), [[0.0; 3]; 2]);
    assert!(
        direct_field(&coincident, &src, 0.0)
            .iter()
            .flatten()
            .all(|x| x.is_nan())
    );
}

#[test]
fn direct_field_preserves_representable_contributions_at_extreme_scales() {
    let wide = direct_field(&[[-f64::MAX], [f64::MAX]], &[2.0, -5.0], 0.0);
    assert!((wide[0][0] / -5.0 - 1.0).abs() < 1e-12);
    assert!((wide[1][0] / -2.0 - 1.0).abs() < 1e-12);

    let check = |actual: f64, expected: f64| {
        assert!(
            (actual / expected - 1.0).abs() < 1e-12,
            "{actual:e} vs {expected:e}"
        );
    };
    let far = direct_field(&[[0.0; 3], [1e160, 0.0, 0.0]], &[0.0, 1e300], 0.0);
    check(far[0][0], 1e-20);
    assert_eq!(far[0][1..], [0.0, 0.0]);
    assert_eq!(far[1], [0.0; 3]);
    let near = direct_field(&[[0.0; 3], [0.0, -1e-170, 0.0]], &[0.0, 1e-300], 0.0);
    check(near[0][1], -1e40);
    let softened = direct_field(&[[0.0; 3], [1e-170, 0.0, 0.0]], &[0.0, 1e-300], 1e-170);
    check(softened[0][0], 1e40 / (2.0 * 2.0_f64.sqrt()));
    let subnormal = direct_field(&[[0.0; 2], [1e308, 0.0]], &[0.0, 1.0], 0.0);
    check(subnormal[0][0], 1e-308);

    // The tiny transverse component needs its weight before rounding to the final scale.
    let mixed = direct_field(&[[0.0; 2], [1e100, -1e-200]], &[0.0, 1e300], 0.0);
    check(mixed[0][0], 1e200);
    check(mixed[0][1], -1e-100);

    // A weighted inverse can overflow or underflow even when the last multiplication fits.
    let overflow = direct_field(&[[0.0; 2], [1e-100, 0.0]], &[0.0, 1e120], 0.0);
    check(overflow[0][0], 1e220);
    let underflow = direct_field(&[[0.0; 3], [1e100, 0.0, 0.0]], &[0.0, 1e-30], 0.0);
    check(underflow[0][0], 1e-230);
    let subnormal_radius_cubed = direct_field(&[[0.0; 3], [1e-107, 0.0, 0.0]], &[0.0, 1e-320], 0.0);
    check(subnormal_radius_cubed[0][0], 1e-320_f64 * 1e214);

    // The exponent exceeds the range of the ordinary inverse at this dimension.
    let mut offset = [0.0; 128];
    offset[0] = 256.0;
    let high_dimensional = direct_field(&[[0.0; 128], offset], &[0.0, 256.0], 0.0);
    check(high_dimensional[0][0], 2.0_f64.powi(-1008));
}

#[test]
fn direct_field_into_reuses_storage_and_checks_shapes_before_writing() {
    let pos = [[0.0, 0.0], [1.0, 2.0]];
    let src = [1.0, -1.0];
    let mut out = [[42.0; 2]; 2];
    assert!(matches!(
        direct_field_into(&pos, &src[..1], 0.0, &mut out),
        Err(NumericError::ShapeMismatch(_))
    ));
    assert!(matches!(
        direct_field_into(&pos, &src, 0.0, &mut out[..1]),
        Err(NumericError::ShapeMismatch(_))
    ));
    assert_eq!(out, [[42.0; 2]; 2]);
    for epsilon in [0.5, 1.0] {
        direct_field_into(&pos, &src, epsilon, &mut out).unwrap();
        assert_eq!(out.as_slice(), direct_field(&pos, &src, epsilon));
    }
    assert!(direct_field_into::<2>(&[], &[], 0.0, &mut []).is_ok());
    assert_eq!(direct_field(&pos[..1], &[], 0.0), [[0.0, 0.0]]);
    assert_eq!(direct_field(&[[0.0; 0]; 2], &[], 0.0), [[0.0_f64; 0]; 2]);
}

#[test]
fn direct_field_is_bit_identical_in_single_and_multi_thread_pools() {
    let pos = random_points::<3>(513, 67);
    let mut rng = Rng::new(90);
    let src: Vec<_> = (0..pos.len()).map(|_| rng.normal()).collect();
    let pool = |threads| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap()
    };
    let serial = pool(1).install(|| direct_field(&pos, &src, 0.05));
    let mut parallel = vec![[0.0; 3]; pos.len()];
    pool(4).install(|| direct_field_into(&pos, &src, 0.05, &mut parallel).unwrap());
    for (a, b) in serial.iter().flatten().zip(parallel.iter().flatten()) {
        assert_eq!(a.to_bits(), b.to_bits());
    }
}

fn bh_vs_direct<const D: usize>(
    n: usize,
    signed: bool,
) -> (f64, f64, f64) {
    let pos = random_points::<D>(n, 10 + D as u64);
    let mut rng = Rng::new(77);
    let src: Vec<f64> = (0..n)
        .map(|_| {
            if signed {
                rng.normal()
            } else {
                0.5 + rng.uniform()
            }
        })
        .collect();
    let exact = direct_field(&pos, &src, 0.0);
    let field = |theta: f64| -> Vec<[f64; D]> {
        let t = BhTree::new(&pos, &src, theta, 0.0).unwrap();
        (0..n).map(|i| t.field(&pos[i], i)).collect()
    };
    (
        rel_error(&field(0.0), &exact),
        rel_error(&field(0.3), &exact),
        rel_error(&field(0.8), &exact),
    )
}

#[test]
fn barnes_hut_converges_to_direct_summation_in_one_two_and_three_dimensions() {
    for signed in [false, true] {
        let (e0, e3, e8) = bh_vs_direct::<1>(400, signed);
        assert!(
            e0 < 1e-12 && e3 < 1e-12 && e8 < 1e-12,
            "D=1 signed={signed}: {e0} {e3} {e8}"
        );
        let (e0, e3, e8) = bh_vs_direct::<2>(600, signed);
        assert!(
            e0 < 1e-12 && e3 < 5e-2 && e3 < e8,
            "D=2 signed={signed}: {e0} {e3} {e8}"
        );
        let (e0, e3, e8) = bh_vs_direct::<3>(600, signed);
        assert!(
            e0 < 1e-12 && e3 < 5e-2 && e3 < e8,
            "D=3 signed={signed}: {e0} {e3} {e8}"
        );
    }
}

#[test]
fn one_dimensional_field_is_the_constant_sheet_force() {
    // sources at 0 and 1: the field at 0.5 cancels, at 2 it is the sum of the strengths towards the sources
    let t = BhTree::new(&[[0.0], [1.0]], &[1.0, 3.0], 0.5, 0.0).unwrap();
    assert!((t.field(&[0.5], usize::MAX)[0] - (3.0 - 1.0)).abs() < 1e-12);
    assert!((t.field(&[2.0], usize::MAX)[0] + 4.0).abs() < 1e-12);
}

#[test]
fn coincident_bodies_and_bad_input_are_handled() {
    let pos = [[0.5, 0.5, 0.5], [0.5, 0.5, 0.5], [2.0, 0.0, 0.0]];
    let t = BhTree::new(&pos, &[1.0, 1.0, 1.0], 0.5, 0.01).unwrap();
    let f = t.field(&pos[2], 2);
    assert!(f.iter().all(|c| c.is_finite()));
    assert!(BhTree::<3>::new(&[], &[], 0.5, 0.0).is_err());
    assert!(BhTree::<4>::new(&[[0.0; 4]], &[1.0], 0.5, 0.0).is_err());
}

#[test]
fn gravity_and_coulomb_helpers_have_the_right_signs() {
    let pos = vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
    let a = gravity_accelerations(&pos, &[1.0, 1.0], 1.0, 0.0, 0.5).unwrap();
    assert!((a[0][0] - 0.25).abs() < 1e-12 && (a[1][0] + 0.25).abs() < 1e-12);
    let f = coulomb_forces(&pos, &[1.0, 1.0], 1.0, 0.0, 0.5).unwrap();
    assert!((f[0][0] + 0.25).abs() < 1e-12 && (f[1][0] - 0.25).abs() < 1e-12);
    let f = coulomb_forces(&pos, &[1.0, -1.0], 1.0, 0.0, 0.5).unwrap();
    assert!((f[0][0] - 0.25).abs() < 1e-12);
}

#[test]
fn fmm_2d_matches_direct_summation_and_improves_with_order() {
    let n = 1500;
    let pos = random_points::<2>(n, 5);
    let mut rng = Rng::new(6);
    let src: Vec<f64> = (0..n).map(|_| rng.normal()).collect();
    let exact = direct_field(&pos, &src, 0.0);
    let mut prev = f64::INFINITY;
    for p in [4, 8, 16, 24] {
        let f = fmm2d_field(&pos, &src, p, 20).unwrap();
        let e = rel_error(&f, &exact);
        assert!(e < prev, "p={p}: error {e} did not improve on {prev}");
        prev = e;
    }
    assert!(prev < 1e-9, "p = 24 error {prev}");
    // an FMM and a Barnes-Hut tree agree on the same problem to the BH accuracy
    let t = BhTree::new(&pos, &src, 0.3, 0.0).unwrap();
    let bh: Vec<[f64; 2]> = (0..n).map(|i| t.field(&pos[i], i)).collect();
    let f = fmm2d_field(&pos, &src, 24, 20).unwrap();
    assert!(rel_error(&bh, &f) < 5e-2);
    assert!(fmm2d_field(&[], &[], 8, 20).is_err());
}

#[test]
fn barnes_hut_gravity_with_a_sixth_order_integrator_keeps_a_kepler_orbit() {
    // central mass 1000 with a light test body on a circular orbit of radius 10 (3-D)
    let m: f64 = 1000.0;
    let r: f64 = 10.0;
    let v = (m / r).sqrt();
    let mut q = vec![0.0, 0.0, 0.0, r, 0.0, 0.0];
    let mut p = vec![0.0, 0.0, 0.0, 0.0, v, 0.0];
    let masses = [m, 1e-6];
    let mut acc = |x: &[f64], a: &mut [f64]| {
        let pos = [[x[0], x[1], x[2]], [x[3], x[4], x[5]]];
        let g = gravity_accelerations(&pos, &masses, 1.0, 0.0, 0.5).unwrap();
        for i in 0..2 {
            for k in 0..3 {
                a[3 * i + k] = g[i][k];
            }
        }
    };
    let period = 2.0 * std::f64::consts::PI * r / v;
    let n = 200;
    integrate(
        Scheme::Yoshida6,
        &mut q,
        &mut p,
        period / n as f64,
        n,
        &mut acc,
    );
    let rr = ((q[3] - q[0]).powi(2) + (q[4] - q[1]).powi(2) + (q[5] - q[2]).powi(2)).sqrt();
    assert!((rr - r).abs() < 1e-6, "radius {rr}");
    assert!(
        (q[3] - q[0] - r).abs() < 1e-4,
        "returns after one period: {}",
        q[3] - q[0]
    );
}

fn clustered_points(
    n: usize,
    seed: u64,
) -> Vec<[f64; 3]> {
    // three tight Gaussian clumps of very different size plus a diffuse background
    let mut rng = Rng::new(seed);
    (0..n)
        .map(|i| {
            let (c, s) = match i % 4 {
                | 0 => ([0.0, 0.0, 0.0], 0.02),
                | 1 => ([3.0, 1.0, -2.0], 0.2),
                | 2 => ([-4.0, 2.0, 1.0], 0.005),
                | _ => ([0.0, 0.0, 0.0], 3.0),
            };
            [
                c[0] + s * rng.normal(),
                c[1] + s * rng.normal(),
                c[2] + s * rng.normal(),
            ]
        })
        .collect()
}

#[test]
fn fmm_3d_error_decreases_geometrically_with_order_on_uniform_and_clustered_sets() {
    for clustered in [false, true] {
        let n = 1200;
        let pos = if clustered {
            clustered_points(n, 21)
        } else {
            random_points::<3>(n, 22)
        };
        let mut rng = Rng::new(8);
        let src: Vec<f64> = (0..n).map(|_| rng.normal()).collect();
        let exact = direct_field(&pos, &src, 0.0);
        let mut prev = f64::INFINITY;
        let mut errs = Vec::new();
        for p in [2usize, 4, 6, 8, 10] {
            let f = fmm3d_field(&pos, &src, p, 16).unwrap();
            let e = rel_error(&f, &exact);
            errs.push(e);
            assert!(e < prev, "clustered={clustered} p={p}: {e} vs {prev}");
            prev = e;
        }
        assert!(
            errs[4] < 1e-4 && errs[4] < errs[0] * 1e-2,
            "clustered={clustered}: {errs:?}"
        );
    }
}

#[test]
fn fmm_3d_potential_and_field_match_direct_sums_and_tightening_theta_helps() {
    use rssn::sim::fmm::FmmOptions;
    use rssn::sim::fmm::fmm3d;
    let n = 500;
    let pos = random_points::<3>(n, 31);
    let mut rng = Rng::new(32);
    let src: Vec<f64> = (0..n).map(|_| 0.5 + rng.uniform()).collect();
    let exact = direct_field(&pos, &src, 0.0);
    let mut pot = vec![0.0; n];
    for i in 0..n {
        for j in 0..n {
            if i != j {
                let r2: f64 = (0..3).map(|k| (pos[i][k] - pos[j][k]).powi(2)).sum();
                pot[i] += src[j] / r2.sqrt();
            }
        }
    }
    let mut last = f64::INFINITY;
    for theta in [0.7, 0.5, 0.3] {
        let r = fmm3d(
            &pos,
            &src,
            &FmmOptions {
                order: 8,
                leaf_size: 12,
                theta,
                eps: 0.0,
            },
        )
        .unwrap();
        let ef = rel_error(&r.field, &exact);
        let ep = (r
            .potential
            .iter()
            .zip(&pot)
            .map(|(a, b)| (a - b).powi(2))
            .sum::<f64>()
            / pot.iter().map(|b| b * b).sum::<f64>())
        .sqrt();
        assert!(
            ef < last && ep < 1e-2,
            "theta={theta}: field {ef} potential {ep}"
        );
        last = ef;
    }
    assert!(last < 1e-6, "theta = 0.3 order 8: {last}");
    // a single leaf (huge leaf size) is exact direct summation
    let r = fmm3d(
        &pos,
        &src,
        &FmmOptions {
            order: 3,
            leaf_size: 10_000,
            theta: 0.5,
            eps: 0.0,
        },
    )
    .unwrap();
    assert!(rel_error(&r.field, &exact) < 1e-13);
    assert!(fmm3d(&[], &[], &FmmOptions::default()).is_err());
    assert!(
        fmm3d(
            &pos,
            &src,
            &FmmOptions {
                theta: 1.5,
                ..FmmOptions::default()
            }
        )
        .is_err()
    );
}

#[test]
fn fmm_2d_is_adaptive_on_strongly_clustered_points() {
    let n = 1500;
    let mut rng = Rng::new(41);
    let pos: Vec<[f64; 2]> = (0..n)
        .map(|i| {
            let s = if i % 3 == 0 { 0.001 } else { 1.0 };
            [s * rng.normal(), s * rng.normal()]
        })
        .collect();
    let src: Vec<f64> = (0..n).map(|_| rng.normal()).collect();
    let exact = direct_field(&pos, &src, 0.0);
    let f = fmm2d_field(&pos, &src, 20, 16).unwrap();
    assert!(rel_error(&f, &exact) < 1e-7);
}

//! Barnes-Hut and 2-D FMM against direct summation; gravity with symplectic integration.

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

fn bh_vs_direct<const D: usize>(
    n: usize,
    signed: bool,
) -> (f64, f64, f64) {
    let pos = random_points::<D>(n, 10 + D as u64);
    let mut rng = Rng::new(77);
    let src: Vec<f64> = (0..n).map(|_| if signed { rng.normal() } else { 0.5 + rng.uniform() }).collect();
    let exact = direct_field(&pos, &src, 0.0);
    let field = |theta: f64| -> Vec<[f64; D]> {
        let t = BhTree::new(&pos, &src, theta, 0.0).unwrap();
        (0..n).map(|i| t.field(&pos[i], i)).collect()
    };
    (rel_error(&field(0.0), &exact), rel_error(&field(0.3), &exact), rel_error(&field(0.8), &exact))
}

#[test]
fn barnes_hut_converges_to_direct_summation_in_one_two_and_three_dimensions() {
    for signed in [false, true] {
        let (e0, e3, e8) = bh_vs_direct::<1>(400, signed);
        assert!(e0 < 1e-12 && e3 < 1e-12 && e8 < 1e-12, "D=1 signed={signed}: {e0} {e3} {e8}");
        let (e0, e3, e8) = bh_vs_direct::<2>(600, signed);
        assert!(e0 < 1e-12 && e3 < 5e-2 && e3 < e8, "D=2 signed={signed}: {e0} {e3} {e8}");
        let (e0, e3, e8) = bh_vs_direct::<3>(600, signed);
        assert!(e0 < 1e-12 && e3 < 5e-2 && e3 < e8, "D=3 signed={signed}: {e0} {e3} {e8}");
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
    integrate(Scheme::Yoshida6, &mut q, &mut p, period / n as f64, n, &mut acc);
    let rr = ((q[3] - q[0]).powi(2) + (q[4] - q[1]).powi(2) + (q[5] - q[2]).powi(2)).sqrt();
    assert!((rr - r).abs() < 1e-6, "radius {rr}");
    assert!((q[3] - q[0] - r).abs() < 1e-4, "returns after one period: {}", q[3] - q[0]);
}

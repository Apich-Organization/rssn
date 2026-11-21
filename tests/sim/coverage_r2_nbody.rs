//! Coverage round 2: the 1PN (EIH) helpers of the Hermite integrator, its block-time bookkeeping
//! and the 2-D fast multipole method.

use rssn::kernels::random::Rng;
use rssn::sim::fmm::FmmOptions;
use rssn::sim::fmm::fmm2d;
use rssn::sim::hermite::*;

const G: f64 = 1.0;

fn newtonian_acc(
    mass: &[f64],
    pos: &[[f64; 3]],
) -> Vec<[f64; 3]> {
    let n = mass.len();
    let mut a = vec![[0.0; 3]; n];
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
            let r = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
            for k in 0..3 {
                a[i][k] += G * mass[j] * d[k] / (r * r * r);
            }
        }
    }
    a
}

fn three_body() -> (Vec<f64>, Vec<[f64; 3]>, Vec<[f64; 3]>) {
    (
        vec![1.0, 0.6, 0.2],
        vec![[0.0, 0.0, 0.0], [1.3, 0.2, -0.1], [-0.4, 1.1, 0.5]],
        vec![[0.0, 0.1, 0.0], [0.0, 0.7, 0.1], [0.3, -0.2, 0.4]],
    )
}

fn norm3(a: [f64; 3]) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}

#[test]
fn eih_accelerations_reduce_to_newton_and_the_correction_scales_as_c_minus_two() {
    let (m, x, v) = three_body();
    let an = newtonian_acc(&m, &x);
    let huge = eih_accelerations(&m, &x, &v, &an, G, 1.0e7);
    for i in 0..3 {
        for k in 0..3 {
            assert!(
                (huge[i][k] - an[i][k]).abs() < 1e-10,
                "c -> infinity gives Newton"
            );
        }
    }
    let corr = |c: f64| -> Vec<[f64; 3]> {
        let a = eih_accelerations(&m, &x, &v, &an, G, c);
        (0..3)
            .map(|i| [a[i][0] - an[i][0], a[i][1] - an[i][1], a[i][2] - an[i][2]])
            .collect()
    };
    let (c1, c2) = (corr(40.0), corr(80.0));
    for i in 0..3 {
        let (n1, n2) = (norm3(c1[i]), norm3(c2[i]));
        assert!(n1 > 0.0);
        assert!(
            (n1 / n2 - 4.0).abs() < 0.02,
            "particle {i}: ratio {}",
            n1 / n2
        );
    }
}

#[test]
fn eih_accelerations_respect_symmetry_and_translation_invariance() {
    // two equal masses, symmetric placement: a0 = -a1 (mirror symmetry through the origin)
    let m = vec![1.0, 1.0];
    let x = vec![[0.5, 0.0, 0.0], [-0.5, 0.0, 0.0]];
    let v = vec![[0.0, 0.4, 0.0], [0.0, -0.4, 0.0]];
    let an = newtonian_acc(&m, &x);
    let a = eih_accelerations(&m, &x, &v, &an, G, 20.0);
    for k in 0..3 {
        assert!(
            (a[0][k] + a[1][k]).abs() < 1e-12,
            "momentum balance, component {k}"
        );
    }
    // translation of every position by the same vector changes nothing
    let (m3, x3, v3) = three_body();
    let an3 = newtonian_acc(&m3, &x3);
    let shift = [3.0, -2.0, 7.5];
    let xs: Vec<[f64; 3]> = x3
        .iter()
        .map(|p| [p[0] + shift[0], p[1] + shift[1], p[2] + shift[2]])
        .collect();
    let a0 = eih_accelerations(&m3, &x3, &v3, &an3, G, 25.0);
    let a1 = eih_accelerations(&m3, &xs, &v3, &an3, G, 25.0);
    for i in 0..3 {
        for k in 0..3 {
            assert!((a0[i][k] - a1[i][k]).abs() < 1e-11);
        }
    }
    // a pure rotation about z of positions, velocities and Newtonian accelerations rotates the result
    let (cth, sth) = (0.8_f64.cos(), 0.8_f64.sin());
    let rot = |p: [f64; 3]| [cth * p[0] - sth * p[1], sth * p[0] + cth * p[1], p[2]];
    let xr: Vec<_> = x3.iter().map(|&p| rot(p)).collect();
    let vr: Vec<_> = v3.iter().map(|&p| rot(p)).collect();
    let anr: Vec<_> = an3.iter().map(|&p| rot(p)).collect();
    let ar = eih_accelerations(&m3, &xr, &vr, &anr, G, 25.0);
    for i in 0..3 {
        let want = rot(a0[i]);
        for k in 0..3 {
            assert!((ar[i][k] - want[k]).abs() < 1e-11, "rotation covariance");
        }
    }
}

#[test]
fn eih_energy_has_the_newtonian_limit_and_is_invariant_under_rigid_motions_and_time_reversal() {
    let (m, x, v) = three_body();
    let newton = {
        let mut e = 0.0;
        for i in 0..3 {
            e += 0.5 * m[i] * (v[i][0].powi(2) + v[i][1].powi(2) + v[i][2].powi(2));
            for j in i + 1..3 {
                let d = [x[i][0] - x[j][0], x[i][1] - x[j][1], x[i][2] - x[j][2]];
                e -= G * m[i] * m[j] / norm3(d);
            }
        }
        e
    };
    let e_inf = eih_energy(&m, &x, &v, G, 1.0e8);
    assert!((e_inf - newton).abs() < 1e-12, "{e_inf} vs {newton}");
    // the 1PN term scales as 1/c^2
    let d1 = eih_energy(&m, &x, &v, G, 40.0) - newton;
    let d2 = eih_energy(&m, &x, &v, G, 80.0) - newton;
    assert!(d1.abs() > 1e-9 && (d1 / d2 - 4.0).abs() < 0.02, "{d1} {d2}");
    // translation, rotation and v -> -v leave it unchanged
    let e0 = eih_energy(&m, &x, &v, G, 30.0);
    let xs: Vec<_> = x
        .iter()
        .map(|p| [p[0] + 1.0, p[1] - 4.0, p[2] + 0.25])
        .collect();
    assert!((eih_energy(&m, &xs, &v, G, 30.0) - e0).abs() < 1e-12);
    let (cth, sth) = (1.1_f64.cos(), 1.1_f64.sin());
    let rot = |p: [f64; 3]| [p[0], cth * p[1] - sth * p[2], sth * p[1] + cth * p[2]];
    let xr: Vec<_> = x.iter().map(|&p| rot(p)).collect();
    let vr: Vec<_> = v.iter().map(|&p| rot(p)).collect();
    assert!((eih_energy(&m, &xr, &vr, G, 30.0) - e0).abs() < 1e-12);
    let vneg: Vec<_> = v.iter().map(|p| [-p[0], -p[1], -p[2]]).collect();
    assert!((eih_energy(&m, &x, &vneg, G, 30.0) - e0).abs() < 1e-12);
    // a single particle at rest carries no energy
    assert_eq!(
        eih_energy(&[2.0], &[[1.0, 2.0, 3.0]], &[[0.0; 3]], G, 10.0),
        0.0
    );
}

#[test]
fn next_time_advances_monotonically_with_the_block_steps() {
    let e: f64 = 0.3;
    let vp = ((1.0 + e) / (1.0 - e)).sqrt();
    let cfg = HermiteConfig {
        g: 1.0,
        eps: 0.0,
        order: HermiteOrder::Fourth,
        policy: StepPolicy::Block {
            eta: 0.01,
            dt_max: 1.0 / 16.0,
            dt_min: 1e-9,
        },
        t0: 0.0,
    };
    let mut s = HermiteSystem::new(
        vec![1.0, 1e-3],
        vec![[0.0; 3], [1.0 - e, 0.0, 0.0]],
        vec![[0.0; 3], [0.0, vp, 0.0]],
        &cfg,
    )
    .unwrap();
    let mut last = s.next_time();
    assert!(last > 0.0 && last.is_finite() && last <= 1.0 / 16.0 + 1e-15);
    for _ in 0..200 {
        s.step().unwrap();
        let t = s.next_time();
        assert!(
            t >= last - 1e-15,
            "block end times never run backwards: {t} < {last}"
        );
        assert!(t.is_finite());
        last = t;
    }
    assert!(last > 0.5);
    // the system clock never passes the next block end
    assert!(s.time <= s.next_time() + 1e-15);
}

// -------------------------------------------------------------- 2-D FMM

fn random_points(
    n: usize,
    seed: u64,
) -> (Vec<[f64; 2]>, Vec<f64>) {
    let mut rng = Rng::new(seed);
    let pos: Vec<[f64; 2]> = (0..n).map(|_| [rng.uniform(), rng.uniform()]).collect();
    let src: Vec<f64> = (0..n).map(|_| rng.uniform() - 0.5).collect();
    (pos, src)
}

fn direct2(
    pos: &[[f64; 2]],
    src: &[f64],
    eps: f64,
) -> Vec<[f64; 2]> {
    let n = pos.len();
    let mut out = vec![[0.0; 2]; n];
    for i in 0..n {
        for j in 0..n {
            if i == j {
                continue;
            }
            let d = [pos[j][0] - pos[i][0], pos[j][1] - pos[i][1]];
            let f = src[j] / (d[0] * d[0] + d[1] * d[1] + eps * eps);
            out[i][0] += f * d[0];
            out[i][1] += f * d[1];
        }
    }
    out
}

fn rel_err(
    a: &[[f64; 2]],
    b: &[[f64; 2]],
) -> f64 {
    let num: f64 = a
        .iter()
        .zip(b)
        .map(|(p, q)| (p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2))
        .sum();
    let den: f64 = b.iter().map(|q| q[0] * q[0] + q[1] * q[1]).sum();
    (num / den).sqrt()
}

#[test]
fn fmm2d_matches_the_direct_sum_and_converges_geometrically_in_the_order() {
    let (pos, src) = random_points(400, 3);
    let exact = direct2(&pos, &src, 0.0);
    let mut prev = f64::MAX;
    for order in [2usize, 4, 8, 12] {
        let opts = FmmOptions {
            order,
            leaf_size: 12,
            theta: 0.5,
            eps: 0.0,
        };
        let f = fmm2d(&pos, &src, &opts).unwrap();
        let e = rel_err(&f, &exact);
        assert!(
            e < prev * 0.7 + 1e-14,
            "order {order}: {e} not below {prev}"
        );
        prev = e;
    }
    assert!(prev < 1e-6, "order 12 error {prev}");
}

#[test]
fn fmm2d_is_linear_in_the_sources_and_translation_invariant() {
    let (pos, src) = random_points(150, 8);
    let opts = FmmOptions {
        order: 10,
        leaf_size: 10,
        theta: 0.5,
        eps: 0.0,
    };
    let f1 = fmm2d(&pos, &src, &opts).unwrap();
    let src2: Vec<f64> = src.iter().map(|s| 2.5 * s).collect();
    let f2 = fmm2d(&pos, &src2, &opts).unwrap();
    for (a, b) in f1.iter().zip(&f2) {
        assert!((2.5 * a[0] - b[0]).abs() < 1e-9 * (1.0 + b[0].abs()));
        assert!((2.5 * a[1] - b[1]).abs() < 1e-9 * (1.0 + b[1].abs()));
    }
    let moved: Vec<[f64; 2]> = pos.iter().map(|p| [p[0] + 100.0, p[1] - 40.0]).collect();
    let f3 = fmm2d(&moved, &src, &opts).unwrap();
    assert!(rel_err(&f3, &f1) < 1e-6, "translation invariance");
    // the net force of a pair of opposite unit sources: the field they exert on each other is
    // equal and opposite in direction times the source strength (Newton's third law, in this kernel
    // the forces are -s_j times the unit vector, so the two have equal magnitude)
    let two = fmm2d(
        &[[0.0, 0.0], [1.0, 0.0]],
        &[1.0, 1.0],
        &FmmOptions::default(),
    )
    .unwrap();
    assert!((two[0][0] - 1.0).abs() < 1e-12 && (two[1][0] + 1.0).abs() < 1e-12);
    assert!(two[0][1].abs() < 1e-15);
}

#[test]
fn fmm2d_rejects_invalid_input() {
    let opts = FmmOptions::default();
    assert!(fmm2d(&[], &[], &opts).is_err());
    assert!(fmm2d(&[[0.0, 0.0], [1.0, 1.0]], &[1.0], &opts).is_err());
    for theta in [0.0, 1.0, -0.5, 2.0, f64::NAN] {
        let bad = FmmOptions { theta, ..opts };
        assert!(
            fmm2d(&[[0.0, 0.0], [1.0, 1.0]], &[1.0, 1.0], &bad).is_err(),
            "theta {theta}"
        );
    }
    // a single particle feels nothing
    let one = fmm2d(&[[0.3, 0.4]], &[5.0], &opts).unwrap();
    assert_eq!(one, vec![[0.0, 0.0]]);
}

#[test]
fn fmm2d_is_deterministic_and_handles_coincident_clusters() {
    let (mut pos, src) = random_points(60, 1);
    // many points in one tight cluster stress the adaptive tree depth
    for p in pos.iter_mut().take(30) {
        p[0] = 0.5 + 1e-6 * p[0];
        p[1] = 0.5 + 1e-6 * p[1];
    }
    let opts = FmmOptions {
        order: 8,
        leaf_size: 8,
        theta: 0.5,
        eps: 0.0,
    };
    let a = fmm2d(&pos, &src, &opts).unwrap();
    let b = fmm2d(&pos, &src, &opts).unwrap();
    assert_eq!(a, b);
    assert!(a.iter().all(|f| f[0].is_finite() && f[1].is_finite()));
    assert!(rel_err(&a, &direct2(&pos, &src, 0.0)) < 1e-3);
}

/// The documented kernel is `d / (|d|^2 + eps^2)`. The expansions are those of the
/// unsoftened `1/z`, so a pair uses them only when the gap between the boxes makes the
/// difference smaller than the truncation error; closer pairs are summed directly.
/// (This used to apply eps to the leaf sums only and was 3% off at eps = 0.05.)
#[test]
fn fmm2d_softening_is_consistent_between_near_and_far_field() {
    let (pos, src) = random_points(400, 3);
    let opts = FmmOptions {
        order: 12,
        leaf_size: 12,
        theta: 0.4,
        eps: 0.05,
    };
    let soft = fmm2d(&pos, &src, &opts).unwrap();
    assert!(rel_err(&soft, &direct2(&pos, &src, 0.05)) < 1e-3);
}

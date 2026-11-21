//! Coverage round 2: Kerr helper functions and symplectic helpers that had no direct test.

use std::f64::consts::PI;

use rssn::sim::kerr::*;
use rssn::sim::symplectic::*;

// ------------------------------------------------------------------- Kerr

#[test]
fn photon_orbit_radius_limits_and_defining_equation() {
    // Schwarzschild: the light ring is at 3M for both senses
    assert!((kerr_photon_orbit_radius(1.0, 0.0, true) - 3.0).abs() < 1e-12);
    assert!((kerr_photon_orbit_radius(1.0, 0.0, false) - 3.0).abs() < 1e-12);
    // extremal limits: prograde -> M, retrograde -> 4M
    assert!((kerr_photon_orbit_radius(1.0, 1.0, true) - 1.0).abs() < 1e-9);
    assert!((kerr_photon_orbit_radius(1.0, 1.0, false) - 4.0).abs() < 1e-9);
    // Bardeen: r^{3/2} - 3 M r^{1/2} +- 2 a M^{1/2} = 0 (upper sign: prograde)
    for &m in &[1.0, 2.5] {
        for k in 0..=9 {
            let a = m * k as f64 / 10.0;
            let rp = kerr_photon_orbit_radius(m, a, true);
            let rr = kerr_photon_orbit_radius(m, a, false);
            let f = |r: f64, s: f64| r.powf(1.5) - 3.0 * m * r.sqrt() + s * 2.0 * a * m.sqrt();
            assert!(f(rp, 1.0).abs() < 1e-9 * m.powf(1.5), "prograde a={a}");
            assert!(f(rr, -1.0).abs() < 1e-9 * m.powf(1.5), "retrograde a={a}");
            assert!(rp <= 3.0 * m + 1e-12 && rr >= 3.0 * m - 1e-12);
            // outside the horizon
            assert!(rp > m + (m * m - a * a).sqrt() - 1e-12);
        }
    }
    // mass scaling
    let r1 = kerr_photon_orbit_radius(1.0, 0.7, true);
    let r3 = kerr_photon_orbit_radius(3.0, 2.1, true);
    assert!((r3 - 3.0 * r1).abs() < 1e-12);
}

#[test]
fn shadow_boundary_is_a_circle_for_schwarzschild_and_mirror_symmetric_for_kerr() {
    let n = 40;
    let pts = shadow_boundary_analytic(1.0, 0.0, PI / 2.0, n).unwrap();
    assert_eq!(pts.len(), 2 * n);
    let rad = 27.0_f64.sqrt();
    for (x, y) in &pts {
        assert!(((x * x + y * y).sqrt() - rad).abs() < 1e-12);
    }
    // Kerr, inclined observer: the curve is symmetric under beta -> -beta
    let kerr = shadow_boundary_analytic(1.0, 0.9, 1.2, n).unwrap();
    assert_eq!(kerr.len(), 2 * n);
    for &(x, y) in &kerr {
        assert!(
            kerr.iter()
                .any(|&(x2, y2)| (x2 - x).abs() < 1e-9 && (y2 + y).abs() < 1e-9),
            "mirror of ({x}, {y}) missing"
        );
    }
    // the shadow of a spinning hole is shifted: the centroid in alpha is nonzero and
    // its size stays between the extremal and Schwarzschild bounds
    let cx: f64 = kerr.iter().map(|p| p.0).sum::<f64>() / kerr.len() as f64;
    assert!(cx.abs() > 0.1, "frame dragging shifts the shadow: {cx}");
    let (xmin, xmax) = kerr
        .iter()
        .fold((f64::MAX, f64::MIN), |(a, b), p| (a.min(p.0), b.max(p.0)));
    assert!(xmax - xmin > 4.0 && xmax - xmin < 2.0 * rad + 1.0);
    // invalid input
    assert!(shadow_boundary_analytic(1.0, 1.0, 1.0, 10).is_err());
    assert!(shadow_boundary_analytic(1.0, 0.5, 0.0, 10).is_err());
    assert!(shadow_boundary_analytic(1.0, 0.5, 1.0, 1).is_err());
    // scaling with mass
    let a = shadow_boundary_analytic(1.0, 0.5, 1.0, 12).unwrap();
    let b = shadow_boundary_analytic(2.0, 1.0, 1.0, 12).unwrap();
    for (p, q) in a.iter().zip(&b) {
        assert!((2.0 * p.0 - q.0).abs() < 1e-8 && (2.0 * p.1 - q.1).abs() < 1e-8);
    }
}

#[test]
fn polar_period_is_none_for_equatorial_and_periodic_for_inclined_orbits() {
    let eq = KerrOrbit::from_apsides(1.0, 0.5, (10.0, 25.0), 0.0, true).unwrap();
    assert!(eq.polar_period().is_none());
    let tilted = KerrOrbit::from_apsides(1.0, 0.5, (10.0, 25.0), 0.5, true).unwrap();
    let tp = tilted
        .polar_period()
        .expect("inclined orbit has a polar period");
    assert!(tp > 0.0 && tp.is_finite());
    for &l in &[0.0, 0.37, 1.1] {
        let a = tilted.cos_theta(l);
        let b = tilted.cos_theta(l + tp);
        assert!(
            (a - b).abs() < 1e-9,
            "cos theta must repeat after the Mino period"
        );
        assert!(a.abs() <= 0.5 + 1e-9, "|cos theta| bounded by z_minus");
    }
    // the half period flips the sign of cos(theta)
    assert!((tilted.cos_theta(0.3) + tilted.cos_theta(0.3 + 0.5 * tp)).abs() < 1e-9);
    // radial motion has its own period and spans exactly [r_p, r_a]
    let tr = tilted.radial_period();
    assert!((tilted.radius(0.0) - 10.0).abs() < 1e-8);
    assert!((tilted.radius(0.5 * tr) - 25.0).abs() < 1e-8);
    assert!((tilted.radius(tr) - 10.0).abs() < 1e-8);
    let (_, omegas) = tilted.frequencies();
    assert!(omegas.iter().all(|w| w.is_finite() && *w > 0.0));
}

#[test]
fn unbound_geodesic_impact_parameter_and_constants() {
    let mu2 = 1.0;
    let e = 1.25;
    let l = 6.0;
    let g = UnboundGeodesic::new(
        1.0,
        0.0,
        mu2,
        (e, l, 0.0),
        (500.0, PI / 2.0, 0.0),
        (false, true),
    )
    .unwrap();
    let b = g.impact_parameter().unwrap();
    assert!((b - l / (e * e - mu2).sqrt()).abs() < 1e-12);
    // photon: b = L / E
    let ph = UnboundGeodesic::new(
        1.0,
        0.0,
        0.0,
        (1.0, 7.0, 0.0),
        (500.0, PI / 2.0, 0.0),
        (false, true),
    )
    .unwrap();
    assert!((ph.impact_parameter().unwrap() - 7.0).abs() < 1e-12);
    // parabolic: no impact parameter
    let par = UnboundGeodesic::new(
        1.0,
        0.0,
        1.0,
        (1.0, 5.0, 0.0),
        (500.0, PI / 2.0, 0.0),
        (false, true),
    )
    .unwrap();
    assert!(par.impact_parameter().is_none());
    // equatorial orbit: the measured Carter constant is zero
    assert!(g.carter_constant().abs() < 1e-12);
    // invalid start data
    let bad = |m: f64, a: f64, c: (f64, f64, f64), pos: (f64, f64, f64)| {
        UnboundGeodesic::new(m, a, mu2, c, pos, (false, true)).is_err()
    };
    assert!(bad(1.0, 1.5, (e, l, 0.0), (500.0, 1.0, 0.0)));
    assert!(bad(1.0, 0.0, (0.5, l, 0.0), (500.0, 1.0, 0.0)));
    assert!(bad(1.0, 0.0, (e, l, -1.0), (500.0, 1.0, 0.0)));
    assert!(bad(1.0, 0.0, (e, l, 0.0), (1.9, 1.0, 0.0)));
    assert!(bad(0.0, 0.0, (e, l, 0.0), (500.0, 1.0, 0.0)));
}

#[test]
fn weak_field_deflection_decreases_with_impact_parameter_and_matches_einstein() {
    // light bending 4M/b for large b (Schwarzschild, first order)
    let mut prev = f64::MAX;
    for &b in &[200.0, 400.0] {
        let mut g = UnboundGeodesic::new(
            1.0,
            0.0,
            0.0,
            (1.0, b, 0.0),
            (1.0e5, PI / 2.0, 0.0),
            (false, true),
        )
        .unwrap();
        match g.scatter(1e-10).unwrap() {
            | ScatterOutcome::Escaped { deflection, .. } => {
                let want = 4.0 / b;
                assert!(
                    ((deflection - want) / want).abs() < 0.03,
                    "b={b}: {deflection} vs {want}"
                );
                assert!(deflection < prev);
                prev = deflection;
            },
            | ScatterOutcome::Captured => panic!("must not be captured at b = {b}"),
        }
    }
    // below the critical impact parameter 3 sqrt(3) M the photon is captured
    let mut g = UnboundGeodesic::new(
        1.0,
        0.0,
        0.0,
        (1.0, 4.5, 0.0),
        (1.0e3, PI / 2.0, 0.0),
        (false, true),
    )
    .unwrap();
    assert_eq!(g.scatter(1e-10).unwrap(), ScatterOutcome::Captured);
    // an outward start cannot be scattered
    let mut out = UnboundGeodesic::new(
        1.0,
        0.0,
        0.0,
        (1.0, 4.5, 0.0),
        (1.0e3, PI / 2.0, 0.0),
        (true, true),
    )
    .unwrap();
    assert!(out.scatter(1e-10).is_err());
}

#[test]
fn trace_ray_rejects_forbidden_rays_and_captures_central_ones() {
    // radial ray (xi = eta = 0) from a far observer is captured
    let o = trace_ray(1.0, 0.0, PI / 2.0, 100.0, (0.0, 0.0), true).unwrap();
    assert_eq!(o, RayOutcome::Captured);
    // an impact parameter far outside the shadow escapes
    let o = trace_ray(1.0, 0.0, PI / 2.0, 100.0, (30.0, 0.0), true).unwrap();
    assert!(matches!(o, RayOutcome::Escaped { .. }));
    // eta < 0 is forbidden at the equator
    assert!(trace_ray(1.0, 0.0, PI / 2.0, 100.0, (0.0, -10.0), true).is_err());
}

// -------------------------------------------------------------- symplectic

#[test]
fn composition_weights_are_symmetric_sum_to_one_and_reject_bad_orders() {
    for bad in [0usize, 1, 3, 5, 7] {
        assert!(composition_weights(bad).is_err(), "order {bad}");
    }
    assert_eq!(composition_weights(2).unwrap(), vec![1.0]);
    let mut prev_len = 1;
    for order in [2usize, 4, 6, 8, 10] {
        let w = composition_weights(order).unwrap();
        let sum: f64 = w.iter().sum();
        assert!((sum - 1.0).abs() < 1e-13, "order {order}: sum {sum}");
        for (a, b) in w.iter().zip(w.iter().rev()) {
            assert!((a - b).abs() < 1e-13, "palindrome");
        }
        assert!(w.len() >= prev_len);
        prev_len = w.len();
        // the cubic moment must vanish for order >= 4
        if order >= 4 {
            let s3: f64 = w.iter().map(|x| x.powi(3)).sum();
            assert!(s3.abs() < 1e-12, "order {order}: sum w^3 = {s3}");
        }
    }
}

#[test]
fn composed_ttl_step_of_order_two_equals_the_plain_step_and_odd_order_errs() {
    let mass = [1.0];
    let mut grad = |q: &[f64], g: &mut [f64]| {
        // V(q) = 1 + q^2/2  (positive), g = -dV/dq
        g[0] = -q[0];
        1.0 + 0.5 * q[0] * q[0]
    };
    let mut a = TtlState {
        q: vec![1.0],
        p: vec![0.3],
        t: 0.0,
    };
    let mut b = a.clone();
    ttl_step(&mut a, &mass, 0.5, 0.01, &mut grad);
    ttl_step_composed(&mut b, &mass, 0.5, 0.01, 2, &mut grad).unwrap();
    assert_eq!(a.q, b.q);
    assert_eq!(a.p, b.p);
    assert_eq!(a.t, b.t);
    assert!(ttl_step_composed(&mut b, &mass, 0.5, 0.01, 3, &mut grad).is_err());
    assert!(ttl_step_composed(&mut b, &mass, 0.5, 0.01, 0, &mut grad).is_err());
}

#[test]
fn triple_jump_lifts_the_order_of_verlet_and_keeps_stage_sums() {
    let verlet = Scheme::Verlet.stages();
    let tj = triple_jump(&verlet, 2);
    let (d, k): (f64, f64) = tj.iter().fold((0.0, 0.0), |s, x| (s.0 + x.0, s.1 + x.1));
    assert!((d - 1.0).abs() < 1e-13 && (k - 1.0).abs() < 1e-13);
    // harmonic oscillator error: order 4 beats order 2 and scales as h^4
    let err = |stages: &[(f64, f64)], n: usize| {
        let h = 2.0 / n as f64;
        let (mut q, mut p) = ([1.0], [0.0]);
        for _ in 0..n {
            step(
                stages,
                &mut q,
                &mut p,
                h,
                &mut |q: &[f64], a: &mut [f64]| a[0] = -q[0],
            );
        }
        let t = 2.0_f64;
        ((q[0] - t.cos()).powi(2) + (p[0] + t.sin()).powi(2)).sqrt()
    };
    let e1 = err(&tj, 20);
    let e2 = err(&tj, 40);
    let rate = (e1 / e2).log2();
    assert!((rate - 4.0).abs() < 0.5, "triple jump order {rate}");
    assert!(e2 < err(&verlet, 40));
}

#[test]
fn ks_matrix_rows_are_orthogonal_with_norm_r_and_generate_the_hopf_map() {
    let us: [[f64; 4]; 3] = [
        [0.3, -1.2, 0.8, 0.5],
        [1.0, 0.0, 0.0, 0.0],
        [-0.4, 0.7, 0.1, -2.0],
    ];
    for u in us {
        let l = ks_matrix(&u);
        let r: f64 = u.iter().map(|c| c * c).sum();
        for i in 0..4 {
            for j in 0..4 {
                let dot: f64 = (0..4).map(|k| l[i][k] * l[j][k]).sum();
                let want = if i == j { r } else { 0.0 };
                assert!((dot - want).abs() < 1e-12, "L L^T [{i}][{j}]");
            }
        }
        // x = L u (first three rows) has |x| = |u|^2
        let x: Vec<f64> = (0..3)
            .map(|i| (0..4).map(|k| l[i][k] * u[k]).sum())
            .collect();
        let nx = (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]).sqrt();
        assert!((nx - r).abs() < 1e-12);
        // fourth component of L u vanishes (the bilinear relation)
        let w: f64 = (0..4).map(|k| l[3][k] * u[k]).sum();
        assert!(w.abs() < 1e-12);
    }
}

#[test]
fn ks_cartesian_round_trip_in_both_hemispheres() {
    let cases = [
        ([1.0, 0.3, -0.2], [0.1, 1.0, 0.2]),
        ([-1.0, 0.3, -0.2], [0.1, 1.0, 0.2]),
        ([-2.0, 0.0, 0.0], [0.0, 0.5, 0.5]),
        ([0.0, 1.5, 0.5], [-0.3, 0.0, 0.9]),
        ([1e-3, 0.0, 1e-3], [0.0, 10.0, 0.0]),
    ];
    for (x, v) in cases {
        let (u, up) = ks_from_cartesian(&x, &v);
        let (x2, v2) = ks_to_cartesian(&u, &up);
        for i in 0..3 {
            assert!((x[i] - x2[i]).abs() < 1e-12 * (1.0 + x[i].abs()), "x {x:?}");
            assert!((v[i] - v2[i]).abs() < 1e-9 * (1.0 + v[i].abs()), "v {x:?}");
        }
        // the KS bilinear constraint: row 4 of L(u) dotted with u' vanishes
        let l = ks_matrix(&u);
        let c: f64 = (0..4).map(|k| l[3][k] * up[k]).sum();
        assert!(c.abs() < 1e-9, "KS bilinear constraint {c}");
    }
}

fn no_perturbation(
    _x: &[f64; 3],
    _v: &[f64; 3],
    _t: f64,
) -> [f64; 3] {
    [0.0; 3]
}

#[test]
fn ks_advance_circular_orbit_matches_the_analytic_motion() {
    for &mu in &[1.0, 4.0] {
        let r = 2.0_f64;
        let vc = (mu / r).sqrt();
        let omega = (mu / r.powi(3)).sqrt();
        let period = 2.0 * PI / omega;
        for &frac in &[0.1, 0.25, 0.5, 0.9] {
            let dt = frac * period;
            let (x, v) = ks_advance(
                &[r, 0.0, 0.0],
                &[0.0, vc, 0.0],
                mu,
                &no_perturbation,
                dt,
                0.05,
            )
            .unwrap();
            let th = omega * dt;
            let want_x = [r * th.cos(), r * th.sin(), 0.0];
            let want_v = [-vc * th.sin(), vc * th.cos(), 0.0];
            for i in 0..3 {
                assert!(
                    (x[i] - want_x[i]).abs() < 1e-6,
                    "mu {mu} frac {frac} x[{i}]"
                );
                assert!(
                    (v[i] - want_v[i]).abs() < 1e-6,
                    "mu {mu} frac {frac} v[{i}]"
                );
            }
        }
    }
}

#[test]
fn ks_advance_eccentric_orbit_conserves_integrals_and_is_reversible() {
    let mu = 1.0;
    let x0 = [1.0, 0.0, 0.0];
    let v0 = [0.0, 0.5, 0.2]; // inclined, strongly eccentric
    let e0 = 0.5 * (0.25 + 0.04) - 1.0;
    let l0 = [0.0, -0.2, 0.5];
    let (x, v) = ks_advance(&x0, &v0, mu, &no_perturbation, 7.0, 0.02).unwrap();
    let r = (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]).sqrt();
    let e = 0.5 * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) - mu / r;
    assert!((e - e0).abs() < 1e-8, "energy drift {}", e - e0);
    let l = [
        x[1] * v[2] - x[2] * v[1],
        x[2] * v[0] - x[0] * v[2],
        x[0] * v[1] - x[1] * v[0],
    ];
    for i in 0..3 {
        assert!((l[i] - l0[i]).abs() < 1e-8, "L[{i}]");
    }
    // time reversal: flip the velocity and advance again
    let (xb, vb) = ks_advance(&x, &[-v[0], -v[1], -v[2]], mu, &no_perturbation, 7.0, 0.02).unwrap();
    for i in 0..3 {
        assert!((xb[i] - x0[i]).abs() < 1e-6, "x back [{i}]");
        assert!((vb[i] + v0[i]).abs() < 1e-6, "v back [{i}]");
    }
    // halving eta converges (agreement between two step controls)
    let (x2, _) = ks_advance(&x0, &v0, mu, &no_perturbation, 7.0, 0.01).unwrap();
    for i in 0..3 {
        assert!((x2[i] - x[i]).abs() < 1e-6);
    }
}

#[test]
fn ks_advance_rejects_bad_arguments() {
    let p = &no_perturbation;
    let (x, v) = ([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
    assert!(ks_advance(&x, &v, 1.0, p, 0.0, 0.1).is_err());
    assert!(ks_advance(&x, &v, 1.0, p, -1.0, 0.1).is_err());
    assert!(ks_advance(&x, &v, 0.0, p, 1.0, 0.1).is_err());
    assert!(ks_advance(&x, &v, 1.0, p, 1.0, 0.0).is_err());
    assert!(ks_advance(&[0.0; 3], &v, 1.0, p, 1.0, 0.1).is_err());
}

#[test]
fn matrix_to_quat_inverts_quat_to_matrix_on_all_four_branches() {
    // unit quaternions, scalar-first; they exercise trace>0 and each dominant-diagonal branch
    let raw: [[f64; 4]; 7] = [
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
        [0.1, 0.9, 0.3, -0.2],
        [-0.3, 0.1, 0.95, 0.2],
        [0.2, -0.3, 0.1, -0.9],
    ];
    for q in raw {
        let n = q.iter().map(|c| c * c).sum::<f64>().sqrt();
        let q = [q[0] / n, q[1] / n, q[2] / n, q[3] / n];
        let m = quat_to_matrix(&q);
        for i in 0..3 {
            for j in 0..3 {
                let d: f64 = (0..3).map(|k| m[i][k] * m[j][k]).sum();
                assert!((d - if i == j { 1.0 } else { 0.0 }).abs() < 1e-12);
            }
        }
        let q2 = matrix_to_quat(&m);
        let norm = q2.iter().map(|c| c * c).sum::<f64>().sqrt();
        assert!((norm - 1.0).abs() < 1e-12);
        let dot: f64 = (0..4).map(|k| q[k] * q2[k]).sum();
        assert!(
            (dot.abs() - 1.0).abs() < 1e-12,
            "q and -q are the same rotation: {q:?} vs {q2:?}"
        );
        let m2 = quat_to_matrix(&q2);
        for i in 0..3 {
            for j in 0..3 {
                assert!((m[i][j] - m2[i][j]).abs() < 1e-12);
            }
        }
    }
}

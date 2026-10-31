//! Exact bound Kerr geodesics in Mino time vs direct integration; frequencies; shadows.

use rssn::sim::kerr::*;

/// Mino-time ODE reference (RK4): state `[r, r', z, z', t, phi]`.
fn mino_rk4(
    o: &KerrOrbit,
    theta0: f64,
    z_dir: f64,
    lambda: f64,
    n: usize,
) -> [f64; 6] {
    let (m, a, e, l, q) = (o.m, o.a, o.energy, o.lz, o.carter);
    let beta = a * a * (1.0 - e * e);
    let rhs = |y: &[f64; 6]| -> [f64; 6] {
        let (r, z) = (y[0], y[2]);
        let k = (l - a * e).powi(2) + q;
        // R'(r)/2
        let c4 = e * e - 1.0;
        let c3 = 2.0 * m;
        let c2 = a * a * (e * e - 1.0) - l * l - q;
        let c1 = 2.0 * m * k;
        let rpp = 0.5 * (4.0 * c4 * r.powi(3) + 3.0 * c3 * r * r + 2.0 * c2 * r + c1);
        let zpp = -(q + l * l + beta) * z + 2.0 * beta * z.powi(3);
        let delta = r * r - 2.0 * m * r + a * a;
        let p = e * (r * r + a * a) - a * l;
        let tdot = (r * r + a * a) * p / delta + a * (l - a * e * (1.0 - z * z));
        let phidot = a * p / delta - a * e + l / (1.0 - z * z);
        [y[1], rpp, y[3], zpp, tdot, phidot]
    };
    let z0 = theta0.cos();
    let zp0 = if o.carter == 0.0 {
        0.0
    } else {
        // (dz/dlambda)^2 = Q - (Q + L^2 + beta) z^2 + beta z^4 ; sign: theta increasing => z decreasing
        let v = (o.carter - (o.carter + l * l + beta) * z0 * z0 + beta * z0.powi(4))
            .max(0.0)
            .sqrt();
        z_dir * v
    };
    let mut y = [o.roots[1], 0.0, z0, zp0, 0.0, 0.0];
    let h = lambda / n as f64;
    for _ in 0..n {
        let k1 = rhs(&y);
        let mut t = y;
        for i in 0..6 {
            t[i] = y[i] + 0.5 * h * k1[i];
        }
        let k2 = rhs(&t);
        for i in 0..6 {
            t[i] = y[i] + 0.5 * h * k2[i];
        }
        let k3 = rhs(&t);
        for i in 0..6 {
            t[i] = y[i] + h * k3[i];
        }
        let k4 = rhs(&t);
        for i in 0..6 {
            y[i] += h / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
        }
    }
    y
}

fn check_against_reference(
    o: &KerrOrbit,
    tol: f64,
) {
    let lam_r = o.radial_period();
    let theta0 = std::f64::consts::FRAC_PI_2;
    for &frac in &[0.3, 1.0, 2.7] {
        let lam = frac * lam_r;
        let n = (lam / 2e-4) as usize + 1;
        let y = mino_rk4(o, theta0, -1.0, lam, n);
        let (t, r, th, ph) = o.state(lam);
        assert!((r - y[0]).abs() < tol * r, "r at {frac}: {r} vs {}", y[0]);
        assert!(
            (th.cos() - y[2]).abs() < tol,
            "z at {frac}: {} vs {}",
            th.cos(),
            y[2]
        );
        assert!(
            (t - y[4]).abs() < tol * t.abs().max(1.0),
            "t at {frac}: {t} vs {}",
            y[4]
        );
        assert!(
            (ph - y[5]).abs() < tol * ph.abs().max(1.0),
            "phi at {frac}: {ph} vs {}",
            y[5]
        );
    }
}

#[test]
fn schwarzschild_orbits_match_direct_mino_integration() {
    let eq = KerrOrbit::from_apsides(1.0, 0.0, (8.0, 20.0), 0.0, true).unwrap();
    check_against_reference(&eq, 1e-9);
    let tilted = KerrOrbit::from_apsides(1.0, 0.0, (8.0, 20.0), 0.5, true).unwrap();
    check_against_reference(&tilted, 1e-9);
}

#[test]
fn kerr_equatorial_and_inclined_orbits_match_direct_mino_integration() {
    for &(a, prograde, zm) in &[
        (0.9, true, 0.0),
        (0.9, false, 0.0),
        (0.7, true, 0.55),
        (0.95, false, 0.8),
    ] {
        let o = KerrOrbit::from_apsides(1.0, a, (9.0, 24.0), zm, prograde).unwrap();
        // E < 1 and the roots bracket the requested apsides
        assert!(
            o.energy < 1.0 && (o.roots[1] - 9.0).abs() < 1e-7 && (o.roots[0] - 24.0).abs() < 1e-7
        );
        check_against_reference(&o, 2e-9);
        if zm > 0.0 {
            assert!((o.z_minus - zm).abs() < 1e-8);
        }
    }
}

#[test]
fn frequencies_give_the_schwarzschild_perihelion_advance() {
    // large semi-latus rectum: advance per radial period ~ 6 pi M / p
    let (rp, ra) = (190.0, 310.0);
    let o = KerrOrbit::from_apsides(1.0, 0.0, (rp, ra), 0.0, true).unwrap();
    let (ups, om) = o.frequencies();
    let p = 2.0 * rp * ra / (rp + ra);
    let advance = 2.0 * std::f64::consts::PI * (om[2] / om[0] - 1.0);
    let expected = 6.0 * std::f64::consts::PI / p;
    assert!(
        (advance - expected).abs() < 0.05 * expected,
        "{advance} vs {expected}"
    );
    // Schwarzschild: spherical symmetry makes the polar and azimuthal frequencies equal
    let tilted = KerrOrbit::from_apsides(1.0, 0.0, (rp, ra), 0.6, true).unwrap();
    let (_, om2) = tilted.frequencies();
    assert!(
        (om2[1] - om2[2]).abs() < 1e-10 * om2[2],
        "{} {}",
        om2[1],
        om2[2]
    );
    assert!(ups[0] > 0.0 && ups[3] > 0.0);
    // Lense-Thirring: for a > 0 the polar and azimuthal frequencies differ
    let kerr = KerrOrbit::from_apsides(1.0, 0.9, (30.0, 60.0), 0.6, true).unwrap();
    let (_, om3) = kerr.frequencies();
    assert!(om3[2] - om3[1] > 1e-6);
    // the frequency definition agrees with the exact state: after one radial period
    // the radius is back and phi advanced by Upsilon_phi/Upsilon_r * 2 pi on average
    let lam = 40.0 * kerr.radial_period();
    let (t1, _, _, ph1) = kerr.state(lam);
    let (t2, _, _, ph2) = kerr.state(2.0 * lam);
    let slope_t = (t2 - t1) / lam;
    let slope_ph = (ph2 - ph1) / lam;
    let (u, _) = kerr.frequencies();
    assert!(
        (slope_t - u[3]).abs() < 1e-4 * u[3],
        "Gamma {slope_t} vs {}",
        u[3]
    );
    assert!(
        (slope_ph - u[2]).abs() < 1e-3 * u[2].abs().max(1e-3),
        "Upsilon_phi {slope_ph} vs {}",
        u[2]
    );
}

#[test]
fn kerr_orbit_constructor_validates_its_input() {
    assert!(
        KerrOrbit::new(
            1.0,
            0.5,
            (1.2, 4.0, 0.0),
            (10.0, std::f64::consts::FRAC_PI_2, 0.0),
            (true, true)
        )
        .is_err()
    );
    assert!(KerrOrbit::from_apsides(1.0, 0.5, (20.0, 10.0), 0.0, true).is_err());
    // unbound / plunging: periapsis inside the separatrix
    assert!(KerrOrbit::from_apsides(1.0, 0.0, (1.5, 30.0), 0.0, true).is_err());
}

fn dist_to_polyline(
    p: (f64, f64),
    poly: &[(f64, f64)],
) -> f64 {
    let mut best = f64::INFINITY;
    for w in poly.windows(2) {
        let (a, b) = (w[0], w[1]);
        let (dx, dy) = (b.0 - a.0, b.1 - a.1);
        let l2 = dx * dx + dy * dy;
        let t = if l2 > 0.0 {
            (((p.0 - a.0) * dx + (p.1 - a.1) * dy) / l2).clamp(0.0, 1.0)
        } else {
            0.0
        };
        best = best.min((p.0 - a.0 - t * dx).hypot(p.1 - a.1 - t * dy));
    }
    best
}

#[test]
fn schwarzschild_shadow_is_the_circle_of_radius_three_root_three() {
    let pts = trace_shadow(1.0, 0.0, 1.2, 1e4, 12).unwrap();
    for (x, y) in pts {
        assert!(
            (x.hypot(y) - 27.0_f64.sqrt()).abs() < 1e-6,
            "{}",
            x.hypot(y)
        );
    }
    let an = shadow_boundary_analytic(1.0, 0.0, 1.2, 50).unwrap();
    assert!(
        an.iter()
            .all(|&(x, y)| (x.hypot(y) - 27.0_f64.sqrt()).abs() < 1e-12)
    );
}

#[test]
fn ray_traced_kerr_shadow_matches_the_bardeen_curve() {
    for &(a, th) in &[
        (0.6, 1.0),
        (0.9, std::f64::consts::FRAC_PI_2),
        (0.998, 1.3),
        (-0.8, 0.9),
    ] {
        let analytic = shadow_boundary_analytic(1.0, a, th, 4000).unwrap();
        // closed polyline
        let mut poly = analytic.clone();
        poly.push(analytic[0]);
        let traced = trace_shadow(1.0, a, th, 1e4, 24).unwrap();
        for p in traced {
            let d = dist_to_polyline(p, &poly);
            assert!(
                d < 2e-4,
                "a={a} theta={th}: traced boundary point {p:?} is {d} from the Bardeen curve"
            );
        }
    }
    // the equatorial prograde / retrograde edges follow from the photon-orbit radii
    let an = shadow_boundary_analytic(1.0, 0.9, std::f64::consts::FRAC_PI_2, 2000).unwrap();
    let (xmin, xmax) = an
        .iter()
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), p| {
            (lo.min(p.0), hi.max(p.0))
        });
    assert!(xmin < -2.0 && xmax > 6.0 && (xmax - xmin) > 8.0);
    assert!(shadow_boundary_analytic(1.0, 1.2, 1.0, 10).is_err());
    assert!(trace_ray(1.0, 0.5, 1.0, 100.0, (0.0, -50.0), true).is_err());
}

#[test]
fn escaping_rays_report_a_celestial_direction() {
    use rssn::sim::kerr::RayOutcome;
    // a ray with a large impact parameter is only mildly deflected
    match trace_ray(
        1.0,
        0.0,
        std::f64::consts::FRAC_PI_2,
        1e4,
        (30.0, 0.0),
        true,
    )
    .unwrap()
    {
        | RayOutcome::Escaped { phi, .. } => {
            // total azimuth swing = pi + deflection angle 4M/b + 15 pi M^2/(4 b^2)
            let deflection = phi.abs() - std::f64::consts::PI;
            assert!(
                (deflection - 0.1464).abs() < 0.01,
                "deflection {deflection}"
            );
        },
        | RayOutcome::Captured => panic!("must escape"),
    }
    assert_eq!(
        trace_ray(1.0, 0.0, std::f64::consts::FRAC_PI_2, 1e4, (4.0, 0.0), true).unwrap(),
        RayOutcome::Captured
    );
}

fn wrap(a: f64) -> f64 {
    let mut d = a % (2.0 * std::f64::consts::PI);
    if d > std::f64::consts::PI {
        d -= 2.0 * std::f64::consts::PI;
    }
    if d < -std::f64::consts::PI {
        d += 2.0 * std::f64::consts::PI;
    }
    d
}

#[test]
fn kerr_schild_integration_reproduces_the_exact_mino_orbit() {
    for &(a, zm, prograde) in &[(0.9, 0.0, true), (0.7, 0.55, true)] {
        let o = KerrOrbit::from_apsides(1.0, a, (9.0, 24.0), zm, prograde).unwrap();
        let lam_r = o.radial_period();
        let theta0 = std::f64::consts::FRAC_PI_2;
        let mut g = KerrSchildGeodesic::from_constants(
            1.0,
            a,
            (o.energy, o.lz, o.carter),
            (o.roots[1], theta0, 0.0),
            1.0,
            (true, true),
        )
        .unwrap();
        assert!(
            (g.hamiltonian() + 0.5).abs() < 1e-10,
            "mass shell {}",
            g.hamiltonian()
        );
        let sigma = |lam: f64| {
            let r = o.radius(lam);
            let c = o.cos_theta(lam);
            r * r + a * a * c * c
        };
        let tau_of = |lam: f64| {
            let n = 4000;
            let h = lam / n as f64;
            let mut s = sigma(0.0) + sigma(lam);
            for k in 1..n {
                s += if k % 2 == 1 { 4.0 } else { 2.0 } * sigma(k as f64 * h);
            }
            s * h / 3.0
        };
        let mut tau_prev = 0.0;
        for k in 1..=5 {
            let lam = 0.45 * k as f64 * lam_r;
            let tau = tau_of(lam);
            assert!(g.advance(tau - tau_prev, 1e-12, 1.0));
            tau_prev = tau;
            let (r, th, ph) = g.boyer_lindquist();
            let (_, r_ex, th_ex, ph_ex) = o.state(lam);
            assert!(
                (r - r_ex).abs() < 2e-7 * r_ex,
                "a={a} k={k}: r {r} vs {r_ex}"
            );
            assert!(
                (th.cos() - th_ex.cos()).abs() < 2e-7,
                "a={a} k={k}: cos theta {} vs {}",
                th.cos(),
                th_ex.cos()
            );
            assert!(
                wrap(ph - ph_ex).abs() < 2e-6,
                "a={a} k={k}: phi {ph} vs {ph_ex}"
            );
        }
        assert!((g.energy() - o.energy).abs() < 1e-9 && (g.angular_momentum() - o.lz).abs() < 1e-8);
    }
}

#[test]
fn radial_plunge_crosses_the_schwarzschild_horizon_with_the_cycloid_proper_time() {
    // fall from rest at r0 = 10 M: r = (r0/2)(1 + cos eta), tau = (r0/2) sqrt(r0/2M) (eta + sin eta)
    let r0: f64 = 10.0;
    let e = (1.0 - 2.0 / r0).sqrt();
    let mut g = KerrSchildGeodesic::from_constants(
        1.0,
        0.0,
        (e, 0.0, 0.0),
        (r0, std::f64::consts::FRAC_PI_2, 0.0),
        1.0,
        (false, true),
    )
    .unwrap();
    for &r_target in &[4.0_f64, 2.0, 1.0, 0.5] {
        let eta = (2.0 * r_target / r0 - 1.0).acos();
        let tau = 0.5 * r0 * (r0 / 2.0).sqrt() * (eta + eta.sin());
        let dtau = tau - g.tau;
        assert!(g.advance(dtau, 1e-12, 0.1));
        assert!(
            (g.radius() - r_target).abs() < 2e-8,
            "r = {} vs {r_target}",
            g.radius()
        );
        assert!((g.energy() - e).abs() < 1e-10);
        assert!((g.hamiltonian() + 0.5).abs() < 1e-9);
    }
}

#[test]
fn kerr_plunge_conserves_energy_angular_momentum_and_carter_constant_across_the_horizon() {
    let a: f64 = 0.9;
    let (e, l, q) = (0.95, 0.8, 1.0);
    let mut g =
        KerrSchildGeodesic::from_constants(1.0, a, (e, l, q), (8.0, 1.2, 0.3), 1.0, (false, true))
            .unwrap();
    let rh = 1.0 + (1.0 - a * a).sqrt();
    let mut crossed = false;
    for _ in 0..400 {
        if !g.advance(0.05, 1e-12, 0.9) {
            break;
        }
        if g.radius() < rh {
            crossed = true;
        }
        assert!((g.energy() - e).abs() < 1e-9);
        assert!(
            (g.angular_momentum() - l).abs() < 1e-8,
            "Lz {}",
            g.angular_momentum()
        );
        assert!((g.hamiltonian() + 0.5).abs() < 1e-9);
        assert!(
            (g.carter_constant() - q).abs() < 1e-6,
            "Q {} at r = {}",
            g.carter_constant(),
            g.radius()
        );
    }
    assert!(crossed, "the geodesic must cross the horizon");
}

#[test]
fn polar_geodesic_passes_through_the_axis_without_coordinate_trouble() {
    // Schwarzschild, L_z = 0: the orbit lies in a plane containing the axis and crosses it
    let (e, q) = (0.98, 20.0);
    let mut g = KerrSchildGeodesic::from_constants(
        1.0,
        0.0,
        (e, 0.0, q),
        (12.0, 1.0, 0.0),
        1.0,
        (false, true),
    )
    .unwrap();
    let j0 = (q).sqrt();
    let (mut crossings, mut last_sign) = (0, g.x[1].signum());
    let mut worst_y = 0.0_f64;
    for _ in 0..4000 {
        assert!(g.advance(0.25, 1e-12, 2.5));
        worst_y = worst_y.max(g.x[2].abs());
        let s = g.x[1].signum();
        if s != last_sign {
            crossings += 1;
            last_sign = s;
        }
        // total angular momentum |x cross p|
        let (x, p) = (&g.x, &g.p);
        let lv = [
            x[2] * p[3] - x[3] * p[2],
            x[3] * p[1] - x[1] * p[3],
            x[1] * p[2] - x[2] * p[1],
        ];
        let jt = (lv[0] * lv[0] + lv[1] * lv[1] + lv[2] * lv[2]).sqrt();
        assert!((jt - j0).abs() < 1e-8, "J = {jt}");
        assert!((g.energy() - e).abs() < 1e-9);
    }
    assert!(worst_y < 1e-9, "orbit left its plane: {worst_y}");
    assert!(
        crossings >= 4,
        "the orbit must cross the axis repeatedly: {crossings}"
    );
}

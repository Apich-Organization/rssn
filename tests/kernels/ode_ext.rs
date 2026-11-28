//! Tests for ETDRK4, exp-Rosenbrock, IMEX, symplectic and variational integrators.

use rssn::kernels::dense::Mat;
use rssn::kernels::ode_adaptive::Event;
use rssn::kernels::ode_adaptive::OdeOptions;
use rssn::kernels::ode_ext::ImexOptions;
use rssn::kernels::ode_ext::ImexScheme;
use rssn::kernels::ode_ext::SplitMethod;
use rssn::kernels::ode_ext::VariationalOptions;
use rssn::kernels::ode_ext::etdrk4;
use rssn::kernels::ode_ext::etdrk4_diag;
use rssn::kernels::ode_ext::exprb32;
use rssn::kernels::ode_ext::imex_rk_with;
use rssn::kernels::ode_ext::symplectic_split;
use rssn::kernels::ode_ext::variational_lobatto_with;
use rssn::kernels::ode_sdirk::SdirkOptions;
use rssn::kernels::ode_sdirk::sdirk3;
use rssn::kernels::ode_sdirk::sdirk3_events;
use rssn::kernels::ode_sdirk::sdirk3_fixed;
use rssn::kernels::ode_sdirk::sdirk3_fixed_events;

fn order(
    e1: f64,
    e2: f64,
) -> f64 {
    (e1 / e2).log2()
}

#[test]
fn etdrk4_order_four_and_stiffness() {
    // Bernoulli: y' = lam y + y^2, y = 1 / ((1/y0 + 1/lam) e^{-lam t} - 1/lam)
    let lam = -3.0;
    let exact = |t: f64| 1.0 / ((1.0 / 0.5 + 1.0 / lam) * (-lam * t).exp() - 1.0 / lam);
    let l = Mat::from_rows(&[vec![lam]]).unwrap();
    let err = |n: usize| {
        let tr = etdrk4(&l, |_t, y, o| o[0] = y[0] * y[0], 0.0, 1.0, &[0.5], n).unwrap();
        (tr[n][0] - exact(1.0)).abs()
    };
    let (e1, e2) = (err(8), err(16));
    assert!(order(e1, e2) > 3.7, "{e1} {e2}");
    let diag_err = |n: usize| {
        let tr = etdrk4_diag(&[lam], |_t, y, o| o[0] = y[0] * y[0], 0.0, 1.0, &[0.5], n).unwrap();
        (tr[n][0] - exact(1.0)).abs()
    };
    assert!((diag_err(8) - e1).abs() < 1e-10 && order(diag_err(8), diag_err(16)) > 3.7);
    // stiff linear part: lambda = -1e4, exact solution of y' = lam (y - cos t) - sin t is y = cos t
    let l = Mat::from_rows(&[vec![-1e4]]).unwrap();
    let tr = etdrk4(
        &l,
        |t, _y, o| o[0] = -1e4 * (-t.cos()) - t.sin(),
        0.0,
        2.0,
        &[1.0],
        10,
    )
    .unwrap();
    assert!((tr[10][0] - 2.0f64.cos()).abs() < 1e-5, "{}", tr[10][0]);
    // linear systems are exact (N = 0): rotation
    let l = Mat::from_rows(&[vec![0.0, -2.0], vec![2.0, 0.0]]).unwrap();
    let tr = etdrk4(&l, |_t, _y, o| o.fill(0.0), 0.0, 3.0, &[1.0, 0.0], 3).unwrap();
    assert!((tr[3][0] - 6.0f64.cos()).abs() < 1e-13 && (tr[3][1] - 6.0f64.sin()).abs() < 1e-13);
    // dense and diagonal agree on a diagonal system with a coupled nonlinearity
    let lam = [-1.0, -50.0, -2000.0];
    let ld = Mat::from_rows(&[
        vec![-1.0, 0.0, 0.0],
        vec![0.0, -50.0, 0.0],
        vec![0.0, 0.0, -2000.0],
    ])
    .unwrap();
    let nl = |_t: f64, y: &[f64], o: &mut [f64]| {
        o[0] = y[1] * y[2];
        o[1] = -y[0] * y[2];
        o[2] = 0.3 * y[0] * y[1];
    };
    let a = etdrk4(&ld, nl, 0.0, 1.0, &[1.0, 0.5, 0.2], 40).unwrap();
    let b = etdrk4_diag(&lam, nl, 0.0, 1.0, &[1.0, 0.5, 0.2], 40).unwrap();
    for i in 0..3 {
        assert!((a[40][i] - b[40][i]).abs() < 1e-8);
    }
    assert!(etdrk4(&l, |_t, _y, _o| {}, 0.0, 1.0, &[1.0], 4).is_err());
}

#[test]
fn exponential_rosenbrock_adaptive() {
    // linear problem is integrated exactly with large steps
    let f = |_t: f64, y: &[f64], o: &mut [f64]| {
        o[0] = -100.0 * y[0] + y[1];
        o[1] = -y[0] - 100.0 * y[1];
    };
    let opts = OdeOptions {
        rtol: 1e-8,
        atol: 1e-10,
        ..OdeOptions::default()
    };
    let s = exprb32(f, None, 0.0, 1.0, &[1.0, 0.0], &opts).unwrap();
    let want = (-100.0f64).exp();
    let y = s.y.last().unwrap();
    assert!((y[0] - want * 1.0f64.cos()).abs() < 1e-8 && (y[1] + want * 1.0f64.sin()).abs() < 1e-8);
    assert!(s.t.len() < 12, "linear steps {}", s.t.len());
    // nonlinear stiff: y' = -50 (y - cos t) - sin t has y = cos t (time dependent)
    let g = |t: f64, y: &[f64], o: &mut [f64]| o[0] = -50.0 * (y[0] - t.cos()) - t.sin();
    let s = exprb32(g, None, 0.0, 3.0, &[1.0], &opts).unwrap();
    assert!((s.y.last().unwrap()[0] - 3.0f64.cos()).abs() < 1e-6);
    // van der Pol (mu = 5) against tight-tolerance reference from the same method
    let vdp = |_t: f64, y: &[f64], o: &mut [f64]| {
        o[0] = y[1];
        o[1] = 5.0 * (1.0 - y[0] * y[0]) * y[1] - y[0];
    };
    let tight = OdeOptions {
        rtol: 1e-11,
        atol: 1e-12,
        ..OdeOptions::default()
    };
    let r = exprb32(vdp, None, 0.0, 5.0, &[2.0, 0.0], &tight).unwrap();
    let mut errs = Vec::new();
    for tol in [1e-4, 1e-6] {
        let o = OdeOptions {
            rtol: tol,
            atol: tol * 1e-2,
            ..OdeOptions::default()
        };
        let s = exprb32(vdp, None, 0.0, 5.0, &[2.0, 0.0], &o).unwrap();
        let (a, b) = (s.y.last().unwrap(), r.y.last().unwrap());
        errs.push(((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2)).sqrt());
    }
    assert!(errs[1] < errs[0] && errs[1] < 1e-3, "{errs:?}");
    let jac = |_t: f64, y: &[f64]| {
        Mat::from_rows(&[
            vec![0.0, 1.0],
            vec![-10.0 * y[0] * y[1] - 1.0, 5.0 * (1.0 - y[0] * y[0])],
        ])
        .unwrap()
    };
    let sj = exprb32(
        vdp,
        Some(&jac),
        0.0,
        5.0,
        &[2.0, 0.0],
        &OdeOptions {
            rtol: 1e-8,
            atol: 1e-10,
            ..OdeOptions::default()
        },
    )
    .unwrap();
    let (a, b) = (sj.y.last().unwrap(), r.y.last().unwrap());
    assert!((a[0] - b[0]).abs() < 1e-5);
}

#[test]
fn imex_orders_and_stability() {
    // y' = -k (y - cos t) [implicit] + (-sin t) [explicit]: y = cos t for any k
    // plus a genuinely coupled test: y' = -2 y (implicit) + sin(y) ... compare with fine solution
    let fe = |_t: f64, y: &[f64], o: &mut [f64]| {
        o[0] = y[0].sin();
        o[1] = -y[0];
    };
    let fi = |_t: f64, y: &[f64], o: &mut [f64]| {
        o[0] = -3.0 * y[0] + y[1];
        o[1] = -y[1];
    };
    let reference = imex_rk_with(
        fe,
        fi,
        0.0,
        1.0,
        &[1.0, 0.5],
        &ImexOptions {
            jac_im: None,
            steps: 4096,
            scheme: ImexScheme::Ars343,
        },
    )
    .unwrap();
    let r = reference.last().unwrap().clone();
    for (scheme, lo) in [
        (ImexScheme::Euler, 0.9),
        (ImexScheme::Ars222, 1.9),
        (ImexScheme::Ars343, 2.8),
    ] {
        let err = |n: usize| {
            let t = imex_rk_with(
                fe,
                fi,
                0.0,
                1.0,
                &[1.0, 0.5],
                &ImexOptions {
                    jac_im: None,
                    steps: n,
                    scheme,
                },
            )
            .unwrap();
            let y = t.last().unwrap();
            ((y[0] - r[0]).powi(2) + (y[1] - r[1]).powi(2)).sqrt()
        };
        let (e1, e2) = (err(20), err(40));
        assert!(
            order(e1, e2) > lo,
            "{scheme:?}: order {} ({e1}, {e2})",
            order(e1, e2)
        );
    }
    // stiff implicit part with a large step: stable, tracks the slow manifold y = cos t
    let fe = |t: f64, _y: &[f64], o: &mut [f64]| o[0] = -t.sin();
    let fi = |t: f64, y: &[f64], o: &mut [f64]| o[0] = -1e5 * (y[0] - t.cos());
    for scheme in [ImexScheme::Euler, ImexScheme::Ars222, ImexScheme::Ars343] {
        let t = imex_rk_with(
            fe,
            fi,
            0.0,
            2.0,
            &[1.0],
            &ImexOptions {
                jac_im: None,
                steps: 20,
                scheme,
            },
        )
        .unwrap();
        let tol = if scheme == ImexScheme::Euler {
            0.2
        } else {
            5e-3
        };
        assert!(
            (t[20][0] - 2.0f64.cos()).abs() < tol,
            "{scheme:?}: {}",
            t[20][0]
        );
    }
    assert!(
        imex_rk_with(
            fe,
            fi,
            0.0,
            1.0,
            &[],
            &ImexOptions {
                jac_im: None,
                steps: 4,
                scheme: ImexScheme::Euler
            }
        )
        .is_err()
    );
}

fn pendulum_accel(q: &[f64]) -> Vec<f64> {
    vec![-q[0].sin()]
}

fn energy(
    q: &[f64],
    p: &[f64],
) -> f64 {
    0.5 * p[0] * p[0] - q[0].cos()
}

#[test]
fn symplectic_composition_orders() {
    let (q0, p0) = ([1.2], [0.0]);
    // high-accuracy reference from the 8th order method with many steps
    let (qr, pr) = symplectic_split(
        pendulum_accel,
        &q0,
        &p0,
        0.0,
        6.0,
        400,
        SplitMethod::Yoshida8,
    )
    .unwrap();
    let err = |m: SplitMethod, n: usize| {
        let (q, p) = symplectic_split(pendulum_accel, &q0, &p0, 0.0, 6.0, n, m).unwrap();
        ((q[0] - qr[0]).powi(2) + (p[0] - pr[0]).powi(2)).sqrt()
    };
    for (m, n1, lo) in [
        (SplitMethod::Verlet, 40, 1.9),
        (SplitMethod::TripleJump(4), 20, 3.8),
        (SplitMethod::TripleJump(6), 10, 5.6),
        (SplitMethod::Yoshida6, 10, 5.6),
    ] {
        let (e1, e2) = (err(m, n1), err(m, 2 * n1));
        assert!(
            order(e1, e2) > lo,
            "{m:?}: order {} ({e1}, {e2})",
            order(e1, e2)
        );
    }
    // eighth order: coarse steps keep the error far below round-off sensitive range
    let e8 = |n: usize| {
        let (q, p) =
            symplectic_split(pendulum_accel, &q0, &p0, 0.0, 6.0, n, SplitMethod::Yoshida8).unwrap();
        let (qq, pp) = symplectic_split(
            pendulum_accel,
            &q0,
            &p0,
            0.0,
            6.0,
            2000,
            SplitMethod::TripleJump(6),
        )
        .unwrap();
        ((q[0] - qq[0]).powi(2) + (p[0] - pp[0]).powi(2)).sqrt()
    };
    assert!(
        e8(48) < 1e-8 && order(e8(24), e8(48)) > 7.0,
        "{} {}",
        e8(24),
        e8(48)
    );
    // bounded energy error over a long run
    let (q, p) = symplectic_split(
        pendulum_accel,
        &q0,
        &p0,
        0.0,
        2000.0,
        20000,
        SplitMethod::Yoshida6,
    )
    .unwrap();
    assert!((energy(&q, &p) - energy(&q0, &p0)).abs() < 1e-6);
    assert!(
        symplectic_split(
            pendulum_accel,
            &q0,
            &p0,
            0.0,
            1.0,
            10,
            SplitMethod::TripleJump(5)
        )
        .is_err()
    );
}

#[test]
fn variational_integrators_order_and_invariants() {
    let m = Mat::from_rows(&[vec![1.0]]).unwrap();
    let gv = |q: &[f64]| vec![q[0].sin()];
    let (q0, p0) = ([1.2], [0.0]);
    let (qr, pr) = symplectic_split(
        pendulum_accel,
        &q0,
        &p0,
        0.0,
        4.0,
        400,
        SplitMethod::Yoshida8,
    )
    .unwrap();
    let err = |s: usize, n: usize| {
        let (q, p) = variational_lobatto_with(
            &m,
            gv,
            &q0,
            &p0,
            0.0,
            4.0,
            &VariationalOptions { steps: n, stages: s },
        )
        .unwrap();
        ((q[0] - qr[0]).powi(2) + (p[0] - pr[0]).powi(2)).sqrt()
    };
    for (s, n1, lo) in [(2usize, 40usize, 1.9), (3, 10, 3.7), (4, 5, 5.4)] {
        let (e1, e2) = (err(s, n1), err(s, 2 * n1));
        assert!(
            order(e1, e2) > lo,
            "stages {s}: order {} ({e1}, {e2})",
            order(e1, e2)
        );
    }
    // two-body problem: angular momentum is conserved exactly (to round-off) for any step
    let mass = Mat::identity(2);
    let kepler = |q: &[f64]| {
        let r3 = (q[0] * q[0] + q[1] * q[1]).powf(1.5);
        vec![q[0] / r3, q[1] / r3]
    };
    let (q, p) = variational_lobatto_with(
        &mass,
        kepler,
        &[1.0, 0.0],
        &[0.0, 0.8],
        0.0,
        20.0,
        &VariationalOptions {
            steps: 100,
            stages: 3,
        },
    )
    .unwrap();
    let l0 = 1.0 * 0.8;
    assert!((q[0] * p[1] - q[1] * p[0] - l0).abs() < 1e-10);
    assert!(
        variational_lobatto_with(
            &mass,
            kepler,
            &[1.0, 0.0],
            &[0.0, 0.8],
            0.0,
            1.0,
            &VariationalOptions { steps: 10, stages: 1 }
        )
        .is_err()
    );
}

// -------------------------------------------------------------------------
// SDIRK3 continuous extension (dense output) and event location
// -------------------------------------------------------------------------

/// `y' = y cos(t)`, `y(0) = 1`, exact `y = exp(sin t)`: smooth and non-stiff.
fn cosine_rhs(
    t: f64,
    y: &[f64],
    dy: &mut [f64],
) {
    dy[0] = y[0] * t.cos();
}

/// Jacobian of [`cosine_rhs`].
fn cosine_jac(
    t: f64,
    _y: &[f64],
) -> Mat {
    Mat::from_rows(&[vec![t.cos()]]).unwrap()
}

#[test]
fn sdirk3_dense_output_is_third_order() {
    let exact = |t: f64| t.sin().exp();
    let point = 0.37_f64;
    let mut prev = f64::NAN;
    let mut mean_order = 0.0;
    let mut count = 0;
    for k in 0..4 {
        let steps = 10usize << k;
        let sol =
            sdirk3_fixed_events(cosine_rhs, cosine_jac, 0.0, 1.0, &[1.0], steps, &[]).unwrap();
        let got = sol.eval(point).unwrap()[0];
        let err = (got - exact(point)).abs();
        if k > 0 {
            let p = (prev / err).log2();
            eprintln!("dense steps = {steps}: err = {err:.3e}, observed order {p:.2}");
            assert!(p > 2.8, "dense-output order {p} at {steps} steps");
            mean_order += p;
            count += 1;
        }
        prev = err;
    }
    let mean_order = mean_order / count as f64;
    assert!(mean_order > 2.9, "mean dense-output order {mean_order}");
}

#[test]
fn sdirk3_dense_output_endpoint_exactness_and_range() {
    let sol = sdirk3_fixed_events(cosine_rhs, cosine_jac, 0.0, 1.0, &[1.0], 20, &[]).unwrap();
    // The Hermite interpolant is the identity at both ends of a step.
    assert!((sol.eval(0.0).unwrap()[0] - 1.0).abs() < 1e-14);
    let (end, want_end) = (sol.eval(1.0).unwrap()[0], sol.y.last().unwrap()[0]);
    assert!((end - want_end).abs() < 1e-14);
    // Outside the integrated range there is no extension.
    assert!(sol.eval(-0.1).is_none());
    assert!(sol.eval(1.1).is_none());
}

#[test]
fn sdirk3_adaptive_solution_exposes_dense_output() {
    let opts = SdirkOptions {
        rtol: 1e-9,
        atol: 1e-11,
        ..SdirkOptions::default()
    };
    let sol = sdirk3(cosine_rhs, cosine_jac, 0.0, 1.0, &[1.0], &opts).unwrap();
    for &t in &[0.0_f64, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0] {
        let got = sol.eval(t).unwrap()[0];
        let want = (t.sin()).exp();
        assert!(
            (got - want).abs() < 1e-6,
            "eval at t = {t}: {got} vs {want}"
        );
    }
    assert!(sol.eval(-0.5).is_none());
    assert!(sol.eval(1.5).is_none());
}

#[test]
fn sdirk3_events_find_the_cosine_zero_crossings() {
    // g(t) = cos(t): crossings at pi/2 and 3 pi/2 in [0, 5].
    let events = [Event {
        g: &|t: f64, _y: &[f64]| t.cos(),
        terminal: false,
        direction: 0,
    }];
    let opts = SdirkOptions {
        rtol: 1e-9,
        atol: 1e-11,
        ..SdirkOptions::default()
    };
    let f = |_t: f64, y: &[f64], dy: &mut [f64]| dy[0] = -y[0];
    let jac = |_t: f64, _y: &[f64]| Mat::from_rows(&[vec![-1.0]]).unwrap();
    let sol = sdirk3_events(f, jac, 0.0, 5.0, &[1.0], &opts, &events).unwrap();
    assert_eq!(sol.events.len(), 2, "expected two crossings");
    let want = [std::f64::consts::FRAC_PI_2, 1.5 * std::f64::consts::PI];
    for (hit, want_t) in sol.events.iter().zip(want) {
        assert!(
            (hit.t - want_t).abs() < 1e-8,
            "event at {:.17}, expected {want_t}",
            hit.t
        );
        assert!((hit.y[0] - (-hit.t).exp()).abs() < 1e-4);
    }
    // The run still reaches t1: a non-terminal event does not stop it.
    assert!((sol.t.last().unwrap() - 5.0).abs() < 1e-12);
}

#[test]
fn sdirk3_terminal_event_stops_at_the_ground() {
    // Projectile: y'' = -g, starting at height 10 with zero velocity. The
    // ground y = 0 is reached at sqrt(2 h / g).
    let g = 9.81_f64;
    let f = |_t: f64, y: &[f64], dy: &mut [f64]| {
        dy[0] = y[1];
        dy[1] = -g;
    };
    let jac = |_t: f64, _y: &[f64]| Mat::from_rows(&[vec![0.0, 1.0], vec![0.0, 0.0]]).unwrap();
    let ground = [Event {
        g: &|_t: f64, y: &[f64]| y[0],
        terminal: true,
        direction: -1,
    }];
    let opts = SdirkOptions {
        rtol: 1e-9,
        atol: 1e-11,
        ..SdirkOptions::default()
    };
    let sol = sdirk3_events(f, jac, 0.0, 2.0, &[10.0, 0.0], &opts, &ground).unwrap();
    let want_t = (2.0 * 10.0 / g).sqrt();
    assert_eq!(sol.events.len(), 1);
    let hit = &sol.events[0];
    assert!(
        (hit.t - want_t).abs() < 1e-6,
        "impact at {:.17}, expected {want_t}",
        hit.t
    );
    assert!(hit.y[0].abs() < 1e-6, "height at impact {}", hit.y[0]);
    // A terminal event ends the trajectory at the event point.
    let last = sol.t.last().unwrap();
    assert!((last - hit.t).abs() < 1e-14);
    assert!((last - 2.0).abs() > 0.1);
}

#[test]
fn sdirk3_fixed_events_matches_fixed_without_events() {
    let fixed = sdirk3_fixed(cosine_rhs, cosine_jac, 0.0, 1.0, &[1.0], 32).unwrap();
    let no_events = sdirk3_fixed_events(cosine_rhs, cosine_jac, 0.0, 1.0, &[1.0], 32, &[]).unwrap();
    assert_eq!(no_events.y.last().unwrap(), &fixed);
    assert_eq!(no_events.t.len(), 33);

    // The adaptive event entry point with an empty list is the plain solver.
    let opts = SdirkOptions::default();
    let plain = sdirk3(cosine_rhs, cosine_jac, 0.0, 1.0, &[1.0], &opts).unwrap();
    let via_events = sdirk3_events(cosine_rhs, cosine_jac, 0.0, 1.0, &[1.0], &opts, &[]).unwrap();
    assert_eq!(plain.t, via_events.t);
    assert_eq!(plain.y, via_events.y);
    assert_eq!(plain.nfev, via_events.nfev);
    assert_eq!(plain.rejected, via_events.rejected);
    assert!((plain.eval(0.5).unwrap()[0] - (0.5_f64.sin()).exp()).abs() < 1e-6);
}

#[test]
fn sdirk3_fixed_events_find_crossings_on_fixed_grid() {
    let events = [Event {
        g: &|t: f64, _y: &[f64]| t.cos(),
        terminal: false,
        direction: 0,
    }];
    let f = |_t: f64, y: &[f64], dy: &mut [f64]| dy[0] = -y[0];
    let jac = |_t: f64, _y: &[f64]| Mat::from_rows(&[vec![-1.0]]).unwrap();
    let sol = sdirk3_fixed_events(f, jac, 0.0, 5.0, &[1.0], 200, &events).unwrap();
    assert_eq!(sol.events.len(), 2);
    let want = [std::f64::consts::FRAC_PI_2, 1.5 * std::f64::consts::PI];
    for (hit, want_t) in sol.events.iter().zip(want) {
        assert!(
            (hit.t - want_t).abs() < 1e-8,
            "fixed-grid event at {}",
            hit.t
        );
    }
}

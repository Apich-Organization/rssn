//! Tests for ETDRK4, exp-Rosenbrock, IMEX, symplectic and variational integrators.

use rssn::kernels::dense::Mat;
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

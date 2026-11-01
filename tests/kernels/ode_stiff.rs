//! Tests for the Radau IIA(5) and variable-order BDF stiff integrators.

use rssn::kernels::dense::Mat;
use rssn::kernels::ode_adaptive::OdeOptions;
use rssn::kernels::ode_adaptive::rosenbrock23;
use rssn::kernels::ode_stiff::bdf;
use rssn::kernels::ode_stiff::BdfOptions;
use rssn::kernels::ode_stiff::bdf_jac;
use rssn::kernels::ode_stiff::bdf_with;
use rssn::kernels::ode_stiff::radau5;
use rssn::kernels::ode_stiff::radau5_fixed;
use rssn::kernels::ode_stiff::radau5_jac;

fn robertson(
    _t: f64,
    y: &[f64],
    dy: &mut [f64],
) {
    dy[0] = -0.04 * y[0] + 1e4 * y[1] * y[2];
    dy[1] = 0.04 * y[0] - 1e4 * y[1] * y[2] - 3e7 * y[1] * y[1];
    dy[2] = 3e7 * y[1] * y[1];
}

fn robertson_jac(
    _t: f64,
    y: &[f64],
) -> Mat {
    Mat::from_rows(&[
        vec![-0.04, 1e4 * y[2], 1e4 * y[1]],
        vec![0.04, -1e4 * y[2] - 6e7 * y[1], -1e4 * y[1]],
        vec![0.0, 6e7 * y[1], 0.0],
    ])
    .unwrap()
}

fn vdp(
    _t: f64,
    y: &[f64],
    dy: &mut [f64],
) {
    dy[0] = y[1];
    dy[1] = 1000.0 * (1.0 - y[0] * y[0]) * y[1] - y[0];
}

// Hairer-Wanner reference values at t = 40.
const ROB40: [f64; 3] = [0.7158270687, 9.185534764e-6, 0.2841637457];

fn check_rob(
    y: &[f64],
    rel: f64,
) {
    for i in 0..3 {
        assert!(
            (y[i] - ROB40[i]).abs() <= rel * ROB40[i],
            "y{i} = {} vs {}",
            y[i],
            ROB40[i]
        );
    }
    assert!((y.iter().sum::<f64>() - 1.0).abs() < 1e-8);
}

#[test]
fn radau5_robertson() {
    let opts = OdeOptions {
        rtol: 1e-9,
        atol: 1e-13,
        ..OdeOptions::default()
    };
    let sol = radau5(robertson, 0.0, 40.0, &[1.0, 0.0, 0.0], &opts).unwrap();
    check_rob(sol.y.last().unwrap(), 1e-5);
    let sol = radau5_jac(robertson, robertson_jac, 0.0, 40.0, &[1.0, 0.0, 0.0], &opts).unwrap();
    check_rob(sol.y.last().unwrap(), 1e-5);
    assert!(sol.t.len() < 1500, "too many steps: {}", sol.t.len());
}

#[test]
fn bdf_robertson() {
    let opts = OdeOptions {
        rtol: 1e-9,
        atol: 1e-13,
        ..OdeOptions::default()
    };
    let sol = bdf(robertson, 0.0, 40.0, &[1.0, 0.0, 0.0], &opts, 5).unwrap();
    check_rob(sol.y.last().unwrap(), 1e-4);
    let sol = bdf_jac(
        robertson,
        robertson_jac,
        0.0,
        40.0,
        &[1.0, 0.0, 0.0],
        &opts,
        5,
    )
    .unwrap();
    check_rob(sol.y.last().unwrap(), 1e-4);
}

#[test]
fn robertson_long_time() {
    let opts = OdeOptions {
        rtol: 1e-8,
        atol: 1e-12,
        ..OdeOptions::default()
    };
    let a = radau5(robertson, 0.0, 1e5, &[1.0, 0.0, 0.0], &opts).unwrap();
    let b = bdf(robertson, 0.0, 1e5, &[1.0, 0.0, 0.0], &opts, 5).unwrap();
    let (ya, yb) = (a.y.last().unwrap(), b.y.last().unwrap());
    assert!(
        (ya[0] - yb[0]).abs() < 1e-4 * ya[0].abs().max(1e-3),
        "{} {}",
        ya[0],
        yb[0]
    );
    assert!((ya.iter().sum::<f64>() - 1.0).abs() < 1e-8);
    assert!(ya[0] > 0.0 && ya[0] < 0.5);
}

#[test]
fn van_der_pol_mu_1000() {
    let opts = OdeOptions {
        rtol: 1e-8,
        atol: 1e-10,
        ..OdeOptions::default()
    };
    let tf = 600.0;
    let r = radau5(vdp, 0.0, tf, &[2.0, 0.0], &opts).unwrap();
    let b = bdf(vdp, 0.0, tf, &[2.0, 0.0], &opts, 5).unwrap();
    let ro = rosenbrock23(
        vdp,
        0.0,
        tf,
        &[2.0, 0.0],
        &OdeOptions {
            rtol: 1e-9,
            atol: 1e-11,
            ..opts
        },
    )
    .unwrap();
    let (yr, yb, yo) = (
        r.y.last().unwrap(),
        b.y.last().unwrap(),
        ro.y.last().unwrap(),
    );
    assert!(
        (yr[0] - yo[0]).abs() < 2e-4,
        "radau {} rosen {}",
        yr[0],
        yo[0]
    );
    assert!(
        (yb[0] - yo[0]).abs() < 2e-3,
        "bdf {} rosen {}",
        yb[0],
        yo[0]
    );
    // slow manifold: t = mu (ln y - y^2/2 - ln 2 + 2); solve for y at tf
    let mut y: f64 = 1.5;
    for _ in 0..60 {
        let g = 1000.0 * (y.ln() - y * y / 2.0 - 2.0_f64.ln() + 2.0) - tf;
        let dg = 1000.0 * (1.0 / y - y);
        y -= g / dg;
    }
    assert!((yr[0] - y).abs() < 0.01, "radau {} manifold {}", yr[0], y);
    // a full period (about 1614) crosses the fast jump to the lower branch
    let long = radau5(vdp, 0.0, 2000.0, &[2.0, 0.0], &opts).unwrap();
    assert!(long.y.iter().any(|v| v[0] < -1.9));
}

#[test]
fn radau5_fixed_step_order_five() {
    // y' = -2y + sin t, y(0) = 1: y = (2 sin t - cos t)/5 + (6/5) e^{-2t}
    let f = |t: f64, y: &[f64], dy: &mut [f64]| dy[0] = -2.0 * y[0] + t.sin();
    let exact = |t: f64| (2.0 * t.sin() - t.cos()) / 5.0 + 1.2 * (-2.0 * t).exp();
    let e = |n: usize| (radau5_fixed(f, 0.0, 2.0, &[1.0], n).unwrap()[0] - exact(2.0)).abs();
    let (e1, e2, e3) = (e(4), e(8), e(16));
    let (o1, o2) = ((e1 / e2).log2(), (e2 / e3).log2());
    assert!(o1 > 4.6 && o1 < 6.5, "order {o1}");
    assert!(o2 > 4.6 && o2 < 6.5, "order {o2}");
}

#[test]
fn adaptive_error_decreases_with_tolerance() {
    let c = 1.0_f64.cos();
    let f = move |_t: f64, y: &[f64], dy: &mut [f64]| dy[0] = -50.0 * (y[0] - c);
    let exact = c + (2.0 - c) * (-50.0_f64).exp();
    let mut last_r = f64::INFINITY;
    let mut last_b = f64::INFINITY;
    for k in [4, 6, 8, 10] {
        let o = OdeOptions {
            rtol: 10f64.powi(-k),
            atol: 10f64.powi(-k - 2),
            ..OdeOptions::default()
        };
        let er = (radau5(f, 0.0, 1.0, &[2.0], &o).unwrap().y.last().unwrap()[0] - exact).abs();
        let eb = (bdf(f, 0.0, 1.0, &[2.0], &o, 5).unwrap().y.last().unwrap()[0] - exact).abs();
        assert!(er < last_r.max(1e-12) * 1.5);
        assert!(eb < last_b.max(1e-12) * 1.5);
        last_r = er;
        last_b = eb;
    }
    assert!(last_r < 1e-9, "radau {last_r}");
    assert!(last_b < 1e-8, "bdf {last_b}");
}

#[test]
fn bdf_order_matters() {
    let f = |_t: f64, y: &[f64], dy: &mut [f64]| {
        dy[0] = y[1];
        dy[1] = -y[0] - 0.1 * y[1];
    };
    let o = OdeOptions {
        rtol: 1e-8,
        atol: 1e-10,
        ..OdeOptions::default()
    };
    let w = (1.0f64 - 0.0025).sqrt();
    let exact = |t: f64| (-0.05 * t).exp() * ((w * t).cos() + 0.05 / w * (w * t).sin());
    let e1 = (bdf(f, 0.0, 5.0, &[1.0, 0.0], &o, 1)
        .unwrap()
        .y
        .last()
        .unwrap()[0]
        - exact(5.0))
    .abs();
    let e5 = (bdf(f, 0.0, 5.0, &[1.0, 0.0], &o, 5)
        .unwrap()
        .y
        .last()
        .unwrap()[0]
        - exact(5.0))
    .abs();
    assert!(e5 < 1e-5, "bdf5 {e5}");
    assert!(e5 * 20.0 < e1, "e1 {e1} e5 {e5}");
}

#[test]
fn backward_integration_and_errors() {
    let f = |_t: f64, y: &[f64], dy: &mut [f64]| dy[0] = -y[0];
    let o = OdeOptions::default();
    let s = radau5(f, 1.0, 0.0, &[(-1.0f64).exp()], &o).unwrap();
    assert!((s.y.last().unwrap()[0] - 1.0).abs() < 1e-6);
    let s = bdf(f, 1.0, 0.0, &[(-1.0f64).exp()], &o, 5).unwrap();
    assert!((s.y.last().unwrap()[0] - 1.0).abs() < 1e-5);
    assert!(radau5(f, 0.0, 0.0, &[1.0], &o).is_err());
    assert!(bdf(f, 0.0, 1.0, &[], &o, 5).is_err());
}

#[test]
fn bdf_accepts_a_time_dependent_mass_matrix() {
    // `(1 + t) y' = -y` with `y(0) = 1` has the exact solution `1/(1 + t)`.
    // The corrector is `M(t_{n+1}) (y - y_pred) = h gamma0 f`, so no `dM/dt` is
    // needed: the mass is simply evaluated at the new time of each step.
    let f = |_t: f64, y: &[f64], dy: &mut [f64]| {
        dy[0] = -y[0];
    };
    let mass_at = |t: f64| Mat::from_rows(&[vec![1.0 + t]]).unwrap_or_else(|e| panic!("{e}"));
    let options = BdfOptions {
        max_order: 5,
        mass_fn: Some(&mass_at),
        ..BdfOptions::default()
    };
    let opts = OdeOptions {
        rtol: 1e-10,
        atol: 1e-12,
        ..OdeOptions::default()
    };
    let solution = bdf_with(f, 0.0, 2.0, &[1.0], &opts, &options)
        .unwrap_or_else(|e| panic!("{e:?}"));
    let last = solution.t.last().copied().unwrap_or(f64::NAN);
    assert!((last - 2.0).abs() < 1e-12);
    let y = solution.y.last().and_then(|row| row.first()).copied().unwrap_or(f64::NAN);
    let exact = 1.0 / (1.0 + last);
    assert!((y - exact).abs() < 1e-8, "{y} vs {exact}");

    // The same problem written with a constant mass and the `1/(1+t)` folded
    // into the right-hand side must agree: the two paths describe one system.
    let folded = |t: f64, y: &[f64], dy: &mut [f64]| {
        dy[0] = -y[0] / (1.0 + t);
    };
    let constant = BdfOptions {
        max_order: 5,
        ..BdfOptions::default()
    };
    let reference = bdf_with(folded, 0.0, 2.0, &[1.0], &opts, &constant)
        .unwrap_or_else(|e| panic!("{e:?}"));
    let reference_y = reference
        .y
        .last()
        .and_then(|row| row.first())
        .copied()
        .unwrap_or(f64::NAN);
    assert!(
        (y - reference_y).abs() < 1e-7,
        "{y} vs {reference_y}"
    );
}

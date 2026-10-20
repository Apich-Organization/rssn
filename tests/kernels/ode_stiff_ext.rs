//! Tests for dense output, events, variable order and DAEs in the stiff integrators.

use rssn::kernels::dense::Mat;
use rssn::kernels::ode_adaptive::{Event, OdeOptions, rosenbrock23};
use rssn::kernels::ode_stiff::{RadauConfig, bdf_ex, radau, radau_fixed};

fn decay(_t: f64, y: &[f64], dy: &mut [f64]) {
    dy[0] = -y[0];
}

fn osc(_t: f64, y: &[f64], dy: &mut [f64]) {
    dy[0] = y[1];
    dy[1] = -y[0];
}

#[test]
fn radau_stage_orders() {
    let f = |t: f64, y: &[f64], dy: &mut [f64]| dy[0] = y[0] * t.cos();
    let exact = 1.0_f64.sin().exp();
    let err = |n: usize, s: usize| (radau_fixed(f, 0.0, 1.0, &[1.0], n, s).unwrap()[0] - exact).abs();
    for (s, lo, n1, n2) in [(3usize, 4.5, 2usize, 4usize), (5, 8.0, 1, 2)] {
        let (e1, e2) = (err(n1, s), err(n2, s));
        let p = (e1 / e2).log2();
        assert!(p > lo, "stages {s}: observed order {p} (errors {e1}, {e2})");
    }
    // order 13: a single step is already accurate to ~1e-11
    assert!(err(1, 7) < 1e-10);
    assert!(err(2, 7) < 1e-12);
}

#[test]
fn radau_dense_output_and_variable_order() {
    let opts = OdeOptions { rtol: 1e-10, atol: 1e-12, ..OdeOptions::default() };
    let sol = radau(decay, None, 0.0, 5.0, &[1.0], &opts, &[], &RadauConfig::default()).unwrap();
    for k in 0..=100 {
        let t = 5.0 * f64::from(k) / 100.0;
        let v = sol.eval(t).unwrap()[0];
        assert!((v - (-t).exp()).abs() < 1e-8, "t={t}: {v}");
    }
    assert!(sol.eval(5.1).is_none());
    assert!(sol.orders.iter().any(|&o| o >= 9), "orders {:?}", sol.orders);
    let fixed = RadauConfig { min_stages: 3, max_stages: 3, ..RadauConfig::default() };
    let s5 = radau(decay, None, 0.0, 5.0, &[1.0], &opts, &[], &fixed).unwrap();
    assert!(s5.orders.iter().all(|&o| o == 5));
    assert!(sol.t.len() < s5.t.len(), "{} vs {}", sol.t.len(), s5.t.len());
    // backward integration
    let back =
        radau(decay, None, 1.0, 0.0, &[(-1.0f64).exp()], &opts, &[], &RadauConfig::default()).unwrap();
    assert!((back.eval(0.3).unwrap()[0] - (-0.3f64).exp()).abs() < 1e-8);
}

#[test]
fn radau_events() {
    let opts = OdeOptions { rtol: 1e-10, atol: 1e-12, ..OdeOptions::default() };
    let g = |_t: f64, y: &[f64]| y[0] - 0.5;
    let ev = [Event { g: &g, terminal: true, direction: -1 }];
    let sol = radau(decay, None, 0.0, 5.0, &[1.0], &opts, &ev, &RadauConfig::default()).unwrap();
    assert_eq!(sol.events.len(), 1);
    assert!((sol.events[0].t - 2.0_f64.ln()).abs() < 1e-8);
    assert!((sol.t.last().unwrap() - 2.0_f64.ln()).abs() < 1e-8);
    // wrong direction: no event, run completes
    let ev = [Event { g: &g, terminal: true, direction: 1 }];
    let sol = radau(decay, None, 0.0, 5.0, &[1.0], &opts, &ev, &RadauConfig::default()).unwrap();
    assert!(sol.events.is_empty());
    assert!((sol.t.last().unwrap() - 5.0).abs() < 1e-12);
    // oscillator zeros of cos: pi/2 (falling), 3pi/2 (rising), 5pi/2 (falling)
    let g0 = |_t: f64, y: &[f64]| y[0];
    let ev = [Event { g: &g0, terminal: false, direction: 0 }];
    let sol = radau(osc, None, 0.0, 10.0, &[1.0, 0.0], &opts, &ev, &RadauConfig::default()).unwrap();
    assert_eq!(sol.events.len(), 3);
    for (k, e) in sol.events.iter().enumerate() {
        let want = std::f64::consts::FRAC_PI_2 + std::f64::consts::PI * k as f64;
        assert!((e.t - want).abs() < 1e-6, "{} vs {want}", e.t);
    }
    let rising = [Event { g: &g0, terminal: true, direction: 1 }];
    let sol =
        radau(osc, None, 0.0, 10.0, &[1.0, 0.0], &opts, &rising, &RadauConfig::default()).unwrap();
    assert!((sol.events[0].t - 1.5 * std::f64::consts::PI).abs() < 1e-6);
}

#[test]
fn bdf_dense_and_events() {
    let opts = OdeOptions { rtol: 1e-9, atol: 1e-11, ..OdeOptions::default() };
    let sol = bdf_ex(decay, None, 0.0, 5.0, &[1.0], &opts, 5, &[]).unwrap();
    for k in 0..=50 {
        let t = 5.0 * f64::from(k) / 50.0;
        let v = sol.eval(t).unwrap()[0];
        assert!((v - (-t).exp()).abs() < 2e-6, "t={t}: {v}");
    }
    let g = |_t: f64, y: &[f64]| y[0] - 0.25;
    let ev = [Event { g: &g, terminal: true, direction: -1 }];
    let sol = bdf_ex(decay, None, 0.0, 5.0, &[1.0], &opts, 5, &ev).unwrap();
    assert_eq!(sol.events.len(), 1);
    assert!((sol.events[0].t - 4.0_f64.ln()).abs() < 1e-5, "{}", sol.events[0].t);
    assert!((sol.t.last().unwrap() - 4.0_f64.ln()).abs() < 1e-5);
    let g0 = |_t: f64, y: &[f64]| y[0];
    let ev = [Event { g: &g0, terminal: false, direction: 1 }];
    let sol = bdf_ex(osc, None, 0.0, 8.0, &[1.0, 0.0], &opts, 5, &ev).unwrap();
    assert_eq!(sol.events.len(), 1);
    assert!((sol.events[0].t - 1.5 * std::f64::consts::PI).abs() < 1e-4);
}

#[test]
fn radau_dae_mass_matrix() {
    // y' = -y + z, 0 = z - y^2  (so y' = y^2 - y, y = 1 / (1 + e^t))
    let f = |_t: f64, u: &[f64], du: &mut [f64]| {
        du[0] = -u[0] + u[1];
        du[1] = u[1] - u[0] * u[0];
    };
    let m = Mat::from_rows(&[vec![1.0, 0.0], vec![0.0, 0.0]]).unwrap();
    let opts = OdeOptions { rtol: 1e-9, atol: 1e-11, ..OdeOptions::default() };
    let cfg = RadauConfig { mass: Some(&m), ..RadauConfig::default() };
    let sol = radau(f, None, 0.0, 4.0, &[0.5, 0.25], &opts, &[], &cfg).unwrap();
    let last = sol.y.last().unwrap();
    let exact = 1.0 / (1.0 + 4.0_f64.exp());
    assert!((last[0] - exact).abs() < 1e-7, "{} vs {exact}", last[0]);
    assert!((last[1] - exact * exact).abs() < 1e-7);
    // mass matrix of the wrong size
    let bad = Mat::identity(3);
    let cfg = RadauConfig { mass: Some(&bad), ..RadauConfig::default() };
    assert!(radau(f, None, 0.0, 1.0, &[0.5, 0.25], &opts, &[], &cfg).is_err());
}

#[test]
fn radau_robertson_dae() {
    // Index-1 form of Robertson's problem: the third equation is the conservation law.
    let f = |_t: f64, y: &[f64], dy: &mut [f64]| {
        dy[0] = -0.04 * y[0] + 1e4 * y[1] * y[2];
        dy[1] = 0.04 * y[0] - 1e4 * y[1] * y[2] - 3e7 * y[1] * y[1];
        dy[2] = y[0] + y[1] + y[2] - 1.0;
    };
    let m = Mat::from_rows(&[vec![1.0, 0.0, 0.0], vec![0.0, 1.0, 0.0], vec![0.0, 0.0, 0.0]]).unwrap();
    let opts = OdeOptions { rtol: 1e-8, atol: 1e-12, ..OdeOptions::default() };
    let cfg = RadauConfig { mass: Some(&m), ..RadauConfig::default() };
    let sol = radau(f, None, 0.0, 40.0, &[1.0, 0.0, 0.0], &opts, &[], &cfg).unwrap();
    let y = sol.y.last().unwrap();
    let want = [0.7158270687, 9.185534764e-6, 0.2841637457];
    for i in 0..3 {
        assert!((y[i] - want[i]).abs() <= 1e-5 * want[i].max(1e-3), "y{i} = {}", y[i]);
    }
}

#[test]
fn rosenbrock_dense_output() {
    let opts = OdeOptions { rtol: 1e-7, atol: 1e-9, ..OdeOptions::default() };
    let sol = rosenbrock23(osc, 0.0, 6.0, &[1.0, 0.0], &opts).unwrap();
    for k in 0..=60 {
        let t = 0.1 * f64::from(k);
        let v = sol.eval(t).unwrap();
        assert!((v[0] - t.cos()).abs() < 2e-4 && (v[1] + t.sin()).abs() < 2e-4, "t={t}");
    }
}

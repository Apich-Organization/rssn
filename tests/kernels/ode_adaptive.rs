//! Tests for the adaptive / stiff / symplectic / BVP ODE kernels.

use rssn::kernels::ode_adaptive::Event;
use rssn::kernels::ode_adaptive::OdeOptions;
use rssn::kernels::ode_adaptive::bvp_fd;
use rssn::kernels::ode_adaptive::dae_index1;
use rssn::kernels::ode_adaptive::dde_rk4;
use rssn::kernels::ode_adaptive::dopri5;
use rssn::kernels::ode_adaptive::leapfrog;
use rssn::kernels::ode_adaptive::rosenbrock23;
use rssn::kernels::ode_adaptive::rosenbrock23_events;
use rssn::kernels::ode_adaptive::shooting;
use rssn::kernels::ode_adaptive::yoshida4;
use rssn::kernels::ode_stiff::bdf_events;
use rssn::kernels::ode_stiff::radau5_events;

fn decay(
    _t: f64,
    y: &[f64],
    dy: &mut [f64],
) {
    dy[0] = -y[0];
}

fn oscillator(
    _t: f64,
    y: &[f64],
    dy: &mut [f64],
) {
    dy[0] = y[1];
    dy[1] = -y[0];
}

#[test]
fn dopri5_accuracy_and_dense_output() {
    let opts = OdeOptions {
        rtol: 1e-10,
        atol: 1e-12,
        ..OdeOptions::default()
    };
    let sol = dopri5(oscillator, 0.0, 10.0, &[0.0, 1.0], &opts, &[]).unwrap();
    let y = sol.y.last().unwrap();
    assert!((y[0] - 10.0_f64.sin()).abs() < 1e-8);
    assert!((y[1] - 10.0_f64.cos()).abs() < 1e-8);
    // dense output between grid points
    for t in [0.37, 2.51, 7.77, 9.99] {
        let v = sol.eval(t).unwrap();
        assert!((v[0] - f64::sin(t)).abs() < 1e-7, "t={t}");
    }
    assert!(sol.eval(11.0).is_none());
}

#[test]
fn dopri5_backward_integration() {
    let sol = dopri5(
        decay,
        1.0,
        0.0,
        &[(-1.0_f64).exp()],
        &OdeOptions::default(),
        &[],
    )
    .unwrap();
    assert!((sol.y.last().unwrap()[0] - 1.0).abs() < 1e-7);
}

#[test]
fn dopri5_events() {
    // ball dropped from h=10 with g=9.81: hits ground at sqrt(2*10/9.81)
    let fall = |_t: f64, y: &[f64], dy: &mut [f64]| {
        dy[0] = y[1];
        dy[1] = -9.81;
    };
    let ground = |_t: f64, y: &[f64]| y[0];
    let events = [Event {
        g: &ground,
        terminal: true,
        direction: -1,
    }];
    let sol = dopri5(
        fall,
        0.0,
        10.0,
        &[10.0, 0.0],
        &OdeOptions::default(),
        &events,
    )
    .unwrap();
    let expect = (20.0_f64 / 9.81).sqrt();
    assert_eq!(sol.events.len(), 1);
    assert!((sol.events[0].t - expect).abs() < 1e-8);
    assert!((sol.t.last().unwrap() - expect).abs() < 1e-8);
    // non-terminal zero crossings of cos: pi/2 + k pi in [0, 10]
    let g = |_t: f64, y: &[f64]| y[1];
    let ev = [Event {
        g: &g,
        terminal: false,
        direction: 0,
    }];
    let sol = dopri5(
        oscillator,
        0.0,
        10.0,
        &[0.0, 1.0],
        &OdeOptions::default(),
        &ev,
    )
    .unwrap();
    assert_eq!(sol.events.len(), 3);
    assert!((sol.events[0].t - std::f64::consts::FRAC_PI_2).abs() < 1e-8);
}

#[test]
fn rosenbrock_stiff() {
    // y' = -1000 (y - cos t) - sin t, solution y = cos t
    let f = |t: f64, y: &[f64], dy: &mut [f64]| {
        dy[0] = -100_000.0 * (y[0] - t.cos()) - t.sin();
    };
    let opts = OdeOptions {
        rtol: 1e-6,
        atol: 1e-9,
        ..OdeOptions::default()
    };
    let sol = rosenbrock23(f, 0.0, 2.0, &[1.0], &opts).unwrap();
    assert!((sol.y.last().unwrap()[0] - 2.0_f64.cos()).abs() < 1e-4);
    // a stiff solver needs far fewer steps than the explicit one on Van der Pol-like decay
    let explicit = dopri5(f, 0.0, 2.0, &[1.0], &opts, &[]).unwrap();
    assert!(sol.t.len() < explicit.t.len());
}

#[test]
fn symplectic_energy_conservation() {
    let acc = |q: &[f64]| vec![-q[0]];
    let (q, p) = leapfrog(acc, &[1.0], &[0.0], 0.01, 100_000);
    let e = 0.5 * (q[0] * q[0] + p[0] * p[0]);
    assert!((e - 0.5).abs() < 1e-4);
    let (q4, p4) = yoshida4(acc, &[1.0], &[0.0], 0.1, 1000);
    assert!((q4[0] - 100.0_f64.cos()).abs() < 2e-3);
    assert!((p4[0] + 100.0_f64.sin()).abs() < 2e-3);
}

#[test]
fn bvp_finite_difference_and_shooting() {
    // y'' = y, y(0)=0, y(1)=1 -> sinh(x)/sinh(1)
    let (x, y) = bvp_fd(|_x, y, _dy| y, 0.0, 1.0, 0.0, 1.0, 199, |s| s).unwrap();
    for (xi, yi) in x.iter().zip(&y).step_by(20) {
        assert!((yi - xi.sinh() / 1.0_f64.sinh()).abs() < 1e-4);
    }
    // same by shooting on slope
    let f = |_t: f64, y: &[f64], dy: &mut [f64]| {
        dy[0] = y[1];
        dy[1] = y[0];
    };
    let p = shooting(
        f,
        0.0,
        1.0,
        |p| vec![0.0, p[0]],
        |e| vec![e[0] - 1.0],
        &[0.5],
        &OdeOptions::default(),
    )
    .unwrap();
    assert!((p[0] - 1.0 / 1.0_f64.sinh()).abs() < 1e-7);
}

#[test]
fn dae_index1_conservation() {
    // y' = -z, 0 = z - y  -> y = exp(-t)
    let traj = dae_index1(
        |_t, _y, z| vec![-z[0]],
        |_t, y, z| vec![z[0] - y[0]],
        &[1.0],
        &[1.0],
        0.0,
        0.001,
        1000,
    )
    .unwrap();
    let last = traj.last().unwrap();
    assert!((last[0] - (-1.0_f64).exp()).abs() < 5e-4);
    assert!((last[1] - last[0]).abs() < 1e-10);
}

#[test]
fn dde_delayed_decay() {
    // y' = -y(t-1), history 1: on [0,1] y = 1 - t; on [1,2] y = 1 - t + (t-1)^2/2
    let (ts, ys) = dde_rk4(
        |_t, _y, yd| vec![-yd[0]],
        |_t| vec![1.0],
        1.0,
        0.0,
        2.0,
        0.01,
    )
    .unwrap();
    let at = |t: f64| {
        ts.iter()
            .position(|&x| (x - t).abs() < 1e-9)
            .map(|i| ys[i][0])
            .unwrap()
    };
    assert!((at(1.0) - 0.0).abs() < 1e-6);
    assert!((at(2.0) - (1.0 - 2.0 + 0.5)).abs() < 1e-5);
    assert!(
        dde_rk4(
            |_t, _y, yd| vec![-yd[0]],
            |_t| vec![1.0],
            0.001,
            0.0,
            1.0,
            0.01
        )
        .is_err()
    );
}

#[test]
fn every_integrator_with_a_dense_output_locates_events() {
    // `y'' = -y` from `(1, 0)` is `y = cos(t)`, so `y[0]` crosses zero at
    // pi/2 and 3 pi/2 in [0, 5]. Each integrator locates the crossing on its
    // own continuous extension, so the times agree to the method's order.
    let osc = |_t: f64, y: &[f64], dy: &mut [f64]| {
        dy[0] = y[1];
        dy[1] = -y[0];
    };
    let zero = |_t: f64, y: &[f64]| y[0];
    let events = [Event {
        g: &zero,
        terminal: false,
        direction: 0,
    }];
    let opts = OdeOptions {
        rtol: 1e-10,
        atol: 1e-12,
        ..OdeOptions::default()
    };
    let want = [
        std::f64::consts::FRAC_PI_2,
        3.0 * std::f64::consts::FRAC_PI_2,
    ];
    let check = |name: &str, hits: &[f64], tolerance: f64| {
        assert_eq!(hits.len(), 2, "{name}: {} events", hits.len());
        for (got, expected) in hits.iter().zip(&want) {
            assert!(
                (got - expected).abs() < tolerance,
                "{name}: {got} vs {expected}"
            );
        }
    };
    let times = |sol: &rssn::kernels::ode_adaptive::OdeSolution| -> Vec<f64> {
        sol.events.iter().map(|hit| hit.t).collect()
    };

    let sol = dopri5(osc, 0.0, 5.0, &[1.0, 0.0], &opts, &events)
        .unwrap_or_else(|e| panic!("dopri5: {e:?}"));
    check("dopri5", &times(&sol), 1e-8);

    // Rosenbrock23's extension is the third-order Hermite interpolant, so its
    // event time is good to the interpolation error, not to the step tolerance.
    let sol = rosenbrock23_events(osc, 0.0, 5.0, &[1.0, 0.0], &opts, &events)
        .unwrap_or_else(|e| panic!("rosenbrock23: {e:?}"));
    check("rosenbrock23", &times(&sol), 1e-6);

    let sol = radau5_events(osc, 0.0, 5.0, &[1.0, 0.0], &opts, &events)
        .unwrap_or_else(|e| panic!("radau5: {e:?}"));
    check("radau5", &times(&sol), 1e-8);

    let sol = bdf_events(osc, 0.0, 5.0, &[1.0, 0.0], &opts, 5, &events)
        .unwrap_or_else(|e| panic!("bdf: {e:?}"));
    check("bdf", &times(&sol), 1e-8);
}

#[test]
fn a_terminal_event_ends_every_integrator_at_the_event() {
    // A ball dropped from rest hits the ground at `sqrt(2 h / g)`; a terminal
    // event must stop the integration there, not at the end of the interval.
    let g = 9.81_f64;
    let fall = move |_t: f64, y: &[f64], dy: &mut [f64]| {
        dy[0] = y[1];
        dy[1] = -g;
    };
    let ground = |_t: f64, y: &[f64]| y[0];
    let events = [Event {
        g: &ground,
        terminal: true,
        direction: -1,
    }];
    let opts = OdeOptions {
        rtol: 1e-10,
        atol: 1e-12,
        ..OdeOptions::default()
    };
    let expected = (2.0 * 10.0 / g).sqrt();
    for (name, sol) in [
        (
            "dopri5",
            dopri5(fall, 0.0, 5.0, &[10.0, 0.0], &opts, &events),
        ),
        (
            "rosenbrock23",
            rosenbrock23_events(fall, 0.0, 5.0, &[10.0, 0.0], &opts, &events),
        ),
        (
            "radau5",
            radau5_events(fall, 0.0, 5.0, &[10.0, 0.0], &opts, &events),
        ),
        (
            "bdf",
            bdf_events(fall, 0.0, 5.0, &[10.0, 0.0], &opts, 5, &events),
        ),
    ] {
        let sol = sol.unwrap_or_else(|e| panic!("{name}: {e:?}"));
        let last = sol.t.last().copied().unwrap_or(f64::NAN);
        assert!(
            (last - expected).abs() < 1e-6,
            "{name}: stopped at {last}, want {expected}"
        );
        assert_eq!(sol.events.len(), 1, "{name}");
        assert!(
            events
                .get(sol.events[0].index)
                .is_some_and(|event| event.terminal),
            "{name}: the hit is not on a terminal event"
        );
    }
}

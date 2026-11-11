//! Robustness of the numerical drivers on adversarial problems: NaN and
//! infinite right-hand sides, finite-time blow-up, degenerate intervals and
//! tolerances, empty states, degenerate polynomials and objectives. A
//! solver may return an error or any value, but must not panic and must
//! return within a bound. Each case runs on its own thread with a watchdog.

use std::sync::mpsc;
use std::time::Duration;

use rssn::kernels::ode::OdeSolverMethod;
use rssn::kernels::ode::solve_adaptive;
use rssn::kernels::ode::solve_fixed;
use rssn::kernels::ode_adaptive::OdeOptions;
use rssn::kernels::ode_adaptive::dopri5;
use rssn::kernels::ode_adaptive::rosenbrock23;
use rssn::kernels::ode_stiff::bdf;
use rssn::kernels::ode_stiff::radau5;
use rssn::kernels::optim;
use rssn::kernels::rootfind;

/// Runs `body` on a fresh thread; reports a panic or a hang as a failure
/// string. Panic messages are suppressed (they are the expected failure
/// mode being tested).
fn guarded<F: FnOnce() + Send + 'static>(
    name: &'static str,
    body: F,
) -> Option<String> {
    let (tx, rx) = mpsc::channel();
    let handle = std::thread::Builder::new()
        .name(name.into())
        .stack_size(32 * 1024 * 1024)
        .spawn(move || {
            let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(body));
            let _ = tx.send(result.is_ok());
        })
        .unwrap();
    match rx.recv_timeout(Duration::from_secs(30)) {
        | Ok(true) => {
            let _ = handle.join();
            None
        },
        | Ok(false) => {
            let _ = handle.join();
            Some(format!("{name}: panicked"))
        },
        | Err(_) => Some(format!("{name}: no result within 30 s")),
    }
}

const NAN: f64 = f64::NAN;
const INF: f64 = f64::INFINITY;

type Rhs = fn(f64, &[f64], &mut [f64]);

fn rhs_nan(
    _: f64,
    _: &[f64],
    d: &mut [f64],
) {
    d.fill(NAN);
}
fn rhs_inf(
    _: f64,
    _: &[f64],
    d: &mut [f64],
) {
    d.fill(INF);
}
fn rhs_blowup(
    _: f64,
    y: &[f64],
    d: &mut [f64],
) {
    for (di, yi) in d.iter_mut().zip(y) {
        *di = yi * yi;
    }
}
fn rhs_stiff(
    _: f64,
    y: &[f64],
    d: &mut [f64],
) {
    for (di, yi) in d.iter_mut().zip(y) {
        *di = -1e9 * yi;
    }
}
fn rhs_discontinuous(
    t: f64,
    y: &[f64],
    d: &mut [f64],
) {
    for (di, yi) in d.iter_mut().zip(y) {
        *di = if t < 0.5 {
            1.0
        } else {
            -yi.signum() * 1e3
        };
    }
}
fn rhs_zero(
    _: f64,
    _: &[f64],
    d: &mut [f64],
) {
    d.fill(0.0);
}
fn rhs_overflow(
    _: f64,
    y: &[f64],
    d: &mut [f64],
) {
    for (di, yi) in d.iter_mut().zip(y) {
        *di = yi * 1e300;
    }
}

fn problems() -> Vec<(&'static str, Rhs, Vec<f64>, (f64, f64))> {
    vec![
        ("nan rhs", rhs_nan, vec![1.0], (0.0, 1.0)),
        ("inf rhs", rhs_inf, vec![1.0], (0.0, 1.0)),
        ("blow-up", rhs_blowup, vec![1.0], (0.0, 2.0)),
        ("stiff", rhs_stiff, vec![1.0, 2.0], (0.0, 1.0)),
        ("discontinuous", rhs_discontinuous, vec![0.0], (0.0, 1.0)),
        ("nan y0", rhs_zero, vec![NAN], (0.0, 1.0)),
        ("inf y0", rhs_zero, vec![INF], (0.0, 1.0)),
        ("empty state", rhs_zero, vec![], (0.0, 1.0)),
        ("empty interval", rhs_zero, vec![1.0], (1.0, 1.0)),
        ("backward interval", rhs_zero, vec![1.0], (1.0, 0.0)),
        ("nan interval", rhs_zero, vec![1.0], (NAN, 1.0)),
        ("infinite interval", rhs_zero, vec![1.0], (0.0, INF)),
        ("huge y0", rhs_overflow, vec![1e300], (0.0, 1.0)),
        ("tiny interval", rhs_zero, vec![1.0], (0.0, 5e-324)),
    ]
}

fn options() -> Vec<OdeOptions> {
    vec![
        OdeOptions {
            max_steps: 5_000,
            ..OdeOptions::default()
        },
        OdeOptions {
            rtol: 0.0,
            atol: 0.0,
            max_steps: 2_000,
            ..OdeOptions::default()
        },
        OdeOptions {
            rtol: NAN,
            atol: -1.0,
            max_steps: 2_000,
            ..OdeOptions::default()
        },
        OdeOptions {
            h0: Some(0.0),
            hmax: 0.0,
            max_steps: 2_000,
            ..OdeOptions::default()
        },
        OdeOptions {
            h0: Some(-1.0),
            rtol: INF,
            max_steps: 2_000,
            ..OdeOptions::default()
        },
        OdeOptions {
            max_steps: 0,
            ..OdeOptions::default()
        },
    ]
}

#[test]
fn ode_solvers_survive_adversarial_problems() {
    let mut failures = Vec::new();
    for (pname, f, y0, span) in problems() {
        for method in [
            OdeSolverMethod::Euler,
            OdeSolverMethod::Heun,
            OdeSolverMethod::RungeKutta4,
        ] {
            for steps in [0_usize, 1, 7, 200] {
                let y0 = y0.clone();
                failures.extend(guarded("solve_fixed", move || {
                    let _ = solve_fixed(f, &y0, span, steps, method);
                }));
            }
        }
        for (rtol, atol, max) in [
            (1e-8, 1e-10, 3000_usize),
            (0.0, 0.0, 1000),
            (NAN, 1.0, 10),
            (1e-8, 1e-10, 0),
        ] {
            let y0 = y0.clone();
            failures.extend(guarded("solve_adaptive", move || {
                let _ = solve_adaptive(f, &y0, span, rtol, atol, max);
            }));
        }
        let y0_outer = y0.clone();
        for opts in options() {
            let y0 = y0_outer.clone();
            failures.extend(guarded("dopri5", move || {
                let _ = dopri5(f, span.0, span.1, &y0, &opts, &[]);
            }));
            let y0 = y0_outer.clone();
            failures.extend(guarded("rosenbrock23", move || {
                let _ = rosenbrock23(f, span.0, span.1, &y0, &opts);
            }));
            let y0 = y0_outer.clone();
            failures.extend(guarded("radau5", move || {
                let _ = radau5(f, span.0, span.1, &y0, &opts);
            }));
            for order in [0_usize, 1, 5, 99] {
                let y0 = y0_outer.clone();
                failures.extend(guarded("bdf", move || {
                    let _ = bdf(f, span.0, span.1, &y0, &opts, order);
                }));
            }
        }
        let _ = pname;
    }
    assert!(
        failures.is_empty(),
        "{} failures:\n{}",
        failures.len(),
        failures.join("\n")
    );
}

#[test]
fn stiff_problem_is_solved_accurately_by_the_stiff_methods() {
    // y' = -1e9 y over [0, 1] must decay to ~0 with an implicit method in
    // a modest number of steps; an explicit method would need 1e9.
    let opts = OdeOptions {
        max_steps: 20_000,
        ..OdeOptions::default()
    };
    let r = radau5(rhs_stiff, 0.0, 1.0, &[1.0, 2.0], &opts).expect("radau5");
    let last = r.y.last().unwrap();
    assert!(last.iter().all(|v| v.abs() < 1e-6), "{last:?}");
    let b = bdf(rhs_stiff, 0.0, 1.0, &[1.0, 2.0], &opts, 5).expect("bdf");
    assert!(b.y.last().unwrap().iter().all(|v| v.abs() < 1e-6));
}

#[test]
fn finite_time_blow_up_is_reported_not_returned_as_a_value() {
    // y' = y², y(0) = 1 blows up at t = 1; integrating to t = 2 must fail
    // (or at least stop short), never return a "solution" at t = 2.
    let opts = OdeOptions {
        max_steps: 100_000,
        ..OdeOptions::default()
    };
    match dopri5(rhs_blowup, 0.0, 2.0, &[1.0], &opts, &[]) {
        | Err(_) => {},
        | Ok(sol) => {
            let t_end = *sol.t.last().unwrap();
            assert!(t_end < 1.0 + 1e-6, "claimed a solution up to t = {t_end}");
        },
    }
}

#[test]
fn root_finders_survive_degenerate_input() {
    let mut failures = Vec::new();
    let polys: Vec<Vec<f64>> = vec![
        vec![],
        vec![0.0],
        vec![5.0],
        vec![0.0, 0.0, 0.0],
        vec![1.0, 0.0],
        vec![0.0, 0.0, 1.0],
        vec![1.0, 2.0, 0.0, 0.0],
        vec![NAN, 1.0],
        vec![1.0, NAN, 1.0],
        vec![INF, 1.0, 1.0],
        vec![1e300, -1e300, 1e300],
        vec![1e-300, 1e300, 1e-300],
        vec![1.0, -3.0, 3.0, -1.0],
        vec![1.0, 0.0, 0.0, 0.0, 0.0, -1.0],
        // Wilkinson-20 style coefficients, repeated roots
        (0..21)
            .map(|k| if k % 2 == 0 { 1.0 } else { -2.0 })
            .collect(),
        vec![1.0; 61],
    ];
    for p in polys {
        let (a, b, c, d, e) = (p.clone(), p.clone(), p.clone(), p.clone(), p);
        failures.extend(guarded("aberth", move || {
            let _ = rootfind::polynomial_roots_aberth(&a);
        }));
        failures.extend(guarded("companion", move || {
            let _ = rootfind::polynomial_roots_companion(&b);
        }));
        failures.extend(guarded("jenkins-traub", move || {
            let _ = rootfind::polynomial_roots_jenkins_traub(&c);
        }));
        failures.extend(guarded("homotopy", move || {
            let _ = rootfind::homotopy_univariate(&d);
        }));
        failures.extend(guarded("polynomial::find_roots", move || {
            let _ = rssn::kernels::polynomial::Polynomial::new(e).find_roots();
        }));
    }
    for (lo, hi, tol, iter) in [
        (0.0, 1.0, 1e-12, 100),
        (1.0, 0.0, 1e-12, 100),
        (NAN, 1.0, 1e-12, 100),
        (0.0, INF, 1e-12, 100),
        (0.0, 1.0, 0.0, 100),
        (0.0, 1.0, NAN, 100),
        (0.0, 1.0, 1e-12, 0),
        (-1e308, 1e308, 1e-12, 100_000),
    ] {
        for f in [
            (|x: f64| x * x + 1.0) as fn(f64) -> f64,
            |x| x - 0.3,
            |_| NAN,
            |x| 1.0 / x,
            |x| x.tan(),
        ] {
            failures.extend(guarded("brent", move || {
                let _ = rootfind::brent(f, lo, hi, tol, iter);
            }));
        }
    }
    assert!(
        failures.is_empty(),
        "{} failures:\n{}",
        failures.len(),
        failures.join("\n")
    );
}

#[test]
fn brent_rejects_a_bracket_without_sign_change() {
    assert!(rootfind::brent(|x| x * x + 1.0, -1.0, 1.0, 1e-12, 100).is_err());
    let r = rootfind::brent(|x| x * x * x - 2.0, 0.0, 2.0, 1e-14, 200).unwrap();
    assert!((r - 2f64.cbrt()).abs() < 1e-12);
}

#[test]
fn optimisers_survive_degenerate_objectives() {
    let mut failures = Vec::new();
    let fns: Vec<(&'static str, fn(&[f64]) -> f64)> = vec![
        ("nan", |_| NAN),
        ("inf", |_| INF),
        ("neg-inf", |x| -x.iter().map(|v| v * v).sum::<f64>() * INF),
        ("flat", |_| 1.0),
        ("quadratic", |x| x.iter().map(|v| v * v).sum()),
        ("discontinuous", |x| {
            if x.first().copied().unwrap_or(0.0) > 0.0 {
                1.0
            } else {
                0.0
            }
        }),
        ("unbounded", |x| x.iter().sum::<f64>()),
    ];
    let starts: Vec<Vec<f64>> = vec![
        vec![],
        vec![0.0],
        vec![NAN],
        vec![1.0, -2.0],
        vec![INF, 1.0],
        vec![1e300; 3],
    ];
    for (name, f) in fns {
        for x0 in &starts {
            for (step, tol, iters) in [
                (0.1, 1e-8, 200_usize),
                (0.0, 1e-8, 50),
                (NAN, 0.0, 50),
                (1.0, -1.0, 0),
            ] {
                let x0 = x0.clone();
                failures.extend(guarded(name, move || {
                    let _ = optim::nelder_mead(f, &x0, step, tol, iters);
                }));
            }
            let x0c = x0.clone();
            failures.extend(guarded(name, move || {
                let _ = optim::bfgs::<_, fn(&[f64]) -> Vec<f64>>(f, None, &x0c, 1e-8, 100);
            }));
            let x0c = x0.clone();
            failures.extend(guarded(name, move || {
                let _ = optim::lbfgs::<_, fn(&[f64]) -> Vec<f64>>(f, None, &x0c, 5, 1e-8, 100);
            }));
            let x0c = x0.clone();
            failures.extend(guarded(name, move || {
                let _ = optim::numerical_gradient(&f, &x0c);
            }));
        }
    }
    assert!(
        failures.is_empty(),
        "{} failures:\n{}",
        failures.len(),
        failures.join("\n")
    );
}

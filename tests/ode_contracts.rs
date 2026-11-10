//! Adaptive solvers share the same host callback and tolerance contracts.
use rssn::kernels::ode_adaptive::Event;
use rssn::kernels::ode_adaptive::OdeError;
use rssn::kernels::ode_adaptive::OdeOptions;
use rssn::kernels::ode_adaptive::OdeSolution;
use rssn::kernels::ode_adaptive::{
    self,
};
use rssn::kernels::ode_stiff;
use std::cell::Cell;

type Rhs<'a> = &'a dyn Fn(f64, &[f64], &mut [f64]);
type Solver =
    fn(Rhs<'_>, f64, f64, &[f64], &OdeOptions, &[Event<'_>]) -> Result<OdeSolution, OdeError>;
fn solvers() -> [(&'static str, Solver); 5] {
    [
        ("dopri5", |f, a, b, y, o, e| {
            ode_adaptive::dopri5(f, a, b, y, o, e)
        }),
        ("rosenbrock23", |f, a, b, y, o, e| {
            ode_adaptive::rosenbrock23_events(f, a, b, y, o, e)
        }),
        ("radau5", |f, a, b, y, o, e| {
            ode_stiff::radau5_events(f, a, b, y, o, e)
        }),
        ("bdf", |f, a, b, y, o, e| {
            ode_stiff::bdf_events(f, a, b, y, o, 5, e)
        }),
        ("sdirk3", |f, a, b, y, o, e| {
            let options = rssn::kernels::ode_sdirk::SdirkOptions {
                rtol: o.rtol,
                atol: o.atol,
                h0: o.h0,
                hmax: o.hmax,
                max_steps: o.max_steps,
                ..Default::default()
            };
            rssn::kernels::ode_sdirk::sdirk3_events(
                f,
                |_, y| rssn::kernels::dense::Mat::zeros(y.len(), y.len()),
                a,
                b,
                y,
                &options,
                e,
            )
        }),
    ]
}

#[test]
fn malformed_problems_never_invoke_host_callbacks() {
    let cases = [
        OdeOptions {
            rtol: f64::NAN,
            ..OdeOptions::default()
        },
        OdeOptions {
            atol: -1.0,
            ..OdeOptions::default()
        },
        OdeOptions {
            rtol: 0.0,
            atol: 0.0,
            ..OdeOptions::default()
        },
        OdeOptions {
            h0: Some(0.0),
            ..OdeOptions::default()
        },
        OdeOptions {
            h0: Some(f64::INFINITY),
            ..OdeOptions::default()
        },
        OdeOptions {
            hmax: f64::NAN,
            ..OdeOptions::default()
        },
        OdeOptions {
            hmax: -1.0,
            ..OdeOptions::default()
        },
        OdeOptions {
            max_steps: 0,
            ..OdeOptions::default()
        },
    ];
    for (name, solve) in solvers() {
        let calls = Cell::new(0);
        let f = |_: f64, _: &[f64], dy: &mut [f64]| {
            calls.set(calls.get() + 1);
            dy.fill(0.0)
        };
        for options in &cases {
            assert_eq!(
                solve(&f, 0.0, 1.0, &[1.0], options, &[]).unwrap_err(),
                OdeError::Invalid,
                "{name}: {options:?}"
            );
        }
        for (a, b, y) in [
            (f64::NAN, 1.0, vec![1.0]),
            (0.0, f64::INFINITY, vec![1.0]),
            (0.0, 1.0, vec![f64::NAN]),
        ] {
            assert_eq!(
                solve(&f, a, b, &y, &OdeOptions::default(), &[]).unwrap_err(),
                OdeError::Invalid,
                "{name}"
            );
        }
        let event = Event {
            g: &|_, _| {
                calls.set(calls.get() + 1);
                0.0
            },
            terminal: false,
            direction: 2,
        };
        assert_eq!(
            solve(&f, 0.0, 1.0, &[1.0], &OdeOptions::default(), &[event]).unwrap_err(),
            OdeError::Invalid,
            "{name}"
        );
        assert_eq!(calls.get(), 0, "{name} evaluated malformed input");
    }
}

#[test]
fn an_incomplete_derivative_is_not_a_zero_or_stale_component() {
    for (name, solve) in solvers() {
        let error = solve(
            &|_, _, dy| dy[0] = 0.0,
            0.0,
            1.0,
            &[0.0, 0.0],
            &OdeOptions::default(),
            &[],
        )
        .unwrap_err();
        assert_eq!(error, OdeError::NonFinite, "{name}");
    }
}

#[test]
fn relative_only_tolerances_accept_the_exact_zero_solution() {
    let options = OdeOptions {
        atol: 0.0,
        ..OdeOptions::default()
    };
    for (name, solve) in solvers() {
        let solution = solve(&|_, _, dy| dy.fill(0.0), 0.0, 0.01, &[0.0], &options, &[])
            .unwrap_or_else(|error| panic!("{name}: {error:?}"));
        assert_eq!(solution.t.last(), Some(&0.01), "{name}");
        assert!(solution.y.iter().all(|state| state == &[0.0]), "{name}");
    }
}

#[test]
fn startup_and_time_derivative_probes_stay_inside_short_intervals() {
    let options = OdeOptions {
        hmax: 1e-12,
        ..OdeOptions::default()
    };
    for (name, solve) in solvers() {
        for (a, b) in [(0.0, 1e-12), (1e-12, 0.0)] {
            let f = |time: f64, _: &[f64], dy: &mut [f64]| {
                assert!(
                    (0.0..=1e-12).contains(&time),
                    "{name}: time={time} outside [{a},{b}]"
                );
                dy.fill(0.0);
            };
            let result = solve(&f, a, b, &[0.0], &options, &[])
                .unwrap_or_else(|error| panic!("{name}: {error:?}"));
            assert_eq!(result.t.last(), Some(&b), "{name}");
        }
    }
}

#[test]
fn public_rhs_counts_include_startup_and_terminal_events() {
    for (name, solve) in solvers() {
        for terminal in [false, true] {
            let calls = Cell::new(0usize);
            let f = |_: f64, _: &[f64], dy: &mut [f64]| {
                calls.set(calls.get() + 1);
                dy.fill(1.0)
            };
            let event = Event {
                g: &|t, _| t - 0.004,
                terminal,
                direction: 1,
            };
            let solution = solve(&f, 0.0, 0.01, &[0.0], &OdeOptions::default(), &[event])
                .unwrap_or_else(|error| panic!("{name}: {error:?}"));
            assert_eq!(solution.nfev, calls.get(), "{name}: terminal={terminal}");
            if terminal {
                assert!((solution.t.last().unwrap() - 0.004).abs() < 1e-12, "{name}");
            }
        }
    }
}

#[test]
fn nonfinite_events_are_errors_at_start_endpoint_and_during_location() {
    let options = OdeOptions {
        h0: Some(0.01),
        ..OdeOptions::default()
    };
    for (name, solve) in solvers() {
        for kind in 0..3 {
            let event_fn = |t: f64, _: &[f64]| {
                match kind {
                    | 0 => f64::NAN,
                    | 1 => {
                        if t > 0.0 {
                            f64::NAN
                        } else {
                            -1.0
                        }
                    },
                    | _ => {
                        if t > 0.0 && t < 0.01 {
                            f64::NAN
                        } else {
                            t - 0.005
                        }
                    },
                }
            };
            let event = Event {
                g: &event_fn,
                terminal: true,
                direction: 0,
            };
            let error = solve(
                &|_, _, dy| dy.fill(0.0),
                0.0,
                0.01,
                &[0.0],
                &options,
                &[event],
            )
            .unwrap_err();
            assert_eq!(error, OdeError::NonFinite, "{name}: event case {kind}");
        }
    }
}

#[test]
fn bdf_validates_before_time_dependent_mass_or_jacobian_callbacks() {
    let calls = Cell::new(0);
    let mass = |_: f64| {
        calls.set(calls.get() + 1);
        rssn::kernels::dense::Mat::identity(1)
    };
    let bo = ode_stiff::BdfOptions {
        mass_fn: Some(&mass),
        ..Default::default()
    };
    let options = OdeOptions {
        atol: f64::NAN,
        ..Default::default()
    };
    assert_eq!(
        ode_stiff::bdf_with(|_, _, dy| dy.fill(0.0), 0.0, 1.0, &[0.0], &options, &bo).unwrap_err(),
        OdeError::Invalid
    );
    assert_eq!(calls.get(), 0);
}

#[test]
fn exponential_rosenbrock_obeys_the_same_tolerance_probe_and_counter_contracts() {
    use rssn::kernels::ode_ext::exprb32;
    let calls = Cell::new(0usize);
    let f = |time: f64, _: &[f64], dy: &mut [f64]| {
        assert!((0.0..=1e-12).contains(&time), "exp-Rosenbrock probe {time}");
        calls.set(calls.get() + 1);
        dy.fill(0.0);
    };
    let bad = OdeOptions {
        atol: f64::NAN,
        ..Default::default()
    };
    assert_eq!(
        exprb32(f, None, 0.0, 1e-12, &[0.0], &bad).unwrap_err(),
        OdeError::Invalid
    );
    assert_eq!(calls.get(), 0);
    let options = OdeOptions {
        atol: 0.0,
        hmax: 1e-12,
        ..Default::default()
    };
    for (a, b) in [(0.0, 1e-12), (1e-12, 0.0)] {
        calls.set(0);
        let solution = exprb32(f, None, a, b, &[0.0], &options).unwrap();
        assert_eq!(solution.t.last(), Some(&b));
        assert_eq!(solution.nfev, calls.get());
    }
}

#[test]
fn rejected_newton_attempts_are_included_in_all_sdirk_work_counters() {
    use rssn::kernels::ode_sdirk::SdirkOptions;
    use rssn::kernels::ode_sdirk::sdirk3;
    let rhs_calls = Cell::new(0usize);
    let jac_calls = Cell::new(0usize);
    let f = |_: f64, _: &[f64], dy: &mut [f64]| {
        rhs_calls.set(rhs_calls.get() + 1);
        dy.fill(1.0)
    };
    let jac = |_: f64, _: &[f64]| {
        let count = jac_calls.get();
        jac_calls.set(count + 1);
        rssn::kernels::dense::Mat::from_rows(&[vec![if count < 2 { -1e6 } else { 0.0 }]]).unwrap()
    };
    let options = SdirkOptions {
        h0: Some(0.1),
        newton_max_iter: 2,
        ..Default::default()
    };
    let solution = sdirk3(f, jac, 0.0, 0.1, &[0.0], &options).unwrap();
    assert!(solution.rejected > 0);
    assert_eq!(solution.nfev, rhs_calls.get());
    assert_eq!(solution.njev, jac_calls.get());
    assert_eq!(solution.newton_iterations, jac_calls.get());
}

#[test]
fn finite_difference_sdirk_counts_the_rhs_inside_its_jacobians() {
    use rssn::kernels::ode_sdirk::SdirkOptions;
    use rssn::kernels::ode_sdirk::sdirk3_fd;
    let calls = Cell::new(0usize);
    let f = |_: f64, y: &[f64], dy: &mut [f64]| {
        calls.set(calls.get() + 1);
        for (out, input) in dy.iter_mut().zip(y) {
            *out = -*input;
        }
    };
    let result = sdirk3_fd(f, 0.0, 0.01, &[1.0, 2.0], &SdirkOptions::default()).unwrap();
    assert!(result.njev > 0);
    assert_eq!(result.nfev, calls.get());
}

#[test]
fn a_fixed_grid_cannot_claim_an_unrepresentable_substep() {
    use rssn::kernels::ode_sdirk::sdirk3_fixed;
    use rssn::kernels::ode_sdirk::sdirk3_fixed_events;
    let upper = f64::from_bits(1.0_f64.to_bits() + 1);
    for (a, b) in [(1.0, upper), (upper, 1.0)] {
        let f = |_: f64, _: &[f64], dy: &mut [f64]| dy.fill(0.0);
        let jac = |_: f64, _: &[f64]| rssn::kernels::dense::Mat::zeros(1, 1);
        assert_eq!(
            sdirk3_fixed(f, jac, a, b, &[0.0], 2).unwrap_err(),
            OdeError::StepLimit
        );
        assert_eq!(
            sdirk3_fixed_events(f, jac, a, b, &[0.0], 2, &[]).unwrap_err(),
            OdeError::StepLimit
        );
        assert_eq!(
            ode_stiff::radau5_fixed(f, a, b, &[0.0], 2).unwrap_err(),
            OdeError::StepLimit
        );
    }
}

#[test]
fn adjacent_float_intervals_are_valid_when_one_step_is_representable() {
    let upper = f64::from_bits(1.0_f64.to_bits() + 1);
    for (name, solve) in solvers() {
        for (a, b) in [(1.0, upper), (upper, 1.0)] {
            let options = OdeOptions {
                hmax: (b - a).abs(),
                ..Default::default()
            };
            let solution = solve(&|_, _, dy| dy.fill(0.0), a, b, &[0.0], &options, &[])
                .unwrap_or_else(|error| panic!("{name}: {error:?}"));
            assert_eq!(solution.t.last(), Some(&b), "{name}");
        }
    }
}

#[test]
fn stiff_solver_budgets_include_rejected_attempts() {
    let options = OdeOptions {
        h0: Some(0.1),
        max_steps: 1,
        ..Default::default()
    };
    for (name, solve) in solvers()
        .into_iter()
        .filter(|(name, _)| matches!(*name, "radau5" | "bdf"))
    {
        let minimum = Cell::new(f64::INFINITY);
        let f = |time: f64, _: &[f64], dy: &mut [f64]| {
            if time > 0.0 {
                minimum.set(minimum.get().min(time));
                dy.fill(f64::NAN)
            } else {
                dy.fill(0.0)
            }
        };
        let error = solve(&f, 0.0, 0.1, &[0.0], &options, &[]).unwrap_err();
        assert_eq!(error, OdeError::StepLimit, "{name}");
        let first_attempt_minimum = if name == "radau5" {
            0.01
        } else {
            0.099
        };
        assert!(
            minimum.get() >= first_attempt_minimum,
            "{name}: a second rejected attempt consumed an exhausted budget, t={}",
            minimum.get()
        );
    }
}

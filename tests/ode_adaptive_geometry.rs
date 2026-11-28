//! Focused adaptive ODE/DAE/DDE accuracy, event and hard-budget contracts.
use rssn::kernels::ode_geom::*;
use std::cell::Cell;

#[test]
fn every_requested_geometric_method_refines_for_accuracy_and_keeps_dense_output() {
    for method in [
        GeometricMethod::Rk4,
        GeometricMethod::StormerVerlet,
        GeometricMethod::Ruth4,
        GeometricMethod::ImplicitMidpoint,
        GeometricMethod::GaussLegendre4,
        GeometricMethod::GaussLegendre6,
    ] {
        let mut errors = Vec::new();
        let mut nodes = Vec::new();
        for tol in [1e-3, 1e-8] {
            let calls = Cell::new(0);
            let mut f = |_: f64, y: &[f64]| {
                calls.set(calls.get() + 1);
                vec![y[1], -y[0]]
            };
            let sol = solve_adaptive_geometric(
                &mut f,
                &[1.0, 0.0],
                (0.0, 2.0),
                method,
                AdaptiveGeomOptions {
                    rtol: tol,
                    atol: tol * 0.01,
                    h0: 1.5,
                    ..Default::default()
                },
                &[],
                &mut |_, _, _| 1.0,
                &Cell::new(false),
            )
            .unwrap();
            assert_eq!(sol.nfev, calls.get());
            assert_eq!(*sol.t.last().unwrap(), 2.0);
            assert!(sol.rejected > 0);
            let dense = evaluate_history(&sol.t, &sol.y, &sol.rates, 0.731).unwrap();
            errors.push(
                (sol.y.last().unwrap()[0] - 2f64.cos())
                    .abs()
                    .max((dense[0] - 0.731f64.cos()).abs()),
            );
            nodes.push(sol.t.len());
        }
        assert!(nodes[1] > nodes[0], "{method:?}: {nodes:?}");
        assert!(errors[1] < errors[0] * 0.03, "{method:?}: {errors:?}");
    }
}
#[test]
fn geometric_adaptation_preserves_terminal_events_and_never_force_accepts() {
    let spec = [EventSpec {
        terminal: true,
        direction: -1,
    }];
    let mut f = |_: f64, y: &[f64]| vec![y[1], -y[0]];
    let sol = solve_adaptive_geometric(
        &mut f,
        &[1.0, 0.0],
        (0.0, 4.0),
        GeometricMethod::GaussLegendre6,
        AdaptiveGeomOptions::default(),
        &spec,
        &mut |_, _, y| y[0],
        &Cell::new(false),
    )
    .unwrap();
    assert!((sol.hits[0].t - std::f64::consts::FRAC_PI_2).abs() < 2e-7);
    assert!(sol.y.last().unwrap()[0].abs() < 1e-12);
    assert!(evaluate_history(&sol.t, &sol.y, &sol.rates, 2.0).is_err());
    let calls = Cell::new(0);
    let mut f = |_: f64, y: &[f64]| {
        calls.set(calls.get() + 1);
        y.to_vec()
    };
    assert!(
        solve_adaptive_geometric(
            &mut f,
            &[1.0],
            (0.0, 1.0),
            GeometricMethod::Rk4,
            AdaptiveGeomOptions {
                max_evaluations: 3,
                ..Default::default()
            },
            &[],
            &mut |_, _, _| 1.0,
            &Cell::new(false)
        )
        .is_err()
    );
    assert_eq!(calls.get(), 3);
    assert!(
        solve_adaptive_geometric(
            &mut f,
            &[1.0],
            (0.0, 1.0),
            GeometricMethod::Rk4,
            AdaptiveGeomOptions {
                h0: 0.5,
                h_min: 0.4,
                rtol: 1e-13,
                atol: 1e-15,
                ..Default::default()
            },
            &[],
            &mut |_, _, _| 1.0,
            &Cell::new(false)
        )
        .is_err()
    );
}
#[test]
fn residual_dae_bdf_uses_accuracy_control_and_validates_every_residual() {
    let mut errors = Vec::new();
    let mut nodes = Vec::new();
    for tol in [1e-4, 1e-7] {
        let calls = Cell::new(0);
        let mut f = |t: f64, y: &[f64], dy: &[f64]| {
            calls.set(calls.get() + 1);
            vec![dy[0] + y[0], y[1] - t]
        };
        let sol = solve_dae_implicit(
            &mut f,
            &[1.0, 0.0],
            &[-1.0, 1.0],
            (0.0, 1.0),
            DaeOptions {
                rtol: tol,
                atol: tol * 0.01,
                h0: 0.5,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(sol.nfev, calls.get());
        assert_eq!(*sol.t.last().unwrap(), 1.0);
        assert!((sol.y.last().unwrap()[1] - 1.0).abs() < 1e-10);
        errors.push((sol.y.last().unwrap()[0] - (-1f64).exp()).abs());
        nodes.push(sol.t.len());
    }
    assert!(errors[1] < errors[0] * 0.03);
    assert!(nodes[1] > nodes[0]);
    assert!(matches!(
        solve_dae_implicit(
            &mut |_, _, _| vec![],
            &[1.0],
            &[-1.0],
            (0.0, 1.0),
            DaeOptions::default()
        ),
        Err(OdeGeomError::ShapeMismatch(_))
    ));
}
#[test]
fn semi_explicit_dae_counts_actual_calls_and_enforces_budgets() {
    let calls = Cell::new(0);
    let mut f = |_: f64, _: &[f64], z: &[f64]| {
        calls.set(calls.get() + 1);
        vec![-z[0]]
    };
    let mut g = |_: f64, y: &[f64], z: &[f64]| {
        calls.set(calls.get() + 1);
        vec![z[0] - y[0]]
    };
    let sol = solve_dae_index1(
        &mut f,
        &mut g,
        &[1.0],
        &[1.0],
        (0.0, 1.0),
        DaeOptions {
            rtol: 1e-4,
            atol: 1e-6,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(sol.nfev, calls.get());
    assert!((sol.y.last().unwrap()[0] - (-1f64).exp()).abs() < 1e-3);
    calls.set(0);
    assert!(
        solve_dae_index1(
            &mut f,
            &mut g,
            &[1.0],
            &[1.0],
            (0.0, 1.0),
            DaeOptions {
                max_evaluations: 3,
                ..Default::default()
            }
        )
        .is_err()
    );
    assert_eq!(calls.get(), 3);
}
#[test]
fn delay_steps_use_accepted_history_and_reject_overlapping_trials() {
    let mut history = |_: f64| vec![1.0];
    let mut f = |_: f64, _: &[f64], past: &[Vec<f64>]| vec![-past[0][0]];
    let sol = solve_dde(
        &mut f,
        &mut history,
        &[1.0],
        &[1.0],
        (0.0, 2.0),
        DdeOptions {
            h0: 1.5,
            ..Default::default()
        },
    )
    .unwrap();
    assert!((sol.y.last().unwrap()[0] + 0.5).abs() < 1e-6);
    let mut f = |t: f64, _: &[f64], lookup: &dyn Fn(f64) -> Vec<f64>| vec![-lookup(t - 0.05)[0]];
    let state = solve_dde_state_dependent(
        &mut f,
        &mut history,
        &[1.0],
        (0.0, 0.2),
        DdeOptions {
            h0: 0.2,
            ..Default::default()
        },
    )
    .unwrap();
    assert!(state.rejected > 0);
    assert_eq!(*state.t.last().unwrap(), 0.2);
    let mut f = |_: f64, y: &[f64], _: &[Vec<f64>]| y.to_vec();
    assert!(
        solve_dde(
            &mut f,
            &mut history,
            &[1.0],
            &[1.0],
            (0.0, 2.0),
            DdeOptions {
                h0: 0.5,
                h_min: 0.4,
                rtol: 1e-13,
                atol: 1e-15,
                ..Default::default()
            }
        )
        .is_err()
    );
}

mod finite_ranges {
    use rssn::kernels::ode_geom::DaeOptions;
    use rssn::kernels::ode_geom::OdeGeomError;
    use rssn::kernels::ode_geom::hermite_segment;
    use rssn::kernels::ode_geom::solve_dae_implicit;

    #[test]
    fn hermite_segments_handle_finite_endpoints_with_overwide_spans() {
        let large = f64::MAX;
        assert_eq!(
            hermite_segment(-large, large, &[1.0], &[3.0], &[0.0], &[0.0], 0.0),
            Ok(vec![2.0])
        );
        assert_eq!(
            hermite_segment(-large, large, &[1.0], &[3.0], &[0.0], &[0.0], -large),
            Ok(vec![1.0])
        );
        assert_eq!(
            hermite_segment(-large, large, &[1.0], &[3.0], &[0.0], &[0.0], large),
            Ok(vec![3.0])
        );
        assert!(matches!(
            hermite_segment(-large, large, &[0.0], &[0.0], &[f64::MAX], &[0.0], 0.0),
            Err(OdeGeomError::Domain(_))
        ));
        let small_slope = 1.0 / large;
        let slope =
            hermite_segment(-large, large, &[0.0], &[0.0], &[small_slope], &[0.0], 0.0).unwrap();
        assert!((slope[0] - 0.25).abs() < f64::EPSILON);
        assert!(matches!(
            hermite_segment(0.0, 1.0, &[1.0], &[2.0], &[f64::NAN], &[0.0], 0.0),
            Err(OdeGeomError::Domain(_))
        ));
    }

    #[test]
    fn implicit_dae_never_calls_the_residual_with_unrepresentable_bdf_states() {
        let mut saw_nonfinite = false;
        let large = f64::MAX;
        let result = solve_dae_implicit(
            &mut |t, y, yp| {
                saw_nonfinite |= !t.is_finite() || y.iter().chain(yp).any(|v| !v.is_finite());
                vec![0.0]
            },
            &[large],
            &[large],
            (0.0, 2.0),
            DaeOptions {
                h0: 2.0,
                max_steps: 1,
                ..Default::default()
            },
        );
        assert!(!saw_nonfinite, "residual evaluated with a nonfinite trial");
        assert!(matches!(result, Err(OdeGeomError::NotConverged)));
    }

    #[test]
    fn implicit_dae_rejects_a_step_without_a_distinct_midpoint() {
        let mut calls = 0;
        let result = solve_dae_implicit(
            &mut |_t, y, yp| {
                calls += 1;
                assert!(y.iter().chain(yp).all(|v| v.is_finite()));
                vec![0.0]
            },
            &[1.0],
            &[0.0],
            (1.0, 1.0_f64.next_up()),
            DaeOptions::default(),
        );
        assert!(matches!(result, Err(OdeGeomError::NotConverged)));
        assert_eq!(calls, 3);
    }
}

mod geometric_time_grid {
    use rssn::kernels::ode::OdeStepMode;
    use rssn::kernels::ode_geom::AdaptiveGeomOptions;
    use rssn::kernels::ode_geom::EventSpec;
    use rssn::kernels::ode_geom::Events;
    use rssn::kernels::ode_geom::FixedGeomSolution;
    use rssn::kernels::ode_geom::GeometricMethod;
    use rssn::kernels::ode_geom::OdeGeomError;
    use rssn::kernels::ode_geom::evaluate_history;
    use rssn::kernels::ode_geom::solve_adaptive_geometric;
    use rssn::kernels::ode_geom::solve_geometric;
    use std::cell::Cell;

    const ULP: f64 = f64::EPSILON;

    fn time(ulps: u32) -> f64 {
        1.0 + f64::from(ulps) * ULP
    }

    fn options(initial: u32) -> AdaptiveGeomOptions {
        AdaptiveGeomOptions {
            h0: f64::from(initial) * ULP,
            h_min: ULP,
            atol: 0.1,
            rtol: 0.0,
            ..Default::default()
        }
    }

    fn solve_linear(
        span: (f64, f64),
        options: AdaptiveGeomOptions,
        specs: &[EventSpec],
        guard: &mut Events<'_>,
        calls: &Cell<usize>,
    ) -> Result<FixedGeomSolution, OdeGeomError> {
        solve_adaptive_geometric(
            &mut |_, _| {
                calls.set(calls.get() + 1);
                vec![1.0 / ULP]
            },
            &[0.0],
            span,
            GeometricMethod::Rk4,
            options,
            specs,
            guard,
            &Cell::new(false),
        )
    }

    #[test]
    fn adaptive_geometric_rejects_an_unsplittable_trial_and_keeps_reverse_contract() {
        let calls = Cell::new(0);
        let result = solve_linear(
            (time(0), time(1)),
            options(1),
            &[],
            &mut |_, _, _| unreachable!(),
            &calls,
        );
        assert!(matches!(result, Err(OdeGeomError::NotConverged)));
        assert_eq!(calls.get(), 1);

        let calls = Cell::new(0);
        let result = solve_linear(
            (time(3), time(0)),
            options(3),
            &[EventSpec {
                direction: 0,
                terminal: true,
            }],
            &mut |_, _, _| unreachable!(),
            &calls,
        );
        assert!(matches!(result, Err(OdeGeomError::Domain(_))));
        assert_eq!(calls.get(), 0);
    }

    #[test]
    fn adaptive_geometric_checks_hermite_at_the_actual_odd_grid_midpoint() {
        for (start, end, middle, middle_state) in [(0, 3, 2, 2.0), (1, 4, 2, 1.0)] {
            let calls = Cell::new(0);
            let answer = solve_linear(
                (time(start), time(end)),
                AdaptiveGeomOptions {
                    max_steps: 1,
                    max_evaluations: 15,
                    ..options(3)
                },
                &[],
                &mut |_, _, _| unreachable!(),
                &calls,
            )
            .unwrap();
            assert_eq!(answer.t, [time(start), time(middle), time(end)]);
            assert_eq!(answer.y, [vec![0.0], vec![middle_state], vec![3.0]]);
            assert_eq!(answer.nfev, 15);
            assert_eq!(answer.nfev, calls.get());
            assert_eq!(answer.rejected, 0);
        }
    }

    #[test]
    fn adaptive_geometric_balances_the_final_clipped_grid_within_step_bounds() {
        let calls = Cell::new(0);
        let joined = solve_linear(
            (time(0), time(3)),
            AdaptiveGeomOptions {
                max_steps: 1,
                ..options(2)
            },
            &[],
            &mut |_, _, _| unreachable!(),
            &calls,
        )
        .unwrap();
        assert_eq!(joined.t, [time(0), time(2), time(3)]);
        assert_eq!(joined.nfev, 15);

        let calls = Cell::new(0);
        let split = solve_linear(
            (time(0), time(5)),
            AdaptiveGeomOptions {
                max_steps: 2,
                h_max: 4.0 * ULP,
                h_min: 3.0 * ULP,
                ..options(4)
            },
            &[],
            &mut |_, _, _| unreachable!(),
            &calls,
        )
        .unwrap();
        assert_eq!(split.t, [time(0), time(2), time(3), time(4), time(5)]);
        assert_eq!(split.nfev, 29);
        assert_eq!(split.nfev, calls.get());

        let calls = Cell::new(0);
        let blocked = solve_linear(
            (time(0), time(3)),
            AdaptiveGeomOptions {
                h_max: 2.0 * ULP,
                ..options(2)
            },
            &[],
            &mut |_, _, _| unreachable!(),
            &calls,
        );
        assert!(matches!(blocked, Err(OdeGeomError::NotConverged)));
        assert_eq!(calls.get(), 1);

        let calls = Cell::new(0);
        let out_of_budget = solve_linear(
            (time(0), time(3)),
            AdaptiveGeomOptions {
                max_evaluations: 14,
                ..options(2)
            },
            &[],
            &mut |_, _, _| unreachable!(),
            &calls,
        );
        assert!(matches!(out_of_budget, Err(OdeGeomError::NotConverged)));
        assert_eq!(calls.get(), 14);
    }

    #[test]
    fn adaptive_geometric_events_report_state_at_a_representable_time_in_order() {
        let calls = Cell::new(0);
        let specs = [
            EventSpec {
                direction: 1,
                terminal: true,
            },
            EventSpec {
                direction: 1,
                terminal: false,
            },
            EventSpec {
                direction: -1,
                terminal: false,
            },
            EventSpec {
                direction: 1,
                terminal: false,
            },
        ];
        let mut seen = Vec::new();
        let mut guard = |index: usize, t: f64, y: &[f64]| {
            assert!((y[0] - (t - time(0)) / ULP).abs() <= 1e-14);
            assert!(!seen.contains(&(index, t)));
            seen.push((index, t));
            match index {
                | 0 => y[0] - 2.25,
                | 1 => y[0] - 0.25,
                | 2 => 1.25 - y[0],
                | 3 => y[0] - 3.0,
                | _ => unreachable!(),
            }
        };
        let answer =
            solve_linear((time(0), time(6)), options(6), &specs, &mut guard, &calls).unwrap();
        let hits: Vec<_> = answer.hits.iter().map(|hit| (hit.index, hit.t)).collect();
        assert_eq!(hits, [(1, time(1)), (2, time(2)), (0, time(3))]);
        for (hit, expected) in answer.hits.iter().zip([1.0, 2.0, 3.0]) {
            assert!((hit.y[0] - expected).abs() <= 1e-14);
        }
        assert_eq!(answer.t, [time(0), time(3)]);
        assert_eq!(answer.rates[1], [1.0 / ULP]);
        assert_eq!(answer.nfev, 15);
        assert_eq!(answer.nfev, calls.get());
        assert_eq!(seen.len(), 13);
    }

    #[test]
    fn adaptive_geometric_terminal_root_in_subnormal_time_uses_the_later_endpoint() {
        let tiny = f64::from_bits(1);
        let specs = [
            EventSpec {
                direction: 0,
                terminal: false,
            },
            EventSpec {
                direction: 1,
                terminal: true,
            },
        ];
        let mut calls = 0;
        let mut guard_calls = 0;
        let answer = solve_adaptive_geometric(
            &mut |_, _| {
                calls += 1;
                vec![0.0]
            },
            &[1.0],
            (-tiny, 2.0 * tiny),
            GeometricMethod::Rk4,
            AdaptiveGeomOptions {
                h0: 2.0 * tiny,
                h_min: tiny,
                max_steps: 1,
                ..Default::default()
            },
            &specs,
            &mut |index, t, y| {
                guard_calls += 1;
                assert_eq!(y, [1.0]);
                t / tiny + if index == 0 { 0.75 } else { 0.25 }
            },
            &Cell::new(false),
        )
        .unwrap();
        assert_eq!(answer.t, [-tiny, 0.0]);
        assert_eq!(answer.y, [vec![1.0], vec![1.0]]);
        assert_eq!(
            answer
                .hits
                .iter()
                .map(|hit| (hit.index, hit.t))
                .collect::<Vec<_>>(),
            [(0, 0.0), (1, 0.0)]
        );
        assert_eq!(guard_calls, 4);
        assert_eq!(answer.nfev, 15);
        assert_eq!(answer.nfev, calls);
    }
    #[test]
    fn adaptive_geometric_reschedules_a_rejected_final_step_after_prior_progress() {
        let mut rhs_times = Vec::new();
        let mut accepted_guards = Vec::new();
        let answer = solve_adaptive_geometric(
            &mut |t, _| {
                rhs_times.push(t);
                let elapsed = ((t - time(2)) / ULP).max(0.0);
                vec![0.002 * elapsed.powi(3) / ULP]
            },
            &[0.0],
            (time(0), time(8)),
            GeometricMethod::Rk4,
            AdaptiveGeomOptions {
                max_steps: 4,
                max_evaluations: 57,
                ..options(2)
            },
            &[EventSpec {
                direction: 0,
                terminal: false,
            }],
            &mut |_, t, _| {
                accepted_guards.push(t);
                1.0
            },
            &Cell::new(false),
        )
        .unwrap();
        // Each RK4 trial makes three four-stage calls and two endpoint calls.
        let trial_bounds: Vec<_> = rhs_times[1..]
            .as_chunks::<14>()
            .0
            .iter()
            .map(|trial| (trial[0], trial[3]))
            .collect();
        assert_eq!(
            trial_bounds,
            [
                (time(0), time(2)),
                (time(2), time(8)),
                (time(2), time(6)),
                (time(6), time(8)),
            ]
        );
        assert_eq!(
            answer.t,
            [
                time(0),
                time(1),
                time(2),
                time(4),
                time(6),
                time(7),
                time(8)
            ]
        );
        assert_eq!(answer.nfev, 57);
        assert_eq!(answer.nfev, rhs_times.len());
        assert_eq!(answer.rejected, 1);
        assert_eq!(accepted_guards, answer.t);
        let exact_endpoint = 0.002 * 6.0_f64.powi(4) / 4.0;
        assert!((answer.y.last().unwrap()[0] - exact_endpoint).abs() < 0.1);
    }

    #[test]
    fn geometric_dispatch_keeps_dense_history_shape_and_exact_terminal_endpoint() {
        let mut calls = 0;
        let answer = solve_geometric(
            &mut |_, _| {
                calls += 1;
                vec![1.0 / ULP]
            },
            &[0.0],
            (time(0), time(3)),
            GeometricMethod::Rk4,
            OdeStepMode::Adaptive {
                rtol: 0.0,
                atol: 0.1,
                dt_init: Some(2.0 * ULP),
                dt_min: Some(ULP),
                dt_max: None,
                max_steps: 1,
            },
            &[EventSpec {
                direction: 1,
                terminal: true,
            }],
            &mut |_, _, y| y[0] - 3.0,
            &Cell::new(false),
        )
        .unwrap();
        assert_eq!(answer.t, [time(0), time(2), time(3)]);
        assert_eq!(answer.y, [vec![0.0], vec![2.0], vec![3.0]]);
        assert_eq!(answer.rates, vec![vec![1.0 / ULP]; 3]);
        assert_eq!(answer.nfev, 15);
        assert_eq!(answer.nfev, calls);
        assert_eq!(answer.hits[0].t, time(3));
        assert_eq!(answer.hits[0].y, [3.0]);
        assert_eq!(
            evaluate_history(&answer.t, &answer.y, &answer.rates, time(1)).unwrap(),
            [1.0],
        );
    }
}

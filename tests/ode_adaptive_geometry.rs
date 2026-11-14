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

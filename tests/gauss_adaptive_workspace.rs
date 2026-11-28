use std::cell::Cell;

use rssn::kernels::ode_geom::AdaptiveGeomOptions;
use rssn::kernels::ode_geom::EventSpec;
use rssn::kernels::ode_geom::FixedGeomSolution;
use rssn::kernels::ode_geom::GeometricMethod;
use rssn::kernels::ode_geom::OdeGeomError;
use rssn::kernels::ode_geom::solve_adaptive_geometric;

#[test]
fn constant_rate_trials_keep_callback_order_and_exact_budget() {
    let offset2 = 3.0_f64.sqrt() / 6.0;
    let offset3 = 15.0_f64.sqrt() / 10.0;
    for (method, nodes) in [
        (GeometricMethod::ImplicitMidpoint, vec![0.5]),
        (
            GeometricMethod::GaussLegendre4,
            vec![0.5 - offset2, 0.5 + offset2],
        ),
        (
            GeometricMethod::GaussLegendre6,
            vec![0.5 - offset3, 0.5, 0.5 + offset3],
        ),
    ] {
        let mut expected_times = vec![0.0];
        for (start, step) in [(0.0, 1.0), (0.0, 0.5), (0.5, 0.5)] {
            expected_times.push(start);
            expected_times.extend(nodes.iter().map(|node| start + step * node));
        }
        expected_times.extend([0.5, 1.0]);

        let mut actual_times = Vec::new();
        let solution = solve_adaptive_geometric(
            &mut |t, y| {
                actual_times.push(t);
                assert!((y[0] - (2.0 + t)).abs() < 4e-15);
                vec![1.0]
            },
            &[2.0],
            (0.0, 1.0),
            method,
            AdaptiveGeomOptions {
                h0: 1.0,
                max_evaluations: expected_times.len(),
                ..Default::default()
            },
            &[],
            &mut |_, _, _| 1.0,
            &Cell::new(false),
        )
        .unwrap();
        assert_eq!(solution.nfev, expected_times.len());
        assert_eq!(solution.rejected, 0);
        assert_eq!(actual_times, expected_times);
        assert_eq!(solution.t, [0.0, 0.5, 1.0]);
        assert!((solution.y[2][0] - 3.0).abs() < 4e-15);
    }
}

fn harmonic(max_evaluations: usize) -> (Result<FixedGeomSolution, OdeGeomError>, Vec<[u64; 3]>) {
    let mut trace = Vec::new();
    let result = solve_adaptive_geometric(
        &mut |t, y| {
            trace.push([t.to_bits(), y[0].to_bits(), y[1].to_bits()]);
            vec![-y[1], y[0]]
        },
        &[1.0, 0.0],
        (0.0, 2.0),
        GeometricMethod::GaussLegendre6,
        AdaptiveGeomOptions {
            rtol: 1e-9,
            atol: 1e-11,
            h0: 1.5,
            max_evaluations,
            ..Default::default()
        },
        &[],
        &mut |_, _, _| 1.0,
        &Cell::new(false),
    );
    (result, trace)
}

#[test]
fn rejected_retries_and_small_budgets_replay_the_same_callback_prefix() {
    let (result, accepted_trace) = harmonic(100_000);
    let result = result.unwrap();
    assert!(result.rejected > 0);
    assert_eq!(result.nfev, accepted_trace.len());
    let fourth_start = accepted_trace
        .iter()
        .enumerate()
        .filter(|(_, observed)| observed[0] == 0.0_f64.to_bits())
        .nth(3)
        .unwrap()
        .0;
    for budget in [1, fourth_start + 1, accepted_trace.len() - 1] {
        let (result, trace) = harmonic(budget);
        assert!(matches!(result, Err(OdeGeomError::NotConverged)));
        assert_eq!(trace, accepted_trace[..budget]);
    }
    let final_state = result.y.last().unwrap();
    assert!((final_state[0] - 2.0_f64.cos()).hypot(final_state[1] - 2.0_f64.sin()) < 1e-9);
}

#[test]
fn terminal_event_in_second_half_preserves_the_hermite_rate_and_budget() {
    let calls = Cell::new(0usize);
    let invalid = Cell::new(false);
    let solution = solve_adaptive_geometric(
        &mut |_, _| {
            calls.set(calls.get() + 1);
            vec![1.0]
        },
        &[0.0],
        (0.0, 1.0),
        GeometricMethod::GaussLegendre4,
        AdaptiveGeomOptions {
            h0: 1.0,
            max_evaluations: 12,
            ..Default::default()
        },
        &[EventSpec {
            terminal: true,
            direction: 1,
        }],
        &mut |_, _, y| y[0] - 0.625,
        &invalid,
    )
    .unwrap();

    assert!(!invalid.get());
    assert_eq!(calls.get(), 12);
    assert_eq!(solution.nfev, calls.get());
    assert_eq!(solution.rejected, 0);
    assert_eq!(solution.hits.len(), 1);
    let hit = &solution.hits[0];
    assert!((hit.t - 0.625).abs() < 2e-15);
    assert!((hit.y[0] - 0.625).abs() < 2e-15);
    assert_eq!(solution.t.last().unwrap().to_bits(), hit.t.to_bits());
    assert!((solution.rates.last().unwrap()[0] - 1.0).abs() < 2e-15);
}

mod fixed_workspace {
    use rssn::kernels::ode_geom::GeometricMethod;
    use rssn::kernels::ode_geom::sample_fixed_geometric;
    use rssn::kernels::ode_geom::solve_fixed_geometric;
    use std::cell::Cell;

    fn gauss_times(
        start: f64,
        h: f64,
    ) -> [f64; 4] {
        let offset = 15.0_f64.sqrt() / 10.0;
        [
            start,
            start + h * (0.5 - offset),
            start + h * 0.5,
            start + h * (0.5 + offset),
        ]
    }

    #[test]
    fn fixed_workspace_preserves_explicit_stage_override_and_evaluations() {
        let mut seen = Vec::new();
        let solution = solve_fixed_geometric(
            &mut |time, _| {
                seen.push(time);
                vec![1.0]
            },
            &[2.0],
            (0.0, 1.0),
            2,
            GeometricMethod::Rk4,
            Some(3),
            &[],
            &mut |_, _, _| 1.0,
            &Cell::new(false),
        )
        .unwrap();
        let mut expected = vec![0.0];
        for start in [0.0, 0.5] {
            expected.extend(gauss_times(start, 0.5));
            expected.push(start + 0.5);
        }
        assert_eq!(seen, expected);
        assert_eq!(solution.nfev, expected.len());
        assert_eq!(solution.t, [0.0, 0.5, 1.0]);
        assert!((solution.y[2][0] - 3.0).abs() < 4e-15);
    }

    #[test]
    fn sampled_workspace_preserves_callback_order_across_output_times() {
        let mut seen = Vec::new();
        let result = sample_fixed_geometric(
            &mut |time, _| {
                seen.push(time);
                vec![1.0]
            },
            &[2.0],
            &[0.0, 0.0, 0.5, 1.0],
            GeometricMethod::GaussLegendre6,
            2,
        )
        .unwrap();
        let expected = [
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.25),
            (0.25, 0.25),
            (0.5, 0.25),
            (0.75, 0.25),
        ]
        .into_iter()
        .flat_map(|(start, h)| gauss_times(start, h))
        .collect::<Vec<_>>();
        assert_eq!(seen, expected);
        for (state, expected) in result.iter().zip([2.0, 2.0, 2.5, 3.0]) {
            assert!((state[0] - expected).abs() < 4e-15);
        }
    }
}

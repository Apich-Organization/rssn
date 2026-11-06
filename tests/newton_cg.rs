use rssn::kernels::newton_cg::{
    NewtonCgError, NewtonCgOptions, NewtonCgStatus, minimize, minimize_tape,
};

#[test]
fn matrix_free_diagonal_problem_uses_vector_actions_in_high_dimension() {
    let n = 512;
    let options = NewtonCgOptions {
        gradient_tolerance: 1e-7,
        max_iterations: 100,
        ..Default::default()
    };
    let result = minimize(
        |x| {
            x.iter()
                .enumerate()
                .map(|(i, x)| 0.5 * (1.0 + (i % 13) as f64) * (x - (i as f64).sin()).powi(2))
                .sum()
        },
        |x, g| {
            for (i, (x, g)) in x.iter().zip(g).enumerate() {
                *g = (1.0 + (i % 13) as f64) * (x - (i as f64).sin());
            }
        },
        |_, v, h| {
            for (i, (v, h)) in v.iter().zip(h).enumerate() {
                *h = (1.0 + (i % 13) as f64) * v;
            }
        },
        &vec![0.0; n],
        &options,
    )
    .unwrap();
    assert_eq!(result.status, NewtonCgStatus::FirstOrderConverged);
    assert!(result.value < 1e-12);
    assert!(result.hessian_vector_evaluations < n);
    assert!(result.gradient_norm <= options.gradient_tolerance);
    for (i, x) in result.x.iter().enumerate() {
        assert!((x - (i as f64).sin()).abs() < 1e-7);
    }
}

#[test]
fn negative_curvature_takes_the_boundary_and_escapes_a_nonstationary_region() {
    let result = minimize(
        |x| x[0].powi(4) - x[0] * x[0],
        |x, g| g[0] = 4.0 * x[0].powi(3) - 2.0 * x[0],
        |x, v, h| h[0] = (12.0 * x[0] * x[0] - 2.0) * v[0],
        &[0.1],
        &NewtonCgOptions::default(),
    )
    .unwrap();
    assert_eq!(result.status, NewtonCgStatus::FirstOrderConverged);
    assert!(result.negative_curvature_steps > 0 && result.boundary_steps > 0);
    assert!((result.x[0] - std::f64::consts::FRAC_1_SQRT_2).abs() < 1e-7);
    assert!((result.value + 0.25).abs() < 1e-12);
}

#[test]
fn domain_failures_are_rejected_and_do_not_replace_the_last_good_iterate() {
    let mut invalid = 0;
    let result = minimize(
        |x| {
            let value = x[0].ln().powi(2);
            if !value.is_finite() {
                invalid += 1;
            }
            value
        },
        |x, g| g[0] = 2.0 * x[0].ln() / x[0],
        |x, v, h| h[0] = 2.0 * (1.0 - x[0].ln()) / x[0].powi(2) * v[0],
        &[10.0],
        &NewtonCgOptions {
            initial_radius: 20.0,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(result.status, NewtonCgStatus::FirstOrderConverged);
    assert!(invalid > 0);
    assert!((result.x[0] - 1.0).abs() < 1e-7);
}

#[test]
fn no_progress_and_budget_exhaustion_are_not_convergence() {
    let options = NewtonCgOptions {
        max_iterations: 0,
        ..Default::default()
    };
    let result = minimize(
        |x| x[0] * x[0],
        |x, g| g[0] = 2.0 * x[0],
        |_, v, h| h[0] = 2.0 * v[0],
        &[2.0],
        &options,
    )
    .unwrap();
    assert_eq!(result.status, NewtonCgStatus::IterationLimit);
    assert_eq!(result.iterations, 0);
    assert_eq!(result.x, vec![2.0]);
    let result = minimize(
        |_| 1.0,
        |_, g| g[0] = 1.0,
        |_, v, h| h[0] = v[0],
        &[2.0],
        &NewtonCgOptions::default(),
    )
    .unwrap();
    assert_eq!(result.status, NewtonCgStatus::RadiusTooSmall);
    assert_eq!(result.x, vec![2.0]);
    assert_eq!(result.accepted_steps, 0);
}

#[test]
fn unwritten_derivative_buffers_are_detected_and_bad_trial_gradients_rejected() {
    let error = minimize(
        |x| x[0] * x[0],
        |_, _| {},
        |_, v, h| h[0] = 2.0 * v[0],
        &[1.0],
        &NewtonCgOptions::default(),
    )
    .unwrap_err();
    assert_eq!(error, NewtonCgError::NonFiniteInitial);
    let error = minimize(
        |x| x[0] * x[0],
        |x, g| g[0] = 2.0 * x[0],
        |_, _, _| {},
        &[1.0],
        &NewtonCgOptions::default(),
    )
    .unwrap_err();
    assert_eq!(error, NewtonCgError::NonFiniteModel);
    let result = minimize(
        |x| x[0] * x[0],
        |x, g| {
            if x[0] == 1.0 {
                g[0] = 2.0;
            }
        },
        |_, v, h| h[0] = 2.0 * v[0],
        &[1.0],
        &NewtonCgOptions {
            min_radius: 1e-6,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(result.status, NewtonCgStatus::RadiusTooSmall);
    assert_eq!(result.x, vec![1.0]);
    assert_eq!(result.accepted_steps, 0);
    assert!(result.gradient_evaluations > 1);
}

#[test]
fn stationary_saddles_are_only_first_order_certificates_and_invalid_options_fail() {
    let result = minimize(
        |x| -x[0] * x[0],
        |x, g| g[0] = -2.0 * x[0],
        |_, v, h| h[0] = -2.0 * v[0],
        &[0.0],
        &NewtonCgOptions::default(),
    )
    .unwrap();
    assert_eq!(result.status, NewtonCgStatus::FirstOrderConverged);
    assert_eq!(result.hessian_vector_evaluations, 0);
    for options in [
        NewtonCgOptions {
            initial_radius: 0.0,
            ..Default::default()
        },
        NewtonCgOptions {
            acceptance_threshold: 0.3,
            ..Default::default()
        },
        NewtonCgOptions {
            max_cg_iterations: 0,
            ..Default::default()
        },
        NewtonCgOptions {
            gradient_tolerance: f64::NAN,
            ..Default::default()
        },
    ] {
        assert_eq!(
            minimize(
                |_| 0.0,
                |_, g| g.fill(0.0),
                |_, _, h| h.fill(0.0),
                &[1.0],
                &options
            )
            .unwrap_err(),
            NewtonCgError::InvalidInput
        );
    }
}

#[test]
fn analytic_tape_reverse_and_hvp_drive_rosenbrock_without_a_dense_hessian() {
    use rssn::{
        backend::tape::lower,
        graph::{Engine, Graph},
    };
    let mut graph = Graph::new();
    Engine::install(&mut graph, &rssn::rules::standard()).unwrap();
    let root = graph.parse("(1-x)^2+100*(y-x^2)^2").unwrap();
    let inputs = [
        graph.interner_mut().symbol("x"),
        graph.interner_mut().symbol("y"),
    ];
    let tape = lower(&graph, &[root], &inputs, &[]).unwrap().optimise();
    let result = minimize_tape(&tape, 0, &[-1.2, 1.0], &[], &NewtonCgOptions::default()).unwrap();
    assert_eq!(
        result.status,
        NewtonCgStatus::FirstOrderConverged,
        "{result:?}"
    );
    assert!(result.value < 1e-14);
    assert!(result.x.iter().all(|x| (x - 1.0).abs() < 1e-6));
    assert_eq!(
        minimize_tape(&tape, 1, &[-1.2, 1.0], &[], &NewtonCgOptions::default()).unwrap_err(),
        NewtonCgError::InvalidInput
    );
}

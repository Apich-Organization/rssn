#[test]
fn overdetermined_fit_and_analytic_jacobian() {
    let f = |x: &[f64]| vec![x[0] + x[1] - 3., 2. * x[0] - x[1], x[0] - 1.];
    let jac = |_: &[f64]| {
        rssn::kernels::dense::Mat::from_rows(&[vec![1., 1.], vec![2., -1.], vec![1., 0.]]).unwrap()
    };
    let r = rssn::kernels::rootfind::least_squares(
        f,
        Some(&jac),
        &[0., 0.],
        rssn::kernels::rootfind::LeastSquaresOptions::default(),
    )
    .unwrap();
    assert!(r.converged(), "{r:?}");
    assert!(
        (r.x[0] - 1.).abs() < 1e-8 && (r.x[1] - 2.).abs() < 1e-8,
        "{r:?}"
    );
}
#[test]
fn rank_deficiency_and_nonzero_optimal_residual() {
    let f = |x: &[f64]| vec![x[0] + x[1] - 1., x[0] + x[1] - 2.];
    let r = rssn::kernels::rootfind::least_squares(
        f,
        None,
        &[0., 0.],
        rssn::kernels::rootfind::LeastSquaresOptions::default(),
    )
    .unwrap();
    assert_eq!(
        r.status,
        rssn::kernels::rootfind::LeastSquaresStatus::GradientConverged,
        "{r:?}"
    );
    assert!((r.x[0] + r.x[1] - 1.5).abs() < 1e-6, "{r:?}");
    assert!((r.residual_norm - 0.5f64.sqrt()).abs() < 1e-8);
}
#[test]
fn tiny_and_huge_scales_and_budget() {
    for scale in [1e-180, 1e180] {
        let f = |x: &[f64]| vec![scale * (x[0] - 2.)];
        let options = rssn::kernels::rootfind::LeastSquaresOptions {
            residual_tolerance: 0.,
            gradient_tolerance: 1e-10,
            ..Default::default()
        };
        let r = rssn::kernels::rootfind::least_squares(f, None, &[0.], options).unwrap();
        assert!(r.converged(), "{scale} {r:?}");
        assert!((r.x[0] - 2.).abs() < 1e-7, "{scale} {r:?}");
    }
    let o = rssn::kernels::rootfind::LeastSquaresOptions {
        max_evaluations: 1,
        ..Default::default()
    };
    let r = rssn::kernels::rootfind::least_squares(|x| vec![x[0] - 1.], None, &[0.], o).unwrap();
    assert_eq!(
        r.status,
        rssn::kernels::rootfind::LeastSquaresStatus::EvaluationLimit
    );
    assert_eq!(r.evaluations, 1);
}
#[test]
fn malformed_shapes_nonfinite_and_stagnation_fail_closed() {
    let o = rssn::kernels::rootfind::LeastSquaresOptions::default();
    assert_eq!(
        rssn::kernels::rootfind::least_squares(|_| vec![f64::NAN], None, &[0.], o).unwrap_err(),
        rssn::kernels::rootfind::LeastSquaresError::InvalidResidual
    );
    let j = |_: &[f64]| {
        rssn::kernels::dense::Mat {
            rows: 1,
            cols: 1,
            data: vec![],
        }
    };
    assert_eq!(
        rssn::kernels::rootfind::least_squares(|x| vec![x[0] + 1.], Some(&j), &[0.], o)
            .unwrap_err(),
        rssn::kernels::rootfind::LeastSquaresError::InvalidJacobian
    );
    let f = |x: &[f64]| {
        if x[0] == 0. {
            vec![1.]
        } else {
            vec![2.]
        }
    };
    let j = |_: &[f64]| rssn::kernels::dense::Mat::from_rows(&[vec![1.]]).unwrap();
    let r = rssn::kernels::rootfind::least_squares(f, Some(&j), &[0.], o).unwrap();
    assert_eq!(
        r.status,
        rssn::kernels::rootfind::LeastSquaresStatus::Stagnated
    );
    assert!(!r.converged());
}

#[test]
fn scaled_gradient_does_not_vanish_from_cancellation_or_product_underflow() {
    use rssn::kernels::dense::Mat;
    use rssn::kernels::rootfind::LeastSquaresOptions;
    use rssn::kernels::rootfind::LeastSquaresStatus;
    use rssn::kernels::rootfind::least_squares;
    let options = LeastSquaresOptions {
        gradient_tolerance: 0.0,
        residual_tolerance: 0.0,
        ..Default::default()
    };
    let residual = |_: &[f64]| vec![1.0, 1e-16, -1.0];
    let jac = |_: &[f64]| Mat::from_rows(&[vec![1.0], vec![1.0], vec![1.0]]).unwrap();
    let result = least_squares(residual, Some(&jac), &[0.0], options).unwrap();
    assert_ne!(
        result.status,
        LeastSquaresStatus::GradientConverged,
        "{result:?}"
    );
    let tiny = f64::from_bits(1);
    let residual = |_: &[f64]| vec![1.0, -1.0, tiny];
    let jac = |_: &[f64]| Mat::from_rows(&[vec![1.0], vec![1.0], vec![tiny]]).unwrap();
    let result = least_squares(residual, Some(&jac), &[0.0], options).unwrap();
    assert_ne!(
        result.status,
        LeastSquaresStatus::GradientConverged,
        "{result:?}"
    );
}

#[test]
fn resolves_ill_conditioned_columns_without_forming_normal_equations() {
    use rssn::kernels::dense::Mat;
    use rssn::kernels::rootfind::LeastSquaresOptions;
    use rssn::kernels::rootfind::least_squares;
    let epsilon = 1e-8;
    let residual = |x: &[f64]| {
        vec![
            x[0] + x[1] - 1.0,
            x[0] + (1.0 + epsilon) * x[1] - (1.0 - epsilon),
            x[0] + (1.0 - epsilon) * x[1] - (1.0 + epsilon),
        ]
    };
    let jacobian = |_: &[f64]| {
        Mat::from_rows(&[
            vec![1.0, 1.0],
            vec![1.0, 1.0 + epsilon],
            vec![1.0, 1.0 - epsilon],
        ])
        .unwrap()
    };
    let options = LeastSquaresOptions {
        max_iterations: 300,
        residual_tolerance: 1e-13,
        gradient_tolerance: 0.0,
        step_tolerance: 0.0,
        ..Default::default()
    };
    let result = least_squares(residual, Some(&jacobian), &[0.0, 0.0], options).unwrap();
    assert!(result.converged(), "{result:?}");
    assert!((result.x[0] - 2.0).abs() < 1e-4, "{result:?}");
}

#[test]
fn finite_difference_retries_the_admissible_side_at_a_domain_boundary() {
    let residual = |x: &[f64]| vec![(-x[0]).sqrt() - 1.0];
    let result =
        rssn::kernels::rootfind::least_squares(residual, None, &[0.0], Default::default()).unwrap();
    assert!(result.converged(), "{result:?}");
    assert!((result.x[0] + 1.0).abs() < 1e-6, "{result:?}");
}

#[test]
fn system_solvers_reject_malformed_values_and_handle_large_residuals() {
    use rssn::kernels::dense::Mat;
    use rssn::kernels::rootfind::RootError;
    use rssn::kernels::rootfind::broyden;
    use rssn::kernels::rootfind::newton_system;
    assert_eq!(
        newton_system(
            |_| vec![f64::NAN],
            None::<&fn(&[f64]) -> Mat>,
            &[0.0],
            1e-9,
            10
        )
        .unwrap_err(),
        RootError::Invalid
    );
    assert_eq!(
        broyden(|_| vec![], &[0.0], 1e-9, 10).unwrap_err(),
        RootError::Invalid
    );
    assert_eq!(
        newton_system(|_| vec![], None::<&fn(&[f64]) -> Mat>, &[], 1e-9, 10).unwrap_err(),
        RootError::Invalid
    );
    let changes_size = |x: &[f64]| {
        if x[0] == 0.0 {
            vec![1.0]
        } else {
            vec![]
        }
    };
    assert_eq!(
        newton_system(changes_size, None::<&fn(&[f64]) -> Mat>, &[0.0], 1e-9, 10).unwrap_err(),
        RootError::Invalid
    );
    assert_eq!(
        broyden(changes_size, &[0.0], 1e-9, 10).unwrap_err(),
        RootError::Invalid
    );
    let f = |x: &[f64]| vec![1e200 * (x[0] - 1.0)];
    let j = |_: &[f64]| Mat::from_rows(&[vec![1e200]]).unwrap();
    let result = newton_system(f, Some(&j), &[0.0], 1e-10, 10).unwrap();
    assert_eq!(result.x[0], 1.0);
}

#[test]
fn standalone_jacobian_reports_changing_residual_shape_without_panicking() {
    use rssn::kernels::rootfind::RootError;
    use rssn::kernels::rootfind::numerical_jacobian;
    use rssn::kernels::rootfind::try_numerical_jacobian;
    let f = |x: &[f64]| {
        if x[0] == 0.0 {
            vec![1.0]
        } else {
            vec![]
        }
    };
    assert_eq!(
        try_numerical_jacobian(&f, &[0.0], &[1.0]).unwrap_err(),
        RootError::Invalid
    );
    let legacy = numerical_jacobian(&f, &[0.0], &[1.0]);
    assert_eq!((legacy.rows, legacy.cols), (1, 1));
    assert!(legacy.data[0].is_nan());
}

#[test]
fn solves_rosenbrock_residuals_with_an_analytic_jacobian() {
    use rssn::kernels::dense::Mat;
    use rssn::kernels::rootfind::LeastSquaresOptions;
    use rssn::kernels::rootfind::least_squares;
    let residual = |x: &[f64]| vec![10.0 * (x[1] - x[0] * x[0]), 1.0 - x[0]];
    let jacobian =
        |x: &[f64]| Mat::from_rows(&[vec![-20.0 * x[0], 10.0], vec![-1.0, 0.0]]).unwrap();
    let options = LeastSquaresOptions {
        max_iterations: 300,
        ..Default::default()
    };
    let result = least_squares(residual, Some(&jacobian), &[-1.2, 1.0], options).unwrap();
    assert!(result.converged(), "{result:?}");
    assert!(
        (result.x[0] - 1.0).abs() < 1e-8 && (result.x[1] - 1.0).abs() < 1e-8,
        "{result:?}"
    );
}

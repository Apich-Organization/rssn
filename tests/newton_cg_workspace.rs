use rssn::kernels::newton_cg::NewtonCgError;
use rssn::kernels::newton_cg::NewtonCgOptions;
use rssn::kernels::newton_cg::NewtonCgStatus;
use rssn::kernels::newton_cg::NewtonCgWorkspace;
use rssn::kernels::newton_cg::minimize;

#[test]
fn oversized_positive_or_negative_directional_curvature_reaches_the_boundary() {
    let options = NewtonCgOptions {
        gradient_tolerance: 0.0,
        initial_radius: 1.0,
        max_radius: 1.0,
        max_iterations: 1,
        ..Default::default()
    };
    let linear = 2.0_f64.powi(700);
    for quadratic in [1.0, -1.0] {
        let result = minimize(
            |x| 0.5 * quadratic * x[0] * x[0] - linear * x[0],
            |x, g| g[0] = quadratic * x[0] - linear,
            |_, direction, action| action[0] = quadratic * direction[0],
            &[0.0],
            &options,
        )
        .unwrap();
        assert_eq!(result.x[0].to_bits(), 1.0_f64.to_bits());
        assert_eq!(result.value.to_bits(), (-linear).to_bits());
        assert_eq!(result.boundary_steps, 1);
        assert_eq!(
            result.negative_curvature_steps,
            usize::from(quadratic < 0.0)
        );
        assert_eq!(result.accepted_steps, 1);
        assert_eq!(result.hessian_vector_evaluations, 2);
    }
}

#[test]
fn underflowed_directional_curvature_finds_the_representable_quadratic_minimum() {
    let coefficient = 2.0_f64.powi(-200);
    let linear = 2.0_f64.powi(-600);
    let result = minimize(
        |x| 0.5 * (coefficient * x[0]) * x[0] - linear * x[0],
        |x, g| g[0] = coefficient * x[0] - linear,
        |_, direction, action| action[0] = coefficient * direction[0],
        &[0.0],
        &NewtonCgOptions {
            gradient_tolerance: 0.0,
            max_iterations: 1,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(result.x[0].to_bits(), 2.0_f64.powi(-400).to_bits());
    assert_eq!(result.value.to_bits(), (-2.0_f64.powi(-1001)).to_bits());
    assert_eq!(result.status, NewtonCgStatus::FirstOrderConverged);
    assert_eq!(result.negative_curvature_steps, 0);
    assert_eq!(result.boundary_steps, 0);
    assert_eq!(result.hessian_vector_evaluations, 2);
}

#[test]
fn workspace_reuses_storage_across_dimensions_and_after_a_failed_action() {
    fn value(x: &[f64]) -> f64 {
        0.5 * x.iter().map(|v| v * v).sum::<f64>()
    }
    fn gradient(
        x: &[f64],
        g: &mut [f64],
    ) {
        g.copy_from_slice(x);
    }
    fn hessian(
        _: &[f64],
        direction: &[f64],
        action: &mut [f64],
    ) {
        action.copy_from_slice(direction);
    }

    let mut workspace = NewtonCgWorkspace::default();
    let options = NewtonCgOptions::default();
    let first = workspace
        .minimize(value, gradient, hessian, &[1.0, 2.0, 3.0], &options)
        .unwrap();
    let failure = workspace.minimize(
        value,
        gradient,
        |_, direction, action| action[0] = direction[0],
        &[1.0, 2.0],
        &options,
    );
    assert_eq!(failure.unwrap_err(), NewtonCgError::NonFiniteModel);
    for point in [&[3.0, -4.0][..], &[0.0][..], &[2.0, 5.5, 3.0, -1.0][..]] {
        let fresh = minimize(value, gradient, hessian, point, &options).unwrap();
        let reused = workspace
            .minimize(value, gradient, hessian, point, &options)
            .unwrap();
        assert_eq!(reused.status, NewtonCgStatus::FirstOrderConverged);
        assert_eq!(fresh.status, reused.status);
        assert_eq!(fresh.value.to_bits(), reused.value.to_bits());
        assert_eq!(
            fresh.gradient_norm.to_bits(),
            reused.gradient_norm.to_bits()
        );
        assert_eq!(
            fresh.hessian_vector_evaluations,
            reused.hessian_vector_evaluations
        );
        assert_eq!(
            fresh.x.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
            reused.x.iter().map(|x| x.to_bits()).collect::<Vec<_>>()
        );
    }
    assert_eq!(first.x.len(), 3);
    assert_eq!(first.status, NewtonCgStatus::FirstOrderConverged);
}

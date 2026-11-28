use std::cell::Cell;

use rssn::kernels::ode_geom::GaussLegendreWorkspace;
use rssn::kernels::ode_geom::OdeGeomError;
use rssn::kernels::ode_geom::gauss_legendre_step;

#[test]
fn implicit_fallback_samples_the_actual_gauss_nodes() {
    for (t0, h, sign) in [(0.0, 1.0, -1.0), (1.0, -1.0, 1.0)] {
        let mut times = Vec::new();
        let (next, count) = gauss_legendre_step(
            &mut |time, state| {
                times.push(time);
                vec![sign * 48.0 * (time - 0.5).powi(2) * state[0]]
            },
            t0,
            &[1.0],
            h,
            2,
        )
        .unwrap();
        // Both collocation nodes have z = h * lambda = -4, so the
        // fourth-order [2/2] Padé factor is (1-2+4/3)/(1+2+4/3).
        assert!((next[0] - 1.0 / 13.0).abs() < 2e-14, "{next:?}");
        assert_eq!(count, times.len());
        assert!(
            count > 1 + 12 * 2,
            "the implicit fallback must be exercised"
        );
        let offset = 3.0_f64.sqrt() / 6.0;
        let nodes = [t0 + h * (0.5 - offset), t0 + h * (0.5 + offset)];
        assert!(
            times[1..]
                .iter()
                .all(|time| { nodes.iter().any(|node| (node - time).abs() < f64::EPSILON) })
        );
    }
}

#[test]
fn different_stage_derivatives_match_the_direct_two_by_two_solve() {
    let r = 3.0_f64.sqrt() / 6.0;
    let (left, right) = (0.5 - r, 0.5 + r);
    let (l, q) = (-4.0 * (1.0 + left), -4.0 * (1.0 + right));
    // Independent inversion of the scalar coupled stage system
    // (I - diag(lambda) A) K = lambda.
    let (a, b, c, d) = (
        1.0 - 0.25 * l,
        -l * (0.25 - r),
        -q * (0.25 + r),
        1.0 - 0.25 * q,
    );
    let denominator = a * d - b * c;
    let expected = 1.0 + 0.5 * ((d * l - b * q) + (a * q - c * l)) / denominator;
    let (next, _) = gauss_legendre_step(
        &mut |t, y| vec![-4.0 * (1.0 + t) * y[0]],
        0.0,
        &[1.0],
        1.0,
        2,
    )
    .unwrap();
    assert!((next[0] - expected).abs() < 2e-14, "{next:?} != {expected}");
}

#[test]
fn homogeneous_linear_steps_resolve_the_same_pade_factor_across_scales() {
    let mut workspace = GaussLegendreWorkspace::default();
    for amplitude in [
        2.0_f64.powi(-700),
        1.0,
        2.0_f64.powi(700),
        -2.0_f64.powi(-700),
    ] {
        for (stages, expected) in [(1, -1.0 / 3.0), (2, 1.0 / 13.0), (3, 1.0 / 77.0)] {
            let count = Cell::new(0);
            let (next, evaluations) = workspace
                .step(
                    &mut |_, state| {
                        count.set(count.get() + 1);
                        state.iter().map(|value| -4.0 * value).collect()
                    },
                    0.0,
                    &[amplitude, amplitude * 0.5],
                    1.0,
                    stages,
                )
                .unwrap();
            assert_eq!(evaluations, count.get());
            assert!((next[0] / amplitude - expected).abs() < 2e-14);
            assert!((next[1] / amplitude - 0.5 * expected).abs() < 2e-14);
        }
    }
}

#[test]
fn harmonic_oscillator_keeps_second_fourth_and_sixth_order() {
    let mut workspace = GaussLegendreWorkspace::default();
    let end = 2.0_f64;
    for stages in 1..=3 {
        let mut errors = Vec::new();
        for steps in [2, 4, 8] {
            let h = end / f64::from(steps);
            let mut state = vec![1.0, 0.0];
            for i in 0..steps {
                state = workspace
                    .step(
                        &mut |_, state| vec![-state[1], state[0]],
                        f64::from(i) * h,
                        &state,
                        h,
                        stages,
                    )
                    .unwrap()
                    .0;
            }
            let error = (state[0] - end.cos()).hypot(state[1] - end.sin());
            errors.push(error);
        }
        let order = (errors[1] / errors[2]).log2();
        assert!(
            (order - (2 * stages) as f64).abs() < 0.2,
            "{order}: {errors:?}"
        );
    }
}

#[test]
fn invalid_inputs_and_staged_errors_leave_the_workspace_reusable() {
    let mut workspace = GaussLegendreWorkspace::default();
    for (t, y, h, stages) in [
        (0.0, f64::NAN, 1.0, 2),
        (f64::INFINITY, 1.0, 1.0, 2),
        (0.0, 1.0, f64::NAN, 2),
        (f64::MAX, 1.0, f64::MAX, 2),
        (0.0, 1.0, 1.0, 4),
    ] {
        assert!(matches!(
            workspace.step(
                &mut |_, _| panic!("invalid input reached RHS"),
                t,
                &[y],
                h,
                stages
            ),
            Err(OdeGeomError::Domain(_))
        ));
    }
    assert!(matches!(
        workspace.step(
            &mut |t, _| {
                if t > 0.0 {
                    vec![]
                } else {
                    vec![1.0; 2]
                }
            },
            0.0,
            &[1.0, 2.0],
            0.5,
            3
        ),
        Err(OdeGeomError::ShapeMismatch(_))
    ));
    let (next, count) = workspace
        .step(&mut |t, _| vec![t.powi(3)], 0.0, &[0.0], 1.0, 2)
        .unwrap();
    assert!(count >= 3);
    assert!((next[0] - 0.25).abs() < 2e-15);
}

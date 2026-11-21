//! Edge cases and pathwise checks for drift-implicit steps and the samplers.

use std::cell::Cell;

use rssn::kernels::random::Rng;
use rssn::kernels::sde::brownian_increments;
use rssn::kernels::sde_ext::DriftImplicitScheme;
use rssn::kernels::sde_ext::SdeExtError;
use rssn::kernels::sde_ext::drift_implicit_core;
use rssn::kernels::sde_ext::drift_implicit_sampled;
use rssn::kernels::sde_ext::drift_implicit_step_with_areas;
use rssn::kernels::sde_ext::sample_alpha_stable_with;
use rssn::kernels::sde_ext::sample_poisson_with;

fn zero(
    _t: f64,
    _state: &[f64],
    out: &mut [f64],
) {
    out.fill(0.0);
}

#[test]
fn sampled_paths_share_every_step_of_the_prescribed_path() {
    let drift = |t: f64, y: &[f64], out: &mut [f64]| out[0] = t - y[0];
    let diffusion = |_t, y: &[f64], out: &mut [f64]| out[0] = 0.5 * y[0];
    let increments = brownian_increments(&mut Rng::new(212), 1, 11, 0.125);
    for scheme in [
        DriftImplicitScheme::BackwardEuler,
        DriftImplicitScheme::Trapezoidal,
        DriftImplicitScheme::Theta(0.0),
        DriftImplicitScheme::Theta(0.7),
        DriftImplicitScheme::DriftImplicitMilstein,
        DriftImplicitScheme::SplitStepBackwardEuler,
    ] {
        let full = drift_implicit_core(
            &drift,
            &diffusion,
            1,
            0.5,
            0.125,
            &[1.0],
            &increments,
            scheme,
        )
        .unwrap();
        for stride in [1, 3, 11, usize::MAX] {
            let sparse = drift_implicit_sampled(
                &drift,
                &diffusion,
                1,
                0.5,
                0.125,
                &[1.0],
                &increments,
                scheme,
                stride,
            )
            .unwrap();
            let expected: Vec<_> = full
                .iter()
                .enumerate()
                .filter_map(|(step, state)| {
                    (step.is_multiple_of(stride) || step == increments.len())
                        .then_some(state[0].to_bits())
                })
                .collect();
            assert_eq!(
                sparse
                    .iter()
                    .map(|state| state[0].to_bits())
                    .collect::<Vec<_>>(),
                expected
            );
        }
    }
}

#[test]
fn invalid_grid_shape_and_retention_fail_before_callbacks() {
    let unexpected =
        |_t, _y: &[f64], _out: &mut [f64]| panic!("invalid request reached a callback");
    let increments = [vec![0.0], vec![0.0]];
    let scheme = DriftImplicitScheme::BackwardEuler;
    for (t, dt, stride) in [
        (0.0, 1.0, 0),
        (f64::MAX, f64::MAX, 1),
        (2.0_f64.powi(53), 0.5, 1),
    ] {
        assert!(matches!(
            drift_implicit_sampled(
                &unexpected,
                &unexpected,
                1,
                t,
                dt,
                &[0.0],
                &increments,
                scheme,
                stride
            ),
            Err(SdeExtError::Invalid(_))
        ));
    }
    assert!(
        drift_implicit_core(
            &unexpected,
            &unexpected,
            1,
            0.0,
            1.0,
            &[0.0],
            &[vec![0.0], vec![f64::NAN]],
            scheme
        )
        .is_err()
    );
    assert!(
        drift_implicit_core(
            &unexpected,
            &unexpected,
            1,
            0.0,
            1.0,
            &[0.0],
            &[],
            DriftImplicitScheme::Theta(f64::NAN)
        )
        .is_err()
    );
    assert!(
        drift_implicit_core(
            &unexpected,
            &unexpected,
            4_096,
            0.0,
            1.0,
            &[0.0],
            &[vec![0.0; 4_096]],
            DriftImplicitScheme::DriftImplicitMilstein
        )
        .is_err()
    );
    assert_eq!(
        drift_implicit_sampled(
            &unexpected,
            &unexpected,
            1,
            f64::MAX,
            1.0,
            &[4.0],
            &[],
            scheme,
            2
        )
        .unwrap(),
        [vec![4.0]]
    );
    let wide = vec![1.0; 4_096];
    assert_eq!(
        drift_implicit_core(&unexpected, &unexpected, 1, 0.0, 1.0, &wide, &[], scheme).unwrap(),
        [wide]
    );
}

#[test]
fn implicit_tolerance_resolves_small_coordinates_beside_large_ones() {
    let drift = |_t, y: &[f64], out: &mut [f64]| out[0] = -y[0];
    let actual = drift_implicit_core(
        &drift,
        &zero,
        1,
        0.0,
        1.0,
        &[1e-190, 1e150],
        &[vec![0.0]],
        DriftImplicitScheme::BackwardEuler,
    )
    .unwrap();
    assert!((actual[1][0] / 5e-191 - 1.0).abs() < 1e-13);
    assert_eq!(actual[1][1], 1e150);
}

#[test]
fn finite_difference_stays_in_range_at_the_top_of_the_state_domain() {
    let drift = |_t, y: &[f64], out: &mut [f64]| {
        assert!(y[0].is_finite());
        out[0] = -1e305;
    };
    let path = drift_implicit_core(
        &drift,
        &zero,
        1,
        0.0,
        1.0,
        &[f64::MAX],
        &[vec![0.0]],
        DriftImplicitScheme::BackwardEuler,
    )
    .unwrap();
    assert_eq!(path[1][0], f64::MAX - 1e305);
}

#[test]
fn split_step_evaluates_diffusion_at_the_implicit_drift_state() {
    let drift = |_t, _y: &[f64], out: &mut [f64]| out[0] = 1.0;
    let diffusion = |_t, y: &[f64], out: &mut [f64]| out[0] = y[0];
    let increments = [vec![1.0]];
    for (scheme, expected) in [
        (DriftImplicitScheme::BackwardEuler, 5.0),
        (DriftImplicitScheme::SplitStepBackwardEuler, 6.0),
    ] {
        let path =
            drift_implicit_core(&drift, &diffusion, 1, 0.0, 1.0, &[2.0], &increments, scheme)
                .unwrap();
        assert_eq!(path[1][0], expected);
    }
}

#[test]
fn zero_diffusion_rows_start_empty_on_each_callback() {
    let diffusion = |t: f64, _y: &[f64], out: &mut [f64]| {
        if t == 0.0 {
            out[0] = 2.0;
        }
    };
    let states = drift_implicit_core(
        &zero,
        &diffusion,
        1,
        0.0,
        0.5,
        &[0.0],
        &[vec![1.0], vec![1.0]],
        DriftImplicitScheme::Theta(0.0),
    )
    .unwrap();
    assert_eq!(states, [vec![0.0], vec![2.0], vec![2.0]]);
}

#[test]
fn final_nonfinite_drift_cannot_be_reported_as_a_converged_root() {
    let calls = Cell::new(0);
    let previous = Cell::new(0.0);
    let drift = |_t, y: &[f64], out: &mut [f64]| {
        let call = calls.get();
        calls.set(call + 1);
        out[0] = match (call, call % 12) {
            | (300, _) => f64::NAN,
            | (_, 0) => {
                previous.set(y[0]);
                y[0] + 1.0
            },
            | (_, 1) => previous.get() + 1.0,
            | _ => y[0] + 1.0,
        };
    };
    let result = drift_implicit_core(
        &drift,
        &zero,
        1,
        0.0,
        1.0,
        &[0.0],
        &[vec![0.0]],
        DriftImplicitScheme::BackwardEuler,
    );
    assert_eq!(calls.get(), 301);
    assert!(matches!(result, Err(SdeExtError::NonFinite(_))));
}

#[test]
fn noncommutative_step_uses_the_supplied_order_of_ito_integrals() {
    let diffusion = |_t, y: &[f64], out: &mut [f64]| out.copy_from_slice(&[1.0, 0.0, 0.0, y[0]]);
    let result = drift_implicit_step_with_areas(
        &zero,
        &diffusion,
        2,
        0.0,
        0.2,
        &[2.0, 0.0],
        &[0.2, 0.3],
        &[0.0, 0.15, 0.01, 0.0],
        DriftImplicitScheme::DriftImplicitMilstein,
    )
    .unwrap();
    assert_eq!(result[0], 2.2);
    assert!((result[1] - 0.75).abs() < 1e-10);
}

#[test]
fn alpha_stable_uses_its_documented_draws_and_rejects_nonfinite_samples() {
    let unit_draws = Cell::new(0);
    let normal_draws = Cell::new(0);
    let mut unit = || {
        unit_draws.set(unit_draws.get() + 1);
        0.75
    };
    let mut normal = || {
        normal_draws.set(normal_draws.get() + 1);
        0.0
    };
    let cauchy = sample_alpha_stable_with(&mut unit, &mut normal, 1.0, 0.0, 3.0, 2.0).unwrap();
    assert!((cauchy - 5.0).abs() < 1e-14);
    assert_eq!(unit_draws.get(), 2);
    assert_eq!(normal_draws.get(), 0);
    let gaussian =
        sample_alpha_stable_with(&mut unit, &mut normal, 2.0, 0.0, f64::MAX, 4.0).unwrap();
    assert_eq!(gaussian, 4.0);
    assert_eq!(unit_draws.get(), 2);
    assert_eq!(normal_draws.get(), 1);
    assert!(matches!(
        sample_alpha_stable_with(&mut || f64::NAN, &mut normal, 1.5, 0.0, 1.0, 0.0),
        Err(SdeExtError::Invalid(_))
    ));
    assert!(matches!(
        sample_alpha_stable_with(&mut unit, &mut || 2.0, 2.0, 0.0, f64::MAX, 0.0),
        Err(SdeExtError::NonFinite(_))
    ));
}

#[test]
fn poisson_has_the_expected_moments_on_both_sides_of_the_sampler_cutoff() {
    for expected in [12.0, 30.0, 31.0, 800.0] {
        let mut rng = Rng::new(1_089);
        let mut unit = || rng.uniform();
        let mut normal = || panic!("Poisson should only consume uniform draws");
        let samples: Vec<_> = (0..8192)
            .map(|_| sample_poisson_with(&mut unit, &mut normal, expected) as f64)
            .collect();
        let mean = samples.iter().sum::<f64>() / samples.len() as f64;
        let variance = samples
            .iter()
            .map(|sample| (sample - mean).powi(2))
            .sum::<f64>()
            / samples.len() as f64;
        assert!(
            (mean - expected).abs() < 0.08 * expected.sqrt(),
            "{expected}: mean {mean}"
        );
        assert!(
            (variance / expected - 1.0).abs() < 0.07,
            "{expected}: variance {variance}"
        );
    }
    assert_eq!(
        sample_poisson_with(&mut || panic!(), &mut || panic!(), 0.0),
        0
    );
}

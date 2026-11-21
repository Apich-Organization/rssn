//! Tests for deterministic admission and sampled work limits in fixed jump diffusion.

use rssn::sim::stochastic::JumpDistribution;
use rssn::sim::stochastic::JumpSpec;
use rssn::sim::stochastic::Method;
use rssn::sim::stochastic::Options;
use rssn::sim::stochastic::SdeError;
use rssn::sim::stochastic::Structure;
use rssn::sim::stochastic::integrate_jump_diffusion;
use rssn::sim::stochastic::integrate_jump_diffusion_compensated;

fn fixed(
    dt: f64,
    steps: usize,
    structure: Structure,
) -> Options {
    Options::fixed(Method::EulerMaruyama, structure, 77, 0.0, dt, steps)
}

fn constant(intensity: f64) -> JumpSpec {
    JumpSpec {
        intensity,
        distribution: JumpDistribution::Constant(vec![1.0]),
        compensated: false,
    }
}

fn reject_before_callbacks(
    initial: &[f64],
    spec: &JumpSpec,
    options: &Options,
) -> SdeError {
    integrate_jump_diffusion(
        |_, _| panic!("drift called for an invalid request"),
        |_, _| panic!("diffusion called for an invalid request"),
        |_, _, _| panic!("jump called for an invalid request"),
        initial,
        spec,
        options,
    )
    .unwrap_err()
}

#[test]
fn nonfinite_poisson_mean_is_an_invalid_schedule_before_callbacks() {
    let initial = [0.0];
    let spec = constant(f64::MAX);
    let grid = Options::grid(
        Method::EulerMaruyama,
        Structure::Diagonal,
        77,
        vec![0.0, 2.0],
    );
    for options in [fixed(2.0, 1, Structure::Diagonal), grid] {
        let error = reject_before_callbacks(&initial, &spec, &options);
        assert!(matches!(&error, SdeError::Invalid(detail) if detail.contains("Poisson mean")));
        assert_eq!(error.code(), "E0802");
    }

    let error = integrate_jump_diffusion_compensated(
        |_, _| panic!("drift called for an invalid request"),
        |_, _| panic!("diffusion called for an invalid request"),
        |_, _, _| panic!("jump called for an invalid request"),
        |_, _| panic!("expected jump called for an invalid request"),
        &initial,
        &spec,
        &fixed(2.0, 1, Structure::Diagonal),
    )
    .unwrap_err();
    assert_eq!(error.code(), "E0802");
}

#[test]
fn finite_means_and_wide_marks_are_checked_against_the_full_simulation_budget() {
    let one_step = fixed(1.0, 1, Structure::Diagonal);
    let two_steps = fixed(1.0, 2, Structure::Diagonal);
    let very_wide = JumpSpec {
        intensity: 0.25,
        distribution: JumpDistribution::DoubleExponential {
            p_up: 0.5,
            eta_up: 1.0,
            eta_down: 1.0,
            dim: usize::MAX,
        },
        compensated: false,
    };
    let width_limited = JumpSpec {
        intensity: 600_000.0,
        distribution: JumpDistribution::Constant(vec![1.0; 16]),
        compensated: false,
    };
    for (spec, options) in [
        (constant(2.0_f64.powi(54)), &one_step),
        (constant(f64::MAX), &one_step),
        (constant(600_000.0), &two_steps),
        (very_wide, &one_step),
        (width_limited, &one_step),
    ] {
        let error = reject_before_callbacks(&[1.0], &spec, options);
        assert!(
            matches!(&error, SdeError::StepLimit(detail) if detail.contains("expected jump count"))
        );
        assert_eq!(error.code(), "E0807");
    }
    let state_limited = reject_before_callbacks(&[1.0; 16], &constant(600_000.0), &one_step);
    assert!(
        matches!(&state_limited, SdeError::StepLimit(detail) if detail.contains("expected jump count"))
    );
}

#[test]
fn fixed_path_and_brownian_storage_are_bounded_before_presampling() {
    let spec = constant(0.0);
    let too_many_steps = fixed(1.0, usize::MAX, Structure::Diagonal);
    let state_limited = fixed(1.0, 524_288, Structure::Diagonal);
    let noise_limited = fixed(
        1.0,
        40_000,
        Structure::General {
            noise: 257,
            commutative: true,
        },
    );
    let error = reject_before_callbacks(&[0.0], &spec, &too_many_steps);
    assert!(matches!(&error, SdeError::StepLimit(detail) if detail.contains("step count")));
    for (initial, options) in [
        (&[0.0; 16][..], &state_limited),
        (&[0.0][..], &noise_limited),
    ] {
        let error = reject_before_callbacks(initial, &spec, options);
        assert!(matches!(&error, SdeError::StepLimit(detail) if detail.contains("value budget")));
    }
}

#[test]
fn an_expected_count_at_the_limit_reaches_the_solver_before_sampling() {
    let callback_error = SdeError::Invalid("drift sentinel".to_owned());
    let error = integrate_jump_diffusion(
        |_, _| Err(callback_error.clone()),
        |_, _| Ok(vec![0.0]),
        |_, _, _| panic!("jump called after a drift error"),
        &[0.0],
        &constant(1_048_576.0),
        &fixed(1.0, 1, Structure::Diagonal),
    )
    .unwrap_err();
    assert_eq!(error, callback_error);
}

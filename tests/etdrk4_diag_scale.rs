//! Scale and callback regressions for diagonal ETDRK4.

use std::cell::Cell;

use rssn::kernels::ode::OdeStepMode;
use rssn::kernels::ode_adaptive::OdeError;
use rssn::kernels::ode_adaptive::OdeOptions;
use rssn::kernels::ode_ext::etdrk4_diag;
use rssn::kernels::ode_ext::etdrk4_diag_adaptive;
use rssn::kernels::ode_ext::etdrk4_diag_with;

#[test]
fn dissipative_linear_rates_keep_representable_forcing() {
    let result = etdrk4_diag(
        &[-1e120, -1e308, -1e308],
        |_, _, rate| rate.copy_from_slice(&[1.0, 1.0, 1e308]),
        0.0,
        2.0,
        &[0.25; 3],
        1,
    )
    .unwrap();
    // The exponential transients are below binary64 range for all three modes.
    for (actual, expected) in result[1].iter().zip([1e-120, 1e-308, 1.0]) {
        assert!(
            (actual / expected - 1.0).abs() < 2e-15,
            "{actual:e} vs {expected:e}"
        );
    }
}

#[test]
fn a_long_step_retains_the_zero_mode_and_the_real_weight_transition() {
    let h = 1e308;
    let linear = [0.0, -1e-322, -1e-308, 1e-308, -2.0001e-308, 2.0001e-308];
    let result = etdrk4_diag(
        &linear,
        |_, _, rate| rate.fill(1e-308),
        0.0,
        h,
        &[0.0; 6],
        1,
    )
    .unwrap();
    for (&lambda, &actual) in linear.iter().zip(&result[1]) {
        let z = h * lambda;
        let want = if z.to_bits() == 0 {
            h * 1e-308
        } else {
            (h * 1e-308) * z.exp_m1() / z
        };
        assert!(
            (actual / want - 1.0).abs() < 3e-15,
            "z={z:e}: {actual:e} vs {want:e}"
        );
    }
}

#[test]
fn the_midpoint_weight_survives_when_the_scaled_eigenvalue_overflows() {
    let h: f64 = 1e308;
    let midpoint = (0.5 * h).to_bits();
    let result = etdrk4_diag(
        &[-2.0],
        |time, _, rate| {
            rate[0] = if time.to_bits() == midpoint {
                1.0
            } else {
                0.0
            };
        },
        0.0,
        h,
        &[0.0],
        1,
    )
    .unwrap();
    // Only the two midpoint samples are nonzero. 4 f_2 is asymptotic to
    // 4 / (h lambda^2), with the omitted relative correction below one ulp.
    assert!((result[1][0] / h.recip() - 1.0).abs() < 3e-15, "{result:?}");
}

#[test]
fn nonlinear_fixed_and_adaptive_paths_follow_the_same_analytic_solution() {
    let exact = |t: f64| 1.0 / ((2.0 - 1.0 / 3.0) * (3.0 * t).exp() + 1.0 / 3.0);
    let evals = Cell::new(0);
    let nonlinear = |_: f64, state: &[f64], rate: &mut [f64]| {
        evals.set(evals.get() + 1);
        rate[0] = state[0] * state[0];
    };
    let fixed = etdrk4_diag_with(
        &[-3.0],
        nonlinear,
        1.0,
        0.0,
        &[exact(1.0)],
        OdeStepMode::fixed(1.0 / 128.0),
    )
    .unwrap();
    assert_eq!(fixed.t.len(), 129);
    assert_eq!(fixed.nfev, 4 * 128);
    assert_eq!(fixed.nfev, evals.get());
    assert!((fixed.y.last().unwrap()[0] - exact(0.0)).abs() < 3e-9);

    evals.set(0);
    let options = OdeOptions {
        rtol: 1e-10,
        atol: 1e-12,
        h0: Some(0.125),
        hmax: 0.125,
        ..OdeOptions::default()
    };
    let adaptive = etdrk4_diag_adaptive(&[-3.0], nonlinear, 0.0, 1.0, &[0.5], &options).unwrap();
    assert!(adaptive.nfev > 12);
    assert_eq!(adaptive.nfev, evals.get());
    for (t, state) in adaptive.t.iter().zip(&adaptive.y) {
        assert!(
            (state[0] - exact(*t)).abs() < 2e-9,
            "state {state:?} at {t}"
        );
    }
}

#[test]
fn invalid_ranges_do_not_call_the_nonlinearity() {
    let never = |_: f64, _: &[f64], _: &mut [f64]| panic!("unexpected callback");
    assert_eq!(
        etdrk4_diag(&[0.0], never, 0.0, 1.0, &[0.0], usize::MAX),
        Err(OdeError::StepLimit)
    );
    assert_eq!(
        etdrk4_diag(&[f64::NAN], never, 0.0, 1.0, &[0.0], 1),
        Err(OdeError::Invalid)
    );
    assert!(matches!(
        etdrk4_diag_adaptive(
            &[f64::INFINITY],
            never,
            0.0,
            1.0,
            &[0.0],
            &OdeOptions::default()
        ),
        Err(OdeError::Invalid)
    ));
}

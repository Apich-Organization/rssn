//! Independent reference tests for prepared Levin quadrature.

use num_complex::Complex64;
use rssn::kernels::oscillatory::LevinPlan;
use rssn::kernels::oscillatory::Oscillatory;
use rssn::kernels::oscillatory::levin;
use rssn::kernels::oscillatory::levin_adaptive;
use std::cell::Cell;

fn complex(value: Oscillatory) -> Complex64 {
    Complex64::new(value.re, value.im)
}

fn exponential_integral(w: f64) -> Complex64 {
    if w == 0.0 {
        return Complex64::ONE;
    }
    let half_sine = (0.5 * w).sin();
    Complex64::new(w.sin() / w, half_sine * (2.0 * half_sine / w))
}

#[test]
fn zero_and_small_frequency_join_the_high_frequency_rule() {
    for n in [2, 4, 8, 16] {
        for w in [
            0.0, 1e-300, 1e-14, 1e-8, 1e-3, 0.07, 0.2, 0.5, 1.0, 2.0, 4.0, 10.0, 50.0,
        ] {
            let actual = levin(&|_| 1.0, &|x| x, &|_| 1.0, 0.0, 1.0, w, n).unwrap();
            let expected = exponential_integral(w);
            assert!(
                (complex(actual) - expected).norm() < 2e-12,
                "n={n} w={w}: {actual:?} vs {expected:?}"
            );
            if w > 0.0 && w < 1e-7 {
                assert!((actual.im / expected.im - 1.0).abs() < 1e-13, "n={n} w={w}");
            }
        }
    }
    let adaptive = levin_adaptive(
        &|x| x.powi(4),
        &|_| unreachable!(),
        &|_| unreachable!(),
        -1.0,
        2.0,
        0.0,
        1e-12,
        32,
    )
    .unwrap();
    assert!((adaptive.re - 6.6).abs() < 1e-13);
    assert_eq!(adaptive.im, 0.0);
    assert!(adaptive.error < 1e-12);
}

#[test]
fn the_prepared_rule_reuses_phase_work_and_preserves_orientation() {
    let phase_calls = Cell::new(0);
    let derivative_calls = Cell::new(0);
    let g = |x: f64| {
        phase_calls.set(phase_calls.get() + 1);
        x * x
    };
    let dg = |x: f64| {
        derivative_calls.set(derivative_calls.get() + 1);
        2.0 * x
    };
    let w = 1.0e5;
    let plan = LevinPlan::new(&g, &dg, 1.0, 2.0, w, 12).unwrap();
    assert_eq!(derivative_calls.get(), 12);
    assert_eq!(phase_calls.get(), 2);
    assert_eq!(plan.nodes().len(), 12);
    assert_eq!(plan.nodes()[0], 2.0);
    assert_eq!(plan.nodes()[11], 1.0);

    let left = Complex64::from_polar(1.0, w);
    let right = Complex64::from_polar(1.0, 4.0 * w);
    let expected = (right - left) / Complex64::new(0.0, w);
    for scale in [-2.25, 0.0, 0.25, 1.0, 3.0] {
        let values: Vec<_> = plan.nodes().iter().map(|x| 2.0 * scale * x).collect();
        let actual = plan.integrate(&|x| 2.0 * scale * x).unwrap();
        assert!((actual.re - plan.integrate_samples(&values).unwrap().re).abs() < 1e-24);
        assert!((complex(actual) - expected * scale).norm() < 5e-19);
        assert!(actual.error.is_nan());
    }
    assert_eq!((phase_calls.get(), derivative_calls.get()), (2, 12));
    let backwards = LevinPlan::new(&|x| x * x, &|x| 2.0 * x, 2.0, 1.0, w, 12).unwrap();
    assert!((complex(backwards.integrate(&|x| 2.0 * x).unwrap()) + expected).norm() < 5e-19);
    let negative = LevinPlan::new(&|x| x * x, &|x| 2.0 * x, 1.0, 2.0, -w, 12).unwrap();
    assert!((complex(negative.integrate(&|x| 2.0 * x).unwrap()) - expected.conj()).norm() < 5e-19);
}

#[test]
fn width_and_phase_derivative_products_handle_extreme_scales() {
    let huge = 2.0_f64.powi(1023);
    let reciprocal_huge = f64::from_bits(1 << 51);
    let constant = |_: f64| unreachable!();
    let large = LevinPlan::new(&constant, &constant, -huge, huge, 0.0, 6).unwrap();
    assert!(large.nodes().iter().all(|x| x.is_finite()));
    let result = large.integrate(&|_| reciprocal_huge).unwrap();
    assert!((result.re - 2.0).abs() < 4e-14);
    assert_eq!(result.im, 0.0);

    let small = LevinPlan::new(&constant, &constant, 0.0, f64::from_bits(1), 0.0, 6).unwrap();
    let result = small.integrate(&|_| huge).unwrap();
    assert!((result.re / 2.0_f64.powi(-51) - 1.0).abs() < 2e-14);
    assert_eq!(result.im, 0.0);

    let (interval_end, derivative, frequency) =
        (2.0_f64.powi(-700), 2.0_f64.powi(150), 2.0_f64.powi(900));
    assert!((derivative * frequency).is_infinite());
    let angle = 2.0_f64.powi(350);
    let plan = LevinPlan::new(
        &|x| derivative * x,
        &|_| derivative,
        0.0,
        interval_end,
        frequency,
        6,
    )
    .unwrap();
    let result = plan.integrate(&|_| 2.0_f64.powi(1000)).unwrap();
    let expected = Complex64::new(angle.sin(), 1.0 - angle.cos());
    assert!((complex(result) * 2.0_f64.powi(50) - expected).norm() < 3e-14);
}

#[test]
fn invalid_samples_and_phases_are_reported_and_empty_intervals_need_no_samples() {
    let id = |x| x;
    let one = |_| 1.0;
    for n in [0, 1, usize::MAX] {
        assert!(LevinPlan::new(&id, &one, 0.0, 1.0, 1.0, n).is_err());
    }
    assert!(LevinPlan::new(&id, &one, 0.0, 1.0, f64::NAN, 4).is_err());
    assert!(LevinPlan::new(&id, &|_| f64::NAN, 0.0, 1.0, 1.0, 4).is_err());
    assert!(LevinPlan::new(&|_| f64::MAX, &|_| 0.0, 0.0, 1.0, 2.0, 4).is_err());
    assert!(LevinPlan::new(&id, &|_| f64::MAX, 0.0, 1.0, f64::MAX, 4).is_err());
    let plan = LevinPlan::new(&id, &one, 0.0, 4.0, 0.0, 4).unwrap();
    assert!(plan.integrate(&|_| f64::MAX).is_err());
    assert!(plan.integrate(&|_| f64::INFINITY).is_err());
    assert!(plan.integrate_samples(&[1.0; 3]).is_err());
    let empty = LevinPlan::new(&|_| unreachable!(), &|_| unreachable!(), 3.0, 3.0, 1.0, 4).unwrap();
    assert!(empty.nodes().is_empty());
    assert_eq!(
        complex(empty.integrate(&|_| unreachable!()).unwrap()),
        Complex64::ZERO
    );
}

#[test]
fn cancelled_scopes_stop_preparation_and_repeated_amplitudes_before_callbacks() {
    use std::sync::Arc;
    use std::sync::atomic::AtomicBool;

    let ready = LevinPlan::new(&|x| x, &|_| 1.0, 0.0, 1.0, 100.0, 8).unwrap();
    let _guard = rssn::cancellation::install(Arc::new(AtomicBool::new(true)));
    assert!(
        ready
            .integrate(&|_| unreachable!())
            .unwrap_err()
            .contains("cancelled")
    );
    assert!(
        LevinPlan::new(&|_| unreachable!(), &|_| unreachable!(), 0.0, 1.0, 1.0, 8)
            .unwrap_err()
            .contains("cancelled")
    );
}

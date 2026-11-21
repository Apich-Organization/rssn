//! Analytic references for diagonal finite differences at binary64 limits.

use std::cell::Cell;

use rssn::kernels::calculus::hessian;
use rssn::kernels::calculus::hessian_diagonal;
use rssn::kernels::calculus::hessian_diagonal_with_steps;
use rssn::kernels::calculus::laplacian;
use rssn::kernels::numeric_common::NumericError;

fn assert_relative(
    actual: f64,
    expected: f64,
    tolerance: f64,
) {
    assert!(
        actual.is_finite() && (actual / expected - 1.0).abs() <= tolerance,
        "{actual:e} != {expected:e}",
    );
}

#[test]
fn diagonal_retains_constant_and_small_curvature_when_scales_are_large() {
    assert_eq!(hessian(|_| f64::MAX, &[1.0]), [0.0]);
    assert_eq!(hessian_diagonal(|_| -f64::MAX, &[1e300]), [0.0]);
    let center = 1e300;
    let f = |point: &[f64]| ((point[0] - center) * 1e-150).powi(2);
    assert_relative(hessian(f, &[center])[0], 2e-300, 4e-15);
    assert_relative(hessian_diagonal(f, &[center])[0], 2e-300, 4e-15);
    assert_relative(laplacian(f, &[center]), 2e-300, 4e-15);
}

#[test]
fn scaled_second_difference_handles_opposite_maximum_ordinates() {
    let f = |point: &[f64]| {
        if point[0] == 0.0 {
            -f64::MAX
        } else {
            f64::MAX
        }
    };
    let result = hessian_diagonal_with_steps(f, &[0.0], &[4.0]).unwrap();
    assert_eq!(result, [f64::MAX * 0.25]);
}

#[test]
fn explicit_steps_resolve_unequal_adjacent_floats_and_a_tiny_axis() {
    let base = [2.0_f64.powi(52), 1e-150];
    let f = |point: &[f64]| {
        let dx = point[0] - base[0];
        let dy = (point[1] - base[1]) * 1e150;
        dx * dx + 7.0 * dx + dy * dy
    };
    let h = hessian_diagonal_with_steps(f, &base, &[1e-20, 1e-155]).unwrap();
    assert_eq!(h[0], 2.0);
    assert_relative(h[1], 2e300, 8e-15);
}

#[test]
fn inward_stencil_is_exact_for_a_local_cubic_at_finite_limits() {
    let unit = f64::MAX - f64::MAX.next_down();
    let amplitude = 2.0_f64.powi(1000);
    for boundary in [-f64::MAX, f64::MAX] {
        let calls = Cell::new(0);
        let f = |point: &[f64]| {
            assert!(point[0].is_finite());
            calls.set(calls.get() + 1);
            let local = (point[0] - boundary) / unit;
            amplitude * (local * local + local * local * local)
        };
        let result = hessian_diagonal_with_steps(f, &[boundary], &[1.0]).unwrap();
        assert_eq!(calls.get(), 4);
        assert_relative(result[0], 2.0_f64.powi(-941), 4e-14);
    }
}

#[test]
fn laplacian_uses_only_axis_samples() {
    const DIMENSION: usize = 32;
    let calls = Cell::new(0);
    let f = |point: &[f64]| {
        calls.set(calls.get() + 1);
        assert!(point.iter().filter(|&&x| x != 0.0).count() <= 1);
        point.iter().map(|x| x * x).sum()
    };
    let point = [0.0; DIMENSION];
    let lap = laplacian(f, &point);
    assert_eq!(calls.get(), 2 * DIMENSION + 1);
    assert_relative(lap, 64.0, 4e-15);
    calls.set(0);
    let diagonal = hessian_diagonal(f, &point);
    assert_eq!(calls.get(), 2 * DIMENSION + 1);
    assert_eq!(lap, diagonal.into_iter().sum::<f64>());
}

#[test]
fn explicit_step_validation_happens_before_function_sampling() {
    let f = |_: &[f64]| -> f64 { panic!("invalid and empty inputs must not use the callback") };
    assert!(matches!(
        hessian_diagonal_with_steps(f, &[0.0], &[]),
        Err(NumericError::ShapeMismatch(_))
    ));
    for step in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            hessian_diagonal_with_steps(f, &[1.0], &[step]),
            Err(NumericError::Domain(_))
        ));
    }
    assert!(matches!(
        hessian_diagonal_with_steps(f, &[f64::INFINITY], &[1.0]),
        Err(NumericError::Domain(_))
    ));
    assert!(hessian_diagonal(f, &[]).is_empty());
    assert_eq!(laplacian(f, &[]), 0.0);
}

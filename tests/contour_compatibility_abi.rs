//! Shared contour panels preserve estimates, sample counts, and raw call behavior.
use num_complex::Complex64;
use rssn::kernels::contour_deform::ContourSegment;
use rssn::kernels::contour_deform::cauchy_principal_value;
use rssn::kernels::contour_deform::integrate_segment;
use rssn::kernels::contour_deform::integrate_segment_checked;
use rssn::kernels::contour_deform::residue_at_infinity_numeric;
use std::cell::Cell;
use std::cell::RefCell;
use std::collections::HashSet;

#[test]
fn legacy_line_and_fixed_residue_keep_values_and_minimum_sample_count() {
    let line = ContourSegment::Line {
        start: Complex64::ZERO,
        end: Complex64::ONE,
    };
    let calls = Cell::new(0);
    let integral = integrate_segment(
        &|z| {
            calls.set(calls.get() + 1);
            z * z
        },
        line,
        1e-11,
    );
    assert!((integral.re - 1.0 / 3.0).abs() < 1e-14);
    assert_eq!(calls.replace(0), 15);
    let residue = residue_at_infinity_numeric(
        |z| {
            calls.set(calls.get() + 1);
            Complex64::ONE / z
        },
        2.0,
        16,
    );
    assert!((residue + Complex64::ONE).norm() < 1e-13);
    assert_eq!(calls.get(), 32);
}

#[test]
fn shared_first_panel_is_reused_if_adaptive_refinement_is_needed() {
    let sampled = RefCell::new(HashSet::new());
    let line = ContourSegment::Line {
        start: Complex64::ZERO,
        end: Complex64::ONE,
    };
    let result = integrate_segment_checked(
        |z| {
            assert!(sampled.borrow_mut().insert(z.re.to_bits()));
            let wave = (Complex64::I * 90.0 * z).exp();
            wave + Complex64::new((-(90.0 * (z.re - 0.29)).powi(2)).exp(), 0.0)
        },
        line,
        1e-10,
        128,
    )
    .unwrap();
    assert!(result.converged && result.evaluations > 15);
    assert_eq!(result.evaluations, sampled.into_inner().len());
}

#[test]
fn legacy_principal_value_keeps_filtered_poles_and_empty_reverse_interval() {
    let f = |x| Complex64::from(1.0 / x);
    let result = cauchy_principal_value(f, -1.0, 2.0, &[0.0, 0.0, 4.0], 1e-11);
    assert!((result.re - 2.0_f64.ln()).abs() < 1e-11);
    assert_eq!(
        cauchy_principal_value(f, 2.0, -1.0, &[0.0], 1e-11),
        Complex64::ZERO
    );
}

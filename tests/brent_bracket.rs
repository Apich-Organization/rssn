//! Analytic and representability references for Brent's safeguarded search.

use rssn::kernels::rootfind::RootError;
use rssn::kernels::rootfind::brent;
use rssn::kernels::rootfind::brent_with_bracket;
use std::cell::Cell;

#[test]
fn finite_normalized_secant_finds_a_root_when_the_ordinary_formula_overflows() {
    let scale = 2.0_f64.powi(1022);
    let calls = Cell::new(0);
    let result = brent_with_bracket(
        |x| {
            calls.set(calls.get() + 1);
            x * scale
        },
        -3.0,
        1.0,
        1e-20,
        1,
    )
    .unwrap();
    assert_eq!(result.root, 0.0);
    assert_eq!(result.bracket, [0.0; 2]);
    assert_eq!(result.evaluations, 3);
    assert_eq!(result.evaluations, calls.get());
    assert_eq!(brent(|x| scale * x, 1.0, -3.0, 1e-20, 1), Ok(0.0));
}

#[test]
fn interpolation_has_the_same_budget_for_widely_scaled_quadratic_values() {
    let exact = 0.3;
    let calls = Cell::new(0);
    let mut expected_evaluations = None;
    for exponent in [0, -1000, 1000] {
        let scale = 2.0_f64.powi(exponent);
        let f = |x| (x - exact) * (x + 1.1) * scale;
        let counted = |x| {
            calls.set(calls.get() + 1);
            f(x)
        };
        let result = brent_with_bracket(counted, 0.0, 1.0, 1e-13, 14).unwrap();
        assert!((result.root - exact).abs() < 1e-13);
        assert!(result.bracket[0] <= exact && exact <= result.bracket[1]);
        assert_eq!(calls.replace(0), result.evaluations);
        if let Some(reference) = expected_evaluations {
            assert_eq!(result.evaluations, reference);
        } else {
            expected_evaluations = Some(result.evaluations);
        }
    }
}

#[test]
fn opposite_extreme_endpoints_and_subnormal_displacements_are_representable() {
    let limit = f64::MAX;
    for (a, b) in [(-limit, limit), (limit, -limit)] {
        let result = brent_with_bracket(|x| x / limit - 0.25, a, b, 1e292, 5).unwrap();
        assert!((result.root / limit - 0.25).abs() <= f64::EPSILON);
        assert!(result.bracket[0] <= limit * 0.25 && limit * 0.25 <= result.bracket[1]);
        assert!(result.evaluations <= 6);
    }

    let shift = f64::from_bits(2);
    let near_zero =
        brent_with_bracket(|x| x + shift, -f64::MAX, 0.0, f64::from_bits(1), 1).unwrap();
    assert_eq!(near_zero.root.to_bits(), (-shift).to_bits());
    assert_eq!(near_zero.evaluations, 3);
}

#[test]
fn the_reported_bracket_is_sorted_and_preserves_the_sign_change() {
    let f = |x: f64| x * x - 2.0;
    let result = brent_with_bracket(f, 4.0, 1.0, 1e-10, 50).unwrap();
    let exact = 2.0_f64.sqrt();
    let [lo, hi] = result.bracket;
    assert!(lo <= exact && exact <= hi);
    assert!(lo <= result.root && result.root <= hi);
    assert!(lo < hi);
    assert!(f(lo) < 0.0 && f(hi) > 0.0);
    assert!(hi - lo < 1e-10);

    let a = 1.0_f64;
    let b = a.next_up();
    let calls = Cell::new(0);
    let result = brent_with_bracket(
        |x| {
            calls.set(calls.get() + 1);
            if x == a { -0.5 } else { 0.125 }
        },
        b,
        a,
        f64::from_bits(1),
        1,
    )
    .unwrap();
    assert_eq!(result.root, b);
    assert_eq!(result.bracket, [a, b]);
    assert_eq!(result.evaluations, 2);
    assert_eq!(calls.get(), 2);
}

#[test]
fn errors_and_endpoint_evaluation_order_follow_the_existing_contract() {
    assert_eq!(
        brent_with_bracket(|_| unreachable!(), f64::NAN, 1.0, 1e-9, 10),
        Err(RootError::Invalid)
    );
    assert_eq!(
        brent_with_bracket(|_| unreachable!(), 0.0, 1.0, 0.0, 10),
        Err(RootError::Invalid)
    );
    assert_eq!(
        brent_with_bracket(|x| x * x + 1.0, -1.0, 1.0, 1e-9, 10),
        Err(RootError::NotBracketed)
    );
    assert_eq!(
        brent_with_bracket(|x| x, -1.0, 2.0, 1e9, 0),
        Err(RootError::NoConvergence)
    );
    let calls = Cell::new(0);
    let invalid = brent_with_bracket(
        |x| {
            calls.set(calls.get() + 1);
            if x == 0.0 { 0.0 } else { f64::NAN }
        },
        0.0,
        1.0,
        1e-9,
        0,
    );
    assert_eq!(invalid, Err(RootError::Invalid));
    assert_eq!(calls.get(), 2);
    let endpoint = brent_with_bracket(|x| x, 0.0, 1.0, 1e-9, 0).unwrap();
    assert_eq!(endpoint.bracket, [0.0; 2]);
    assert_eq!(endpoint.evaluations, 2);
}

//! Scale-aware tests of the same-sign fallback in nearby root search.

use rssn::kernels::rootfind::RootError;
use rssn::kernels::rootfind::root_near;

#[test]
fn double_roots_remain_resolvable_when_squares_underflow_or_overflow() {
    let expected = 2.0_f64.sqrt();
    for scale in [-1e200, -1e-200, 1e-200, 1.0, 1e200] {
        let mut first = Vec::new();
        let mut second = Vec::new();
        let a = root_near(
            &mut |x| {
                first.push(x.to_bits());
                ((x * x - 2.0).powi(2)) * scale
            },
            0.0,
        )
        .unwrap();
        let b = root_near(
            &mut |x| {
                second.push(x.to_bits());
                ((x * x - 2.0).powi(2)) * scale
            },
            0.0,
        )
        .unwrap();
        assert!((a.abs() - expected).abs() < 2e-14, "{scale}: {a}");
        assert_eq!(a.to_bits(), b.to_bits());
        assert_eq!(first, second);
    }
}

#[test]
fn the_reported_root_stays_within_a_subnormal_callbacks_zero_plateau() {
    let mut residual = |x: f64| (x * x - 2.0).powi(2) * 1e-310;
    let found = root_near(&mut residual, 0.0).unwrap();
    assert_eq!(residual(found).to_bits(), 0);
    assert!((found.abs() - 2.0_f64.sqrt()).abs() < 1e-7, "{found}");
}

#[test]
fn same_sign_constants_do_not_count_as_touching_roots() {
    for constant in [
        -1e200,
        -1e-200,
        -f64::from_bits(1),
        f64::from_bits(1),
        1e-200,
        1e200,
    ] {
        for guess in [0.0, 1e20] {
            let mut evaluations = 0;
            let result = root_near(
                &mut |_| {
                    evaluations += 1;
                    constant
                },
                guess,
            );
            assert_eq!(result, Err(RootError::NoConvergence), "{constant}, {guess}");
            // One Newton probe set, at most 30 sign brackets and five
            // golden attempts with no more than 203 evaluations apiece.
            assert!(evaluations <= 5 + 2 * 30 + 5 * 203);
        }
    }
}

#[test]
fn positive_offsets_are_rejected_using_the_local_residual_scale() {
    for scale in [1e-200, 1.0, 1e200] {
        for offset in [1e-24, 1e-12] {
            let mut f = |x: f64| ((x * x - 2.0).powi(2) + offset) * scale;
            assert_eq!(
                root_near(&mut f, 0.0),
                Err(RootError::NoConvergence),
                "{scale:e}, {offset:e}",
            );
        }
    }
}

#[test]
fn touching_search_expands_past_an_edge_and_requires_a_resolved_local_minimum() {
    let expected = 2.0_f64.sqrt();
    let mut tiny = |x: f64| (x - expected).powi(2) * 1e-200;
    let found = root_near(&mut tiny, 0.0).unwrap();
    assert!((found - expected).abs() < 1e-14, "{found}");

    let root = 1e12 + 0.25;
    let mut far = |x: f64| (x - root).powi(2) * 1e-200;
    let found = root_near(&mut far, 1e12).unwrap();
    assert!((found - root).abs() < 0.002, "{found}");

    let mut open_edge = |x: f64| {
        if x <= 0.25 {
            f64::NAN
        } else {
            (x - 0.25).powi(2) * 1e-200
        }
    };
    assert_eq!(
        root_near(&mut open_edge, 0.0),
        Err(RootError::NoConvergence)
    );
}

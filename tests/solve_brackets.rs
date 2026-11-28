//! IEEE boundaries and callback budgets for scalar interval root solving.

use rssn::kernels::rootfind::root_near;
use rssn::kernels::solve::solve_root;
use rssn::kernels::solve::solve_root_bisection;
use std::cell::Cell;
use std::cell::RefCell;

#[test]
fn bisection_brackets_opposite_subnormal_values_without_a_product() {
    let least = f64::from_bits(1);
    let residual = |x: f64| {
        if x.to_bits() == 0.25_f64.to_bits() {
            0.0
        } else {
            least.copysign(x - 0.25)
        }
    };
    for interval in [(0.0, 1.0), (1.0, 0.0), (-f64::MAX, f64::MAX)] {
        let root = solve_root_bisection(residual, interval, 0.0, 1_200).unwrap();
        assert_eq!(root.to_bits(), 0.25_f64.to_bits(), "{interval:?}");
    }
    for same_sign in [least, -least] {
        let error = solve_root_bisection(|_| same_sign, (0.0, 1.0), 0.0, 12).unwrap_err();
        assert!(error.contains("does not bracket a root"), "{error}");
    }
}

#[test]
fn bisection_does_not_discard_the_root_when_the_midpoint_product_underflows() {
    let residual = |x: f64| (x - 0.25) * if x <= 0.75 { 1e-200 } else { 1e200 };
    for interval in [(0.0, 1.0), (1.0, 0.0)] {
        let root = solve_root_bisection(residual, interval, 1e-250, 32).unwrap();
        assert_eq!(root.to_bits(), 0.25_f64.to_bits(), "{interval:?}");
    }
}

#[test]
fn bisection_checks_coordinates_tolerances_and_sampled_values() {
    let calls = Cell::new(0);
    let f = |x| {
        calls.set(calls.get() + 1);
        x
    };
    for endpoint in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(solve_root_bisection(f, (-1.0, endpoint), 1e-12, 8).is_err());
        assert!(solve_root(f, endpoint, 1e-12, 8).is_err());
    }
    for tolerance in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(solve_root_bisection(f, (-1.0, 1.0), tolerance, 8).is_err());
        assert!(solve_root(f, 0.0, tolerance, 8).is_err());
    }
    assert_eq!(calls.get(), 0);

    for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let bad_endpoint = |x: f64| if x > 0.0 { value } else { 0.0 };
        let bad_midpoint = |x: f64| {
            if x.abs().to_bits() == 0 {
                value
            } else {
                x
            }
        };
        for residual in [&bad_endpoint as &dyn Fn(f64) -> f64, &bad_midpoint] {
            let error = solve_root_bisection(residual, (-1.0, 1.0), 1e-12, 8).unwrap_err();
            assert_eq!(error, "Function value is not finite.");
        }
    }
}

#[test]
fn bisection_retains_endpoint_precedence_and_the_midpoint_budget() {
    let calls = RefCell::new(Vec::<f64>::new());
    let residual = |x| {
        calls.borrow_mut().push(x);
        x - 0.75
    };
    for (interval, expected_order) in [((0.0, 2.0), [0.0, 2.0]), ((2.0, 0.0), [2.0, 0.0])] {
        assert_eq!(
            solve_root_bisection(residual, interval, 0.0, 0)
                .unwrap()
                .to_bits(),
            1.0_f64.to_bits()
        );
        assert_eq!(calls.take(), expected_order);
        assert_eq!(
            solve_root_bisection(residual, interval, 0.0, 1)
                .unwrap()
                .to_bits(),
            0.5_f64.to_bits()
        );
        assert_eq!(calls.take(), [expected_order[0], expected_order[1], 1.0]);
    }
    let first = solve_root_bisection(|_| -0.0, (1.0, -1.0), 0.0, 0).unwrap();
    assert_eq!(first.to_bits(), 1.0_f64.to_bits());
    assert!(solve_root_bisection(|_| 1.0, (1.0, 1.0), 0.0, 20).is_err());
    let single = solve_root_bisection(|x| x - 1.0, (1.0, 1.0), 0.0, 20).unwrap();
    assert_eq!(single.to_bits(), 1.0_f64.to_bits());
}

#[test]
fn bisection_stops_without_reevaluating_adjacent_endpoints() {
    let a = 1.0_f64;
    let b = f64::from_bits(a.to_bits() + 1);
    for interval in [(a, b), (b, a)] {
        let calls = Cell::new(0);
        let residual = |x: f64| {
            calls.set(calls.get() + 1);
            if x.to_bits() == a.to_bits() {
                -1.0
            } else {
                1.0
            }
        };
        let result = solve_root_bisection(residual, interval, 0.0, 16_000).unwrap();
        assert_eq!(result.to_bits(), f64::midpoint(a, b).to_bits());
        assert_eq!(calls.get(), 2);
    }
}

#[test]
fn automatic_root_search_uses_the_next_valid_subnormal_bracket() {
    let least = f64::from_bits(1);
    let residual = |x: f64| {
        if x.to_bits() == (-2.0_f64).to_bits() {
            f64::NAN
        } else if x.to_bits() == 2.5_f64.to_bits() {
            0.0
        } else {
            least.copysign(x - 2.5)
        }
    };
    let result = solve_root(residual, 0.0, 0.0, 32).unwrap();
    assert_eq!(result.to_bits(), 2.5_f64.to_bits());
}

#[test]
fn automatic_root_search_retries_after_a_nonfinite_derivative_sample() {
    let probed = Cell::new(false);
    let residual = |x: f64| {
        assert!(x.is_finite());
        if x.to_bits() == 1e-6_f64.to_bits() {
            probed.set(true);
            f64::INFINITY
        } else {
            x - 0.25
        }
    };
    let root = solve_root(residual, 0.0, 1e-12, 40).unwrap();
    assert!(probed.get());
    assert!((root - 0.25).abs() < 1e-12, "{root}");

    let no_root = |x| {
        if x > 0.0 {
            f64::INFINITY
        } else {
            1.0
        }
    };
    assert!(solve_root(no_root, 0.0, 1e-12, 40).is_err());
}

#[test]
fn nearby_root_search_shares_the_bracket_for_tiny_measured_residual_scales() {
    let root = root_near(&mut |x| (x - 2.5) * 1e-200, 0.0).unwrap();
    assert!((root - 2.5).abs() <= 3e-14, "{root}");
}

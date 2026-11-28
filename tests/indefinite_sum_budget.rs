use std::cell::Cell;
use std::cell::RefCell;

use rssn::kernels::indefinite_sum::IndefiniteSumConfig;
use rssn::kernels::indefinite_sum::indefinite_product_with_budget;
use rssn::kernels::indefinite_sum::indefinite_sum;
use rssn::kernels::indefinite_sum::indefinite_sum_with_budget;

#[test]
fn invalid_coordinates_and_unrequestable_spans_fail_before_callbacks() {
    let unused = |_: f64| -> f64 { panic!("validation should precede callbacks") };
    let default = IndefiniteSumConfig::default();
    for x in [f64::INFINITY, f64::NEG_INFINITY, f64::NAN, 2.0_f64.powi(60)] {
        assert!(indefinite_sum_with_budget(unused, x, &default, usize::MAX).is_err());
    }
    for h in [f64::INFINITY, f64::NEG_INFINITY, f64::NAN] {
        let cfg = IndefiniteSumConfig { h, ..default.clone() };
        assert!(indefinite_sum_with_budget(unused, 0.0, &cfg, 10).is_err());
    }
    for step in [0.0, -0.0, f64::INFINITY, f64::NAN, f64::from_bits(1)] {
        let cfg = IndefiniteSumConfig {
            step,
            ..default.clone()
        };
        assert!(indefinite_sum_with_budget(unused, 1.0, &cfg, 10).is_err());
    }
    for x in [5_000_001.0, -5_000_001.0, 100_000_000.25] {
        assert!(
            indefinite_sum(unused, x, &default)
                .unwrap_err()
                .contains("budget")
        );
    }
    let cfg = IndefiniteSumConfig {
        h: -2.0_f64.powi(971),
        step: f64::MAX,
        ..default
    };
    // A nearly integral backward offset rounds to -1. Its physical lattice
    // sample is below -MAX even though the requested endpoint is finite.
    assert!(
        indefinite_sum_with_budget(unused, -f64::MAX, &cfg, 1)
            .unwrap_err()
            .contains("lattice sample")
    );
}

#[test]
fn exact_lattice_uses_half_open_endpoints_and_counts_each_term() {
    let points = RefCell::new(Vec::new());
    let f = |coordinate| {
        points.borrow_mut().push(coordinate);
        coordinate
    };
    let cfg = IndefiniteSumConfig {
        h: 1.0,
        step: 0.5,
        ..Default::default()
    };
    assert!(indefinite_sum_with_budget(f, 3.0, &cfg, 3).is_err());
    assert!(points.borrow().is_empty());
    let total = indefinite_sum_with_budget(f, 3.0, &cfg, 4).unwrap();
    assert_eq!(total.to_bits(), 7.0_f64.to_bits());
    assert_eq!(*points.borrow(), [1.0, 1.5, 2.0, 2.5]);

    points.borrow_mut().clear();
    let reverse = IndefiniteSumConfig {
        h: 2.0,
        step: -0.5,
        ..Default::default()
    };
    let total = indefinite_sum_with_budget(f, 3.0, &reverse, 2).unwrap();
    assert_eq!(total.to_bits(), (-5.5_f64).to_bits());
    assert_eq!(*points.borrow(), [3.0, 2.5]);

    let unused = |_: f64| -> f64 { panic!("normalization point needs zero evaluations") };
    let normalized = indefinite_sum_with_budget(unused, 2.0, &reverse, 0).unwrap();
    assert!(normalized.abs() < f64::EPSILON);
    let product = indefinite_product_with_budget(unused, 2.0, &reverse, 0).unwrap();
    assert_eq!(product.to_bits(), 1.0_f64.to_bits());
}

#[test]
fn scaled_finite_span_recovers_and_nonfinite_probes_are_never_passed_through() {
    let cfg = IndefiniteSumConfig {
        h: -f64::MAX,
        step: f64::MAX,
        ..Default::default()
    };
    let points = RefCell::new(Vec::new());
    let result = indefinite_sum_with_budget(
        |coordinate| {
            points.borrow_mut().push(coordinate);
            1.0
        },
        f64::MAX,
        &cfg,
        2,
    )
    .unwrap();
    assert_eq!(result.to_bits(), 2.0_f64.to_bits());
    assert_eq!(*points.borrow(), [-f64::MAX, 0.0]);

    let cfg = IndefiniteSumConfig {
        step: f64::MAX,
        ..Default::default()
    };
    let result = indefinite_sum(
        |coordinate| {
            assert!(coordinate.is_finite());
            1.0
        },
        f64::MAX * 0.5,
        &cfg,
    )
    .unwrap();
    assert!((result - 0.5).abs() < 1e-12);
}

#[test]
fn fractional_budget_reserves_a_complete_strategy_and_charges_the_fallback() {
    let cfg = IndefiniteSumConfig::default();
    let count = Cell::new(0);
    let evaluated = |coordinate: f64| {
        count.set(count.get() + 1);
        coordinate.exp2()
    };
    let x = 2.7_f64;
    let expected = x.exp2() - 1.0;

    assert!(indefinite_sum_with_budget(evaluated, x, &cfg, 15).is_err());
    assert_eq!(count.replace(0), 0);
    let taylor_only = indefinite_sum_with_budget(evaluated, x, &cfg, 16).unwrap();
    assert!((taylor_only - expected).abs() < 1e-8);
    assert_eq!(count.replace(0), 16);

    assert!(indefinite_sum_with_budget(evaluated, x, &cfg, 301).is_err());
    assert_eq!(count.replace(0), 286);
    let both = indefinite_sum_with_budget(evaluated, x, &cfg, 302).unwrap();
    assert_eq!(both.to_bits(), taylor_only.to_bits());
    assert_eq!(count.replace(0), 302);

    let result = indefinite_sum_with_budget(
        |_| {
            count.set(count.get() + 1);
            3.0
        },
        0.25,
        &cfg,
        186,
    )
    .unwrap();
    assert!((result - 0.75).abs() < 1e-10);
    assert_eq!(count.get(), 186);
}

#[test]
fn impossible_shift_uses_available_taylor_and_oversized_interpolation_fails() {
    let cfg = IndefiniteSumConfig {
        shift: usize::MAX,
        ..Default::default()
    };
    let count = Cell::new(0);
    let result = indefinite_sum_with_budget(
        |_| {
            count.set(count.get() + 1);
            4.0
        },
        0.25,
        &cfg,
        14,
    )
    .unwrap();
    assert!((result - 1.0).abs() < 1e-12);
    assert_eq!(count.get(), 14);

    let huge_taylor = IndefiniteSumConfig {
        taylor_terms: usize::MAX,
        ..cfg
    };
    let unused = |_: f64| -> f64 { panic!("neither strategy fits the request") };
    assert!(indefinite_sum_with_budget(unused, 0.25, &huge_taylor, usize::MAX).is_err());
    // An integer request does not use either interpolation parameter.
    assert!(indefinite_sum_with_budget(unused, 0.0, &huge_taylor, 0).is_ok());
}

#[test]
fn product_uses_the_same_budget_and_remains_positive_factor_only() {
    let cfg = IndefiniteSumConfig::default();
    let count = Cell::new(0);
    let factor = |_| {
        count.set(count.get() + 1);
        2.0
    };
    assert!(indefinite_product_with_budget(factor, 3.0, &cfg, 2).is_err());
    assert_eq!(count.get(), 0);
    let product = indefinite_product_with_budget(factor, 3.0, &cfg, 3).unwrap();
    assert!((product - 8.0).abs() < 1e-12);
    assert_eq!(count.get(), 3);
    assert!(indefinite_product_with_budget(|_| -2.0, 1.0, &cfg, 1).is_err());
}

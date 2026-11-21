//! Independent reference tests for bounded series acceleration.

use rssn::kernels::series::sum_terms_to_infinity;
use rssn::kernels::series::sum_to_infinity;
use std::cell::Cell;

#[test]
fn checkpoints_cannot_request_terms_beyond_the_budget() {
    for budget in [0, 1, 3, 7, 15, 16, 17, 23, 31, 32, 33, 49, 63, 64, 65, 127] {
        let calls = Cell::new(0);
        let actual = sum_to_infinity(
            |k| {
                calls.set(calls.get() + 1);
                assert!(calls.get() <= budget, "term requested beyond cap {budget}");
                1.0 / k
            },
            1,
            0.0,
            budget,
        );
        assert_eq!(actual.terms, calls.get());
        assert_eq!(actual.terms, budget);
    }
}

#[test]
fn finite_prefixes_keep_the_partial_sum_when_acceleration_is_underdetermined() {
    let calls = Cell::new(0);
    let values = [1.25, -0.5, 2.0];
    let source = values.into_iter().inspect(|_| calls.set(calls.get() + 1));
    let result = sum_terms_to_infinity(source, 1e-14, 100);
    assert_eq!(result.terms, 3);
    assert_eq!(calls.get(), 3);
    assert_eq!(result.value, 2.75);
    assert_eq!(result.error, f64::INFINITY);

    let no_terms = sum_terms_to_infinity(std::iter::empty(), 1e-14, 100);
    assert!(no_terms.value.is_nan());
    assert_eq!(no_terms.error, f64::INFINITY);
    assert_eq!(no_terms.terms, 0);
}

#[test]
fn the_integer_index_stops_at_the_last_representable_i64() {
    let calls = Cell::new(0_u32);
    let result = sum_to_infinity(
        |index| {
            assert_eq!(index, i64::MAX as f64);
            calls.set(calls.get() + 1);
            f64::from(calls.get())
        },
        i64::MAX - 2,
        1e-14,
        8,
    );
    assert_eq!(calls.get(), 3);
    assert_eq!(result.terms, 3);
    assert_eq!(result.value, 6.0);
    assert_eq!(result.error, f64::INFINITY);
}

#[test]
fn the_nonfinite_term_or_overflow_counts_as_the_last_evaluation() {
    for values in [
        [1.0, f64::NAN, 0.5],
        [1.0, f64::INFINITY, 0.5],
        [f64::MAX, f64::MAX, -f64::MAX],
    ] {
        let calls = Cell::new(0);
        let source = values.into_iter().inspect(|_| calls.set(calls.get() + 1));
        let result = sum_terms_to_infinity(source, 1e-13, 99);
        assert_eq!(calls.get(), 2);
        assert_eq!(result.terms, 2);
        assert!(result.value.is_nan());
        assert_eq!(result.error, f64::INFINITY);
    }
}

#[test]
fn lazy_recurrences_and_indexed_series_approach_analytic_sums() {
    let geometric = std::iter::successors(Some(1.0_f64), |last| Some(0.8 * last));
    let result = sum_terms_to_infinity(geometric, 1e-13, 123);
    assert!((result.value - 5.0).abs() < 5e-12, "{result:?}");
    assert!(result.terms < 123);

    let telescoping = |k: f64| (k * (k + 1.0)).recip();
    let streamed =
        sum_terms_to_infinity((1..=4097).map(|k| telescoping(f64::from(k))), 1e-13, 4097);
    let indexed = sum_to_infinity(telescoping, 1, 1e-13, 4097);
    assert_eq!(streamed, indexed);
    assert!((streamed.value - 1.0).abs() < 5e-11, "{streamed:?}");
}

//! Boundary and numerical reference checks for the inclusive finite sum.

use rssn::kernels::series::sum_range;
use std::cell::Cell;
use std::cell::RefCell;

#[test]
fn the_maximum_index_is_evaluated_once_with_the_final_compensation() {
    let calls = Cell::new(0);
    let sequence = [2.0_f64.powi(52), 0.5, 0.5];
    let value = sum_range(
        |index| {
            let position = calls.get();
            assert!(
                position < sequence.len(),
                "evaluated past the maximum index"
            );
            assert_eq!(index, (i64::MAX - 2 + position as i64) as f64);
            calls.set(position + 1);
            sequence[position]
        },
        i64::MAX - 2,
        i64::MAX,
    );
    assert_eq!(calls.get(), sequence.len());
    assert_eq!(value, 2.0_f64.powi(52) + 1.0);
}

#[test]
fn empty_and_singleton_ranges_at_both_integer_limits() {
    let untouched = sum_range(|_| unreachable!(), i64::MAX, i64::MIN);
    assert_eq!(untouched.to_bits(), 0.0_f64.to_bits());
    for bound in [i64::MIN, i64::MAX] {
        let calls = Cell::new(0);
        assert_eq!(
            sum_range(
                |index| {
                    calls.set(calls.get() + 1);
                    assert_eq!(index, bound as f64);
                    -9.25
                },
                bound,
                bound
            ),
            -9.25
        );
        assert_eq!(calls.get(), 1);
    }
}

#[test]
fn term_order_is_inclusive_across_zero_and_across_the_precision_boundary() {
    for (start, count) in [(-5_i64, 12_i64), ((1_i64 << 53) - 3, 8), (i64::MIN, 2050)] {
        let seen = RefCell::new(Vec::new());
        let total = sum_range(
            |index| {
                seen.borrow_mut().push(index);
                0.25
            },
            start,
            start + (count - 1),
        );
        let expected: Vec<_> = (0..count).map(|offset| (start + offset) as f64).collect();
        assert_eq!(seen.into_inner(), expected);
        assert_eq!(total, 0.25 * count as f64);
    }
}

#[test]
fn a_long_mixed_sign_range_matches_the_exact_integer_polynomial_sum() {
    let (from, to) = (-35_789_i64, 40_319_i64);
    let denominator = 32.0;
    let exact = (from..to + 1).map(|k| k * k - 75 * k + 310).sum::<i64>();
    let actual = sum_range(|k| (k * k - 75.0 * k + 310.0) / denominator, from, to);
    assert_eq!(actual, exact as f64 / denominator);
}

#[test]
fn a_nonfinite_partial_sum_preserves_callback_order_and_ieee_result() {
    assert_eq!(sum_range(|_| f64::INFINITY, 0, 0), f64::INFINITY);
    let calls = Cell::new(0);
    let result = sum_range(
        |index| {
            calls.set(calls.get() + 1);
            if index == 0.0 {
                f64::INFINITY
            } else {
                index
            }
        },
        0,
        2,
    );
    assert!(result.is_nan());
    assert_eq!(calls.get(), 3);
}

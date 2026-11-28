//! Big-integer checks for numerical counts of set partitions.

use num_bigint::BigInt;
use num_bigint::BigUint;
use num_traits::ToPrimitive;
use rssn::kernels::combinatorics::bell;
use rssn::kernels::combinatorics::bell_batch;
use rssn::kernels::combinatorics::stirling_second;

fn choose(
    n: u64,
    k: u64,
) -> BigUint {
    (1..=k).fold(BigUint::from(1_u64), |value, i| value * (n - i + 1) / i)
}

fn stirling_by_surjections(
    n: u32,
    k: u64,
) -> BigUint {
    let mut numerator = BigInt::from(0);
    for j in 0..=k {
        let term = BigInt::from(choose(k, j)) * BigInt::from(j).pow(n);
        if (k - j).is_multiple_of(2) {
            numerator += term;
        } else {
            numerator -= term;
        }
    }
    let denominator = (1..=k).fold(BigInt::from(1), |value, j| value * j);
    (numerator / denominator).to_biguint().unwrap()
}

fn exact_bell(last: usize) -> Vec<BigUint> {
    // Bell's triangle sums adjacent values and obtains a new Bell number at
    // the start of each row. This avoids the production Stirling recurrence.
    let mut previous = vec![BigUint::from(1_u64)];
    let mut numbers = previous.clone();
    for n in 1..=last {
        let mut carry = previous[n - 1].clone();
        let mut row = vec![carry.clone()];
        for value in &previous {
            carry += value;
            row.push(carry.clone());
        }
        numbers.push(row[0].clone());
        previous = row;
    }
    numbers
}

fn assert_count(
    actual: f64,
    expected: &BigUint,
) {
    let expected = expected.to_f64().unwrap_or(f64::INFINITY);
    if expected <= 9_007_199_254_740_992.0 || expected.is_infinite() {
        assert_eq!(actual, expected);
    } else {
        assert!(
            actual.is_finite() && (actual / expected - 1.0).abs() < 1e-13,
            "{actual:e} != {expected:e}"
        );
    }
}

#[test]
fn stirling_matches_exact_surjection_counts() {
    for n in [0, 1, 12, 30, 32, 33, 34, 40, 41, 42, 45, 80, 128] {
        for k in 0..=u64::from(n) {
            assert_count(
                stirling_second(u64::from(n), k),
                &stirling_by_surjections(n, k),
            );
        }
    }
    for (n, k) in [
        (160, 20),
        (180, 30),
        (218, 100),
        (300, 274),
        (383, 256),
        (644, 3),
        (1024, 2),
        (1025, 2),
    ] {
        assert_count(
            stirling_second(u64::from(n), k),
            &stirling_by_surjections(n, k),
        );
    }
}

#[test]
fn stirling_keeps_large_near_diagonal_counts_and_rejects_overflow() {
    for n in [1_000_000, 1_u64 << 55, u64::MAX] {
        assert_count(stirling_second(n, n - 1), &choose(n, 2));
        assert_count(
            stirling_second(n, n - 2),
            &(choose(n, 3) + choose(n, 4) * 3_u32),
        );
        assert_count(
            stirling_second(n, n - 3),
            &(choose(n, 4) + choose(n, 5) * 10_u32 + choose(n, 6) * 15_u32),
        );
        assert_eq!(stirling_second(n, n), 1.0);
        assert_eq!(stirling_second(n, 0), 0.0);
        assert_eq!(stirling_second(n, 1), 1.0);
        assert_eq!(stirling_second(n, 2), f64::INFINITY);
    }
    assert_eq!(stirling_second(u64::MAX - 1, u64::MAX), 0.0);
    assert_eq!(stirling_second(u64::MAX, u64::MAX - 128), f64::INFINITY);
}

#[test]
fn bell_matches_big_integer_triangle_through_first_overflow() {
    let reference = exact_bell(220);
    for (n, value) in reference.iter().enumerate() {
        assert_count(bell(n as u64), value);
    }
    assert!(bell(218).is_finite());
    assert_eq!(bell(219), f64::INFINITY);
    assert_eq!(bell(u64::MAX), f64::INFINITY);
}

#[test]
fn bell_batch_preserves_order_repetitions_and_overflow() {
    assert!(bell_batch(&[]).is_empty());
    let indices = [218, 0, 9, 219, 2, 9, u64::MAX, 1];
    let expected = indices.map(bell);
    assert_eq!(bell_batch(&indices), expected);
}

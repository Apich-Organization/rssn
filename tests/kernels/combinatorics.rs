//! Counting functions and recurrences (ported from `numerical_combinatorics_test.rs`).

use rssn::kernels::combinatorics::*;

#[test]

fn test_factorial() {
    assert_eq!(factorial(0), 1.0);

    assert_eq!(factorial(1), 1.0);

    assert_eq!(factorial(5), 120.0);

    assert!((factorial(10) - 3628800.0).abs() < 1e-9);
}

#[test]

fn test_permutations() {
    assert_eq!(permutations(5, 2), 20.0);

    assert_eq!(permutations(5, 5), 120.0);

    assert_eq!(permutations(5, 6), 0.0);
}

#[test]

fn test_combinations() {
    assert_eq!(combinations(5, 2), 10.0);

    assert_eq!(combinations(5, 5), 1.0);

    assert_eq!(combinations(5, 0), 1.0);

    assert_eq!(combinations(5, 6), 0.0);
}

#[test]

fn test_solve_recurrence_numerical() {
    // Fibonacci: a_n = a_{n-1} + a_{n-2}
    // coeffs = [1.0, 1.0], initial = [0.0, 1.0]
    let coeffs = vec![1.0, 1.0];

    let initial = vec![0.0, 1.0];

    // F(0)=0, F(1)=1, F(2)=1, F(3)=2, F(4)=3, F(5)=5
    assert_eq!(
        solve_recurrence_numerical(&coeffs, &initial, 0).unwrap(),
        0.0
    );

    assert_eq!(
        solve_recurrence_numerical(&coeffs, &initial, 1).unwrap(),
        1.0
    );

    assert_eq!(
        solve_recurrence_numerical(&coeffs, &initial, 5).unwrap(),
        5.0
    );
}

#[test]

fn test_stirling_second() {
    // S(0, 0) = 1
    assert_eq!(stirling_second(0, 0), 1.0);

    // S(n, n) = 1
    assert_eq!(stirling_second(5, 5), 1.0);

    // S(n, 1) = 1
    assert_eq!(stirling_second(5, 1), 1.0);

    // S(3, 2) = 3 ({1,2}, {3}; {1,3}, {2}; {2,3}, {1})
    assert_eq!(stirling_second(3, 2), 3.0);

    // S(4, 2) = 7
    assert_eq!(stirling_second(4, 2), 7.0);
}

#[test]

fn test_bell() {
    // B(0) = 1
    assert_eq!(bell(0), 1.0);

    // B(1) = 1
    assert_eq!(bell(1), 1.0);

    // B(3) = 5 (1+3+1 = 5)
    assert_eq!(bell(3), 5.0);
}

#[test]

fn test_catalan() {
    // C_0 = 1
    assert_eq!(catalan(0), 1.0);

    // C_1 = 1
    assert_eq!(catalan(1), 1.0);

    // C_2 = 2
    assert_eq!(catalan(2), 2.0);

    // C_3 = 5
    assert_eq!(catalan(3), 5.0);
}

#[test]

fn test_rising_factorial() {
    // x^(0) = 1
    assert_eq!(rising_factorial(2.0, 0), 1.0);

    // 2^(3) = 2 * 3 * 4 = 24
    assert_eq!(rising_factorial(2.0, 3), 24.0);
}

#[test]

fn test_falling_factorial() {
    // x_0 = 1
    assert_eq!(falling_factorial(2.0, 0), 1.0);

    // 4_2 = 4 * 3 = 12
    assert_eq!(falling_factorial(4.0, 2), 12.0);

    // 2_3 = 2 * 1 * 0 = 0
    assert_eq!(falling_factorial(2.0, 3), 0.0);
}

#[cfg(test)]
mod proptests {

    use proptest::prelude::*;
    use proptest::test_runner::RngSeed;

    use super::*;

    proptest! {
        #![proptest_config(ProptestConfig { rng_seed: RngSeed::Fixed(0x5EED), failure_persistence: None, ..ProptestConfig::default() })]

        #[test]
        fn prop_factorial_increasing(n in 1..20u64) {
            prop_assert!(factorial(n) >= factorial(n-1));
        }

        #[test]
        fn prop_combinations_symmetry(n in 0..20u64, k in 0..20u64) {
             if k <= n {
                 prop_assert_eq!(combinations(n, k), combinations(n, n - k));
             }
        }

        #[test]
        fn prop_n_choose_k_le_2_pow_n(n in 0..20u64, k in 0..20u64) {
            let n_f64 = n as f64;
            let expected_max = 2.0f64.powf(n_f64);
            if k <= n {
               prop_assert!(combinations(n, k) <= expected_max);
            }
        }

        #[test]
        fn prop_stirling_le_bell(n in 0..10u64, k in 0..10u64) {
             if k <= n {
                 prop_assert!(stirling_second(n, k) <= bell(n));
             }
        }

        #[test]
        fn prop_rising_falling_relationship(x in -10.0..10.0f64, n in 0..5u64) {
            // x^(n) = (-1)^n * (-x)_n ? No
            // x_n = x(x-1)...(x-n+1)
            // x^(n) = x(x+1)...(x+n-1)
            // (-x)_n = (-x)(-x-1)...(-x-n+1) = (-1)^n * x(x+1)...(x+n-1) = (-1)^n * x^(n)

            let term1 = falling_factorial(-x, n);
            let term2 = if n % 2 == 0 { 1.0 } else { -1.0 } * rising_factorial(x, n);
            prop_assert!((term1 - term2).abs() < 1e-9);
        }
    }
}

// ============================================================================
// Added: exact reference values and edge cases
// ============================================================================

#[test]
fn factorial_limits() {
    assert!(factorial(170).is_finite());
    assert!((factorial(170) / 7.257_415_615_307_994e306 - 1.0).abs() < 1e-12);
    assert_eq!(factorial(171), f64::INFINITY);
    assert_eq!(factorial(20), 2_432_902_008_176_640_000.0);
}

#[test]
fn permutations_and_combinations_reference_values() {
    assert_eq!(permutations(10, 3), 720.0);
    assert_eq!(permutations(7, 0), 1.0);
    assert_eq!(combinations(52, 5), 2_598_960.0);
    assert_eq!(combinations(30, 15), 155_117_520.0);
    // Pascal's rule
    for n in 1..25u64 {
        for k in 1..n {
            assert_eq!(combinations(n, k), combinations(n - 1, k - 1) + combinations(n - 1, k), "n={n} k={k}");
        }
    }
}

#[test]
fn stirling_numbers_reference_values() {
    assert_eq!(stirling_second(10, 3), 9330.0);
    assert_eq!(stirling_second(6, 3), 90.0);
    assert_eq!(stirling_second(5, 0), 0.0);
    assert_eq!(stirling_second(3, 5), 0.0);
    // S(n, 2) = 2^(n-1) - 1
    for n in 2..15u64 {
        assert_eq!(stirling_second(n, 2), 2f64.powi(n as i32 - 1) - 1.0);
    }
}

#[test]
fn bell_and_catalan_reference_values() {
    let bell_seq = [1.0, 1.0, 2.0, 5.0, 15.0, 52.0, 203.0, 877.0, 4140.0, 21147.0, 115_975.0];
    for (n, b) in bell_seq.iter().enumerate() {
        assert_eq!(bell(n as u64), *b, "B({n})");
    }
    let catalan_seq = [1.0, 1.0, 2.0, 5.0, 14.0, 42.0, 132.0, 429.0, 1430.0, 4862.0, 16_796.0];
    for (n, c) in catalan_seq.iter().enumerate() {
        assert_eq!(catalan(n as u64), *c, "C({n})");
    }
}

#[test]
fn recurrence_solver_more_cases() {
    // a_n = 2 a_{n-1}, a_0 = 3  =>  3 * 2^n
    assert_eq!(solve_recurrence_numerical(&[2.0], &[3.0], 10), Ok(3072.0));
    // Fibonacci F(30) = 832040
    assert_eq!(solve_recurrence_numerical(&[1.0, 1.0], &[0.0, 1.0], 30), Ok(832_040.0));
    // Lucas numbers L(10) = 123
    assert_eq!(solve_recurrence_numerical(&[1.0, 1.0], &[2.0, 1.0], 10), Ok(123.0));
    // Mismatched initial conditions.
    assert!(solve_recurrence_numerical(&[1.0, 1.0], &[0.0], 5).is_err());
}

#[test]
fn rising_and_falling_factorials_match_gamma_and_binomials() {
    assert_eq!(rising_factorial(1.0, 5), 120.0);
    assert_eq!(falling_factorial(5.0, 5), 120.0);
    assert_eq!(falling_factorial(7.0, 3), permutations(7, 3));
    assert!((rising_factorial(0.5, 3) - 0.5 * 1.5 * 2.5).abs() < 1e-12);
}

proptest::proptest! {
    #![proptest_config(proptest::prelude::ProptestConfig { rng_seed: proptest::test_runner::RngSeed::Fixed(0x5EED), failure_persistence: None, ..proptest::prelude::ProptestConfig::default() })]

    /// Row sums of Pascal's triangle are 2^n.
    #[test]
    fn prop_binomial_row_sums(n in 0..30u64) {
        let s: f64 = (0..=n).map(|k| combinations(n, k)).sum();
        proptest::prop_assert_eq!(s, 2f64.powi(n as i32));
    }

    /// P(n, k) = C(n, k) * k!
    #[test]
    fn prop_permutations_relate_to_combinations(n in 0..20u64, k in 0..20u64) {
        proptest::prop_assume!(k <= n);
        let lhs = permutations(n, k);
        let rhs = combinations(n, k) * factorial(k);
        proptest::prop_assert!((lhs - rhs).abs() <= 1e-9 * lhs.max(1.0));
    }

    /// Bell numbers satisfy B(n+1) = sum_k C(n, k) B(k).
    #[test]
    fn prop_bell_recurrence(n in 0..12u64) {
        let rhs: f64 = (0..=n).map(|k| combinations(n, k) * bell(k)).sum();
        proptest::prop_assert_eq!(bell(n + 1), rhs);
    }

    /// Catalan recurrence C(n+1) = sum C(i) C(n-i).
    #[test]
    fn prop_catalan_convolution(n in 0..12u64) {
        let rhs: f64 = (0..=n).map(|i| catalan(i) * catalan(n - i)).sum();
        proptest::prop_assert_eq!(catalan(n + 1), rhs);
    }

    /// Fibonacci solution via the recurrence solver satisfies F(n+2) = F(n+1) + F(n).
    #[test]
    fn prop_fibonacci_relation(n in 0..60usize) {
        let f = |k| solve_recurrence_numerical(&[1.0, 1.0], &[0.0, 1.0], k).unwrap_or(f64::NAN);
        proptest::prop_assert_eq!(f(n + 2), f(n + 1) + f(n));
    }
}

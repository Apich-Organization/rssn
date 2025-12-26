use rssn::numerical::combinatorics::*;

#[test]

fn test_factorial() {

    assert_eq!(factorial(0), 1.0);

    assert_eq!(factorial(1), 1.0);

    assert_eq!(factorial(5), 120.0);

    assert!(
        (factorial(10) - 3628800.0)
            .abs()
            < 1e-9
    );
}

#[test]

fn test_permutations() {

    assert_eq!(
        permutations(5, 2),
        20.0
    );

    assert_eq!(
        permutations(5, 5),
        120.0
    );

    assert_eq!(
        permutations(5, 6),
        0.0
    );
}

#[test]

fn test_combinations() {

    assert_eq!(
        combinations(5, 2),
        10.0
    );

    assert_eq!(
        combinations(5, 5),
        1.0
    );

    assert_eq!(
        combinations(5, 0),
        1.0
    );

    assert_eq!(
        combinations(5, 6),
        0.0
    );
}

#[test]

fn test_solve_recurrence_numerical() {

    // Fibonacci: a_n = a_{n-1} + a_{n-2}
    // coeffs = [1.0, 1.0], initial = [0.0, 1.0]
    let coeffs = vec![1.0, 1.0];

    let initial = vec![0.0, 1.0];

    // F(0)=0, F(1)=1, F(2)=1, F(3)=2, F(4)=3, F(5)=5
    assert_eq!(
        solve_recurrence_numerical(
            &coeffs,
            &initial,
            0
        )
        .unwrap(),
        0.0
    );

    assert_eq!(
        solve_recurrence_numerical(
            &coeffs,
            &initial,
            1
        )
        .unwrap(),
        1.0
    );

    assert_eq!(
        solve_recurrence_numerical(
            &coeffs,
            &initial,
            5
        )
        .unwrap(),
        5.0
    );
}

#[test]

fn test_stirling_second() {

    // S(0, 0) = 1
    assert_eq!(
        stirling_second(0, 0),
        1.0
    );

    // S(n, n) = 1
    assert_eq!(
        stirling_second(5, 5),
        1.0
    );

    // S(n, 1) = 1
    assert_eq!(
        stirling_second(5, 1),
        1.0
    );

    // S(3, 2) = 3 ({1,2}, {3}; {1,3}, {2}; {2,3}, {1})
    assert_eq!(
        stirling_second(3, 2),
        3.0
    );

    // S(4, 2) = 7
    assert_eq!(
        stirling_second(4, 2),
        7.0
    );
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
    assert_eq!(
        rising_factorial(2.0, 0),
        1.0
    );

    // 2^(3) = 2 * 3 * 4 = 24
    assert_eq!(
        rising_factorial(2.0, 3),
        24.0
    );
}

#[test]

fn test_falling_factorial() {

    // x_0 = 1
    assert_eq!(
        falling_factorial(2.0, 0),
        1.0
    );

    // 4_2 = 4 * 3 = 12
    assert_eq!(
        falling_factorial(4.0, 2),
        12.0
    );

    // 2_3 = 2 * 1 * 0 = 0
    assert_eq!(
        falling_factorial(2.0, 3),
        0.0
    );
}

#[cfg(test)]

mod proptests {

    use proptest::prelude::*;

    use super::*;

    proptest! {
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

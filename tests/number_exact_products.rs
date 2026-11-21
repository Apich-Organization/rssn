//! Exact combinatorial products at word and big-integer boundaries.

use num_bigint::BigInt;
use num_traits::One;
use num_traits::Zero;
use rssn::kernels::number_exact::binomial;
use rssn::kernels::number_exact::falling_factorial;
use rssn::kernels::number_exact::permutations;

#[test]
fn exact_products_keep_empty_and_impossible_selections_distinct() {
    assert_eq!(permutations(-1), BigInt::zero());
    for n in [0, 1, 9, i64::MAX] {
        assert_eq!(falling_factorial(n, 0), BigInt::one());
        assert_eq!(binomial(n, 0), BigInt::one());
    }
    assert_eq!(falling_factorial(-1, 0), BigInt::zero());
    assert_eq!(falling_factorial(2, 3), BigInt::zero());
    assert_eq!(binomial(1, -1), BigInt::zero());
    assert_eq!(falling_factorial(6, 3), BigInt::from(120));
    assert_eq!(permutations(10), BigInt::from(3_628_800));
}

#[test]
fn balanced_products_match_a_sequential_reference() {
    for n in [0, 1, 31, 32, 33, 65, 129, 256, 1024] {
        let expected: BigInt = (2..=n).map(BigInt::from).product();
        assert_eq!(permutations(n), expected);
        assert_eq!(falling_factorial(n, n), expected);
    }
    for (n, k) in [
        (200, 64),
        (200, 65),
        (256, 128),
        (1000, 500),
        (i64::MAX, 87),
    ] {
        let mut expected = BigInt::one();
        for i in 0..k {
            expected = expected * BigInt::from(n - i) / BigInt::from(i + 1);
        }
        assert_eq!(binomial(n, k), expected);
        assert_eq!(binomial(n, n - k), expected);
        assert_eq!(falling_factorial(n, k), expected * permutations(k));
    }
}

//! Finite fields GF(p) and GF(2^8) (ported from `numerical_finite_field_test.rs`).
//! The `matrix::Field` impl for `PrimeFieldElement` no longer exists.

use proptest::prelude::*;
use proptest::test_runner::RngSeed;
use rssn::kernels::finite_field::{
    PrimeFieldElement, gf256_add, gf256_div, gf256_inv, gf256_mul, gf256_pow,
};

fn cfg() -> ProptestConfig {
    ProptestConfig {
        rng_seed: RngSeed::Fixed(0x5EED),
        failure_persistence: None,
        ..ProptestConfig::default()
    }
}

#[test]
fn prime_field_arithmetic() {
    let p = 11;
    let a = PrimeFieldElement::new(7, p);
    let b = PrimeFieldElement::new(5, p);
    assert_eq!((a + b).value, 1); // 12 mod 11
    assert_eq!((a * b).value, 2); // 35 mod 11
    assert_eq!((a - b).value, 2);
    assert_eq!((b - a).value, 9); // -2 mod 11
    assert_eq!((-a).value, 4);
    assert_eq!((a + b).modulus, p);
}

#[test]
fn prime_field_construction_reduces() {
    assert_eq!(PrimeFieldElement::new(25, 7).value, 4);
    assert_eq!(PrimeFieldElement::new(7, 7).value, 0);
}

#[test]
fn prime_field_inverse() {
    let a = PrimeFieldElement::new(7, 11);
    let inv = a.inverse().unwrap_or_else(|| panic!("7 is invertible mod 11"));
    assert_eq!(inv.value, 8); // 7 * 8 = 56 = 1 mod 11
    assert_eq!((a * inv).value, 1);
    assert!(PrimeFieldElement::new(0, 11).inverse().is_none());
    // Composite modulus: 4 has no inverse mod 8, 3 does.
    assert!(PrimeFieldElement::new(4, 8).inverse().is_none());
    assert_eq!(PrimeFieldElement::new(3, 8).inverse().map(|e| e.value), Some(3));
}

#[test]
fn prime_field_division_and_assign_ops() {
    let a = PrimeFieldElement::new(3, 13);
    let b = PrimeFieldElement::new(5, 13);
    let q = a / b;
    assert_eq!((q * b).value, a.value);
    let mut x = a;
    x += b;
    assert_eq!(x.value, 8);
    x -= b;
    assert_eq!(x, a);
    x *= b;
    assert_eq!(x.value, 2); // 15 mod 13
    x /= b;
    assert_eq!(x, a);
}

#[test]
fn prime_field_pow() {
    let a = PrimeFieldElement::new(2, 11);
    assert_eq!(a.pow(5).value, 10); // 32 mod 11
    assert_eq!(a.pow(10).value, 1); // Fermat
    assert_eq!(a.pow(0).value, 1);
    // Large exponent stays exact via u128 arithmetic.
    let p = 1_000_000_007u64;
    assert_eq!(PrimeFieldElement::new(123_456, p).pow(p - 1).value, 1);
}

#[test]
fn gf256_addition_is_xor() {
    assert_eq!(gf256_add(0x57, 0x83), 0x57 ^ 0x83);
    assert_eq!(gf256_add(0xAB, 0xAB), 0);
}

#[test]
fn gf256_multiplication_and_inverse() {
    let a = 0x57;
    let inv = gf256_inv(a).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(gf256_mul(a, inv), 1);
    assert!(gf256_inv(0).is_err());
    // Polynomial 0x11d: x^8 = x^4 + x^3 + x^2 + 1, so 2 * 0x80 = 0x1d.
    assert_eq!(gf256_mul(2, 0x80), 0x1d);
    assert_eq!(gf256_mul(0, 0x53), 0);
    assert_eq!(gf256_mul(1, 0x53), 0x53);
}

#[test]
fn gf256_division() {
    assert_eq!(gf256_div(gf256_mul(0x1c, 0x53), 0x53), Ok(0x1c));
    assert_eq!(gf256_div(0, 5), Ok(0));
    assert!(gf256_div(5, 0).is_err());
}

#[test]
fn gf256_power() {
    assert_eq!(gf256_pow(2, 0), 1);
    assert_eq!(gf256_pow(2, 1), 2);
    assert_eq!(gf256_pow(2, 2), 4);
    assert_eq!(gf256_pow(2, 3), 8);
    assert_eq!(gf256_pow(2, 8), 0x1d);
    // 2 generates the multiplicative group of order 255.
    assert_eq!(gf256_pow(2, 255), 1);
    assert_eq!(gf256_pow(0, 5), 0);
}

#[test]
fn gf256_two_is_a_primitive_element() {
    let mut seen = [false; 256];
    for e in 0..255u64 {
        let v = gf256_pow(2, e) as usize;
        assert!(!seen[v], "2^{e} repeated value {v}");
        seen[v] = true;
    }
    assert!(!seen[0]);
}

proptest! {
    #![proptest_config(cfg())]

    #[test]
    fn prop_prime_field_mul_inverse(v in 1..100u64, p_idx in 0..10usize) {
        let primes = [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29];
        let p = primes[p_idx];
        let val = v % p;
        prop_assume!(val != 0);
        let a = PrimeFieldElement::new(val, p);
        let inv = a.inverse().ok_or_else(|| TestCaseError::fail("no inverse in a prime field"))?;
        prop_assert_eq!((a * inv).value, 1);
    }

    #[test]
    fn prop_prime_field_distributes(a in 0..1000u64, b in 0..1000u64, c in 0..1000u64) {
        let p = 10_007;
        let (a, b, c) = (PrimeFieldElement::new(a, p), PrimeFieldElement::new(b, p), PrimeFieldElement::new(c, p));
        prop_assert_eq!(a * (b + c), a * b + a * c);
        prop_assert_eq!((a + b) - b, a);
    }

    #[test]
    fn prop_gf256_mul_inverse(a in 1..=255u8) {
        let inv = gf256_inv(a).map_err(TestCaseError::fail)?;
        prop_assert_eq!(gf256_mul(a, inv), 1);
    }

    #[test]
    fn prop_gf256_field_axioms(a in any::<u8>(), b in any::<u8>(), c in any::<u8>()) {
        prop_assert_eq!(gf256_mul(a, b), gf256_mul(b, a));
        prop_assert_eq!(gf256_mul(a, gf256_mul(b, c)), gf256_mul(gf256_mul(a, b), c));
        prop_assert_eq!(gf256_mul(a, gf256_add(b, c)), gf256_add(gf256_mul(a, b), gf256_mul(a, c)));
    }
}

mod modulus_agnostic {
    use num_traits::{One, Zero};
    use rssn::kernels::finite_field::PrimeFieldElement;

    #[test]
    fn zero_and_one_adopt_the_other_operands_modulus() {
        let x = PrimeFieldElement::new(5, 7);
        for (a, b) in [(x + PrimeFieldElement::zero(), x), (PrimeFieldElement::zero() + x, x)] {
            assert_eq!(a, b);
        }
        assert_eq!(x * PrimeFieldElement::one(), x);
        assert_eq!(PrimeFieldElement::one() * x, x);
        assert_eq!(x - PrimeFieldElement::zero(), x);
        assert_eq!((PrimeFieldElement::zero() - x).value, 2);
        assert_eq!((PrimeFieldElement::zero() - x).modulus, 7);
        assert_eq!((x / PrimeFieldElement::one()), x);
        assert_eq!((PrimeFieldElement::one() / x).value, 3); // 5 * 3 = 15 = 1 mod 7
    }

    #[test]
    fn constants_are_flagged_and_zero_one_predicates_work() {
        assert!(PrimeFieldElement::zero().is_agnostic());
        assert!(PrimeFieldElement::one().is_agnostic());
        assert!(PrimeFieldElement::zero().is_zero());
        assert!(PrimeFieldElement::one().is_one());
        assert!(!PrimeFieldElement::new(3, 7).is_agnostic());
        assert!(PrimeFieldElement::new(7, 7).is_zero());
    }

    #[test]
    fn sums_and_products_fold_from_the_identity() {
        let xs: Vec<_> = (1..=6).map(|v| PrimeFieldElement::new(v, 7)).collect();
        let sum = xs.iter().fold(PrimeFieldElement::zero(), |a, &b| a + b);
        assert_eq!((sum.value, sum.modulus), (0, 7)); // 21 mod 7
        let prod = xs.iter().fold(PrimeFieldElement::one(), |a, &b| a * b);
        assert_eq!((prod.value, prod.modulus), (6, 7)); // Wilson: 6! = -1 mod 7
    }

    #[test]
    fn agnostic_values_reduce_when_bound() {
        let two = PrimeFieldElement::one() + PrimeFieldElement::one();
        assert_eq!(two.value, 2);
        assert!(two.is_agnostic());
        let y = PrimeFieldElement::new(6, 7);
        assert_eq!((y + two).value, 1);
        assert_eq!((two * y).value, 5);
        // an agnostic value larger than the modulus is reduced, not wrapped
        let big = PrimeFieldElement::new(100, 0);
        assert_eq!((PrimeFieldElement::new(0, 7) + big).value, 2);
    }

    #[test]
    fn agnostic_arithmetic_never_panics() {
        let z = PrimeFieldElement::zero();
        let o = PrimeFieldElement::one();
        assert_eq!((z - o).value, 0); // saturates
        assert_eq!(o.inverse(), Some(o));
        assert_eq!(z.inverse(), None);
        assert_eq!((o / z).value, 0);
        assert_eq!(o.pow(1000).value, 1);
        assert_eq!((-z).value, 0);
    }

    #[test]
    fn large_modulus_does_not_overflow() {
        let p = 18_446_744_073_709_551_557_u64; // largest prime below 2^64
        let a = PrimeFieldElement::new(p - 1, p);
        assert_eq!((a + a).value, p - 2);
        assert_eq!((PrimeFieldElement::new(1, p) - a).value, 2);
        assert_eq!((a * a).value, 1);
    }
}

//! Exact reference cases for integer powers of truncated holonomic series.

#![allow(clippy::unwrap_used)]

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::One;
use num_traits::Zero;
use rssn::kernels::holonomic::Atom;
use rssn::kernels::holonomic::Expr;
use rssn::kernels::holonomic::RatFn;

type Q = BigRational;

fn q(
    n: i64,
    d: i64,
) -> Q {
    Q::new(BigInt::from(n), BigInt::from(d))
}

fn integer_power(
    atom: Atom,
    power: u32,
) -> Expr {
    Expr::Pow(Box::new(Expr::Atom(atom)), power)
}

#[test]
fn maximum_integer_exponent_matches_signed_binomial_coefficients() {
    let exponent = u32::MAX;
    let expr = integer_power(
        Atom::Rational(RatFn::poly(vec![q(-1, 1), q(2, 1)])),
        exponent,
    );
    let mut binomial = BigInt::one();
    let expected: Vec<_> = (0..24)
        .map(|k| {
            // (-1 + 2x)^e has coefficient (-1)^(e-k) 2^k C(e,k).
            let sign = if k % 2 == 0 { -1 } else { 1 };
            let result = Q::from_integer((&binomial << k) * sign);
            binomial = (&binomial * (exponent - k as u32)) / (k + 1);
            result
        })
        .collect();
    assert_eq!(expr.series(expected.len()).unwrap(), expected);
}

#[test]
fn factors_positive_order_and_preserves_gaps_and_the_leading_term() {
    let exponent = 17_u32;
    let rat = RatFn::new(vec![q(0, 1), q(0, 1), q(-2, 5)], vec![q(1, 1), q(-1, 3)]).unwrap();
    let expr = integer_power(Atom::Rational(rat), exponent);
    let order = 2 * exponent as usize;
    let scale = Q::new(
        BigInt::from(-2).pow(exponent),
        BigInt::from(5).pow(exponent),
    );
    assert_eq!(expr.series(order).unwrap(), vec![Q::zero(); order]);
    assert_eq!(expr.series(order + 1).unwrap()[order], scale);

    // [x^k](1 - x/3)^(-e) = C(e + k - 1, k) / 3^k.
    let mut expected = vec![Q::zero(); order];
    let mut term = scale;
    for k in 0..24 {
        expected.push(term.clone());
        term *= q(i64::from(exponent) + k, 3 * (k + 1));
    }
    assert_eq!(expr.series(expected.len()).unwrap(), expected);

    let monomial = Atom::Rational(RatFn::poly(vec![q(0, 1), q(0, 1), q(-2, 5)]));
    let huge = integer_power(monomial, u32::MAX);
    assert_eq!(huge.series(64).unwrap(), vec![Q::zero(); 64]);

    let gapped = Atom::Rational(RatFn::poly(vec![q(1, 1), q(0, 1), q(1, 1)]));
    let expected: Vec<_> = [1, 0, 3, 0, 3, 0, 1, 0, 0, 0, 0].map(|c| q(c, 1)).to_vec();
    assert_eq!(
        integer_power(gapped, 3).series(expected.len()).unwrap(),
        expected
    );
}

#[test]
fn dense_exponential_power_agrees_with_its_known_generating_function() {
    let exponent = 256_u32;
    let expr = integer_power(Atom::Exp(RatFn::poly(vec![q(0, 1), q(1, 1)])), exponent);
    let mut term = Q::one();
    let expected: Vec<_> = (1..=24)
        .map(|k| {
            let value = term.clone();
            term *= q(i64::from(exponent), k);
            value
        })
        .collect();
    assert_eq!(expr.series(expected.len()).unwrap(), expected);
}

#[test]
fn low_powers_zero_and_invalid_atoms_keep_their_series_contract() {
    let polynomial = [q(1, 2), q(-2, 3), q(1, 5)];
    let atom = Atom::Rational(RatFn::poly(polynomial.to_vec()));
    assert_eq!(
        integer_power(atom.clone(), 1).series(3).unwrap(),
        polynomial
    );
    assert_eq!(
        integer_power(atom, 2).series(7).unwrap(),
        [
            q(1, 4),
            q(-2, 3),
            q(29, 45),
            q(-4, 15),
            q(1, 25),
            q(0, 1),
            q(0, 1)
        ]
    );

    let zero = Atom::Rational(RatFn::zero());
    assert_eq!(
        integer_power(zero.clone(), 0).series(3).unwrap(),
        [q(1, 1), q(0, 1), q(0, 1)]
    );
    assert_eq!(
        integer_power(zero.clone(), u32::MAX).series(3).unwrap(),
        vec![Q::zero(); 3]
    );
    assert!(integer_power(zero, u32::MAX).series(0).unwrap().is_empty());

    let invalid = Atom::Sin(vec![q(1, 1)]);
    assert!(integer_power(invalid.clone(), 0).series(4).is_none());
    assert!(integer_power(invalid, 2).series(0).is_none());
}

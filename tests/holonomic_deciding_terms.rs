//! Exact coefficient and search bounds for holonomic zero decisions.

#![allow(clippy::unwrap_used)]

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::One;
use num_traits::Zero;
use rssn::kernels::holonomic::Limits;
use rssn::kernels::holonomic::Operator;
use rssn::kernels::holonomic::deciding_terms;

type Q = BigRational;

fn q(
    n: i64,
    d: i64,
) -> Q {
    Q::new(BigInt::from(n), BigInt::from(d))
}

const fn limits(max_terms: usize) -> Limits {
    Limits {
        max_order: 128,
        max_degree: 128,
        max_terms,
    }
}

fn falling_operator(
    coefficients: &[Q],
    first_order: usize,
    shift: i64,
) -> Operator {
    let mut op = vec![vec![]; first_order + coefficients.len()];
    for (i, c) in coefficients.iter().enumerate() {
        let order = first_order + i;
        let degree = usize::try_from(i128::try_from(order).unwrap() - i128::from(shift)).unwrap();
        op[order] = vec![Q::zero(); degree];
        op[order].push(c.clone());
    }
    op
}

#[test]
fn no_nonnegative_roots_can_certify_zero_with_a_zero_term_budget() {
    let huge = Q::from_integer(BigInt::one() << 192_usize);
    // x f' + huge f has recurrence multiplier n + huge.
    let linear = falling_operator(&[huge.clone(), q(1, 1)], 0, 0);
    assert_eq!(deciding_terms(&linear, &limits(0)), Some(0));

    // -huge - n(n-1)...(n-63) is strictly negative at nonnegative integers.
    let mut same_sign = vec![Q::zero(); 65];
    same_sign[0] = -huge;
    same_sign[64] = q(-1, 1);
    assert_eq!(
        deciding_terms(&falling_operator(&same_sign, 0, 0), &limits(0)),
        Some(0)
    );

    assert_eq!(deciding_terms(&[vec![Q::one()]], &limits(0)), Some(0));
    assert_eq!(deciding_terms(&[], &limits(usize::MAX)), None);
    assert_eq!(deciding_terms(&[vec![], vec![Q::zero()]], &limits(8)), None);
}

#[test]
fn common_falling_factor_counts_origin_and_high_derivative_orders() {
    // x f'=0 leaves just the constant coefficient free (the m=0 root).
    let origin = falling_operator(&[Q::one()], 1, 0);
    assert_eq!(deciding_terms(&origin, &limits(0)), None);
    assert_eq!(deciding_terms(&origin, &limits(1)), Some(1));

    // f^(64)=0 has arbitrary coefficients below degree 64.
    let derivative = falling_operator(&[Q::one()], 64, 64);
    assert_eq!(deciding_terms(&derivative, &limits(64)), Some(64));
    assert_eq!(deciding_terms(&derivative, &limits(63)), None);
}

#[test]
fn higher_orders_scan_the_inclusive_root_bound_in_integer_coordinates() {
    // The residual P(z) = (z-4)(z+1/2) has its only positive integer root
    // at the inclusive Cauchy bound 4. The leading order supplies two more
    // free coefficients, and the negative recurrence shift changes neither.
    let coefficients = [q(-2, 1), q(-5, 2), q(1, 1)];
    let mut op = falling_operator(&coefficients, 2, -3);
    // This term has a lower shift and cannot contribute to the leading term.
    op[4].push(q(1000, 7));
    assert_eq!(deciding_terms(&op, &limits(7)), Some(7));
    assert_eq!(deciding_terms(&op, &limits(6)), None);

    // Independently convert a factored polynomial to falling-factorial
    // coordinates by forward differences. Its largest integer root is 9.
    let mut values: Vec<BigInt> = (0..=4)
        .map(|x| BigInt::from((x - 1) * (x - 4) * (x - 9) * (x + 2)))
        .collect();
    let mut factorial = BigInt::one();
    let mut coefficients = Vec::new();
    for i in 1..=5 {
        coefficients.push(Q::new(values[0].clone(), factorial.clone()));
        values = values.windows(2).map(|w| &w[1] - &w[0]).collect();
        factorial *= i;
    }
    let op = falling_operator(&coefficients, 0, 0);
    assert_eq!(deciding_terms(&op, &limits(512)), Some(10));
}

#[test]
fn linear_integer_root_and_term_count_are_exact_at_usize_endpoints() {
    let max = BigInt::from(usize::MAX);
    let root = Q::from_integer(&max - 1_u8);
    let op = falling_operator(&[-root, Q::one()], 0, -1);
    assert_eq!(deciding_terms(&op, &limits(usize::MAX)), Some(usize::MAX));
    assert_eq!(deciding_terms(&op, &limits(usize::MAX - 1)), None);

    let root_too_high = Q::from_integer(max.clone());
    let op = falling_operator(&[-root_too_high, Q::one()], 0, -1);
    assert_eq!(deciding_terms(&op, &limits(usize::MAX)), None);

    let positive_noninteger = Q::new(2 * max + 1, BigInt::from(2));
    let op = falling_operator(&[-positive_noninteger, Q::one()], 0, 0);
    assert_eq!(deciding_terms(&op, &limits(0)), Some(0));
}

//! Exact enclosure regressions for arbitrary-precision real ball division.

use num_bigint::BigInt;
use num_rational::BigRational;
use rssn::kernels::mp::Ball;
use rssn::kernels::mp::Real;

fn rational(
    numerator: i64,
    denominator: i64,
) -> BigRational {
    BigRational::new(BigInt::from(numerator), BigInt::from(denominator))
}

fn exact(value: &Real) -> BigRational {
    let exponent = value.exponent();
    let shift = exponent.unsigned_abs();
    if exponent >= 0 {
        BigRational::from_integer(value.mantissa() << shift)
    } else {
        BigRational::new(value.mantissa().clone(), BigInt::from(1_u32) << shift)
    }
}

fn assert_encloses(
    ball: &Ball,
    minimum: BigRational,
    maximum: BigRational,
) {
    let midpoint = exact(ball.midpoint());
    let radius = exact(ball.radius());
    let lower = &midpoint - &radius;
    let upper = &midpoint + &radius;
    assert!(
        lower <= minimum,
        "lower endpoint {lower} excludes {minimum}"
    );
    assert!(
        upper >= maximum,
        "upper endpoint {upper} excludes {maximum}"
    );
}

#[test]
fn point_division_encloses_the_exact_rational_quotient() {
    let one = Ball::point(Real::one(1));
    let three = Ball::point(Real::from_i64(3, 16).unwrap());
    for sign in [1, -1] {
        let numerator = Ball::point(Real::from_i64(sign, 16).unwrap());
        let quotient = numerator.div(&three, 4).unwrap();
        assert_eq!(exact(quotient.midpoint()), rational(11 * sign, 32));
        assert_encloses(&quotient, rational(sign, 3), rational(sign, 3));
    }

    let exact_quotient = three.div(&one, 2).unwrap();
    assert!(exact_quotient.radius().is_zero());
    let wide = Ball::point(Real::from_i64(1025, 11).unwrap());
    assert_encloses(
        &wide.div(&one, 1).unwrap(),
        rational(1025, 1),
        rational(1025, 1),
    );
}

#[test]
fn quotient_encloses_both_input_radii_with_either_divisor_sign() {
    let radius = Real::one(16).scale_pow2(-10).unwrap();
    let numerator = Ball::new(Real::one(16), radius.clone()).unwrap();
    let positive = Ball::new(Real::from_i64(3, 16).unwrap(), radius).unwrap();
    assert_encloses(
        &numerator.div(&positive, 4).unwrap(),
        rational(1023, 3073),
        rational(1025, 3071),
    );
    assert_encloses(
        &numerator.div(&positive.neg(), 4).unwrap(),
        rational(-1025, 3071),
        rational(-1023, 3073),
    );
}

#[test]
fn divisor_zero_is_rejected_only_when_in_the_closed_interval() {
    let positive = Real::one(16);
    let numerator = Ball::point(Real::zero(16));
    for midpoint in [Real::zero(16), positive.clone(), positive.neg()] {
        let overlaps_zero = Ball::new(midpoint, positive.clone()).unwrap();
        assert_eq!(
            numerator.div(&overlaps_zero, 4).unwrap_err(),
            "division by an interval containing zero"
        );
    }
    let radius = Real::from_i64(1023, 16).unwrap().scale_pow2(-10).unwrap();
    let excludes_zero = Ball::new(positive.clone(), radius).unwrap();
    let numerator = Ball::point(positive);
    assert_encloses(
        &numerator.div(&excludes_zero, 4).unwrap(),
        rational(1024, 2047),
        rational(1024, 1),
    );
}

#[test]
fn midpoint_rounding_is_enclosed_at_the_exponent_floor() {
    let exponent = rssn::kernels::mp::MAX_EXP.saturating_neg();
    let narrow = Ball::point(Real::from_i64(17, 5).unwrap().scale_pow2(exponent).unwrap());
    let two = Ball::point(Real::from_i64(2, 1).unwrap());
    let quotient = narrow.div(&two, 4).unwrap();
    let expected = exact(narrow.midpoint()) * rational(1, 2);
    assert_encloses(&quotient, expected.clone(), expected);
}

#[test]
fn quotient_bounds_do_not_square_the_divisor_scale() {
    for scale in [-170_000, 170_000] {
        let mid = Real::one(16).scale_pow2(scale).unwrap();
        let rad = Real::one(16).scale_pow2(scale - 1).unwrap();
        let input = Ball::new(mid, rad).unwrap();
        let quotient = input.div(&input, 16).unwrap();
        assert_encloses(&quotient, rational(1, 3), rational(3, 1));
    }
}

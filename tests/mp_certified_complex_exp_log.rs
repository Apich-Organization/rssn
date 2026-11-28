//! Independent enclosure checks for MP complex modulus, exp and principal ln.

use num_bigint::BigInt;
use num_rational::BigRational;
use rssn::kernels::mp::Ball;
use rssn::kernels::mp::MAX_BITS;
use rssn::kernels::mp::MAX_EXP;
use rssn::kernels::mp::MpComplex;
use rssn::kernels::mp::Real;
use rssn::kernels::mp::Round;

fn int(n: i64) -> BigRational {
    BigRational::from_integer(n.into())
}

fn dyadic(
    mantissa: i64,
    exponent: i64,
    precision: u32,
) -> Real {
    Real::rounded(mantissa.into(), exponent, precision, Round::NearestEven).unwrap()
}

fn exact(x: &Real) -> BigRational {
    let shift = x.exponent().unsigned_abs();
    if x.exponent() < 0 {
        BigRational::new(x.mantissa().clone(), BigInt::from(1) << shift)
    } else {
        BigRational::from_integer(x.mantissa() << shift)
    }
}

fn endpoints(ball: &Ball) -> (BigRational, BigRational) {
    let mid = exact(&ball.mid);
    let rad = exact(&ball.rad);
    (&mid - &rad, mid + rad)
}

fn encloses(
    ball: &Ball,
    lo: &BigRational,
    hi: &BigRational,
) {
    let (lower, upper) = endpoints(ball);
    assert!(
        &lower <= lo && &upper >= hi,
        "{ball:?}: lower covered {}, upper covered {}",
        &lower <= lo,
        &upper >= hi
    );
}

fn ln2_reference() -> (BigRational, BigRational) {
    let mut total = int(0);
    for k in 1..=150 {
        total += BigRational::new(1.into(), BigInt::from(k) << k as u32);
    }
    let upper = &total + BigRational::new(1.into(), BigInt::from(151) << 150_u32);
    (total, upper)
}

fn quarter_pi_reference() -> (BigRational, BigRational) {
    let atan = |n| {
        let x = int(n).recip();
        let square = &x * &x;
        let mut power = x;
        let mut sum = int(0);
        for k in 0..100 {
            let term = &power / int(2 * k + 1);
            if k % 2 == 0 {
                sum += &term;
            } else {
                sum -= &term;
            }
            power *= &square;
        }
        let hi = &sum + power / int(201);
        (sum, hi)
    };
    let (a, b) = atan(2);
    let (c, d) = atan(3);
    (a + c, b + d)
}

#[test]
fn rectangle_modulus_uses_nearest_and_farthest_corners() {
    let point = MpComplex::point(dyadic(3, 0, 2), dyadic(1, 2, 1));
    let abs = point.abs(64).unwrap();
    assert_eq!(exact(&abs.mid), int(5));
    assert!(abs.rad.is_zero());
    let re = Ball::new(dyadic(1, 1, 1), Real::one(1)).unwrap();
    let im = Ball::new(dyadic(1, 2, 1), Real::one(1)).unwrap();
    let bounds = MpComplex::new(re, im).abs(64).unwrap();
    let (low, high) = endpoints(&bounds);
    assert!(low > int(3) && low.pow(2) <= int(10));
    assert!(high < int(6) && high.pow(2) >= int(34));
}

#[test]
fn complex_exponential_contains_interior_trigonometric_extrema() {
    let zero = Real::zero(1);
    let span = Ball::new(zero.clone(), dyadic(3, 0, 2)).unwrap();
    let result = MpComplex::new(Ball::point(zero.clone()), span)
        .exp(32)
        .unwrap();
    // The interval [-3, 3] contains ±π/2, so its sine attains both extrema.
    encloses(&result.im, &int(-1), &int(1));
    let (_, real_upper) = endpoints(&result.re);
    assert!(real_upper >= int(1));

    let unit = MpComplex::point(zero.clone(), zero.clone()).exp(8).unwrap();
    assert_eq!(exact(&unit.re.mid), int(1));
    assert!(unit.re.rad.is_zero() && unit.im.mid.is_zero() && unit.im.rad.is_zero());

    let positive = MpComplex::point(Real::one(1), zero).exp(40).unwrap();
    let mut factorial = BigInt::from(1);
    let mut lower = int(1);
    for k in 1..=36 {
        factorial *= k;
        lower += BigRational::new(1.into(), factorial.clone());
    }
    let upper = &lower + BigRational::new(2.into(), factorial * 37);
    encloses(&positive.re, &lower, &upper);
}

#[test]
fn principal_log_encloses_separated_rectangles_and_rejects_the_cut() {
    let one = Real::one(1);
    let zero = Real::zero(1);
    let rectangle = MpComplex::new(
        Ball::new(zero.clone(), one.clone()).unwrap(),
        Ball::point(one.clone()),
    );
    let result = rectangle.ln(64).unwrap();
    let (_, ln2_hi) = ln2_reference();
    let (pi_lo, pi_hi) = quarter_pi_reference();
    encloses(&result.re, &int(0), &(ln2_hi / int(2)));
    encloses(&result.im, &pi_lo, &(int(3) * pi_hi));

    let cut = MpComplex::new(
        Ball::point(one.neg()),
        Ball::new(zero.clone(), dyadic(1, -3, 1)).unwrap(),
    );
    assert!(cut.ln(32).is_err());
    assert!(MpComplex::point(zero.clone(), zero.clone()).ln(32).is_err());
    let log_one = MpComplex::point(one, zero).ln(32).unwrap();
    assert!(log_one.re.mid.is_zero() && log_one.re.rad.is_zero());
    assert!(log_one.im.mid.is_zero() && log_one.im.rad.is_zero());
}

#[test]
fn logarithm_keeps_norm_scaling_outside_the_real_exponent_range() {
    let exact_tiny = dyadic(1, -MAX_EXP, 1);
    let tiny_log = MpComplex::point(exact_tiny, Real::zero(1)).ln(64).unwrap();
    let (ln2_low, ln2_high) = ln2_reference();
    encloses(
        &tiny_log.re,
        &(-int(MAX_EXP) * &ln2_high),
        &(-int(MAX_EXP) * &ln2_low),
    );
    assert!(tiny_log.im.mid.is_zero() && tiny_log.im.rad.is_zero());

    let biggest = Real::rounded(
        (BigInt::from(1) << MAX_BITS) - 1_u32,
        MAX_EXP,
        MAX_BITS,
        Round::NearestEven,
    )
    .unwrap();
    let huge = MpComplex::point(biggest.clone(), biggest);
    assert!(huge.abs(64).is_err());
    let logarithm = huge.ln(64).unwrap();
    // The exact log modulus is (MAX_EXP + MAX_BITS + 1/2) ln(2)
    // plus ln(1 - 2^-MAX_BITS), bounded below by -2 * 2^-MAX_BITS.
    let scale = int(MAX_EXP + i64::from(MAX_BITS)) + int(2).recip();
    let correction = BigRational::new(2.into(), BigInt::from(1) << MAX_BITS);
    encloses(
        &logarithm.re,
        &(&scale * ln2_low - correction),
        &(scale * ln2_high),
    );
    let (lo, hi) = quarter_pi_reference();
    encloses(&logarithm.im, &lo, &hi);
}

#[test]
fn session_exposes_declared_complex_rectangle_contracts() {
    use rssn::api::Session;
    use rssn::kernel::CallCtx;
    use rssn::kernel::Value;
    use rssn::kernel::manifest::SignatureSource;
    use rssn::kernel::signature::DimRuleInfo;
    use rssn::kernel::value::RecordValue;

    let record = |fields: &[(&str, Value)]| {
        Value::Record(RecordValue {
            fields: fields
                .iter()
                .map(|(name, value)| ((*name).into(), value.clone()))
                .collect(),
        })
    };
    let dyadic = |mantissa: i64, exponent: i64, precision: u32| {
        record(&[
            ("mantissa", Value::int(mantissa)),
            ("exponent", Value::int(exponent)),
            ("precision", Value::int(precision)),
        ])
    };
    let ball = |midpoint, radius| record(&[("midpoint", midpoint), ("radius", radius)]);
    let rectangle = |real, imag| record(&[("real", real), ("imag", imag)]);
    let session = Session::new();
    let ctx = CallCtx::new();
    let call = |op: &str, input: Value, precision: u32| {
        session.call_value(op, &[input, Value::int(precision)], &ctx)
    };
    let zero = dyadic(0, 0, 1);
    let one = dyadic(1, 0, 1);
    let zero_ball = ball(zero.clone(), zero.clone());
    let one_ball = ball(one.clone(), zero.clone());
    assert_eq!(
        call(
            "num_mp_complex_abs",
            rectangle(
                ball(dyadic(3, 0, 2), zero.clone()),
                ball(dyadic(1, 2, 1), zero.clone())
            ),
            3,
        )
        .unwrap(),
        ball(dyadic(5, 0, 3), dyadic(0, 0, 3)),
    );
    let origin = rectangle(zero_ball.clone(), zero_ball.clone());
    assert_eq!(
        call("num_mp_complex_exp", origin.clone(), 1).unwrap(),
        rectangle(one_ball.clone(), zero_ball),
    );
    assert!(call("num_mp_complex_ln", origin, 32).is_err());
    assert!(
        call(
            "num_mp_complex_ln",
            rectangle(ball(zero, one), one_ball),
            32
        )
        .is_ok()
    );

    let manifest = session.manifest();
    for method in ["abs", "exp", "ln"] {
        let op = manifest.op(&format!("num_mp_complex_{method}")).unwrap();
        assert_eq!(op.signature_source, SignatureSource::Declared);
        assert_eq!(op.dim_rule, Some(DimRuleInfo::Dimensionless));
        assert_eq!(
            op.rust_path,
            Some(format!("crate::kernels::mp::complex::MpComplex::{method}")),
        );
    }
}

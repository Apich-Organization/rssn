//! Rational references and branch cases for principal complex inverse functions.

use num_bigint::BigInt;
use num_rational::BigRational;
use rssn::kernels::mp::Ball;
use rssn::kernels::mp::MpComplex;
use rssn::kernels::mp::Real;
use rssn::kernels::mp::Round;

type Bounds = (BigRational, BigRational);

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

fn contains(
    ball: &Ball,
    (lower, upper): &Bounds,
) {
    let midpoint = exact(&ball.mid);
    let radius = exact(&ball.rad);
    assert!(
        &midpoint - &radius <= *lower && midpoint + radius >= *upper,
        "{ball:?}"
    );
}

fn exact_zero(ball: &Ball) {
    assert!(ball.mid.is_zero() && ball.rad.is_zero(), "{ball:?}");
}

/// Bounds independent atan/atanh series for |x| <= 2/3. The remainder
/// is at most (2/3)^513 / (1 - (2/3)^2), below 2^-160.
fn angle_series(
    x: BigRational,
    alternating: bool,
) -> Bounds {
    let mut power = x.clone();
    let mut square = &x * &x;
    if alternating {
        square = -square;
    }
    let mut sum = int(0);
    for k in 0..256 {
        sum += &power / int(2 * k + 1);
        power *= &square;
    }
    let tail = BigRational::new(1.into(), BigInt::from(1) << 160_u32);
    (&sum - &tail, sum + tail)
}

fn quarter_pi() -> Bounds {
    let (lo1, hi1) = angle_series(int(2).recip(), true);
    let (lo2, hi2) = angle_series(int(3).recip(), true);
    (lo1 + lo2, hi1 + hi2)
}

fn ln_reference(n: i64) -> Bounds {
    let atanh = angle_series(int(n - 1) / int(n + 1), false);
    (atanh.0 * int(2), atanh.1 * int(2))
}

fn scale(
    (lo, hi): Bounds,
    numerator: i64,
    denominator: i64,
) -> Bounds {
    let factor = int(numerator) / int(denominator);
    (lo * &factor, hi * factor)
}

#[test]
fn exact_axes_use_certified_real_angles() {
    let half = dyadic(1, -1, 1);
    let one = Real::one(1);
    let zero = Real::zero(1);
    let real = MpComplex::point(half.clone(), zero.clone());
    let asin = real.asin(64).unwrap();
    contains(&asin.re, &scale(quarter_pi(), 2, 3));
    exact_zero(&asin.im);
    let acos = real.acos(64).unwrap();
    contains(&acos.re, &scale(quarter_pi(), 4, 3));
    exact_zero(&acos.im);
    let atanh = real.atanh(64).unwrap();
    contains(&atanh.re, &scale(ln_reference(3), 1, 2));
    exact_zero(&atanh.im);

    let asinh = MpComplex::point(zero.clone(), half).asinh(64).unwrap();
    contains(&asinh.im, &scale(quarter_pi(), 2, 3));
    exact_zero(&asinh.re);
    let atanh_i = MpComplex::point(zero.clone(), one.clone())
        .atanh(64)
        .unwrap();
    contains(&atanh_i.im, &quarter_pi());
    exact_zero(&atanh_i.re);
    let atan = MpComplex::point(one, zero.clone()).atan(64).unwrap();
    contains(&atan.re, &quarter_pi());
    exact_zero(&atan.im);

    let origin = MpComplex::zero(1).atan(1).unwrap();
    exact_zero(&origin.re);
    exact_zero(&origin.im);
    let one = MpComplex::one(1).acosh(1).unwrap();
    exact_zero(&one.re);
    exact_zero(&one.im);
}

#[test]
fn off_axis_values_and_rectangle_images_use_correct_principal_branches() {
    let z = MpComplex::point(Real::one(1), Real::one(1));
    let atan = z.atan(64).unwrap();
    let (pi_lo, pi_hi) = quarter_pi();
    let (angle_lo, angle_hi) = scale(angle_series(int(2).recip(), true), 1, 2);
    contains(&atan.re, &(pi_lo + angle_lo, pi_hi + angle_hi));
    contains(&atan.im, &scale(ln_reference(5), 1, 4));

    // This rectangle contains 5/8 + i*3*sqrt(3)/8, which equals both
    // sin(pi/6 + i ln(2)) and cosh(ln(2) + i*pi/3).
    let sqrt3_floor = (BigInt::from(3) << 220_u32).sqrt();
    let imag_lo = Real::rounded(sqrt3_floor * 3, -113, 128, Round::NearestEven).unwrap();
    let imag = Ball::new(imag_lo, dyadic(3, -113, 2)).unwrap();
    let rectangle = MpComplex::new(Ball::point(dyadic(5, -3, 3)), imag);
    let asin = rectangle.asin(64).unwrap();
    contains(&asin.re, &scale(quarter_pi(), 2, 3));
    contains(&asin.im, &ln_reference(2));
    let acosh = rectangle.acosh(64).unwrap();
    contains(&acosh.re, &ln_reference(2));
    contains(&acosh.im, &scale(quarter_pi(), 4, 3));
}

#[test]
fn inverse_sine_uses_the_cancellation_free_logarithm() {
    let x = dyadic(1, 5000, 1);
    let value = MpComplex::point(Real::zero(1), x.clone()).asin(64).unwrap();
    let (base_lo, base_hi) = scale(ln_reference(2), 5001, 1);
    // 0 < asinh(x) - ln(2x) < x^-2 / 4.
    let correction = exact(&dyadic(1, -10002, 1));
    let upper = base_hi + correction;
    contains(&value.im, &(base_lo.clone(), upper.clone()));
    assert!(value.re.contains_zero());
    let reflected = MpComplex::point(x.neg(), Real::zero(1)).asinh(64).unwrap();
    contains(&reflected.re, &(-upper, -base_lo));
    assert!(reflected.im.contains_zero());
}

#[test]
fn branch_points_and_rectangles_have_deliberate_domain_behavior() {
    let one = MpComplex::one(1);
    assert!(one.asin(64).is_ok());
    assert!(one.acos(64).is_ok());
    assert!(one.atanh(64).is_err());
    assert!(MpComplex::i(1).atan(64).is_err());
    assert!(MpComplex::i(1).asinh(64).is_ok());

    let radius = dyadic(1, -3, 1);
    let strip = Ball::new(Real::zero(1), radius.clone()).unwrap();
    let right_cut = MpComplex::new(Ball::point(dyadic(1, 1, 1)), strip.clone());
    assert!(right_cut.asin(64).is_err());
    assert!(right_cut.acos(64).is_err());
    assert!(right_cut.atanh(64).is_err());
    let imaginary_cut = MpComplex::new(strip.clone(), Ball::point(dyadic(1, 1, 1)));
    assert!(imaginary_cut.atan(64).is_err());
    assert!(imaginary_cut.asinh(64).is_err());
    let crosses_acosh_cut = MpComplex::new(Ball::point(Real::zero(1)), strip);
    assert!(crosses_acosh_cut.acosh(64).is_err());
}

#[test]
fn session_declares_all_principal_inverse_rectangle_contracts() {
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
    let dyadic = |mantissa: i64| {
        record(&[
            ("mantissa", Value::int(mantissa)),
            ("exponent", Value::int(0)),
            ("precision", Value::int(1)),
        ])
    };
    let zero = dyadic(0);
    let ball = |midpoint| record(&[("midpoint", midpoint), ("radius", zero.clone())]);
    let zero_ball = ball(zero.clone());
    let origin = record(&[("real", zero_ball.clone()), ("imag", zero_ball.clone())]);
    let one = record(&[("real", ball(dyadic(1))), ("imag", zero_ball.clone())]);
    let imaginary_unit = record(&[("real", zero_ball), ("imag", ball(dyadic(1)))]);
    let session = Session::new();
    let ctx = CallCtx::new();
    let manifest = session.manifest();
    for method in ["asin", "acos", "atan", "asinh", "acosh", "atanh"] {
        let name = format!("num_mp_complex_{method}");
        let op = manifest.op(&name).unwrap();
        assert_eq!(op.signature_source, SignatureSource::Declared);
        assert_eq!(op.dim_rule, Some(DimRuleInfo::Dimensionless));
        assert_eq!(
            op.rust_path,
            Some(format!("crate::kernels::mp::complex::MpComplex::{method}")),
        );
        let operand = if matches!(method, "acos" | "acosh") {
            &one
        } else {
            &origin
        };
        let actual = session
            .call_value(&name, &[operand.clone(), Value::int(1)], &ctx)
            .unwrap();
        assert_eq!(actual, origin, "{method}");
    }
    assert!(
        session
            .call_value(
                "num_mp_complex_atan",
                &[imaginary_unit, Value::int(32)],
                &ctx
            )
            .is_err()
    );
}

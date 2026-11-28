//! Exact and rational reference checks for principal complex roots and powers.

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

fn bounds(ball: &Ball) -> (BigRational, BigRational) {
    let mid = exact(&ball.mid);
    let rad = exact(&ball.rad);
    (&mid - &rad, mid + rad)
}

fn contains(
    ball: &Ball,
    left: &BigRational,
    right: &BigRational,
) {
    let (lower, upper) = bounds(ball);
    assert!(&lower <= left && &upper >= right, "{ball:?}");
}

fn point(
    real: i64,
    imag: i64,
) -> MpComplex {
    MpComplex::point(
        Real::from_i64(real, 8).unwrap(),
        Real::from_i64(imag, 8).unwrap(),
    )
}

fn is_exact(
    ball: &Ball,
    expected: &BigRational,
) {
    assert_eq!(&exact(&ball.mid), expected);
    assert!(ball.rad.is_zero());
}

fn sqrt2_bounds() -> (BigRational, BigRational) {
    let scale = BigInt::from(1) << 128_u32;
    let floor = (BigInt::from(2) << 256_u32).sqrt();
    (
        BigRational::new(floor.clone(), scale.clone()),
        BigRational::new(floor + 1, scale),
    )
}

#[test]
fn principal_square_root_matches_exact_squares_in_both_half_planes() {
    for (x, y, re, im) in [(3, 4, 2, 1), (-3, 4, 1, 2), (-3, -4, 1, -2), (-4, 0, 0, 2)] {
        let actual = point(x, y).sqrt(64).unwrap();
        is_exact(&actual.re, &int(re));
        is_exact(&actual.im, &int(im));
    }
    let centered_on_cut = MpComplex::new(
        Ball::new(dyadic(-5, -1, 3), dyadic(3, -1, 2)).unwrap(),
        Ball::new(Real::zero(1), dyadic(1, -2, 1)).unwrap(),
    );
    let result = centered_on_cut.sqrt(32).unwrap();
    contains(&result.im, &int(-2), &int(2));
    let starts_on_cut = MpComplex::new(
        Ball::new(dyadic(-5, -1, 3), dyadic(3, -1, 2)).unwrap(),
        Ball::new(dyadic(1, -3, 1), dyadic(1, -3, 1)).unwrap(),
    );
    let one_sided = starts_on_cut.sqrt(32).unwrap();
    let (low, high) = bounds(&one_sided.im);
    assert!(low > int(0) && low <= int(1) && high >= int(2));
}

#[test]
fn square_root_preserves_scaled_norms_outside_the_real_exponent_range() {
    let tiny = MpComplex::point(dyadic(1, -MAX_EXP, 1), Real::zero(1));
    let root = tiny.sqrt(64).unwrap();
    is_exact(&root.re, &exact(&dyadic(1, -MAX_EXP / 2, 1)));
    is_exact(&root.im, &int(0));

    let largest = Real::rounded(
        (BigInt::from(1) << MAX_BITS) - 1_u32,
        MAX_EXP,
        MAX_BITS,
        Round::NearestEven,
    )
    .unwrap();
    let huge = MpComplex::point(largest.clone(), largest.clone());
    assert!(huge.abs(64).is_err());
    let root = huge.sqrt(64).unwrap();
    let (sqrt2_lower, sqrt2_upper) = sqrt2_bounds();
    let magnitude = exact(&largest);
    for (component, sign) in [(&root.re, int(1)), (&root.im, int(-1))] {
        let (lower, upper) = bounds(component);
        assert!(lower > int(0));
        assert!(lower.pow(2) * int(2) <= (&sqrt2_lower + &sign) * &magnitude);
        assert!(upper.pow(2) * int(2) >= (&sqrt2_upper + sign) * &magnitude);
    }
}

#[test]
fn small_coordinate_avoids_subtracting_nearly_equal_norms() {
    let negative = MpComplex::point(dyadic(-1, 4000, 1), dyadic(1, -4000, 1));
    let actual = negative.sqrt(64).unwrap();
    // u = y / (2v), where 2^2000 <= v < 2^2000 * (1 + 2^-256).
    let power_of_two = exact(&dyadic(1, -6001, 1));
    let epsilon = exact(&dyadic(1, -256, 1));
    let lower = &power_of_two * (int(1) - epsilon);
    let (actual_low, _) = bounds(&actual.re);
    assert!(actual_low > int(0));
    contains(&actual.re, &lower, &power_of_two);
}

#[test]
fn complex_power_checks_exact_integer_range_and_preserves_principal_values() {
    let too_large = MpComplex::point(dyadic(1, 64, 63), Real::zero(1));
    let error = point(2, 0).pow(&too_large, 32).unwrap_err();
    assert!(error.contains("signed 64-bit range"));
    let smallest = MpComplex::point(Real::from_i64(i64::MIN, 1).unwrap(), Real::zero(1));
    let identity = point(-1, 0).pow(&smallest, 32).unwrap();
    is_exact(&identity.re, &int(1));
    is_exact(&identity.im, &int(0));

    let half = MpComplex::point(dyadic(1, -1, 1), Real::zero(1));
    let result = point(0, 1).pow(&half, 64).unwrap();
    let (lo, hi) = sqrt2_bounds();
    let (lo, hi) = (lo / int(2), hi / int(2));
    contains(&result.re, &lo, &hi);
    contains(&result.im, &lo, &hi);
    assert!(point(-1, 0).pow(&half, 64).is_err());
}

#[test]
fn session_declares_principal_root_and_power_contracts() {
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
    let zero = dyadic(0, 0, 1);
    let ball = |midpoint| record(&[("midpoint", midpoint), ("radius", zero.clone())]);
    let rectangle = |real, imag| record(&[("real", ball(real)), ("imag", ball(imag))]);
    let session = Session::new();
    let ctx = CallCtx::new();
    let neg_four = rectangle(dyadic(-1, 2, 1), zero.clone());
    assert_eq!(
        session
            .call_value("num_mp_complex_sqrt", &[neg_four, Value::int(1)], &ctx)
            .unwrap(),
        rectangle(zero.clone(), dyadic(1, 1, 1)),
    );
    let two = rectangle(dyadic(1, 1, 1), zero.clone());
    let huge_integer = rectangle(dyadic(1 << 62, 2, 63), zero.clone());
    assert!(
        session
            .call_value(
                "num_mp_complex_pow",
                &[two, huge_integer, Value::int(16)],
                &ctx
            )
            .is_err()
    );
    let manifest = session.manifest();
    for method in ["sqrt", "pow"] {
        let op = manifest.op(&format!("num_mp_complex_{method}")).unwrap();
        assert_eq!(op.signature_source, SignatureSource::Declared);
        assert_eq!(op.dim_rule, Some(DimRuleInfo::Dimensionless));
        assert_eq!(
            op.rust_path,
            Some(format!("crate::kernels::mp::complex::MpComplex::{method}")),
        );
    }
}

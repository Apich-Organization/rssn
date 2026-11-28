//! High-precision rational references for complex rectangle trigonometry.

use num_bigint::BigInt;
use num_rational::BigRational;
use rssn::kernels::mp::Ball;
use rssn::kernels::mp::MpComplex;
use rssn::kernels::mp::Real;
use rssn::kernels::mp::Round;

type Bounds = (BigRational, BigRational);
type Operation = fn(&MpComplex, u32) -> Result<MpComplex, String>;

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

fn encloses(
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

fn sum(
    (left_lo, left_hi): Bounds,
    (right_lo, right_hi): Bounds,
) -> Bounds {
    (left_lo + right_lo, left_hi + right_hi)
}

fn negate((lower, upper): Bounds) -> Bounds {
    (-upper, -lower)
}

fn product(
    left: &Bounds,
    right: &Bounds,
) -> Bounds {
    let values = [
        &left.0 * &right.0,
        &left.0 * &right.1,
        &left.1 * &right.0,
        &left.1 * &right.1,
    ];
    (
        values.iter().min().unwrap().clone(),
        values.iter().max().unwrap().clone(),
    )
}

fn divide(
    numerator: &Bounds,
    (lower, upper): Bounds,
) -> Bounds {
    assert!(lower > int(0));
    product(numerator, &(upper.recip(), lower.recip()))
}

/// Rational Taylor reference for |x| <= 1. The first omitted term is
/// below 1/64!, and every later magnitude is at most half its predecessor.
fn series(
    x: &BigRational,
    odd: bool,
    alternating: bool,
) -> Bounds {
    let mut term = if odd { x.clone() } else { int(1) };
    let mut sum = term.clone();
    let square = if alternating { -x * x } else { x * x };
    let offset = i64::from(odd);
    for k in 1..32 {
        term = term * &square / int((2 * k + offset - 1) * (2 * k + offset));
        sum += &term;
    }
    let remainder = BigRational::new(1.into(), BigInt::from(1) << 150_u32);
    (&sum - &remainder, sum + remainder)
}

fn reference(
    x: &BigRational,
    y: &BigRational,
) -> [(Bounds, Bounds); 6] {
    let (sin_x, cos_x) = (series(x, true, true), series(x, false, true));
    let (sinh_x, cosh_x) = (series(x, true, false), series(x, false, false));
    let (sin_y, cos_y) = (series(y, true, true), series(y, false, true));
    let (sinh_y, cosh_y) = (series(y, true, false), series(y, false, false));
    let dx = x * int(2);
    let dy = y * int(2);
    let tan_denom = sum(series(&dx, false, true), series(&dy, false, false));
    let tanh_denom = sum(series(&dx, false, false), series(&dy, false, true));
    [
        (product(&sin_x, &cosh_y), product(&cos_x, &sinh_y)),
        (product(&cos_x, &cosh_y), negate(product(&sin_x, &sinh_y))),
        (
            divide(&series(&dx, true, true), tan_denom.clone()),
            divide(&series(&dy, true, false), tan_denom),
        ),
        (product(&sinh_x, &cos_y), product(&cosh_x, &sin_y)),
        (product(&cosh_x, &cos_y), product(&sinh_x, &sin_y)),
        (
            divide(&series(&dx, true, false), tanh_denom.clone()),
            divide(&series(&dy, true, true), tanh_denom),
        ),
    ]
}

fn operations() -> [Operation; 6] {
    [
        MpComplex::sin,
        MpComplex::cos,
        MpComplex::tan,
        MpComplex::sinh,
        MpComplex::cosh,
        MpComplex::tanh,
    ]
}

#[test]
fn complex_trigonometry_encloses_independent_point_and_rectangle_references() {
    let x = dyadic(1, -2, 1);
    let y = dyadic(-1, -1, 1);
    let pt = MpComplex::point(x.clone(), y.clone());
    for (operation, (real, imag)) in operations()
        .into_iter()
        .zip(reference(&exact(&x), &exact(&y)))
    {
        let actual = operation(&pt, 64).unwrap();
        encloses(&actual.re, &real);
        encloses(&actual.im, &imag);
    }
    let rad = dyadic(1, -5, 1);
    let rectangle = MpComplex::new(
        Ball::new(x.clone(), rad.clone()).unwrap(),
        Ball::new(y.clone(), rad.clone()).unwrap(),
    );
    let (x, y, rad) = (exact(&x), exact(&y), exact(&rad));
    let corners = [
        reference(&(&x - &rad), &(&y - &rad)),
        reference(&(&x - &rad), &(&y + &rad)),
        reference(&(&x + &rad), &(&y - &rad)),
        reference(&(&x + &rad), &(&y + &rad)),
    ];
    for (index, operation) in operations().into_iter().enumerate() {
        let actual = operation(&rectangle, 64).unwrap();
        for corner in &corners {
            encloses(&actual.re, &corner[index].0);
            encloses(&actual.im, &corner[index].1);
        }
    }
}

#[test]
fn origin_maps_to_exact_integer_points_at_one_bit() {
    for (index, operation) in operations().into_iter().enumerate() {
        let value = operation(&MpComplex::zero(1), 1).unwrap();
        let expected = i64::from(matches!(index, 1 | 4));
        assert_eq!(exact(&value.re.mid), int(expected));
        assert!(value.re.rad.is_zero() && value.im.mid.is_zero() && value.im.rad.is_zero());
    }
}

#[test]
fn tangent_axes_saturate_without_overflow_and_reject_poles() {
    let large = dyadic(1, 20, 1);
    let near_one = (
        int(1) - BigRational::new(1.into(), BigInt::from(1) << 1024_u32),
        int(1),
    );
    let tanh = MpComplex::point(large.clone(), Real::zero(1))
        .tanh(32)
        .unwrap();
    encloses(&tanh.re, &near_one);
    assert!(tanh.im.mid.is_zero() && tanh.im.rad.is_zero());
    let tan = MpComplex::point(Real::zero(1), large).tan(32).unwrap();
    encloses(&tan.im, &near_one);
    assert!(tan.re.mid.is_zero() && tan.re.rad.is_zero());

    let pole = Ball::new(dyadic(3, -1, 2), dyadic(1, -2, 1)).unwrap();
    let zero = Ball::point(Real::zero(1));
    assert!(MpComplex::from_real(pole.clone()).tan(32).is_err());
    assert!(MpComplex::new(zero, pole).tanh(32).is_err());
}

#[test]
fn session_declares_all_six_complex_trig_contracts() {
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
    let dyadic = |mantissa: i64, exponent: i64| {
        record(&[
            ("mantissa", Value::int(mantissa)),
            ("exponent", Value::int(exponent)),
            ("precision", Value::int(1)),
        ])
    };
    let zero = dyadic(0, 0);
    let ball = |midpoint| record(&[("midpoint", midpoint), ("radius", zero.clone())]);
    let rectangle = |real, imag| record(&[("real", real), ("imag", imag)]);
    let zero_ball = ball(zero.clone());
    let one_ball = ball(dyadic(1, 0));
    let origin = rectangle(zero_ball.clone(), zero_ball.clone());
    let unit = rectangle(one_ball, zero_ball.clone());
    let session = Session::new();
    let ctx = CallCtx::new();
    let manifest = session.manifest();
    for method in ["sin", "cos", "tan", "sinh", "cosh", "tanh"] {
        let name = format!("num_mp_complex_{method}");
        let op = manifest.op(&name).unwrap();
        assert_eq!(op.signature_source, SignatureSource::Declared);
        assert_eq!(op.dim_rule, Some(DimRuleInfo::Dimensionless));
        assert_eq!(
            op.rust_path,
            Some(format!("crate::kernels::mp::complex::MpComplex::{method}")),
        );
        let expected = if matches!(method, "cos" | "cosh") {
            &unit
        } else {
            &origin
        };
        let actual = session
            .call_value(&name, &[origin.clone(), Value::int(1)], &ctx)
            .unwrap();
        assert_eq!(&actual, expected, "{method}");
    }
    let tall = rectangle(zero_ball, ball(dyadic(1, 20)));
    assert!(
        session
            .call_value("num_mp_complex_tan", &[tall, Value::int(32)], &ctx)
            .is_ok()
    );
}

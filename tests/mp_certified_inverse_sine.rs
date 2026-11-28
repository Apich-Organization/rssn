//! Independent reference checks for certified real inverse sine and cosine.

use num_bigint::BigInt;
use num_rational::BigRational;
use rssn::kernels::mp::Ball;
use rssn::kernels::mp::MAX_BITS;
use rssn::kernels::mp::Real;
use rssn::kernels::mp::Round;
use rssn::kernels::mp::funcs;
use rssn::kernels::mp::parse;

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
    lower: &BigRational,
    upper: &BigRational,
) {
    let (lo, hi) = endpoints(ball);
    assert!(
        &lo <= lower && &hi >= upper,
        "{ball:?}: lower covered {}, upper covered {}",
        &lo <= lower,
        &hi >= upper
    );
}

fn quarter_pi() -> (BigRational, BigRational) {
    let reference = |denominator| {
        let ratio = int(denominator).recip();
        let square = &ratio * &ratio;
        let mut power = ratio;
        let mut sum = int(0);
        for k in 0..120 {
            let term = &power / int(2 * k + 1);
            if k % 2 == 0 {
                sum += &term;
            } else {
                sum -= &term;
            }
            power *= &square;
        }
        let high = &sum + power / int(241);
        (sum, high)
    };
    // atan(1/2) + atan(1/3) = π/4, independent of Machin's formula.
    let (a, b) = reference(2);
    let (c, d) = reference(3);
    (a + c, b + d)
}

#[test]
fn certified_endpoints_and_signs_contain_independent_pi_values() {
    let (pi_low, pi_high) = quarter_pi();
    let half = dyadic(1, -1, 1);
    let one = Real::one(1);
    let zero = Real::zero(1);
    for prec in [1, 4, 32, 128] {
        encloses(
            &funcs::asin(&one, prec).unwrap(),
            &(int(2) * &pi_low),
            &(int(2) * &pi_high),
        );
        encloses(
            &funcs::asin(&one.neg(), prec).unwrap(),
            &(-int(2) * &pi_high),
            &(-int(2) * &pi_low),
        );
        encloses(
            &funcs::acos(&one.neg(), prec).unwrap(),
            &(int(4) * &pi_low),
            &(int(4) * &pi_high),
        );
        encloses(
            &funcs::acos(&zero, prec).unwrap(),
            &(int(2) * &pi_low),
            &(int(2) * &pi_high),
        );
        encloses(
            &funcs::asin(&half, prec).unwrap(),
            &(int(2) * &pi_low / int(3)),
            &(int(2) * &pi_high / int(3)),
        );
        encloses(
            &funcs::acos(&half, prec).unwrap(),
            &(int(4) * &pi_low / int(3)),
            &(int(4) * &pi_high / int(3)),
        );
        encloses(
            &funcs::acos(&half.neg(), prec).unwrap(),
            &(int(8) * &pi_low / int(3)),
            &(int(8) * &pi_high / int(3)),
        );
        let at_one = funcs::acos(&one, prec).unwrap();
        assert!(at_one.mid.is_zero() && at_one.rad.is_zero());
        let at_zero = funcs::asin(&zero, prec).unwrap();
        assert!(at_zero.mid.is_zero() && at_zero.rad.is_zero());
    }
}

#[test]
fn nearly_unit_dyadics_retain_tiny_acos_and_accept_valid_inputs() {
    let near = Real::rounded(
        (BigInt::from(1) << 128_u32) - 1_u32,
        -128,
        128,
        Round::NearestEven,
    )
    .unwrap();
    let out = funcs::acos(&near, 64).unwrap();
    assert_eq!(out.mid.floor_log2(), Some(-64));
    let (low, high) = endpoints(&out);
    assert!(low > int(0) && high < int(1));
    // Independent Taylor bounds, using monotonicity of cos near the origin.
    let cos_bounds = |theta: &BigRational| {
        let square = theta * theta;
        let mut power = int(1);
        let mut sum = power.clone();
        for k in 1..=6 {
            power = power * &square / int((2 * k - 1) * (2 * k));
            if k % 2 == 0 {
                sum += &power;
            } else {
                sum -= &power;
            }
        }
        let lower = &sum - power * square / int(13 * 14);
        (lower, sum)
    };
    assert!(cos_bounds(&low).0 >= exact(&near));
    assert!(cos_bounds(&high).1 <= exact(&near));

    let at_max = Real::rounded(
        (BigInt::from(1) << MAX_BITS) - 1_u32,
        -i64::from(MAX_BITS),
        MAX_BITS,
        Round::NearestEven,
    )
    .unwrap();
    assert_eq!(
        funcs::acos(&at_max, 48).unwrap().mid.floor_log2(),
        Some(-i64::from(MAX_BITS / 2))
    );
    assert!(funcs::asin(&at_max, 4).is_ok());
    assert!(funcs::acos(&at_max.neg(), 4).is_ok());
}

#[test]
fn closed_intervals_preserve_boundary_order_and_reject_outside_domain() {
    let value = Ball::new(dyadic(1, -2, 1), dyadic(3, -2, 2)).unwrap();
    let (pi_low, pi_high) = quarter_pi();
    encloses(
        &funcs::asin_ball(&value, 64).unwrap(),
        &(-int(2) * &pi_high / int(3)),
        &(int(2) * &pi_high),
    );
    encloses(
        &funcs::acos_ball(&value, 64).unwrap(),
        &int(0),
        &(int(8) * &pi_high / int(3)),
    );
    let outside = Ball::new(dyadic(3, -2, 2), dyadic(1, -1, 1)).unwrap();
    assert!(funcs::asin_ball(&outside, 64).is_err());
    assert!(funcs::acos_ball(&outside.neg(), 64).is_err());
    assert!(funcs::asin(&dyadic(5, -2, 3), 64).is_err());
    assert!(funcs::acos(&dyadic(-5, -2, 3), 64).is_err());
    let text = parse::evaluate("acos(-1)", 64).unwrap();
    encloses(&text, &(int(4) * pi_low), &(int(4) * pi_high));
}

#[test]
fn native_session_keeps_declared_endpoint_and_interval_contracts() {
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
                .map(|(n, v)| ((*n).into(), v.clone()))
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
    let session = Session::new();
    let context = CallCtx::new();
    let call = |name: &str, args: &[Value]| session.call_value(name, args, &context);
    let one = dyadic(1, 0, 1);
    let zero = dyadic(0, 0, 1);
    let entire = record(&[("midpoint", zero.clone()), ("radius", one.clone())]);
    let outside = record(&[("midpoint", one.clone()), ("radius", one.clone())]);
    let Value::Record(asin_one) = call("num_mp_real_asin", &[one.clone(), Value::int(8)]).unwrap()
    else {
        panic!("expected a dyadic ball");
    };
    assert_ne!(asin_one.get("radius"), Some(&dyadic(0, 0, 8)));
    assert_eq!(
        call("num_mp_real_acos", &[one, Value::int(1)]).unwrap(),
        record(&[("midpoint", zero.clone()), ("radius", zero)]),
    );
    assert!(call("num_mp_ball_asin", &[entire.clone(), Value::int(16)]).is_ok());
    assert!(call("num_mp_ball_acos", &[entire, Value::int(16)]).is_ok());
    assert!(call("num_mp_ball_asin", &[outside.clone(), Value::int(16)]).is_err());
    assert!(call("num_mp_ball_acos", &[outside, Value::int(16)]).is_err());

    let manifest = session.manifest();
    for (name, source) in [
        ("real_asin", "asin"),
        ("ball_asin", "asin_ball"),
        ("real_acos", "acos"),
        ("ball_acos", "acos_ball"),
    ] {
        let op = manifest.op(&format!("num_mp_{name}")).unwrap();
        assert_eq!(op.signature_source, SignatureSource::Declared);
        assert_eq!(op.dim_rule, Some(DimRuleInfo::Dimensionless));
        assert_eq!(
            op.rust_path,
            Some(format!("crate::kernels::mp::funcs::{source}"))
        );
    }
}
